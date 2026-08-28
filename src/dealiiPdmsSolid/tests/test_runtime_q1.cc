#define DEALII_PDMS_NO_MAIN
#include "../dealiiPdmsSolid.cc"

#include <deal.II/base/mpi.h>
#include <cstdlib>
#include <fstream>
#include <random>

using namespace dealii;

static void declare_parameters(ParameterHandler &prm)
{
  prm.declare_entry("mesh","",Patterns::Anything());
  prm.declare_entry("interface transfer","legacyNearestNeighbour",Patterns::Anything());
  prm.declare_entry("interface boundary","4",Patterns::Integer(0));
  prm.declare_entry("clamped boundary 1","1",Patterns::Integer(0));
  prm.declare_entry("clamped boundary 2","2",Patterns::Integer(0));
  prm.declare_entry("clamped boundary 3","3",Patterns::Integer(0));
  prm.declare_entry("state input","",Patterns::Anything());
}

static std::vector<runtime_q1::Query<3>> parse_queries(const std::string &contents)
{
  std::istringstream input(contents);
  std::vector<runtime_q1::Query<3>> queries; std::string line;
  bool canonical_manifest=false;
  while (std::getline(input,line))
    {
      if (line.rfind("#",0)==0) continue;
      canonical_manifest=line.find("globalFluidPointId")!=std::string::npos;
      break;
    }
  std::uint64_t sequential_id=0;
  while (std::getline(input,line))
    {
      std::replace(line.begin(),line.end(),',',' '); std::istringstream row(line);
      runtime_q1::Query<3> query;
      if (canonical_manifest)
        { unsigned int owner; if (row>>query.id>>query.point[0]>>query.point[1]>>query.point[2]>>owner) queries.push_back(query); }
      else
        { query.id=sequential_id++; if (row>>query.point[0]>>query.point[1]>>query.point[2]) queries.push_back(query); }
    }
  return queries;
}

static std::vector<runtime_q1::Query<3>> read_queries_collectively(const std::string &path)
{
  const unsigned int rank=Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  std::string contents;
  if (rank==0)
    {
      std::ifstream input(path,std::ios::binary); AssertThrow(input,ExcMessage("cannot open query fixture"));
      contents.assign(std::istreambuf_iterator<char>(input),std::istreambuf_iterator<char>());
    }
  std::uint64_t size=contents.size();
  MPI_Bcast(&size,1,MPI_UINT64_T,0,MPI_COMM_WORLD);
  AssertThrow(size<=256ULL*1024ULL*1024ULL,ExcMessage("query manifest exceeds diagnostic limit"));
  if (rank!=0) contents.resize(size);
  if (size) MPI_Bcast(contents.data(),static_cast<int>(size),MPI_CHAR,0,MPI_COMM_WORLD);
  return parse_queries(contents);
}

int main(int argc,char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc,argv,1);
  try
    {
      AssertThrow(argc==3,ExcMessage("usage: testRuntimeQ1 solid.msh query.csv"));
      ParameterHandler prm; declare_parameters(prm); prm.set("mesh",argv[1]);
      PdmsSolid solid(prm); solid.initialize_runtime();
      std::string dof_manifest_hash; solid.dof_manifest_payload(dof_manifest_hash);
      const auto queries=read_queries_collectively(argv[2]);
      const char *tolerance_text=std::getenv("RUNTIME_Q1_REL_TOL");
      const double tolerance=tolerance_text ? std::stod(tolerance_text) : 1e-10;
      const auto op=solid.build_canonical_interface_q1_operator(queries,tolerance);
      AssertThrow(op.rows.size()==3*queries.size(),ExcMessage("RuntimeQ1 row count"));
      AssertThrow(op.maximum_residual<=1e-12,ExcMessage("RuntimeQ1 geometric residual"));

      std::map<types::global_dof_index,Point<3>> support;
      std::map<std::uint64_t,Point<3>> query_points;
      for (const auto &query:queries) query_points.emplace(query.id,query.point);
      DoFTools::map_dofs_to_support_points(MappingQ1<3>(),solid.runtime_dofs(),support);
      double partition=0,linear=0,rotation=0,work_difference=0;
      Tensor<1,3> omega; omega[0]=.7; omega[1]=-.2; omega[2]=.4;
      std::mt19937_64 generator(20260826); std::uniform_real_distribution<double> random(-1,1);
      Vector<double> virtual_dof(solid.runtime_dof_count());
      std::vector<double> transpose(solid.runtime_dof_count(),0.0),transpose_repeat(solid.runtime_dof_count(),0.0);
      for (auto &value:virtual_dof) value=random(generator);
      solid.runtime_constraints().distribute(virtual_dof);
      double point_work=0,dof_work=0; Tensor<1,3> point_force,solid_force,point_moment,solid_moment;
      std::size_t nonzeros=0;
      for (const auto &row:op.rows)
        {
          double sum=0,coordinate=0,rigid=0,virtual_value=0;
          const auto &point=query_points.at(row.query_id);
          for (const auto &entry:row.entries)
            {
              sum+=entry.weight; coordinate+=entry.weight*support.at(entry.dof)[row.component];
              Tensor<1,3> x; for (unsigned int d=0;d<3;++d) x[d]=support.at(entry.dof)[d];
              rigid+=entry.weight*cross_product_3d(omega,x)[row.component];
              virtual_value+=entry.weight*virtual_dof[entry.dof];
            }
          partition=std::max(partition,std::abs(sum-1));
          linear=std::max(linear,std::abs(coordinate-point[row.component]));
          Tensor<1,3> p; for (unsigned int d=0;d<3;++d) p[d]=point[d];
          rotation=std::max(rotation,std::abs(rigid-cross_product_3d(omega,p)[row.component]));
          const double force=random(generator); point_work+=virtual_value*force;
          nonzeros+=row.entries.size(); point_force[row.component]+=force;
          Tensor<1,3> point_load; point_load[row.component]=force;
          point_moment+=cross_product_3d(p,point_load);
          for (const auto &entry:row.entries)
            {
              const double contribution=entry.weight*force;
              dof_work+=virtual_dof[entry.dof]*contribution; transpose[entry.dof]+=contribution;
              solid_force[row.component]+=contribution;
              Tensor<1,3> x,load; for (unsigned int d=0;d<3;++d) x[d]=support.at(entry.dof)[d];
              load[row.component]=contribution; solid_moment+=cross_product_3d(x,load);
            }
        }
      transpose_repeat=transpose;
      work_difference=std::abs(point_work-dof_work)/std::max({std::abs(point_work),std::abs(dof_work),1e-30});
      const double force_error=(point_force-solid_force).norm()/std::max(point_force.norm(),1e-30);
      const double moment_error=(point_moment-solid_moment).norm()/std::max(point_moment.norm(),1e-30);
      AssertThrow(partition<=1e-14 && linear<=1e-12 && rotation<=1e-12,
                  ExcMessage("RuntimeQ1 reproduction failure"));
      AssertThrow(work_difference<=1e-12,ExcMessage("RuntimeQ1 dual work failure"));
      AssertThrow(force_error<=1e-12 && moment_error<=1e-12,ExcMessage("RuntimeQ1 force/moment failure"));
      AssertThrow(transpose==transpose_repeat,ExcMessage("RuntimeQ1 transpose repetition failure"));
      const auto min_graph=Utilities::MPI::min(runtime_q1::fnv(op.graph_hash),MPI_COMM_WORLD);
      const auto max_graph=Utilities::MPI::max(runtime_q1::fnv(op.graph_hash),MPI_COMM_WORLD);
      AssertThrow(min_graph==max_graph,ExcMessage("RuntimeQ1 MPI graph mismatch"));
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
        std::cout << std::setprecision(17) << "RUNTIME_Q1 PASS rows=" << op.rows.size()
          << " graphHash=" << op.graph_hash << " weightsHash=" << op.weights_hash
          << " queryHash=" << op.query_manifest_hash << " supportHash=" << op.support_hash
          << " dofManifestHash=" << dof_manifest_hash
          << " geometryHash=" << op.geometry_hash << " nonzeros=" << nonzeros
          << " maxResidual=" << op.maximum_residual
          << " partition=" << partition << " linear=" << linear << " rotation=" << rotation
          << " work=" << work_difference << " force=" << force_error << " moment=" << moment_error
          << " ambiguous=" << op.multiple_candidates << " snapped=" << op.snapped_coordinates
          << " maxSnap=" << op.maximum_reference_snap
          << " faceInterior=" << op.face_interior_points << " edge=" << op.edge_points
          << " vertex=" << op.vertex_points << " equivalentCandidates=" << op.equivalent_candidate_count
          << " maxPrunedMass=" << op.maximum_pruned_weight_mass << std::endl;
    }
  catch (const std::exception &error)
    { std::cerr << "RUNTIME_Q1 FAIL: " << error.what() << std::endl; return 1; }
}
