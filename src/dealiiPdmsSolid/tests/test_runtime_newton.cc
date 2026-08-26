#define DEALII_PDMS_NO_MAIN
#include "../dealiiPdmsSolid.cc"

#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <iostream>

using namespace dealii;

static void declare_runtime_parameters(ParameterHandler &prm)
{
  prm.declare_entry("mesh","",Patterns::Anything());
  prm.declare_entry("state input","",Patterns::Anything());
  prm.declare_entry("interface transfer","dualConservative",Patterns::Anything());
  prm.declare_entry("interface manifest","",Patterns::Anything());
  prm.declare_entry("interface graph hash","",Patterns::Anything());
  prm.declare_entry("interface weights hash","weights-v1",Patterns::Anything());
  prm.declare_entry("interface dof manifest hash","dofs-v1",Patterns::Anything());
  prm.declare_entry("interface time index","7",Patterns::Integer(0));
  prm.declare_entry("interface outer corrector","3",Patterns::Integer(0));
  prm.declare_entry("interface operator version","11",Patterns::Integer(0));
  prm.declare_entry("interface z version","5",Patterns::Integer(0));
  prm.declare_entry("clamped boundary 1","1",Patterns::Integer(0));
  prm.declare_entry("clamped boundary 2","2",Patterns::Integer(0));
  prm.declare_entry("clamped boundary 3","3",Patterns::Integer(0));
  prm.declare_entry("rho","970",Patterns::Double(0));
  prm.declare_entry("mu","778200",Patterns::Double(0));
  prm.declare_entry("bulk modulus","3631600",Patterns::Double(0));
  prm.declare_entry("delta t","1e-7",Patterns::Double(0));
  prm.declare_entry("solid impedance","67350",Patterns::Double(0));
  prm.declare_entry("time integration","backwardEuler",Patterns::Anything());
  prm.declare_entry("minimum jacobian","1e-8",Patterns::Double(0));
  prm.declare_entry("interface boundary","4",Patterns::Integer(0));
}

int main(int argc,char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc,argv,1);
  try
    {
      AssertThrow(argc==3,ExcMessage("usage: testRuntimeNewton solid.msh manifest"));
      const auto manifest=interface_sparsity::read_manifest(argv[2]);
      ParameterHandler prm; declare_runtime_parameters(prm);
      prm.set("mesh",argv[1]); prm.set("interface manifest",argv[2]);
      prm.set("interface graph hash",manifest.hash_graph);
      PdmsSolid solid(prm); solid.initialize_runtime();
      const auto &node_dofs=solid.runtime_node_dofs();
      const auto &constraints=solid.runtime_constraints();
      types::global_dof_index free_dof=numbers::invalid_dof_index;
      for (const auto &node:node_dofs)
        for (unsigned int c=0;c<3 && free_dof==numbers::invalid_dof_index;++c)
          if (!constraints.is_constrained(node[c])) free_dof=node[c];
      AssertThrow(free_dof!=numbers::invalid_dof_index,ExcMessage("no free interface DoF"));
      unsigned int component=0;
      for (const auto &node:node_dofs) for (unsigned int c=0;c<3;++c)
        if (node[c]==free_dof) component=c;

      generalized_interface::Stamp stamp;
      stamp.time_index=7; stamp.outer_corrector=3; stamp.operator_version=11; stamp.z_version=5;
      stamp.hash_graph=manifest.hash_graph; stamp.hash_weights="weights-v1";
      stamp.dof_manifest_hash="dofs-v1"; stamp.valid=true;
      auto messages=[&](double displacement,double force0,double tangent_value)
        {
          generalized_interface::ForceMessage force; force.stamp=stamp;
          force.values.push_back({free_dof,component,force0-tangent_value*displacement});
          force.checksum=generalized_interface::checksum(force);
          generalized_interface::TangentMessage tangent; tangent.stamp=stamp;
          tangent.values.push_back({free_dof,free_dof,tangent_value});
          tangent.checksum=generalized_interface::checksum(tangent);
          return std::make_pair(force,tangent);
        };
      const double force0=1.0,tangent_value=100.0;
      auto data=messages(0,force0,tangent_value); solid.set_generalized_interface_data(data.first,data.second);
      solid.assemble_runtime({});
      Vector<double> rhs0=solid.runtime_rhs(); SparseMatrix<double> matrix0;
      matrix0.reinit(solid.runtime_sparsity()); matrix0.copy_from(solid.runtime_matrix());
      solid.assemble_runtime({});
      Vector<double> repeat=solid.runtime_rhs(); repeat-=rhs0;
      AssertThrow(repeat.l2_norm()==0,ExcMessage("repeated runtime assembly accumulated force"));

      Vector<double> direction(rhs0.size()); direction[free_dof]=1.0;
      Vector<double> jd(rhs0.size()); matrix0.vmult(jd,direction);
      const double epsilons[]={1e-4,1e-6,1e-8,1e-10,1e-12};
      double best=1e300,wrong_sign_best=1e300;
      for (const double epsilon:epsilons)
        {
          solid.runtime_solution()=0; solid.runtime_solution()[free_dof]=epsilon;
          auto shifted=messages(epsilon,force0,tangent_value);
          solid.set_generalized_interface_data(shifted.first,shifted.second);
          solid.assemble_runtime({});
          Vector<double> plus=solid.runtime_rhs(); plus*=-1;
          solid.runtime_solution()=0; solid.runtime_solution()[free_dof]=-epsilon;
          shifted=messages(-epsilon,force0,tangent_value);
          solid.set_generalized_interface_data(shifted.first,shifted.second);
          solid.assemble_runtime({});
          Vector<double> minus=solid.runtime_rhs(); minus*=-1;
          Vector<double> fd(plus); fd-=minus; fd/=(2*epsilon);
          Vector<double> error(fd); error-=jd;
          Vector<double> wrong_sign(fd); wrong_sign+=jd;
          best=std::min(best,error.l2_norm()/std::max(jd.l2_norm(),1e-300));
          wrong_sign_best=std::min(wrong_sign_best,
                                   wrong_sign.l2_norm()/std::max(jd.l2_norm(),1e-300));
        }
      std::cerr << "runtime_fd_best=" << best << std::endl;
      AssertThrow(best<1e-8,ExcMessage("runtime assemble_newton finite difference failed"));
      AssertThrow(wrong_sign_best>1.0,
                  ExcMessage("runtime finite difference did not detect inverted tangent sign"));
      solid.runtime_solution()=0;
      auto zero=messages(0,0,0);
      solid.set_generalized_interface_data(zero.first,zero.second);
      solid.assemble_runtime({});
      AssertThrow(solid.runtime_rhs()[free_dof]==0,
                  ExcMessage("exactly zero valid generalized force was not accepted"));
      solid.clear_provisional_generalized_interface_state();
      solid.set_generalized_interface_load(data.first);
      bool rejected_missing_tangent=false;
      try { solid.assemble_runtime({}); }
      catch (const std::exception &) { rejected_missing_tangent=true; }
      AssertThrow(rejected_missing_tangent,ExcMessage("runtime accepted missing tangent"));
      auto stale=data.first; --stale.stamp.outer_corrector;
      stale.checksum=generalized_interface::checksum(stale);
      bool rejected_stale=false;
      try { solid.set_generalized_interface_load(stale); }
      catch (const std::exception &) { rejected_stale=true; }
      AssertThrow(rejected_stale,ExcMessage("runtime accepted stale load"));
      auto wrong_hash=data.first; wrong_hash.stamp.hash_weights="wrong";
      wrong_hash.checksum=generalized_interface::checksum(wrong_hash);
      bool rejected_hash=false;
      try { solid.set_generalized_interface_load(wrong_hash); }
      catch (const std::exception &) { rejected_hash=true; }
      AssertThrow(rejected_hash,ExcMessage("runtime accepted incorrect hash"));
      generalized_interface::Expected next_expected{7,4,11,5,manifest.hash_graph,
                                                      "weights-v1","dofs-v1"};
      solid.begin_generalized_interface_corrector(next_expected);
      stamp.outer_corrector=4;
      auto next=messages(0,force0,tangent_value);
      solid.set_generalized_interface_data(next.first,next.second);
      solid.assemble_runtime({});
      const unsigned int ranks=Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
      AssertThrow(Utilities::MPI::min(rhs0.l2_norm(),MPI_COMM_WORLD)==
                  Utilities::MPI::max(rhs0.l2_norm(),MPI_COMM_WORLD),
                  ExcMessage("replicated runtime force differs"));
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
        std::cout << "{\"ranks\":" << ranks << ",\"free_dof\":" << free_dof
                  << ",\"unit_rhs\":" << rhs0[free_dof]
                  << ",\"repeat_error\":" << repeat.l2_norm()
                  << ",\"runtime_fd_best\":" << best
                  << ",\"wrong_sign_error\":" << wrong_sign_best << "}" << std::endl;
    }
  catch (const std::exception &e)
    { std::cerr << "testRuntimeNewton failed: " << e.what() << std::endl; return 1; }
}
