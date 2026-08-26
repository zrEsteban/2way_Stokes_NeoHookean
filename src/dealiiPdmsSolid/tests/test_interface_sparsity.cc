#include "InterfaceSparsityExtension.H"
#include "GeneralizedInterfaceLoad.H"

#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <map>

using namespace dealii;
using namespace interface_sparsity;

int main(int argc,char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc,argv,1);
  try
    {
      AssertThrow(argc==3,ExcMessage("usage: testInterfaceSparsity solid.msh manifest"));
      Triangulation<3> mesh; GridIn<3> reader; reader.attach_triangulation(mesh);
      std::ifstream mesh_input(argv[1]); AssertThrow(mesh_input,ExcMessage("mesh absent"));
      reader.read_msh(mesh_input);
      FESystem<3> fe(FE_Q<3>(1),3); DoFHandler<3> dofs(mesh); dofs.distribute_dofs(fe);
      AffineConstraints<double> constraints;
      for (const unsigned int boundary:{1u,2u,3u})
        VectorTools::interpolate_boundary_values(
          dofs,boundary,Functions::ZeroFunction<3>(3),constraints);
      constraints.close();
      const auto manifest=read_manifest(argv[2]);
      AssertThrow(manifest.hash_graph==
        "6ee94c285c4e5d8619f9f08c3a666eb57dd46ddfe47b2b02be25c1248b4abc39",
        ExcMessage("hashGraph does not match canonical H graph"));
      const auto node_dofs=map_nodes_to_dofs<3>(manifest,dofs);
      DynamicSparsityPattern legacy(dofs.n_dofs());
      DoFTools::make_sparsity_pattern(dofs,legacy,constraints,false);
      const std::size_t legacy_nnz=legacy.n_nonzero_elements();
      DynamicSparsityPattern extended(dofs.n_dofs());
      DoFTools::make_sparsity_pattern(dofs,extended,constraints,false);
      const auto stats=augment_interface_sparsity<3>(extended,manifest,node_dofs,constraints);
      AssertThrow(stats.added>0,ExcMessage("legacy pattern unexpectedly contains dual graph"));
      AssertThrow(stats.final_missing==0,ExcMessage("extended graph is incomplete"));
      SparsityPattern pattern; pattern.copy_from(extended);
      SparseMatrix<double> matrix; matrix.reinit(pattern);

      auto expansion=[&](unsigned int node,unsigned int c)
        { return constraint_expansion(node_dofs.at(node)[c],constraints); };
      for (std::size_t r=0;r<manifest.rows.size();++r)
        {
          const double z=20000.0+80000.0*(double(r)+1.0)/manifest.rows.size();
          for (unsigned int c=0;c<3;++c)
            for (const auto &a:manifest.rows[r].entries)
              for (const auto &b:manifest.rows[r].entries)
                for (const auto &i:expansion(a.node,c)) for (const auto &j:expansion(b.node,c))
                  { AssertThrow(pattern.exists(i.first,j.first),ExcMessage("required connection absent"));
                    matrix.add(i.first,j.first,
                               i.second*j.second*a.weight*b.weight*manifest.rows[r].area*z); }
        }
      matrix.compress(VectorOperation::add);
      std::vector<unsigned int> components(dofs.n_dofs(),numbers::invalid_unsigned_int);
      for (const auto &node:node_dofs) for (unsigned int c=0;c<3;++c) components[node[c]]=c;
      generalized_interface::Expected expected{7,3,11,5,manifest.hash_graph,
                                                "weights-canonical-v1","dofs-canonical-v1"};
      generalized_interface::Stamp stamp;
      stamp.time_index=expected.time_index; stamp.outer_corrector=expected.outer_corrector;
      stamp.operator_version=expected.operator_version; stamp.z_version=expected.z_version;
      stamp.hash_graph=expected.hash_graph; stamp.hash_weights=expected.hash_weights;
      stamp.dof_manifest_hash=expected.dof_manifest_hash; stamp.valid=true;
      std::map<types::global_dof_index,double> force_map;
      std::map<std::pair<types::global_dof_index,types::global_dof_index>,double> tangent_map;
      std::map<std::pair<types::global_dof_index,types::global_dof_index>,double> uniform_tangent_map;
      std::mt19937_64 protocol_generator(20260825);
      std::normal_distribution<double> protocol_normal;
      double fluid_work=0;
      Vector<double> virtual_displacement(dofs.n_dofs());
      for (unsigned int i=0;i<virtual_displacement.size();++i) virtual_displacement[i]=protocol_normal(protocol_generator);
      constraints.set_zero(virtual_displacement);
      for (std::size_t r=0;r<manifest.rows.size();++r)
        {
          const double z=20000.0+80000.0*(double(r)+1.0)/manifest.rows.size();
          for (unsigned int c=0;c<3;++c)
            {
              const double traction=(c==0 ? 1200.0 : protocol_normal(protocol_generator)*300.0);
              double hd=0;
              for (const auto &a:manifest.rows[r].entries)
                { const auto dof=node_dofs[a.node][c];
                  force_map[dof]+=-a.weight*manifest.rows[r].area*traction;
                  hd+=a.weight*virtual_displacement[dof]; }
              fluid_work+=hd*manifest.rows[r].area*traction;
              for (const auto &a:manifest.rows[r].entries)
                for (const auto &b:manifest.rows[r].entries)
                  {
                    tangent_map[{node_dofs[a.node][c],node_dofs[b.node][c]}]+=
                      a.weight*b.weight*manifest.rows[r].area*z;
                    uniform_tangent_map[{node_dofs[a.node][c],node_dofs[b.node][c]}]+=
                      a.weight*b.weight*manifest.rows[r].area*67350.0;
                  }
            }
        }
      generalized_interface::ForceMessage force; force.stamp=stamp;
      for (const auto &v:force_map) force.values.push_back({v.first,components[v.first],v.second});
      force.checksum=generalized_interface::checksum(force);
      generalized_interface::TangentMessage tangent; tangent.stamp=stamp;
      for (const auto &v:tangent_map) tangent.values.push_back({v.first.first,v.first.second,v.second});
      tangent.checksum=generalized_interface::checksum(tangent);
      SparseMatrix<double> protocol_matrix; protocol_matrix.reinit(pattern);
      Vector<double> protocol_rhs(dofs.n_dofs());
      generalized_interface::GeneralizedLoad api(constraints,pattern,components,expected);
      api.receive(force,tangent); api.add_to_newton(protocol_rhs,protocol_matrix);
      protocol_matrix.compress(VectorOperation::add);
      Vector<double> expected_force(dofs.n_dofs());
      for (const auto &v:force_map)
        for (const auto &i:constraint_expansion(v.first,constraints)) expected_force[i.first]+=i.second*v.second;
      Vector<double> force_difference(protocol_rhs); force_difference-=expected_force;
      const double force_error=force_difference.l2_norm()/std::max(expected_force.l2_norm(),1e-300);
      AssertThrow(force_error<1e-12,ExcMessage("generalized force API mismatch"));
      Vector<double> protocol_kx(dofs.n_dofs()); protocol_matrix.vmult(protocol_kx,virtual_displacement);
      Vector<double> protocol_reference(dofs.n_dofs()); matrix.vmult(protocol_reference,virtual_displacement);
      Vector<double> matrix_difference(protocol_kx); matrix_difference-=protocol_reference;
      const double protocol_kx_error=matrix_difference.l2_norm()/std::max(protocol_reference.l2_norm(),1e-300);
      AssertThrow(protocol_kx_error<1e-12,ExcMessage("generalized tangent API mismatch"));
      const double solid_work=virtual_displacement*protocol_rhs;
      const double work_error=std::abs(fluid_work+solid_work)/
        std::max({std::abs(fluid_work),std::abs(solid_work),1e-30});
      AssertThrow(work_error<1e-12,ExcMessage("constrained virtual work mismatch"));

      // R_Gamma(d)=-f0+J_Gamma*d; the API convention is J_Gamma=-df/dd.
      const double epsilons[]={1e-2,1e-4,1e-6,1e-8,1e-10,1e-12};
      double best_fd=1e300, inverted_fd=0;
      Vector<double> base_r(protocol_rhs); base_r*=-1.0;
      for (const double epsilon:epsilons)
        { Vector<double> perturbed(base_r); perturbed.add(epsilon,protocol_kx);
          Vector<double> quotient(perturbed); quotient-=base_r; quotient/=epsilon;
          quotient-=protocol_kx;
          best_fd=std::min(best_fd,quotient.l2_norm()/std::max(protocol_kx.l2_norm(),1e-300));
          Vector<double> wrong(protocol_kx); wrong*=-1; wrong-=protocol_kx;
          inverted_fd=wrong.l2_norm()/std::max(protocol_kx.l2_norm(),1e-300); }
      AssertThrow(best_fd<1e-8 && inverted_fd>1.9,
                  ExcMessage("finite difference did not validate/detect tangent sign"));

      // Complete synthetic structural residual: linear structural stiffness,
      // cubic internal force, and the real assembled interface contribution.
      Vector<double> state(dofs.n_dofs()),direction(virtual_displacement);
      for (unsigned int i=0;i<state.size();++i) state[i]=1e-3*protocol_normal(protocol_generator);
      constraints.set_zero(state); direction/=std::max(direction.l2_norm(),1e-300);
      auto full_residual=[&](const Vector<double> &d)
        { Vector<double> result(dofs.n_dofs()),interface_part(dofs.n_dofs());
          protocol_matrix.vmult(interface_part,d); result=interface_part; result-=protocol_rhs;
          for (unsigned int i=0;i<d.size();++i) result[i]+=2.5e4*d[i]+1e8*d[i]*d[i]*d[i];
          return result; };
      Vector<double> full_jd(dofs.n_dofs()); protocol_matrix.vmult(full_jd,direction);
      for (unsigned int i=0;i<state.size();++i)
        full_jd[i]+=(2.5e4+3e8*state[i]*state[i])*direction[i];
      const Vector<double> full_base=full_residual(state);
      std::vector<double> full_fd_errors;
      for (const double epsilon:epsilons)
        { Vector<double> shifted(state); shifted.add(epsilon,direction);
          Vector<double> quotient=full_residual(shifted); quotient-=full_base; quotient/=epsilon;
          quotient-=full_jd;
          full_fd_errors.push_back(quotient.l2_norm()/std::max(full_jd.l2_norm(),1e-300)); }
      const double full_best=*std::min_element(full_fd_errors.begin(),full_fd_errors.end());
      AssertThrow(full_best<1e-8 && full_fd_errors.front()>full_best,
                  ExcMessage("complete residual finite-difference convergence missing"));

      auto rejected=[&](auto bad_force,auto bad_tangent)
        { try { generalized_interface::GeneralizedLoad candidate(constraints,pattern,components,expected);
                candidate.receive(bad_force,bad_tangent); }
          catch (const std::exception &) { return true; } return false; };
      auto stale=force; --stale.stamp.outer_corrector; stale.checksum=generalized_interface::checksum(stale);
      AssertThrow(rejected(stale,tangent),ExcMessage("stale force accepted"));
      auto bad_hash=force; bad_hash.stamp.hash_graph="bad"; bad_hash.checksum=generalized_interface::checksum(bad_hash);
      AssertThrow(rejected(bad_hash,tangent),ExcMessage("bad hash accepted"));
      auto duplicate=force; duplicate.values.push_back(duplicate.values.front()); duplicate.checksum=generalized_interface::checksum(duplicate);
      AssertThrow(rejected(duplicate,tangent),ExcMessage("duplicate force accepted"));
      auto bad_units=tangent; bad_units.units="Pa"; bad_units.checksum=generalized_interface::checksum(bad_units);
      AssertThrow(rejected(force,bad_units),ExcMessage("bad tangent units accepted"));
      auto zero=force; for (auto &v:zero.values) v.value=0; zero.checksum=generalized_interface::checksum(zero);
      generalized_interface::GeneralizedLoad zero_api(constraints,pattern,components,expected);
      zero_api.receive(zero,tangent); AssertThrow(zero_api.is_received(),ExcMessage("valid zero force rejected"));
      generalized_interface::TangentMessage uniform_tangent; uniform_tangent.stamp=stamp;
      for (const auto &v:uniform_tangent_map)
        uniform_tangent.values.push_back({v.first.first,v.first.second,v.second});
      uniform_tangent.checksum=generalized_interface::checksum(uniform_tangent);
      auto free_dof=force.values.begin();
      while (free_dof!=force.values.end() && constraints.is_constrained(free_dof->dof)) ++free_dof;
      AssertThrow(free_dof!=force.values.end(),ExcMessage("no free interface DoF"));
      auto unit=zero; unit.values.clear(); unit.values.push_back({free_dof->dof,free_dof->component,1.0});
      unit.checksum=generalized_interface::checksum(unit);
      generalized_interface::GeneralizedLoad unit_api(constraints,pattern,components,expected);
      unit_api.receive(unit,uniform_tangent);
      SparseMatrix<double> unit_matrix; unit_matrix.reinit(pattern); Vector<double> unit_rhs(dofs.n_dofs());
      unit_api.add_to_newton(unit_rhs,unit_matrix);
      AssertThrow(unit_rhs[free_dof->dof]==1.0 && unit_rhs.l1_norm()==1.0,
                  ExcMessage("unit generalized nodal force changed or integrated"));
      Vector<double> x(dofs.n_dofs()),kx(dofs.n_dofs()),reference(dofs.n_dofs());
      std::mt19937_64 generator(20260825); std::normal_distribution<double> normal;
      for (unsigned int i=0;i<x.size();++i) x[i]=normal(generator);
      constraints.set_zero(x); matrix.vmult(kx,x);
      for (std::size_t r=0;r<manifest.rows.size();++r)
        {
          const double z=20000.0+80000.0*(double(r)+1.0)/manifest.rows.size();
          for (unsigned int c=0;c<3;++c)
            { double hx=0;
              for (const auto &a:manifest.rows[r].entries)
                for (const auto &i:expansion(a.node,c)) hx+=a.weight*i.second*x[i.first];
              const double q=manifest.rows[r].area*z*hx;
              for (const auto &a:manifest.rows[r].entries)
                for (const auto &i:expansion(a.node,c)) reference[i.first]+=a.weight*i.second*q; }
        }
      Vector<double> difference(kx); difference-=reference;
      const double action_error=difference.l2_norm()/std::max(reference.l2_norm(),1e-300);
      AssertThrow(action_error<1e-12,ExcMessage("real SparseMatrix Kx mismatch"));
      const double energy=matrix.matrix_scalar_product(x,x);
      AssertThrow(energy>=-1e-13*std::max(std::abs(energy),1.0),
                  ExcMessage("K is not positive semidefinite"));
      double symmetry2=0,norm2=0;
      for (unsigned int i=0;i<matrix.m();++i) for (auto e=matrix.begin(i);e!=matrix.end(i);++e)
        { const double d=e->value()-matrix.el(e->column(),i); symmetry2+=d*d; norm2+=e->value()*e->value(); }
      const double symmetry=std::sqrt(symmetry2)/std::max(std::sqrt(norm2),1e-300);
      AssertThrow(symmetry<1e-12,ExcMessage("K is not symmetric"));
      const unsigned int ranks=Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
      const auto min_added=Utilities::MPI::min(stats.added,MPI_COMM_WORLD);
      const auto max_added=Utilities::MPI::max(stats.added,MPI_COMM_WORLD);
      AssertThrow(min_added==max_added,ExcMessage("partition-dependent graph"));
      const auto min_force_checksum=Utilities::MPI::min(force.checksum,MPI_COMM_WORLD);
      const auto max_force_checksum=Utilities::MPI::max(force.checksum,MPI_COMM_WORLD);
      const auto min_tangent_checksum=Utilities::MPI::min(tangent.checksum,MPI_COMM_WORLD);
      const auto max_tangent_checksum=Utilities::MPI::max(tangent.checksum,MPI_COMM_WORLD);
      const double min_rhs_norm=Utilities::MPI::min(protocol_rhs.l2_norm(),MPI_COMM_WORLD);
      const double max_rhs_norm=Utilities::MPI::max(protocol_rhs.l2_norm(),MPI_COMM_WORLD);
      const double min_matrix_norm=Utilities::MPI::min(protocol_matrix.frobenius_norm(),MPI_COMM_WORLD);
      const double max_matrix_norm=Utilities::MPI::max(protocol_matrix.frobenius_norm(),MPI_COMM_WORLD);
      AssertThrow(min_force_checksum==max_force_checksum && min_tangent_checksum==max_tangent_checksum &&
                  min_rhs_norm==max_rhs_norm && min_matrix_norm==max_matrix_norm,
                  ExcMessage("replicated interface assembly differs between ranks"));
      bool rejected_bad_id=false;
      try
        { auto invalid=node_dofs;
          auto target=invalid.begin();
          while (target!=invalid.end() && constraints.is_constrained((*target)[0])) ++target;
          AssertThrow(target!=invalid.end(),ExcMessage("no free DoF for negative test"));
          (*target)[0]=dofs.n_dofs();
          DynamicSparsityPattern negative(dofs.n_dofs());
          DoFTools::make_sparsity_pattern(dofs,negative,constraints,false);
          augment_interface_sparsity<3>(negative,manifest,invalid,constraints); }
      catch (const std::exception &) { rejected_bad_id=true; }
      AssertThrow(rejected_bad_id,ExcMessage("out-of-range DoF was accepted"));
      bool rejected_finalized=false;
      try
        { DynamicSparsityPattern finalized;
          augment_interface_sparsity<3>(finalized,manifest,node_dofs,constraints); }
      catch (const std::exception &) { rejected_finalized=true; }
      AssertThrow(rejected_finalized,ExcMessage("late/finalized augmentation was accepted"));
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
        std::cout << std::setprecision(17)
                  << "{\"ranks\":" << ranks << ",\"backend\":\"dealii::SparseMatrix<double> replicated\""
                  << ",\"legacy_nnz\":" << legacy_nnz
                  << ",\"required\":" << stats.required << ",\"added\":" << stats.added
                  << ",\"final_nnz\":" << pattern.n_nonzero_elements()
                  << ",\"missing\":" << stats.final_missing
                  << ",\"effective_hash_fnv64\":\"" << std::hex << stats.effective_hash << std::dec << "\""
                  << ",\"action_error\":" << action_error << ",\"symmetry_error\":" << symmetry
                  << ",\"energy\":" << energy
                  << ",\"negative_bad_id\":true,\"negative_finalized\":true"
                  << ",\"force_error\":" << force_error
                  << ",\"protocol_kx_error\":" << protocol_kx_error
                  << ",\"work_error\":" << work_error << ",\"best_fd_error\":" << best_fd
                  << ",\"inverted_sign_fd_error\":" << inverted_fd
                  << ",\"full_fd_best\":" << full_best
                  << ",\"full_fd_large_epsilon\":" << full_fd_errors.front()
                  << ",\"full_fd_small_epsilon\":" << full_fd_errors.back()
                  << "}" << std::endl;
    }
  catch (const std::exception &e)
    { std::cerr << "testInterfaceSparsity failed: " << e.what() << std::endl; return 1; }
}
