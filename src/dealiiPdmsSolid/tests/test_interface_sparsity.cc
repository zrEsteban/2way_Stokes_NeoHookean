#include "InterfaceSparsityExtension.H"

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
                  << "}" << std::endl;
    }
  catch (const std::exception &e)
    { std::cerr << "testInterfaceSparsity failed: " << e.what() << std::endl; return 1; }
}
