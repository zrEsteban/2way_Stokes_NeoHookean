#include "../RuntimeQ1InterfaceOperator.H"
#include <deal.II/base/mpi.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <iostream>

using namespace dealii;

template <class Callable> static void expect_failure(Callable callable)
{
  bool failed=false; try { callable(); } catch (const std::exception &) { failed=true; }
  AssertThrow(failed,ExcMessage("negative RuntimeQ1 edge test unexpectedly passed"));
}

int main(int argc,char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc,argv,1);
  try
    {
      Triangulation<3> tria; GridGenerator::hyper_cube(tria,0,1); tria.refine_global(1);
      for (const auto &cell:tria.active_cell_iterators())
        for (unsigned int f=0;f<GeometryInfo<3>::faces_per_cell;++f)
          if (cell->face(f)->at_boundary()) cell->face(f)->set_boundary_id(4);
      FESystem<3> fe(FE_Q<3>(1),3); DoFHandler<3> dofs(tria); dofs.distribute_dofs(fe);
      AffineConstraints<double> constraints; constraints.close(); MappingQ1<3> mapping;
      runtime_q1::Builder<3> builder(mapping,dofs,constraints,4,"reference",1,1e-10);
      std::vector<runtime_q1::Query<3>> queries={
        {1,Point<3>(0,.31,.42)}, {2,Point<3>(0,0,.37)}, {3,Point<3>(0,0,0)},
        {4,Point<3>(4e-15,0,.63)}};
      const auto first=builder.build(queries),second=builder.build(queries);
      AssertThrow(first.graph_hash==second.graph_hash && first.weights_hash==second.weights_hash,
                  ExcMessage("edge/vertex canonicalization is not deterministic"));
      AssertThrow(first.face_interior_points==1 && first.edge_points==2 && first.vertex_points==1,
                  ExcMessage("edge/vertex codimension classification failed"));
      AssertThrow(first.multiple_candidates==3,ExcMessage("shared edge/vertex candidates not audited"));
      for (const auto &row:first.rows)
        {
          double sum=0; for (const auto &entry:row.entries) sum+=entry.weight;
          AssertThrow(std::abs(sum-1)<=1e-14,ExcMessage("edge/vertex row sum"));
          const std::size_t expected=row.query_id==3 ? 1 : (row.query_id==1 ? 4 : 2);
          AssertThrow(row.entries.size()==expected,ExcMessage("edge/vertex Q1 support"));
        }
      expect_failure([&]{ builder.build({{1,Point<3>(0,0,0)},{1,Point<3>(0,0,0)}}); });
      expect_failure([&]{ builder.build({{5,Point<3>(.5,.5,.5)}}); });
      expect_failure([&]{ runtime_q1::Builder<3> wrong(mapping,dofs,constraints,99); wrong.build(queries); });
      expect_failure([&]{ runtime_q1::Builder<3> loose(mapping,dofs,constraints,4,"reference",1,1e-7); });
      expect_failure([&]{ runtime_q1::Builder<3> frame(mapping,dofs,constraints,4,"current"); });
      std::cout << "RUNTIME_Q1_EDGE PASS face=1 edge=2 vertex=1 negative=5" << std::endl;
    }
  catch (const std::exception &error)
    { std::cerr << "RUNTIME_Q1_EDGE FAIL: " << error.what() << std::endl; return 1; }
}
