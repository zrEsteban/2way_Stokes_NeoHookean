#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "NeoHookeanMaterial.H"
#include "TimeIntegration.H"
#include "InterfaceSparsityExtension.H"
#include "GeneralizedInterfaceLoad.H"

using namespace dealii;

namespace
{
constexpr unsigned int dim = 3;

struct Sample
{
  Point<dim> x;
  Tensor<1, dim> traction;
  Tensor<1, dim> velocity;
};

std::vector<double> split_numbers(const std::string &line)
{
  std::vector<double> values;
  std::stringstream stream(line);
  std::string token;
  while (std::getline(stream, token, ','))
    {
      try { values.push_back(std::stod(token)); }
      catch (const std::exception &) { return {}; }
    }
  return values;
}

std::vector<Sample> read_samples(const std::string &name)
{
  std::ifstream input(name);
  AssertThrow(input, ExcMessage("Cannot open Robin input: " + name));
  std::vector<Sample> samples;
  std::string line;
  while (std::getline(input, line))
    {
      const auto a = split_numbers(line);
      if (a.size() != 9) continue;
      Sample s;
      for (unsigned int d=0; d<dim; ++d)
        {
          s.x[d] = a[d];
          s.traction[d] = a[3+d];
          s.velocity[d] = a[6+d];
        }
      samples.push_back(s);
    }
  AssertThrow(!samples.empty(), ExcMessage("Robin input has no numeric samples"));
  return samples;
}

std::vector<Point<dim>> read_query_points(const std::string &name)
{
  std::vector<Point<dim>> points;
  if (name.empty()) return points;
  std::ifstream input(name);
  AssertThrow(input, ExcMessage("Cannot open query input: " + name));
  std::string line;
  while (std::getline(input, line))
    {
      const auto a = split_numbers(line);
      if (a.size() != dim) continue;
      Point<dim> point;
      for (unsigned int d=0; d<dim; ++d) point[d]=a[d];
      points.push_back(point);
    }
  AssertThrow(!points.empty(), ExcMessage("Query input has no numeric points"));
  return points;
}

const Sample &nearest(const std::vector<Sample> &samples, const Point<dim> &p)
{
  std::size_t best = 0;
  double distance = std::numeric_limits<double>::max();
  for (std::size_t i=0; i<samples.size(); ++i)
    if (p.distance_square(samples[i].x) < distance)
      {
        best = i;
        distance = p.distance_square(samples[i].x);
      }
  return samples[best];
}
}

class PdmsSolid
{
public:
  explicit PdmsSolid(ParameterHandler &parameters)
    : prm(parameters), fe(FE_Q<dim>(1), dim), dofs(mesh) {}
  void run();
  void initialize_runtime() { read_mesh(); setup(); }
  void set_generalized_interface_data(
    const generalized_interface::ForceMessage &force,
    const generalized_interface::TangentMessage &tangent)
  {
    set_generalized_interface_load(force);
    set_generalized_interface_tangent(tangent);
  }
  void set_generalized_interface_load(const generalized_interface::ForceMessage &force)
  {
    AssertThrow(generalized_load && runtime_expected,
                ExcMessage("load received before generalized manifest/runtime"));
    generalized_interface::validate_stamp(force.stamp,*runtime_expected);
    AssertThrow(force.units=="N" && force.checksum==generalized_interface::checksum(force),
                ExcMessage("invalid generalized force units/checksum"));
    provisional_force=std::make_unique<generalized_interface::ForceMessage>(force);
    receive_generalized_pair_if_complete();
  }
  void set_generalized_interface_tangent(const generalized_interface::TangentMessage &tangent)
  {
    AssertThrow(generalized_load && runtime_expected,
                ExcMessage("tangent received before generalized manifest/runtime"));
    generalized_interface::validate_stamp(tangent.stamp,*runtime_expected);
    AssertThrow(tangent.units=="N/m" && tangent.checksum==generalized_interface::checksum(tangent),
                ExcMessage("invalid generalized tangent units/checksum"));
    provisional_tangent=std::make_unique<generalized_interface::TangentMessage>(tangent);
    receive_generalized_pair_if_complete();
  }
  void clear_provisional_generalized_interface_state()
  {
    provisional_force.reset(); provisional_tangent.reset();
    if (generalized_load) generalized_load->clear_provisional();
  }
  void begin_generalized_interface_corrector(const generalized_interface::Expected &expected)
  {
    AssertThrow(runtime_expected && expected.hash_graph==runtime_expected->hash_graph,
                ExcMessage("hashGraph change after SparseMatrix::reinit is forbidden"));
    clear_provisional_generalized_interface_state();
    runtime_expected=std::make_unique<generalized_interface::Expected>(expected);
    generalized_load=std::make_unique<generalized_interface::GeneralizedLoad>(
      constraints,sparsity,interface_components,*runtime_expected);
  }
  double assemble_runtime(const std::vector<Sample> &samples)
  { return assemble_newton(samples); }
  const Vector<double> &runtime_rhs() const { return rhs; }
  const SparseMatrix<double> &runtime_matrix() const { return matrix; }
  Vector<double> &runtime_solution() { return solution; }
  const std::vector<std::array<types::global_dof_index,dim>> &runtime_node_dofs() const
  { return interface_node_dofs; }
  const AffineConstraints<double> &runtime_constraints() const { return constraints; }
  const SparsityPattern &runtime_sparsity() const { return sparsity; }

private:
  void receive_generalized_pair_if_complete()
  {
    if (provisional_force && provisional_tangent)
      generalized_load->receive(*provisional_force,*provisional_tangent);
  }
  void read_mesh();
  void setup();
  double assemble_newton(const std::vector<Sample> &samples);
  void solve_newton(const std::vector<Sample> &samples);
  void read_state();
  void write_state() const;
  void write_results(const std::vector<Sample> &samples,
                     const std::vector<Point<dim>> &queries) const;

  ParameterHandler &prm;
  Triangulation<dim> mesh;
  FESystem<dim> fe;
  DoFHandler<dim> dofs;
  AffineConstraints<double> constraints;
  SparsityPattern sparsity;
  SparseMatrix<double> matrix;
  Vector<double> solution, old_solution, old_velocity;
  Vector<double> older_solution, older_velocity, rhs;
  unsigned int history_depth = 0;
  std::vector<unsigned int> interface_components;
  std::vector<std::array<types::global_dof_index,dim>> interface_node_dofs;
  std::unique_ptr<generalized_interface::GeneralizedLoad> generalized_load;
  std::unique_ptr<generalized_interface::Expected> runtime_expected;
  std::unique_ptr<generalized_interface::ForceMessage> provisional_force;
  std::unique_ptr<generalized_interface::TangentMessage> provisional_tangent;
};

void PdmsSolid::read_mesh()
{
  GridIn<dim> reader;
  reader.attach_triangulation(mesh);
  std::ifstream input(prm.get("mesh"));
  AssertThrow(input, ExcMessage("Cannot open mesh: " + prm.get("mesh")));
  reader.read_msh(input);
}

void PdmsSolid::setup()
{
  dofs.distribute_dofs(fe);
  constraints.clear();
  for (const std::string key : {"clamped boundary 1", "clamped boundary 2", "clamped boundary 3"})
    VectorTools::interpolate_boundary_values(
      dofs, prm.get_integer(key), Functions::ZeroFunction<dim>(dim), constraints);
  constraints.close();
  DynamicSparsityPattern dsp(dofs.n_dofs());
  DoFTools::make_sparsity_pattern(dofs, dsp, constraints, false);
  const std::string transfer=prm.get("interface transfer");
  if (transfer=="dualConservative")
    {
      const std::string manifest_path=prm.get("interface manifest");
      AssertThrow(!manifest_path.empty(),ExcMessage(
        "dualConservative requires interface manifest before matrix reinit"));
      const auto manifest=interface_sparsity::read_manifest(manifest_path);
      const std::string expected_hash=prm.get("interface graph hash");
      AssertThrow(!expected_hash.empty() && manifest.hash_graph==expected_hash,
                  ExcMessage("interface manifest hashGraph mismatch"));
      const auto node_dofs=interface_sparsity::map_nodes_to_dofs<dim>(manifest,dofs);
      interface_node_dofs=node_dofs;
      const auto stats=interface_sparsity::augment_interface_sparsity<dim>(
        dsp,manifest,node_dofs,constraints);
      AssertThrow(stats.final_missing==0,ExcMessage("Incomplete dual interface sparsity"));
      std::cout << "Dual interface graph prepared before SparseMatrix::reinit: hashGraph="
                << manifest.hash_graph << ", required=" << stats.required
                << ", added=" << stats.added << std::endl;
      interface_components.assign(dofs.n_dofs(),numbers::invalid_unsigned_int);
      for (const auto &node:node_dofs)
        for (unsigned int component=0;component<dim;++component)
          interface_components[node[component]]=component;
    }
  sparsity.copy_from(dsp);
  matrix.reinit(sparsity);
  if (transfer=="dualConservative")
    {
      generalized_interface::Expected expected{
        static_cast<unsigned int>(prm.get_integer("interface time index")),
        static_cast<unsigned int>(prm.get_integer("interface outer corrector")),
        static_cast<std::uint64_t>(prm.get_integer("interface operator version")),
        static_cast<std::uint64_t>(prm.get_integer("interface z version")),
        prm.get("interface graph hash"),prm.get("interface weights hash"),
        prm.get("interface dof manifest hash")};
      AssertThrow(!expected.hash_weights.empty() && !expected.dof_manifest_hash.empty(),
                  ExcMessage("dualConservative requires all runtime hashes"));
      runtime_expected=std::make_unique<generalized_interface::Expected>(expected);
      generalized_load=std::make_unique<generalized_interface::GeneralizedLoad>(
        constraints,sparsity,interface_components,*runtime_expected);
    }
  solution.reinit(dofs.n_dofs());
  old_solution.reinit(dofs.n_dofs());
  old_velocity.reinit(dofs.n_dofs());
  older_solution.reinit(dofs.n_dofs());
  older_velocity.reinit(dofs.n_dofs());
  rhs.reinit(dofs.n_dofs());
  read_state();
}

void PdmsSolid::read_state()
{
  const std::string state = prm.get("state input");
  if (state.empty()) return;
  std::ifstream input(state);
  if (!input) return;
  old_solution.block_read(input);
  old_velocity.block_read(input);
  AssertThrow(old_solution.size() == dofs.n_dofs(),
              ExcMessage("State vector does not match the current FEM mesh"));
  solution = old_solution;
  // Legacy states contain only (d^n,v^n).  Extended states append
  // (d^{n-1},v^{n-1},history_depth), so existing cases restart with one BE
  // step before BDF2 is enabled.
  input >> std::ws;
  if (input.peek() != std::char_traits<char>::eof())
    {
      older_solution.block_read(input);
      older_velocity.block_read(input);
      input >> history_depth;
      history_depth = std::min(history_depth,2u);
    }
  else
    {
      older_solution = old_solution;
      older_velocity = old_velocity;
      history_depth = 0;
    }
}

void PdmsSolid::write_state() const
{
  const std::string state = prm.get("state output");
  if (state.empty()) return;
  std::ofstream output(state);
  solution.block_write(output);
  const double dt = prm.get_double("delta t");
  const bool bdf2 = prm.get("time integration") == "bdf2" && history_depth >= 1;
  Vector<double> velocity(solution);
  if (bdf2)
    {
      velocity *= 3.0;
      velocity.add(-4.0,old_solution);
      velocity.add(1.0,older_solution);
      velocity /= 2.0*dt;
    }
  else
    {
      velocity -= old_solution;
      velocity /= dt;
    }
  velocity.block_write(output);
  old_solution.block_write(output);
  old_velocity.block_write(output);
  output << std::min(history_depth+1,2u) << '\n';
}

double PdmsSolid::assemble_newton(const std::vector<Sample> &samples)
{
  matrix = 0;
  rhs = 0;
  const double rho = prm.get_double("rho");
  const double mu = prm.get_double("mu");
  const double bulk = prm.get_double("bulk modulus");
  const double dt = prm.get_double("delta t");
  const double impedance = prm.get_double("solid impedance");
  const bool bdf2 = prm.get("time integration") == "bdf2" && history_depth >= 1;
  const double velocity_coefficient =
    dealii_pdms::time_integration::velocity_coefficient(dt,bdf2);
  const double acceleration_coefficient =
    dealii_pdms::time_integration::acceleration_coefficient(dt,bdf2);
  const QGauss<dim> quadrature(2);
  const QGauss<dim-1> face_quadrature(2);
  FEValues<dim> values(fe, quadrature,
                       update_values | update_gradients | update_JxW_values);
  FEFaceValues<dim> face_values(fe, face_quadrature,
                                update_values | update_quadrature_points |
                                update_JxW_values);
  const FEValuesExtractors::Vector displacement(0);
  FullMatrix<double> cell_matrix(fe.n_dofs_per_cell(), fe.n_dofs_per_cell());
  Vector<double> cell_rhs(fe.n_dofs_per_cell());
  std::vector<types::global_dof_index> local(fe.n_dofs_per_cell());
  std::vector<Tensor<2,dim>> gradients(quadrature.size());
  std::vector<Tensor<1,dim>> displacements(quadrature.size());
  std::vector<Tensor<1,dim>> old_displacements(quadrature.size());
  std::vector<Tensor<1,dim>> old_velocities(quadrature.size());
  std::vector<Tensor<1,dim>> older_displacements(quadrature.size());
  std::vector<Tensor<1,dim>> older_velocities(quadrature.size());
  double minimum_j = std::numeric_limits<double>::max();

  for (const auto &cell : dofs.active_cell_iterators())
    {
      cell_matrix = 0;
      cell_rhs = 0;
      cell->get_dof_indices(local);
      values.reinit(cell);
      values[displacement].get_function_gradients(solution, gradients);
      values[displacement].get_function_values(solution, displacements);
      values[displacement].get_function_values(old_solution, old_displacements);
      values[displacement].get_function_values(old_velocity, old_velocities);
      if (bdf2)
        {
          values[displacement].get_function_values
            (older_solution, older_displacements);
          values[displacement].get_function_values
            (older_velocity, older_velocities);
        }
      for (unsigned int q=0; q<quadrature.size(); ++q)
        {
          Tensor<2,dim> F = unit_symmetric_tensor<dim>();
          F += gradients[q];
          const double J = determinant(F);
          minimum_j = std::min(minimum_j,J);
          AssertThrow(J > prm.get_double("minimum jacobian"),
                      ExcMessage("neo-Hookean deformation rejected: det(F)="+
                                 std::to_string(J)));
          const Tensor<2,dim> P = dealii_pdms::first_piola<dim>(F,mu,bulk);
          const Tensor<4,dim> A = dealii_pdms::consistent_tangent<dim>(F,mu,bulk);
          Tensor<1,dim> velocity;
          Tensor<1,dim> acceleration;
          velocity = dealii_pdms::time_integration::velocity
            (displacements[q],old_displacements[q],older_displacements[q],dt,bdf2);
          acceleration = dealii_pdms::time_integration::acceleration
            (velocity,old_velocities[q],older_velocities[q],dt,bdf2);
          for (unsigned int i=0; i<fe.n_dofs_per_cell(); ++i)
            {
              const auto grad_i = values[displacement].gradient(i,q);
              const auto phi_i = values[displacement].value(i,q);
              cell_rhs(i) -= (scalar_product(grad_i,P)+rho*(phi_i*acceleration))*values.JxW(q);
              for (unsigned int j=0; j<fe.n_dofs_per_cell(); ++j)
                {
                  const auto grad_j = values[displacement].gradient(j,q);
                  double material = 0;
                  for (unsigned int a=0; a<dim; ++a)
                    for (unsigned int b=0; b<dim; ++b)
                      for (unsigned int c=0; c<dim; ++c)
                        for (unsigned int d=0; d<dim; ++d)
                          material += grad_i[a][b]*A[a][b][c][d]*grad_j[c][d];
                  cell_matrix(i,j) +=
                    (material+rho*acceleration_coefficient
                     *(phi_i*values[displacement].value(j,q)))
                    *values.JxW(q);
                }
            }
        }

      if (prm.get("interface transfer")=="legacyNearestNeighbour")
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->face(f)->at_boundary() &&
            cell->face(f)->boundary_id() == prm.get_integer("interface boundary"))
          {
            face_values.reinit(cell, f);
            for (unsigned int q=0; q<face_quadrature.size(); ++q)
              {
                const Sample &s = nearest(samples, face_values.quadrature_point(q));
                const Tensor<1,dim> robin = s.traction + impedance*s.velocity;
                for (unsigned int i=0; i<fe.n_dofs_per_cell(); ++i)
                  {
                    const auto phi_i = face_values[displacement].value(i,q);
                    Tensor<1,dim> u, u_old, u_older;
                    for (unsigned int j=0; j<fe.n_dofs_per_cell(); ++j)
                      {
                        const auto phi_j=face_values[displacement].value(j,q);
                        u += solution[local[j]]*phi_j;
                        u_old += old_solution[local[j]]*phi_j;
                        if (bdf2)
                          u_older += older_solution[local[j]]*phi_j;
                      }
                    const Tensor<1,dim> solid_velocity =
                      dealii_pdms::time_integration::velocity
                        (u,u_old,u_older,dt,bdf2);
                    cell_rhs(i) +=
                      (phi_i*(robin-impedance*solid_velocity))
                      *face_values.JxW(q);
                    for (unsigned int j=0; j<fe.n_dofs_per_cell(); ++j)
                      cell_matrix(i,j) += impedance*velocity_coefficient *
                        (phi_i*face_values[displacement].value(j,q))*face_values.JxW(q);
                  }
              }
          }

      constraints.distribute_local_to_global(cell_matrix, cell_rhs, local, matrix, rhs);
    }
  if (prm.get("interface transfer")=="dualConservative")
    {
      AssertThrow(generalized_load && generalized_load->is_received(),
                  ExcMessage("dualConservative requires valid force and tangent before assemble_newton"));
      generalized_load->add_to_newton(rhs,matrix);
    }
  return minimum_j;
}

void PdmsSolid::solve_newton(const std::vector<Sample> &samples)
{
  const unsigned int maximum = prm.get_integer("maximum newton iterations");
  const double absolute_tolerance = prm.get_double("newton tolerance");
  const double relative_tolerance = prm.get_double("newton relative tolerance");
  double initial_residual = -1;
  for (unsigned int iteration=0; iteration<maximum; ++iteration)
    {
      const double minimum_j = assemble_newton(samples);
      const double residual = rhs.l2_norm();
      if (initial_residual < 0) initial_residual = residual;
      std::cout << "Newton " << iteration << ": |R|=" << residual
                << ", minJ=" << minimum_j << std::endl;
      if (residual <= std::max(absolute_tolerance,
                               relative_tolerance*initial_residual)) return;
      Vector<double> increment(solution.size());
      if (prm.get("linear solver") == "direct")
        {
          SparseDirectUMFPACK direct;
          direct.initialize(matrix);
          direct.vmult(increment,rhs);
        }
      else
        {
          SolverControl control(prm.get_integer("linear maximum iterations"),
                                prm.get_double("linear relative tolerance")*
                                std::max(rhs.l2_norm(),1e-30));
          SolverCG<Vector<double>> cg(control);
          PreconditionSSOR<SparseMatrix<double>> preconditioner;
          preconditioner.initialize(matrix,prm.get_double("ssor relaxation"));
          cg.solve(matrix,increment,rhs,preconditioner);
          std::cout << "  linear CG iterations=" << control.last_step()
                    << ", residual=" << control.last_value() << std::endl;
        }
      constraints.distribute(increment);
      const Vector<double> base(solution);
      double damping = prm.get_double("newton damping");
      bool accepted = false;
      for (unsigned int line_search=0; line_search<10; ++line_search)
        {
          solution = base;
          solution.add(damping,increment);
          try
            {
              assemble_newton(samples);
              const double trial_residual=rhs.l2_norm();
              if (trial_residual < residual || damping <= 1.0/512.0)
                {
                  std::cout << "  line search: alpha=" << damping
                            << ", |R_trial|=" << trial_residual << std::endl;
                  accepted=true;
                  break;
                }
            }
          catch (const std::exception &)
            {}
          damping *= 0.5;
        }
      AssertThrow(accepted, ExcMessage("Newton line search failed to find a valid state"));
    }
  AssertThrow(false, ExcMessage("neo-Hookean Newton iteration did not converge"));
}

void PdmsSolid::write_results(const std::vector<Sample> &samples,
                              const std::vector<Point<dim>> &queries) const
{
  const double dt = prm.get_double("delta t");
  const double impedance = prm.get_double("solid impedance");
  const bool bdf2 = prm.get("time integration") == "bdf2" && history_depth >= 1;
  std::ofstream output(prm.get("output"));
  output << "x,y,z,ux,uy,uz,vx,vy,vz,tx,ty,tz,ax,ay,az\n"
         << std::setprecision(17);
  const auto write_point = [&](const Point<dim> &point, const Sample &s)
    {
      Vector<double> value(dim), old_value(dim), old_velocity_value(dim);
      Vector<double> older_value(dim), older_velocity_value(dim);
      VectorTools::point_value(dofs,solution,point,value);
      VectorTools::point_value(dofs,old_solution,point,old_value);
      VectorTools::point_value(dofs,old_velocity,point,old_velocity_value);
      if (bdf2)
        {
          VectorTools::point_value(dofs,older_solution,point,older_value);
          VectorTools::point_value(dofs,older_velocity,point,older_velocity_value);
        }
      Tensor<1,dim> u, velocity, acceleration;
      for (unsigned int d=0; d<dim; ++d)
        {
          u[d]=value[d];
          velocity[d]=dealii_pdms::time_integration::velocity
            (value[d],old_value[d],older_value[d],dt,bdf2);
          acceleration[d]=dealii_pdms::time_integration::acceleration
            (velocity[d],old_velocity_value[d],older_velocity_value[d],dt,bdf2);
        }
      const Tensor<1,dim> traction =
        s.traction + impedance*(s.velocity-velocity);
      output << point[0] << ',' << point[1] << ',' << point[2];
      for (unsigned int d=0; d<dim; ++d) output << ',' << u[d];
      for (unsigned int d=0; d<dim; ++d) output << ',' << velocity[d];
      for (unsigned int d=0; d<dim; ++d) output << ',' << traction[d];
      for (unsigned int d=0; d<dim; ++d) output << ',' << acceleration[d];
      output << '\n';
    };
  if (queries.empty())
    for (const Sample &s : samples) write_point(s.x,s);
  else
    for (const Point<dim> &point : queries) write_point(point,nearest(samples,point));

  DataOut<dim> data;
  data.attach_dof_handler(dofs);
  std::vector<std::string> names(dim, "displacement");
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);
  data.add_data_vector(solution, names, DataOut<dim>::type_dof_data, interpretation);
  data.build_patches();
  std::ofstream vtk(prm.get("vtk output"));
  data.write_vtu(vtk);
}

void PdmsSolid::run()
{
  read_mesh();
  setup();
  AssertThrow(prm.get("interface transfer")=="legacyNearestNeighbour",
              ExcMessage("dualConservative runtime API is available, but the standalone "
                         "executable has no live transport; inject data before assemble_newton"));
  const auto samples = read_samples(prm.get("input"));
  const auto queries = read_query_points(prm.get("query input"));
  solve_newton(samples);
  write_results(samples,queries);
  write_state();
  std::cout << "deal.II solid: cells=" << mesh.n_active_cells()
            << ", dofs=" << dofs.n_dofs()
            << ", time integration="
            << (prm.get("time integration") == "bdf2" && history_depth >= 1
                  ? "BDF2" : "backwardEuler startup")
            << std::endl;
}

#ifndef DEALII_PDMS_NO_MAIN
int main(int argc, char **argv)
{
  try
    {
      ParameterHandler prm;
      prm.declare_entry("mesh", "solid.msh", Patterns::Anything());
      prm.declare_entry("input", "robin-in.csv", Patterns::Anything());
      prm.declare_entry("output", "robin-out.csv", Patterns::Anything());
      prm.declare_entry("query input", "", Patterns::Anything());
      prm.declare_entry("vtk output", "solid.vtu", Patterns::Anything());
      prm.declare_entry("rho", "970", Patterns::Double(0));
      prm.declare_entry("mu", "778200", Patterns::Double(0));
      prm.declare_entry("bulk modulus", "3631600", Patterns::Double(0));
      prm.declare_entry("delta t", "1e-7", Patterns::Double(0));
      prm.declare_entry("time integration", "backwardEuler",
                        Patterns::Selection("backwardEuler|bdf2"));
      prm.declare_entry("solid impedance", "67350", Patterns::Double(0));
      prm.declare_entry("interface transfer", "legacyNearestNeighbour",
                        Patterns::Selection("legacyNearestNeighbour|dualConservative"));
      prm.declare_entry("interface manifest", "", Patterns::Anything());
      prm.declare_entry("interface graph hash", "", Patterns::Anything());
      prm.declare_entry("interface weights hash", "", Patterns::Anything());
      prm.declare_entry("interface dof manifest hash", "", Patterns::Anything());
      prm.declare_entry("interface time index", "0", Patterns::Integer(0));
      prm.declare_entry("interface outer corrector", "0", Patterns::Integer(0));
      prm.declare_entry("interface operator version", "0", Patterns::Integer(0));
      prm.declare_entry("interface z version", "0", Patterns::Integer(0));
      prm.declare_entry("interface boundary", "4", Patterns::Integer(0));
      prm.declare_entry("clamped boundary 1", "1", Patterns::Integer(0));
      prm.declare_entry("clamped boundary 2", "2", Patterns::Integer(0));
      prm.declare_entry("clamped boundary 3", "3", Patterns::Integer(0));
      prm.declare_entry("state input", "", Patterns::Anything());
      prm.declare_entry("state output", "", Patterns::Anything());
      prm.declare_entry("minimum jacobian", "1e-8", Patterns::Double(0));
      prm.declare_entry("maximum newton iterations", "20", Patterns::Integer(1));
      prm.declare_entry("newton tolerance", "1e-18", Patterns::Double(0));
      prm.declare_entry("newton relative tolerance", "1e-3", Patterns::Double(0));
      prm.declare_entry("newton damping", "1", Patterns::Double(0,1));
      prm.declare_entry("linear solver", "iterative",
                        Patterns::Selection("iterative|direct"));
      prm.declare_entry("linear maximum iterations", "5000", Patterns::Integer(1));
      prm.declare_entry("linear relative tolerance", "1e-10", Patterns::Double(0));
      prm.declare_entry("ssor relaxation", "1.2", Patterns::Double(0,2));
      AssertThrow(argc == 2, ExcMessage("usage: dealiiPdmsSolid parameters.prm"));
      prm.parse_input(argv[1]);
      PdmsSolid(prm).run();
    }
  catch (const std::exception &e)
    {
      std::cerr << "deal.II PDMS solid failed: " << e.what() << std::endl;
      return 1;
    }
}
#endif
