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
#include "RuntimeQ1InterfaceOperator.H"
#include "GeneralizedInterfaceLoad.H"
#include "ExecutableProtocol.H"

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
  void run_executable_protocol();
  void initialize_runtime() { read_mesh(); setup(); }
  void prepare_executable_protocol() { read_mesh(); prepare_setup(); }
  void finish_executable_protocol(const interface_sparsity::Manifest &manifest,
                                  const generalized_interface::Expected &expected)
  { finish_setup(&manifest,&expected); }
  std::string dof_manifest_payload(std::string &hash) const;
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
  const DoFHandler<dim> &runtime_dofs() const { return dofs; }
  types::global_dof_index runtime_dof_count() const { return dofs.n_dofs(); }
  runtime_q1::Operator<dim> build_canonical_interface_q1_operator(
    const std::vector<runtime_q1::Query<dim>> &queries,
    const double relative_tolerance=1e-10) const
  {
    runtime_q1::Builder<dim> builder(mapping,dofs,constraints,
      prm.get_integer("interface boundary"),"reference",1,relative_tolerance);
    return builder.build(queries);
  }
  std::string runtime_assembly_diagnostic();
  std::string structural_state_payload(const std::string &kind) const;
  void accept_protocol_time_step();
  void rollback_protocol_time_step();

private:
  void receive_generalized_pair_if_complete()
  {
    if (provisional_force && provisional_tangent)
      generalized_load->receive(*provisional_force,*provisional_tangent);
  }
  void read_mesh();
  void prepare_setup();
  void finish_setup(const interface_sparsity::Manifest *manifest_override=nullptr,
                    const generalized_interface::Expected *expected_override=nullptr);
  void setup();
  double assemble_newton(const std::vector<Sample> &samples);
  void solve_newton(const std::vector<Sample> &samples);
  void read_state();
  void write_state() const;
  void write_results(const std::vector<Sample> &samples,
                     const std::vector<Point<dim>> &queries) const;

  ParameterHandler &prm;
  Triangulation<dim> mesh;
  MappingQ1<dim> mapping;
  FESystem<dim> fe;
  DoFHandler<dim> dofs;
  AffineConstraints<double> constraints;
  SparsityPattern sparsity;
  SparseMatrix<double> matrix;
  Vector<double> solution, old_solution, old_velocity;
  Vector<double> older_solution, older_velocity, accepted_acceleration, rhs;
  unsigned int history_depth = 0;
  std::vector<unsigned int> interface_components;
  std::vector<std::array<types::global_dof_index,dim>> interface_node_dofs;
  std::unique_ptr<generalized_interface::GeneralizedLoad> generalized_load;
  std::unique_ptr<generalized_interface::Expected> runtime_expected;
  std::unique_ptr<generalized_interface::ForceMessage> provisional_force;
  std::unique_ptr<generalized_interface::TangentMessage> provisional_tangent;
  std::unique_ptr<DynamicSparsityPattern> pending_sparsity;
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
  prepare_setup();
  finish_setup();
}

void PdmsSolid::prepare_setup()
{
  dofs.distribute_dofs(fe);
  constraints.clear();
  for (const std::string key : {"clamped boundary 1", "clamped boundary 2", "clamped boundary 3"})
    VectorTools::interpolate_boundary_values(
      dofs, prm.get_integer(key), Functions::ZeroFunction<dim>(dim), constraints);
  constraints.close();
  pending_sparsity=std::make_unique<DynamicSparsityPattern>(dofs.n_dofs());
  DoFTools::make_sparsity_pattern(dofs, *pending_sparsity, constraints, false);
}

void PdmsSolid::finish_setup(const interface_sparsity::Manifest *manifest_override,
                             const generalized_interface::Expected *expected_override)
{
  AssertThrow(pending_sparsity,ExcMessage("setup finalized without prepared sparsity"));
  DynamicSparsityPattern &dsp=*pending_sparsity;
  const std::string transfer=prm.get("interface transfer");
  if (transfer=="dualConservative")
    {
      interface_sparsity::Manifest loaded;
      if (!manifest_override)
        { const std::string path=prm.get("interface manifest");
          AssertThrow(!path.empty(),ExcMessage(
            "dualConservative requires interface manifest before matrix reinit"));
          loaded=interface_sparsity::read_manifest(path); }
      const auto &manifest=manifest_override ? *manifest_override : loaded;
      const std::string expected_hash=expected_override ? expected_override->hash_graph
                                                        : prm.get("interface graph hash");
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
  pending_sparsity.reset();
  matrix.reinit(sparsity);
  if (transfer=="dualConservative")
    {
      generalized_interface::Expected configured{
        static_cast<unsigned int>(prm.get_integer("interface time index")),
        static_cast<unsigned int>(prm.get_integer("interface outer corrector")),
        static_cast<std::uint64_t>(prm.get_integer("interface operator version")),
        static_cast<std::uint64_t>(prm.get_integer("interface z version")),
        prm.get("interface graph hash"),prm.get("interface weights hash"),
        prm.get("interface dof manifest hash")};
      const auto &expected=expected_override ? *expected_override : configured;
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
  accepted_acceleration.reinit(dofs.n_dofs());
  rhs.reinit(dofs.n_dofs());
  read_state();
}

std::string PdmsSolid::dof_manifest_payload(std::string &hash) const
{
  std::map<types::global_dof_index,Point<dim>> support;
  DoFTools::map_dofs_to_support_points(MappingQ1<dim>(),dofs,support);
  std::vector<unsigned int> components(dofs.n_dofs(),numbers::invalid_unsigned_int);
  std::vector<types::global_dof_index> local(fe.n_dofs_per_cell());
  for (const auto &cell:dofs.active_cell_iterators())
    { cell->get_dof_indices(local); for (unsigned int i=0;i<local.size();++i)
        components[local[i]]=fe.system_to_component_index(i).first; }
  std::ostringstream body; body << std::setprecision(17)
    << "schemaVersion 1\ndimension " << dim << "\nscalarBytes " << sizeof(double)
    << "\nreplicatedMatrix 1\nentries " << dofs.n_dofs() << '\n';
  for (types::global_dof_index id=0;id<dofs.n_dofs();++id)
    { const auto &p=support.at(id); body << id << ' ' << components[id] << ' '
        << p[0] << ' ' << p[1] << ' ' << p[2] << " 0 "
        << (constraints.is_constrained(id) ? 1 : 0) << '\n'; }
  body << "end\n";
  std::ostringstream value; value << std::hex << executable_protocol::fnv(body.str());
  hash=value.str(); return "dofManifestHash "+hash+"\n"+body.str();
}

std::string PdmsSolid::runtime_assembly_diagnostic()
{
  Vector<double> interface_rhs(rhs.size());
  SparseMatrix<double> interface_matrix; interface_matrix.reinit(sparsity);
  generalized_load->add_to_newton(interface_rhs,interface_matrix);
  double interface_matrix_norm=0; std::size_t interface_entries=0;
  for (unsigned int i=0;i<interface_matrix.m();++i)
    for (auto entry=interface_matrix.begin(i);entry!=interface_matrix.end(i);++entry)
      if (entry->value()!=0)
        { interface_matrix_norm+=entry->value()*entry->value(); ++interface_entries; }
  assemble_newton({});
  const Vector<double> first_rhs(rhs);
  SparseMatrix<double> first_matrix; first_matrix.reinit(sparsity); first_matrix.copy_from(matrix);
  assemble_newton({});
  Vector<double> rhs_delta(rhs); rhs_delta-=first_rhs;
  double matrix_delta=0,matrix_norm=0;
  std::ostringstream canonical; canonical << std::hexfloat;
  for (unsigned int i=0;i<matrix.m();++i)
    for (auto entry=matrix.begin(i);entry!=matrix.end(i);++entry)
      { const double prior=first_matrix.el(i,entry->column());
        matrix_delta=std::max(matrix_delta,std::abs(entry->value()-prior));
        matrix_norm+=entry->value()*entry->value();
        if (entry->value()!=0) canonical << i << ' ' << entry->column() << ' '
                                         << entry->value() << '\n'; }
  std::ostringstream rhs_text; rhs_text << std::hexfloat;
  for (unsigned int i=0;i<rhs.size();++i) if (rhs[i]!=0) rhs_text << i << ' ' << rhs[i] << '\n';
  const std::uint64_t rhs_hash=executable_protocol::fnv(rhs_text.str());
  const std::uint64_t jacobian_hash=executable_protocol::fnv(canonical.str());
  const unsigned int ranks=Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  AssertThrow(Utilities::MPI::min(rhs.l2_norm(),MPI_COMM_WORLD)==
              Utilities::MPI::max(rhs.l2_norm(),MPI_COMM_WORLD),
              ExcMessage("replicated residual norm mismatch"));
  AssertThrow(Utilities::MPI::min(std::sqrt(matrix_norm),MPI_COMM_WORLD)==
              Utilities::MPI::max(std::sqrt(matrix_norm),MPI_COMM_WORLD),
              ExcMessage("replicated Jacobian norm mismatch"));
  AssertThrow(Utilities::MPI::min(rhs_hash,MPI_COMM_WORLD)==
              Utilities::MPI::max(rhs_hash,MPI_COMM_WORLD) &&
              Utilities::MPI::min(jacobian_hash,MPI_COMM_WORLD)==
              Utilities::MPI::max(jacobian_hash,MPI_COMM_WORLD),
              ExcMessage("replicated residual/Jacobian checksum mismatch"));
  std::ostringstream out; out << std::setprecision(17)
    << "ranks " << ranks << "\nrhsNorm " << rhs.l2_norm()
    << "\njacobianNorm " << std::sqrt(matrix_norm)
    << "\nrhsChecksum " << rhs_hash << "\njacobianChecksum " << jacobian_hash
    << "\ninterfaceRhsNorm " << interface_rhs.l2_norm()
    << "\ninterfaceJacobianNorm " << std::sqrt(interface_matrix_norm)
    << "\ninterfaceJacobianEntries " << interface_entries
    << "\nrhsRepeatError " << rhs_delta.l2_norm()
    << "\njacobianRepeatError " << matrix_delta << "\nassemblies 2\nend\n";
  return out.str();
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
      input >> std::ws;
      if (input.peek()!=std::char_traits<char>::eof())
        accepted_acceleration.block_read(input);
      else
        {
          accepted_acceleration=old_velocity;
          accepted_acceleration-=older_velocity;
          accepted_acceleration/=prm.get_double("delta t");
        }
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
  Vector<double> acceleration(velocity);
  acceleration.add(-1.0,old_velocity);
  if (bdf2)
    {
      acceleration *= 3.0;
      acceleration.add(-1.0,old_velocity);
      acceleration.add(1.0,older_velocity);
      acceleration /= 2.0*dt;
    }
  else
    acceleration /= dt;
  acceleration.block_write(output);
}

std::string PdmsSolid::structural_state_payload(const std::string &kind) const
{
  AssertThrow(kind=="accepted" || kind=="provisional",
              ExcMessage("invalid structural state kind"));
  const bool provisional=kind=="provisional";
  const double dt=prm.get_double("delta t");
  const bool bdf2=prm.get("time integration")=="bdf2" && history_depth>=1;
  Vector<double> velocity(dofs.n_dofs()),acceleration(dofs.n_dofs());
  const Vector<double> &displacement=provisional ? solution : old_solution;
  if (provisional)
    for (types::global_dof_index id=0;id<dofs.n_dofs();++id)
      {
        velocity[id]=dealii_pdms::time_integration::velocity
          (solution[id],old_solution[id],older_solution[id],dt,bdf2);
        acceleration[id]=dealii_pdms::time_integration::acceleration
          (velocity[id],old_velocity[id],older_velocity[id],dt,bdf2);
      }
  else
    { velocity=old_velocity; acceleration=accepted_acceleration; }

  std::map<types::global_dof_index,unsigned int> interface_dofs;
  for (const auto &node:interface_node_dofs)
    for (unsigned int component=0;component<dim;++component)
      interface_dofs.emplace(node[component],component);
  std::ostringstream out; out << std::setprecision(17)
    << "stateKind " << kind << "\nintegrationScheme " << (bdf2 ? "BDF2" : "BE")
    << "\ndisplacementUnits m\nvelocityUnits m/s\naccelerationUnits m/s2\nvalid 1\nentries "
    << interface_dofs.size() << '\n';
  for (const auto &entry:interface_dofs)
    out << entry.first << ' ' << entry.second << ' ' << displacement[entry.first]
        << ' ' << velocity[entry.first] << ' ' << acceleration[entry.first] << '\n';
  out << "end\n";
  return out.str();
}

void PdmsSolid::accept_protocol_time_step()
{
  const double dt=prm.get_double("delta t");
  const bool bdf2=prm.get("time integration")=="bdf2" && history_depth>=1;
  Vector<double> velocity(dofs.n_dofs()),acceleration(dofs.n_dofs());
  for (types::global_dof_index id=0;id<dofs.n_dofs();++id)
    {
      velocity[id]=dealii_pdms::time_integration::velocity
        (solution[id],old_solution[id],older_solution[id],dt,bdf2);
      acceleration[id]=dealii_pdms::time_integration::acceleration
        (velocity[id],old_velocity[id],older_velocity[id],dt,bdf2);
    }
  older_solution=old_solution; older_velocity=old_velocity;
  old_solution=solution; old_velocity=velocity; accepted_acceleration=acceleration;
  history_depth=std::min(history_depth+1,2u);
  clear_provisional_generalized_interface_state();
  write_state();
}

void PdmsSolid::rollback_protocol_time_step()
{
  solution=old_solution;
  clear_provisional_generalized_interface_state();
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

void PdmsSolid::run_executable_protocol()
{
  const unsigned int rank=Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  std::unique_ptr<executable_protocol::UnixChannel> channel;
  if (rank==0) channel=std::make_unique<executable_protocol::UnixChannel>(
    prm.get("protocol endpoint"),prm.get_integer("protocol timeout seconds"));
  const auto receive=[&]()
    {
      std::string wire,error;
      if (rank==0)
        try { wire=executable_protocol::encode(channel->receive()); }
        catch (const std::exception &failure) { error=failure.what(); }
      error=Utilities::MPI::broadcast(MPI_COMM_WORLD,error,0);
      AssertThrow(error.empty(),ExcMessage("external frame rejected: "+error));
      wire=Utilities::MPI::broadcast(MPI_COMM_WORLD,wire,0);
      return executable_protocol::decode(wire);
    };
  const auto send=[&](const executable_protocol::Type type,const executable_protocol::Frame &request,
                      const std::string &payload)
    { if (rank==0) { executable_protocol::Frame reply; reply.message_type=type;
        reply.sequence=request.sequence; reply.producer="dealiiPdmsSolid";
        reply.time_index=request.time_index; reply.outer_corrector=request.outer_corrector;
        reply.operator_version=request.operator_version; reply.payload=payload; channel->send(reply); } };
  executable_protocol::Frame current;
  try
    {
      current=receive();
      AssertThrow(current.message_type==executable_protocol::Type::Hello,
                  ExcMessage("AwaitingHello: expected Hello"));
      AssertThrow(current.payload==
        "dimension 3\nscalarBytes 8\nsignConvention R=Rs-fGamma;J=Js+JGamma\n"
        "forceUnits N\ntangentUnits N/m\ntimeIntegration BE,BDF2\nreplicatedMatrix 1\nend\n",
        ExcMessage("Hello capabilities mismatch"));
      send(executable_protocol::Type::Capabilities,current,current.payload);

      prepare_executable_protocol();
      std::string dof_hash; const std::string dof_payload=dof_manifest_payload(dof_hash);
      send(executable_protocol::Type::DofManifest,current,dof_payload);

      current=receive();
      AssertThrow(current.message_type==executable_protocol::Type::OperatorManifest,
                  ExcMessage("AwaitingOperatorManifest: expected OperatorManifest"));
      std::istringstream metadata(current.payload); std::string key,graph,weights,units,sign;
      std::uint64_t operator_version=0,z_version=0; std::size_t manifest_length=0;
      metadata >> key >> graph; AssertThrow(key=="hashGraph",ExcMessage("missing hashGraph"));
      metadata >> key >> weights; AssertThrow(key=="hashWeights",ExcMessage("missing hashWeights"));
      metadata >> key >> operator_version; AssertThrow(key=="operatorVersion",ExcMessage("missing operatorVersion"));
      metadata >> key >> z_version; AssertThrow(key=="zVersion",ExcMessage("missing zVersion"));
      metadata >> key >> units; AssertThrow(key=="units" && units=="H:1,W:m2",ExcMessage("operator units mismatch"));
      metadata >> key; std::getline(metadata,sign); AssertThrow(key=="sign" && sign==" fGamma=-HtWtf",ExcMessage("operator sign mismatch"));
      metadata >> key >> manifest_length; AssertThrow(key=="manifestBytes",ExcMessage("missing manifestBytes"));
      metadata.get(); const std::streampos position=metadata.tellg();
      AssertThrow(position>=0 && current.payload.size()==static_cast<std::size_t>(position)+manifest_length,
                  ExcMessage("operator manifest length/trailing bytes mismatch"));
      const std::string manifest_text=current.payload.substr(static_cast<std::size_t>(position));
      const auto manifest=interface_sparsity::read_manifest_text(manifest_text);
      AssertThrow(manifest.hash_graph==graph && current.operator_version==operator_version,
                  ExcMessage("operator manifest hashes/version mismatch"));
      generalized_interface::Expected expected{
        static_cast<unsigned int>(current.time_index),static_cast<unsigned int>(current.outer_corrector),
        operator_version,z_version,graph,weights,dof_hash};
      finish_executable_protocol(manifest,expected);
      send(executable_protocol::Type::Ready,current,
           "hashGraph "+graph+"\nhashWeights "+weights+"\ndofManifestHash "+dof_hash+"\nend\n");

      std::unique_ptr<executable_protocol::Frame> pending_force,pending_tangent;
      std::map<std::uint64_t,std::uint64_t> seen;
      std::uint64_t expected_time=current.time_index,expected_outer=current.outer_corrector;
      for (;;)
        {
          current=receive();
          const auto prior=seen.find(current.sequence);
          if (prior!=seen.end())
            { AssertThrow(prior->second==current.message_checksum,
                          ExcMessage("sequenceNumber replay with different content"));
              send(executable_protocol::Type::Ack,current,"idempotent 1\nend\n"); continue; }
          seen.emplace(current.sequence,current.message_checksum);
          if (current.message_type==executable_protocol::Type::ForceMessage)
            { (void)executable_protocol::parse_force(current);
              pending_force=std::make_unique<executable_protocol::Frame>(current);
              send(executable_protocol::Type::Ack,current,"pending ForceMessage\nend\n"); }
          else if (current.message_type==executable_protocol::Type::TangentMessage)
            { (void)executable_protocol::parse_tangent(current);
              pending_tangent=std::make_unique<executable_protocol::Frame>(current);
              send(executable_protocol::Type::Ack,current,"pending TangentMessage\nend\n"); }
          else if (current.message_type==executable_protocol::Type::ClearProvisionalState)
            { AssertThrow((current.time_index==expected_time && current.outer_corrector==expected_outer+1) ||
                          (current.time_index==expected_time+1 && current.outer_corrector==1),
                          ExcMessage("invalid corrector transition"));
              expected_time=current.time_index; expected_outer=current.outer_corrector;
              pending_force.reset(); pending_tangent.reset();
              clear_provisional_generalized_interface_state();
              send(executable_protocol::Type::Ack,current,"cleared 1\nend\n"); }
          else if (current.message_type==executable_protocol::Type::ActivateCorrectorState)
            {
              AssertThrow(pending_force && pending_tangent,ExcMessage("incomplete corrector state"));
              AssertThrow(current.time_index==expected_time && current.outer_corrector==expected_outer,
                          ExcMessage("stale or future corrector activation"));
              std::istringstream activation(current.payload); std::uint64_t force_checksum=0,tangent_checksum=0;
              std::string activation_graph,activation_weights;
              activation >> key >> force_checksum; AssertThrow(key=="forceMessageChecksum",ExcMessage("missing force checksum"));
              activation >> key >> tangent_checksum; AssertThrow(key=="tangentMessageChecksum",ExcMessage("missing tangent checksum"));
              activation >> key >> activation_graph; AssertThrow(key=="hashGraph",ExcMessage("missing activation graph"));
              activation >> key >> activation_weights; AssertThrow(key=="hashWeights",ExcMessage("missing activation weights"));
              activation >> key; AssertThrow(key=="end",ExcMessage("truncated activation")); activation >> std::ws;
              AssertThrow(activation.eof() && force_checksum==pending_force->message_checksum &&
                tangent_checksum==pending_tangent->message_checksum && activation_graph==graph &&
                activation_weights==weights,ExcMessage("ActivateCorrectorState mismatch"));
              AssertThrow(current.time_index==pending_force->time_index &&
                current.outer_corrector==pending_force->outer_corrector &&
                current.operator_version==pending_force->operator_version &&
                current.time_index==pending_tangent->time_index &&
                current.outer_corrector==pending_tangent->outer_corrector &&
                current.operator_version==pending_tangent->operator_version,
                ExcMessage("corrector stamp mismatch"));
              const generalized_interface::Expected active{
                static_cast<unsigned int>(current.time_index),static_cast<unsigned int>(current.outer_corrector),
                current.operator_version,z_version,graph,weights,dof_hash};
              begin_generalized_interface_corrector(active);
              set_generalized_interface_data(executable_protocol::parse_force(*pending_force),
                                             executable_protocol::parse_tangent(*pending_tangent));
              send(executable_protocol::Type::Ack,current,"activated 1\nend\n");
              send(executable_protocol::Type::AssemblyResult,current,runtime_assembly_diagnostic());
            }
          else if (current.message_type==executable_protocol::Type::RequestStructuralState)
            {
              std::istringstream request(current.payload); std::string state_kind;
              request >> key >> state_kind; AssertThrow(key=="stateKind",ExcMessage("missing stateKind"));
              request >> key; AssertThrow(key=="end",ExcMessage("truncated state request")); request >> std::ws;
              AssertThrow(request.eof(),ExcMessage("trailing state request bytes"));
              if (state_kind=="provisional")
                {
                  AssertThrow(generalized_load && generalized_load->is_received(),
                              ExcMessage("provisional state requested before corrector activation"));
                  solve_newton({});
                }
              send(executable_protocol::Type::StructuralStateMessage,current,
                   structural_state_payload(state_kind));
            }
          else if (current.message_type==executable_protocol::Type::AcceptTimeStep)
            {
              AssertThrow(generalized_load && generalized_load->is_received(),
                          ExcMessage("cannot accept without active corrector state"));
              accept_protocol_time_step(); pending_force.reset(); pending_tangent.reset();
              send(executable_protocol::Type::Ack,current,"accepted 1\nend\n");
            }
          else if (current.message_type==executable_protocol::Type::RejectTimeStep ||
                   current.message_type==executable_protocol::Type::RollbackToAcceptedState)
            {
              rollback_protocol_time_step(); pending_force.reset(); pending_tangent.reset();
              send(executable_protocol::Type::Ack,current,"rolledBack 1\nend\n");
            }
          else if (current.message_type==executable_protocol::Type::Shutdown)
            { pending_force.reset(); pending_tangent.reset();
              clear_provisional_generalized_interface_state();
              send(executable_protocol::Type::Ack,current,"shutdown 1\nend\n"); break; }
          else AssertThrow(false,ExcMessage("unexpected message in Ready state"));
        }
    }
  catch (const std::exception &error)
    {
      if (rank==0 && channel)
        try { send(executable_protocol::Type::Nack,current,
                   std::string("error ")+error.what()+"\nend\n"); } catch (...) {}
      throw;
    }
}

void PdmsSolid::run()
{
  if (prm.get("interface transfer")=="dualConservative")
    { run_executable_protocol(); return; }
  read_mesh();
  setup();
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
      prm.declare_entry("protocol endpoint", "", Patterns::Anything());
      prm.declare_entry("protocol timeout seconds", "10", Patterns::Integer(1,3600));
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
      if (prm.get("interface transfer")=="dualConservative")
        {
          AssertThrow(!prm.get("protocol endpoint").empty(),
                      ExcMessage("dualConservative requires protocol endpoint"));
          Utilities::MPI::MPI_InitFinalize mpi(argc,argv,1);
          PdmsSolid(prm).run();
        }
      else
        PdmsSolid(prm).run();
    }
  catch (const std::exception &e)
    {
      std::cerr << "deal.II PDMS solid failed: " << e.what() << std::endl;
      return 1;
    }
}
#endif
