#pragma once

#include <deal.II/base/point.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <array>
#include <cstdint>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace interface_sparsity
{
using dealii::types::global_dof_index;

struct Entry { unsigned int node; double weight; };
struct Row { unsigned int face; double area; std::vector<Entry> entries; };
struct Manifest
{
  std::string hash_graph;
  std::vector<dealii::Point<3>> nodes;
  std::vector<Row> rows;
};

inline Manifest read_manifest(std::istream &in)
{
  std::string word; unsigned int version=0;
  in >> word >> version;
  AssertThrow(word=="G2B1A_INTERFACE_GRAPH" && version==1,
              dealii::ExcMessage("Unsupported interface manifest schema"));
  Manifest m; std::size_t count=0;
  in >> word >> m.hash_graph;
  AssertThrow(word=="hashGraph" && m.hash_graph.size()==64,
              dealii::ExcMessage("Invalid hashGraph"));
  in >> word >> count;
  AssertThrow(word=="nodes", dealii::ExcMessage("Missing manifest nodes"));
  m.nodes.resize(count);
  for (std::size_t i=0;i<count;++i)
    { unsigned int id; in >> word >> id >> m.nodes[i][0] >> m.nodes[i][1] >> m.nodes[i][2];
      AssertThrow(word=="node" && id==i, dealii::ExcMessage("Invalid node ID/order")); }
  in >> word >> count;
  AssertThrow(word=="rows", dealii::ExcMessage("Missing manifest rows"));
  m.rows.resize(count);
  for (std::size_t r=0;r<count;++r)
    { std::size_t n; in >> word >> m.rows[r].face >> m.rows[r].area >> n;
      AssertThrow(word=="row" && m.rows[r].face==r && n>0,
                  dealii::ExcMessage("Invalid row ID/order/support"));
      m.rows[r].entries.resize(n);
      for (auto &e:m.rows[r].entries)
        { in >> e.node >> e.weight;
          AssertThrow(e.node<m.nodes.size(),dealii::ExcMessage("Manifest node ID out of range")); }
    }
  in >> word;
  AssertThrow(in && word=="end",dealii::ExcMessage("Truncated interface manifest"));
  return m;
}

inline Manifest read_manifest(const std::string &path)
{
  std::ifstream in(path);
  AssertThrow(in, dealii::ExcMessage("Cannot open interface manifest: " + path));
  return read_manifest(in);
}

inline Manifest read_manifest_text(const std::string &text)
{
  std::istringstream in(text);
  const auto manifest=read_manifest(in);
  in >> std::ws;
  AssertThrow(in.eof(),dealii::ExcMessage("Trailing bytes in interface manifest"));
  return manifest;
}

template <int dim>
std::vector<std::array<global_dof_index,dim>>
map_nodes_to_dofs(const Manifest &m, const dealii::DoFHandler<dim> &dofs)
{
  std::map<global_dof_index,dealii::Point<dim>> support;
  dealii::DoFTools::map_dofs_to_support_points(dealii::MappingQ1<dim>(),dofs,support);
  std::map<global_dof_index,unsigned int> component;
  std::vector<global_dof_index> local(dofs.get_fe().n_dofs_per_cell());
  for (const auto &cell:dofs.active_cell_iterators())
    { cell->get_dof_indices(local);
      for (unsigned int i=0;i<local.size();++i)
        component[local[i]]=dofs.get_fe().system_to_component_index(i).first; }
  std::vector<std::array<global_dof_index,dim>> result(m.nodes.size());
  for (auto &a:result) a.fill(dealii::numbers::invalid_dof_index);
  for (const auto &item:support)
    { std::size_t best=m.nodes.size(); double distance=1e300;
      for (std::size_t n=0;n<m.nodes.size();++n)
        if (item.second.distance_square(m.nodes[n])<distance)
          { best=n; distance=item.second.distance_square(m.nodes[n]); }
      AssertThrow(best<m.nodes.size() && distance<1e-24,
                  dealii::ExcMessage("DoF support point absent from manifest"));
      const unsigned int c=component.at(item.first);
      AssertThrow(c<dim && result[best][c]==dealii::numbers::invalid_dof_index,
                  dealii::ExcMessage("Invalid or duplicate component mapping"));
      result[best][c]=item.first;
    }
  for (const auto &a:result) for (const auto id:a)
    AssertThrow(id!=dealii::numbers::invalid_dof_index,
                dealii::ExcMessage("Incomplete manifest-to-DoF mapping"));
  return result;
}

inline std::vector<std::pair<global_dof_index,double>>
constraint_expansion(const global_dof_index id,
                     const dealii::AffineConstraints<double> &constraints)
{
  if (!constraints.is_constrained(id)) return {{id,1.0}};
  const auto *entries=constraints.get_constraint_entries(id);
  if (entries==nullptr || entries->empty()) return {};
  return *entries;
}

struct Statistics
{
  std::size_t required=0, already_present=0, added=0, final_missing=0;
  std::uint64_t effective_hash=1469598103934665603ULL;
};

template <int dim>
Statistics augment_interface_sparsity(
  dealii::DynamicSparsityPattern &dsp, const Manifest &manifest,
  const std::vector<std::array<global_dof_index,dim>> &node_dofs,
  const dealii::AffineConstraints<double> &constraints)
{
  AssertThrow(dsp.n_rows()>0,dealii::ExcMessage("Cannot augment a finalized/empty pattern"));
  std::set<std::pair<global_dof_index,global_dof_index>> pairs;
  for (const auto &row:manifest.rows)
    for (unsigned int c=0;c<dim;++c)
      for (const auto &a:row.entries) for (const auto &b:row.entries)
        for (const auto &i:constraint_expansion(node_dofs.at(a.node)[c],constraints))
          for (const auto &j:constraint_expansion(node_dofs.at(b.node)[c],constraints))
            { AssertThrow(i.first<dsp.n_rows() && j.first<dsp.n_cols(),
                          dealii::ExcMessage("Expanded DoF ID out of range"));
              pairs.emplace(i.first,j.first); pairs.emplace(j.first,i.first); }
  Statistics s; s.required=pairs.size();
  for (const auto &p:pairs)
    { const bool existed=dsp.exists(p.first,p.second); s.already_present+=existed;
      if (!existed) { dsp.add(p.first,p.second); ++s.added; }
      s.effective_hash^=p.first; s.effective_hash*=1099511628211ULL;
      s.effective_hash^=p.second; s.effective_hash*=1099511628211ULL; }
  for (const auto &p:pairs) s.final_missing+=!dsp.exists(p.first,p.second);
  return s;
}
}
