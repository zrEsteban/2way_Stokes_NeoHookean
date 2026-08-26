#pragma once

#include "InterfaceSparsityExtension.H"

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <set>
#include <sstream>

namespace generalized_interface
{
using dealii::types::global_dof_index;

struct Stamp
{
  unsigned int schema_version=1;
  std::string transfer_mode="dualConservative";
  unsigned int time_index=0, outer_corrector=0;
  std::uint64_t operator_version=0, z_version=0;
  std::string hash_graph, hash_weights, dof_manifest_hash;
  bool valid=false;
};
struct DofValue { global_dof_index dof; unsigned int component; double value; };
struct Triplet { global_dof_index row,column; double value; };
struct ForceMessage { Stamp stamp; std::string units="N"; std::vector<DofValue> values; std::uint64_t checksum=0; };
struct TangentMessage { Stamp stamp; std::string units="N/m"; std::vector<Triplet> values; std::uint64_t checksum=0; };
struct Expected
{
  unsigned int time_index,outer_corrector;
  std::uint64_t operator_version,z_version;
  std::string hash_graph,hash_weights,dof_manifest_hash;
};

inline void fnv_bytes(std::uint64_t &hash,const void *data,const std::size_t size)
{
  const auto *bytes=static_cast<const unsigned char *>(data);
  for (std::size_t i=0;i<size;++i) { hash^=bytes[i]; hash*=1099511628211ULL; }
}
inline void fnv_string(std::uint64_t &hash,const std::string &value)
{ fnv_bytes(hash,value.data(),value.size()); const unsigned char zero=0; fnv_bytes(hash,&zero,1); }
inline std::string canonical_stamp(const Stamp &s)
{
  std::ostringstream out;
  out << s.schema_version << '|' << s.transfer_mode << '|' << s.time_index << '|'
      << s.outer_corrector << '|' << s.operator_version << '|' << s.z_version << '|'
      << s.hash_graph << '|' << s.hash_weights << '|' << s.dof_manifest_hash << '|'
      << (s.valid ? 1 : 0) << '\n';
  return out.str();
}
inline std::uint64_t checksum(const ForceMessage &m)
{
  std::ostringstream out; out << canonical_stamp(m.stamp) << m.units << '\n' << std::hexfloat;
  for (const auto &v:m.values) out << v.dof << '|' << v.component << '|' << v.value << '\n';
  std::uint64_t h=1469598103934665603ULL; const auto text=out.str(); fnv_string(h,text); return h;
}
inline std::uint64_t checksum(const TangentMessage &m)
{
  std::ostringstream out; out << canonical_stamp(m.stamp) << m.units << '\n' << std::hexfloat;
  for (const auto &v:m.values) out << v.row << '|' << v.column << '|' << v.value << '\n';
  std::uint64_t h=1469598103934665603ULL; const auto text=out.str(); fnv_string(h,text); return h;
}

inline void validate_stamp(const Stamp &s,const Expected &e)
{
  AssertThrow(s.valid,dealii::ExcMessage("generalized interface message is not valid"));
  AssertThrow(s.schema_version==1,dealii::ExcMessage("schemaVersion mismatch"));
  AssertThrow(s.transfer_mode=="dualConservative",dealii::ExcMessage("transferMode mismatch"));
  AssertThrow(s.time_index==e.time_index && s.outer_corrector==e.outer_corrector,
              dealii::ExcMessage("stale or future corrector message"));
  AssertThrow(s.operator_version==e.operator_version,dealii::ExcMessage("operatorVersion mismatch"));
  AssertThrow(s.z_version==e.z_version,dealii::ExcMessage("Z version mismatch"));
  AssertThrow(s.hash_graph==e.hash_graph && s.hash_weights==e.hash_weights &&
              s.dof_manifest_hash==e.dof_manifest_hash,dealii::ExcMessage("interface hash mismatch"));
}

class GeneralizedLoad
{
public:
  GeneralizedLoad(const dealii::AffineConstraints<double> &constraints,
                  const dealii::SparsityPattern &pattern,
                  const std::vector<unsigned int> &components,
                  const Expected &expected)
    : constraints(constraints),pattern(pattern),components(components),expected(expected) {}

  void receive(const ForceMessage &force,const TangentMessage &tangent)
  {
    validate_stamp(force.stamp,expected); validate_stamp(tangent.stamp,expected);
    AssertThrow(force.stamp.operator_version==tangent.stamp.operator_version &&
                force.stamp.time_index==tangent.stamp.time_index &&
                force.stamp.outer_corrector==tangent.stamp.outer_corrector,
                dealii::ExcMessage("force/tangent stamp mismatch"));
    AssertThrow(force.units=="N",dealii::ExcMessage("force units must be N"));
    AssertThrow(tangent.units=="N/m",dealii::ExcMessage("tangent units must be N/m"));
    AssertThrow(force.checksum==checksum(force),dealii::ExcMessage("force checksum mismatch"));
    AssertThrow(tangent.checksum==checksum(tangent),dealii::ExcMessage("tangent checksum mismatch"));
    force_values.clear(); tangent_values.clear();
    std::set<global_dof_index> force_ids;
    for (const auto &v:force.values)
      { validate_dof(v.dof,v.component); AssertThrow(force_ids.insert(v.dof).second,
          dealii::ExcMessage("duplicate generalized force DoF"));
        for (const auto &i:interface_sparsity::constraint_expansion(v.dof,constraints))
          force_values.emplace_back(i.first,i.second*v.value); }
    std::map<std::pair<global_dof_index,global_dof_index>,double> tangent_ids;
    for (const auto &v:tangent.values)
      { validate_dof(v.row,components.at(v.row)); validate_dof(v.column,components.at(v.column));
        AssertThrow(tangent_ids.emplace(std::make_pair(v.row,v.column),v.value).second,
                    dealii::ExcMessage("duplicate generalized tangent entry"));
        for (const auto &i:interface_sparsity::constraint_expansion(v.row,constraints))
          for (const auto &j:interface_sparsity::constraint_expansion(v.column,constraints))
            { AssertThrow(pattern.exists(i.first,j.first),
                          dealii::ExcMessage("generalized tangent outside sparsity pattern"));
              tangent_values.push_back({i.first,j.first,i.second*j.second*v.value}); }
      }
    for (const auto &entry:tangent_ids)
      {
        const auto reverse=tangent_ids.find({entry.first.second,entry.first.first});
        AssertThrow(reverse!=tangent_ids.end() &&
                    std::abs(reverse->second-entry.second)<=
                      1e-13*std::max({std::abs(reverse->second),std::abs(entry.second),1.0}),
                    dealii::ExcMessage("generalized tangent is not symmetric"));
      }
    received=true;
  }

  void add_to_newton(dealii::Vector<double> &rhs,dealii::SparseMatrix<double> &jacobian) const
  {
    AssertThrow(received,dealii::ExcMessage("generalized interface data not received"));
    for (const auto &v:force_values) rhs[v.first]+=v.second; // rhs=-R; R=Rs-f_Gamma
    for (const auto &v:tangent_values) jacobian.add(v.row,v.column,v.value); // J=Js+J_Gamma
  }
  bool is_received() const { return received; }
  void clear_provisional()
  { received=false; force_values.clear(); tangent_values.clear(); }

private:
  void validate_dof(const global_dof_index id,const unsigned int component) const
  {
    AssertThrow(id<components.size(),dealii::ExcMessage("global DoF ID out of range"));
    AssertThrow(component<3 && components[id]==component,
                dealii::ExcMessage("generalized interface component mismatch"));
  }
  const dealii::AffineConstraints<double> &constraints;
  const dealii::SparsityPattern &pattern;
  const std::vector<unsigned int> &components;
  Expected expected;
  bool received=false;
  std::vector<std::pair<global_dof_index,double>> force_values;
  std::vector<Triplet> tangent_values;
};
}
