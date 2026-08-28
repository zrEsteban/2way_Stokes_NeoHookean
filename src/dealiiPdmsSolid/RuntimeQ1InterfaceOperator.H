#pragma once

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/bounding_box.h>
#include <deal.II/base/point.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/numerics/rtree.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace runtime_q1
{
using global_dof_index=dealii::types::global_dof_index;

template <int dim> struct Query
{
  std::uint64_t id=0;
  dealii::Point<dim> point;
};

struct Entry { global_dof_index dof=0; double weight=0; };

template <int dim> struct Row
{
  std::uint64_t query_id=0;
  unsigned int component=0;
  std::string cell_id;
  unsigned int face_index=0;
  dealii::Point<dim> reference;
  double residual=0,row_sum=0;
  std::vector<Entry> entries;
};

template <int dim> struct Operator
{
  std::string coordinate_frame="reference",coordinate_units="m";
  std::uint64_t geometry_version=1;
  std::string query_manifest_hash,support_hash,graph_hash,weights_hash,geometry_hash;
  std::vector<Row<dim>> rows;
  double maximum_residual=0;
  std::size_t multiple_candidates=0;
  std::size_t snapped_coordinates=0;
  double maximum_reference_snap=0;
  std::size_t face_interior_points=0,edge_points=0,vertex_points=0;
  std::size_t equivalent_candidate_count=0;
  double maximum_pruned_weight_mass=0;
};

inline std::uint64_t fnv(const std::string &value)
{
  std::uint64_t hash=1469598103934665603ULL;
  for (const unsigned char byte:value) { hash^=byte; hash*=1099511628211ULL; }
  return hash;
}

inline std::string hexadecimal(const std::uint64_t value)
{ std::ostringstream out; out << std::hex << value; return out.str(); }

template <int dim> class Builder
{
public:
  Builder(const dealii::MappingQ1<dim> &mapping,
          const dealii::DoFHandler<dim> &dofs,
          const dealii::AffineConstraints<double> &constraints,
          const dealii::types::boundary_id boundary,
          const std::string &frame="reference",
          const std::uint64_t geometry_version=1,
          const double relative_tolerance=1e-10)
    : mapping(mapping),dofs(dofs),constraints(constraints),boundary(boundary),
      frame(frame),geometry_version(geometry_version),relative_tolerance(relative_tolerance)
  {
    if (frame!="reference") throw std::runtime_error("RuntimeQ1 requires coordinateFrame=reference");
    if (!(relative_tolerance>0) || relative_tolerance>1e-8)
      throw std::runtime_error("RuntimeQ1 relative tolerance must be in (0,1e-8]");
    const auto &fe=dofs.get_fe();
    if (fe.n_components()!=dim || fe.degree!=1 || fe.get_name().find("FE_Q<3>(1)")==std::string::npos)
      throw std::runtime_error("RuntimeQ1 requires vector FESystem(FE_Q(1),dim)");
  }

  Operator<dim> build(std::vector<Query<dim>> queries) const
  {
    if (queries.empty()) throw std::runtime_error("RuntimeQ1 query list is empty");
    std::sort(queries.begin(),queries.end(),[](const auto &a,const auto &b){return a.id<b.id;});
    for (std::size_t i=0;i<queries.size();++i)
      {
        if (i && queries[i-1].id==queries[i].id) throw std::runtime_error("duplicate queryPointId");
        for (unsigned int d=0;d<dim;++d)
          if (!std::isfinite(queries[i].point[d])) throw std::runtime_error("non-finite query coordinate");
      }

    const auto faces=interface_faces();
    if (faces.empty()) throw std::runtime_error("configured interface boundary has no faces");
    std::vector<std::pair<dealii::BoundingBox<dim>,std::size_t>> indexed_faces;
    double maximum_face_tolerance=0;
    for (std::size_t i=0;i<faces.size();++i)
      {
        indexed_faces.emplace_back(dealii::BoundingBox<dim>({faces[i].lower,faces[i].upper}),i);
        maximum_face_tolerance=std::max(maximum_face_tolerance,scaled_tolerance(faces[i].diameter));
      }
    const auto face_index=dealii::pack_rtree(indexed_faces);
    Operator<dim> result; result.coordinate_frame=frame; result.geometry_version=geometry_version;
    std::ostringstream query_text,geometry_text,support_text,graph_text,weights_text;
    query_text << std::hexfloat; geometry_text << std::hexfloat;
    for (const auto &q:queries) { query_text << q.id; for (unsigned int d=0;d<dim;++d) query_text << ' ' << q.point[d]; query_text << '\n'; }
    for (const auto &face:faces)
      { geometry_text << face.cell->id().to_string() << ' ' << face.face << ' ' << face.diameter;
        for (unsigned int v=0;v<face.cell->face(face.face)->n_vertices();++v)
          for (unsigned int d=0;d<dim;++d) geometry_text << ' ' << face.cell->face(face.face)->vertex(v)[d];
        geometry_text << '\n'; }

    for (const auto &query:queries)
      {
        std::vector<Candidate> candidates;
        double minimum_distance=std::numeric_limits<double>::max();
        std::string closest_cell="none"; unsigned int closest_face=0;
        double closest_h=0,closest_tolerance=0;
        dealii::Point<dim> query_lower=query.point,query_upper=query.point;
        for (unsigned int d=0;d<dim;++d)
          { query_lower[d]-=maximum_face_tolerance; query_upper[d]+=maximum_face_tolerance; }
        std::vector<std::pair<dealii::BoundingBox<dim>,std::size_t>> spatial_candidates;
        face_index.query(boost::geometry::index::intersects(dealii::BoundingBox<dim>({query_lower,query_upper})),
                         std::back_inserter(spatial_candidates));
        std::sort(spatial_candidates.begin(),spatial_candidates.end(),[](const auto &a,const auto &b)
          { return a.second<b.second; });
        for (const auto &indexed_face:spatial_candidates)
          {
            const auto &face=faces[indexed_face.second];
            const double tolerance=scaled_tolerance(face.diameter);
            double box_distance=0;
            for (unsigned int d=0;d<dim;++d)
              box_distance+=std::pow(std::max({face.lower[d]-query.point[d],0.0,query.point[d]-face.upper[d]}),2);
            if (std::sqrt(box_distance)<minimum_distance)
              { minimum_distance=std::sqrt(box_distance); closest_cell=face.cell->id().to_string();
                closest_face=face.face; closest_h=face.diameter; closest_tolerance=tolerance; }
            if (box_distance>tolerance*tolerance) continue;
            try
              {
                const auto raw_unit=mapping.transform_real_to_unit_cell(face.cell,query.point);
                auto unit=raw_unit;
                if (dealii::GeometryInfo<dim>::distance_to_unit_cell(unit)>tolerance/std::max(face.diameter,1e-300)) continue;
                const unsigned int normal=dealii::GeometryInfo<dim>::unit_normal_direction[face.face];
                const double target=dealii::GeometryInfo<dim>::unit_normal_orientation[face.face]<0 ? 0.0 : 1.0;
                if (std::abs(unit[normal]-target)>tolerance/std::max(face.diameter,1e-300)) continue;
                // The inverse MappingQ1 on highly anisotropic cells can leave
                // coordinates of an already accepted edge/vertex a few ulps
                // times the local mapping condition away from 0 or 1. Snap in
                // reference space only; physical membership is rechecked.
                const double condition=face.diameter/std::max(face.minimum_edge,1e-300);
                const double snap_tolerance=tolerance/std::max(face.minimum_edge,1e-300)
                  +512.0*std::numeric_limits<double>::epsilon()*std::max(1.0,condition);
                for (unsigned int d=0;d<dim;++d)
                  {
                    const double endpoint=std::abs(unit[d])<=std::abs(unit[d]-1.0) ? 0.0 : 1.0;
                    const double change=std::abs(unit[d]-endpoint);
                    if (change<=snap_tolerance)
                      {
                        auto trial=unit; trial[d]=endpoint;
                        if (mapping.transform_unit_to_real_cell(face.cell,trial).distance(query.point)<=tolerance)
                          { unit=trial; if (change>0) { ++result.snapped_coordinates;
                              result.maximum_reference_snap=std::max(result.maximum_reference_snap,change); } }
                      }
                  }
                const auto mapped=mapping.transform_unit_to_real_cell(face.cell,unit);
                const double residual=mapped.distance(query.point);
                if (residual>tolerance) continue;
                candidates.push_back(make_candidate(query,face,raw_unit,unit,residual));
              }
            catch (const std::exception &) {}
          }
        if (candidates.empty())
          { std::ostringstream error; error << std::setprecision(17)
              << "queryPointId " << query.id << " not on boundary " << boundary
              << "; x=" << query.point << "; minBoxDistance=" << minimum_distance;
            error << "; tolerance=" << closest_tolerance << "; hFace=" << closest_h
                  << "; relativeDistance=" << minimum_distance/std::max(closest_h,1e-300)
                  << "; closestCellId=" << closest_cell << "; face=" << closest_face;
            throw std::runtime_error(error.str()); }
        std::sort(candidates.begin(),candidates.end(),[](const auto &a,const auto &b)
          { return a.cell_id<b.cell_id || (a.cell_id==b.cell_id && a.face_index<b.face_index); });
        const char *diagnostic_query=std::getenv("RUNTIME_Q1_DIAGNOSTIC_QUERY");
        if (diagnostic_query && query.id==std::stoull(diagnostic_query))
          {
            std::cerr << std::setprecision(17) << "RUNTIME_Q1_CANDIDATES query=" << query.id
                      << " point=" << query.point << " count=" << candidates.size() << '\n';
            for (std::size_t k=0;k<candidates.size();++k)
              {
                const auto &candidate=candidates[k];
                std::cerr << "candidate=" << k << " cell=" << candidate.cell_id
                          << " face=" << candidate.face_index << " rawXi=" << candidate.raw_reference
                          << " canonicalXi=" << candidate.reference
                          << " residual=" << candidate.residual
                          << " diameter=" << candidate.face_diameter
                          << " minEdge=" << candidate.minimum_edge << '\n';
                for (unsigned int component=0;component<dim;++component)
                  {
                    double sum=0; std::cerr << " component=" << component << " row=";
                    for (const auto &entry:candidate.component_entries[component])
                      { sum+=entry.weight; std::cerr << '(' << entry.dof << ',' << entry.weight << ')'; }
                    std::cerr << " sum=" << sum << '\n';
                  }
              }
          }
        if (candidates.size()>1)
          {
            ++result.multiple_candidates;
            for (std::size_t i=1;i<candidates.size();++i)
              {
                if (!equivalent(candidates[0],candidates[i]))
                  throw std::runtime_error("incompatible interface rows at edge/vertex for queryPointId "+std::to_string(query.id));
                ++result.equivalent_candidate_count;
              }
          }
        unsigned int endpoint_coordinates=0;
        for (unsigned int d=0;d<dim;++d)
          if (candidates[0].reference[d]==0.0 || candidates[0].reference[d]==1.0) ++endpoint_coordinates;
        if (endpoint_coordinates==1) ++result.face_interior_points;
        else if (endpoint_coordinates==2) ++result.edge_points;
        else if (endpoint_coordinates==3) ++result.vertex_points;
        for (unsigned int component=0;component<dim;++component)
          {
            Row<dim> row; row.query_id=query.id; row.component=component;
            row.cell_id=candidates[0].cell_id; row.face_index=candidates[0].face_index;
            row.reference=candidates[0].reference; row.residual=candidates[0].residual;
            row.entries=candidates[0].component_entries[component];
            for (const auto &entry:row.entries) row.row_sum+=entry.weight;
            if (std::abs(row.row_sum-1.0)>1e-14)
              { std::ostringstream error; error << std::setprecision(17)
                  << "RuntimeQ1 row does not form partition of unity: query=" << row.query_id
                  << " component=" << component << " sum=" << row.row_sum
                  << " cell=" << row.cell_id << " face=" << row.face_index
                  << " xi=" << row.reference; throw std::runtime_error(error.str()); }
            result.maximum_residual=std::max(result.maximum_residual,row.residual);
            support_text << row.query_id << ' ' << row.component << ' ' << row.cell_id << ' ' << row.face_index << '\n';
            graph_text << row.query_id << ' ' << row.component;
            weights_text << std::hexfloat << row.query_id << ' ' << row.component;
            for (const auto &entry:row.entries) { graph_text << ' ' << entry.dof; weights_text << ' ' << entry.dof << ' ' << entry.weight; }
            graph_text << '\n'; weights_text << '\n'; result.rows.push_back(std::move(row));
          }
      }
    result.query_manifest_hash=hexadecimal(fnv(query_text.str()));
    result.geometry_hash=hexadecimal(fnv(geometry_text.str()));
    result.support_hash=hexadecimal(fnv(support_text.str()));
    result.graph_hash=hexadecimal(fnv(graph_text.str()));
    result.weights_hash=hexadecimal(fnv(weights_text.str()));
    return result;
  }

private:
  using Cell=typename dealii::DoFHandler<dim>::active_cell_iterator;
  struct Face { Cell cell; unsigned int face; dealii::Point<dim> lower,upper; double diameter,minimum_edge; };
  struct Candidate
  {
    std::string cell_id; unsigned int face_index=0; dealii::Point<dim> raw_reference,reference;
    double residual=0,face_diameter=0,minimum_edge=0;
    std::vector<std::vector<Entry>> component_entries;
  };

  double scaled_tolerance(const double h) const
  { return std::max(128.0*std::numeric_limits<double>::epsilon()*std::max(1.0,h),relative_tolerance*h); }

  std::vector<Face> interface_faces() const
  {
    std::vector<Face> faces;
    for (const auto &cell:dofs.active_cell_iterators())
      for (unsigned int f=0;f<dealii::GeometryInfo<dim>::faces_per_cell;++f)
        if (cell->face(f)->at_boundary() && cell->face(f)->boundary_id()==boundary)
          {
            double minimum_edge=std::numeric_limits<double>::max();
            for (unsigned int e=0;e<cell->face(f)->n_lines();++e)
              minimum_edge=std::min(minimum_edge,cell->face(f)->line(e)->diameter());
            Face face{cell,f,cell->face(f)->vertex(0),cell->face(f)->vertex(0),
                      cell->face(f)->diameter(),minimum_edge};
            for (unsigned int v=1;v<cell->face(f)->n_vertices();++v)
              for (unsigned int d=0;d<dim;++d)
                { face.lower[d]=std::min(face.lower[d],cell->face(f)->vertex(v)[d]);
                  face.upper[d]=std::max(face.upper[d],cell->face(f)->vertex(v)[d]); }
            faces.push_back(face);
          }
    std::sort(faces.begin(),faces.end(),[](const auto &a,const auto &b)
      { return a.cell->id().to_string()<b.cell->id().to_string() ||
          (a.cell->id().to_string()==b.cell->id().to_string() && a.face<b.face); });
    return faces;
  }

  Candidate make_candidate(const Query<dim> &,const Face &face,
                           const dealii::Point<dim> &raw_unit,
                           const dealii::Point<dim> &unit,const double residual) const
  {
    Candidate candidate; candidate.cell_id=face.cell->id().to_string();
    candidate.face_index=face.face; candidate.raw_reference=raw_unit;
    candidate.reference=unit; candidate.residual=residual;
    candidate.face_diameter=face.diameter; candidate.minimum_edge=face.minimum_edge;
    candidate.component_entries.resize(dim);
    std::vector<global_dof_index> local(dofs.get_fe().n_dofs_per_cell());
    face.cell->get_dof_indices(local);
    // Keep every numerically non-zero FE value.  Inverse MappingQ1 can place
    // an on-face point a few ulps off the reference face; pruning those small
    // but legitimate Q1 values breaks partition of unity.
    const double prune=0.0;
    for (unsigned int i=0;i<local.size();++i)
      for (unsigned int component=0;component<dim;++component)
        {
          const double value=dofs.get_fe().shape_value_component(i,unit,component);
          if (std::abs(value)>prune) candidate.component_entries[component].push_back({local[i],value});
        }
    for (auto &entries:candidate.component_entries)
      std::sort(entries.begin(),entries.end(),[](const auto &a,const auto &b){return a.dof<b.dof;});
    return candidate;
  }

  void expand_constraint(const global_dof_index id,const double weight,
                         std::map<global_dof_index,double> &linear,double &offset,
                         std::vector<global_dof_index> &stack) const
  {
    if (!constraints.is_constrained(id)) { linear[id]+=weight; return; }
    if (std::find(stack.begin(),stack.end(),id)!=stack.end())
      throw std::runtime_error("cyclic AffineConstraints in RuntimeQ1 equivalence");
    offset+=weight*constraints.get_inhomogeneity(id);
    const auto *entries=constraints.get_constraint_entries(id);
    if (!entries) return;
    stack.push_back(id);
    for (const auto &entry:*entries)
      expand_constraint(entry.first,weight*entry.second,linear,offset,stack);
    stack.pop_back();
  }

  std::pair<std::map<global_dof_index,double>,double>
  constrained_functional(const std::vector<Entry> &row) const
  {
    std::map<global_dof_index,double> linear; double offset=0;
    std::vector<global_dof_index> stack;
    for (const auto &entry:row) expand_constraint(entry.dof,entry.weight,linear,offset,stack);
    for (auto it=linear.begin();it!=linear.end();)
      if (it->second==0) it=linear.erase(it); else ++it;
    return {linear,offset};
  }

  bool equivalent(const Candidate &a,const Candidate &b) const
  {
    const double condition=std::max(a.face_diameter/std::max(a.minimum_edge,1e-300),
                                    b.face_diameter/std::max(b.minimum_edge,1e-300));
    const double roundoff_bound=512.0*std::numeric_limits<double>::epsilon()*std::max(1.0,condition);
    for (unsigned int c=0;c<a.component_entries.size();++c)
      {
        if (a.component_entries[c].size()!=b.component_entries[c].size()) return false;
        for (std::size_t i=0;i<a.component_entries[c].size();++i)
          if (a.component_entries[c][i].dof!=b.component_entries[c][i].dof) return false;
        std::map<global_dof_index,double> difference;
        for (const auto &e:a.component_entries[c]) difference[e.dof]+=e.weight;
        for (const auto &e:b.component_entries[c]) difference[e.dof]-=e.weight;
        for (const auto &e:difference) if (std::abs(e.second)>roundoff_bound) return false;
        const auto ca=constrained_functional(a.component_entries[c]);
        const auto cb=constrained_functional(b.component_entries[c]);
        std::map<global_dof_index,double> constrained_difference=ca.first;
        for (const auto &entry:cb.first) constrained_difference[entry.first]-=entry.second;
        for (const auto &entry:constrained_difference)
          if (std::abs(entry.second)>roundoff_bound) return false;
        if (std::abs(ca.second-cb.second)>roundoff_bound) return false;
      }
    return true;
  }

  const dealii::MappingQ1<dim> &mapping;
  const dealii::DoFHandler<dim> &dofs;
  const dealii::AffineConstraints<double> &constraints;
  const dealii::types::boundary_id boundary;
  const std::string frame;
  const std::uint64_t geometry_version;
  const double relative_tolerance;
};
}
