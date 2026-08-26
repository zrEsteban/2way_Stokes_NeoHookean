#!/usr/bin/env python3
"""Diagnostic-only explicit H/H^T for the current deal.II--OpenFOAM interface."""

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
from mpi4py import MPI


def foam_points(path):
    text = path.read_text()
    return np.array([
        tuple(map(float, m)) for m in re.findall(
            r"\(([-+0-9.eE]+) ([-+0-9.eE]+) ([-+0-9.eE]+)\)", text
        )
    ])


def foam_faces(path):
    text = path.read_text()
    return [list(map(int, value.split())) for value in
            re.findall(r"\d+\(([^()]*)\)", text)]


def gmsh_volume(path, physical_id=10):
    lines = path.read_text().splitlines()
    begin = lines.index("$Nodes")
    count = int(lines[begin + 1])
    nodes = {int(a[0]): np.array(tuple(map(float, a[1:4])))
             for a in (line.split() for line in lines[begin + 2:begin + 2 + count])}
    begin = lines.index("$Elements")
    count = int(lines[begin + 1])
    cells = []
    for line in lines[begin + 2:begin + 2 + count]:
        a = list(map(int, line.split()))
        element_type, n_tags = a[1], a[2]
        tags, connectivity = a[3:3 + n_tags], a[3 + n_tags:]
        if element_type == 5 and tags and tags[0] == physical_id:
            cells.append(connectivity)
    ids = sorted(nodes)
    stable = {node: i for i, node in enumerate(ids)}
    return nodes, cells, stable


def shape(xi, eta, zeta):
    return np.array(((1-xi)*(1-eta)*(1-zeta), xi*(1-eta)*(1-zeta),
                     xi*eta*(1-zeta), (1-xi)*eta*(1-zeta),
                     (1-xi)*(1-eta)*zeta, xi*(1-eta)*zeta,
                     xi*eta*zeta, (1-xi)*eta*zeta))


def shape_derivatives(xi, eta, zeta):
    return (np.array((-(1-eta)*(1-zeta), (1-eta)*(1-zeta), eta*(1-zeta),
                      -eta*(1-zeta), -(1-eta)*zeta, (1-eta)*zeta, eta*zeta,
                      -eta*zeta)),
            np.array((-(1-xi)*(1-zeta), -xi*(1-zeta), xi*(1-zeta),
                      (1-xi)*(1-zeta), -(1-xi)*zeta, -xi*zeta, xi*zeta,
                      (1-xi)*zeta)),
            np.array((-(1-xi)*(1-eta), -xi*(1-eta), -xi*eta,
                      -(1-xi)*eta, (1-xi)*(1-eta), xi*(1-eta), xi*eta,
                      (1-xi)*eta)))


def locate_in_cell(point, coordinates):
    xi = eta = zeta = 0.5
    for _ in range(20):
        weights = shape(xi, eta, zeta)
        residual = weights @ coordinates - point
        dxi, deta, dzeta = shape_derivatives(xi, eta, zeta)
        jacobian = np.column_stack((dxi @ coordinates, deta @ coordinates,
                                    dzeta @ coordinates))
        step = np.linalg.solve(jacobian, residual)
        xi -= step[0]
        eta -= step[1]
        zeta -= step[2]
        if np.linalg.norm(step) < 1e-14:
            break
    weights = shape(xi, eta, zeta)
    distance = np.linalg.norm(weights @ coordinates - point)
    outside = max(0.0, -xi, xi-1.0, -eta, eta-1.0, -zeta, zeta-1.0)
    return distance + outside, weights


def polygon_geometry(coordinates):
    centre = coordinates.mean(axis=0)
    area_vector = np.zeros(3)
    for i in range(len(coordinates)):
        area_vector += np.cross(coordinates[i]-centre,
                                coordinates[(i+1) % len(coordinates)]-centre)/2
    return centre, np.linalg.norm(area_vector)


def build_operator(case):
    mesh = case / "constant/fluid/polyMesh"
    points = foam_points(mesh / "points")
    faces = foam_faces(mesh / "faces")
    boundary = (mesh / "boundary").read_text()
    match = re.search(r"interface\s*\{.*?nFaces\s+(\d+);.*?startFace\s+(\d+);",
                      boundary, re.S)
    if not match:
        raise RuntimeError("interface patch not found")
    n_faces, start_face = map(int, match.groups())
    patch_faces = faces[start_face:start_face+n_faces]
    patch_point_ids = sorted(set(sum(patch_faces, [])))

    nodes, cells, stable = gmsh_volume(case / "dealiiSolid/solid.msh")
    cell_data = []
    for cell in cells:
        # Gmsh hex ordering matches the Q1 vertex ordering used by deal.II.
        coordinates = np.array([nodes[node] for node in cell])
        cell_data.append((cell, coordinates, coordinates.mean(axis=0)))

    point_rows = {}
    max_projection = 0.0
    for point_id in patch_point_ids:
        point = points[point_id]
        candidates = sorted(cell_data, key=lambda q: np.dot(q[2]-point, q[2]-point))[:24]
        best = None
        for cell, coordinates, _ in candidates:
            score, weights = locate_in_cell(point, coordinates)
            if best is None or score < best[0]:
                best = score, cell, weights
        score, cell, weights = best
        max_projection = max(max_projection, score)
        point_rows[point_id] = {stable[node]: float(weight)
                               for node, weight in zip(cell, weights)}

    rows, centres, areas = [], [], []
    for face in patch_faces:
        row = {}
        for point_id in face:
            for column, weight in point_rows[point_id].items():
                row[column] = row.get(column, 0.0) + weight/len(face)
        rows.append(row)
        centre, area = polygon_geometry(points[face])
        centres.append(centre)
        areas.append(area)
    structural_points = np.array([nodes[node] for node, _ in sorted(stable.items(), key=lambda x:x[1])])
    return rows, np.array(centres), np.array(areas), structural_points, start_face, max_projection


def apply_h(rows, values, owned):
    result = np.zeros((len(rows), 3))
    for face in owned:
        for column, weight in rows[face].items():
            result[face] += weight*values[column]
    return result


def apply_ht(rows, values, n_columns, owned):
    result = np.zeros((n_columns, 3))
    for face in owned:
        for column, weight in rows[face].items():
            result[column] += weight*values[face]
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("case", type=Path)
    parser.add_argument("--seed", type=int, default=20260825)
    args = parser.parse_args()
    comm = MPI.COMM_WORLD
    rank, size = comm.rank, comm.size
    rows, centres, areas, structural_points, start_face, max_projection = build_operator(args.case)
    owned = range(rank, len(rows), size)
    n_columns = len(structural_points)
    rng = np.random.default_rng(args.seed)
    virtual = rng.standard_normal((n_columns, 3))
    traction = rng.standard_normal((len(rows), 3))*1.0e3
    face_virtual_local = apply_h(rows, virtual, owned)
    weighted_traction = areas[:, None]*traction
    nodal_local = -apply_ht(rows, weighted_traction, n_columns, owned)
    face_virtual = np.zeros_like(face_virtual_local)
    nodal = np.zeros_like(nodal_local)
    comm.Allreduce(face_virtual_local, face_virtual, op=MPI.SUM)
    comm.Allreduce(nodal_local, nodal, op=MPI.SUM)

    fluid_work_local = sum(np.dot(face_virtual[i], weighted_traction[i]) for i in owned)
    fluid_work = comm.allreduce(fluid_work_local, op=MPI.SUM)
    solid_work = float(np.sum(virtual*nodal))
    floor = 1e-30
    work_error = abs(fluid_work+solid_work)/max(abs(fluid_work), abs(solid_work), floor)
    row_errors = [abs(sum(row.values())-1.0) for row in rows]
    fluid_force = weighted_traction.sum(axis=0)
    force_error = np.linalg.norm(nodal.sum(axis=0)+fluid_force)/max(np.linalg.norm(fluid_force), floor)
    fluid_moment = np.cross(centres, weighted_traction).sum(axis=0)
    solid_moment = np.cross(structural_points, nodal).sum(axis=0)
    moment_scale = max(np.linalg.norm(fluid_moment), np.linalg.norm(solid_moment), floor)
    moment_error = np.linalg.norm(fluid_moment+solid_moment)/moment_scale

    ones = np.ones((n_columns, 3))
    rigid = np.zeros((len(rows), 3))
    comm.Allreduce(apply_h(rows, ones, owned), rigid, op=MPI.SUM)
    partition_error = float(np.max(np.abs(rigid-1.0)))

    linear = structural_points.copy()
    linear_faces = np.zeros((len(rows), 3))
    comm.Allreduce(apply_h(rows, linear, owned), linear_faces, op=MPI.SUM)
    linear_error = float(np.max(np.linalg.norm(linear_faces-centres, axis=1)))
    omega = np.array((0.7, -0.2, 0.4))
    rotation = np.cross(np.broadcast_to(omega, structural_points.shape), structural_points)
    rotation_faces = np.zeros((len(rows), 3))
    comm.Allreduce(apply_h(rows, rotation, owned), rotation_faces, op=MPI.SUM)
    expected_rotation = np.cross(np.broadcast_to(omega, centres.shape), centres)
    rotation_error = float(np.max(np.linalg.norm(rotation_faces-expected_rotation, axis=1)))

    uniform = np.tile(np.array((1200.0, -300.0, 450.0)), (len(rows), 1))
    uniform_force = -apply_ht(rows, areas[:, None]*uniform, n_columns, owned)
    uniform_global = np.zeros_like(uniform_force)
    comm.Allreduce(uniform_force, uniform_global, op=MPI.SUM)
    uniform_resultant_error = np.linalg.norm(
        uniform_global.sum(axis=0)+(areas[:, None]*uniform).sum(axis=0))

    if rank == 0:
        print(json.dumps({
            "ranks": size,
            "rows_scalar": len(rows),
            "columns_scalar": n_columns,
            "rows_vector": 3*len(rows),
            "columns_vector": 3*n_columns,
            "nonzeros_scalar": sum(len(row) for row in rows),
            "patch_start_face": start_face,
            "area_sum": float(areas.sum()),
            "row_sum_max_error": max(row_errors),
            "partition_unity_error": partition_error,
            "work_error": work_error,
            "force_error": float(force_error),
            "uniform_resultant_abs_error": float(uniform_resultant_error),
            "moment_error": float(moment_error),
            "linear_max_abs_error_m": linear_error,
            "rotation_max_abs_error_m": rotation_error,
            "max_surface_projection_error_m": max_projection,
            "ownership": "globalFaceIndex modulo ranks; H^T Allreduce to stable structural-node owner",
        }, sort_keys=True))


if __name__ == "__main__":
    main()
