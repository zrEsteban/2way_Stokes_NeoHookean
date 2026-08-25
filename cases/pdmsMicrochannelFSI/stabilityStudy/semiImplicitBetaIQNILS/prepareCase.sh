#!/usr/bin/env bash
set -euo pipefail

case_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

blockMesh -case "$case_dir" -region fluid
blockMesh -case "$case_dir" -region solid
(
    cd "$case_dir/dealiiSolid"
    gmsh solid.geo -3 -format msh2 -o solid.msh
)

test -s "$case_dir/constant/fluid/polyMesh/faces"
test -s "$case_dir/constant/solid/polyMesh/faces"
test -s "$case_dir/dealiiSolid/solid.msh"
echo "Prepared OpenFOAM and deal.II meshes under $case_dir"
