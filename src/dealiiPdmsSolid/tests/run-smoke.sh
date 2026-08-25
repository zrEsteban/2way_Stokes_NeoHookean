#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
build="${DEALII_PDMS_BUILD:-/tmp/dealii-pdms-build}"
work="$(mktemp -d /tmp/dealii-pdms-smoke.XXXXXX)"
trap 'rm -rf "$work"' EXIT

cmake -S "$root" -B "$build" -DCMAKE_BUILD_TYPE=Release
cmake --build "$build" -j"${BUILD_JOBS:-2}"
ctest --test-dir "$build" --output-on-failure
gmsh "$root/tests/box.geo" -3 -format msh2 -o "$work/solid.msh" -v 1
cp "$root/tests/parameters.prm" "$work/parameters.prm"
{
  echo 'x,y,z,tx,ty,tz,vx,vy,vz'
  echo '0.0005,0.0001,0.0001,0,0,1000,0,0,0'
} > "$work/robin-in.csv"
(cd "$work" && "$build/dealiiPdmsSolid" parameters.prm)

test -s "$work/robin-out.csv"
test -s "$work/solid.vtu"
test -s "$work/solid-state.bin"
awk -F, 'NR>1 {if ($6 > max) max=$6} END {exit !(max>0)}' "$work/robin-out.csv"
echo "PASS: deal.II Robin load produces positive interface displacement"

# Advance one more accepted structural step from the serialized state. This
# exercises both vectors in the restart (displacement and velocity).
cp "$work/solid-state.bin" "$work/accepted-state.bin"
cp "$work/robin-out.csv" "$work/first-step.csv"
{
  echo 'set state input = accepted-state.bin'
  echo 'set state output = second-state.bin'
  echo 'set output = second-step.csv'
  echo 'set vtk output = second-step.vtu'
} >> "$work/parameters.prm"
(cd "$work" && "$build/dealiiPdmsSolid" parameters.prm)
test -s "$work/second-state.bin"
awk -F, 'NR>1 {for(i=4;i<=12;i++) if ($i!=$i) exit 1; n++} END {exit !(n>0)}' \
  "$work/second-step.csv"
cmp -s "$work/accepted-state.bin" "$work/second-state.bin" && {
  echo 'restart did not advance the structural state' >&2
  exit 1
}
echo "PASS: accepted displacement/velocity state restarts and advances"

# Exercise the same rounded PDMS topology and physical boundary tags used by
# the production mesh, at a deliberately small resolution.
(cd "$root/tests" && gmsh rounded_gate.geo -3 -format msh2 \
  -o "$work/rounded-solid.msh" -v 1)
cp "$root/tests/parameters.prm" "$work/rounded.prm"
{
  echo 'set mesh = rounded-solid.msh'
  echo 'set input = rounded-in.csv'
  echo 'set output = rounded-out.csv'
  echo 'set vtk output = rounded.vtu'
  echo 'set interface boundary = 4'
  echo 'set clamped boundary 1 = 1'
  echo 'set clamped boundary 2 = 2'
  echo 'set clamped boundary 3 = 3'
  echo 'set state input ='
  echo 'set state output = rounded-state.bin'
} >> "$work/rounded.prm"
{
  echo 'x,y,z,tx,ty,tz,vx,vy,vz'
  echo '0.0001,0.00005,0.00003,0,0,1000,0,0,0'
} > "$work/rounded-in.csv"
(cd "$work" && "$build/dealiiPdmsSolid" rounded.prm)
awk -F, 'NR>1 {if ($6 > max) max=$6} END {exit !(max>0)}' "$work/rounded-out.csv"
echo "PASS: rounded PDMS topology and boundary ids solve with Robin loading"

# Verify the two limiting behaviours of the solid Robin operator on the same
# mesh and load.  Z=0 is the exact Neumann reference; increasing Z with a zero
# target velocity must drive the interface velocity toward the Dirichlet
# limit.  The reported traction must also close the Robin identity pointwise.
for label in neumann matched dirichlet; do
  case "$label" in
    neumann) impedance=0 ;;
    matched) impedance=67350 ;;
    dirichlet) impedance=6.735e10 ;;
  esac
  cp "$root/tests/parameters.prm" "$work/limit-$label.prm"
  {
    echo "set solid impedance = $impedance"
    echo "set output = limit-$label.csv"
    echo "set vtk output = limit-$label.vtu"
    echo 'set state output ='
  } >> "$work/limit-$label.prm"
  (cd "$work" && "$build/dealiiPdmsSolid" "limit-$label.prm")
done

awk -F, '
  FNR==1 {next}
  {v=sqrt($7*$7+$8*$8+$9*$9); if (v>max[FILENAME]) max[FILENAME]=v}
  END {
    n=max[ARGV[1]]; m=max[ARGV[2]]; d=max[ARGV[3]];
    if (!(n>m && m>d && d/n<1e-5)) exit 1
  }
' "$work/limit-neumann.csv" "$work/limit-matched.csv" \
  "$work/limit-dirichlet.csv"

awk -F, -v Z=67350 '
  NR==1 {next}
  {
    for (d=0; d<3; ++d) {
      # Input traction is (0,0,1000) and target velocity is zero.
      expected=(d==2 ? 1000 : 0)-Z*$(7+d);
      error=$(10+d)-expected;
      scale=1+sqrt(expected*expected);
      if (sqrt(error*error)/scale>1e-10) exit 1
    }
    rows++
  }
  END {exit !(rows>0)}
' "$work/limit-matched.csv"
echo "PASS: Robin operator approaches Neumann/Dirichlet limits and closes pointwise"
