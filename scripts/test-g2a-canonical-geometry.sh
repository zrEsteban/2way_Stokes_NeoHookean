#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
work=${G2A_CANONICAL_WORK:-$(mktemp -d /tmp/g2a-canonical.XXXXXX)}
mkdir -p "$work"
python3 "$root/scripts/generate_canonical_interface_geometry.py" --output "$work/run1" >"$work/run1.log"
python3 "$root/scripts/generate_canonical_interface_geometry.py" --output "$work/run2" >"$work/run2.log"
python3 - "$work" <<'PY'
import json,sys
from pathlib import Path
w=Path(sys.argv[1]); a=json.loads((w/'run1/generation-manifest.json').read_text()); b=json.loads((w/'run2/generation-manifest.json').read_text())
for key in ('geometryContractHash','solid','fluid','masterCrossSectionPoints','manifestHash'):
    assert a[key]==b[key],(key,a[key],b[key])
print('deterministic generation PASS',a['manifestHash'])
PY
python3 "$root/scripts/validate_canonical_interface_geometry.py" "$work/run1" --json "$work/geometry-validation.json"
checkMesh -no-libs -case "$work/run1/convert" -allGeometry -allTopology >"$work/checkMesh.log" 2>&1
if grep -q 'Mesh OK' "$work/checkMesh.log"; then
  echo "unexpected Mesh OK: update this blocked-gate regression" >&2; exit 1
fi
grep -q 'Cells with small determinant.*number of cells: 652' "$work/checkMesh.log"
grep -q 'Failed 1 mesh checks' "$work/checkMesh.log"
runtime=${RUNTIME_Q1_TEST:-/tmp/g2a-runtime-q1-build/testRuntimeQ1}
if "$runtime" "$work/run1/solid.msh" "$work/run1/robin-query.csv" >"$work/runtime-serial.log" 2>&1; then
  echo "unexpected RuntimeQ1 PASS after quality blocker" >&2; exit 1
fi
grep -q 'incompatible interface rows at edge/vertex' "$work/runtime-serial.log"
echo "G2-A canonicalization audit PASS; gate BLOCKED_MESH_QUALITY: $work"
