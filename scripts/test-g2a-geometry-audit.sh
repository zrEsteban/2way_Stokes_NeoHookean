#!/usr/bin/env bash
set -euo pipefail
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
case_dir="$repo_root/cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS"
work=${G2A_GEOMETRY_WORK:-$(mktemp -d /tmp/g2a-geometry.XXXXXX)}
mkdir -p "$work"
run()
{
  local name=$1; shift
  "$@" python3 "$repo_root/scripts/audit_g2a_interface_geometry.py" "$case_dir" \
    --json "$work/$name.json" --csv "$work/$name.csv" >"$work/$name.log"
}
run serial env
run repeat env
run mpi2 mpirun -np 2
run mpi4 mpirun -np 4
python3 - "$work" <<'PY'
import json,sys
from pathlib import Path
w=Path(sys.argv[1]); values=[json.loads((w/f'{n}.json').read_text()) for n in ('serial','repeat','mpi2','mpi4')]
assert {v['classification'] for v in values}=={'CAD_MISMATCH'}
assert len({v['auditHash'] for v in values})==1
v=values[0]
assert v['provenance']['fixtureExactSetMatch'] and not v['provenance']['fixtureGlobalOrderMatch']
assert v['provenance']['fixturePatchLocalOrderMatch']
assert v['searchComparison']['differentDistance']==0
assert v['outsideTolerance']['1e-8']==162
assert v['query20']['global_fluid_point_id']==20251
assert abs(v['query20']['distance']-4.1667767476470444e-7)<1e-18
print('G2-A geometry audit PASS:',v['classification'],v['auditHash'])
PY
echo "artifacts: $work"
