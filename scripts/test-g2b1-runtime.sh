#!/usr/bin/env bash
set -euo pipefail
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
case_dir="$repo_root/cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS"
if [[ -n ${G2B1_RUNTIME_WORK:-} ]]; then
  work=$G2B1_RUNTIME_WORK
  [[ ! -e "$work" ]] || { echo "ERROR: work directory exists: $work" >&2; exit 2; }
  mkdir -p "$work"
else
  work=$(mktemp -d /tmp/g2b1-runtime.XXXXXX)
fi
python3 "$repo_root/scripts/audit-g2b1-sparsity.py" \
  --export-manifest "$work/interface.manifest" "$case_dir" >"$work/reference.json"
cmake -S "$repo_root/src/dealiiPdmsSolid" -B "$work/build" -DCMAKE_BUILD_TYPE=Release
cmake --build "$work/build" --target testRuntimeNewton dealiiPdmsSolid -j2
mesh="$case_dir/dealiiSolid/solid.msh"
"$work/build/testRuntimeNewton" "$mesh" "$work/interface.manifest" >"$work/serial.log" 2>&1
mpirun -np 2 "$work/build/testRuntimeNewton" "$mesh" "$work/interface.manifest" >"$work/mpi2.log" 2>&1
mpirun -np 4 "$work/build/testRuntimeNewton" "$mesh" "$work/interface.manifest" >"$work/mpi4.log" 2>&1
python3 - "$work" <<'PY'
import json,pathlib,sys
root=pathlib.Path(sys.argv[1]); values=[]
for name in ('serial.log','mpi2.log','mpi4.log'):
    line=next(line for line in reversed((root/name).read_text().splitlines()) if line.startswith('{'))
    values.append(json.loads(line))
for ranks,value in zip((1,2,4),values):
    assert value['ranks']==ranks and value['unit_rhs']==1
    assert value['repeat_error']==0 and value['runtime_fd_best']<=1e-8
print('G2-B.1-RUNTIME PdmsSolid::assemble_newton PASS')
for value in values: print(json.dumps(value,sort_keys=True))
PY
echo "artifacts: $work"
