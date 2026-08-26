#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
case_dir="$repo_root/cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS"
work=${G2A_WORK:-$(mktemp -d /tmp/g2a-dual-transfer.XXXXXX)}
[[ ! -e "$work" ]] || { echo "ERROR: work directory exists: $work" >&2; exit 2; }
mkdir -p "$work"

python3 "$repo_root/scripts/g2a_dual_transfer.py" "$case_dir" >"$work/serial.json"
mpirun -np 2 python3 "$repo_root/scripts/g2a_dual_transfer.py" "$case_dir" >"$work/mpi2.json"
mpirun -np 4 python3 "$repo_root/scripts/g2a_dual_transfer.py" "$case_dir" >"$work/mpi4.json"

python3 - "$work" <<'PY'
import json, pathlib, sys
root=pathlib.Path(sys.argv[1])
data=[json.loads((root/name).read_text()) for name in ('serial.json','mpi2.json','mpi4.json')]
for result in data:
    limit=1e-12 if result['ranks']==1 else 1e-10
    assert result['work_error'] <= limit, result
    assert result['partition_unity_error'] <= 5e-13, result
    assert result['force_error'] <= 1e-12, result
    assert result['uniform_resultant_abs_error'] <= 1e-12, result
keys=('rows_scalar','columns_scalar','nonzeros_scalar','area_sum','row_sum_max_error',
      'work_error','force_error','moment_error','linear_max_abs_error_m',
      'rotation_max_abs_error_m','max_surface_projection_error_m')
reference=data[0]
for result in data[1:]:
    for key in keys:
        scale=max(abs(reference[key]),1.0)
        assert abs(result[key]-reference[key]) <= 1e-12*scale, (key,reference,result)
print('G2-A dual transfer PASS')
for result in data: print(json.dumps(result,sort_keys=True))
PY
