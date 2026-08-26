#!/usr/bin/env bash
set -euo pipefail
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
case_dir="$repo_root/cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS"
work=${G2B1A_WORK:-$(mktemp -d /tmp/g2b1a-sparsity.XXXXXX)}
[[ ! -e "$work" ]] || { echo "ERROR: work directory exists: $work" >&2; exit 2; }
mkdir -p "$work"
python3 "$repo_root/scripts/audit-g2b1-sparsity.py" --expect-legacy-missing "$case_dir" >"$work/serial.json"
mpirun -np 2 python3 "$repo_root/scripts/audit-g2b1-sparsity.py" --expect-legacy-missing "$case_dir" >"$work/mpi2.json"
mpirun -np 4 python3 "$repo_root/scripts/audit-g2b1-sparsity.py" --expect-legacy-missing "$case_dir" >"$work/mpi4.json"
python3 - "$work" <<'PY'
import json,pathlib,sys
root=pathlib.Path(sys.argv[1]); values=[json.loads((root/n).read_text()) for n in ('serial.json','mpi2.json','mpi4.json')]
reference=values[0]; assert reference['raw_missing_scalar_pairs']==13124; assert reference['missing_after_extension']==0
for value in values:
    assert value['hashGraph']==reference['hashGraph']
    assert value['effective_missing_scalar_pairs']==reference['effective_missing_scalar_pairs']
    assert all(value['negative_tests'].values())
    limit=1e-12 if value['ranks']==1 else 1e-10
    for result in value['results']:
        assert result['action_error']<=limit and result['symmetry_error']<=limit and result['energy']>=-1e-13
print('G2-B.1a sparse graph PASS')
for value in values: print(json.dumps(value,sort_keys=True))
PY
