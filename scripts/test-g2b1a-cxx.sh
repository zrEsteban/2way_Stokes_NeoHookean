#!/usr/bin/env bash
set -euo pipefail
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
case_dir="$repo_root/cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS"
if [[ -n ${G2B1A_CXX_WORK:-} ]]; then
  work=$G2B1A_CXX_WORK
  [[ ! -e "$work" ]] || { echo "ERROR: work directory exists: $work" >&2; exit 2; }
  mkdir -p "$work"
else
  work=$(mktemp -d /tmp/g2b1a-cxx.XXXXXX)
fi
mkdir -p "$work/build"
python3 "$repo_root/scripts/audit-g2b1-sparsity.py" \
  --expect-legacy-missing --export-manifest "$work/interface.manifest" "$case_dir" \
  >"$work/python-reference.json"
cmake -S "$repo_root/src/dealiiPdmsSolid" -B "$work/build" -DCMAKE_BUILD_TYPE=Release
cmake --build "$work/build" --target testInterfaceSparsity dealiiPdmsSolid -j2
test_exe="$work/build/testInterfaceSparsity"
mesh="$case_dir/dealiiSolid/solid.msh"
"$test_exe" "$mesh" "$work/interface.manifest" >"$work/serial.json"
mpirun -np 2 "$test_exe" "$mesh" "$work/interface.manifest" >"$work/mpi2.json"
mpirun -np 4 "$test_exe" "$mesh" "$work/interface.manifest" >"$work/mpi4.json"

# Negative protocol/lifecycle checks: the production dual branch requires its
# manifest before SparseMatrix::reinit, and the parser rejects identity/ID drift.
sed -e "s|^set mesh.*|set mesh = $mesh|" \
    -e 's|^set input.*|set input = absent.csv|' \
    "$repo_root/src/dealiiPdmsSolid/tests/parameters.prm" >"$work/dual-missing.prm"
printf '\nset interface transfer = dualConservative\n' >>"$work/dual-missing.prm"
if "$work/build/dealiiPdmsSolid" "$work/dual-missing.prm" >"$work/missing.out" 2>&1; then
  echo 'ERROR: dual setup accepted an absent manifest' >&2; exit 1
fi
rg -q 'requires interface manifest before matrix reinit' "$work/missing.out"
sed 's/6ee94c285c4e5d8619f9f08c3a666eb57dd46ddfe47b2b02be25c1248b4abc39/0000000000000000000000000000000000000000000000000000000000000000/' \
  "$work/interface.manifest" >"$work/bad-hash.manifest"
if "$test_exe" "$mesh" "$work/bad-hash.manifest" >"$work/bad-hash.out" 2>&1; then
  echo 'ERROR: C++ test accepted an incompatible hashGraph' >&2; exit 1
fi
rg -q 'hashGraph does not match canonical H graph' "$work/bad-hash.out"
sed '0,/^node 0 /s//node 999999 /' "$work/interface.manifest" >"$work/bad-id.manifest"
if "$test_exe" "$mesh" "$work/bad-id.manifest" >"$work/bad-id.out" 2>&1; then
  echo 'ERROR: C++ parser accepted an invalid node ID' >&2; exit 1
fi
rg -q 'Invalid node ID/order' "$work/bad-id.out"
python3 - "$work" <<'PY'
import json,pathlib,sys
root=pathlib.Path(sys.argv[1])
values=[json.loads((root/name).read_text()) for name in ('serial.json','mpi2.json','mpi4.json')]
reference=values[0]
symbolic=json.loads((root/'python-reference.json').read_text())
for value in values:
    assert value['backend']=='dealii::SparseMatrix<double> replicated'
    assert value['legacy_nnz']==reference['legacy_nnz']
    assert value['required']==reference['required']
    assert value['added']==symbolic['added_vector_nnz_scalar_or_diagonal_Z']
    assert value['missing']==0
    assert value['final_nnz']==reference['final_nnz']
    assert value['effective_hash_fnv64']==reference['effective_hash_fnv64']
    limit=1e-12 if value['ranks']==1 else 1e-10
    assert value['action_error']<=limit and value['symmetry_error']<=limit
    assert value['force_error']<=limit and value['protocol_kx_error']<=limit
    assert value['work_error']<=limit
    assert value['best_fd_error'] <= (1e-8 if value['ranks']==1 else 1e-7)
    assert value['full_fd_best'] <= (1e-8 if value['ranks']==1 else 1e-7)
    assert value['inverted_sign_fd_error'] > 1.9
print('G2-B.1a-CXX real DynamicSparsityPattern/SparseMatrix PASS')
for value in values: print(json.dumps(value,sort_keys=True))
PY
echo "artifacts: $work"
