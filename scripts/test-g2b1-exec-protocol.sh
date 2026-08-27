#!/usr/bin/env bash
set -euo pipefail
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
case_dir="$repo_root/cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS"
work=${G2B1_EXEC_WORK:-$(mktemp -d /tmp/g2b1-exec.XXXXXX)}
[[ ! -e "$work" ]] || [[ -d "$work" && -z "$(find "$work" -mindepth 1 -print -quit)" ]] || {
  echo "ERROR: work directory is not empty: $work" >&2; exit 2; }
mkdir -p "$work"
python3 "$repo_root/scripts/audit-g2b1-sparsity.py" \
  --export-manifest "$work/interface.manifest" "$case_dir" >"$work/reference.json"
cmake -S "$repo_root/src/dealiiPdmsSolid" -B "$work/build" -DCMAKE_BUILD_TYPE=Release
cmake --build "$work/build" --target dealiiPdmsSolid testExecutableProtocol testRuntimeNewton -j2
"$work/build/testExecutableProtocol"
participant=(python3 "$repo_root/scripts/g2b1_exec_participant.py"
  --binary "$work/build/dealiiPdmsSolid" --mesh "$case_dir/dealiiSolid/solid.msh"
  --manifest "$work/interface.manifest")
for ranks in 1 2 4; do
  "${participant[@]}" --work "$work/valid-r${ranks}" --ranks "$ranks"
done
negative=(bad_magic bad_schema unknown_type oversize bad_checksum truncated timeout
  nan inf duplicate_id invalid_id invalid_tangent_id outside_pattern force_pa tangent_units
  wrong_graph wrong_weights wrong_dof_hash stale future sequence_collision force_only
  tangent_only version_mismatch bad_activation partial_disconnect)
for test_name in "${negative[@]}"; do
  python3 "$repo_root/scripts/g2b1_exec_negative.py" "$test_name" \
    --binary "$work/build/dealiiPdmsSolid" --mesh "$case_dir/dealiiSolid/solid.msh" \
    --manifest "$work/interface.manifest" --work "$work/negative"
done
echo "G2-B.1-EXEC-PROTOCOL PASS"
echo "artifacts: $work"
