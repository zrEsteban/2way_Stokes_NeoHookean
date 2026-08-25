#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source_case_rel=cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS
tmp_root=$(mktemp -d /tmp/g0-constantHs-regression.XXXXXX)
trap 'rm -rf "$tmp_root"' EXIT
case_dir="$tmp_root/case"

command -v foamDictionary >/dev/null || {
    echo "ERROR: initialise OpenFOAM before running this test" >&2
    exit 2
}

mkdir -p "$case_dir"
git -C "$repo_root" archive HEAD:"$source_case_rel" | tar -x -C "$case_dir"
cp "$repo_root/$source_case_rel/decomposeCase.sh" "$case_dir/decomposeCase.sh"

"$case_dir/prepareCase.sh" >/dev/null
TWO_WAY_FSI_ROOT="$repo_root" "$case_dir/decomposeCase.sh" >/dev/null
"$repo_root/scripts/validate-decomposition-constantHs.sh" "$case_dir"

if TWO_WAY_FSI_ROOT="$repo_root" "$case_dir/decomposeCase.sh" >/dev/null 2>&1; then
    echo "ERROR: decomposition accepted pre-existing processor directories" >&2
    exit 1
fi

field="$case_dir/processor0/0/fluid/p"
entry=boundaryField.interface.constantHs
expected=$(foamDictionary -no-libs -precision 17 "$field" -entry "$entry" -value)

foamDictionary -no-libs "$field" -entry "$entry" -remove >/dev/null
if "$repo_root/scripts/validate-decomposition-constantHs.sh" "$case_dir" >/dev/null 2>&1; then
    echo "ERROR: validator accepted a missing constantHs" >&2
    exit 1
fi

foamDictionary -no-libs "$field" -entry "$entry" -set 1e-12 >/dev/null
if "$repo_root/scripts/validate-decomposition-constantHs.sh" "$case_dir" >/dev/null 2>&1; then
    echo "ERROR: validator accepted an altered constantHs" >&2
    exit 1
fi

foamDictionary -no-libs -precision 17 "$field" -entry "$entry" \
    -set "$expected" >/dev/null
"$repo_root/scripts/validate-decomposition-constantHs.sh" "$case_dir" >/dev/null

echo "decomposition constantHs regression PASS"
