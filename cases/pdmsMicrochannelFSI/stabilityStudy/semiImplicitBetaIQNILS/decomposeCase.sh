#!/usr/bin/env bash
set -euo pipefail

case_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=${TWO_WAY_FSI_ROOT:-$(cd "$case_dir/../../../.." && pwd)}
[[ -x "$repo_root/scripts/validate-decomposition-constantHs.sh" ]] || {
    echo "ERROR: set TWO_WAY_FSI_ROOT to the project repository root" >&2
    exit 2
}

command -v decomposePar >/dev/null || {
    echo "ERROR: decomposePar is unavailable; initialise OpenFOAM first" >&2
    exit 2
}

shopt -s nullglob
existing=("$case_dir"/processor*)
if (( ${#existing[@]} > 0 )); then
    echo "ERROR: processor directories already exist; use a clean case copy" >&2
    printf '  %s\n' "${existing[@]}" >&2
    exit 2
fi

decomposePar -no-libs -case "$case_dir" -region fluid
decomposePar -no-libs -case "$case_dir" -region solid
"$repo_root/scripts/validate-decomposition-constantHs.sh" "$case_dir"
