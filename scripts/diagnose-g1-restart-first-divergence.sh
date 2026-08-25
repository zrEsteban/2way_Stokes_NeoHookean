#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
case_rel=cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS
work=${G1R_WORK:-$(mktemp -d /tmp/g1r-first-divergence.XXXXXX)}
[[ ! -e "$work" ]] || { echo "ERROR: work directory exists: $work" >&2; exit 2; }
mkdir -p "$work"

new_case()
{
    local name=$1 end_time=$2 case_dir="$work/$1"
    mkdir -p "$case_dir"
    git -C "$repo_root" archive HEAD:"$case_rel" | tar -x -C "$case_dir"
    foamDictionary -no-libs "$case_dir/system/controlDict" -entry endTime -set "$end_time" >/dev/null
    foamDictionary -no-libs "$case_dir/system/controlDict" -entry writeInterval -set 1e-7 >/dev/null
    foamDictionary -no-libs "$case_dir/system/controlDict" -entry writePrecision -set 17 >/dev/null
    foamDictionary -no-libs "$case_dir/constant/solid/solidProperties" \
        -entry nonLinearGeometryUpdatedLagrangianCoeffs.restart -set true >/dev/null
    "$case_dir/prepareCase.sh" >"$work/$name.prepare.log" 2>&1
}

run_case()
{
    local name=$1 suffix=$2
    G1_RESTART_DIAGNOSTICS=1 "$work/$name/runCase.sh" \
        >"$work/$name.$suffix.log" 2>&1
}

new_case continuous 6e-7
new_case split 5e-7
run_case continuous full
run_case split first

find "$work/split/5e-07" -type f -print0 | sort -z | xargs -0 sha256sum \
    >"$work/split.before-read.sha256"
foamDictionary -no-libs "$work/split/system/controlDict" -entry startFrom -set latestTime >/dev/null
foamDictionary -no-libs "$work/split/system/controlDict" -entry endTime -set 6e-7 >/dev/null
run_case split restart

rg '^G1R ' "$work/continuous.full.log" >"$work/continuous.g1r"
rg '^G1R ' "$work/split.restart.log" >"$work/restart.g1r"
diff -u "$work/continuous.g1r" "$work/restart.g1r" >"$work/first-divergence.diff" || true
diff -qr "$work/continuous/5e-07" "$work/split/5e-07" \
    >"$work/snapshot.diff" || true

echo "G1-R evidence: $work"
cat "$work/first-divergence.diff"
