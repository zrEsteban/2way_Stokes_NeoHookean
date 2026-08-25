#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
case_rel=cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS
if [[ -n ${G1_RESTART_WORK:-} ]]; then
    work=$G1_RESTART_WORK
    [[ ! -e "$work" ]] || { echo "ERROR: work directory exists: $work" >&2; exit 2; }
    mkdir -p "$work"
else
    work=$(mktemp -d /tmp/g1-restart.XXXXXX)
    trap 'rm -rf "$work"' EXIT
fi

command -v foamDictionary >/dev/null || {
    echo "ERROR: initialise OpenFOAM and the project environment" >&2
    exit 2
}
: "${DEALII_PDMS_SOLID_BIN:=$(command -v dealiiPdmsSolid)}"
export DEALII_PDMS_SOLID_BIN

new_case()
{
    local name=$1 end_time=$2
    local case_dir="$work/$name"
    mkdir -p "$case_dir"
    git -C "$repo_root" archive HEAD:"$case_rel" | tar -x -C "$case_dir"
    foamDictionary -no-libs "$case_dir/system/controlDict" \
        -entry endTime -set "$end_time" >/dev/null
    foamDictionary -no-libs "$case_dir/system/controlDict" \
        -entry writeInterval -set 1e-7 >/dev/null
    foamDictionary -no-libs "$case_dir/system/controlDict" \
        -entry writePrecision -set 17 >/dev/null
    foamDictionary -no-libs "$case_dir/constant/fsiProperties" \
        -entry robinRobinCoeffs.stateAudit -set yes >/dev/null
    foamDictionary -no-libs "$case_dir/constant/solid/solidProperties" \
        -entry nonLinearGeometryUpdatedLagrangianCoeffs.restart \
        -set true >/dev/null
    "$case_dir/prepareCase.sh" >"$work/$name.prepare.log" 2>&1
}

run_case()
{
    local name=$1 parallel=$2 suffix=$3
    local case_dir="$work/$name"
    local args=()
    [[ "$parallel" == yes ]] && args=(-parallel)
    set +e
    "$case_dir/runCase.sh" "${args[@]}" >"$work/$name.$suffix.log" 2>&1
    local rc=$?
    set -e
    echo "$rc" >"$work/$name.$suffix.exit"
    (( rc == 0 )) || { echo "ERROR: $name/$suffix exit $rc" >&2; exit 1; }
}

restart_modes=${G1_RESTART_MODES:-"serial mpi"}
for mode in $restart_modes; do
    continuous="$mode-continuous"
    split="$mode-split"
    new_case "$continuous" 1e-6
    new_case "$split" 5e-7
    if [[ "$mode" == mpi ]]; then
        "$work/$continuous/decomposeCase.sh" >"$work/$continuous.decompose.log" 2>&1
        "$work/$split/decomposeCase.sh" >"$work/$split.decompose.log" 2>&1
        parallel=yes
    else
        parallel=no
    fi

    run_case "$continuous" "$parallel" full
    run_case "$split" "$parallel" first
    cp "$work/$split/postProcessing/fsiResiduals.dat" \
        "$work/$split.first-fsiResiduals.dat"
    if [[ "$parallel" == yes ]]; then
        while IFS= read -r processor_dir; do
            test -s "$processor_dir/4e-07/fluid/U"
            test -s "$processor_dir/5e-07/solid/DD"
        done < <(find "$work/$split" -maxdepth 1 -type d \
            -name 'processor*' -print | sort)
    else
        test -s "$work/$split/4e-07/fluid/U"
        test -s "$work/$split/5e-07/solid/DD"
    fi
    foamDictionary -no-libs "$work/$split/system/controlDict" \
        -entry startFrom -set latestTime >/dev/null
    foamDictionary -no-libs "$work/$split/system/controlDict" \
        -entry endTime -set 1e-6 >/dev/null
    run_case "$split" "$parallel" restart
    {
        head -1 "$work/$split.first-fsiResiduals.dat"
        tail -n +2 "$work/$split.first-fsiResiduals.dat"
        tail -n +2 "$work/$split/postProcessing/fsiResiduals.dat"
    } >"$work/$split.combined-fsiResiduals.dat"

    test "$(rg -c 'G1 state audit: accepted time step' \
        "$work/$continuous.full.log")" -eq 10
    test "$(( $(rg -c 'G1 state audit: accepted time step' "$work/$split.first.log") \
        + $(rg -c 'G1 state audit: accepted time step' "$work/$split.restart.log") ))" -eq 10
    rg -q 'time integration=BDF2' "$work/$split.restart.log"
    rg -q '^End$' "$work/$continuous.full.log"
    rg -q '^End$' "$work/$split.restart.log"

    cmp "$work/$continuous/dealiiSolid/robin-out.csv" \
        "$work/$split/dealiiSolid/robin-out.csv"
    cmp "$work/$continuous/postProcessing/fsiResiduals.dat" \
        "$work/$split.combined-fsiResiduals.dat"

    sha256sum "$work/$continuous/postProcessing/fsiResiduals.dat" \
        "$work/$continuous/dealiiSolid/robin-out.csv"
done

if [[ " $restart_modes " == *" serial "* && " $restart_modes " == *" mpi "* ]]; then
for name in serial-continuous mpi-continuous; do
    file="$work/$name/dealiiSolid/robin-out.csv"
    { head -1 "$file"; tail -n +2 "$file" | sort -t, -k1,1g -k2,2g -k3,3g; } \
        >"$work/$name.sorted.csv"
done

paste -d, "$work/serial-continuous.sorted.csv" "$work/mpi-continuous.sorted.csv" |
awk -F, '
  NR>1 {
    for (i=4;i<=15;i++) {
      d=$(i)-$(i+15); if (d<0) d=-d;
      a=$(i); if (a<0) a=-a;
      if (d>maxDiff[i]) maxDiff[i]=d;
      if (a>maxRef[i]) maxRef[i]=a;
    }
  }
  END {
    for (i=4;i<=15;i++)
      if (maxDiff[i]/(maxRef[i]+1e-300)>5e-5) exit 1
  }
'
fi

if rg -n -i '(^|[^[:alpha:]])(nan|[-+]?inf(inity)?)([^[:alpha:]]|$)' \
    "$work"/*.log; then
    echo "ERROR: NaN/Inf found" >&2
    exit 1
fi

echo "G1 restart PASS: continuous/split serial+MPI and serial/MPI rtol=5e-5"
