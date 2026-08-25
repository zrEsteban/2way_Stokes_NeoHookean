#!/bin/bash
set -euo pipefail

case_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root=${TWO_WAY_FSI_ROOT:-$(cd "$case_dir/../../../.." && pwd)}
state_dir="$case_dir/dealiiSolid"
export PDMS_CASE_DIR="$case_dir"

if ! command -v solids4Foam >/dev/null; then
    echo "ERROR: solids4Foam is not available; source env/activate-project.sh" >&2
    exit 2
fi
if [[ -z ${DEALII_PDMS_SOLID_BIN:-} ]]; then
    if ! DEALII_PDMS_SOLID_BIN=$(command -v dealiiPdmsSolid); then
        echo "ERROR: dealiiPdmsSolid is not available" >&2
        exit 2
    fi
    export DEALII_PDMS_SOLID_BIN
fi

if [[ -e "$state_dir/accepted-state.bin" || -e "$state_dir/trial-state.bin" ]]; then
    if ! grep -Eq 'startFrom[[:space:]]+latestTime' "$case_dir/system/controlDict"; then
        echo "ERROR: hay un estado deal.II previo pero OpenFOAM no usa latestTime" >&2
        exit 2
    fi
    echo "Continuando OpenFOAM/deal.II desde el último estado aceptado"
fi

parallel=false
for arg in "$@"; do
    if [[ "$arg" == "-parallel" ]]; then
        parallel=true
        [[ -x "$repo_root/scripts/validate-decomposition-constantHs.sh" ]] || {
            echo "ERROR: set TWO_WAY_FSI_ROOT to the project repository root" >&2
            exit 2
        }
        "$repo_root/scripts/validate-decomposition-constantHs.sh" "$case_dir"
        break
    fi
done

cd "$case_dir"
if [[ "$parallel" == true ]]; then
    command -v mpirun >/dev/null || {
        echo "ERROR: mpirun is unavailable" >&2
        exit 2
    }
    mpi_procs=${FSI_MPI_PROCS:-$(foamDictionary -no-libs \
        "$case_dir/system/fluid/decomposeParDict" \
        -entry numberOfSubdomains -value)}
    exec mpirun -np "$mpi_procs" solids4Foam -case "$case_dir" "$@"
fi
exec solids4Foam -case "$case_dir" "$@"
