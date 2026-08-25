#!/bin/bash
set -euo pipefail

case_dir="$(cd "$(dirname "$0")" && pwd)"
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

cd "$case_dir"
exec solids4Foam -case "$case_dir" "$@"
