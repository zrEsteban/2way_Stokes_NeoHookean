#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
work=${G2A_FILLET_QUALITY_WORK:-/tmp/g2a-fillet-mesh-quality}
mkdir -p "$work"

run_candidate() {
    local id=$1
    shift
    local output="$work/$id"
    test ! -e "$output" || { echo "candidate output already exists: $output" >&2; return 2; }
    local start=$SECONDS
    python3 "$root/scripts/generate_canonical_interface_geometry.py" --output "$output" "$@" >"$work/$id.manifest.log"
    checkMesh -no-libs -case "$output/convert" -allGeometry -allTopology -writeAllFields >"$work/$id.checkMesh.log" 2>&1
    printf '%s\n' "$((SECONDS-start))" >"$work/$id.runtime-seconds"
    python3 "$root/scripts/validate_canonical_interface_geometry.py" "$output" --json "$work/$id.geometry.json" >"$work/$id.geometry.log"
}

run_candidate original
run_candidate axial160 --fluid-axial-cells 160
run_candidate axial320 --fluid-axial-cells 320
run_candidate grading080 --fluid-side-progression 0.8
run_candidate grading100 --fluid-side-progression 1.0
run_candidate transverse2 --fluid-width-cells 48 --fluid-height-cells 20
run_candidate combo --fluid-axial-cells 320 --fluid-side-progression 0.8
run_candidate axial640 --fluid-axial-cells 640
run_candidate axial1280 --fluid-axial-cells 1280
run_candidate axial1920 --fluid-axial-cells 1920
