#!/bin/bash
set -euo pipefail

case_dir="$(cd "$(dirname "$0")" && pwd)"
find "$case_dir" -maxdepth 1 -mindepth 1 -type d \
    ! -name 0 \
    -regextype posix-extended -regex '.*/[0-9]+(\.[0-9]+)?(e[-+]?[0-9]+)?' \
    -exec rm -rf -- {} +
rm -rf -- "$case_dir/postProcessing"
rm -f -- "$case_dir"/log.* "$case_dir/run.pid"
rm -f -- "$case_dir/dealiiSolid/accepted-state.bin" \
    "$case_dir/dealiiSolid/trial-state.bin" \
    "$case_dir/dealiiSolid/robin-in.csv" \
    "$case_dir/dealiiSolid/robin-out.csv" \
    "$case_dir/dealiiSolid/robin-query.csv" \
    "$case_dir/dealiiSolid/solid-trial.vtu"
