#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
case_rel=cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS
tmp_root=$(mktemp -d /tmp/g1-robin-serialization.XXXXXX)
trap 'rm -rf "$tmp_root"' EXIT

command -v foamDictionary >/dev/null || {
    echo "ERROR: initialise OpenFOAM and the isolated project environment" >&2
    exit 2
}

copy_case()
{
    local destination=$1
    mkdir -p "$destination"
    git -C "$repo_root" archive HEAD:"$case_rel" | tar -x -C "$destination"
    "$destination/prepareCase.sh" >/dev/null
}

semantic_manifest()
{
    local field=$1 output=$2
    : > "$output"
    for entry in \
        dimensions \
        boundaryField.interface.type \
        boundaryField.interface.constantHs \
        boundaryField.interface.value \
        boundaryField.interface.prevPressure \
        boundaryField.interface.prevAcceleration \
        boundaryField.interface.coeff0 \
        boundaryField.interface.coeff1 \
        boundaryField.interface.rhs
    do
        printf '%s=' "$entry" >> "$output"
        foamDictionary -no-libs -precision 17 "$field" -entry "$entry" -value \
            >> "$output"
    done
}

copy_case "$tmp_root/plugin-a"
copy_case "$tmp_root/plugin-b"
copy_case "$tmp_root/no-libs"

for label in plugin-a plugin-b; do
    decomposePar -case "$tmp_root/$label" -region fluid -force >/dev/null
    "$repo_root/scripts/validate-decomposition-constantHs.sh" "$tmp_root/$label" \
        >/dev/null
    semantic_manifest "$tmp_root/$label/processor0/0/fluid/p" \
        "$tmp_root/$label.manifest"
done
cmp "$tmp_root/plugin-a.manifest" "$tmp_root/plugin-b.manifest"

decomposePar -no-libs -case "$tmp_root/no-libs" -region fluid -force >/dev/null
"$repo_root/scripts/validate-decomposition-constantHs.sh" "$tmp_root/no-libs" \
    >/dev/null

original_value=$(foamDictionary -no-libs -precision 17 \
    "$tmp_root/plugin-a/0/fluid/p" \
    -entry boundaryField.interface.constantHs -value)
written_value=$(foamDictionary -no-libs -precision 17 \
    "$tmp_root/plugin-a/processor0/0/fluid/p" \
    -entry boundaryField.interface.constantHs -value)
[[ "$written_value" == "$original_value" ]]

# Reuse G0's positive and negative validator controls.
"$repo_root/scripts/test-decomposition-preserves-constantHs.sh" >/dev/null

echo "Robin serialization PASS: plugin/no-libs, round-trip and idempotence"
