#!/usr/bin/env bash
set -euo pipefail

usage()
{
    echo "Usage: $0 CASE_DIR" >&2
    exit 2
}

[[ $# -eq 1 ]] || usage
case_dir=$(cd "$1" && pwd)
original="$case_dir/0/fluid/p"

command -v foamDictionary >/dev/null || {
    echo "ERROR: foamDictionary is unavailable; initialise OpenFOAM first" >&2
    exit 2
}
[[ -f "$original" ]] || {
    echo "ERROR: missing original field: $original" >&2
    exit 2
}

mapfile -t fsi_patches < <(
    foamDictionary -no-libs "$original" -entry boundaryField -keywords |
    while IFS= read -r patch; do
        patch_type=$(foamDictionary -no-libs "$original" \
            -entry "boundaryField.$patch.type" -value 2>/dev/null || true)
        if [[ "$patch_type" == "pdmsElasticWallPressure" ]]; then
            printf '%s\n' "$patch"
        fi
    done
)

(( ${#fsi_patches[@]} > 0 )) || {
    echo "ERROR: no pdmsElasticWallPressure patch in $original" >&2
    exit 1
}

original_dimensions=$(foamDictionary -no-libs -precision 17 "$original" \
    -entry dimensions -value)

shopt -s nullglob
processor_dirs=("$case_dir"/processor*)
(( ${#processor_dirs[@]} > 0 )) || {
    echo "ERROR: no processor directories under $case_dir" >&2
    exit 1
}

checked=0
for processor_dir in "${processor_dirs[@]}"; do
    [[ -d "$processor_dir" ]] || continue
    field="$processor_dir/0/fluid/p"
    [[ -f "$field" ]] || {
        echo "ERROR: missing decomposed field: $field" >&2
        exit 1
    }
    dimensions=$(foamDictionary -no-libs -precision 17 "$field" \
        -entry dimensions -value)
    [[ "$dimensions" == "$original_dimensions" ]] || {
        echo "ERROR: dimensions differ in $field: '$dimensions' != '$original_dimensions'" >&2
        exit 1
    }

    for patch in "${fsi_patches[@]}"; do
        patch_type=$(foamDictionary -no-libs "$field" \
            -entry "boundaryField.$patch.type" -value 2>/dev/null) || {
            echo "ERROR: FSI patch $patch missing in $field" >&2
            exit 1
        }
        [[ "$patch_type" == "pdmsElasticWallPressure" ]] || {
            echo "ERROR: FSI patch $patch has type '$patch_type' in $field" >&2
            exit 1
        }
        original_value=$(foamDictionary -no-libs -precision 17 "$original" \
            -entry "boundaryField.$patch.constantHs" -value 2>/dev/null) || {
            echo "ERROR: constantHs missing in original patch $patch" >&2
            exit 1
        }
        value=$(foamDictionary -no-libs -precision 17 "$field" \
            -entry "boundaryField.$patch.constantHs" -value 2>/dev/null) || {
            echo "ERROR: constantHs missing in $field patch $patch" >&2
            exit 1
        }
        [[ "$value" == "$original_value" ]] || {
            echo "ERROR: constantHs differs in $field patch $patch: '$value' != '$original_value'" >&2
            exit 1
        }
        ((checked += 1))
    done
done

expected=$(( ${#processor_dirs[@]} * ${#fsi_patches[@]} ))
(( checked == expected )) || {
    echo "ERROR: checked $checked constantHs entries, expected $expected" >&2
    exit 1
}

echo "constantHs validation PASS: $checked entries, dimensions $original_dimensions"
