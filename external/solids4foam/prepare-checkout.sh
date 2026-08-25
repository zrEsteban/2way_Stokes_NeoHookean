#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=manifest.env
source "$script_dir/manifest.env"

if [[ $# -ne 1 ]]; then
    echo "usage: $0 NEW_CHECKOUT_DIRECTORY" >&2
    exit 2
fi

destination=$1
if [[ -e "$destination" ]]; then
    echo "destination already exists: $destination" >&2
    exit 2
fi

clone_source=${SOLIDS4FOAM_MIRROR:-$SOLIDS4FOAM_URL}
git clone --no-checkout "$clone_source" "$destination"
git -C "$destination" checkout --detach "$SOLIDS4FOAM_COMMIT"

patch_file="$script_dir/$SOLIDS4FOAM_PATCH"
printf '%s  %s\n' "$SOLIDS4FOAM_PATCH_SHA256" "$patch_file" | sha256sum -c -
git -C "$destination" apply --check "$patch_file"
git -C "$destination" apply "$patch_file"

actual_diff_sha=$(
    git -C "$destination" diff --binary | sha256sum | awk '{print $1}'
)
if [[ "$actual_diff_sha" != "$SOLIDS4FOAM_RESULT_DIFF_SHA256" ]]; then
    echo "unexpected patched-tree diff SHA-256: $actual_diff_sha" >&2
    exit 1
fi

test "$(git -C "$destination" rev-parse HEAD)" = "$SOLIDS4FOAM_COMMIT"
test -z "$(git -C "$destination" ls-files --others --exclude-standard)"
echo "Prepared solids4foam $SOLIDS4FOAM_COMMIT with diff $actual_diff_sha"
