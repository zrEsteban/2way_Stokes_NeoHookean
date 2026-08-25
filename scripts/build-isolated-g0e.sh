#!/usr/bin/env bash
set -euo pipefail

if [[ -z ${WM_PROJECT_DIR:-} || -z ${WM_OPTIONS:-} ]]; then
    echo "source OpenFOAM v2512 before running this script" >&2
    exit 1
fi
if [[ -z ${PETSC_DIR:-} || -z ${PETSC_ARCH:-} ]]; then
    echo "PETSC_DIR and PETSC_ARCH are required" >&2
    exit 1
fi

project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
build_id=${PROJECT_BUILD_ID:-g0e-846a5f0}
prefix="$project_root/platforms/$build_id/$WM_OPTIONS"
s4f_root=${PROJECT_S4F_ROOT:-$project_root/external/solids4foam/source}

test -d "$s4f_root/.git"
test "$(git -C "$s4f_root" rev-parse HEAD)" = \
    4b254fa5260e0ae94640d7404089bde73907fc2d
test "$(git -C "$s4f_root" diff --binary | sha256sum | awk '{print $1}')" = \
    9099bbfcc247235ab9d736e5ed3e7587aceae5e1e4ad359e936db04a3ed207d0

export WM_PROJECT_USER_DIR="$project_root"
export FOAM_USER_APPBIN="$prefix/bin"
export FOAM_USER_LIBBIN="$prefix/lib"
export FOAM_MODULE_APPBIN="$FOAM_USER_APPBIN"
export FOAM_MODULE_LIBBIN="$FOAM_USER_LIBBIN"
export S4F_ROOT="$s4f_root"

filter_legacy_platform()
{
    local value=$1 entry filtered=
    local IFS=:
    for entry in $value; do
        case "$entry" in
            */OpenFOAM/ezamora-v2512/platforms/*) continue ;;
        esac
        filtered=${filtered:+$filtered:}$entry
    done
    printf '%s' "$filtered"
}

PATH=$(filter_legacy_platform "$PATH")
LD_LIBRARY_PATH=$(filter_legacy_platform "${LD_LIBRARY_PATH:-}")
export PATH="$FOAM_USER_APPBIN:$PATH"
export LD_LIBRARY_PATH="$FOAM_USER_LIBBIN:$PETSC_DIR/$PETSC_ARCH/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
unset -f filter_legacy_platform
mkdir -p "$FOAM_USER_APPBIN" "$FOAM_USER_LIBBIN" "$prefix/build-dealii"

# The checkout is disposable and fixed; build products never target the legacy
# OpenFOAM user tree.
S4F_NO_FILE_FIXES=1 "$s4f_root/Allwmake" -j "${BUILD_JOBS:-2}"

wclean libso "$project_root/src/robinRobinCoupling"
wmake libso "$project_root/src/robinRobinCoupling"
wclean libso "$project_root/src/fiveParameterMooneyRivlinElastic"
wmake libso "$project_root/src/fiveParameterMooneyRivlinElastic"

cmake -S "$project_root/src/dealiiPdmsSolid" -B "$prefix/build-dealii" \
    -DCMAKE_BUILD_TYPE=Release
cmake --build "$prefix/build-dealii" --parallel "${BUILD_JOBS:-2}"
install -m 755 "$prefix/build-dealii/dealiiPdmsSolid" \
    "$FOAM_USER_APPBIN/dealiiPdmsSolid"

sha256sum \
    "$FOAM_USER_APPBIN/solids4Foam" \
    "$FOAM_USER_APPBIN/dealiiPdmsSolid" \
    "$FOAM_USER_LIBBIN/libsolids4FoamModels.so" \
    "$FOAM_USER_LIBBIN/librobinRobinCoupling.so" \
    "$FOAM_USER_LIBBIN/libfiveParameterMooneyRivlinElastic.so"
