#!/usr/bin/env bash

# Source this file after OpenFOAM v2512 has been initialised.
if [[ -z ${WM_PROJECT_DIR:-} || -z ${WM_OPTIONS:-} ]]; then
    echo "OpenFOAM is not initialised: WM_PROJECT_DIR/WM_OPTIONS missing" >&2
    return 1 2>/dev/null || exit 1
fi

project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export TWO_WAY_FSI_ROOT="$project_root"
build_id=${PROJECT_BUILD_ID:-g0e-846a5f0}
project_prefix="$project_root/platforms/$build_id/$WM_OPTIONS"
export S4F_ROOT=${PROJECT_S4F_ROOT:-$project_root/external/solids4foam/source}
export FOAM_USER_APPBIN="$project_prefix/bin"
export FOAM_USER_LIBBIN="$project_prefix/lib"
export WM_PROJECT_USER_DIR="$project_root"

if [[ ! -d "$S4F_ROOT/.git" ]]; then
    echo "fixed solids4foam checkout not found: $S4F_ROOT" >&2
    return 1 2>/dev/null || exit 1
fi
if [[ -z ${PETSC_DIR:-} || -z ${PETSC_ARCH:-} ]]; then
    echo "PETSC_DIR and PETSC_ARCH must identify the fixed PETSc build" >&2
    return 1 2>/dev/null || exit 1
fi

petsc_lib="$PETSC_DIR/$PETSC_ARCH/lib"

# OpenFOAM's user configuration may inject the legacy user platform. Remove
# those entries instead of leaving an accidental fallback behind the isolated
# prefix.
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
export LD_LIBRARY_PATH="$FOAM_USER_LIBBIN:$petsc_lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
unset -f filter_legacy_platform

expected_solver="$FOAM_USER_APPBIN/solids4Foam"
expected_models="$FOAM_USER_LIBBIN/libsolids4FoamModels.so"
expected_dealii="$FOAM_USER_APPBIN/dealiiPdmsSolid"
if [[ ! -x "$expected_solver" || ! -x "$expected_dealii" || ! -f "$expected_models" ]]; then
    echo "isolated solids4foam build is incomplete under $project_prefix" >&2
    return 1 2>/dev/null || exit 1
fi

resolved_solver=$(readlink -f "$(command -v solids4Foam)")
if [[ "$resolved_solver" != "$expected_solver" ]]; then
    echo "unexpected solids4Foam precedence: $resolved_solver" >&2
    return 1 2>/dev/null || exit 1
fi

for artifact in \
    "$expected_solver" \
    "$expected_dealii" \
    "$expected_models" \
    "$FOAM_USER_LIBBIN/librobinRobinCoupling.so" \
    "$FOAM_USER_LIBBIN/libfiveParameterMooneyRivlinElastic.so"
do
    if [[ ! -e "$artifact" ]]; then
        echo "required isolated artifact is missing: $artifact" >&2
        return 1 2>/dev/null || exit 1
    fi
    linkage=$(ldd "$artifact" 2>&1)
    if grep -q 'ezamora-v2512' <<<"$linkage"; then
        echo "forbidden runtime dependency in $artifact" >&2
        grep 'ezamora-v2512' <<<"$linkage" >&2
        return 1 2>/dev/null || exit 1
    fi
    if grep -q 'not found' <<<"$linkage"; then
        echo "unresolved runtime dependency in $artifact" >&2
        grep 'not found' <<<"$linkage" >&2
        return 1 2>/dev/null || exit 1
    fi
done

echo "Project root:       $project_root"
echo "OpenFOAM:           $WM_PROJECT_DIR ($WM_PROJECT_VERSION, $WM_OPTIONS)"
echo "solids4foam source: $S4F_ROOT @ $(git -C "$S4F_ROOT" rev-parse HEAD)"
echo "isolated prefix:    $project_prefix"
echo "PETSc:              $PETSC_DIR/$PETSC_ARCH"
echo "solids4Foam:        $resolved_solver"
