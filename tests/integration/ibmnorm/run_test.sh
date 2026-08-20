#!/bin/bash
# Test script for TEST_IBMNORM - host IBM normal-velocity validation
#
# Runs the runmode 1011 checks on the host routines modibm::ibmnorm is built
# from: solid, which pins the field inside the solid, and the two second-order
# advection corrections that remove the flux the scheme sent across solid
# faces.
#
# These are deliberately not a comparison against the device port. That
# comparison lives in tests_cuda.f90, it needs a GPU, and being symmetric it
# passes when both sides share a mistake. What runs here instead are properties
# the host routines must hold on their own:
#   - a solid velocity point ends at exactly zero, and nothing else moves;
#   - a solid scalar point walled in on all six sides takes the interior value,
#     while one with fluid neighbours takes their mean - so a constant field
#     stays constant across the boundary, and the interior value cannot reach
#     it;
#   - the liberal correction vanishes on a constant field for any velocity,
#     which is the flux it exists to cancel;
#   - both corrections vanish at zero velocity, are linear in the scalar, and
#     move only boundary points.
#
# The point lists are also checked to hold no cell twice: the device port runs
# one thread per point, so a repeat that a sequential loop tolerates would be
# two threads writing one cell.
#
# CPU only, and no GPU needed. Running this runmode against a GPU build fails
# with an explanatory message rather than passing vacuously.

set -eu

# Configuration
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
# shellcheck disable=SC1091
source "${REPO_ROOT}/tests/integration/common/runtime_modules.sh"
load_udales_runtime_modules

UDALES_BUILD="${UDALES_BUILD:-${REPO_ROOT}/build/cpu/debug/u-dales}"
CASE_SOURCE="${CASE_SOURCE:-${REPO_ROOT}/tests/cases/101}"
NAMELIST_SOURCE="${NAMELIST_SOURCE:-${SCRIPT_DIR}/namoptions.1011}"
NAMELIST_SOURCE_MPI="${NAMELIST_SOURCE_MPI:-${SCRIPT_DIR}/namoptions.1011.tworank}"
NAMELIST="${NAMELIST:-namoptions.101}"
TMPDIR_PARENT="${TMPDIR_PARENT:-}"
if [ -z "${MPIEXEC:-}" ] && command -v mpiifort >/dev/null 2>&1; then
    MPIEXEC="$(dirname "$(command -v mpiifort)")/mpiexec"
else
    MPIEXEC="${MPIEXEC:-mpiexec}"
fi
MPI_LAUNCH_EXTRA_ARGS="${MPI_LAUNCH_EXTRA_ARGS:-}"
MPI_VERSION_OUTPUT="$("$MPIEXEC" --version 2>/dev/null || true)"

if printf '%s\n' "$MPI_VERSION_OUTPUT" | grep -Eqi "Open MPI|OpenRTE"; then
    MPI_LAUNCH_EXTRA_ARGS="--oversubscribe ${MPI_LAUNCH_EXTRA_ARGS}"
    TMPDIR="${TMPDIR:-/tmp}"
    export TMPDIR
    export OMPI_MCA_prte_tmpdir_base="${OMPI_MCA_prte_tmpdir_base:-$TMPDIR}"
    export PRTE_MCA_prte_tmpdir_base="${PRTE_MCA_prte_tmpdir_base:-$TMPDIR}"
    if [ -z "$TMPDIR_PARENT" ]; then
        TMPDIR_PARENT="$TMPDIR"
    fi
fi

# Check if executable exists
if [ ! -f "$UDALES_BUILD" ]; then
    echo "ERROR: u-dales executable not found at: $UDALES_BUILD"
    echo "Please build u-dales in debug mode first:"
    echo "  ${REPO_ROOT}/tools/build_executable.sh common debug"
    exit 1
fi

# Check if case source exists
if [ ! -d "$CASE_SOURCE" ]; then
    echo "ERROR: Case source not found: $CASE_SOURCE"
    exit 1
fi

# $1 label, $2 namelist source
run_one () {
    variant="$1"
    namelist_src="$2"

    if [ ! -f "$namelist_src" ]; then
        echo "ERROR: Namelist source not found: $namelist_src"
        return 1
    fi

    nprocx="$(awk -F= '/^[[:space:]]*nprocx[[:space:]]*=/ {gsub(/[[:space:]]/,"",$2); print $2; exit}' "$namelist_src")"
    nprocy="$(awk -F= '/^[[:space:]]*nprocy[[:space:]]*=/ {gsub(/[[:space:]]/,"",$2); print $2; exit}' "$namelist_src")"
    if [ -z "$nprocx" ] || [ -z "$nprocy" ]; then
        echo "ERROR: Failed to determine nprocx/nprocy from $namelist_src"
        return 1
    fi
    nprocs="${NPROCS:-$((nprocx * nprocy))}"

    if [ -n "$TMPDIR_PARENT" ]; then
        RUN_DIR="$(mktemp -d "${TMPDIR_PARENT%/}/ibmnorm-XXXXXX")"
    else
        RUN_DIR="$(mktemp -d)"
    fi

    cp -r "$CASE_SOURCE"/. "$RUN_DIR"/
    cp "$namelist_src" "$RUN_DIR/$NAMELIST"

    echo "=========================================="
    echo "Running TEST_IBMNORM (${variant})"
    echo "=========================================="
    echo "MPI processes: $nprocs"
    echo "Executable:    $UDALES_BUILD"
    echo "Case source:   $(basename "$CASE_SOURCE")"
    echo "Run directory: $RUN_DIR"
    echo "Namelist:      $NAMELIST"
    echo "=========================================="
    echo ""

    ( cd "$RUN_DIR" && "$MPIEXEC" $MPI_LAUNCH_EXTRA_ARGS -n "$nprocs" "$UDALES_BUILD" "$NAMELIST" )
    rc=$?

    rm -rf "$RUN_DIR"
    return $rc
}

EXIT_CODE=0
run_one "single rank" "$NAMELIST_SOURCE" || EXIT_CODE=$?
if [ "$EXIT_CODE" -eq 0 ]; then
    run_one "two ranks, split domain" "$NAMELIST_SOURCE_MPI" || EXIT_CODE=$?
fi

echo ""
echo "=========================================="
if [ "$EXIT_CODE" -eq 0 ]; then
    echo "Test completed successfully"
else
    echo "Test failed with exit code: $EXIT_CODE"
fi
echo "=========================================="

exit $EXIT_CODE
