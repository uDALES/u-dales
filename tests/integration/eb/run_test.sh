#!/bin/bash
# Test script for TEST_EB - facet energy balance validation
#
# Runs the runmode 1012 checks on modEB: the energy balance itself, the flux
# time integration in intqH, and the longwave exchange in calclw.
#
# EB is the one routine in this part of the time loop that deliberately stays
# on the host, so unlike the other runmode tests here there is no device branch
# to compare against. What runs instead are properties the routine has to hold
# on its own:
#   - intqH keeps the facet fluxes on the third Runge-Kutta stage only, and
#     clears them on every stage. That is the contract letting modcuda's
#     updateFacFluxHost copy once a step instead of three times, so it is
#     checked against the routine rather than assumed;
#   - the energy balance fires only on the third stage, and only once its next
#     time has arrived, and marks the facet properties dirty exactly when it
#     rewrites them - the flag the GPU mirror now follows;
#   - a facet in equilibrium does not move: uniform layer temperatures with the
#     net surface flux balancing the emitted longwave is a fixed point, and
#     almost any slip in the matrix assembly, the inverse or the matmuls
#     destroys it;
#   - the surface energy balance closes against faclam, facTdash and the
#     forcing the test set, which is the physics rather than a second copy of
#     the algebra;
#   - the innermost layer is held fixed;
#   - the sparse and dense view-factor paths of calclw agree.
#
# Every check has a negative control that deliberately breaks it, and the
# controls are asserted to fire.
#
# The two-rank pass is what separates the reduced facet flux from the local
# one, and shows the time integral is kept on rank 0 alone.
#
# CPU only, and no GPU needed - though nothing in the routine is device-side,
# so it runs on either build.

set -eu

# Configuration
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
# shellcheck disable=SC1091
source "${REPO_ROOT}/tests/integration/common/runtime_modules.sh"
load_udales_runtime_modules

UDALES_BUILD="${UDALES_BUILD:-${REPO_ROOT}/build/cpu/debug/u-dales}"
CASE_SOURCE="${CASE_SOURCE:-${REPO_ROOT}/tests/cases/064}"
NAMELIST_SOURCE="${NAMELIST_SOURCE:-${SCRIPT_DIR}/namoptions.1012}"
NAMELIST_SOURCE_MPI="${NAMELIST_SOURCE_MPI:-${SCRIPT_DIR}/namoptions.1012.tworank}"
NAMELIST="${NAMELIST:-namoptions.064}"
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
        RUN_DIR="$(mktemp -d "${TMPDIR_PARENT%/}/eb-XXXXXX")"
    else
        RUN_DIR="$(mktemp -d)"
    fi

    cp -r "$CASE_SOURCE"/. "$RUN_DIR"/
    cp "$namelist_src" "$RUN_DIR/$NAMELIST"

    echo "=========================================="
    echo "Running TEST_EB (${variant})"
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
