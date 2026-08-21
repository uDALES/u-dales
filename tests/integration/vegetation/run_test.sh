#!/bin/bash
# Test script for TEST_VEGETATION - canopy drag and energy balance validation
#
# Runs the runmode 1013 checks on vegetation.f90: the momentum drag on the
# staggered face lists, the canopy energy balance, and the scalar deposition.
#
# vegetation_forcing was ported to the GPU by moving every loop onto the device
# and folding the time-independent part of the canopy energy balance into a
# startup cache, which also collapsed the two radiation modes into one routine.
# Neither change may move the physics, so the reference this test compares
# against is not the new code rearranged but the original expressions written
# out from the raw vegetation properties - both Beer-Lambert exponentials, the
# saturation-curve slope, r_a = 130*sqrt(lsize/sqrt(wind2)), and the decoupling
# factor in its original 1/(1 + 2*(gam/(s + 2*gam))*(rs/r_a)) form.
#
# On top of that it pins:
#   - the drag on all three face lists, and that nothing outside those lists
#     moves - checked over the whole declared extent, halos included, because
#     a list built or indexed one cell off still reproduces every value on it;
#   - qt = qtR + qtA, so the diagnostic split stays exhaustive;
#   - qe + qh = q_av_leaf, the absorbed radiation going somewhere;
#   - the mode dispatch, by running legacy and sveg over the same fields and
#     requiring each to match its own reference and to differ from the other;
#   - drag-only mode leaving heat and moisture alone;
#   - scalar deposition, per component.
#
# The shipped vegetation cases give every point the same properties, taken from
# one set of namelist values, so a per-point mix-up would reproduce the right
# answer everywhere. The test therefore imposes a bounded per-point profile on
# lad, rs, lsize, ud, dec and sveg, and asserts it is not flat before trusting
# any comparison. That is what separates m1-m7 and m13-m15 in the table below
# from a vacuous pass.
#
# The two-rank pass is what keeps the routine collective-safe: with the domain
# split in x the single tree lands on one rank, and the other has no vegetation
# points at all. Every rank still has to run init_vegetation the same number of
# times, so the test reduces the point count instead of returning early.
#
# CPU only. The device kernels are covered by
# tests_cuda.f90::test_vegetation_forcing and by the vegetation,
# vegetation-sveg and vegetation-deposition GPU parity cases.

set -eu

# Configuration
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
# shellcheck disable=SC1091
source "${REPO_ROOT}/tests/integration/common/runtime_modules.sh"
load_udales_runtime_modules

UDALES_BUILD="${UDALES_BUILD:-${REPO_ROOT}/build/cpu/debug/u-dales}"
CASE_SOURCE="${CASE_SOURCE:-${REPO_ROOT}/tests/cases/526}"
NAMELIST_SOURCE="${NAMELIST_SOURCE:-${SCRIPT_DIR}/namoptions.1013}"
NAMELIST_SOURCE_MPI="${NAMELIST_SOURCE_MPI:-${SCRIPT_DIR}/namoptions.1013.tworank}"
NAMELIST="${NAMELIST:-namoptions.526}"
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
        RUN_DIR="$(mktemp -d "${TMPDIR_PARENT%/}/veg-XXXXXX")"
    else
        RUN_DIR="$(mktemp -d)"
    fi

    cp -r "$CASE_SOURCE"/. "$RUN_DIR"/
    cp "$namelist_src" "$RUN_DIR/$NAMELIST"

    echo "=========================================="
    echo "Running TEST_VEGETATION (${variant})"
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
