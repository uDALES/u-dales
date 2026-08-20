#!/bin/bash
# Test script for TEST_MASSCORR - host flow-rate correction validation
#
# Runs the runmode 1010 checks on modforces::masscorr, which corrects the
# predicted velocity so the domain meets a prescribed flow rate. It has four
# branches - volume-averaged and outlet-plane, for u and for v - and each
# diagnoses a defect from reductions over the fluid cells, then adds that
# defect uniformly to the tendency.
#
# Checking the three scalars each branch forms internally would restate the
# code. What runs here instead measures, with its own reduction, the rate the
# routine claims to set:
#   - after the call the achieved rate equals uflowrate (or vflowrate), for
#     every branch;
#   - a second call is a no-op, since the target is already met. A defect that
#     is right once but inconsistent with the correction actually applied
#     fails this;
#   - the rate is the global one, so each rank gets a different field and all
#     ranks are checked against the same answer;
#   - u and v are independent, the correction is uniform over ib:ie, jb:je,
#     kb:ke and reaches nothing else, and linoutflow or all four switches off
#     is a genuine no-op.
#
# Two passes. On a single rank the local reduction and the global one are the
# same number, so a dropped MPI_ALLREDUCE is invisible; the second pass splits
# the domain in y and gives every rank a different slab of one global field.
#
# The v outlet branch is skipped at two ranks: its outlet is a row of constant
# j, which lives on one rank once the domain is split in y, while the branch
# all-reduces over every rank. That is a property of the scheme, not of this
# test, so the check runs only where the row is whole.
#
# This covers the host branch only. The OpenACC kernels are covered by
# tests_cuda.f90::test_masscorr, which runs under UDALES_RUN_CUDA_SELFTEST on a
# Debug GPU build. Running this runmode against a GPU build fails with an
# explanatory message rather than passing vacuously.

set -eu

# Configuration
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
# shellcheck disable=SC1091
source "${REPO_ROOT}/tests/integration/common/runtime_modules.sh"
load_udales_runtime_modules

UDALES_BUILD="${UDALES_BUILD:-${REPO_ROOT}/build/cpu/debug/u-dales}"
CASE_SOURCE="${CASE_SOURCE:-${REPO_ROOT}/tests/cases/101}"
NAMELIST_SOURCE="${NAMELIST_SOURCE:-${SCRIPT_DIR}/namoptions.1010}"
NAMELIST_SOURCE_MPI="${NAMELIST_SOURCE_MPI:-${SCRIPT_DIR}/namoptions.1010.tworank}"
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
        RUN_DIR="$(mktemp -d "${TMPDIR_PARENT%/}/masscorr-XXXXXX")"
    else
        RUN_DIR="$(mktemp -d)"
    fi

    cp -r "$CASE_SOURCE"/. "$RUN_DIR"/
    cp "$namelist_src" "$RUN_DIR/$NAMELIST"

    echo "=========================================="
    echo "Running TEST_MASSCORR (${variant})"
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
