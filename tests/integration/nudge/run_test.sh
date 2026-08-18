#!/bin/bash
# Test script for TEST_NUDGE - profile-nudging tendency validation
#
# This script runs the host-branch test for modforces::nudge.
#
# nudge relaxes the upper part of the domain towards prescribed profiles at rate
# 1/tnudge, starting at level kb+nnudge. The test pins two properties that the
# CPU and GPU branches have to agree on:
#   - the relaxation starts at kb+nnudge, so levels below it stay untouched;
#   - it is applied over each array's full declared extent, halos included,
#     because the CPU branch uses whole-array assignment. An explicit-loop GPU
#     port that covered only ib..ie would diverge, so the halo cells are
#     asserted separately.
#
# This covers the host branch only. The OpenACC kernels are covered by
# tests_cuda.f90::test_nudge_profiles, which runs under UDALES_RUN_CUDA_SELFTEST
# on a Debug GPU build. Running this runmode against a GPU build fails with an
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
NAMELIST_SOURCE="${NAMELIST_SOURCE:-${SCRIPT_DIR}/namoptions.1007}"
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

# Check if namelist exists
if [ ! -f "$NAMELIST_SOURCE" ]; then
    echo "ERROR: Namelist source not found: $NAMELIST_SOURCE"
    exit 1
fi

if [ -z "${NPROCS:-}" ]; then
    NPROCX="$(awk -F= '/^[[:space:]]*nprocx[[:space:]]*=/ {gsub(/[[:space:]]/,"",$2); print $2; exit}' "$NAMELIST_SOURCE")"
    NPROCY="$(awk -F= '/^[[:space:]]*nprocy[[:space:]]*=/ {gsub(/[[:space:]]/,"",$2); print $2; exit}' "$NAMELIST_SOURCE")"
    if [ -z "$NPROCX" ] || [ -z "$NPROCY" ]; then
        echo "ERROR: Failed to determine nprocx/nprocy from $NAMELIST_SOURCE"
        exit 1
    fi
    NPROCS="$((NPROCX * NPROCY))"
fi

if [ -n "$TMPDIR_PARENT" ]; then
    RUN_DIR="$(mktemp -d "${TMPDIR_PARENT%/}/nudge-XXXXXX")"
else
    RUN_DIR="$(mktemp -d)"
fi
trap 'rm -rf "$RUN_DIR"' EXIT

cp -r "$CASE_SOURCE"/. "$RUN_DIR"/
cp "$NAMELIST_SOURCE" "$RUN_DIR/$NAMELIST"

cd "$RUN_DIR" || exit 1

echo "=========================================="
echo "Running TEST_NUDGE"
echo "=========================================="
echo "MPI processes: $NPROCS"
echo "Executable:    $UDALES_BUILD"
echo "Case source:   $(basename "$CASE_SOURCE")"
echo "Run directory: $RUN_DIR"
echo "Namelist:      $NAMELIST"
echo "=========================================="
echo ""

# Run the test
"$MPIEXEC" $MPI_LAUNCH_EXTRA_ARGS -n "$NPROCS" "$UDALES_BUILD" "$NAMELIST"

# Capture exit code
EXIT_CODE=$?

echo ""
echo "=========================================="
if [ "$EXIT_CODE" -eq 0 ]; then
    echo "Test completed successfully"
else
    echo "Test failed with exit code: $EXIT_CODE"
fi
echo "=========================================="

exit $EXIT_CODE
