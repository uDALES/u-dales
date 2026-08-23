#!/bin/bash
# Test script for TEST_DRIVER_PLANES - driver inlet plane generation
#
# This script runs the host-branch test for moddriver::drivergen with idriver
# 2: that one call leaves all twelve driver inlet planes current, in each of
# the four time branches it can take, and that the m planes are written on
# rk3step 0 and 3 only.
#
# That last fact is why the test exists. The GPU driver inlet mirrors the
# twelve planes on the device, and initCUDA fills them there and then, because
# boundary_device only uploads on the stage drivergen runs on - so nothing
# uploads a complete set until the end of the first step. Allocating without
# filling left the first two stages of the first step reading whatever the
# allocation landed on, which cost case 452 a divergence of 4.33 on its first
# step against 2.7E-13 for the same case on the host.
#
# The store is built in the test rather than read from a driver file, so the
# expected values follow from the interpolation arithmetic rather than from a
# fixture that would have to be regenerated with the code it checks.
#
# The device side of the same invariant - that the mirror is live when
# initCUDA returns - is tests_cuda.f90::test_driver_plane_seed, which runs
# under UDALES_RUN_CUDA_SELFTEST on a Debug GPU build.

set -eu

# Configuration
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
# shellcheck disable=SC1091
source "${REPO_ROOT}/tests/integration/common/runtime_modules.sh"
load_udales_runtime_modules

UDALES_BUILD="${UDALES_BUILD:-${REPO_ROOT}/build/cpu/debug/u-dales}"
CASE_SOURCE="${CASE_SOURCE:-${REPO_ROOT}/tests/cases/101}"
NAMELIST_SOURCE="${NAMELIST_SOURCE:-${SCRIPT_DIR}/namoptions.1015}"
MPI_NAMELIST_SOURCE="${MPI_NAMELIST_SOURCE:-${SCRIPT_DIR}/namoptions.1015.mpi}"
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

# $1 label, $2 namelist source, $3 vertical grid override or empty.
run_one () {
    variant="$1"
    namelist_src="$2"
    prof_override="$3"

    if [ ! -f "$namelist_src" ]; then
        echo "ERROR: Namelist source not found: $namelist_src"
        return 1
    fi

    expnr="$(awk -F= '/^[[:space:]]*iexpnr[[:space:]]*=/ {gsub(/[[:space:]]/,"",$2); print $2; exit}' "$namelist_src")"
    if [ -z "$expnr" ]; then
        echo "ERROR: Failed to determine iexpnr from $namelist_src"
        return 1
    fi
    namelist="namoptions.${expnr}"

    nprocs="${NPROCS:-}"
    if [ -z "$nprocs" ]; then
        nprocx="$(awk -F= '/^[[:space:]]*nprocx[[:space:]]*=/ {gsub(/[[:space:]]/,"",$2); print $2; exit}' "$namelist_src")"
        nprocy="$(awk -F= '/^[[:space:]]*nprocy[[:space:]]*=/ {gsub(/[[:space:]]/,"",$2); print $2; exit}' "$namelist_src")"
        if [ -z "$nprocx" ] || [ -z "$nprocy" ]; then
            echo "ERROR: Failed to determine nprocx/nprocy from $namelist_src"
            return 1
        fi
        nprocs="$((nprocx * nprocy))"
    fi

    if [ -n "$TMPDIR_PARENT" ]; then
        RUN_DIR="$(mktemp -d "${TMPDIR_PARENT%/}/driver-planes-XXXXXX")"
    else
        RUN_DIR="$(mktemp -d)"
    fi

    cp -r "$CASE_SOURCE"/. "$RUN_DIR"/
    cp "$namelist_src" "$RUN_DIR/$namelist"
    if [ -n "$prof_override" ]; then
        if [ ! -f "$prof_override" ]; then
            echo "ERROR: override not found: $prof_override"
            rm -rf "$RUN_DIR"
            return 1
        fi
        cp "$prof_override" "$RUN_DIR/prof.inp.${expnr}"
    fi

    echo "=========================================="
    echo "Running TEST_DRIVER_PLANES (${variant})"
    echo "=========================================="
    echo "MPI processes: $nprocs"
    echo "Executable:    $UDALES_BUILD"
    echo "Case source:   $(basename "$CASE_SOURCE")"
    echo "Run directory: $RUN_DIR"
    echo "Namelist:      $namelist"
    echo "=========================================="
    echo ""

    ( cd "$RUN_DIR" && "$MPIEXEC" $MPI_LAUNCH_EXTRA_ARGS -n "$nprocs" "$UDALES_BUILD" "$namelist" )
    rc=$?

    rm -rf "$RUN_DIR"
    return $rc
}

EXIT_CODE=0

run_one "one rank" "$NAMELIST_SOURCE" "" || EXIT_CODE=$?
# The driver planes are indexed jb..je, which on one rank spans the whole
# domain - so a plane built over the global extent is indistinguishable from
# one built over the rank's own. The expected values carry zstart(2), so the
# two-rank pass is what makes that difference able to fail.
if [ "$EXIT_CODE" -eq 0 ]; then
    run_one "two ranks in y" "$MPI_NAMELIST_SOURCE" "" || EXIT_CODE=$?
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
