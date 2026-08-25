#!/bin/bash
# Test script for TEST_TIMEDEP - the prescribed time-dependent forcings
#
# This script runs modtimedep::timedepsurf and modtimedep::timedepnudge against
# interpolation tables installed in the solver rather than read from a file, so
# that every entry can be made unique in level, in time and across the five
# fluxes and four profiles.
#
# Both routines are host code in either build, so unlike the tstep runmode this
# one is worth running against a GPU build too: it covers the bracket search and
# the freeze past the end of the table that timedep_device shares with timedep,
# and the surface fluxes, which the GPU path does not touch at all - the wall
# function kernels take those five scalars as implicit firstprivate at launch.
#
# What it does not cover is the OpenACC branch of the profile interpolation.
# That is tests_cuda.f90::test_timedep_nudge_device, which runs under
# UDALES_RUN_CUDA_SELFTEST on a Debug GPU build, and the timedep-nudging case in
# the GPU parity matrix.

set -eu

# Configuration
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
# shellcheck disable=SC1091
source "${REPO_ROOT}/tests/integration/common/runtime_modules.sh"
load_udales_runtime_modules

UDALES_BUILD="${UDALES_BUILD:-${REPO_ROOT}/build/cpu/debug/u-dales}"
CASE_SOURCE="${CASE_SOURCE:-${REPO_ROOT}/tests/cases/101}"
NAMELIST_SOURCE="${NAMELIST_SOURCE:-${SCRIPT_DIR}/namoptions.1018}"
MPI_NAMELIST_SOURCE="${MPI_NAMELIST_SOURCE:-${SCRIPT_DIR}/namoptions.1018.mpi}"
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
        RUN_DIR="$(mktemp -d "${TMPDIR_PARENT%/}/timedep-XXXXXX")"
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
    echo "Running TEST_TIMEDEP (${variant})"
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

# One rank first: the tables and the profiles are global, so a single rank is
# where the interpolation itself is easiest to read a failure out of.
run_one "single rank" "$NAMELIST_SOURCE" "" || EXIT_CODE=$?
# Then two, because the profiles are not decomposed and every rank interpolates
# its own copy from its own broadcast table. A rank that brackets differently
# from its neighbours nudges towards a different target on each side of the
# decomposition, and nothing downstream would say so.
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
