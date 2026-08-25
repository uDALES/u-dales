#!/bin/bash
# Test script for TEST_TSTEP - the adaptive time step limits
#
# This script runs the host-branch test for modtstep::tstep_limits, the fused
# reduction that produces the Courant and diffusion numbers tstep_update picks
# dt from, plus the geometry cache inittstep builds for it.
#
# Not modchecksim's pair of the same name. Those are diagnostics, printed and
# discarded; these choose the time step, and they are a different expression -
# absolute velocities, the cell-centred dxi and dyi rather than the half-level
# dxhi, and a division by dzh rather than a multiply by dzhi.
#
# The reduction is driven with a state whose answer follows from the grid
# rather than from re-running the loop, so a transcription error in the port
# cannot be reproduced by the expectation. What the passes below add on top of
# that is the grid and the decomposition the loop runs over.
#
# This covers the host branch only. The OpenACC kernel is covered by
# tests_cuda.f90::test_tstep_limits, which runs under UDALES_RUN_CUDA_SELFTEST
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
NAMELIST_SOURCE="${NAMELIST_SOURCE:-${SCRIPT_DIR}/namoptions.1017}"
MPI_NAMELIST_SOURCE="${MPI_NAMELIST_SOURCE:-${SCRIPT_DIR}/namoptions.1017.mpi}"
PROF_STRETCHED="${PROF_STRETCHED:-${SCRIPT_DIR}/prof.inp.stretched}"
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
        RUN_DIR="$(mktemp -d "${TMPDIR_PARENT%/}/tstep-XXXXXX")"
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
    echo "Running TEST_TSTEP (${variant})"
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

# The case ships a uniform vertical grid, where dzh and dzh2i hold the same
# value at every level. The test substitutes a dzh of its own before it drives
# the reduction, so the loop's vertical indexing is exercised on either grid -
# but the cache is checked against the grid inittstep actually saw, and on a
# uniform one dzh2i(k) + dx2i + dy2i is the same number at every k. The
# stretched pass is what makes that check able to fail, and what puts a real
# spread of magnitudes through the exact-equality comparison.
run_one "uniform vertical grid" "$NAMELIST_SOURCE" "" || EXIT_CODE=$?
if [ "$EXIT_CODE" -eq 0 ]; then
    run_one "stretched vertical grid" "$NAMELIST_SOURCE" "$PROF_STRETCHED" || EXIT_CODE=$?
fi
# On one rank jb..je spans the whole domain, so a loop written over the global
# extent is indistinguishable from one written over the rank's own. The
# two-rank pass separates them, and asserts the expectations are rank-local.
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
