#!/bin/bash
# Test script for TEST_THERMODYNAMICS - host thermodynamics contract
#
# This script pins the facts about modthermodynamics::thermodynamics that the
# GPU port reproduces rather than corrects.
#
# The device self-test in tests_cuda.f90 runs the host routine and the device
# routine over one seed and requires them to agree, so it catches any
# disagreement between the two and nothing that both get wrong together. This
# test is that blind spot. The largest item in it: thermo writes ql0 one k
# level low, because its ql dummy is declared kb:ke+kh while its qt and thl
# dummies are declared kb-kh:ke+kh and the first call passes the whole of ql0.
# Sequence association then stores the saturation computed from level k at
# level k-kh and never writes ql0's top level at all. thermo_device declares
# its dummies the same way on purpose, so if this test starts failing because
# the host was corrected, the device declaration has to be corrected with it.
#
# The other checks are the surface override in calc_halflev, the eps1 clamp on
# dthvdz and its range, and the two halo behaviours of calthv - the saturated
# branch writes the interior only, the unsaturated one assigns thv0h whole.

set -eu

# Configuration
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
# shellcheck disable=SC1091
source "${REPO_ROOT}/tests/integration/common/runtime_modules.sh"
load_udales_runtime_modules

UDALES_BUILD="${UDALES_BUILD:-${REPO_ROOT}/build/cpu/debug/u-dales}"
CASE_SOURCE="${CASE_SOURCE:-${REPO_ROOT}/tests/cases/101}"
NAMELIST_SOURCE="${NAMELIST_SOURCE:-${SCRIPT_DIR}/namoptions.1016}"
MPI_NAMELIST_SOURCE="${MPI_NAMELIST_SOURCE:-${SCRIPT_DIR}/namoptions.1016.mpi}"
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
        RUN_DIR="$(mktemp -d "${TMPDIR_PARENT%/}/thermodynamics-XXXXXX")"
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
    echo "Running TEST_THERMODYNAMICS (${variant})"
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
# The slab averages inside diagfld reduce across ranks, and the fluid-cell
# masks the checks run over are per-rank. One rank owns the whole y extent and
# cannot tell a rank-local index from a global one; the two-rank pass can.
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
