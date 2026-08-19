#!/bin/bash
# Test script for TEST_IBM_WALLFUN - host IBM wall function validation
#
# Runs the runmode 1008 checks on the host implementations that modibm's GPU
# port mirrors: diffc_corr, wallfunmom and wallfunheat.
#
# These are deliberately not a comparison against the device port. That
# comparison lives in tests_cuda.f90, it needs a GPU, and being symmetric it
# passes when both sides share a mistake. What runs here instead are properties
# the host routines must hold on their own:
#   - no solid neighbours, or a constant field, means no correction at all;
#   - one isolated solid face has a value worked out by hand;
#   - the correction accumulates into the tendency rather than replacing it;
#   - only cells in the boundary point list are written;
#   - the momentum and heat reported per facet equal what left the field.
#
# CPU only, and no GPU needed.

set -eu

# Configuration
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
# shellcheck disable=SC1091
source "${REPO_ROOT}/tests/integration/common/runtime_modules.sh"
load_udales_runtime_modules

UDALES_BUILD="${UDALES_BUILD:-${REPO_ROOT}/build/cpu/debug/u-dales}"
CASE_SOURCE="${CASE_SOURCE:-${REPO_ROOT}/tests/cases/101}"
NAMELIST_SOURCE="${NAMELIST_SOURCE:-${SCRIPT_DIR}/namoptions.1008}"
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

PROF_STRETCHED="${PROF_STRETCHED:-${SCRIPT_DIR}/prof.inp.stretched}"
REC_CASE_SOURCE="${REC_CASE_SOURCE:-${REPO_ROOT}/tests/cases/100}"
REC_NAMELIST_SOURCE="${REC_NAMELIST_SOURCE:-${SCRIPT_DIR}/namoptions.1008.reconstruction}"
REC_TFACINIT="${REC_TFACINIT:-${SCRIPT_DIR}/Tfacinit.inp.reconstruction}"

# Run once per vertical grid. The case ships a uniform one, where dzf and dzh
# are equal at every level - so the cell-volume reciprocals dxdydzfi and
# dxdydzhi hold identical values and a swap between them, or either one built
# from the wrong spacing, is invisible. Every other case in the suite is uniform
# too, including the CPU-GPU parity matrix, so without the stretched pass
# nothing anywhere would catch it. The stretched pass is what makes the
# cell-volume checks, and the momentum and heat conservation cross-checks,
# able to fail.
# $1 label, $2 case directory, $3 namelist source, $4 vertical grid override
# or empty, $5 Tfacinit override or empty.
run_one () {
    variant="$1"
    case_dir="$2"
    namelist_src="$3"
    prof_override="$4"
    tfac_override="$5"

    expnr="${namelist_src##*.}"
    case "$expnr" in
        [0-9][0-9][0-9]) ;;
        *) expnr="$(awk -F= '/^[[:space:]]*iexpnr[[:space:]]*=/ {gsub(/[[:space:]]/,"",$2); print $2; exit}' "$namelist_src")" ;;
    esac
    namelist="namoptions.${expnr}"

    if [ -n "$TMPDIR_PARENT" ]; then
        RUN_DIR="$(mktemp -d "${TMPDIR_PARENT%/}/ibm_wallfun-XXXXXX")"
    else
        RUN_DIR="$(mktemp -d)"
    fi

    cp -r "$case_dir"/. "$RUN_DIR"/
    cp "$namelist_src" "$RUN_DIR/$namelist"
    for pair in "prof.inp.${expnr}|$prof_override" "Tfacinit.inp.${expnr}|$tfac_override"; do
        dest="${pair%%|*}"
        src="${pair#*|}"
        [ -n "$src" ] || continue
        if [ ! -f "$src" ]; then
            echo "ERROR: override not found: $src"
            rm -rf "$RUN_DIR"
            return 1
        fi
        cp "$src" "$RUN_DIR/$dest"
    done

    echo "=========================================="
    echo "Running TEST_IBM_WALLFUN (${variant})"
    echo "=========================================="
    echo "MPI processes: $NPROCS"
    echo "Executable:    $UDALES_BUILD"
    echo "Case source:   $(basename "$case_dir")"
    echo "Run directory: $RUN_DIR"
    echo "Namelist:      $namelist"
    echo "=========================================="
    echo ""

    ( cd "$RUN_DIR" && "$MPIEXEC" $MPI_LAUNCH_EXTRA_ARGS -n "$NPROCS" "$UDALES_BUILD" "$namelist" )
    rc=$?

    rm -rf "$RUN_DIR"
    return $rc
}

EXIT_CODE=0
run_one "uniform vertical grid" "$CASE_SOURCE" "$NAMELIST_SOURCE" "" "" || EXIT_CODE=$?
if [ "$EXIT_CODE" -eq 0 ]; then
    run_one "stretched vertical grid" "$CASE_SOURCE" "$NAMELIST_SOURCE" "$PROF_STRETCHED" "" || EXIT_CODE=$?
fi
if [ "$EXIT_CODE" -eq 0 ]; then
    run_one "reconstruction" "$REC_CASE_SOURCE" "$REC_NAMELIST_SOURCE" "" "$REC_TFACINIT" || EXIT_CODE=$?
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
