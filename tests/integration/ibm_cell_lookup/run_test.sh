#!/bin/bash
# Test script for TEST_IBM_CELL_LOOKUP - IBM grid cell lookup validation
#
# This script runs the cell-lookup test for the Immersed Boundary Method
# reconstruction path.
#
# modibm::cell_index locates the cell containing a reconstruction point with
# count(p >= grid). That is only equivalent to the backward search it replaced,
# findloc(p >= grid, .true., 1, back=.true.), when the grid array is
# monotonically increasing and fully initialised. (findloc is not used because
# the NVHPC runtime aborts on logical arrays with
# "FINDLOC: unimplemented for data type".)
#
# The test validates that precondition, checks cell_index against an explicit
# backward search over a probe sweep, and re-derives the reconstruction indices
# initibm computed for the case.
#
# It runs against tests/cases/100, whose geometry has facet sections close
# enough to a wall to require full reconstruction. Grid-aligned cases never
# execute the lookup at all, so the test fails rather than passing vacuously if
# it finds no reconstruction sections.

set -eu

# Configuration
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
# shellcheck disable=SC1091
source "${REPO_ROOT}/tests/integration/common/runtime_modules.sh"
load_udales_runtime_modules

UDALES_BUILD="${UDALES_BUILD:-${REPO_ROOT}/build/cpu/debug/u-dales}"
CASE_SOURCE="${CASE_SOURCE:-${REPO_ROOT}/tests/cases/100}"
NAMELIST_SOURCE="${NAMELIST_SOURCE:-${SCRIPT_DIR}/namoptions.1006}"
NAMELIST="${NAMELIST:-namoptions.100}"
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
    RUN_DIR="$(mktemp -d "${TMPDIR_PARENT%/}/ibm-cell-lookup-XXXXXX")"
else
    RUN_DIR="$(mktemp -d)"
fi
trap 'rm -rf "$RUN_DIR"' EXIT

cp -r "$CASE_SOURCE"/. "$RUN_DIR"/
cp "$NAMELIST_SOURCE" "$RUN_DIR/$NAMELIST"

cd "$RUN_DIR" || exit 1

echo "=========================================="
echo "Running TEST_IBM_CELL_LOOKUP"
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
