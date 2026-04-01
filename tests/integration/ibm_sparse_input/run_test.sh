#!/bin/bash
# Test script for TEST_SPARSE_IJK - IBM file reading validation
#
# This script runs the sparse i,j,k coordinate file reading test
# for the Immersed Boundary Method implementation.
#
# The test validates that read_sparse_ijk() correctly distributes
# solid and fluid boundary point data across MPI ranks.

set -eu

# Configuration
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
UDALES_BUILD="${UDALES_BUILD:-${REPO_ROOT}/build/debug/u-dales}"
CASE_SOURCE="${CASE_SOURCE:-${REPO_ROOT}/tests/cases/101}"
NAMELIST="${NAMELIST:-namoptions.101}"
TMPDIR_PARENT="${TMPDIR_PARENT:-}"
MPIEXEC="${MPIEXEC:-mpiexec}"

# Check if executable exists
if [ ! -f "$UDALES_BUILD" ]; then
    echo "ERROR: u-dales executable not found at: $UDALES_BUILD"
    echo "Please build u-dales in debug mode first:"
    echo "  cd ${REPO_ROOT}/build/debug && cmake -DCMAKE_BUILD_TYPE=Debug ../.. && make"
    exit 1
fi

# Check if case source exists
if [ ! -d "$CASE_SOURCE" ]; then
    echo "ERROR: Case source not found: $CASE_SOURCE"
    exit 1
fi

# Check if namelist exists in source
if [ ! -f "$CASE_SOURCE/$NAMELIST" ]; then
    echo "ERROR: Namelist file not found: $CASE_SOURCE/$NAMELIST"
    exit 1
fi

if [ -z "${NPROCS:-}" ]; then
    NPROCX="$(awk -F= '/^[[:space:]]*nprocx[[:space:]]*=/ {gsub(/[[:space:]]/,"",$2); print $2; exit}' "$CASE_SOURCE/$NAMELIST")"
    NPROCY="$(awk -F= '/^[[:space:]]*nprocy[[:space:]]*=/ {gsub(/[[:space:]]/,"",$2); print $2; exit}' "$CASE_SOURCE/$NAMELIST")"
    if [ -z "$NPROCX" ] || [ -z "$NPROCY" ]; then
        echo "ERROR: Failed to determine nprocx/nprocy from $CASE_SOURCE/$NAMELIST"
        exit 1
    fi
    NPROCS="$((NPROCX * NPROCY))"
fi

if [ -n "$TMPDIR_PARENT" ]; then
    RUN_DIR="$(mktemp -d "${TMPDIR_PARENT%/}/sparse-ijk-XXXXXX")"
else
    RUN_DIR="$(mktemp -d)"
fi
trap 'rm -rf "$RUN_DIR"' EXIT

cp -r "$CASE_SOURCE"/. "$RUN_DIR"/
cp "$SCRIPT_DIR"/solid_*_[0-9]_[0-9].txt "$RUN_DIR"/
cp "$SCRIPT_DIR"/fluid_boundary_*_[0-9]_[0-9].txt "$RUN_DIR"/

cd "$RUN_DIR" || exit 1

echo "=========================================="
echo "Running TEST_SPARSE_IJK"
echo "=========================================="
echo "MPI processes: $NPROCS"
echo "Executable:    $UDALES_BUILD"
echo "Case source:   $(basename "$CASE_SOURCE")"
echo "Run directory: $RUN_DIR"
echo "Namelist:      $NAMELIST"
echo "=========================================="
echo ""

# Run the test
"$MPIEXEC" -n "$NPROCS" "$UDALES_BUILD" "$NAMELIST"

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
