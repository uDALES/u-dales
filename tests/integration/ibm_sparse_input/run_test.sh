#!/bin/bash
# Test script for TEST_SPARSE_IJK - IBM file reading validation
#
# This script runs the sparse i,j,k coordinate file reading test
# for the Immersed Boundary Method implementation.
#
# The test validates that read_sparse_ijk() correctly distributes
# solid and fluid boundary point data across MPI ranks.

set -eu

if ! command -v module >/dev/null 2>&1 && [ -f /etc/profile.d/modules.sh ]; then
    # Ensure the cluster module command is available in non-login shells.
    # The supported solver build uses the Intel/iimpi stack below.
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
fi

if command -v module >/dev/null 2>&1; then
    module load intel/2025a netCDF/4.9.2-iimpi-2023a netCDF-Fortran/4.6.1-iimpi-2023a FFTW/3.3.9-intel-2021a CMake/3.29.3-GCCcore-13.3.0 git/2.45.1-GCCcore-13.3.0
fi

# Configuration
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
UDALES_BUILD="${UDALES_BUILD:-${REPO_ROOT}/build/debug/u-dales}"
CASE_SOURCE="${CASE_SOURCE:-${REPO_ROOT}/tests/cases/101}"
NAMELIST_SOURCE="${NAMELIST_SOURCE:-${SCRIPT_DIR}/namoptions.1004}"
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
    echo "  cd ${REPO_ROOT}/build/debug && cmake -DCMAKE_BUILD_TYPE=Debug ../.. && make"
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
    RUN_DIR="$(mktemp -d "${TMPDIR_PARENT%/}/sparse-ijk-XXXXXX")"
else
    RUN_DIR="$(mktemp -d)"
fi
trap 'rm -rf "$RUN_DIR"' EXIT

cp -r "$CASE_SOURCE"/. "$RUN_DIR"/
cp "$NAMELIST_SOURCE" "$RUN_DIR/$NAMELIST"
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
