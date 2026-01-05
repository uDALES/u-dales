#!/bin/bash
# Test script for TEST_SPARSE_IJK - IBM file reading validation
#
# This script runs the sparse i,j,k coordinate file reading test
# for the Immersed Boundary Method implementation.
#
# The test validates that read_sparse_ijk() correctly distributes
# solid and fluid boundary point data across MPI ranks.

# Configuration
NPROCS=4
UDALES_BUILD="../build/debug/u-dales"
CASE_DIR="../examples/101"
NAMELIST="namoptions.101"

# Check if executable exists
if [ ! -f "$UDALES_BUILD" ]; then
    echo "ERROR: u-dales executable not found at: $UDALES_BUILD"
    echo "Please build u-dales in debug mode first:"
    echo "  cd ../build/debug && cmake -DCMAKE_BUILD_TYPE=Debug ../.. && make"
    exit 1
fi

# Check if case directory exists
if [ ! -d "$CASE_DIR" ]; then
    echo "ERROR: Case directory not found: $CASE_DIR"
    exit 1
fi

# Check if namelist exists
if [ ! -f "$CASE_DIR/$NAMELIST" ]; then
    echo "ERROR: Namelist file not found: $CASE_DIR/$NAMELIST"
    exit 1
fi

# Change to case directory
cd "$CASE_DIR" || exit 1

echo "=========================================="
echo "Running TEST_SPARSE_IJK"
echo "=========================================="
echo "MPI processes: $NPROCS"
echo "Executable:    $UDALES_BUILD"
echo "Case:          $(basename $CASE_DIR)"
echo "Namelist:      $NAMELIST"
echo "=========================================="
echo ""

# Run the test
mpiexec -n $NPROCS "$UDALES_BUILD" "$NAMELIST"

# Capture exit code
EXIT_CODE=$?

echo ""
echo "=========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "Test completed successfully"
else
    echo "Test failed with exit code: $EXIT_CODE"
fi
echo "=========================================="

exit $EXIT_CODE
