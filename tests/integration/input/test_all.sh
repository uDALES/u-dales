#!/usr/bin/env bash
# Simplified run_all script for tests/integration/input
# Assumptions:
# - This script is located in tests/integration/input
# - The two python test scripts are in the same directory
# - mpirun and other runtime dependencies are available

# Do not exit on first failure; capture each test's exit code and report at the end
set -uo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

PYTHON=python3
ROUNDTRIP_TEST=./test_roundtrip_input.py
SOURCE_TEST=./test_sourcecode.py

echo "Running round-trip test (default)..."
${PYTHON} ${ROUNDTRIP_TEST} --mode default
rc1=$?
if [ ${rc1} -eq 0 ]; then
  echo "test_roundtrip_input (default): PASSED"
else
  echo "test_roundtrip_input (default): FAILED"
fi

echo "Running round-trip test (random)..."
${PYTHON} ${ROUNDTRIP_TEST} --mode random
rc2=$?
if [ ${rc2} -eq 0 ]; then
  echo "test_roundtrip_input (random): PASSED"
else
  echo "test_roundtrip_input (random): FAILED"
fi

echo "Running sourcecode test..."
${PYTHON} ${SOURCE_TEST}
rc3=$?
if [ ${rc3} -eq 0 ]; then
  echo "test_sourcecode: PASSED"
else
  echo "test_sourcecode: FAILED"
fi

if [ ${rc1} -eq 0 ] && [ ${rc2} -eq 0 ] && [ ${rc3} -eq 0 ]; then
  echo
  echo "=== ALL TESTS PASSED ==="
  exit 0
else
  echo
  echo "=== SOME TESTS FAILED ==="
  echo "test_roundtrip_input (default):  $( [ ${rc1} -eq 0 ] && echo PASSED || echo FAILED )"
  echo "test_roundtrip_input (random):   $( [ ${rc2} -eq 0 ] && echo PASSED || echo FAILED )"
  echo "test_sourcecode:            $( [ ${rc3} -eq 0 ] && echo PASSED || echo FAILED )"
  exit 1
fi
