#!/bin/bash
# Test script for TEST_TREES_SPARSE_INPUT - sparse tree forcing validation
#
# This script runs the sparse tree input test, which compares block-based tree
# forcing against the sparse vegetation path on the committed 525 case.

set -eu

NPROCS="${NPROCS:-4}"
NPROCX="${NPROCX:-2}"
NPROCY="${NPROCY:-2}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
UDALES_BUILD="${UDALES_BUILD:-${REPO_ROOT}/build/debug/u-dales}"
CASE_SOURCE="${CASE_SOURCE:-${REPO_ROOT}/tests/cases/525}"
NAMELIST="${NAMELIST:-namoptions.525}"
TMPDIR_PARENT="${TMPDIR_PARENT:-}"

if [ ! -f "$UDALES_BUILD" ]; then
    echo "ERROR: u-dales executable not found at: $UDALES_BUILD"
    echo "Please build u-dales in debug mode first:"
    echo "  cd ${REPO_ROOT}/build/debug && cmake -DCMAKE_BUILD_TYPE=Debug ../.. && make"
    exit 1
fi

if [ ! -d "$CASE_SOURCE" ]; then
    echo "ERROR: Case source not found: $CASE_SOURCE"
    exit 1
fi

if [ ! -f "$CASE_SOURCE/$NAMELIST" ]; then
    echo "ERROR: Namelist file not found: $CASE_SOURCE/$NAMELIST"
    exit 1
fi

if [ -n "$TMPDIR_PARENT" ]; then
    RUN_DIR="$(mktemp -d "${TMPDIR_PARENT%/}/tree-sparse-XXXXXX")"
else
    RUN_DIR="$(mktemp -d)"
fi
trap 'rm -rf "$RUN_DIR"' EXIT

cp -r "$CASE_SOURCE"/. "$RUN_DIR"/
cd "$RUN_DIR" || exit 1

if [ "$NPROCS" -gt 4 ]; then
    echo "ERROR: Refusing to run tree_sparse_compare with more than 4 MPI ranks (requested $NPROCS)"
    exit 1
fi

if [ "$((NPROCX * NPROCY))" -ne "$NPROCS" ]; then
    echo "ERROR: NPROCS ($NPROCS) must equal NPROCX * NPROCY ($NPROCX * $NPROCY)"
    exit 1
fi

python3 - "$NAMELIST" "$NPROCX" "$NPROCY" <<'PY'
from pathlib import Path
import re
import sys

namelist = Path(sys.argv[1])
nprocx = sys.argv[2]
nprocy = sys.argv[3]
text = namelist.read_text(encoding="utf-8")
text = re.sub(r"(?m)^(\s*nprocx\s*=\s*)\d+(\s*)$", rf"\g<1>{nprocx}\2", text)
text = re.sub(r"(?m)^(\s*nprocy\s*=\s*)\d+(\s*)$", rf"\g<1>{nprocy}\2", text)
namelist.write_text(text, encoding="utf-8")
PY

echo "=========================================="
echo "Running TEST_TREES_SPARSE_COMPARE"
echo "=========================================="
echo "MPI processes: $NPROCS"
echo "Executable:    $UDALES_BUILD"
echo "Case source:   $(basename "$CASE_SOURCE")"
echo "Run directory: $RUN_DIR"
echo "Namelist:      $NAMELIST"
echo "=========================================="
echo ""

mpiexec -n "$NPROCS" "$UDALES_BUILD" "$NAMELIST"
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
