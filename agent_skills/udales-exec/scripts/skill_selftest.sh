#!/usr/bin/env bash
set -euo pipefail

# Self-test for the udales-exec skill.
# Default: environment fingerprint + wrapper presence + Python tests.
# --full: run execution-side tests only (no builds).

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
MODE="${1:-}"

echo "udales-exec self-test"
echo "repo_root: ${ROOT_DIR}"

echo "== detect_env =="
"${ROOT_DIR}/skills/udales-exec/scripts/detect_env.sh"

echo "== required scripts =="
required=(
  "${ROOT_DIR}/tools/build_executable.sh"
  "${ROOT_DIR}/tools/build_preprocessing.sh"
  "${ROOT_DIR}/tools/hpc_execute.sh"
  "${ROOT_DIR}/tools/hpc_gather.sh"
  "${ROOT_DIR}/tests/run_tests.py"
)

missing=0
for path in "${required[@]}"; do
  if [ -e "$path" ]; then
    echo "ok: $path"
  else
    echo "missing: $path"
    missing=$((missing + 1))
  fi
done

echo "== optional integration scripts =="
optional=(
  "${ROOT_DIR}/tests/integration/mpi_operators/run_test.sh"
  "${ROOT_DIR}/tests/integration/ibm_sparse_input/run_test.sh"
  "${ROOT_DIR}/tests/integration/processor_boundaries/test_processor_boundaries.py"
)

for path in "${optional[@]}"; do
  if [ -e "$path" ]; then
    echo "ok: $path"
  else
    echo "missing: $path"
  fi
done

if [ "$missing" -ne 0 ]; then
  echo "self-test failed: missing ${missing} required script(s)"
  exit 1
fi

echo "== python tests =="
if command -v module >/dev/null 2>&1; then
  module load Python/3.9.6-GCCcore-11.2.0 >/dev/null 2>&1 || true
fi

VENV_PATH="${UD_VENV:-/rds/general/user/${USER}/home/udales/.venv}"
PYTHON_BIN="${UD_PYTHON:-python3}"

if [ -f "${VENV_PATH}/bin/activate" ]; then
  # shellcheck disable=SC1091
  source "${VENV_PATH}/bin/activate"
fi

if command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  export PYTHONPATH="${ROOT_DIR}/tools/python:${PYTHONPATH:-}"
  if ! "$PYTHON_BIN" -c "import udprep.directshortwave_f2py" >/dev/null 2>&1; then
    echo "missing: udprep.directshortwave_f2py (expected after preprocessing build)"
    echo "hint: run ./tools/build_preprocessing.sh <system> to build f2py extensions"
    exit 1
  fi
  "$PYTHON_BIN" -m unittest discover -s "${ROOT_DIR}/tools/python/tests" -p 'test_*.py'
else
  echo "missing: ${PYTHON_BIN}"
  exit 1
fi

echo "self-test passed"

if [ "$MODE" != "--full" ]; then
  exit 0
fi

echo "== full mode enabled (execution only, no builds) =="

if [ -f /etc/profile.d/modules.sh ]; then
  # shellcheck disable=SC1091
  source /etc/profile.d/modules.sh >/dev/null 2>&1 || true
fi

if command -v module >/dev/null 2>&1; then
  module load intel/2025a netCDF/4.9.2-iimpi-2023a netCDF-Fortran/4.6.1-iimpi-2023a \
    FFTW/3.3.9-intel-2021a CMake/3.29.3-GCCcore-13.3.0 git/2.45.1-GCCcore-13.3.0
fi

UD_BUILD="${UD_BUILD:-${ROOT_DIR}/build/debug/u-dales}"
if [ ! -x "$UD_BUILD" ]; then
  echo "missing: UD_BUILD executable at ${UD_BUILD}"
  echo "set UD_BUILD to a built solver binary before running --full"
  exit 1
fi

echo "== mpi operator integration =="
if [ -x "${ROOT_DIR}/tests/integration/mpi_operators/run_test.sh" ]; then
  UD_BUILD="$UD_BUILD" bash "${ROOT_DIR}/tests/integration/mpi_operators/run_test.sh"
fi

echo "full self-test complete"
