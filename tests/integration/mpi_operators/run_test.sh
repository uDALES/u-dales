#!/bin/bash
# Direct MPI operator and architecture validation using runmodes 1005/1006 on case 100.

set -eu

if ! command -v module >/dev/null 2>&1 && [ -f /etc/profile.d/modules.sh ]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
fi

if command -v module >/dev/null 2>&1; then
    module load intel/2025a netCDF/4.9.2-iimpi-2023a netCDF-Fortran/4.6.1-iimpi-2023a FFTW/3.3.9-intel-2021a CMake/3.29.3-GCCcore-13.3.0 git/2.45.1-GCCcore-13.3.0
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
UDALES_BUILD="${UDALES_BUILD:-${REPO_ROOT}/build/debug/u-dales}"
CASE_SOURCE="${CASE_SOURCE:-${REPO_ROOT}/tests/cases/100}"
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

if [ ! -f "$UDALES_BUILD" ]; then
    echo "ERROR: u-dales executable not found at: $UDALES_BUILD"
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
    ROOT_RUN_DIR="$(mktemp -d "${TMPDIR_PARENT%/}/mpi-operators-XXXXXX")"
else
    ROOT_RUN_DIR="$(mktemp -d)"
fi
trap 'rm -rf "$ROOT_RUN_DIR"' EXIT

failures=0

run_mode() {
    local test_name="$1"
    local runmode="$2"
    local success_pattern="$3"
    local label="$4"
    local npx="$5"
    local npy="$6"
    local np="$7"
    local run_dir="${ROOT_RUN_DIR}/${test_name}/${label}"

    mkdir -p "$run_dir"
    cp -r "$CASE_SOURCE"/. "$run_dir"/

    python3 - <<PY
from pathlib import Path
path = Path("${run_dir}/${NAMELIST}")
text = path.read_text()
if "runmode" not in text:
    text = text.replace("&RUN\n", "&RUN\n  runmode      = ${runmode}\n", 1)
else:
    import re
    text = re.sub(r"runmode\\s*=\\s*\\d+", "runmode      = ${runmode}", text, count=1)
text = text.replace("nprocx       = 16", "nprocx       = ${npx}")
text = text.replace("nprocy       = 16", "nprocy       = ${npy}")
path.write_text(text)
PY

    echo "=========================================="
    echo "Running ${test_name} [$label]"
    echo "MPI processes: $np"
    echo "nprocx/nprocy: $npx/$npy"
    echo "Run directory: $run_dir"
    echo "=========================================="

    (
        cd "$run_dir" || exit 1
        "$MPIEXEC" $MPI_LAUNCH_EXTRA_ARGS -n "$np" "$UDALES_BUILD" "$NAMELIST" > run.log 2>&1
    )

    if ! grep -q "$success_pattern" "${run_dir}/run.log"; then
        echo "FAIL: ${test_name} [$label]"
        tail -n 120 "${run_dir}/run.log" || true
        failures=$((failures + 1))
    fi
}

for test_name in TEST_MPI_OPERATORS TEST_MPI_ARCHITECTURE; do
    if [ "$test_name" = "TEST_MPI_OPERATORS" ]; then
        runmode=1005
        success_pattern="ALL TESTS PASSED: tests_mpi_fieldops"
    else
        runmode=1006
        success_pattern="ALL TESTS PASSED: tests_mpi_architecture"
    fi

    run_mode "$test_name" "$runmode" "$success_pattern" serial 1 1 1
    run_mode "$test_name" "$runmode" "$success_pattern" xsplit 2 1 2
    run_mode "$test_name" "$runmode" "$success_pattern" ysplit 1 2 2
    run_mode "$test_name" "$runmode" "$success_pattern" xysplit 2 2 4
done

echo "=========================================="
if [ "$failures" -eq 0 ]; then
    echo "All TEST_MPI_OPERATORS and TEST_MPI_ARCHITECTURE runs passed"
    echo "=========================================="
    exit 0
else
    echo "TEST_MPI_OPERATORS failures: $failures"
    echo "Preserving run directory: $ROOT_RUN_DIR"
    trap - EXIT
    echo "=========================================="
    exit 1
fi
