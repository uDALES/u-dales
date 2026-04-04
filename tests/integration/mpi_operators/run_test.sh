#!/bin/bash
# Direct MPI operator validation using runmode 1005 on case 100.

set -u

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
    local label="$1"
    local npx="$2"
    local npy="$3"
    local np="$4"
    local run_dir="${ROOT_RUN_DIR}/${label}"

    mkdir -p "$run_dir"
    cp -r "$CASE_SOURCE"/. "$run_dir"/

    python3 - <<PY
from pathlib import Path
import re

path = Path("${run_dir}/${NAMELIST}")
text = path.read_text()
if "runmode" not in text:
    text = text.replace("&RUN\n", "&RUN\n  runmode      = 1005\n", 1)
else:
    text = re.sub(r"runmode\\s*=\\s*\\d+", "runmode      = 1005", text, count=1)
text = re.sub(r"(?m)^(\\s*nprocx\\s*=\\s*)\\d+(\\s*)$", r"\\g<1>${npx}\\2", text)
text = re.sub(r"(?m)^(\\s*nprocy\\s*=\\s*)\\d+(\\s*)$", r"\\g<1>${npy}\\2", text)
path.write_text(text)
PY

    echo "=========================================="
    echo "Running TEST_MPI_OPERATORS [$label]"
    echo "MPI processes: $np"
    echo "nprocx/nprocy: $npx/$npy"
    echo "Run directory: $run_dir"
    echo "=========================================="

    local run_rc=0
    (
        cd "$run_dir" || exit 1
        "$MPIEXEC" $MPI_LAUNCH_EXTRA_ARGS -n "$np" "$UDALES_BUILD" "$NAMELIST" > run.log 2>&1
    ) || run_rc=$?

    if [ "$run_rc" -ne 0 ]; then
        echo "FAIL: TEST_MPI_OPERATORS [$label] exited with code $run_rc"
        tail -n 120 "${run_dir}/run.log" || true
        failures=$((failures + 1))
        return
    fi

    if ! grep -q "ALL TESTS PASSED: tests_mpi_operators" "${run_dir}/run.log"; then
        echo "FAIL: TEST_MPI_OPERATORS [$label]"
        tail -n 120 "${run_dir}/run.log" || true
        failures=$((failures + 1))
    fi
}

run_mode serial 1 1 1
run_mode xsplit 2 1 2
run_mode ysplit 1 2 2
run_mode xysplit 2 2 4

echo "=========================================="
if [ "$failures" -eq 0 ]; then
    echo "All TEST_MPI_OPERATORS runs passed"
    echo "=========================================="
    exit 0
else
    echo "TEST_MPI_OPERATORS failures: $failures"
    echo "Preserving run directory: $ROOT_RUN_DIR"
    trap - EXIT
    echo "=========================================="
    exit 1
fi
