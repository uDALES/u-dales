#!/bin/bash
# SGS statistics validation using runmode 1006 on the stretched-grid case 300.

set -u

if ! command -v module >/dev/null 2>&1 && [ -f /etc/profile.d/modules.sh ]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
fi

# NOTE: this default matches the other solver suites, but it does NOT match the
# stack tools/build_executable.sh actually builds with on CX3
# (intel/2021a + netCDF-Fortran/4.5.3). Running a netCDF-writing case with the
# 4.6.1 runtime against a 4.5.3-linked binary segfaults inside nf90_def_var.
# This test never writes NetCDF so it is unaffected; override
# UDALES_RUNTIME_MODULES if you need the build-matching stack.
UDALES_RUNTIME_MODULES="${UDALES_RUNTIME_MODULES:-intel/2025a netCDF/4.9.2-iimpi-2023a netCDF-Fortran/4.6.1-iimpi-2023a FFTW/3.3.9-intel-2021a CMake/3.29.3-GCCcore-13.3.0 git/2.45.1-GCCcore-13.3.0}"

if command -v module >/dev/null 2>&1; then
    module load $UDALES_RUNTIME_MODULES
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
UDALES_BUILD="${UDALES_BUILD:-${REPO_ROOT}/build/debug/u-dales}"
CASE_SOURCE="${CASE_SOURCE:-${REPO_ROOT}/tests/cases/300}"
NAMELIST_DIR="${NAMELIST_DIR:-${SCRIPT_DIR}}"
NAMELIST="${NAMELIST:-namoptions.300}"
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

if [ -n "$TMPDIR_PARENT" ]; then
    ROOT_RUN_DIR="$(mktemp -d "${TMPDIR_PARENT%/}/sgs-statistics-XXXXXX")"
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
    local namelist_source="$5"
    local run_dir="${ROOT_RUN_DIR}/${label}"

    mkdir -p "$run_dir"
    cp -r "$CASE_SOURCE"/. "$run_dir"/
    if [ ! -f "$namelist_source" ]; then
        echo "ERROR: Namelist source not found: $namelist_source"
        exit 1
    fi
    cp "$namelist_source" "$run_dir/$NAMELIST"

    echo "=========================================="
    echo "Running TEST_SGS_STATS [$label]"
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
        echo "FAIL: TEST_SGS_STATS [$label] exited with code $run_rc"
        tail -n 120 "${run_dir}/run.log" || true
        failures=$((failures + 1))
        return
    fi

    if ! grep -q "ALL TESTS PASSED: tests_sgs_statistics" "${run_dir}/run.log"; then
        echo "FAIL: TEST_SGS_STATS [$label]"
        tail -n 120 "${run_dir}/run.log" || true
        failures=$((failures + 1))
    fi
}

run_mode serial 1 1 1 "${NAMELIST_DIR}/namoptions.1006.serial"
run_mode xsplit 2 1 2 "${NAMELIST_DIR}/namoptions.1006.xsplit"
run_mode ysplit 1 2 2 "${NAMELIST_DIR}/namoptions.1006.ysplit"
run_mode xysplit 2 2 4 "${NAMELIST_DIR}/namoptions.1006.xysplit"

echo "=========================================="
if [ "$failures" -eq 0 ]; then
    echo "All TEST_SGS_STATS runs passed"
    echo "=========================================="
    exit 0
else
    echo "TEST_SGS_STATS failures: $failures"
    echo "Preserving run directory: $ROOT_RUN_DIR"
    trap - EXIT
    echo "=========================================="
    exit 1
fi
