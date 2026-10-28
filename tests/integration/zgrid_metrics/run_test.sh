#!/bin/bash
# Vertical grid-metric validation using runmode 1006 on case 100.
#
# Runs TEST_ZGRID_METRICS on the uniform case-100 grid and on a stretched
# profile (prof_stretched.inp.100, generated with the udprep lstretchexpcheck
# method: dzlin=0.5 up to hlin=32 m, stretched to 100 m, dz ratio 3.8). The
# test asserts the identities of the solver's face recursion and the derived
# metric arrays, and that a stretched grid is rejected when combined with a
# Poisson configuration whose vertical eigenvalues assume equidistant spacing
# (here BCzp=2), which must make the run fail.

set -u

if ! command -v module >/dev/null 2>&1 && [ -f /etc/profile.d/modules.sh ]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
fi

if command -v module >/dev/null 2>&1 && module avail intel/2025a 2>&1 | grep -q intel; then
    module load intel/2025a netCDF/4.9.2-iimpi-2023a netCDF-Fortran/4.6.1-iimpi-2023a FFTW/3.3.9-intel-2021a CMake/3.29.3-GCCcore-13.3.0 git/2.45.1-GCCcore-13.3.0
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
UDALES_BUILD="${UDALES_BUILD:-${REPO_ROOT}/build/debug/u-dales}"
CASE_SOURCE="${CASE_SOURCE:-${REPO_ROOT}/tests/cases/100}"
NAMELIST_DIR="${NAMELIST_DIR:-${SCRIPT_DIR}}"
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

if [ -n "$TMPDIR_PARENT" ]; then
    ROOT_RUN_DIR="$(mktemp -d "${TMPDIR_PARENT%/}/zgrid-metrics-XXXXXX")"
else
    ROOT_RUN_DIR="$(mktemp -d)"
fi
trap 'rm -rf "$ROOT_RUN_DIR"' EXIT

failures=0

run_mode() {
    local label="$1"
    local np="$2"
    local namelist_source="$3"
    local profile_source="$4"   # empty = keep case profile
    local expect="$5"           # pass | fail
    local run_dir="${ROOT_RUN_DIR}/${label}"

    mkdir -p "$run_dir"
    cp -r "$CASE_SOURCE"/. "$run_dir"/
    if [ ! -f "$namelist_source" ]; then
        echo "ERROR: Namelist source not found: $namelist_source"
        exit 1
    fi
    cp "$namelist_source" "$run_dir/$NAMELIST"
    if [ -n "$profile_source" ]; then
        cp "$profile_source" "$run_dir/prof.inp.100"
    fi

    echo "=========================================="
    echo "Running TEST_ZGRID_METRICS [$label] (expect $expect)"
    echo "MPI processes: $np"
    echo "Run directory: $run_dir"
    echo "=========================================="

    local run_rc=0
    (
        cd "$run_dir" || exit 1
        "$MPIEXEC" $MPI_LAUNCH_EXTRA_ARGS -n "$np" "$UDALES_BUILD" "$NAMELIST" > run.log 2>&1
    ) || run_rc=$?

    if [ "$expect" = "pass" ]; then
        if [ "$run_rc" -ne 0 ] || ! grep -q "PASS: all vertical grid metric checks passed" "${run_dir}/run.log"; then
            echo "FAIL: TEST_ZGRID_METRICS [$label] (rc=$run_rc)"
            tail -n 60 "${run_dir}/run.log" || true
            failures=$((failures + 1))
        fi
    else
        if [ "$run_rc" -eq 0 ] || ! grep -q "stretched grid requires ipoiss=0" "${run_dir}/run.log"; then
            echo "FAIL: TEST_ZGRID_METRICS [$label] was expected to reject the configuration (rc=$run_rc)"
            tail -n 60 "${run_dir}/run.log" || true
            failures=$((failures + 1))
        fi
    fi
}

run_mode uniform           1 "${NAMELIST_DIR}/namoptions.1006.uniform"           ""                                        pass
run_mode stretched         1 "${NAMELIST_DIR}/namoptions.1006.stretched"         "${NAMELIST_DIR}/prof_stretched.inp.100" pass
run_mode stretched_xysplit 4 "${NAMELIST_DIR}/namoptions.1006.stretched_xysplit" "${NAMELIST_DIR}/prof_stretched.inp.100" pass
run_mode stretched_badpois 1 "${NAMELIST_DIR}/namoptions.1006.stretched_badpois" "${NAMELIST_DIR}/prof_stretched.inp.100" fail

echo "=========================================="
if [ "$failures" -eq 0 ]; then
    echo "ALL TESTS PASSED: zgrid_metrics"
    exit 0
else
    echo "FAILURES: $failures zgrid_metrics mode(s) failed"
    exit 1
fi
