#!/bin/bash
# Prescribed facet flux on sloped facets (#409) using runmode 1008 on case 100,
# run on four decompositions whose totals must agree.

set -u

if ! command -v module >/dev/null 2>&1 && [ -f /etc/profile.d/modules.sh ]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
fi

# Load the Imperial (ICL) toolchain only where it exists. On other clusters
# (e.g. ECMWF) the caller provides the runtime via MPIEXEC and the module set
# used to build UDALES_BUILD; an unknown module name must not abort the test.
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

# launcher, extra arguments and the one-rank-per-CPU pre-flight check (never oversubscribe)
# shellcheck disable=SC1091
source "${REPO_ROOT}/tools/mpi_env.sh"

if [ ! -f "$UDALES_BUILD" ]; then
    echo "ERROR: u-dales executable not found at: $UDALES_BUILD"
    exit 1
fi

if [ ! -d "$CASE_SOURCE" ]; then
    echo "ERROR: Case source not found: $CASE_SOURCE"
    exit 1
fi

if [ -n "$TMPDIR_PARENT" ]; then
    ROOT_RUN_DIR="$(mktemp -d "${TMPDIR_PARENT%/}/facet-flux-XXXXXX")"
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
    echo "Running TEST_FACET_FLUX [$label]"
    echo "MPI processes: $np"
    echo "nprocx/nprocy: $npx/$npy"
    echo "Run directory: $run_dir"
    echo "=========================================="

    local run_rc=0
    (
        cd "$run_dir" || exit 1
        mpi_preflight "$np" || exit 1
        "$MPIEXEC" $MPI_LAUNCH_EXTRA_ARGS -n "$np" "$UDALES_BUILD" "$NAMELIST" > run.log 2>&1
    ) || run_rc=$?

    if [ "$run_rc" -ne 0 ]; then
        echo "FAIL: TEST_FACET_FLUX [$label] exited with code $run_rc"
        tail -n 120 "${run_dir}/run.log" || true
        failures=$((failures + 1))
        return
    fi

    if ! grep -q "ALL TESTS PASSED: tests_facet_flux" "${run_dir}/run.log"; then
        echo "FAIL: TEST_FACET_FLUX [$label]"
        tail -n 120 "${run_dir}/run.log" || true
        failures=$((failures + 1))
        return
    fi

    grep -h "FACETFLUX\|Prescribed facet fluxes" "${run_dir}/run.log"
    grep "^ *FACETFLUX" "${run_dir}/run.log" | sed 's/^[^:]*://' > "${ROOT_RUN_DIR}/${label}.totals"
}

run_mode serial 1 1 1 "${NAMELIST_DIR}/namoptions.1008.serial"
run_mode xsplit 2 1 2 "${NAMELIST_DIR}/namoptions.1008.xsplit"
run_mode ysplit 1 2 2 "${NAMELIST_DIR}/namoptions.1008.ysplit"
run_mode xysplit 2 2 4 "${NAMELIST_DIR}/namoptions.1008.xysplit"

# The totals must not depend on the decomposition (summation order only: relative 1e-10).
if [ "$failures" -eq 0 ]; then
    for label in ysplit xsplit xysplit; do
        if ! python3 - "${ROOT_RUN_DIR}/serial.totals" "${ROOT_RUN_DIR}/${label}.totals" <<'PYEOF'
import sys
a, b = [[float(x) for x in open(p).read().split()] for p in sys.argv[1:3]]
sys.exit(0 if len(a) == len(b) == 8 and all(abs(x - y) <= 1e-10 * max(abs(x), abs(y), 1e-300) for x, y in zip(a, b)) else 1)
PYEOF
        then
            echo "FAIL: TEST_FACET_FLUX totals differ between serial and $label"
            cat "${ROOT_RUN_DIR}/serial.totals" "${ROOT_RUN_DIR}/${label}.totals"
            failures=$((failures + 1))
        fi
    done
fi

echo "=========================================="
if [ "$failures" -eq 0 ]; then
    echo "All TEST_FACET_FLUX runs passed"
    echo "=========================================="
    exit 0
else
    echo "TEST_FACET_FLUX failures: $failures"
    echo "Preserving run directory: $ROOT_RUN_DIR"
    trap - EXIT
    echo "=========================================="
    exit 1
fi
