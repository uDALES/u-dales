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
    local family="$6"
    local label="$1"
    local npx="$2"
    local npy="$3"
    local np="$4"
    local namelist_source="$5"
    local run_dir="${ROOT_RUN_DIR}/${family}-${label}"

    mkdir -p "$run_dir"
    cp -r "$CASE_SOURCE"/. "$run_dir"/
    if [ ! -f "$namelist_source" ]; then
        echo "ERROR: Namelist source not found: $namelist_source"
        exit 1
    fi
    cp "$namelist_source" "$run_dir/$NAMELIST"
    if [ "$family" = mixed ]; then # iwalltemp = 2 reads the initial facet temperatures
        { echo "# initial facet temperatures"; awk 'NR > 1 { print "288.000000" }' "$run_dir/facets.inp.100"; } > "$run_dir/Tfacinit.inp.100"
    fi
    if [ "$family" = pertype ]; then
        cp "${SCRIPT_DIR}/factypeflux.inp.100" "$run_dir/"
    fi

    echo "=========================================="
    echo "Running TEST_FACET_FLUX [$family $label]"
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
        echo "FAIL: TEST_FACET_FLUX [$family $label] exited with code $run_rc"
        tail -n 120 "${run_dir}/run.log" || true
        failures=$((failures + 1))
        return
    fi

    if ! grep -q "ALL TESTS PASSED: tests_facet_flux" "${run_dir}/run.log"; then
        echo "FAIL: TEST_FACET_FLUX [$family $label]"
        tail -n 120 "${run_dir}/run.log" || true
        failures=$((failures + 1))
        return
    fi

    local nlines
    nlines=$(grep -c "^Prescribed " "${run_dir}/run.log" || true)
    if [ "$family" = direct ]; then
        if ! grep -q "^Prescribed heat flux (iwalltemp = 1): 287 of 1102 facets" "${run_dir}/run.log" || \
           ! grep -q "^Prescribed moisture flux (iwallmoist = 1): 287 of 1102 facets" "${run_dir}/run.log"; then
            echo "FAIL: TEST_FACET_FLUX [$family $label] start-up lines missing"
            failures=$((failures + 1))
        fi
    elif [ "$nlines" -ne 0 ]; then
        echo "FAIL: TEST_FACET_FLUX [$family $label] unexpected start-up line"
        failures=$((failures + 1))
    fi
    if [ "$family" = pertype ] && ! grep -q "^Read factypeflux.inp.100: prescribed flux for 1 facet types (1102 of 1102 facets)" "${run_dir}/run.log"; then
        echo "FAIL: TEST_FACET_FLUX [$family $label] factypeflux.inp start-up line missing"
        failures=$((failures + 1))
    fi
    grep -h "FACETFLUX\|^Prescribed \|^Read factypeflux" "${run_dir}/run.log"
    grep "^ *FACETFLUX" "${run_dir}/run.log" | sed 's/^[^:]*://' > "${ROOT_RUN_DIR}/${family}-${label}.totals"
}

for family in direct pertype mixed; do
    nlfam=$family
    [ "$family" = pertype ] && nlfam=direct # same namelist; only the copied factypeflux.inp differs
    run_mode serial 1 1 1 "${NAMELIST_DIR}/namoptions.1008.${nlfam}.serial" "$family"
    run_mode xsplit 2 1 2 "${NAMELIST_DIR}/namoptions.1008.${nlfam}.xsplit" "$family"
    run_mode ysplit 1 2 2 "${NAMELIST_DIR}/namoptions.1008.${nlfam}.ysplit" "$family"
    run_mode xysplit 2 2 4 "${NAMELIST_DIR}/namoptions.1008.${nlfam}.xysplit" "$family"
done

# Unsupported wall treatments and a malformed factypeflux.inp stop the run at start-up.
check_abort() {
    local label="$1" expect="$2" sedexpr="$3" fluxfile="$4"
    local run_dir="${ROOT_RUN_DIR}/abort-${label}"
    mkdir -p "$run_dir"
    cp -r "$CASE_SOURCE"/. "$run_dir"/
    sed "$sedexpr" "${NAMELIST_DIR}/namoptions.1008.direct.serial" > "$run_dir/$NAMELIST"
    [ -n "$fluxfile" ] && printf '%s\n' "$fluxfile" > "$run_dir/factypeflux.inp.100"
    local rc=0
    (cd "$run_dir" && mpi_preflight 1 && "$MPIEXEC" $MPI_LAUNCH_EXTRA_ARGS -n 1 "$UDALES_BUILD" "$NAMELIST" > run.log 2>&1) || rc=$?
    if [ "$rc" -eq 0 ] || ! grep -q "$expect" "$run_dir/run.log"; then
        echo "FAIL: TEST_FACET_FLUX abort check [$label] (exit code $rc, expected message: $expect)"
        failures=$((failures + 1))
    fi
}
check_abort iwalltemp3 "iwalltemp must be 1 or 2" 's/^nfcts = 1102/iwalltemp = 3\nnfcts = 1102/' ""
check_abort iwallmoist3 "iwallmoist must be 1 or 2" 's/^nfcts = 1102/iwallmoist = 3\nnfcts = 1102/' ""
check_abort iwallmoist2 "iwallmoist = 2 needs the heat transfer coefficient" 's/^nfcts = 1102/iwallmoist = 2\nnfcts = 1102/' ""
check_abort badfluxfile "cannot read line of factypeflux.inp" 's/^nfcts = 1102/nfcts = 1102/' "1  abc  def"
check_abort unknowntype "not in factypes.inp" 's/^nfcts = 1102/nfcts = 1102/' "77  1.0  1.0"

# The totals must not depend on the decomposition (summation order only: relative 1e-10).
if [ "$failures" -eq 0 ]; then
    for family in direct pertype mixed; do
        for label in ysplit xsplit xysplit; do
            if ! python3 - "${ROOT_RUN_DIR}/${family}-serial.totals" "${ROOT_RUN_DIR}/${family}-${label}.totals" <<'PYEOF'
import sys
a, b = [[float(x) for x in open(p).read().split()] for p in sys.argv[1:3]]
sys.exit(0 if len(a) == len(b) == 10 and all(abs(x - y) <= 1e-10 * max(abs(x), abs(y)) for x, y in zip(a, b)) else 1)
PYEOF
            then
                echo "FAIL: TEST_FACET_FLUX [$family] totals differ between serial and $label"
                cat "${ROOT_RUN_DIR}/${family}-serial.totals" "${ROOT_RUN_DIR}/${family}-${label}.totals"
                failures=$((failures + 1))
            fi
        done
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
