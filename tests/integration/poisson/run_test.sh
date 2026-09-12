#!/bin/bash
# Poisson solver validation using runmode 1019 on case 103.
#
# Runs the in-solver test tests_poisson (src/tests.f90) for every combination
# of Poisson solver (ipoiss 0 and 3) and decomposition (1x1, 2x1, 1x2, 2x2).
# The test feeds the solver the right-hand side obtained by applying its own
# discrete operator to a known field and requires that field back, on every
# rank, so it is independent of the decomposition and of halo exchange.
#
# Environment:
#   UDALES_BUILD       executable (default build/cpu/debug/u-dales)
#   UDALES_GPU         1 to launch through tools/bind.sh, one GPU per rank
#   UDALES_MAX_RANKS   skip decompositions needing more ranks (default 4;
#                      set 2 on a two-GPU node, 1 on a single GPU)
#   UDALES_POISSON_SOLVERS  ipoiss values to run (default "0 3")
#   MPIEXEC, MPI_LAUNCH_EXTRA_ARGS, UDALES_RUNTIME_MODULES as for the other
#   integration runners

set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
# shellcheck disable=SC1091
source "${REPO_ROOT}/tests/integration/common/runtime_modules.sh"
load_udales_runtime_modules

UDALES_GPU="${UDALES_GPU:-0}"
UDALES_BUILD="${UDALES_BUILD:-${REPO_ROOT}/build/cpu/debug/u-dales}"
UDALES_MAX_RANKS="${UDALES_MAX_RANKS:-4}"
UDALES_POISSON_SOLVERS="${UDALES_POISSON_SOLVERS:-0 3}"
CASE_SOURCE="${CASE_SOURCE:-${REPO_ROOT}/tests/regression/david_tests/cases/103}"
NAMELIST_DIR="${NAMELIST_DIR:-${SCRIPT_DIR}}"
NAMELIST="namoptions.103"
TMPDIR_PARENT="${TMPDIR_PARENT:-}"

# The launcher has to belong to the MPI the executable was built with. The
# GPU harness names them UDALES_GPU_MPIEXEC / UDALES_CPU_MPIEXEC; honour those
# first, then MPIEXEC, then the Intel wrapper next to mpiifort.
if [ "$UDALES_GPU" = "1" ] && [ -n "${UDALES_GPU_MPIEXEC:-}" ]; then
    MPIEXEC="$UDALES_GPU_MPIEXEC"
elif [ "$UDALES_GPU" != "1" ] && [ -n "${UDALES_CPU_MPIEXEC:-}" ]; then
    MPIEXEC="$UDALES_CPU_MPIEXEC"
elif [ -z "${MPIEXEC:-}" ] && command -v mpiifort >/dev/null 2>&1; then
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
# Runs happen inside temporary directories, so a relative path would not
# survive the cd.
UDALES_BUILD="$(readlink -f "$UDALES_BUILD")"
if [ ! -d "$CASE_SOURCE" ]; then
    echo "ERROR: Case source not found: $CASE_SOURCE"
    exit 1
fi

BIND=()
if [ "$UDALES_GPU" = "1" ]; then
    BIND=(bash "${REPO_ROOT}/tools/bind.sh")
fi

if [ -n "$TMPDIR_PARENT" ]; then
    ROOT_RUN_DIR="$(mktemp -d "${TMPDIR_PARENT%/}/poisson-XXXXXX")"
else
    ROOT_RUN_DIR="$(mktemp -d)"
fi
trap 'rm -rf "$ROOT_RUN_DIR"' EXIT

failures=0
runs=0
skipped=0

run_mode() {
    local ipoiss="$1"
    local label="$2"
    local npx="$3"
    local npy="$4"
    local np=$((npx * npy))
    local namelist_source="${NAMELIST_DIR}/namoptions.1019.${label}"
    local run_dir="${ROOT_RUN_DIR}/ipoiss${ipoiss}-${label}"

    if [ "$np" -gt "$UDALES_MAX_RANKS" ]; then
        echo "SKIP: TEST_POISSON [ipoiss=$ipoiss $label] needs $np ranks, UDALES_MAX_RANKS=$UDALES_MAX_RANKS"
        skipped=$((skipped + 1))
        return
    fi
    if [ ! -f "$namelist_source" ]; then
        echo "ERROR: Namelist source not found: $namelist_source"
        exit 1
    fi

    mkdir -p "$run_dir"
    cp -r "$CASE_SOURCE"/. "$run_dir"/
    # The committed namelists carry ipoiss = 3; the solver under test is a
    # parameter of this runner.
    sed "s/^ipoiss .*/ipoiss       = ${ipoiss}/" "$namelist_source" > "$run_dir/$NAMELIST"

    echo "=========================================="
    echo "Running TEST_POISSON [ipoiss=$ipoiss $label]"
    echo "MPI processes: $np   nprocx/nprocy: $npx/$npy   executable: $UDALES_BUILD"
    echo "Run directory: $run_dir"
    echo "=========================================="
    runs=$((runs + 1))

    local run_rc=0
    (
        cd "$run_dir" || exit 1
        # ${BIND[@]+"${BIND[@]}"}: an empty array under set -u is an error in
        # bash 3.2 (macOS); this expands to nothing there and to the array
        # elsewhere.
        "$MPIEXEC" $MPI_LAUNCH_EXTRA_ARGS -n "$np" ${BIND[@]+"${BIND[@]}"} "$UDALES_BUILD" "$NAMELIST" > run.log 2>&1
    ) || run_rc=$?

    grep -E '^rank |^max \|rhs\|' "${run_dir}/run.log" || true
    if [ "$run_rc" -ne 0 ]; then
        echo "FAIL: TEST_POISSON [ipoiss=$ipoiss $label] exited with code $run_rc"
        tail -n 60 "${run_dir}/run.log" || true
        failures=$((failures + 1))
        return
    fi
    if ! grep -q "ALL TESTS PASSED: tests_poisson" "${run_dir}/run.log"; then
        echo "FAIL: TEST_POISSON [ipoiss=$ipoiss $label]"
        tail -n 60 "${run_dir}/run.log" || true
        failures=$((failures + 1))
        return
    fi
    echo "PASS: TEST_POISSON [ipoiss=$ipoiss $label]"
}

for ipoiss in $UDALES_POISSON_SOLVERS; do
    run_mode "$ipoiss" serial  1 1
    run_mode "$ipoiss" xsplit  2 1
    run_mode "$ipoiss" ysplit  1 2
    run_mode "$ipoiss" xysplit 2 2
done

echo "=========================================="
if [ "$failures" -eq 0 ]; then
    echo "All TEST_POISSON runs passed ($runs run, $skipped skipped)"
    echo "=========================================="
    exit 0
else
    echo "TEST_POISSON failures: $failures of $runs ($skipped skipped)"
    echo "Preserving run directory: $ROOT_RUN_DIR"
    trap - EXIT
    echo "=========================================="
    exit 1
fi
