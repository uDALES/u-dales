#!/bin/bash
# Analytic checks of the derived hydrostatic base state (src/modbasestate.f90,
# src/modthermodynamics.f90; issue #302), via two runmodes:
#
#   1006 TEST_BASESTATE   -- initbasestate vs the closed-form hydrostatic column,
#                            in uniform, linearly stratified and moist regimes.
#   1007 TEST_BURIED      -- fromztop fed the base profiles (what diagfld's
#                            fully-solid-slab fallback hands it) must reproduce
#                            that reference column: the section 6.6 continuation.
#
# Both dispatch straight after initglobal, so the vertical grid exists but
# readinitfiles has not run. The tests supply their own profiles and call the
# real routines directly, comparing with an independent expectation rather than
# with a previous run. That is the part cases 090/091/092 cannot do: a regression
# comparison only shows the base state has not *changed*, and the 091 assert only
# shows it was *printed*.
#
# Single rank, no geometry, no timesteps: each exits at the dispatch point.

set -u

if ! command -v module >/dev/null 2>&1 && [ -f /etc/profile.d/modules.sh ]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
fi

if command -v module >/dev/null 2>&1; then
    # Keep consistent with the `icl` stack in tools/build_executable.sh.
    module load intel/2021a netCDF/4.8.0-iimpi-2021a netCDF-Fortran/4.5.3-iimpi-2021a FFTW/3.3.9-intel-2021a CMake/3.29.3-GCCcore-13.3.0 git/2.45.1-GCCcore-13.3.0 >/dev/null 2>&1 || true
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
UDALES_BUILD="${UDALES_BUILD:-${REPO_ROOT}/build/debug/u-dales}"

if [ ! -x "$UDALES_BUILD" ]; then
    echo "FAIL: no u-dales executable at $UDALES_BUILD" >&2
    exit 1
fi

# Deliberately no mpiexec: these are single-rank checks that exit at the runmode
# dispatch, and MPI singleton init covers them. Going through mpiexec would add a
# hydra bootstrap that fails on some CX3 nodes, producing a failure that looks
# exactly like a real test failure.

WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/udales-basestate-XXXXXX")"
cleanup() { [ "${KEEP_WORKDIR:-0}" = "1" ] || rm -rf "$WORKDIR"; }
trap cleanup EXIT
export HDF5_USE_FILE_LOCKING=FALSE

# runmode -> (namelist, prof.inp grid file, PASS marker printed by the routine)
run_one() {
    local namelist="$1" prof="$2" marker="$3"
    # initglobal reads zf from prof.inp before the runmode dispatch, so the grid
    # file must be present even though these tests never use the profile itself.
    cp "${SCRIPT_DIR}/${namelist}" "${SCRIPT_DIR}/${prof}" "$WORKDIR/"
    ( cd "$WORKDIR" && "$UDALES_BUILD" "$namelist" > "out.${namelist}.log" 2>&1 )
    local rc=$?
    if [ "$rc" -ne 0 ]; then
        echo "FAIL: ${marker} run exited $rc" >&2
        grep -E "FAIL|${marker}" "$WORKDIR/out.${namelist}.log" >&2 || tail -25 "$WORKDIR/out.${namelist}.log" >&2
        return 1
    fi
    # The routine must actually have been dispatched. Without this a namelist typo
    # that fell through to a normal 1 s run would exit 0 and "pass".
    if ! grep -q "${marker}: PASS" "$WORKDIR/out.${namelist}.log"; then
        echo "FAIL: ${marker} did not run (runmode not dispatched?)" >&2
        tail -25 "$WORKDIR/out.${namelist}.log" >&2
        return 1
    fi
    return 0
}

run_one namoptions.1006 prof.inp.090 tests_basestate          || exit 1
run_one namoptions.1007 prof.inp.091 tests_buried_continuation || exit 1

echo "basestate: PASS (closed-form base state + buried-slab continuation)"
