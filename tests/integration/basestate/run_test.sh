#!/bin/bash
# Unit check of the derived hydrostatic base state (src/modbasestate.f90, #302),
# via runmode 1006 (TEST_BASESTATE).
#
# The solver dispatches this straight after initglobal, so the vertical grid
# exists but readinitfiles has not run. tests_basestate therefore supplies its
# own profiles and calls initbasestate directly, and compares the result with the
# closed-form hydrostatic solution rather than with a previous run. That is the
# part cases 090/091/092 cannot do: a regression comparison only shows the base
# state has not *changed*, and the 091 assert only shows it was *printed*.
#
# Single rank, no geometry, no timesteps: it exits at the dispatch point.

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
NAMELIST="${NAMELIST:-namoptions.1006}"

if [ ! -x "$UDALES_BUILD" ]; then
    echo "FAIL: no u-dales executable at $UDALES_BUILD" >&2
    exit 1
fi

# Deliberately no mpiexec: this is a single-rank check that exits at the runmode
# dispatch, and MPI singleton init covers it. Going through mpiexec would add a
# hydra bootstrap that fails on some CX3 nodes, producing a failure that looks
# exactly like a real test failure.

WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/udales-basestate-XXXXXX")"
cleanup() { [ "${KEEP_WORKDIR:-0}" = "1" ] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

# initglobal reads zf from prof.inp before the runmode dispatch, so the grid
# file must be present even though this test never uses the profile itself.
cp "${SCRIPT_DIR}/${NAMELIST}" "${SCRIPT_DIR}/prof.inp.090" "$WORKDIR/"

export HDF5_USE_FILE_LOCKING=FALSE

( cd "$WORKDIR" && "$UDALES_BUILD" "$NAMELIST" > out.log 2>&1 )
rc=$?

if [ "$rc" -ne 0 ]; then
    echo "FAIL: tests_basestate reported failure (exit $rc)" >&2
    grep -E "FAIL|tests_basestate" "$WORKDIR/out.log" >&2 || tail -25 "$WORKDIR/out.log" >&2
    exit 1
fi

# The runmode must actually have been dispatched. Without this, a namelist typo
# that silently fell through to a normal 1 s run would exit 0 and "pass".
if ! grep -q "tests_basestate: PASS" "$WORKDIR/out.log"; then
    echo "FAIL: tests_basestate did not run (runmode not dispatched?)" >&2
    tail -25 "$WORKDIR/out.log" >&2
    exit 1
fi

echo "basestate: PASS (hydrostatic base state matches closed form)"
