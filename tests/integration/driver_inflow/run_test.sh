#!/bin/bash
# Precursor -> driver smoke test for the inflow module (src/inflow.f90).
#
# Runs case 501 (idriver = 1) to write inlet planes at iplane, feeds them to
# case 502 (idriver = 2), and asserts both complete and that 502 really took its
# inflow from the precursor. Between them this covers initinflow, drivergen,
# writedriverfile, readdriverfile and both branches of exitinflow -- the driver
# lifecycle that examples/949, examples/950 and tests/cases/525 exercise only at
# production size (949 is 256x128x128; 525 asks for 691200 s on 4096 ranks), so
# nothing ever ran it in CI.
#
# Both cases are 8x8x8 on 2 ranks and take seconds, which is what makes this
# runnable as part of the default suite. Vreman closure, matching the SGS model
# every case in experiments/ actually uses. ltempeq/lmoist are on deliberately:
# they are what make BCxT = 3 / BCxq = 3 set lhdriver/lqdriver, so the thl and qt
# driver files get written, read and deallocated too. A neutral case (as in
# experiments/525) would leave those branches untested.

set -u

if ! command -v module >/dev/null 2>&1 && [ -f /etc/profile.d/modules.sh ]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
fi

if command -v module >/dev/null 2>&1; then
    # Keep consistent with the `icl` stack in tools/build_executable.sh: mixing
    # toolchain years does not co-load on CX3.
    module load intel/2021a netCDF/4.8.0-iimpi-2021a netCDF-Fortran/4.5.3-iimpi-2021a FFTW/3.3.9-intel-2021a CMake/3.29.3-GCCcore-13.3.0 git/2.45.1-GCCcore-13.3.0 >/dev/null 2>&1 || true
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
UDALES_BUILD="${UDALES_BUILD:-${REPO_ROOT}/build/debug/u-dales}"
PRECURSOR_SOURCE="${REPO_ROOT}/tests/cases/501"
DRIVER_SOURCE="${REPO_ROOT}/tests/cases/502"

if [ ! -x "$UDALES_BUILD" ]; then
    echo "FAIL: no u-dales executable at $UDALES_BUILD" >&2
    exit 1
fi

if [ -z "${MPIEXEC:-}" ] && command -v mpiifort >/dev/null 2>&1; then
    MPIEXEC="$(dirname "$(command -v mpiifort)")/mpiexec"
else
    MPIEXEC="${MPIEXEC:-mpiexec}"
fi
MPI_LAUNCH_EXTRA_ARGS="${MPI_LAUNCH_EXTRA_ARGS:-}"
MPI_VERSION_OUTPUT="$("$MPIEXEC" --version 2>/dev/null || true)"
if printf '%s\n' "$MPI_VERSION_OUTPUT" | grep -Eqi "Open MPI|OpenRTE"; then
    MPI_LAUNCH_EXTRA_ARGS="--oversubscribe ${MPI_LAUNCH_EXTRA_ARGS}"
fi

WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/udales-driver-inflow-XXXXXX")"
cleanup() { [ "${KEEP_WORKDIR:-0}" = "1" ] || rm -rf "$WORKDIR"; }
trap cleanup EXIT

export HDF5_USE_FILE_LOCKING=FALSE
NPROCS=2

run_case() {
    # $1 = case number, $2 = run directory
    ( cd "$2" && "$MPIEXEC" $MPI_LAUNCH_EXTRA_ARGS -n "$NPROCS" "$UDALES_BUILD" "namoptions.$1" > solver.log 2>&1 )
}

# ---- precursor: writes the inlet planes ----
cp -r "$PRECURSOR_SOURCE" "$WORKDIR/501"
if ! run_case 501 "$WORKDIR/501"; then
    echo "FAIL: precursor 501 (idriver = 1) did not complete" >&2
    tail -25 "$WORKDIR/501/solver.log" >&2
    exit 1
fi

# tdriver is written by driverid 0 only; u/v/w (+ thl/qt here) are per y-rank.
missing=0
for f in tdriver_000.501 udriver_000.501 udriver_001.501 vdriver_000.501 wdriver_000.501 hdriver_000.501 qdriver_000.501; do
    [ -s "$WORKDIR/501/$f" ] || { echo "FAIL: precursor did not write $f" >&2; missing=1; }
done
[ "$missing" -eq 0 ] || exit 1

# ---- driver: consumes them ----
cp -r "$DRIVER_SOURCE" "$WORKDIR/502"
cp "$WORKDIR"/501/*driver_*.501 "$WORKDIR/502/"
if ! run_case 502 "$WORKDIR/502"; then
    echo "FAIL: driver 502 (idriver = 2) did not complete" >&2
    tail -25 "$WORKDIR/502/solver.log" >&2
    exit 1
fi

# BCxT = 3 / BCxq = 3 must have engaged lhdriver/lqdriver. If they had not, the
# solver warns and silently falls back to a non-precursor inflow, which would
# make this test pass while covering less than it claims.
if grep -qi "not given by precursor" "$WORKDIR/502/solver.log"; then
    echo "FAIL: 502 did not take thl/qt inflow from the precursor:" >&2
    grep -i "not given by precursor" "$WORKDIR/502/solver.log" >&2
    exit 1
fi

# The run must have actually advanced, not just initialised and stopped.
for c in 501 502; do
    steps=$(grep -c "Time of Simulation" "$WORKDIR/$c/solver.log" || true)
    [ "${steps:-0}" -ge 1 ] || { echo "FAIL: case $c produced no timesteps" >&2; exit 1; }
done

echo "driver_inflow: PASS (501 idriver=1 -> 502 idriver=2, thl/qt inflow from precursor)"
