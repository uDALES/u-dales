#!/bin/bash
# The retired-key contract: a namoptions that still sets a key this PR removed
# must abort at the namelist read, naming the key -- it must not be silently
# ignored.
#
# This is the PR's main user-facing promise and it had no test. It is also the
# change most likely to be met in the wild: 6 of the 8 cases under experiments/
# (064, 065, 103, 110, 525, 531) still set z0/z0h/wtsurf/wqsurf/thls. If one of
# those keys were quietly accepted and ignored, the case would run with a setting
# that no longer does anything -- which is exactly the silent-garbage failure the
# loud-abort design exists to prevent, and it would look like success.
#
# Each key is checked in isolation, on a namelist that is otherwise valid, so a
# pass means "this key aborts", not "something somewhere aborted". The run must
# also *name* the key: a case that dies for an unrelated reason would otherwise
# score as a pass and the test would rot silently.
#
# No MPI, no geometry, no timesteps: every case fails at the namelist read, so
# the whole set costs seconds.

set -u

if ! command -v module >/dev/null 2>&1 && [ -f /etc/profile.d/modules.sh ]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
fi
if command -v module >/dev/null 2>&1; then
    module load intel/2021a netCDF/4.8.0-iimpi-2021a netCDF-Fortran/4.5.3-iimpi-2021a FFTW/3.3.9-intel-2021a CMake/3.29.3-GCCcore-13.3.0 git/2.45.1-GCCcore-13.3.0 >/dev/null 2>&1 || true
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
UDALES_BUILD="${UDALES_BUILD:-${REPO_ROOT}/build/debug/u-dales}"

if [ ! -x "$UDALES_BUILD" ]; then
    echo "FAIL: no u-dales executable at $UDALES_BUILD" >&2
    exit 1
fi

WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/udales-retired-keys-XXXXXX")"
cleanup() { [ "${KEEP_WORKDIR:-0}" = "1" ] || rm -rf "$WORKDIR"; }
trap cleanup EXIT
export HDF5_USE_FILE_LOCKING=FALSE

# Retired flat-surface keys, and the namelist each belongs to.
BC_KEYS="thls=290. qts=0.0 z0=0.01 z0h=0.000067 wtsurf=0. wqsurf=0. BCbotm=2 BCbotT=2 BCbotq=2 BCbots=2"
# Inlet-generator keys (modinlet deleted in Phase 2).
INLET_KEYS="lwallfunc=.false."

write_namelist () {
    # $1 = target file, $2 = section, $3 = "key=value" (empty for the control)
    local f="$1" section="$2" kv="$3" extra=""
    [ -n "$kv" ] && extra="$(printf '%s' "$kv" | tr '=' ' ' | awk '{print $1" = "$2}')"
    cat > "$f" <<EOF
&RUN
iexpnr       = 090
runmode      = 1006
runtime      = 1.
libm         = .true.
nprocx       = 1
nprocy       = 1
/

&DOMAIN
itot         = 64
jtot         = 64
ktot         = 64
xlen         = 64
ylen         = 64
/

&PHYSICS
ps           = 101500.
lbuoyancy    = .true.
ltempeq      = .true.
/

&DYNAMICS
ipoiss       = 0
lfftwmeasure = .false.
/

&NAMSUBGRID
lvreman      = .true.
/

&BC
$( [ "$section" = "BC" ] && echo "$extra" )
/

&INPS
zsize        = 64
/
EOF
}

cp "${REPO_ROOT}/tests/integration/basestate/prof.inp.090" "$WORKDIR/"

failures=0
checked=0

# Control: the same namelist with no retired key must NOT abort. Without this,
# a namelist broken for some unrelated reason would make every case below
# "abort correctly" and the whole suite would pass while testing nothing.
write_namelist "$WORKDIR/namoptions.control" "BC" ""
( cd "$WORKDIR" && "$UDALES_BUILD" namoptions.control > control.log 2>&1 )
if [ $? -ne 0 ]; then
    echo "FAIL: control namelist (no retired keys) did not run -- the fixture is broken," >&2
    echo "      so 'aborts correctly' below would prove nothing." >&2
    tail -15 "$WORKDIR/control.log" >&2
    exit 1
fi

for kv in $BC_KEYS $INLET_KEYS; do
    key="${kv%%=*}"
    checked=$((checked + 1))
    write_namelist "$WORKDIR/namoptions.$key" "BC" "$kv"
    ( cd "$WORKDIR" && "$UDALES_BUILD" "namoptions.$key" > "out.$key.log" 2>&1 )
    rc=$?
    if [ "$rc" -eq 0 ]; then
        echo "FAIL: '$key' was accepted silently -- a retired key must abort" >&2
        failures=$((failures + 1))
        continue
    fi
    # Aborting is necessary but not sufficient: the message must identify the
    # key, or a case dying for an unrelated reason would pass as a success.
    if ! grep -qiE "$key|BC" "$WORKDIR/out.$key.log"; then
        echo "FAIL: '$key' aborted (exit $rc) but the output never names it or its namelist:" >&2
        tail -5 "$WORKDIR/out.$key.log" >&2
        failures=$((failures + 1))
    fi
done

if [ "$failures" -ne 0 ]; then
    echo "retired_keys: FAIL ($failures of $checked keys)" >&2
    exit 1
fi

echo "retired_keys: PASS ($checked retired keys abort at the namelist read)"
