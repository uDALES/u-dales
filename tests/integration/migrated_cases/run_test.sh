#!/bin/bash
# Smoke test for the two cases Phase 1 migrated off the flat-surface (lbottom)
# scheme onto ground facets.
#
#   103  fixed-temperature ground: the old BCbotT=2 / thls=288 case, now a ground
#        facet held at 288 K (Tfacinit=288, iwalltemp=2). 8x8x8 on 2 ranks.
#
# Until now neither had any automated coverage: 103 lived under
# tests/regression/david_tests, whose run_and_compare was commented out and
# pointed at directories that do not exist, so it never ran a case; the harness
# has since been removed. This asserts the migration's intent rather than
# statistical equivalence against the retired scheme (which would need a pre-PR
# build and a criterion nobody has written down):
#
#   - the case runs to completion and advances in time;
#   - the base state is derived, and for 103 it is the 288 K the fixed-temperature
#     ground migrated to (this is the #302 value that the old thls sentinel used
#     to carry);
#   - no dumped field is NaN or the nodata marker (-999.): the ground-facet floor
#     produces finite fields, not the sentinels a broken bottom BC would leave.
#
# examples/999 (the neutral 2-facet floor, the other migrated case) is 128^3 and
# is exercised as a full-size example rather than here; see tests/cases/103.

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

if [ -z "${MPIEXEC:-}" ] && command -v mpiifort >/dev/null 2>&1; then
    MPIEXEC="$(dirname "$(command -v mpiifort)")/mpiexec"
else
    MPIEXEC="${MPIEXEC:-mpiexec}"
fi
MPI_LAUNCH_EXTRA_ARGS="${MPI_LAUNCH_EXTRA_ARGS:-}"
if printf '%s\n' "$("$MPIEXEC" --version 2>/dev/null || true)" | grep -Eqi "Open MPI|OpenRTE"; then
    MPI_LAUNCH_EXTRA_ARGS="--oversubscribe ${MPI_LAUNCH_EXTRA_ARGS}"
fi

WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/udales-migrated-XXXXXX")"
cleanup() { [ "${KEEP_WORKDIR:-0}" = "1" ] || rm -rf "$WORKDIR"; }
trap cleanup EXIT
export HDF5_USE_FILE_LOCKING=FALSE

run_case_103() {
    local src="${REPO_ROOT}/tests/cases/103" run="$WORKDIR/103"
    cp -r "$src" "$run"
    local np
    np=$(( $(grep -oiE "nprocx *= *[0-9]+" "$run/namoptions.103" | grep -oE "[0-9]+") * \
           $(grep -oiE "nprocy *= *[0-9]+" "$run/namoptions.103" | grep -oE "[0-9]+") ))
    ( cd "$run" && "$MPIEXEC" $MPI_LAUNCH_EXTRA_ARGS -n "$np" "$UDALES_BUILD" namoptions.103 > solver.log 2>&1 )
    if [ $? -ne 0 ]; then
        echo "FAIL: case 103 did not complete" >&2; tail -20 "$run/solver.log" >&2; return 1
    fi
    if [ "$(grep -c 'Time of Simulation' "$run/solver.log")" -lt 1 ]; then
        echo "FAIL: case 103 produced no timesteps" >&2; return 1
    fi
    # The migration replaced a fixed 288 K bottom BC with a 288 K ground facet.
    # The derived base state must report that value.
    if ! grep -Eq "Base state:.*288\.0" "$run/solver.log"; then
        echo "FAIL: case 103 base state is not 288 K (fixed-temperature migration):" >&2
        grep -i "Base state" "$run/solver.log" >&2; return 1
    fi
    # No NaN in any dump; and no nodata marker in the per-cell prognostic/thermo
    # dumps (tdump/fielddump), which the ground-facet fallback keeps finite. The
    # averaged dumps (xydump/ydump/...) may legitimately carry the marker in rows
    # or columns that hold no fluid, so they are checked for NaN only.
    python3 - "$run" <<'PY'
import sys, glob, os
import numpy as np
try:
    import netCDF4 as nc
except Exception as e:
    print(f"FAIL: cannot import netCDF4: {e}", file=sys.stderr); sys.exit(1)
run = sys.argv[1]
bad = []
for path in sorted(glob.glob(os.path.join(run, "*.nc"))):
    base = os.path.basename(path)
    per_cell = base.startswith("tdump") or base.startswith("fielddump")
    with nc.Dataset(path) as ds:
        for name, var in ds.variables.items():
            if var.dtype.kind != "f":
                continue
            a = np.asarray(var[:], dtype=np.float64)
            if a.size == 0:
                continue
            if np.isnan(a).any():
                bad.append(f"{base}:{name} has NaN")
            if per_cell and np.any(np.isclose(a, -999.0, atol=1e-6)):
                bad.append(f"{base}:{name} has the nodata marker (per-cell dump must be finite)")
if bad:
    print("FAIL: case 103 dumps contain sentinels/NaN:", file=sys.stderr)
    for b in bad[:10]:
        print("  " + b, file=sys.stderr)
    sys.exit(1)
PY
    return $?
}

run_case_103 || exit 1
echo "migrated_cases: PASS (103 fixed-temperature ground facet, 288 K base state, no sentinels)"
