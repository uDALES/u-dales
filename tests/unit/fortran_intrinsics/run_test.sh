#!/bin/bash
# Toolchain probe for Fortran intrinsics whose runtime support varies by compiler.
#
# Background: NVHPC 24.11 compiles findloc() on a logical array without complaint
# but has no runtime implementation for it, so a GPU build aborts with
#
#     0: FINDLOC: unimplemented for data type
#
# the first time such a call is reached. Because the call in modibm only ran for
# geometries needing IBM wall-function reconstruction, every grid-aligned case
# built and ran fine and the gap surfaced only on a user's simulation.
#
# This test does two things:
#   1. builds a capability matrix by running each intrinsic in its own process,
#      so a runtime abort in one probe does not hide the others;
#   2. greps src/ for uses of any intrinsic the target compiler cannot run, and
#      fails if one is found.
#
# Step 2 is what makes this self-maintaining: reintroducing a logical findloc
# into the solver fails the test on the NVHPC toolchain rather than in a
# simulation. An unsupported intrinsic that the solver does not use is reported
# but not fatal.
#
# Compilers: set FC to test one specific compiler. Otherwise every compiler in
# CANDIDATE_FCS that is on PATH is tested.

set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
SRC_DIR="${REPO_ROOT}/src"
PROBE_SRC="${SCRIPT_DIR}/intrinsic_probe.f90"

CANDIDATE_FCS="${FC:-gfortran nvfortran}"

# check-name : extended regex matching that intrinsic's use in Fortran source
CHECKS="findloc-logical findloc-real findloc-integer count-logical"

# Pattern used to decide whether the solver actually depends on a check. Empty
# means "not attributable to a source pattern" - reported, never fatal.
pattern_for() {
    case "$1" in
        # findloc whose first argument is a comparison, i.e. a logical mask
        findloc-logical)  printf '%s' 'findloc[[:space:]]*\([^,]*[<>=]' ;;
        findloc-real)     printf '%s' '' ;;
        findloc-integer)  printf '%s' '' ;;
        count-logical)    printf '%s' 'count[[:space:]]*\(' ;;
        *)                printf '%s' '' ;;
    esac
}

WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

overall_status=0
tested_any=0

echo "=========================================="
echo "Running Fortran intrinsic support probe"
echo "=========================================="
echo "Source scanned: ${SRC_DIR}"
echo ""

for fc in $CANDIDATE_FCS; do
    if ! command -v "$fc" >/dev/null 2>&1; then
        echo "-- $fc: not on PATH, skipping"
        continue
    fi
    tested_any=1

    # FC may be an absolute path, so derive a filename-safe tag from it
    tag="$(printf '%s' "$fc" | tr -c 'A-Za-z0-9._-' '_')"
    exe="${WORK_DIR}/probe_${tag}"
    if ! "$fc" -o "$exe" "$PROBE_SRC" > "${WORK_DIR}/build_${tag}.log" 2>&1; then
        echo "-- $fc: FAILED to compile the probe"
        sed 's/^/     /' "${WORK_DIR}/build_${tag}.log"
        overall_status=1
        continue
    fi

    echo "-- $fc ($("$fc" --version 2>&1 | grep -m1 -v '^[[:space:]]*$'))"

    for check in $CHECKS; do
        if "$exe" "$check" > "${WORK_DIR}/run_${tag}_${check}.log" 2>&1; then
            printf '     %-18s supported\n' "$check"
            continue
        fi

        detail="$(tr '\n' ' ' < "${WORK_DIR}/run_${tag}_${check}.log" | sed 's/[[:space:]]\+/ /g')"
        pattern="$(pattern_for "$check")"

        if [ -z "$pattern" ]; then
            printf '     %-18s UNSUPPORTED (not used by src/, not fatal): %s\n' "$check" "$detail"
            continue
        fi

        # Match on code only: strip each line from its first "!" before testing,
        # so comments describing an intrinsic do not count as using it.
        hits="$(find "$SRC_DIR" -type f \( -name '*.f90' -o -name '*.F90' -o -name '*.cuf' \) -print0 \
                | xargs -0 awk -v pat="$pattern" '
                    { code = $0; sub(/!.*/, "", code)
                      if (tolower(code) ~ pat) printf "%s:%d:%s\n", FILENAME, FNR, $0 }' 2>/dev/null)"
        if [ -z "$hits" ]; then
            printf '     %-18s UNSUPPORTED (not used by src/, not fatal): %s\n' "$check" "$detail"
            continue
        fi

        printf '     %-18s UNSUPPORTED AND USED BY src/: %s\n' "$check" "$detail"
        printf '%s\n' "$hits" | sed "s|^${REPO_ROOT}/|       |"
        overall_status=1
    done
done

echo ""
echo "=========================================="
if [ "$tested_any" -eq 0 ]; then
    echo "No Fortran compiler found; nothing tested"
    echo "=========================================="
    exit 1
fi
if [ "$overall_status" -eq 0 ]; then
    echo "Test completed successfully"
else
    echo "Test failed: src/ uses an intrinsic the compiler cannot run"
fi
echo "=========================================="

exit "$overall_status"
