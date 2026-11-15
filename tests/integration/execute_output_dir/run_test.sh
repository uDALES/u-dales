#!/bin/bash
# The execute scripts must never run into an existing output directory, and
# copy_inputs.sh must not create one. uDALES never removes, modifies or
# overwrites existing results; that decision stays with the user.
#
# The run itself is replaced by stubs so this needs neither the solver nor MPI:
# MPIEXEC points at a launcher stub that drops "-n N" and executes the command,
# and DA_BUILD at a script that records that it ran.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
TOOLS="${REPO_ROOT}/tools"

ROOT="$(mktemp -d "${TMPDIR:-/tmp}/udales-execute-output-dir-XXXXXX")"
trap 'rm -rf "$ROOT"' EXIT

failures=0
fail() { echo "FAIL: $*"; failures=$((failures + 1)); }
pass() { echo "PASS: $*"; }

# --- fixtures ---------------------------------------------------------------
EXPDIR="$ROOT/experiments"; WORKDIR="$ROOT/outputs"; STUBS="$ROOT/stubs"
mkdir -p "$EXPDIR/001" "$WORKDIR" "$STUBS"

cat > "$STUBS/u-dales" <<'STUB'
#!/bin/bash
echo "stub u-dales ran with $1"
touch ran.marker
STUB
cat > "$STUBS/mpiexec" <<'STUB'
#!/bin/bash
while [ $# -gt 0 ]; do
    case "$1" in
        -n) shift 2 ;;
        *) break ;;
    esac
done
exec "$@"
STUB
chmod +x "$STUBS/u-dales" "$STUBS/mpiexec"

cat > "$EXPDIR/001/config.sh" <<CFG
export DA_EXPDIR="$EXPDIR"
export DA_WORKDIR="$WORKDIR"
export DA_TOOLSDIR="$TOOLS"
export DA_BUILD="$STUBS/u-dales"
export NCPU=1
CFG
printf '&RUN\niexpnr = 001\nlwarmstart = .false.\n/\n' > "$EXPDIR/001/namoptions.001"
export MPIEXEC="$STUBS/mpiexec"

# --- 1. a fresh case runs ---------------------------------------------------
if bash "$TOOLS/local_execute.sh" "$EXPDIR/001" > "$ROOT/run1.log" 2>&1 \
   && [ -f "$WORKDIR/001/ran.marker" ] && [ -f "$WORKDIR/001/output.001.log" ]; then
    pass "local_execute.sh creates the output directory and runs the case"
else
    fail "local_execute.sh did not run a fresh case"; tail -20 "$ROOT/run1.log"
fi

# --- 2. an existing output directory is refused, untouched --------------------
echo "precious results" > "$WORKDIR/001/results.nc"
rm -f "$WORKDIR/001/ran.marker"
before="$(cd "$WORKDIR/001" && ls | sort)"
bash "$TOOLS/local_execute.sh" "$EXPDIR/001" > "$ROOT/run2.out" 2> "$ROOT/run2.err"; rc=$?
after="$(cd "$WORKDIR/001" && ls | sort)"
if [ $rc -ne 0 ] && grep -q "already exists" "$ROOT/run2.err" && grep -q "rerun case 001 on purpose" "$ROOT/run2.err"; then
    pass "local_execute.sh refuses an existing output directory (rc=$rc) with the message"
else
    fail "local_execute.sh did not refuse the existing output directory (rc=$rc)"; cat "$ROOT/run2.err"
fi
if [ "$before" == "$after" ] && [ "$(cat "$WORKDIR/001/results.nc")" == "precious results" ] && [ ! -f "$WORKDIR/001/ran.marker" ]; then
    pass "existing results left untouched and the solver was not started"
else
    fail "the existing output directory was modified"
fi

# --- 3. the same guard in the cluster scripts -------------------------------
for script in hpc_execute.sh archer_execute.sh ecmwf_execute.sh; do
    if grep -q 'if \[ -e "\$outdir" \]; then' "$TOOLS/$script" && grep -q 'already exists' "$TOOLS/$script"; then
        pass "$script carries the same guard"
    else
        fail "$script lacks the existing-output-directory guard"
    fi
done

# --- 4. copy_inputs.sh must not create the output directory -----------------
if bash "$TOOLS/copy_inputs.sh" "$EXPDIR/001" 2 > "$ROOT/copy.log" 2>&1 \
   && [ -f "$EXPDIR/002/namoptions.002" ] && [ ! -e "$WORKDIR/002" ]; then
    pass "copy_inputs.sh prepares the case without creating $WORKDIR/002"
else
    fail "copy_inputs.sh failed or created the output directory"; tail -20 "$ROOT/copy.log"
fi

echo "=========================================="
if [ "$failures" -eq 0 ]; then
    echo "ALL TESTS PASSED: execute scripts protect existing output directories"
    exit 0
fi
echo "TESTS FAILED: $failures"
exit 1
