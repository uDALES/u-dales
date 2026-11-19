#!/usr/bin/env bash

set -e

# Resolve paths relative to this script's location, regardless of where it is called from
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Script usage function
usage() {
    echo "Usage: $0 <-m|-p> [case1 case2 ...] <ref_data_path> [--tolerance <val>] [--system <common|icl|ecmwf>] [--foreground]"
    echo ""
    echo "Arguments:"
    echo "  -m             : Generate inputs with the MATLAB preprocessing route (common and icl only)"
    echo "  -p             : Generate inputs with the Python preprocessing route"
    echo "  case1 ...      : Optional list of experiment cases (e.g., 100 201 224)"
    echo "                   If not specified, uses default cases: ${DEFAULT_TEST_CASES[*]}"
    echo "  ref_data_path  : Path to the reference data directory (required; last positional argument)"
    echo "  --tolerance    : Max absolute error for numeric files (default: 1e-10)"
    echo "  --system       : System the suite runs on: common, icl or ecmwf (default: common)."
    echo "                   The bare word common, icl or ecmwf among the arguments is accepted too."
    echo "  --foreground   : Run in the current shell and keep the output on the terminal"
    echo "                   (the run ends with the terminal; on ecmwf it runs on the current node)."
    echo ""
    echo "How the suite runs:"
    echo "  common, icl : detached from the terminal (survives closing it); console output goes to"
    echo "                logdir/test_inputs.console.log and the detailed log to logdir/test_inputs.log."
    echo "                When the output is already redirected (not a terminal) it runs in place."
    echo "  ecmwf       : submitted as one Slurm job (qos nf) that builds, writes and compares;"
    echo "                the job output is logdir/test_inputs.console.log."
    echo "                Settings: PREPROC_ACCOUNT (default spnlthee), PREPROC_QOS (nf),"
    echo "                PREPROC_WALLTIME (04:00:00), PREPROC_MEM (32gb)."
    echo ""
    echo "Phases:"
    echo "  PHASE 1 : Build preprocessing applications"
    echo "  PHASE 2 : Run write_inputs.sh for each case to (re-)generate input files"
    echo "  PHASE 3 : Compare generated input files against the reference directory"
    echo ""
    echo "Examples:"
    echo "  $0 -m ref_data/                          # Run default cases with MATLAB route"
    echo "  $0 -m 990 991 ref_data/                  # Run specific cases with MATLAB route"
    echo "  $0 -p 995 ref_data/                      # Generate inputs with Python route"
    echo "  $0 -m 992 ref_data/ --tolerance 1e-8"
    echo "  $0 -m 993 ref_data/ --system icl"
    echo "  $0 -p 993 /path/to/ref_data ecmwf        # Submit the suite as a Slurm job on ECMWF"
    exit 1
}

DEFAULT_TEST_CASES=(990 991 992 993 994 995 997 998 999 807)
PREPROCESS_ROUTE=""
PREPROCESS_ROUTE_NAME=""
PREPROCESS_ROUTE_SET=false

set_preprocess_route() {
    if [ "$PREPROCESS_ROUTE_SET" = true ]; then
        echo "[ERROR] Choose only one preprocessing route: -m or -p"
        usage
    fi

    if [ "$1" = "-m" ]; then
        PREPROCESS_ROUTE="-m"
        PREPROCESS_ROUTE_NAME="MATLAB"
    else
        PREPROCESS_ROUTE="-p"
        PREPROCESS_ROUTE_NAME="Python"
    fi
    PREPROCESS_ROUTE_SET=true
}

# Parse arguments
if [ $# -eq 0 ]; then
    usage
fi

# The arguments as given, for re-launching the suite detached or inside a job.
ORIG_ARGS=("$@")

TEST_CASES=()
POSITIONAL_ARGS=()
TOLERANCE="1e-10"
SYSTEM="common"
FOREGROUND=false
IN_RUN=false        # internal: set when re-launched detached or inside the job

while [ $# -gt 0 ]; do
    case "$1" in
        -m)
            set_preprocess_route "$1"
            ;;
        -p)
            set_preprocess_route "$1"
            ;;
        --tolerance)
            shift
            if [ $# -eq 0 ]; then
                echo "[ERROR] --tolerance requires a value"
                usage
            fi
            TOLERANCE="$1"
            ;;
        --system)
            shift
            if [ $# -eq 0 ]; then
                echo "[ERROR] --system requires a value"
                usage
            fi
            SYSTEM="$1"
            ;;
        --foreground)
            FOREGROUND=true
            ;;
        --in-run)
            IN_RUN=true
            ;;
        --*)
            echo "Unknown option: $1"
            usage
            ;;
        common|icl|ecmwf)
            # A bare system word, as for build_executable.sh, unless it names a
            # directory (then it is a reference data path).
            if [ -d "$1" ]; then
                POSITIONAL_ARGS+=("$1")
            else
                SYSTEM="$1"
            fi
            ;;
        *)
            POSITIONAL_ARGS+=("$1")
            ;;
    esac
    shift
done

if [ "$PREPROCESS_ROUTE_SET" = false ]; then
    echo "[ERROR] Preprocessing route must be set with either -m or -p"
    usage
fi

if [ ${#POSITIONAL_ARGS[@]} -eq 0 ]; then
    echo "[ERROR] Reference data path is required"
    usage
fi

last_positional_index=$((${#POSITIONAL_ARGS[@]} - 1))
REF_DATA_PATH="${POSITIONAL_ARGS[$last_positional_index]}"
TEST_CASES=("${POSITIONAL_ARGS[@]:0:$last_positional_index}")

# Validate --system value
if [ "$SYSTEM" != "common" ] && [ "$SYSTEM" != "icl" ] && [ "$SYSTEM" != "ecmwf" ]; then
    echo "[ERROR] --system must be 'common', 'icl' or 'ecmwf', got: '$SYSTEM'"
    usage
fi

if [ "$SYSTEM" = "ecmwf" ] && [ "$PREPROCESS_ROUTE" = "-m" ]; then
    echo "[ERROR] The MATLAB preprocessing route is not available on ecmwf; use -p"
    usage
fi

if [ ${#TEST_CASES[@]} -eq 0 ]; then
    TEST_CASES=("${DEFAULT_TEST_CASES[@]}")
fi

# The reference path must survive a change of working directory (detached
# re-launch, Slurm job), so make it absolute when it exists.
if [ -d "$REF_DATA_PATH" ]; then
    REF_DATA_PATH="$(cd "$REF_DATA_PATH" && pwd)"
fi

# Color codes (only on a terminal; logs stay free of escape codes)
if [ -t 1 ]; then
    RED='\033[0;31m'
    GREEN='\033[0;32m'
    YELLOW='\033[1;33m'
    NC='\033[0m'
    IS_TTY=true
else
    RED=''
    GREEN=''
    YELLOW=''
    NC=''
    IS_TTY=false
fi

print_green()  { echo -e "${GREEN}$1${NC}"; }
print_red()    { echo -e "${RED}$1${NC}"; }
print_yellow() { echo -e "${YELLOW}$1${NC}"; }

LOGDIR="$REPO_ROOT/logdir"
mkdir -p "$LOGDIR"
LOG_FILE="$LOGDIR/test_inputs.log"
CONSOLE_LOG="$LOGDIR/test_inputs.console.log"

# Default Slurm settings for the ecmwf job (override in the environment).
PREPROC_ACCOUNT="${PREPROC_ACCOUNT:-spnlthee}"
PREPROC_QOS="${PREPROC_QOS:-nf}"
PREPROC_WALLTIME="${PREPROC_WALLTIME:-04:00:00}"
PREPROC_MEM="${PREPROC_MEM:-32gb}"

# Largest nompthreads over the selected cases: the CPU request of the job.
# write_inputs.sh re-reads each case's own value at run time.
max_nompthreads() {
    local max=8
    local case_num value
    for case_num in "${TEST_CASES[@]}"; do
        value=$(grep -ih "nompthreads" "$SCRIPT_DIR/experiments/$case_num/namoptions.$case_num" 2>/dev/null \
            | head -n 1 | sed -e 's/.*=//' -e 's/[!,/].*//' -e 's/[[:space:]]//g')
        if [[ "$value" =~ ^[0-9]+$ ]] && [ "$value" -gt "$max" ]; then
            max=$value
        fi
    done
    echo "$max"
}

# ----------------------------------------------------------------------------
# Launch: detach from the terminal (common, icl) or submit a Slurm job (ecmwf),
# unless asked to run in the foreground or already re-launched.
# ----------------------------------------------------------------------------
if [ "$IN_RUN" = false ] && [ "$FOREGROUND" = false ]; then
    if [ "$SYSTEM" = "ecmwf" ]; then
        if ! command -v sbatch >/dev/null 2>&1; then
            echo "[ERROR] sbatch not found; --system ecmwf submits a Slurm job. Use --foreground inside an allocation."
            exit 1
        fi
        if [[ "$REPO_ROOT" =~ [[:space:]] ]]; then
            echo "[ERROR] The repository path must not contain spaces: ECMWF's sbatch cannot submit from such a directory"
            exit 1
        fi
        if [[ ! "$PREPROC_MEM" =~ ^[0-9]+gb$ ]]; then
            echo "[ERROR] PREPROC_MEM must be set like 32gb (number followed by lowercase gb)"
            exit 1
        fi
        ncpus=$(max_nompthreads)
        mem_directive=""
        if [ "$PREPROC_QOS" != "np" ]; then
            mem_directive="#SBATCH --mem=${PREPROC_MEM%gb}G"
        fi
        # Re-launch inside the job with the original arguments, quoted for bash.
        quoted_args=""
        for arg in "${ORIG_ARGS[@]}"; do
            quoted_args+=" $(printf '%q' "$arg")"
        done
        # The job re-launches from the directory the suite was called from, so
        # a relative reference path keeps its meaning.
        LAUNCH_DIR="$(pwd)"
        JOB_FILE="$LOGDIR/test_inputs.slurm"
        cat <<EOF > "$JOB_FILE"
#!/bin/bash
#SBATCH --job-name=ud_test_preprocess
#SBATCH --account=$PREPROC_ACCOUNT
#SBATCH --qos=$PREPROC_QOS
#SBATCH --time=$PREPROC_WALLTIME
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$ncpus
#SBATCH --hint=nomultithread
#SBATCH --output=test_inputs.console.log
${mem_directive}
set -o pipefail

## Build and runtime environment of the preprocessing tools; keep in step with
## the ecmwf block of tools/build_preprocessing.sh.
module load prgenv/gnu gcc/14.2.0 python3/3.13.16-01 cmake/4.2.4

export PREPROC_WALLTIME="$PREPROC_WALLTIME"
export PREPROC_MEM="$PREPROC_MEM"
export PREPROC_ACCOUNT="$PREPROC_ACCOUNT"
export PREPROC_QOS="$PREPROC_QOS"

cd "$LAUNCH_DIR"
bash "$SCRIPT_DIR/ud_test_preprocess.sh"${quoted_args} --in-run
EOF
        # Submitted from logdir so that the path-free --output lands there
        # (ECMWF's sbatch wrapper splits directives on whitespace).
        rm -f "$CONSOLE_LOG"
        cd "$LOGDIR"
        if ! sbatch "$JOB_FILE"; then
            echo "[ERROR] sbatch failed; job script kept at $JOB_FILE" >&2
            exit 1
        fi
        echo "Preprocessing test suite submitted as a Slurm job (account $PREPROC_ACCOUNT, qos $PREPROC_QOS, $ncpus cpus, $PREPROC_MEM, $PREPROC_WALLTIME)."
        echo "  Job script : $JOB_FILE"
        echo "  Console log: $CONSOLE_LOG   (follow with: tail -f $CONSOLE_LOG)"
        echo "  Detail log : $LOG_FILE"
        exit 0
    elif [ "$IS_TTY" = true ]; then
        # Detach so that the suite survives the terminal closing. Output that
        # would have gone to the terminal goes to the console log instead.
        # setsid -f (util-linux) puts the child in its own session before it
        # execs anything; without it (macOS) nohup alone has to do.
        if setsid -f true 2>/dev/null; then
            setsid -f nohup bash "$SCRIPT_DIR/ud_test_preprocess.sh" "${ORIG_ARGS[@]}" --in-run > "$CONSOLE_LOG" 2>&1 < /dev/null
        else
            nohup bash "$SCRIPT_DIR/ud_test_preprocess.sh" "${ORIG_ARGS[@]}" --in-run > "$CONSOLE_LOG" 2>&1 < /dev/null &
            disown 2>/dev/null || true
        fi
        # Give the child time to detach before this shell exits: a hangup that
        # arrives before nohup/setsid have run would still kill it.
        sleep 1
        echo "Preprocessing test suite running detached; it continues if this terminal is closed."
        echo "  Console log: $CONSOLE_LOG   (follow with: tail -f $CONSOLE_LOG; the PID is on its first line)"
        echo "  Detail log : $LOG_FILE"
        echo "  Pass --foreground to run in this shell instead."
        exit 0
    fi
    # Output is not a terminal (already redirected or nohup'ed): run in place.
fi

# ----------------------------------------------------------------------------
# The suite itself
# ----------------------------------------------------------------------------

# Spinner (terminal only)
show_spinner() {
    local pid=$1
    local message=$2
    local delay=0.1
    local spinstr='|/-\'
    while kill -0 "$pid" 2>/dev/null; do
        local temp=${spinstr#?}
        printf "\r[%c] %s" "$spinstr" "$message"
        local spinstr=$temp${spinstr%"$temp"}
        sleep $delay
    done
    printf "\r"
}

# Wait for a background pid; spinner on a terminal, quiet otherwise. Sets WAIT_STATUS.
wait_quietly() {
    local pid=$1
    local message=$2
    if [ "$IS_TTY" = true ]; then
        show_spinner "$pid" "$message"
    fi
    local errexit_was_set=false
    case $- in
        *e*) errexit_was_set=true ;;
    esac
    set +e
    wait "$pid"
    WAIT_STATUS=$?
    if [ "$errexit_was_set" = true ]; then
        set -e
    fi
}

# Run a command; exits on failure
# Optional 4th arg: suppress_ok=true skips the "[OK]...Done" message
run_with_status_strict() {
    local cmd="$1"
    local description="$2"
    local log_file="$3"
    local suppress_ok="${4:-false}"
    local start_s=$SECONDS

    echo "Starting: $description"

    bash -c "$cmd" >> "$log_file" 2>&1 &
    wait_quietly $! "$description"
    local exit_code=$WAIT_STATUS

    if [ $exit_code -eq 0 ]; then
        [ "$suppress_ok" != "true" ] && printf "[OK] %s - Done (%ds)\n" "$description" "$((SECONDS - start_s))"
        return 0
    else
        printf "[FAILED] %s - Exit Code: %d\n" "$description" "$exit_code"
        print_red "FAILED: $description (exit code: $exit_code)"
        echo "Check $log_file for details."
        exit $exit_code
    fi
}

# Run a command; does not exit on failure — returns exit code
run_with_status() {
    local cmd="$1"
    local description="$2"
    local log_file="$3"
    local start_s=$SECONDS

    echo "Starting: $description"

    bash -c "$cmd" >> "$log_file" 2>&1 &
    wait_quietly $! "$description"
    local exit_code=$WAIT_STATUS

    if [ $exit_code -eq 0 ]; then
        printf "[OK] %s - Done (%ds)\n" "$description" "$((SECONDS - start_s))"
    else
        printf "[FAILED] %s - Exit Code: %d (%ds)\n" "$description" "$exit_code" "$((SECONDS - start_s))"
    fi
    return $exit_code
}

# Strip terminal control sequences from a log before grepping it.
strip_log() {
    if command -v col >/dev/null 2>&1; then
        col -b < "$1"
    else
        cat "$1"
    fi
}

ERROR_PATTERN="^Error |^Error:|Undefined function|Undefined variable|not enough input|out of memory|Traceback|Exception|RuntimeError|ValueError|ModuleNotFoundError"

# Log file
echo "u-dales Input Write & Comparison Started at $(date)" > "$LOG_FILE"
echo "========================================" >> "$LOG_FILE"
echo "  Test experiment path : $SCRIPT_DIR/experiments/" >> "$LOG_FILE"
echo "  Reference path       : $REF_DATA_PATH" >> "$LOG_FILE"
echo "  Cases                : ${TEST_CASES[*]}" >> "$LOG_FILE"
echo "  Tolerance            : $TOLERANCE" >> "$LOG_FILE"
echo "  System               : $SYSTEM" >> "$LOG_FILE"
echo "  Host                 : $(hostname)" >> "$LOG_FILE"
echo "  Preprocessing route  : $PREPROCESS_ROUTE_NAME ($PREPROCESS_ROUTE)" >> "$LOG_FILE"
echo "========================================" >> "$LOG_FILE"

echo "Starting u-dales input write & comparison (PID $$)..."
echo ""
echo "========================================"
echo "  Run settings"
echo "========================================"
echo "  Test experiment path : $SCRIPT_DIR/experiments/"
echo "  Reference path       : $REF_DATA_PATH"
echo "  Cases                : ${TEST_CASES[*]}"
echo "  Tolerance            : $TOLERANCE"
echo "  System               : $SYSTEM"
echo "  Host                 : $(hostname)"
echo "  Preprocessing route  : $PREPROCESS_ROUTE_NAME ($PREPROCESS_ROUTE)"
echo "  Log file             : $LOG_FILE"
echo "  Started              : $(date)"
echo "========================================"
echo ""

# Pre-flight checks
echo "========================================"
echo "  Pre-flight checks"
echo "========================================"

preflight_ok=true

if [ ! -f "$REPO_ROOT/tools/build_preprocessing.sh" ]; then
    print_red "[MISSING] build_preprocessing.sh not found in: $REPO_ROOT/tools"
    preflight_ok=false
fi

if [ ! -f "$REPO_ROOT/tools/write_inputs.sh" ]; then
    print_red "[MISSING] write_inputs.sh not found in: $REPO_ROOT/tools"
    preflight_ok=false
fi

if [ ! -f "$REPO_ROOT/tools/ud_compare_inputs.py" ]; then
    print_red "[MISSING] ud_compare_inputs.py not found in: $REPO_ROOT/tools"
    preflight_ok=false
fi

VENV_DIR="$REPO_ROOT/tools/python/.venv"

if [ ! -d "$VENV_DIR" ]; then
    print_red "[MISSING] Python virtual environment not found: $VENV_DIR"
    print_red "         Run 'bash $REPO_ROOT/tools/python/setup_venv.sh $SYSTEM preprocessing_tools' to set it up."
    preflight_ok=false
elif [ ! -f "$VENV_DIR/bin/python3" ] && [ ! -f "$VENV_DIR/bin/python" ]; then
    print_red "[MISSING] Python interpreter not found inside virtual environment: $VENV_DIR"
    preflight_ok=false
fi

for case_num in "${TEST_CASES[@]}"; do
    if [ ! -d "$SCRIPT_DIR/experiments/$case_num" ]; then
        print_red "[MISSING] Experiment input directory not found: $SCRIPT_DIR/experiments/$case_num"
        preflight_ok=false
    fi
    if [ ! -d "$REF_DATA_PATH/$case_num" ]; then
        print_red "[MISSING] Reference directory not found: $REF_DATA_PATH/$case_num"
        preflight_ok=false
    fi
done

if [ "$preflight_ok" = false ]; then
    echo ""
    print_red "Pre-flight checks FAILED. Resolve the issues above before running."
    exit 1
fi

# Resolve Python interpreter inside the venv
if [ -f "$VENV_DIR/bin/python3" ]; then
    VENV_PYTHON="$VENV_DIR/bin/python3"
elif [ -f "$VENV_DIR/bin/python" ]; then
    VENV_PYTHON="$VENV_DIR/bin/python"
fi

print_green "  All pre-flight checks passed."
echo ""

# Phase 1: Build preprocessing applications
echo "========================================"
echo "PHASE 1: Build preprocessing applications"
echo "========================================"
if [ "$PREPROCESS_ROUTE" = "-p" ]; then
    BUILD_PREPROC_TARGET="preprocessing_tools"
else
    BUILD_PREPROC_TARGET="view3d"
fi
BUILD_PREPROC_CMD="cd '$REPO_ROOT' && tools/build_preprocessing.sh '$SYSTEM' '$BUILD_PREPROC_TARGET'"
run_with_status_strict "$BUILD_PREPROC_CMD" "Building preprocessing applications" "$LOG_FILE"
print_green "[SUCCESS] Preprocessing applications built successfully."
echo ""

echo "========================================"
echo "PHASE 2: Generating input files and comparing with reference data"
echo "========================================"

# Per-case Phase 2 + Phase 3
overall_success=true
case_results=()

for case_num in "${TEST_CASES[@]}"; do
    echo "========================================"
    echo "PROCESSING CASE $case_num"
    echo "========================================"

    # Clean experiment directory — remove generated files, keep source/config files
    find "$SCRIPT_DIR/experiments/$case_num" -maxdepth 1 -type f \
        ! -iname "*.stl"             \
        ! -name "trees.inp.*"        \
        ! -name "heatpump.inp.*"     \
        ! -name "scalarsource*"      \
        ! -name "config*"            \
        ! -name "namoptions.*"       \
        ! -name "timedepnudge.inp.*" \
        ! -name "timedeplw.inp.*"    \
        -delete
    echo ""

    # Step 1: Write inputs. write_inputs.sh runs in the foreground (start mode f)
    # on the current node, so its exit code is the preprocessing's.
    echo "STEP 1: Write inputs for case $case_num"
    echo "=========================================="
    case_dir="$SCRIPT_DIR/experiments/$case_num"
    WRITE_INPUTS_LOG="$case_dir/write_inputs.$case_num.log"
    WRITE_CMD="cd '$REPO_ROOT' && tools/write_inputs.sh '$PREPROCESS_ROUTE' '$case_dir' f '$SYSTEM'"
    set +e
    run_with_status "$WRITE_CMD" "Writing inputs for case $case_num" "$LOG_FILE"
    write_exit_code=$?
    set -e

    if [ ! -f "$WRITE_INPUTS_LOG" ]; then
        print_red "[FAILED] write_inputs log not found: $WRITE_INPUTS_LOG"
        print_red "         write_inputs may not have run at all for case $case_num."
        case_results+=("$case_num: FAILED (write_inputs log missing)")
        overall_success=false
        echo ""
        continue
    fi
    if [ $write_exit_code -ne 0 ] || strip_log "$WRITE_INPUTS_LOG" | grep -qiE "$ERROR_PATTERN"; then
        print_red "[FAILED] write_inputs failed for case $case_num (exit code $write_exit_code):"
        strip_log "$WRITE_INPUTS_LOG" | grep -iE "$ERROR_PATTERN" | head -10 | while read -r line; do
            print_red "  $line"
        done
        echo "Full write_inputs log: $WRITE_INPUTS_LOG"
        case_results+=("$case_num: FAILED (write_inputs errors)")
        overall_success=false
        echo ""
        continue
    fi
    print_green "[SUCCESS] Case $case_num inputs written successfully"
    echo ""

    # Step 2: Compare inputs
    echo "STEP 2: Compare inputs for case $case_num"
    echo "============================================"

    COMPARE_CMD="$VENV_PYTHON $REPO_ROOT/tools/ud_compare_inputs.py $case_num $SCRIPT_DIR/experiments/ $REF_DATA_PATH $TOLERANCE"

    set +e

    CASE_LOG_FILE="${LOG_FILE}.case_${case_num}"
    run_with_status "$COMPARE_CMD" "Comparing input files for case $case_num" "$CASE_LOG_FILE"
    exit_code=$?

    echo "=== INPUT COMPARISON RESULTS FOR CASE $case_num ===" >> "$LOG_FILE"
    cat "$CASE_LOG_FILE" >> "$LOG_FILE"
    echo "=== END CASE $case_num ===" >> "$LOG_FILE"
    echo "" >> "$LOG_FILE"

    set -e

    if [ $exit_code -eq 0 ]; then
        case_results+=("$case_num: PASSED")
        print_green "[SUCCESS] Case $case_num: All input file comparisons PASSED"
    else
        case_results+=("$case_num: FAILED")
        overall_success=false
        print_red "[FAILED] Case $case_num: Some input file comparisons FAILED"
        echo "  See $LOG_FILE for details."
    fi

    rm -f "$CASE_LOG_FILE"
    echo ""
done

# Summary
echo ""
echo "========================================" >> "$LOG_FILE"
echo "Input write & comparison completed at $(date)" >> "$LOG_FILE"

echo "========================================"
echo "  u-dales Input Write & Comparison Summary"
echo "========================================"
echo "  Reference path : $REF_DATA_PATH"
echo "  Cases tested   : ${TEST_CASES[*]}"
echo "  Preprocess via : $PREPROCESS_ROUTE_NAME ($PREPROCESS_ROUTE)"
echo "  Finished       : $(date)"
echo ""
echo "  Individual case results:"
for result in "${case_results[@]}"; do
    case_num=$(echo "$result" | cut -d: -f1)
    status=$(echo "$result" | cut -d: -f2- | xargs)
    if [ "$status" = "PASSED" ]; then
        echo -e "    - $case_num: ${GREEN}PASSED${NC}"
    else
        echo -e "    - $case_num: ${RED}${status}${NC}"
    fi
done
echo ""

if [ "$overall_success" = true ]; then
    print_green "[SUCCESS] ALL INPUT COMPARISONS PASSED"
    echo "Check $LOG_FILE for detailed logs."
    exit 0
else
    print_red "[FAILED] SOME INPUT COMPARISONS FAILED"
    echo "Check $LOG_FILE for detailed logs."
    exit 1
fi
