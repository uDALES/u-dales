#!/usr/bin/env bash

set -e

# Resolve paths relative to this script's location, regardless of where it is called from
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Script usage function
usage() {
    echo "Usage: $0 <ref_data_path> [case1 case2 ...] [--tolerance <val>] [--system <common|icl>]"
    echo ""
    echo "Arguments:"
    echo "  ref_data_path  : Path to the reference data directory (required)"
    echo "  case1 ...      : Optional list of experiment cases (e.g., 100 201 224)"
    echo "                   If not specified, uses default cases: ${DEFAULT_TEST_CASES[*]}"
    echo "  --tolerance    : Max absolute error for numeric files (default: 1e-10)"
    echo "  --system       : Build system type for preprocessing: common or icl (default: common)"
    echo ""
    echo "Phases:"
    echo "  PHASE 1 : Build preprocessing applications (View3D)"
    echo "  PHASE 2 : Run write_inputs.sh for each case to (re-)generate input files"
    echo "  PHASE 3 : Compare generated input files against the reference directory"
    echo ""
    echo "Examples:"
    echo "  $0 ref_data/                             # Run default cases"
    echo "  $0 ref_data/ 990 991                     # Run specific cases"
    echo "  $0 ref_data/ 992 --tolerance 1e-8"
    echo "  $0 ref_data/ 993 --system icl"
    exit 1
}

DEFAULT_TEST_CASES=(990 991 992 993 994 995 997 998 999 807)

# Parse arguments
if [ $# -eq 0 ]; then
    usage
fi

REF_DATA_PATH="$1"
shift

TEST_CASES=()
TOLERANCE="1e-10"
SYSTEM="common"

while [ $# -gt 0 ]; do
    case "$1" in
        --tolerance)
            shift
            TOLERANCE="$1"
            ;;
        --system)
            shift
            SYSTEM="$1"
            ;;
        --*)
            echo "Unknown option: $1"
            usage
            ;;
        *)
            TEST_CASES+=("$1")
            ;;
    esac
    shift
done

# Validate --system value
if [ "$SYSTEM" != "common" ] && [ "$SYSTEM" != "icl" ]; then
    echo "[ERROR] --system must be 'common' or 'icl', got: '$SYSTEM'"
    usage
fi

if [ ${#TEST_CASES[@]} -eq 0 ]; then
    TEST_CASES=("${DEFAULT_TEST_CASES[@]}")
fi

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

print_green()  { echo -e "${GREEN}$1${NC}"; }
print_red()    { echo -e "${RED}$1${NC}"; }
print_yellow() { echo -e "${YELLOW}$1${NC}"; }

# Spinner
show_spinner() {
    local pid=$1
    local message=$2
    local delay=0.1
    local spinstr='|/-\'
    while [ "$(ps a | awk '{print $1}' | grep $pid)" ]; do
        local temp=${spinstr#?}
        printf "\r[%c] %s" "$spinstr" "$message"
        local spinstr=$temp${spinstr%"$temp"}
        sleep $delay
    done
    printf "\r"
}

# Run a command with spinner; exits on failure
# Optional 4th arg: suppress_ok=true skips the "[OK]...Done" message
run_with_status_strict() {
    local cmd="$1"
    local description="$2"
    local log_file="$3"
    local suppress_ok="${4:-false}"

    echo "Starting: $description"

    bash -c "$cmd" >> "$log_file" 2>&1 &
    local pid=$!
    show_spinner $pid "$description"
    wait $pid
    local exit_code=$?

    if [ $exit_code -eq 0 ]; then
        [ "$suppress_ok" != "true" ] && printf "[OK] %s - Done\n" "$description"
        return 0
    else
        printf "[FAILED] %s - Exit Code: %d\n" "$description" "$exit_code"
        print_red "FAILED: $description (exit code: $exit_code)"
        echo "Check $log_file for details."
        exit $exit_code
    fi
}

# Wait for MATLAB launched by write_inputs.sh to finish.
# write_inputs.sh backgrounds MATLAB with: nohup matlab ... -r "expnr=N; write_inputs; quit"
# We detect completion by polling for a running MATLAB process matching that expnr.
# Usage: wait_for_matlab <case_num>
wait_for_matlab() {
    local case_num="$1"
    local elapsed=0

    # Brief initial pause to allow MATLAB to be spawned
    sleep 3

    while pgrep -f "expnr=$case_num" > /dev/null 2>&1; do
        elapsed=$((elapsed + 5))
        printf "\r\033[K  write_inputs running (expnr=%s) — elapsed: %ds" "$case_num" "$elapsed"
        sleep 5
    done
    printf "\r\033[K  write_inputs finished (expnr=%s) after ~%ds\n" "$case_num" "$elapsed"

    return 0
}

# Run a command with spinner; does not exit on failure — returns exit code
run_with_status() {
    local cmd="$1"
    local description="$2"
    local log_file="$3"

    echo "Starting: $description"

    bash -c "$cmd" >> "$log_file" 2>&1 &
    local pid=$!
    show_spinner $pid "$description"
    wait $pid
    local exit_code=$?

    printf "[OK] %s - Done\n" "$description"
    return $exit_code
}

# Log file
LOGDIR="$REPO_ROOT/logdir"
mkdir -p "$LOGDIR"
LOG_FILE="$LOGDIR/test_inputs.log"
echo "u-dales Input Write & Comparison Started at $(date)" > "$LOG_FILE"
echo "========================================" >> "$LOG_FILE"
echo "  Test experiment path : $REF_DATA_PATH" >> "$LOG_FILE"
echo "  Reference path       : $SCRIPT_DIR/experiments/" >> "$LOG_FILE"
echo "  Cases                : ${TEST_CASES[*]}" >> "$LOG_FILE"
echo "  Tolerance            : $TOLERANCE" >> "$LOG_FILE"
echo "  System               : $SYSTEM" >> "$LOG_FILE"
echo "========================================" >> "$LOG_FILE"

echo "Starting u-dales input write & comparison..."
echo ""
echo "========================================"
echo "  Run settings"
echo "========================================"
echo "  Test experiment path : $REF_DATA_PATH"
echo "  Reference path       : $SCRIPT_DIR/experiments/"
echo "  Cases                : ${TEST_CASES[*]}"
echo "  Tolerance            : $TOLERANCE"
echo "  System               : $SYSTEM"
echo "  Log file             : $LOG_FILE"
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

if [ ! -d "$REPO_ROOT/tools/.venv_netcdf" ]; then
    print_red "[MISSING] NetCDF virtual environment not found: $REPO_ROOT/tools/.venv_netcdf"
    print_red "         Run '$REPO_ROOT/tools/ud_set_nc_venv.sh' to set it up."
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
if [ -f "$REPO_ROOT/tools/.venv_netcdf/bin/python3" ]; then
    VENV_PYTHON="$REPO_ROOT/tools/.venv_netcdf/bin/python3"
elif [ -f "$REPO_ROOT/tools/.venv_netcdf/bin/python" ]; then
    VENV_PYTHON="$REPO_ROOT/tools/.venv_netcdf/bin/python"
fi

print_green "  All pre-flight checks passed."
echo ""

# Phase 1: Build preprocessing applications
echo "========================================"
echo "PHASE 1: Build preprocessing applications"
echo "========================================" 
BUILD_PREPROC_CMD="cd '$REPO_ROOT' && tools/build_preprocessing.sh '$SYSTEM'"
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
        -delete
    echo ""

    # Step 1: Write inputs
    echo "STEP 1: Write inputs for case $case_num"
    echo "=========================================="
    WRITE_CMD="cd '$REPO_ROOT' && tools/write_inputs.sh '$SCRIPT_DIR/experiments/$case_num'"
    run_with_status_strict "$WRITE_CMD" "Writing inputs for case $case_num" "$LOG_FILE" true

    # write_inputs.sh launches MATLAB in the background and exits immediately.
    # Poll until the MATLAB process for this case finishes before comparing.
    wait_for_matlab "$case_num"

    # Check MATLAB log for errors
    MATLAB_LOG="$SCRIPT_DIR/experiments/$case_num/write_inputs.$case_num.log"
    if [ ! -f "$MATLAB_LOG" ]; then
        print_red "[FAILED] write_inputs log not found: $MATLAB_LOG"
        print_red "         write_inputs may not have run at all for case $case_num."
        case_results+=("$case_num: FAILED (write_inputs log missing)")
        overall_success=false
        echo ""
        continue
    fi
    if col -b < "$MATLAB_LOG" | grep -qiE "^Error |^Error:|Undefined function|Undefined variable|not enough input|out of memory"; then
        print_red "[FAILED] write_inputs reported errors for case $case_num:"
        col -b < "$MATLAB_LOG" | grep -iE "^Error |^Error:|Undefined function|Undefined variable|not enough input|out of memory" | head -10 | while read -r line; do
            print_red "  $line"
        done
        echo "Full write_inputs log: $MATLAB_LOG"
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
