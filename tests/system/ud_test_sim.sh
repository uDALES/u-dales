#!/usr/bin/env bash

set -e

# Resolve paths relative to this script's location, regardless of where it is called from
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Script usage function
usage() {
    echo "Usage: $0 <ref_data_path> [case1 case2 ...] [--tolerance <val>] [--tol-thl <val>] [--system <common|gpu>]"
    echo ""
    echo "Arguments:"
    echo "  ref_data_path  : Path to the reference data directory (required)"
    echo "  case1 ...      : Optional list of test cases to run (e.g., 100 201 402)"
    echo "                   If not specified, uses default cases for the system"
    echo "  --tolerance    : Max absolute error for NetCDF output comparisons (default: 1e-6)"
    echo "  --tol-thl      : Tolerance for temperature variables (default: same as --tolerance)"
    echo "  --system       : Build system type: common (CPU) or gpu (default: common)"
    echo ""
    echo "Examples:"
    echo "  $0 ref_data/                                    # Run default CPU test cases"
    echo "  $0 ref_data/ --system common                   # Run default CPU test cases"
    echo "  $0 ref_data/ --system gpu                      # Run default GPU test cases"
    echo "  $0 ref_data/ 100 201                           # Run specific CPU test cases"
    echo "  $0 ref_data/ 402 502 452 --system gpu          # Run specific GPU test cases"
    echo "  $0 ref_data/ 224 --tolerance 1e-8"
    echo "  $0 ref_data/ 224 --tolerance 1e-8 --tol-thl 1e-7"
    exit 1
}

# Parse command line arguments
if [ $# -eq 0 ]; then
    usage
fi

REF_DATA_PATH="$1"
shift

# Remaining args: test cases and optional flags
TEST_CASES=()
TOLERANCE="1e-6"
TOL_THL=""
SYSTEM="common"

while [ $# -gt 0 ]; do
    case "$1" in
        --tolerance)
            shift
            TOLERANCE="$1"
            ;;
        --tol-thl)
            shift
            TOL_THL="$1"
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

# Default tol_thl to tolerance if not explicitly set
if [ -z "$TOL_THL" ]; then
    TOL_THL="$TOLERANCE"
fi

# Validate --system value
if [ "$SYSTEM" != "common" ] && [ "$SYSTEM" != "gpu" ]; then
    echo "[ERROR] --system must be 'common' or 'gpu', got: '$SYSTEM'"
    usage
fi

# Set default test cases based on system
if [ "$SYSTEM" = "common" ]; then
    DEFAULT_TEST_CASES=(100 218 224 242 807)
elif [ "$SYSTEM" = "gpu" ]; then
    DEFAULT_TEST_CASES=(402 502 452 410 411)
else
    echo "[ERROR] Invalid --system value: '$SYSTEM'. Must be 'common' or 'gpu'."
    exit 1
fi

# Use provided test cases or defaults
if [ ${#TEST_CASES[@]} -eq 0 ]; then
    TEST_CASES=("${DEFAULT_TEST_CASES[@]}")
fi

# Color codes for terminal output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored text
print_green() {
    echo -e "${GREEN}$1${NC}"
}

print_red() {
    echo -e "${RED}$1${NC}"
}

print_yellow() {
    echo -e "${YELLOW}$1${NC}"
}

# Function to show a spinner while command is running
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
    printf "\r"  # Clear the spinner line
}

# Function to run a command with spinner and status messages
run_with_status() {
    local cmd="$1"
    local description="$2"
    local success_msg="$3"
    local log_file="$4"
    
    echo "Starting: $description"
    
    # Run command in background and capture output to log
    bash -c "$cmd" >> "$log_file" 2>&1 &
    local pid=$!
    
    # Show spinner
    show_spinner $pid "$description"
    
    # Wait for completion and check exit code
    wait $pid
    local exit_code=$?
    
    if [ $exit_code -eq 0 ]; then
        printf "[OK] %s - Done\n" "$description"
        if [ -n "$success_msg" ]; then
            print_green "$success_msg"
        fi
        return 0
    else
        printf "[FAILED] %s - Exit Code: %d\n" "$description" "$exit_code"
        print_red "FAILED: $description (exit code: $exit_code)"
        echo "Check $log_file for details."
        exit $exit_code
    fi
}

# Function to run comparison tests without exiting on failure
run_comparison_with_status() {
    local cmd="$1"
    local description="$2"
    local success_msg="$3"
    local log_file="$4"
    
    echo "Starting: $description"
    
    # Run command in background and capture output to log
    bash -c "$cmd" >> "$log_file" 2>&1 &
    local pid=$!
    
    # Show spinner
    show_spinner $pid "$description"
    
    # Wait for completion and check exit code
    wait $pid
    local exit_code=$?
    
    # Always show OK since the script executed (regardless of test results)
    printf "[OK] %s - Done\n" "$description"
    
    if [ -n "$success_msg" ] && [ $exit_code -eq 0 ]; then
        print_green "$success_msg"
    fi
    
    return $exit_code
}

# Log file
LOGDIR="$REPO_ROOT/logdir"
mkdir -p "$LOGDIR"
LOG_FILE="$LOGDIR/test_${SYSTEM}.log"
echo "u-dales Simulation Test Suite Started at $(date)" > "$LOG_FILE"
echo "========================================" >> "$LOG_FILE"
echo "  Experiments path : $SCRIPT_DIR/experiments/" >> "$LOG_FILE"
echo "  Reference path   : $REF_DATA_PATH" >> "$LOG_FILE"
echo "  Cases            : ${TEST_CASES[*]}" >> "$LOG_FILE"
echo "  Tolerance        : $TOLERANCE" >> "$LOG_FILE"
echo "  Tol THL          : $TOL_THL" >> "$LOG_FILE"
echo "  System           : $SYSTEM" >> "$LOG_FILE"
echo "========================================" >> "$LOG_FILE"

echo "Starting u-dales simulation test suite..."
echo ""
echo "========================================"
echo "  Run settings"
echo "========================================"
echo "  Experiments path : $SCRIPT_DIR/experiments/"
echo "  Reference path   : $REF_DATA_PATH"
echo "  Cases            : ${TEST_CASES[*]}"
echo "  Tolerance        : $TOLERANCE"echo "  Tol THL          : $TOL_THL"echo "  System           : $SYSTEM"
echo "  Log file         : $LOG_FILE"
echo "========================================"
echo ""

# Pre-flight checks: verify all required directories and files exist before starting
echo "========================================"
echo "  Pre-flight checks"
echo "========================================"

preflight_ok=true

# Required repo-level scripts
for required_script in \
    "$REPO_ROOT/tools/build_executable.sh" \
    "$REPO_ROOT/tools/local_execute.sh"
do
    if [ ! -f "$required_script" ]; then
        print_red "[MISSING] Required script not found: $required_script"
        preflight_ok=false
    fi
done

# ud_compare_outputs.py must be present in tools/
if [ ! -f "$REPO_ROOT/tools/ud_compare_outputs.py" ]; then
    print_red "[MISSING] ud_compare_outputs.py not found in: $REPO_ROOT/tools"
    preflight_ok=false
fi

VENV_DIR="$REPO_ROOT/tools/python/.venv"

# Python virtual environment must exist for the comparison phase
if [ ! -d "$VENV_DIR" ]; then
    print_red "[MISSING] Python virtual environment not found: $VENV_DIR"
    print_red "         Run '$REPO_ROOT/tools/python/setup_venv.sh' to set it up."
    preflight_ok=false
elif [ ! -f "$VENV_DIR/bin/python3" ] && [ ! -f "$VENV_DIR/bin/python" ]; then
    print_red "[MISSING] Python interpreter not found inside virtual environment: $VENV_DIR"
    preflight_ok=false
fi

# Each requested test case must have an experiments input directory and a ref_data directory
for case_num in "${TEST_CASES[@]}"; do
    if [ ! -d "$SCRIPT_DIR/experiments/$case_num" ]; then
        print_red "[MISSING] Experiment input directory not found: $SCRIPT_DIR/experiments/$case_num"
        preflight_ok=false
    fi

    if [ ! -d "$REF_DATA_PATH/$case_num" ]; then
        print_red "[MISSING] Reference data directory not found: $REF_DATA_PATH/$case_num"
        preflight_ok=false
    fi
done

if [ "$preflight_ok" = false ]; then
    echo ""
    print_red "Pre-flight checks FAILED. Resolve the issues above before running the test suite."
    exit 1
fi

# Resolve the Python interpreter inside the venv (prefer python3, fall back to python)
if [ -f "$VENV_DIR/bin/python3" ]; then
    VENV_PYTHON="$VENV_DIR/bin/python3"
elif [ -f "$VENV_DIR/bin/python" ]; then
    VENV_PYTHON="$VENV_DIR/bin/python"
fi

print_green "  All pre-flight checks passed."
echo ""

# Phase 1: Build executable (only once)
echo "PHASE 1: Build u-dales executable ($SYSTEM)"
echo "======================================"

# Set build command based on system
if [ "$SYSTEM" = "common" ]; then
    BUILD_CMD="cd '$REPO_ROOT' && echo 'Building from: $REPO_ROOT' && rm -rf build && echo 'REMOVED old u-dales build directory.' && tools/build_executable.sh common release"
elif [ "$SYSTEM" = "gpu" ]; then
    BUILD_CMD="cd '$REPO_ROOT' && echo 'Building from: $REPO_ROOT' && rm -rf build && echo 'REMOVED old u-dales build directory.' && module purge && module use /opt/nvidia/hpc_sdk/modulefiles && module load nvhpc/24.11 && module list && tools/build_executable.sh gpu release"
else
    echo "[ERROR] Invalid --system value: '$SYSTEM'. Must be 'common' or 'gpu'."
    exit 1
fi

run_with_status "$BUILD_CMD" "Building u-dales executable" "[SUCCESS] Executable built successfully" "$LOG_FILE"
echo ""

# Phase 2 & 3: Run simulation and comparison for each test case
overall_success=true
case_results=()

for case_num in "${TEST_CASES[@]}"; do
    echo "========================================"
    echo "PROCESSING CASE $case_num"
    echo "========================================"
    
    # Phase 2: Run simulation
    echo "PHASE 2: Simulation for case $case_num"
    echo "============================================="
    # Inline simulation execution
    SIM_CMD="cd '$REPO_ROOT' && rm -rf '$SCRIPT_DIR/outputs/$case_num' && echo 'REMOVED old outputs/$case_num.' && tools/local_execute.sh '$SCRIPT_DIR/experiments/$case_num'"
    run_with_status "$SIM_CMD" "Running simulation for case $case_num" "" "$LOG_FILE"
    simulation_exit_code=$?
    
    # Check simulation results
    if [ $simulation_exit_code -eq 0 ]; then
        print_green "[SUCCESS] Case $case_num simulation completed successfully"
    else
        print_red "[FAILED] Case $case_num simulation failed"
        case_results+=("$case_num: FAILED")
        overall_success=false
        echo ""
        continue  # Skip comparison if simulation failed
    fi
    echo ""
    
    # Phase 3: Compare results
    echo "PHASE 3: Compare results for case $case_num"
    echo "==========================================="

    COMPARE_CMD="$VENV_PYTHON $REPO_ROOT/tools/ud_compare_outputs.py $case_num $SCRIPT_DIR/outputs/ $REF_DATA_PATH $TOLERANCE $TOL_THL"
    
    # Disable exit-on-error for comparison (we want to handle failure gracefully)
    set +e
    
    # Create a temporary log file for this specific case comparison
    CASE_LOG_FILE="${LOG_FILE}.case_${case_num}"
    
    # Use the existing run_comparison_with_status function
    run_comparison_with_status "$COMPARE_CMD" "Running comparison tests for case $case_num" "" "$CASE_LOG_FILE"
    comparison_exit_code=$?
    
    # Append the case-specific log to the main log file
    echo "=== COMPARISON RESULTS FOR CASE $case_num ===" >> "$LOG_FILE"
    cat "$CASE_LOG_FILE" >> "$LOG_FILE"
    echo "=== END CASE $case_num ===" >> "$LOG_FILE"
    echo "" >> "$LOG_FILE"
    
    set -e  # Re-enable exit-on-error
    
    # Track results for this case based on exit code
    if [ $comparison_exit_code -eq 0 ]; then
        case_results+=("$case_num: PASSED")
        print_green "[SUCCESS] Case $case_num: All comparison tests PASSED"
    else
        case_results+=("$case_num: FAILED") 
        overall_success=false
        
        # Check if this was a script execution failure or test result failure
        # Use the case-specific log file for accurate detection
        if grep -q "=== SUMMARY ===" "$CASE_LOG_FILE"; then
            print_yellow "[MIXED] Case $case_num: Comparison completed, but some tests FAILED tolerance"
        else
            print_red "[FAILED] Case $case_num: Comparison script failed to execute properly"
        fi

        # Phase 4: Compare input files to help diagnose why outputs differ
        echo ""
        echo "PHASE 4: Compare input files for case $case_num (triggered by Phase 3 failure)"
        echo "==============================================================================="

        COMPARE_INPUTS_CMD="$VENV_PYTHON $REPO_ROOT/tools/ud_compare_inputs.py $case_num $SCRIPT_DIR/experiments/ $REF_DATA_PATH"

        set +e
        CASE_INPUTS_LOG_FILE="${LOG_FILE}.inputs_${case_num}"
        run_comparison_with_status "$COMPARE_INPUTS_CMD" "Comparing input files for case $case_num" "" "$CASE_INPUTS_LOG_FILE"
        inputs_exit_code=$?

        echo "=== INPUT COMPARISON RESULTS FOR CASE $case_num ===" >> "$LOG_FILE"
        cat "$CASE_INPUTS_LOG_FILE" >> "$LOG_FILE"
        echo "=== END INPUT COMPARISON CASE $case_num ===" >> "$LOG_FILE"
        echo "" >> "$LOG_FILE"

        set -e

        if [ $inputs_exit_code -eq 0 ]; then
            print_green "[INFO] Case $case_num: Input files match reference — output differences are not due to input changes"
        else
            print_yellow "[INFO] Case $case_num: Input files differ from reference — this may explain the output differences"
        fi

        rm -f "$CASE_INPUTS_LOG_FILE"
    fi
    
    # Clean up the temporary case log file
    rm -f "$CASE_LOG_FILE"
    echo ""
done

# Final summary
echo ""
echo "========================================" >> $LOG_FILE
echo "Test suite completed at $(date)" >> $LOG_FILE

echo "========================================"
echo "u-dales Test Suite Summary"
echo "========================================"
echo "System: $SYSTEM"
echo "Cases tested: ${TEST_CASES[*]}"
echo ""
echo "Individual case results:"
for result in "${case_results[@]}"; do
    case_num=$(echo "$result" | cut -d: -f1)
    status=$(echo "$result" | cut -d: -f2- | xargs)
    if [ "$status" = "PASSED" ]; then
        echo -e "  - $case_num: ${GREEN}PASSED${NC}"
    elif [[ "$status" == SKIPPED* ]]; then
        echo -e "  - $case_num: ${YELLOW}${status}${NC}"
    else
        echo -e "  - $case_num: ${RED}${status}${NC}"
    fi
done
echo ""

if [ "$overall_success" = true ]; then
    print_green "[SUCCESS] ALL PHASES COMPLETED SUCCESSFULLY"
    echo "  - Build: SUCCESS"
    echo "  - All simulations: SUCCESS"
    echo "  - All comparisons: SUCCESS"
    echo ""
    echo "Check $LOG_FILE for detailed logs."
    exit 0
else
    print_red "[FAILED] TEST SUITE FAILED"
    echo "  - Build: SUCCESS"
    echo "  - Some simulations/comparisons: FAILED"
    echo ""
    echo "Check $LOG_FILE for detailed logs."
    exit 1
fi
