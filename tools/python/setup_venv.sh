#!/bin/bash
# Setup script for Linux/WSL.
# Creates a Python virtual environment, installs dependencies, and builds
# the required preprocessing artifacts (View3D and f2py extension modules)
# through the standalone preprocessing CMake entry point.
#
# Usage (run from the repository root):
#   bash tools/python/setup_venv.sh <build_system> [build_target]
#
# Arguments:
#   build_system   Build environment to use. Required.
#                  Allowed values : common | icl
#                    common - local Linux / WSL system (no module loading)
#                    icl    - Imperial College London HPC cluster
#                             (loads CMake and Python environment modules)
#
#   build_target   CMake preprocessing target to build.
#                  Allowed values : view3d | preprocessing_tools
#                  Default        : preprocessing_tools
#                    view3d               - View3D executable only
#                    preprocessing_tools  - View3D + directshortwave and
#                                          IBM f2py extension modules
#
# Environment variables (optional overrides):
#   PYTHON_BIN     Python interpreter to use for creating the venv.
#                  Default: python3
#                  Override when the system python3 lacks development headers:
#                    PYTHON_BIN=/opt/pbs/python/bin/python3 bash tools/python/setup_venv.sh common
#
#   VENV_DIR       Explicit path for the virtual environment directory.
#                  Default: tools/python/.venv (relative to the repository root)
#                  Falls back to .venv (repo root) if that legacy directory exists.
#
# Virtual environment behaviour:
#   - If the target venv directory does not exist it is created from scratch.
#   - If it already exists the script prompts whether to recreate it (default: N).
#     Answering N skips installation/rebuild and validates the existing venv.
#
# Examples:
#   # Local setup
#   bash tools/python/setup_venv.sh common
#
#   # Local setup, View3D only (no f2py modules)
#   bash tools/python/setup_venv.sh common view3d
#
#   # ICL cluster, full build
#   bash tools/python/setup_venv.sh icl preprocessing_tools
#
#   # Use a specific Python interpreter
#   PYTHON_BIN=/opt/pbs/python/bin/python3 bash tools/python/setup_venv.sh common

set -euo pipefail

if [ -t 1 ]; then
    RED='\033[0;31m'
    GREEN='\033[0;32m'
    YELLOW='\033[1;33m'
    NC='\033[0m'
else
    RED=''
    GREEN=''
    YELLOW=''
    NC=''
fi

print_green() { echo -e "${GREEN}$1${NC}"; }
print_red() { echo -e "${RED}$1${NC}"; }
print_yellow() { echo -e "${YELLOW}$1${NC}"; }

usage() {
    echo "Usage: bash tools/python/setup_venv.sh <build_system> [build_target]"
    echo ""
    echo "  build_system   Build environment to use (required)"
    echo "                   common  - local Linux / WSL system"
    echo "                   icl     - Imperial College London HPC cluster"
    echo ""
    echo "  build_target   CMake target to build (default: preprocessing_tools)"
    echo "                   view3d               - View3D executable only"
    echo "                   preprocessing_tools  - View3D + f2py extension modules"
}

die() {
    print_red "[ERROR] $1"
    exit 1
}

check_python_version() {
    local python_cmd="$1"
    if ! "$python_cmd" -c "import sys; sys.exit(0 if sys.version_info >= (3, 9) else 1)" >/dev/null 2>&1; then
        die "Python 3.9 or newer is required (found $("$python_cmd" --version 2>&1))."
    fi
}

check_python_headers() {
    local python_cmd="$1"
    local python_include
    python_include="$("$python_cmd" -c 'import sysconfig; print(sysconfig.get_paths().get("include", ""))')"
    if [ ! -f "$python_include/Python.h" ]; then
        print_red "[ERROR] Python development headers were not found for $python_cmd"
        echo "Expected: $python_include/Python.h"
        echo "The directshortwave f2py wrapper is mandatory."
        echo "Set PYTHON_BIN to an interpreter with headers, e.g."
        echo "  PYTHON_BIN=/opt/pbs/python/bin/python3 ./tools/python/setup_venv.sh common"
        exit 1
    fi
}

load_build_environment() {
    if [ "$BUILD_SYSTEM" = "icl" ]; then
        print_yellow "Loading ICL environment modules..."
        if ! command -v module >/dev/null 2>&1; then
            die "Environment modules are not available; cannot use build system 'icl'."
        fi
        module load CMake/3.31.3-GCCcore-14.2.0
        module load Python/3.13.1-GCCcore-14.2.0
    fi
}

run_import_checks() {
    print_yellow "Testing Python package imports..."
    PYTHONPATH="${SCRIPT_DIR}${PYTHONPATH:+:$PYTHONPATH}" python - <<'PY'
import importlib
import sys

modules = [
    ("numpy", "numpy"),
    ("scipy", "scipy"),
    ("xarray", "xarray"),
    ("netCDF4", "netCDF4"),
    ("h5netcdf", "h5netcdf"),
    ("trimesh", "trimesh"),
    ("shapely", "shapely"),
    ("triangle", "triangle"),
    ("matplotlib", "matplotlib"),
    ("f90nml", "f90nml"),
    ("numba", "numba"),
    ("pvlib", "pvlib"),
    # PyVista is the default 3-D backend and a core dependency; Plotly is an
    # optional alternative (requirements-plotly.txt) and is deliberately not
    # required here.
    ("pyvista", "pyvista"),
    # uDALES tools themselves — verifies the editable install resolves without
    # any PYTHONPATH / sys.path manipulation.
    ("udbase", "udbase"),
    ("udgeom", "udgeom"),
    ("udprep", "udprep"),
    ("udvis", "udvis"),
]

failed = False
for label, module_name in modules:
    try:
        module = importlib.import_module(module_name)
        version = getattr(module, "__version__", "available")
        print(f"[OK] {label} {version}")
    except Exception as exc:
        print(f"[ERROR] {label} import failed: {exc}")
        failed = True

sys.exit(1 if failed else 0)
PY
}

run_f2py_checks() {
    if [ "$BUILD_TARGET" != "preprocessing_tools" ]; then
        return
    fi

    print_yellow "Testing f2py preprocessing modules..."
    local missing=0
    if ! compgen -G "${SCRIPT_DIR}/udprep/directshortwave_f2py"*.so >/dev/null; then
        print_red "[ERROR] directshortwave_f2py module missing in ${SCRIPT_DIR}/udprep"
        missing=1
    fi
    if ! compgen -G "${SCRIPT_DIR}/udprep/ibm_preproc_f2py"*.so >/dev/null; then
        print_red "[ERROR] ibm_preproc_f2py module missing in ${SCRIPT_DIR}/udprep"
        missing=1
    fi
    if [ "$missing" -ne 0 ]; then
        echo "Recreate the virtual environment or run:"
        echo "  PREPROCESSING_PYTHON_EXECUTABLE=\"$VENV_DIR/bin/python\" tools/build_preprocessing.sh $BUILD_SYSTEM preprocessing_tools"
        exit 1
    fi

    PYTHONPATH="${SCRIPT_DIR}${PYTHONPATH:+:$PYTHONPATH}" python - <<'PY'
import importlib
import sys

modules = [
    "udprep.directshortwave_f2py",
    "udprep.ibm_preproc_f2py",
]

failed = False
for module_name in modules:
    try:
        importlib.import_module(module_name)
        print(f"[OK] {module_name} imported successfully")
    except Exception as exc:
        print(f"[ERROR] {module_name} import failed: {exc}")
        failed = True

sys.exit(1 if failed else 0)
PY
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UDALES_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
BUILD_SYSTEM="${1:-}"
BUILD_TARGET="${2:-preprocessing_tools}"

if [ -z "$BUILD_SYSTEM" ]; then
    print_red "[ERROR] missing required build system."
    echo ""
    usage
    exit 1
fi

if [[ "$BUILD_SYSTEM" != "common" && "$BUILD_SYSTEM" != "icl" ]]; then
    print_red "[ERROR] invalid build system '${BUILD_SYSTEM}'."
    echo ""
    usage
    exit 1
fi

if [[ "$BUILD_TARGET" != "view3d" && "$BUILD_TARGET" != "preprocessing_tools" ]]; then
    print_red "[ERROR] invalid build target '${BUILD_TARGET}'."
    echo ""
    usage
    exit 1
fi

load_build_environment

DEFAULT_VENV_DIR="${UDALES_ROOT}/tools/python/.venv"

if [ -n "${VENV_DIR:-}" ]; then
    VENV_DIR="$VENV_DIR"
elif [ -d "$DEFAULT_VENV_DIR" ]; then
    VENV_DIR="$DEFAULT_VENV_DIR"
else
    VENV_DIR="$DEFAULT_VENV_DIR"
fi

echo "=========================================="
echo "Setting up Python virtual environment"
echo "=========================================="
echo "u-dales repository: $UDALES_ROOT"
echo "Python interpreter: $PYTHON_BIN"
echo "Virtual environment: $VENV_DIR"
echo "Build system:        $BUILD_SYSTEM"
echo "Build target:        $BUILD_TARGET"
echo ""

if ! command -v "$PYTHON_BIN" &> /dev/null; then
    die "$PYTHON_BIN is not installed"
fi

echo "Python version: $($PYTHON_BIN --version)"
check_python_version "$PYTHON_BIN"

# Create virtual environment if it doesn't exist
SKIP_INSTALL=false
if [ -d "$VENV_DIR" ]; then
    print_yellow "Virtual environment already exists at: $VENV_DIR"
    if ! read -p "Do you want to recreate it? (y/N): " -n 1 -r; then
        REPLY=""
    fi
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        print_yellow "Removing existing virtual environment..."
        rm -rf "$VENV_DIR"
    else
        print_yellow "Using existing virtual environment and running validation checks."
        SKIP_INSTALL=true
    fi
fi

if [ "$SKIP_INSTALL" = false ]; then
    check_python_headers "$PYTHON_BIN"

    echo "Creating virtual environment at: $VENV_DIR"
    "$PYTHON_BIN" -m venv "$VENV_DIR"
fi

if [ ! -f "$VENV_DIR/bin/activate" ]; then
    die "Virtual environment activation script not found: $VENV_DIR/bin/activate"
fi

echo "Activating virtual environment..."
source "$VENV_DIR/bin/activate"
check_python_version python

if [ -z "${MPLCONFIGDIR:-}" ]; then
    export MPLCONFIGDIR="${TMPDIR:-/tmp}/u-dales-matplotlib-${USER:-user}"
    mkdir -p "$MPLCONFIGDIR"
fi

if [ "$SKIP_INSTALL" = false ]; then
    echo "Upgrading pip..."
    python -m pip install --upgrade pip

    echo "Installing dependencies from requirements.txt..."
    if [ -f "${SCRIPT_DIR}/requirements.txt" ]; then
        python -m pip install -r "${SCRIPT_DIR}/requirements.txt"
    else
        print_yellow "[WARN] requirements.txt not found at ${SCRIPT_DIR}/requirements.txt"
    fi

    echo "Installing build dependencies from requirements-build.txt..."
    if [ -f "${SCRIPT_DIR}/requirements-build.txt" ]; then
        python -m pip install -r "${SCRIPT_DIR}/requirements-build.txt"
    else
        print_yellow "[WARN] requirements-build.txt not found at ${SCRIPT_DIR}/requirements-build.txt"
    fi

    # Register the uDALES Python tools as an editable install so `udbase`,
    # `udgeom`, `udprep`, `udvis` (and the flat helper modules) import from any
    # working directory without PYTHONPATH / sys.path edits. This does not copy
    # or build anything, so the in-place f2py extensions in udprep/ keep working.
    echo "Installing uDALES Python tools (editable: pip install -e)..."
    # --no-deps: the runtime dependencies were already installed from
    # requirements.txt above; this step only registers the source dir on the path.
    python -m pip install -e "${SCRIPT_DIR}" --no-deps
else
    echo "Installation steps skipped."
fi

run_import_checks

echo "Building preprocessing tools..."
PREPROCESSING_PYTHON_EXECUTABLE="$(command -v python)" \
    "${UDALES_ROOT}/tools/build_preprocessing.sh" "${BUILD_SYSTEM}" "${BUILD_TARGET}"

run_f2py_checks

echo ""
print_green "=========================================="
print_green "Setup complete!"
print_green "=========================================="
echo ""
echo "To use the virtual environment:"
echo "  1. Activate:   source $VENV_DIR/bin/activate"
echo "  2. Run script: python tools/write_inputs.py [config_dir]"
echo "  3. Deactivate: deactivate"
echo ""
echo "Example workflow:"
echo "  cd $UDALES_ROOT"
echo "  source $VENV_DIR/bin/activate"
echo "  python tools/write_inputs.py"
echo "  deactivate"
echo ""
echo "System tests and comparison scripts also use this environment:"
echo "  $VENV_DIR"
echo ""
