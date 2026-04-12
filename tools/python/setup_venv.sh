#!/bin/bash
# Setup script for Linux/WSL.
# Creates a Python virtual environment, installs dependencies, and builds
# the required preprocessing artifacts (View3D and f2py extension modules)
# through the standalone preprocessing CMake entry point.
#
# Usage (run from the repository root):
#   bash tools/python/setup_venv.sh [build_system] [build_target]
#
# Arguments:
#   build_system   Build environment to use.
#                  Allowed values : common | icl
#                  Default        : common
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
#                    PYTHON_BIN=/opt/pbs/python/bin/python3 bash tools/python/setup_venv.sh
#
#   VENV_DIR       Explicit path for the virtual environment directory.
#                  Default: tools/python/.venv (relative to the repository root)
#                  Falls back to .venv (repo root) if that legacy directory exists.
#
# Virtual environment behaviour:
#   - If the target venv directory does not exist it is created from scratch.
#   - If it already exists the script prompts whether to recreate it (default: N).
#     Answering N exits immediately without reinstalling packages or rebuilding.
#
# Examples:
#   # Default local setup
#   bash tools/python/setup_venv.sh
#
#   # Local setup, View3D only (no f2py modules)
#   bash tools/python/setup_venv.sh common view3d
#
#   # ICL cluster, full build
#   bash tools/python/setup_venv.sh icl preprocessing_tools
#
#   # Use a specific Python interpreter
#   PYTHON_BIN=/opt/pbs/python/bin/python3 bash tools/python/setup_venv.sh

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UDALES_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
BUILD_SYSTEM="${1:-common}"
BUILD_TARGET="${2:-preprocessing_tools}"

if [[ "$BUILD_SYSTEM" != "common" && "$BUILD_SYSTEM" != "icl" ]]; then
    echo "Error: invalid build system '${BUILD_SYSTEM}'."
    echo ""
    echo "Usage: bash tools/python/setup_venv.sh [build_system] [build_target]"
    echo ""
    echo "  build_system   Build environment to use (default: common)"
    echo "                   common  - local Linux / WSL system"
    echo "                   icl     - Imperial College London HPC cluster"
    echo ""
    echo "  build_target   CMake target to build (default: preprocessing_tools)"
    echo "                   view3d               - View3D executable only"
    echo "                   preprocessing_tools  - View3D + f2py extension modules"
    exit 1
fi

if [[ "$BUILD_TARGET" != "view3d" && "$BUILD_TARGET" != "preprocessing_tools" ]]; then
    echo "Error: invalid build target '${BUILD_TARGET}'."
    echo ""
    echo "Usage: bash tools/python/setup_venv.sh [build_system] [build_target]"
    echo ""
    echo "  build_system   Build environment to use (default: common)"
    echo "                   common  - local Linux / WSL system"
    echo "                   icl     - Imperial College London HPC cluster"
    echo ""
    echo "  build_target   CMake target to build (default: preprocessing_tools)"
    echo "                   view3d               - View3D executable only"
    echo "                   preprocessing_tools  - View3D + f2py extension modules"
    exit 1
fi

DEFAULT_VENV_DIR="${UDALES_ROOT}/tools/python/.venv"
LEGACY_VENV_DIR="${UDALES_ROOT}/.venv"

if [ -n "${VENV_DIR:-}" ]; then
    VENV_DIR="$VENV_DIR"
elif [ -d "$DEFAULT_VENV_DIR" ]; then
    VENV_DIR="$DEFAULT_VENV_DIR"
elif [ -d "$LEGACY_VENV_DIR" ]; then
    VENV_DIR="$LEGACY_VENV_DIR"
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
if [ "$VENV_DIR" = "$LEGACY_VENV_DIR" ]; then
    echo "Note: using legacy repo-local virtual environment. New setups default to $DEFAULT_VENV_DIR"
fi
echo ""

if ! command -v "$PYTHON_BIN" &> /dev/null; then
    echo "Error: $PYTHON_BIN is not installed"
    exit 1
fi

echo "Python version: $($PYTHON_BIN --version)"

PYTHON_INCLUDE="$($PYTHON_BIN -c 'import sysconfig; print(sysconfig.get_paths().get("include", ""))')"
if [ ! -f "$PYTHON_INCLUDE/Python.h" ]; then
    echo "Error: Python development headers were not found for $PYTHON_BIN"
    echo "Expected: $PYTHON_INCLUDE/Python.h"
    echo "The directshortwave f2py wrapper is mandatory."
    echo "Set PYTHON_BIN to an interpreter with headers, e.g."
    echo "  PYTHON_BIN=/opt/pbs/python/bin/python3 ./tools/python/setup_venv.sh"
    exit 1
fi

# Create virtual environment if it doesn't exist
if [ -d "$VENV_DIR" ]; then
    echo "Virtual environment already exists at: $VENV_DIR"
    read -p "Do you want to recreate it? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing virtual environment..."
        rm -rf "$VENV_DIR"
    else
        echo "Using existing virtual environment."
        echo "To activate: source \"$VENV_DIR/bin/activate\""
        exit 0
    fi
fi

echo "Creating virtual environment at: $VENV_DIR"
"$PYTHON_BIN" -m venv "$VENV_DIR"

echo "Activating virtual environment..."
source "$VENV_DIR/bin/activate"

echo "Upgrading pip..."
pip install --upgrade pip

echo "Installing dependencies from requirements.txt..."
if [ -f "${SCRIPT_DIR}/requirements.txt" ]; then
    pip install -r "${SCRIPT_DIR}/requirements.txt"
else
    echo "Warning: requirements.txt not found at ${SCRIPT_DIR}/requirements.txt"
fi

echo "Installing build dependencies from requirements-build.txt..."
if [ -f "${SCRIPT_DIR}/requirements-build.txt" ]; then
    pip install -r "${SCRIPT_DIR}/requirements-build.txt"
else
    echo "Warning: requirements-build.txt not found at ${SCRIPT_DIR}/requirements-build.txt"
fi

echo "Building preprocessing tools..."
PREPROCESSING_PYTHON_EXECUTABLE="$(command -v python)" \
    "${UDALES_ROOT}/tools/build_preprocessing.sh" "${BUILD_SYSTEM}" "${BUILD_TARGET}"

echo ""
echo "=========================================="
echo "Setup complete!"
echo "=========================================="
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
