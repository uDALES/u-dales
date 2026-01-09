#!/bin/bash
# Setup script for Linux/WSL
# Creates a Python virtual environment and installs all dependencies

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Go to u-dales root (tools/python -> u-dales)
UDALES_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
# Create venv in parent directory of u-dales
VENV_DIR="$(cd "${UDALES_ROOT}/.." && pwd)/venv-udales"

echo "=========================================="
echo "Setting up Python virtual environment"
echo "=========================================="
echo "u-dales repository: $UDALES_ROOT"
echo "Virtual environment: $VENV_DIR"
echo ""

# Check if Python 3 is available
if ! command -v python3 &> /dev/null; then
    echo "Error: python3 is not installed"
    exit 1
fi

echo "Python version: $(python3 --version)"

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
        echo "To activate: source ../venv-udales/bin/activate (from u-dales root)"
        exit 0
    fi
fi

echo "Creating virtual environment at: $VENV_DIR"
python3 -m venv "$VENV_DIR"

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

echo ""
echo "=========================================="
echo "Setup complete!"
echo "=========================================="
echo ""
echo "To use the virtual environment:"
echo "  1. Activate:   source ../venv-udales/bin/activate (from u-dales root)"
echo "  2. Run script: python tools/write_inputs.py [config_dir]"
echo "  3. Deactivate: deactivate"
echo ""
echo "Example workflow:"
echo "  cd $UDALES_ROOT"
echo "  source ../venv-udales/bin/activate"
echo "  python tools/write_inputs.py"
echo "  deactivate"
echo ""
