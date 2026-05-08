#!/bin/bash
# Setup script to create NetCDF virtual environment for u-dales testing
# Run this once before using ud_test_sim.sh

set -e

# Resolve the directory where this script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="$SCRIPT_DIR/.venv_netcdf"

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

echo "========================================"
echo "u-dales NetCDF Environment Setup"
echo "========================================"
echo ""

# ud_compare_outputs.py requires Python >= 3.9 (f-strings + netCDF4/xarray dependency floors)
if command -v python3 &>/dev/null; then
    PYTHON=python3
elif command -v python &>/dev/null && python -c "import sys; sys.exit(0 if sys.version_info.major == 3 else 1)" &>/dev/null; then
    PYTHON=python
else
    print_red "Error: neither 'python3' nor 'python' (Python 3) found. Please install Python 3.9 or newer and try again."
    exit 1
fi

# Verify minimum version
if ! $PYTHON -c "import sys; sys.exit(0 if sys.version_info >= (3, 9) else 1)" &>/dev/null; then
    print_red "Error: Python 3.9 or newer is required (found $($PYTHON --version)). Please upgrade and try again."
    exit 1
fi
print_green "Using $($PYTHON --version)"

# Check if virtual environment already exists
if [ -d "$VENV_DIR" ]; then
    print_yellow "NetCDF virtual environment already exists."
    read -p "Do you want to recreate it? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        print_yellow "Removing existing virtual environment..."
        rm -rf "$VENV_DIR"
    else
        print_yellow "Using existing virtual environment."
        source "$VENV_DIR/bin/activate"
        SKIP_INSTALL=true
    fi
fi

SKIP_INSTALL=${SKIP_INSTALL:-false}

if [ "$SKIP_INSTALL" = false ]; then
    echo "1. Creating Python virtual environment..."
    $PYTHON -m venv "$VENV_DIR"

    echo "2. Activating virtual environment..."
    source "$VENV_DIR/bin/activate"

    echo "3. Upgrading pip..."
    pip install --upgrade pip --quiet

    echo "4. Installing NumPy..."
    pip install numpy --quiet

    echo "5. Installing NetCDF4..."
    pip install netcdf4 --quiet
fi  # end SKIP_INSTALL

if [ "$SKIP_INSTALL" = true ]; then
    echo "Installation steps 1-5 skipped..."
    echo "Running import checks on existing environment..."
fi

echo "6. Testing NetCDF imports..."
if $PYTHON -c "
import sys
try:
    import numpy as np
    print(f'✓ NumPy {np.__version__} imported successfully')
except ImportError as e:
    print(f'✗ NumPy import failed: {e}')
    sys.exit(1)

try:
    import netCDF4 as nc
    print(f'✓ netCDF4 {nc.__version__} imported successfully')
except ImportError as e:
    print(f'✗ netCDF4 import failed: {e}')
    sys.exit(1)

"; then
    print_green "✓ All NetCDF libraries installed and tested successfully!"
    echo "7. Deactivating virtual environment..."
    deactivate
else
    print_red "✗ Some NetCDF libraries failed to install or import."
    deactivate
    exit 1
fi

echo ""
print_green "========================================"
print_green "NetCDF Environment Setup Complete!"
print_green "========================================"
echo ""
VENV_PATH="$VENV_DIR"
echo "Virtual environment created: $VENV_PATH"
echo ""
echo "Usage:"
echo "  source $VENV_PATH/bin/activate      # Manual activation"
echo "  Run your Python script using $PYTHON"
echo "  deactivate                          # Manual deactivation"
echo ""
echo "To remove the environment later:"
echo "  rm -rf $VENV_PATH"