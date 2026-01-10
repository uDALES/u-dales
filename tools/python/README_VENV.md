# Python Virtual Environment Setup

This directory contains scripts to set up a Python virtual environment for running uDALES preprocessing tools.

**Note:** The virtual environment is created **outside** the u-dales repository (in the parent directory) as `venv-udales/`. This keeps the repository clean and prevents accidentally committing large dependency files.

## One-Time Setup

### Windows (PowerShell)

```powershell
# Navigate to the u-dales root directory
cd D:\Postdoc1\simulation\ud2dev\u-dales

# Run the setup script
.\tools\python\setup_venv.ps1
```

### Linux/WSL (Bash)

```bash
# Navigate to the u-dales root directory
cd /mnt/d/Postdoc1/simulation/ud2dev/u-dales

# Run the setup script
./tools/python/setup_venv.sh
```

## Daily Usage

### Activate Virtual Environment

**Windows:**
```powershell
# From u-dales root directory
..\venv-udales\Scripts\Activate.ps1
```

**Linux/WSL:**
```bash
# From u-dales root directory
source ../venv-udales/bin/activate
```

### Run Your Scripts

Once activated, you'll see `(venv-udales)` in your terminal prompt.

```bash
# Example 1: Run with default hardcoded paths
python tools/write_inputs.py

# Example 2: Run with config directory
python tools/write_inputs.py /path/to/config/directory

# Example 3: Run test script
python tools/python/tests/init_udprep.py
```

#### Config Directory Requirements

When running with a config directory, it must contain:

1. **config.sh** - Shell script with required variables:
   - `DA_EXPDIR` - Path to experiment directory (or parent experiments directory)
   - `DA_TOOLSDIR` - Path to tools directory

2. **namoptions.XXX** - Fortran namelist file where:
   - XXX is a 3-digit experiment number (e.g., 001, 105, 999)
   - Contains `iexpnr` variable matching the XXX in the filename
   - Directory name should also match the experiment number

#### Supported config.sh Formats

All of the following formats are supported:

**Format 1: Command Substitutions (Dynamic Paths)**
```bash
# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Set paths relative to the script location
export DA_EXPDIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
export DA_TOOLSDIR="$(cd "${SCRIPT_DIR}/../../u-dales/tools" && pwd)"
export DA_BUILD="$(cd "${SCRIPT_DIR}/../../u-dales/build/release" && pwd)/u-dales"
export DA_WORKDIR="$(cd "${SCRIPT_DIR}/../../outputs" && pwd)"
export NCPU=8
```

**Format 2: Simple Paths with Tilde Expansion**
```bash
export DA_EXPDIR=~/simulation/testu2/experiments
export DA_TOOLSDIR=~/simulation/testu2/u-dales/tools
export DA_BUILD=~/simulation/testu2/u-dales/build/release/u-dales
export DA_WORKDIR=~/simulation/testu2/outputs
export NCPU=16
```

**Format 3: Mixed with Environment Variables**
```bash
# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Set paths relative to the script location
export DA_EXPDIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
export DA_TOOLSDIR="$(cd "${SCRIPT_DIR}/../../u-dales/tools" && pwd)"
export DA_BUILD="$(cd "${SCRIPT_DIR}/../../u-dales/build/release" && pwd)/u-dales"
export DA_WORKDIR=$EPHEMERAL/antwerp
export NCPU=128
export NNODE=1
export WALLTIME="72:00:00"
export MEM="128gb"
```

**Note:** Additional variables (DA_BUILD, DA_WORKDIR, NCPU, NNODE, WALLTIME, MEM, etc.) in config.sh are allowed and won't cause any issues. Only `DA_EXPDIR` and `DA_TOOLSDIR` are required for the preprocessing scripts.

### Deactivate Virtual Environment

When you're done:
```bash
deactivate
```

## What Gets Installed

The setup script installs all dependencies from `requirements.txt`:

- **Core dependencies:**
  - numpy >= 1.20.0
  - scipy >= 1.7.0
  - xarray >= 0.19.0
  - netCDF4 >= 1.5.0
  - h5netcdf >= 1.0.0
  - trimesh >= 3.9.0
  - matplotlib >= 3.3.0
  - f90nml >= 1.4.4
  - plotly >= 5.0.0

- **Optional (for notebooks):**
  - ipykernel >= 6.0.0
  - nbformat >= 4.2.0

## Troubleshooting

### Windows: Script Execution Policy

If you get an error about script execution being disabled:

```powershell
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
```

### Missing Python

Make sure Python 3 is installed:
- Windows: Download from [python.org](https://www.python.org/downloads/)
- Linux/WSL: `sudo apt-get update && sudo apt-get install python3 python3-venv python3-pip`

### Recreating the Environment

To completely recreate the virtual environment:

**Windows:**
```powershell
# From parent directory of u-dales
Remove-Item -Recurse -Force venv-udales
# From u-dales root
.\tools\python\setup_venv.ps1
```

**Linux/WSL:**
```bash
# From parent directory of u-dales
rm -rf venv-udales
# From u-dales root
./tools/python/setup_venv.sh
```

## Cross-Platform Notes

The virtual environment is **platform-specific** and is created **outside the repository** in the parent directory. If you work on both Windows and WSL:

- The same `venv-udales/` directory works for both if they share the same filesystem (e.g., Windows accessing WSL files)
- For separate filesystems, create the venv separately in each environment
- The venv location (`../venv-udales/`) keeps it out of version control automatically

The Python scripts themselves (`write_inputs.py`, `udprep/*.py`) are fully cross-platform and work identically on both Windows and WSL.
