# Python Environment Setup

This document describes how to create the Python virtual environment required
for uDALES preprocessing, post-processing, and tests.

Python 3.9 or newer is required. Older interpreters are not supported.

---

## Quick Start

### Linux / WSL

Run the setup script from the repository root. It creates the virtual
environment, installs all dependencies, and builds the preprocessing
tools (View3D and f2py extension modules):

```bash
# For a local machine
bash tools/python/setup_venv.sh common

# For the Imperial HPC machine
bash tools/python/setup_venv.sh icl
```

Then activate and use:

```bash
source tools/python/.venv/bin/activate
python tools/write_inputs.py
deactivate
```

### Windows (PowerShell)

Requires Python 3.9+, CMake, Ninja, GCC, and GFortran on `PATH`
(e.g. from conda-forge). Run from the repository root:

```powershell
.\tools\python\setup_venv.ps1
```

Then activate and use:

```powershell
tools\python\.venv\Scripts\Activate.ps1
python tools\write_inputs.py
deactivate
```

> **Path with spaces:** If the repository path contains spaces (e.g. OneDrive),
> map the parent directory to a drive letter first to avoid f2py failures:
> ```powershell
> subst U: "C:\Users\<user>\...\uDALES"
> cd /d U:\u-dales
> $env:PYTHON_BIN = "U:\.venv\Scripts\python.exe"
> .\tools\python\setup_venv.ps1
> ```

---

## Setup Script Reference

Both scripts perform the same logical setup steps. The Linux script requires the build system as its first positional argument; the PowerShell script uses named parameters with Windows/local defaults. On Windows, `tools/build_preprocessing.sh` is bash-only, so the PowerShell script drives CMake directly instead.

### Linux / WSL

```bash
bash tools/python/setup_venv.sh <build_system> [build_target]
```

### Windows (PowerShell)

```
.\tools\python\setup_venv.ps1 [-BuildSystem <system>] [-BuildTarget <target>]
```

### Arguments

| Argument | Allowed values | Default | Description |
|---|---|---|---|
| `build_system` | `common`, `icl` (Linux only) | Linux: required; Windows: `common` | `common` — local Linux/WSL/Windows; `icl` — Imperial College London HPC cluster (Linux only, loads environment modules) |
| `build_target` | `view3d`, `preprocessing_tools` | `preprocessing_tools` | `view3d` — View3D executable only; `preprocessing_tools` — View3D + f2py extension modules |

### Environment variable overrides

| Variable | Default (Linux) | Default (Windows) | Description |
|---|---|---|---|
| `PYTHON_BIN` | `python3` | `python` | Python interpreter for creating the venv. Override when the default lacks development headers. |
| `VENV_DIR` | `tools/python/.venv` | `tools\python\.venv` | Explicit venv path. Falls back to `.venv` at the repo root if that legacy directory exists. |

### Examples (Linux / WSL)

```bash
# Local setup
bash tools/python/setup_venv.sh common

# View3D only
bash tools/python/setup_venv.sh common view3d

# ICL cluster, full build
bash tools/python/setup_venv.sh icl preprocessing_tools

# Use a specific Python interpreter
PYTHON_BIN=/opt/pbs/python/bin/python3 bash tools/python/setup_venv.sh common
```

### Examples (Windows)

```powershell
# Default local setup
.\tools\python\setup_venv.ps1

# View3D only
.\tools\python\setup_venv.ps1 -BuildTarget view3d

# Use a specific Python interpreter
$env:PYTHON_BIN = "C:\Python312\python.exe"
.\tools\python\setup_venv.ps1
```

### Virtual environment behaviour

- If the venv directory does not exist it is created from scratch.
- If it already exists the script prompts whether to recreate it (default: **N**).
  Answering N skips package installation/rebuild and runs validation checks.

---

## What the Setup Scripts Do

Both scripts perform the same logical steps:

1. Validate `build_system` and `build_target` inputs.
2. Resolve the virtual environment directory using the priority chain:
   explicit override → existing default venv → existing legacy venv → new default venv.
3. Check that the chosen Python interpreter is available.
4. Create the virtual environment with `python -m venv`.
5. Install runtime dependencies from `tools/python/requirements.txt`.
6. Install build-time dependencies from `tools/python/requirements-build.txt`
   (`meson` and `ninja`, used to compile the Fortran extensions).
7. Build the preprocessing tools:

   | Artifact | Linux path | Windows path |
   |---|---|---|
   | View3D executable | `tools/preprocessing/build/bin/view3d` | `tools\preprocessing\build\bin\view3d.exe` |
   | MATLAB compat symlink/copy | `tools/View3D/build/src/view3d` | `tools\View3D\src\View3D.exe` |
   | directshortwave f2py | `tools/python/udprep/directshortwave_f2py*.so` | `tools\python\udprep\directshortwave_f2py*.pyd` |
   | IBM preproc f2py | `tools/python/udprep/ibm_preproc_f2py*.so` | `tools\python\udprep\ibm_preproc_f2py*.pyd` |

   On Linux/WSL this is done via `tools/build_preprocessing.sh`.
   On Windows this is done by calling CMake/Ninja directly (with `-G Ninja`,
   `-DCMAKE_C_COMPILER=gcc`, `-DCMAKE_Fortran_COMPILER=gfortran`).

---

## Manual Rebuilds

### Linux / WSL

```bash
source tools/python/.venv/bin/activate
tools/build_preprocessing.sh common preprocessing_tools   # full rebuild
tools/build_preprocessing.sh common view3d                # View3D only
```

### Windows

```powershell
tools\python\.venv\Scripts\Activate.ps1
$py = (Get-Command python).Source
cmake -G Ninja -S tools\preprocessing -B tools\preprocessing\build `
    -DCMAKE_POLICY_VERSION_MINIMUM=3.5 `
    -DPREPROCESSING_PYTHON_EXECUTABLE="$py" `
    -DCMAKE_C_COMPILER=gcc -DCMAKE_Fortran_COMPILER=gfortran `
    -DBUILD_PREPROCESSING_DIRECTSHORTWAVE_F2PY=ON `
    -DBUILD_PREPROCESSING_IBM_F2PY=ON
cmake --build tools\preprocessing\build --target preprocessing_tools
```

---

## Jupyter Notebooks

Register the venv as a Jupyter kernel to use the tutorial notebooks.

**Linux / WSL:**
```bash
source tools/python/.venv/bin/activate
python -m ipykernel install --user --name udales-py --display-name "uDALES Python"
```

**Windows:**
```powershell
tools\python\.venv\Scripts\Activate.ps1
python -m ipykernel install --user --name udales-py --display-name "uDALES Python"
```

---

## Notes

- `directshortwave_f2py`, `ibm_preproc_f2py`, and View3D are required parts of
  the Python preprocessing workflow and are not optional.
- `requirements-build.txt` covers Python-side build tools only. A system
  Fortran compiler (`gfortran`) and C compiler (`gcc`) must also be on `PATH`.
- On Windows, View3D additionally requires a C++ compiler (e.g. MinGW g++ or
  MSVC Build Tools). If not available, the f2py-only build still works.
- The venv directory (`tools/python/.venv`) is git-ignored. Do not commit it.
- After pulling changes, rerun the setup script (or at minimum
  `pip install -r tools/python/requirements.txt`) to pick up new dependencies.
