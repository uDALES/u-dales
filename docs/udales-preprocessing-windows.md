# Windows setup for preprocessing libraries

If you use the repo skills (`udales-detect` and `udales-exec`), a coding agent can do most of this setup and build workflow for you automatically.

The Python preprocessing path on Windows is still experimental.

This page explains how to set up and build the uDALES preprocessing libraries on native Windows using PowerShell.

It covers:

- Python environment setup
- build dependencies for f2py modules
- building preprocessing targets with CMake
- common Windows pitfalls (paths with spaces, missing C++ compiler for View3D)

## Scope

The preprocessing stack contains:

- Python dependencies from `tools/python/requirements.txt`
- Python build dependencies from `tools/python/requirements-build.txt`
- Fortran/C preprocessing targets under `tools/preprocessing`
- f2py extension modules:
  - `tools/python/udprep/directshortwave_f2py*.pyd`
  - `tools/python/udprep/ibm_preproc_f2py*.pyd`

## Requirements (read first)

Install and make available on PATH:

- Python 3.9+
- CMake
- Ninja
- GCC/GFortran (for example from conda-forge)

Optional but required for full preprocessing build (View3D):

- A C++ compiler on Windows (for example MSVC C++ Build Tools or MinGW g++)

### Reference installation on this computer

The commands on this page were validated with the following local setup:

- Python 3.12.7
- CMake 4.1.2
- GCC 15.2.0 (`gcc` from conda-forge)
- GFortran 15.2.0 (`gfortran` from conda-forge)
- Ninja 1.13.0
- Python packages from `tools/python/requirements.txt` (including `numpy`, `scipy`, `pandas`, `xarray`, `netCDF4`, `h5netcdf`, `trimesh`, `shapely`, `matplotlib`, `f90nml`, `numba`, `pvlib`, `pyvista`, `ipykernel`). `plotly` is an optional visualization backend in `requirements-plotly.txt`.
- Python build packages from `tools/python/requirements-build.txt` (`meson`, `ninja`)

You can inspect your local tool versions with:

    python --version
    cmake --version
    gcc --version
    gfortran --version
    ninja --version

## 1. Create or update the preprocessing environment

Run the PowerShell setup script from the repository root:

    .\tools\python\setup_venv.ps1

This creates or updates `tools\python\.venv`, installs the Python dependencies,
and builds the preprocessing artifacts. To build only View3D:

    .\tools\python\setup_venv.ps1 -BuildTarget view3d

If the venv already exists, the setup script asks whether to recreate it.

## 2. Activate the environment

From the repository root:

    tools\python\.venv\Scripts\Activate.ps1

Then run preprocessing tools/tests as needed, for example:

    $env:PREPROC_NCPU = "8"
    python tools\write_inputs.py examples\999

## 3. Important Windows path note (spaces)

If your repository path contains spaces (for example OneDrive paths), `numpy.f2py` can fail during module builds.

Recommended workaround: map the parent folder to a drive letter without spaces.

From the project top-level directory:

    subst U: "C:\Users\<user>\...\uDALES"
    cd /d U:\u-dales

Run the setup script from the mapped path:

    .\tools\python\setup_venv.ps1

If you need an explicit interpreter, set `PYTHON_BIN` before running the setup
script.

## 4. Manual preprocessing rebuilds (Windows)

The setup script normally handles this. Use the manual CMake commands below only
when debugging the preprocessing build directly. Activate the venv first:

    tools\python\.venv\Scripts\Activate.ps1
    $py = (Get-Command python).Source

### Option A: Build Python/Fortran preprocessing libraries (no View3D)

Use this when you do not yet have a Windows C++ compiler installed.

    cmake -G Ninja -S tools/preprocessing -B tools/preprocessing/build ^
      -DCMAKE_POLICY_VERSION_MINIMUM=3.5 ^
      -DPREPROCESSING_PYTHON_EXECUTABLE="$py" ^
      -DCMAKE_C_COMPILER=gcc ^
      -DCMAKE_Fortran_COMPILER=gfortran ^
      -DBUILD_PREPROCESSING_VIEW3D=OFF

    cmake --build tools/preprocessing/build --target preprocessing_tools

This builds:

- `tools/preprocessing/build/bin/IBM_preproc.exe`
- `tools/python/udprep/directshortwave_f2py*.pyd`
- `tools/python/udprep/ibm_preproc_f2py*.pyd`

### Option B: Build full preprocessing stack (including View3D)

Requires a Windows C++ compiler in PATH.

    cmake -G Ninja -S tools/preprocessing -B tools/preprocessing/build ^
      -DCMAKE_POLICY_VERSION_MINIMUM=3.5 ^
      -DPREPROCESSING_PYTHON_EXECUTABLE="$py" ^
      -DCMAKE_C_COMPILER=gcc ^
      -DCMAKE_Fortran_COMPILER=gfortran

    cmake --build tools/preprocessing/build --target preprocessing_tools

## 5. Verify artifacts

Run these checks from `u-dales` root:

    Test-Path tools\preprocessing\build\bin\IBM_preproc.exe
    Get-ChildItem tools\python\udprep -Filter "directshortwave_f2py*.pyd"
    Get-ChildItem tools\python\udprep -Filter "ibm_preproc_f2py*.pyd"

For full builds with View3D enabled:

    Test-Path tools\preprocessing\build\bin\view3d.exe

## Troubleshooting

### CMake selects NMake and fails

Symptom:

- CMake reports NMake generator errors or missing `nmake`

Fix:

- Use `-G Ninja` explicitly.

### f2py fails with "Permission denied" near a path fragment

Symptom:

- Errors referencing a truncated path segment such as `C:/Users/.../OneDrive`

Cause:

- Path with spaces during f2py invocation.

Fix:

- Build from a no-space alias path using `subst` (for example `U:`).

### View3D fails because no C++ compiler is found

Symptom:

- CMake error: no CMAKE_CXX_COMPILER could be found.

Fix:

- Install a Windows C++ compiler and ensure it is in PATH.
- Reconfigure and rebuild without `-DBUILD_PREPROCESSING_VIEW3D=OFF`.

## Daily usage

From project top-level:

    cd .\u-dales
    tools\python\.venv\Scripts\Activate.ps1
    $env:PREPROC_NCPU = "8"
    python tools\write_inputs.py examples\999
