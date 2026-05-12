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
- Python packages from `tools/python/requirements.txt` (including `numpy`, `scipy`, `xarray`, `netCDF4`, `matplotlib`, `f90nml`, `numba`, `pvlib`, `plotly`, `ipykernel`)
- Python build packages from `tools/python/requirements-build.txt` (`meson`, `ninja`)

You can inspect your local tool versions with:

    python --version
    cmake --version
    gcc --version
    gfortran --version
    ninja --version

## 1. Create or activate the shared virtual environment

From the project top-level directory (the parent of `u-dales`):

    cd "C:\Users\<user>\...\uDALES"
    python -m venv .venv
    .\.venv\Scripts\Activate.ps1

If the venv already exists, only activation is needed:

    .\.venv\Scripts\Activate.ps1

## 2. Install Python dependencies

Move to the repository root and install runtime + build packages:

    cd .\u-dales
    python -m pip install --upgrade pip
    python -m pip install -r tools\python\requirements.txt
    python -m pip install -r tools\python\requirements-build.txt

## 3. Important Windows path note (spaces)

If your repository path contains spaces (for example OneDrive paths), `numpy.f2py` can fail during module builds.

Recommended workaround: map the parent folder to a drive letter without spaces.

From the project top-level directory:

    subst U: "C:\Users\<user>\...\uDALES"
    cd /d U:\u-dales

Use `U:\.venv\Scripts\python.exe` as the Python executable for preprocessing builds.

## 4. Build preprocessing targets (Windows)

### Option A: Build Python/Fortran preprocessing libraries (no View3D)

Use this when you do not yet have a Windows C++ compiler installed.

    cmake -G Ninja -S tools/preprocessing -B tools/preprocessing/build ^
      -DCMAKE_POLICY_VERSION_MINIMUM=3.5 ^
      -DPREPROCESSING_PYTHON_EXECUTABLE="U:\.venv\Scripts\python.exe" ^
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
      -DPREPROCESSING_PYTHON_EXECUTABLE="U:\.venv\Scripts\python.exe" ^
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

    .\.venv\Scripts\Activate.ps1
    subst U: "C:\Users\<user>\...\uDALES"
    cd /d U:\u-dales

Then run preprocessing tools/tests as needed.
