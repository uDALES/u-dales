# Fortran f2py wrappers

This folder contains f2py wrapper sources for Fortran routines used by the
Python tooling. The authoritative Fortran implementations live under
`tools/SEB`.

## Build (PowerShell)

From the `u-dales` repo root:

```powershell
./tools/python/fortran/build_f2py.f90
```

Install build dependencies (and Fortran compiler). From within the virtual environment run:
```powershell
conda install -c conda-forge gfortran meson ninja --solver=classic
```

Build dependencies are also listed in `tools/python/requirements-build.txt`.
