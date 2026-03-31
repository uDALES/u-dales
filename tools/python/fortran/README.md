# Fortran f2py wrappers

This folder contains f2py wrapper sources for Fortran routines used by the
Python tooling. The authoritative Fortran implementations live under
`tools/SEB`.

## Linux / macOS

From the `u-dales` repo root:

```bash
source .venv/bin/activate
./tools/python/fortran/build_f2py.sh
```

## Windows (PowerShell)

From the `u-dales` repo root:

```powershell
.\tools\python\fortran\build_f2py.ps1
```

## Requirements

- Python development headers for the interpreter used to build the extension
- `numpy`
- Python-side build helpers from `tools/python/requirements-build.txt`
- a Fortran compiler such as `gfortran`

When building with `gfortran`, the wrapper source may require the
`-ffree-line-length-none` flag. The provided shell and PowerShell build scripts
set the needed flags.
