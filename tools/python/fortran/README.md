# Fortran f2py wrappers

This folder contains f2py wrapper sources for Fortran routines used by the
Python tooling. The authoritative Fortran implementations live under
`tools/SEB`.

## Build (PowerShell)

From the `u-dales` repo root:

```powershell
python -m numpy.f2py -c -m directshortwave_f2py tools/python/fortran/directShortwave_f2py.f90
```

Optional optimization:

```powershell
python -m numpy.f2py -c -m directshortwave_f2py tools/python/fortran/directShortwave_f2py.f90 --opt="-O3"
```
