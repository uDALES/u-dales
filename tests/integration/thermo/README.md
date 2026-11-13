# Thermo integration tests

Solver tests for the moist thermodynamics, Monin-Obukhov (MO) surface and
parametrised radiation. Each test writes a small flat case (16 x 16 x 64,
`thermo_case.py`), runs `u-dales` and compares the netCDF output with an
independent Python reference. Output is single precision, so tolerances allow
float32 rounding.

Run: `python tests/integration/thermo/test_thermo.py -v [TestClass]`.

## Requirements

- `UDALES_BUILD`: path to the `u-dales` executable (default `build/debug/u-dales`);
  tests skip if it is missing.
- `UDALES_RUNTIME_MODULES`: module list of the stack the executable was built
  with (netCDF, netCDF-Fortran, FFTW, compiler runtime).
- Python with `numpy` and `netCDF4`.
- Up to 4 MPI ranks (`TestDecomposition`). uDALES never oversubscribes; with
  fewer CPUs the run is refused by `tests/integration/mpi_env.py`.

## Test classes

- `TestInitChecks`: an unset `thls` sentinel aborts at start-up; a dry case
  without `qts` runs. Reference: the documented start-up checks.
- `TestSurfaceSign`: a positive `wtsurf`/`wqsurf` is a downward flux (decision D1).
  Reference: sign of the resulting surface-layer tendency.
- `TestMO`: `ustar` and the Obukhov length from the `BCbotm = 4` surface, and a
  prescribed `ustar`. Reference: a Python re-implementation of the MO iteration
  with the Businger-Dyer functions.
- `TestCalthv`, `TestCalthvGuard`: the Deardorff surface-level `thv` from MO
  gradients; unsaturated and no-MO cases; a uniform saturated column must not
  divide by zero. Reference: analytic `thv` from `thl`, `qt`, `ql`.
- `TestThermostat`: cloudy-column diagnostics (`ql`, `lwp_bar`) and a dry
  column without moist diagnostics. Reference: Python saturation adjustment of
  the RF01 profiles.
- `TestRadiation`: `thlprad` for liquid, smoke and the DYCOMS large-scale term
  against the closed-form longwave expressions, and the balance of the
  free-troposphere radiative cooling against subsidence warming.
- `TestDecomposition`: radiation, MO surface and moist SGS in a randomised
  10-step run on 1, 2 (x) and 4 (2 x 2) ranks. Reference: the serial run.
