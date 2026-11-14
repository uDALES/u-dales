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
  the RF01 profiles. The record's `ql` is also recomputed from the same record's
  `thl`, `qt` and `exnf` with the uDALES adjustment (checks `rlv`, the
  thermostat timing and that `ql0(k)` belongs to level k, #403).
- `TestRadiation`: `thlprad` for liquid and smoke (`rad_absorber = 'sv'`, also bit-identical to the
  pre-Task-16 `rad_smoke` reference `smoke_pre_task16.txt`, see below), `zi_sv`, and the DYCOMS free-troposphere
  term of the `dycoms_rf01` case hook (flux form vs pointwise form, column integral, subsidence balance),
  against the closed-form longwave expressions
- `TestCaseHooks`: the `&NAMCASES` registry (unknown names abort, `case_hooks = ''` is a no-op), the
  aborts for required parameters that are not set (including `tthermostat`), a namelist read error aborting
  once on all ranks, and the warning for a case group whose hook is not listed
- `TestDecomposition`: radiation, MO surface and moist SGS in a randomised
  short run (`runtime` 2.5 s, `dtmax` 1 s, about 3 steps) on 1, 2 (x) and 4 (2 x 2) ranks. Reference: the serial run.
- `TestSGSOneEqn`: the one-equation SGS model (Task 10) through the in-solver
  `runmode = 1008` (`tests.f90:tests_sgs_oneeqn`), serial, on a synthetic
  stratified state. Kh/Km = 1/Prandtl = 3 where `dthvdz <= 0` and
  `ch1 + ch2*zlt/delta` where `dthvdz > 0` (zlt with `grav/thvf`), interior
  dissipation without the old factor 2, no `e12p(kb)` source without the MO
  surface (work arrays pre-filled with garbage), and the DALES k = 1 sources with
  it, and a synthetic fully solid level (`thvf = -999`), which must get the
  neutral length scale and no buoyancy source. Reference: the closed-form
  expressions, constants recomputed independently.
- `TestMOKernels`: in-solver `runmode = 1012` (`tests.f90:tests_mo_kernels`): `psim`/`psih` against Businger-Dyer
  and Beljaars-Holtslag reference values, `phim`/`phih`, the `dudz_mo`/`dvdz_mo` formula
  `u* phim(z1/L) / (kappa z1) * u/|U|`, and the `calc_obl_iter` round trip from the bulk Richardson number of a known `L`.
- `TestThermoLevels`: in-solver `runmode = 1010` (`tests.f90:tests_thermo_levels`):
  after `thermodynamics`, `ql0(k)` and `ql0h(k)` are the saturation adjustment of
  `thl0`/`qt0` and `thl0h`/`qt0h` at the same level, cell by cell (#403).
- `TestLstend`: the large-scale subsidence through the in-solver `runmode = 1009`
  (`tests.f90:tests_lstend`), serial. `thl`, `qt`, `sv`, `u` and `v` vary in x,
  `whls` changes sign with height and `lmomsubs` is on. Every cell is compared with
  the DALES dev 075ac4e `lstend` upwind formula, including both kb rules.

`smoke_pre_task16.txt` pins the smoke longwave path exactly (`assert_array_equal`) to
the output of the removed `rad_smoke = .true.` code (commit 5c58b2e3, Debug build). It
proves that the Task 16 refactor (`rad_absorber = 'sv'`) changed nothing. Any later,
intended change to the smoke path, the laminar dry column or the thermostat output
breaks it, and a Release `UDALES_BUILD` may differ in the last bit: then retire the test
(or regenerate the fixture from the new code and say so in the commit). Retire it at the
latest when the `thermo` branch is merged.
