# Synthetic inflow generator

`modSyntheticInflow.f90` is a standalone program (Xie & Castro, 2008; Kim,
Castro & Xie, 2013 -- see the references at the end of the source file) that
writes the `*driver_*` files a "driven" target simulation (`idriver = 2`,
see `docs/udales-driver-simulations.md`) reads as its inlet condition,
without having to run a full precursor simulation first. It needs four (or,
with temperature/moisture, up to ten) small text input files, produced here
by `tools/python/udprep/synthetic_inflow.py` (the modern, Python replacement
for `write_Reynolds_stress.m`, which is being retired along with the rest of
the MATLAB preprocessing).

## Building

```bash
gfortran -O2 -fopenmp tools/syntheticInflow/modSyntheticInflow.f90 -o synInflow_executable
```

No other flags or libraries are required: the program only reads plain-text
files (`namoptions.<expnr>`, `prof.inp.<expnr>`, and the
`syntheticInflow_inputs/*.txt` files below) and links only `OMP_LIB` --
there is no netCDF dependency (unlike the MATLAB script it replaces, which
reads `tdump.<expnr>.nc` directly with `ncread`). This is also the exact
command used by the existing wrapper script,
`tools/generate_synthetic_inflow.sh <path-to-case>`, which additionally
copies the input files into a scratch work directory, sets
`OMP_NUM_THREADS`/`OMP_PLACES`, runs the executable, and cleans up
afterwards.

On Imperial's CX3 (or any cluster using environment modules), load a
matching compiler module before *and when running* the binary --
`module load tools/prod foss/2023a` was used to build and verify this; the
compiled binary needs the same `libgfortran.so` at runtime, so a bare shell
without that module loaded will fail with e.g.
`libgfortran.so.5: version 'GFORTRAN_10' not found`.

## Running

```
./synInflow_executable <iexpnr>          # e.g. ./synInflow_executable 995
```

`<iexpnr>` (a 3-digit code, e.g. `995`) selects `namoptions.<iexpnr>` in the
current working directory; it is **not** optional in practice -- without it
the program tries to open a file literally named `namoptions.` (blank-padded)
and fails. `namoptions.<iexpnr>`, `prof.inp.<iexpnr>`, and the
`syntheticInflow_inputs/*.txt` input files (copied or symlinked into the
same directory, alongside the executable, as `tools/generate_synthetic_inflow.sh`
does) must all be in the working directory the program is *run* from -- all
its file opens use bare relative filenames, not a configurable input path.

The program reads (via `read_namelist`, `tools/syntheticInflow/modSyntheticInflow.f90:154`
onwards) the `&RUN`, `&DOMAIN`, `&PHYSICS`, `&DYNAMICS` and `&STG` namelist
blocks of `namoptions.<iexpnr>` -- it does **not** read `&DRIVER` (that
namelist only matters to the solver itself, for `idriver`/`tdriverstart`/
`dtdriver`/`driverstore`); only `iexpnr, dtmax, runtime, nprocy` (`&RUN`),
`itot, jtot, ktot, xlen, ylen` (`&DOMAIN`), `ltempeq, lmoist` (`&PHYSICS`)
and, optionally, `lcalc_time_and_length_scale, lcheck_driver_outputs`
(`&STG`) are actually used; any other namelist variables can be omitted
(missing namelist groups/variables are tolerated -- the Fortran only treats
`iostat > 0`, a genuine parse error, as fatal, not `iostat < 0` for a
group/file that simply isn't present).

## Input-file contract

All four (velocity-only) files live in `experiments/<expnr>/syntheticInflow_inputs/`
and are copied into the run directory before invoking the executable. Each
has one header line (its exact text does not matter -- the Fortran skips it
with `READ (unit, '(a80)')`) followed by `ktot + 1` data rows for
`k = 0..nz` (`nz = ktot`), at cell **edges** running from `z = 0` to
`z = zsize`.

`Reynolds_stress_profiles_velocity.txt`
: `z umean R11 R21 R22 R31 R32 R33`, all in SI units (`umean` m/s; the six
  second moments in m^2/s^2). `R11 = u'u'`, `R21 = u'v'`, `R22 = v'v'`,
  `R31 = u'w'`, `R32 = v'w'`, `R33 = w'w'`.

`length_time_scales_u.txt` / `_v.txt` / `_w.txt`
: `z nLy nLz T` for each velocity component. **`nLy` and `nLz` are integer
  grid-point counts, not metres** -- the Fortran declares them `INTEGER` and
  later does `NUY = 2*nluy` etc. to size the digital filter stencil in grid
  points (`calc_filterCoeff_b`). `T` is an integral **time scale in
  seconds** -- it feeds the Xie & Castro (2008) Eq. (14) AR(1) filter
  coefficients `EXP(-(pi/2)*(dtmax/T))`, so `dtmax` (the target
  simulation's timestep) and `T` must be in the same units.

  These three files are only opened at all when `&STG
  lcalc_time_and_length_scale = .FALSE.`. The Fortran default -- set in the
  variable declaration, *before* `read_namelist` even runs, so it applies
  even if `&STG` is omitted entirely -- is `.TRUE.`, in which case
  `calc_time_and_length_scale` (ported from PALM4U's
  `synthetic_turbulence_generator_mod.f90`) computes `nLy`, `nLz` and `T`
  itself from the local grid spacing (`8 * MIN(dx, dy, dz(k))`) and the
  profile's bulk velocity, and never reads these files. Both modes were
  exercised during verification (see below): the tiny smoke-test case relied
  on the self-computed scales.

Both input paths in `udprep/synthetic_inflow.py` (from a `tdump.<expnr>.nc`
precursor statistics file, or from a precursor's own `*driver_*` plane
records) write all profiles on this same `ktot + 1`-point edge grid, using
the same interpolation convention as the MATLAB `obtain_profile` it
replaces: the bottom edge is forced to `0` (a no-slip / zero-variance-at-
the-wall convention, for the mean and stress file only -- not for the
length/time-scale files, where the bottom edge instead holds the same value
as the first interior level, mirroring `calc_time_and_length_scale`'s own
`nluy(0) = nluy(1)` etc.), interior edges are linearly interpolated, and the
top edge is linearly extrapolated from the last two source points. See the
module docstring and each function's docstring in
`tools/python/udprep/synthetic_inflow.py` for the exact formulas and which
line of `write_Reynolds_stress.m` each one replaces.

## Verification performed

- Built with the command above (`foss/2023a`: GCC 12.3.0 / gfortran),
  0 warnings.
- Unit tests: `tools/python/tests/test_synthetic_inflow.py` (20 cases)
  exercise both input paths against synthetic data with known statistics
  (an AR(1) time series of known correlation time, a sinusoidal-in-y
  structure of known wavelength, a known mean profile), the file
  read/write round trip, and the physical-realizability validation
  (`ReynoldsStressProfile.validate` / `LengthTimeScales.validate`).
- End-to-end smoke test: a tiny case (`itot=8, jtot=16, ktot=8, nprocy=2,
  dtmax=1, runtime=3`, velocity-only) with inputs generated by
  `python -m udprep.synthetic_inflow from-tdump` from a synthetic
  `tdump.nc`. The generator produced:
  - `tdriver_000.<expnr>`: 32 bytes = `(runtime/dtmax + 1) * 8` = 4 records
    of one `real(8)` each, with strictly increasing time stamps
    `0, 1, 2, 3`.
  - `{u,v,w}driver_{000,001}.<expnr>`: 3200 bytes = 4 records of
    `(jtot/nprocy + 2*jh) * (ktot + 2*kh) * 8` = `(8+2)*(8+2)*8` = 800 bytes
    each (`jh = kh = 1` for the default `iadv_mom = iadv_cd2` scheme).
  - Reading these back with `udprep.synthetic_inflow.assemble_driver_plane`
    reproduced the expected `(nt, jtot, ktot)` shape and a plausible mean
    profile.

### What was not obvious

- The runtime linker needs the *same* compiler module loaded as the build
  (see "Building" above) -- easy to miss on a cluster where the build and
  run happen in different shells.
- On this login node, running with `OMP_NUM_THREADS=2` reliably hung (0%
  CPU, no crash, no error) partway through the very first `calc_psi`
  parallel region for this tiny case; `OMP_NUM_THREADS=1` ran the same case
  in ~0.01 s. `tools/generate_synthetic_inflow.sh` sets `OMP_NUM_THREADS=8`
  and `OMP_PLACES=cores` for real (much larger) cases -- worth setting
  `OMP_PLACES=cores` even for small/interactive runs, and treating an
  apparently-stuck run at low CPU usage as a threading issue to retry with
  fewer threads, rather than a hang to wait out.
- `tdriver_000.<expnr>` is the only timestamp file ever written, regardless
  of `nprocy` -- both this generator's `write_driver` and the solver's own
  `writedriverfile` (`src/moddriver.f90`) only write it when `driverid == 0`
  (equivalently, only rank/column 0's `cdriverid` is ever `'000'` for this
  purpose); do not expect `tdriver_001.<expnr>` etc. to appear.
- The command-line argument to the executable *only* selects which
  `namoptions.<N>` file to open; the experiment number used for output
  filenames and looked up for `prof.inp.<N>` comes from `iexpnr` inside that
  file, so keep the two consistent.
