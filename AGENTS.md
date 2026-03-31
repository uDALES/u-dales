# uDALES Agent Notes

## What This Repo Is

uDALES is an urban large-eddy simulation codebase centered on a Fortran MPI
solver, with Python and MATLAB tooling for preprocessing, geometry generation,
 radiation, vegetation, and postprocessing.

## Fast Mental Map

- `src/`: core Fortran solver
- `tools/python/`: current Python tooling and tests
- `tools/matlab/`: legacy MATLAB postprocessing and geometry helpers
- `tools/IBM/`: immersed-boundary preprocessing utilities
- `tools/SEB/`: surface energy balance utilities
- `tools/syntheticInflow/`: inflow-generation tooling
- `tools/View3D/`: external view-factor calculator used by preprocessing
- `examples/`: runnable example cases and reference inputs
- `tests/`: branch-specific test assets organized by scope
- `docs/`: user guides, tutorials, and reference docs
- `2decomp-fft/`: parallel decomposition library submodule/dependency

## Solver Layout

Key Fortran files under `src/`:

- `program.f90`: main entry point and runmode dispatch
- `modglobal.f90`: global parameters/state
- `modstartup.f90`, `readinput.f90`: startup and input parsing
- `modfields.f90`, `modsave.f90`, `modfielddump.f90`, `modstatsdump.f90`, `modstat_nc.f90`: field/state/output handling
- `modboundary.f90`, `modinlet*.f90`, `moddriver.f90`: boundaries, inflow, driver runs
- `modforces.f90`, `modthermodynamics.f90`, `modsubgrid*.f90`, `modpois.f90`: core physics/numerics
- `modibm*.f90`, `initfac.f90`: immersed boundary method
- `vegetation.f90`, `modtrees.f90`, `modEB.f90`, `heatpump.f90`, `modchem.f90`, `modpurifiers.f90`: branch/domain-specific physics
- `tests.f90`: in-solver test routines triggered by runmode

## Tooling

### Python

`tools/python/` is the modern tooling stack.

- `udbase.py`: main analysis/postprocessing interface
- `udprep/`: preprocessing modules for grid, radiation, vegetation, scalars, IBM, SEB, and initialization
- `udgeom/`: geometry generation, manipulation, and viewing
- `tests/`: Python unit tests
- `examples/`: notebooks and tutorial scripts
- `namelists.json`: tooling-facing namelist metadata

### MATLAB

`tools/matlab/` still matters for some workflows.

- `udbase.m`: legacy postprocessing entry point
- `merge_stat.m`, `time_average.m`, `coarsegrain_field.m`: common postprocessing routines
- `+udgeom/`: MATLAB geometry generation/manipulation helpers
- `tools/SEB/`: older MATLAB/Fortran shortwave and SEB preprocessing code, including the original `directShortwave` implementation

### Other Tools

- `tools/IBM/`: preprocessing and geometry tools for IBM workflows
- `tools/examples/`: assorted helper examples
- `tools/singularity/`: container-related utilities

## Examples

`examples/` contains case directories like `001`, `002`, `024`, `101`, `102`,
`201`, `949`, `950`, and `999`.

Typical case contents:

- `namoptions.*`: main simulation configuration
- `lscale.inp.*`, `prof.inp.*`: forcing/profile inputs
- `facets.*`, `factypes.*`, `facet_sections_*`: geometry/facet inputs
- `solid_*`, `fluid_boundary_*`: IBM sparse boundary data
- optional STL, vegetation, radiation, or scalar inputs depending on the case

## Tests

There are two automated test homes in this repo.

- `tests/`: repo-level tests and shared case fixtures, mainly for compiled
  uDALES behaviour, MPI runs, integration cases, and branch regression
- `tools/python/tests/`: Python package tests for `udbase`, `udprep`,
  `udgeom`, solar utilities, and related tooling

Do not assume all tests belong under `tests/`. Put a new test where the thing
being validated actually lives.

`tests/` itself is organized by scope.

- `tests/regression/`: branch-to-branch comparison harness
- `tests/integration/sparse_ijk/`: MPI sparse-input integration test
- `tests/integration/tree_input/`: vegetation/tree integration assets used by solver and Python tests
- `tests/system/`: placeholder for heavier whole-code validation
- `tests/unit/`: placeholder for future isolated non-Python tests
- `tools/python/tests/`: Python unit tests

Primary test entry points:

- `tests/regression/run_tests.py`
- `tests/integration/sparse_ijk/run_test.sh`
- `python -m unittest discover -s tools/python/tests`

Placement guidance:

- Put new tests in `tests/` when they validate executable behaviour, MPI runs,
  shared case assets, or cross-branch/regression outputs.
- Put new tests in `tools/python/tests/` when they validate Python APIs or
  internal Python implementations.
- Put exploratory, visual, or solver-development scripts in
  `tools/python/examples/` or a dedicated dev area, not in the automated test suites.

## Build And Run

Typical debug build:

```bash
cmake -S . -B build/debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build/debug -j 4
```

Typical release build:

```bash
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release
cmake --build build/release -j 4
```

On the cluster, MPI-enabled builds and tests may require loading the project
module stack first.

## Python Environment

- Use the repo-local virtual environment at `.venv/` for Python tooling and tests.
- `tools/python/setup_venv.sh` expects a Python interpreter with development
  headers because the `directshortwave_f2py` extension is mandatory.
- On this cluster, `/opt/pbs/python/bin/python3` is a known-good interpreter
  for creating a build-capable venv when `/usr/bin/python3` lacks `Python.h`.
- View3D and `directshortwave_f2py` are required parts of the Python radiation
  workflow, not optional extras.
- The current `.venv` is intended to be the single canonical environment; avoid
  maintaining parallel repo-local virtual environments.

## Branch-Specific Gotchas

- `tools/View3D` is a submodule. After switching branches, run:
  `git submodule update --init --recursive tools/View3D`
- Build View3D after switching or cloning if the executable is missing:
  `./tools/build_preprocessing.sh common`
- Untracked test artifacts can survive branch switches; Git will not remove them automatically.
- This branch and others may share the same high-level test layout while keeping different test contents.
- CI and scripts should use the refactored regression path:
  `tests/regression/run_tests.py`
