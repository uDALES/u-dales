# Direct Shortwave Integration Tests

This directory contains repo-level direct shortwave integration tests.

Contents:

- `test_directshortwave.py`: direct shortwave integration tests on shared
  committed cases under `tests/cases/`

Covered cases:

- `tests/cases/100/`
  - legacy Fortran `tools/SEB/directShortwave.f90`
  - Python `scanline` via the mandatory `directshortwave_f2py` wrapper
  - Python `facsec`
  - Python `moller`
- `tests/cases/525/`
  - Python `facsec` with vegetation
  - Python `moller` with vegetation

Why it lives here:

- it is an integration/reference test, not a pure Python unit test
- it relies on shared case data under `tests/cases/`
- it checks agreement across multiple implementations rather than one Python API

Current status:

- this is a preprocessing/tooling integration suite
- it is useful development coverage, but it is not intended to define the
  primary supported GitHub Actions merge gate yet
- when prerequisites such as `gfortran` or `directshortwave_f2py` are absent,
  the test should skip cleanly rather than fail opaquely

The committed fixtures include geometry, facet sections, fluid boundary files,
solid masks, and vegetation inputs where needed because the ray-based methods
depend on them.

To run it:

```bash
source ../.venv/bin/activate
python tests/integration/directshortwave/test_directshortwave.py
```

Current practical test settings:

- `ray_density = 12` on case `100`
- `ray_density = 12` on case `525`

For case `100`, this value was chosen from a local convergence study as a good
compromise for `facsec`: close to converged on total energy and facet-level
agreement without the cost of very high ray densities.
