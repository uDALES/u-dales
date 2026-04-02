# Shared Test Cases

This directory contains committed case fixtures shared by repo-level tests.

Current cases:

- `100/`: no-tree direct shortwave reference case
- `101/`: IBM sparse-input reference case
- `525/`: flat-terrain tree case for direct shortwave integration checks
- `526/`: reduced tree case for legacy vegetation regression against `v2.2.0`

Current consumers:

- `integration/ibm_sparse_input/run_test.sh`
  - `101/` as the staged base case for the `read_sparse_ijk()` IBM test
- `integration/directshortwave/test_directshortwave.py`
  - `100/` for the legacy Fortran vs f2py vs Python no-tree reference test
  - `525/` for the tree-enabled `facsec` vs `moller` comparison
- `regression/new_vegetation_module_against_v2.2/run_test.py`
  - `526/` as the staged base case for the reduced legacy vegetation regression

Guidance:

- Put reusable committed fixtures here when more than one test may need them.
- Keep test-specific scripts and assertions under `tests/integration/`,
  `tests/regression/`, `tests/system/`, or `tests/unit/`.
- Keep each consumer suite responsible for documenting its own runtime class:
  baseline CI, CI integration, CI regression, or heavy/manual.
- Do not store generated solver outputs or local binaries here.
