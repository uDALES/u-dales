# Integration Tests

This directory contains end-to-end tests that exercise multiple interacting
components of the model and its tooling.

Current contents:

- `directshortwave/` for committed direct shortwave cases, including no-tree reference case `100` and tree case `525`
- `ibm_sparse_input/` for MPI validation of the sparse IBM input reader on committed case `101`
- `tree_sparse_compare/` for MPI validation of sparse tree forcing on committed case `525`

These tests are branch-specific and complement, rather than replace, the
input-interface integration tests present on other branches.

`directshortwave/` is Python-driven, but it belongs here because it is anchored
to committed repo fixtures in `tests/cases/` and validates agreement across
multiple implementations rather than a single Python API.

`ibm_sparse_input/` and `tree_sparse_compare/` are executable-driven solver tests. They
stage shared inputs from `tests/cases/` into temporary run directories before
execution.

Run it with:

```bash
source .venv/bin/activate
python tests/integration/directshortwave/test_directshortwave.py
```
