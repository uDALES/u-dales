# Integration Tests

This directory contains end-to-end tests that exercise multiple interacting
components of the model and its tooling.

Current contents:

- `sparse_ijk/` for MPI validation of the sparse IBM input reader
- `tree_input/` for vegetation/tree input assets used by solver and Python tests

These tests are branch-specific and complement, rather than replace, the
input-interface integration tests present on other branches.

`tree_input/` currently acts mainly as a shared integration fixture rather than
as a standalone automated test entry point. That is still appropriate here
because the case is used by both repo-level and Python-level checks.
