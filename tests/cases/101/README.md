# Case 101

Shared committed fixture for the IBM sparse-input integration test.

Current consumers:

- `tests/integration/ibm_sparse_input/run_test.sh`

This directory contains the base case inputs only. Rank-local golden files used
to validate `read_sparse_ijk()` stay in `tests/integration/ibm_sparse_input/`
because they are test-specific outputs, not reusable case data.
