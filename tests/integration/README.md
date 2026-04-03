# Integration Tests

This directory contains end-to-end tests that exercise multiple interacting
components of the model and its tooling.

Integration tests are not all the same. In this repository they generally fall
into two groups:

- CI-eligible integration tests: moderate-cost checks that fit GitHub Actions
  and may depend on compilers, MPI, committed fixtures, or preprocessing
  binaries
- heavy/manual integration tests: larger or more resource-intensive checks that
  are intended for local or cluster execution rather than routine GitHub
  Actions runs

Current contents:

- `directshortwave/` for committed direct shortwave cases, including no-tree
  reference case `100` and tree case `525`
- `udprep/` for preprocessing integration checks that exercise real `UDPrep`
  workflows against committed cases or preprocessing binaries
- `ibm_sparse_input/` for MPI validation of the sparse IBM input reader on committed case `101`

`directshortwave/` is Python-driven, but it belongs here because it is anchored
to committed repo fixtures in `tests/cases/` and validates agreement across
multiple implementations rather than a single Python API. It should be treated
as a preprocessing integration suite, not as a pure solver test.

`udprep/` is also Python-driven, but it exercises end-to-end preprocessing
paths that depend on committed cases, compiled wrappers, or external
preprocessing executables. It should be treated as integration coverage rather
than as pure `tools/python` unit coverage.

`ibm_sparse_input/` is an executable-driven solver test. It
stage shared inputs from `tests/cases/` into temporary run directories before
execution.

Current supported status:

- `ibm_sparse_input/`: supported solver-facing integration coverage
- `udprep/`: supported preprocessing/tooling integration coverage
- `directshortwave/`: experimental preprocessing/tooling integration coverage

This directory now feeds two different curated paths through
`tests/test_suites.yml`:

- `supported`: solver-facing integration suites that are part of the supported
  path
- `experimental`: Python/preprocessing suites that are useful development
  coverage but not part of the primary supported merge gate

Each suite here should state whether it is intended for:

- the supported path
- an experimental, non-blocking lane
- local developer use
- cluster/manual execution only

At present, Python preprocessing/tooling integration tests should generally be
treated as experimental and non-blocking unless the interface is explicitly
supported and stable enough to become part of the required merge gate.

Direct shortwave prerequisites:

- Python environment with the project tooling installed
- compiled `directshortwave_f2py` wrapper
- `gfortran` for the legacy reference executable used by the no-tree case

Current status:

- useful development coverage for preprocessing
- not intended to define the primary supported GitHub Actions path yet
- should skip cleanly when prerequisites are absent rather than fail opaquely

Run the direct shortwave integration test with:

```bash
source ../.venv/bin/activate
python tests/integration/directshortwave/test_directshortwave.py
```

Or via the top-level dispatcher:

```bash
python tests/run_tests.py experimental
```

The `experimental` selection is defined in `tests/test_suites.yml`.

Supported solver-facing integration prerequisites:

- MPI launcher/runtime available in `PATH`
- built solver executable available to the integration harness
- project module stack loaded where required on the cluster

Run the supported integration path via:

```bash
python tests/run_tests.py supported --branch-a <branch_a> --branch-b <branch_b> --build-type <build_type>
```

In the current supported selection, the solver-facing integration suites are:

- `integration/ibm_sparse_input/run_test.sh`
- `integration/udprep/` preprocessing integration tests

Experimental integration coverage currently includes:

- `integration/directshortwave/test_directshortwave.py`
