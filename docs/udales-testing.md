# Testing and CI

This page is a developer-facing guide to how uDALES is tested, how to run
those tests locally, and what continuous integration (CI) checks on a pull
request. For the full test philosophy, layer definitions, and manifest
schema, see [`tests/README.md`](https://github.com/uDALES/u-dales/blob/master/tests/README.md);
this page summarises the parts a contributor needs day to day.

## Test layout

Tests live in two places, split by what they validate rather than by
language:

- **`tests/`** — repo-level tests: solver/MPI-driven checks, case-fixture
  integration tests, and branch-to-branch regressions. Organised into
  `unit/`, `integration/`, `system/`, and `regression/` subdirectories, plus
  `cases/` for shared committed case fixtures (e.g. cases `100`, `101`,
  `525`, `526`).
- **`tools/python/tests/`** — unit tests for the Python package (`udbase`,
  `udgeom`, `udprep`, `udvis`, and related utilities), discovered by
  `unittest discover` on the `test_*.py` naming pattern.

Exploratory or plotting-heavy scripts that aren't stable automated tests
belong in `tools/python/examples/`, not either test tree. Directory location
alone doesn't decide what runs where or under what policy — that's
controlled by the manifest described below.

## Running tests locally

All curated test selections are dispatched through
[`tests/run_tests.py`](https://github.com/uDALES/u-dales/blob/master/tests/run_tests.py),
which reads groups from [`tests/test_suites.yml`](https://github.com/uDALES/u-dales/blob/master/tests/test_suites.yml).
Each suite in the manifest is tagged with a `class` (`supported` or
`experimental`), a `kind` (scope), a `component`, a `platform`, and a rough
`cost` tier; these are printed as the suite runs.

List the available selections by passing an invalid one (or reading the
manifest directly — there's no separate `--list` flag):

```bash
python tests/run_tests.py --help
```

Run a named selection:

```bash
python tests/run_tests.py python-library
python tests/run_tests.py supported --branch-a master --branch-b HEAD --build-type Release
python tests/run_tests.py experimental
```

- `python-library` — `tools/python` unit tests plus the Python-driven
  integration/reference suites (directshortwave, udprep integration,
  udbase-vs-MATLAB parity, Python-vs-MATLAB preprocessing parity). Needs the
  conda `udales` environment (or equivalent) and `tools/python` installed
  (`pip install -e tools/python --no-deps`); no compiled solver, no MPI.
- `supported` — the curated merge-gating selection: `python-library` plus
  the branch-comparison regression harness and the solver/MPI-driven
  integration suites (IBM sparse input, MPI operators, processor
  boundaries). These need a compiled uDALES build; the MPI suites read the
  build path from `UDALES_BUILD` (set per suite by the manifest as
  `build/<build-type-lower>/u-dales`) and expect `mpiexec` on `PATH`. The
  regression suite additionally needs `--branch-a`/`--branch-b` (default
  `master`/`HEAD`) and builds both branches itself.
- `supported-macos` — a temporary macOS compatibility selection:
  `python-library` plus only the IBM sparse-input suite, excluding the
  branch-comparison regression until an unrelated macOS/Homebrew CMake
  incompatibility on `master` is fixed.
- `experimental` — coverage not yet in the merge gate: slower
  directshortwave unit/periodic checks, the vegetation-module-vs-`v2.2.0`
  regression, and the MPI averaging regression; also solver/branch-build
  dependent.
- `all` — `supported` + `experimental`.

`--python` (or `UDALES_TEST_PYTHON`) selects the interpreter used to launch
child suites; by default it's whatever interpreter you invoked
`run_tests.py` with, so activate your conda environment or
`tools/python/.venv` first.

To run only the Python package tests directly, without the dispatcher, use
`python -m unittest discover -s tools/python/tests -p "test_*.py"`. This
needs `tools/python` on the Python path (`pip install -e tools/python
--no-deps`, as CI does) and `environment.yml`'s dependencies. The
visualisation tests (`test_udbase_vis.py`) additionally need
`PYVISTA_OFF_SCREEN=true` and a headless GL/Xvfb setup, and can pull in the
optional Plotly backend (`tools/python/requirements-plotly.txt`).

For suites that need a compiled solver (the solver/MPI-driven suites under
`supported`/`experimental`), build first — see the [development
notes](udales-development-notes.md) for a `Debug` build, or the
[installation guide](udales-installation.md) for a `Release` build.

## Continuous integration

`.github/workflows/ci.yml` runs on every push and pull request (a
`[skip ci]` marker in the head commit message skips it on `push`). It has
four jobs:

- **`build-and-supported-tests`** — matrix over `os: [ubuntu-latest,
  macos-latest]` × `build-type: [Debug, Release]`. Installs system
  dependencies (`gfortran`, MPI, NetCDF, Graphviz, FFTW) and the conda
  `udales` environment, builds the preprocessing tools and installs
  `tools/python` editable, builds uDALES with CMake, summarises compiler
  warnings (see below), then runs `tests/run_tests.py supported` — except on
  `macos-latest`, which currently runs `supported-macos` instead.
- **`python-viz`** — Linux-only; installs headless OpenGL/Xvfb and the
  optional Plotly backend, then runs `test_udbase_vis.py` off-screen with
  both the PyVista and Plotly backends via `xvfb-run`.
- **`docs`** — builds the MkDocs site (`mkdocs build --site-dir build/html`)
  and the FORD Fortran API docs (`ford ford.md`), uploads the combined HTML
  as an artifact.
- **`publish-docs`** — only on `push` to `master` (and not `[skip ci]`);
  downloads the `docs` artifact and publishes it to GitHub Pages.

Compiler warnings from each build are summarised by
`.github/scripts/summarise_warnings.sh` into the GitHub Actions step
summary: it counts warnings by `-W` class, flags known-benign classes
(`-Wcompare-reals`, `-Wunused-value`, reviewed in #334) versus everything
else as "review", and lists the actionable warning sites. This step is
reporting-only and never fails the build — CI doesn't pin compiler versions,
so warning sets legitimately vary across runner images.

`.github/workflows/ci-rerun-on-cancel.yml` watches for CI runs that GitHub
itself cancelled (infra/runner loss, not a genuine test failure) and
automatically re-triggers them, capped at 3 attempts, so contributors don't
need to manually re-run infra-flaky jobs.

## Expectations for contributions

Before requesting review, a pull request should have a green
`build-and-supported-tests` matrix (all four OS/build-type combinations) and
a green `python-viz` job — together these run the `supported` (or
`supported-macos`) test selection plus both visualisation backends. The
`docs` job should also build cleanly if you touched documentation. Compiler
warnings are reported but do not gate merging, so a new warning won't fail
CI, but check the step summary and avoid introducing avoidable ones.
`experimental` and `heavy` suites are not part of the required gate, but are
worth running locally if your change touches the areas they cover (see
`tests/test_suites.yml` for what each suite exercises).
