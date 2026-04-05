---
name: udales-exec
description: Run uDALES code and tests across heterogeneous clusters; detect the environment, load the right modules, and choose the correct build/run/test workflow.
---

# Cluster Runbook Skill

This skill helps Codex run uDALES on multiple machines with different module stacks and toolchains. It prioritizes the repo wrappers and the documented workflows, and it records new cluster-specific instructions as they are discovered.

## When to use

Use this skill when the user asks to build, run, or test uDALES on a cluster or any system where the environment varies (modules, MPI launcher, Python stack).

## Note on multi-system use

This skill is designed so an LLM can run uDALES quickly on different systems
despite local tooling quirks. It relies on the existing wrapper scripts and
records new environment details in `references/clusters.md` rather than
changing tooling in place.

## Quick start

1. Identify the system.
2. Load the correct modules.
3. Choose the correct path for:
   - solver build
   - preprocessing build
   - Python tooling/tests
   - MPI/executable tests
4. Record any new cluster-specific steps in `references/clusters.md`.

## Identify the system

Prefer a lightweight probe, in this order:

1. `hostname -f`
2. `uname -a`
3. `module --version` (if modules exist)
4. `command -v mpiexec` and `mpiexec --version`
5. `command -v python3` and `python3 --version`

If a known hostname or module stack matches a known cluster in `references/clusters.md`, use that profile directly. If not, treat it as a new cluster.

## Default repo workflows (use these unless the cluster profile says otherwise)

Solver build:

```bash
./tools/build_executable.sh <system> <debug|release>
```

Preprocessing build (View3D etc.):

```bash
./tools/build_preprocessing.sh <system>
```

Supported tests:

```bash
python tests/run_tests.py supported --branch-a <branch_a> --branch-b <branch_b> --build-type <Debug|Release>
```

MPI integration tests (solver):

```bash
bash tests/integration/ibm_sparse_input/run_test.sh
bash tests/integration/mpi_operators/run_test.sh
python tests/integration/processor_boundaries/test_processor_boundaries.py
```

Python tests:

```bash
python -m unittest discover -s tools/python/tests -p 'test_*.py'
```

## When the cluster is unknown

1. Probe modules and toolchain:
   - `module avail` for `intel`, `gcc`, `mpi`, `netcdf`, `fftw`, `cmake`
   - record exact module names and versions
2. Probe Python:
   - if `Python/3.9.*` module exists, prefer it
3. Attempt a minimal solver build using the repo wrapper.
4. Update `references/clusters.md` with:
   - hostname patterns
   - module load lines
   - build command(s)
   - test command(s)

## Self-update rule

Whenever a new system is encountered or a module stack changes, append a new
entry to `references/clusters.md`. Do not overwrite existing entries unless
the user explicitly approves a cleanup.

## References to consult

Read these files when you need guidance:

- `docs/cluster_workflows.md`
- `tests/README.md`
- `tests/test_suites.yml`
- `tools/build_executable.sh`
- `tools/build_preprocessing.sh`

## Scripts

Use `scripts/detect_env.sh` to gather a consistent fingerprint for unknown
clusters before writing a new profile entry.

Use `scripts/skill_selftest.sh` to run a lightweight check that the repo
wrappers and integration scripts exist, the environment probe runs, and the
Python unit tests execute. Pass `--full` to run execution-side checks only
(no builds).

The Python self-test uses `UD_VENV` and `UD_PYTHON` if set. Example:

```bash
UD_VENV=/rds/general/user/$USER/home/udales/.venv
UD_PYTHON=python3
UD_VENV="$UD_VENV" UD_PYTHON="$UD_PYTHON" \
  skills/udales-exec/scripts/skill_selftest.sh
```

If `UD_VENV` is not set, it falls back to:
`/rds/general/user/$USER/home/udales/.venv`.

For execution-only validation with a prebuilt solver binary:

```bash
UD_VENV=/rds/general/user/$USER/home/udales/.venv
UD_PYTHON=python3
UD_BUILD=/path/to/u-dales
UD_VENV="$UD_VENV" UD_PYTHON="$UD_PYTHON" UD_BUILD="$UD_BUILD" \
  skills/udales-exec/scripts/skill_selftest.sh --full
```
