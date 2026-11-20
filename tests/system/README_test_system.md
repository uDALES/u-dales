# Full Simulation Test Suite

This directory contains scripts to build u-dales, run test simulations, compare their outputs against stored reference results, and verify input files.

For documentation on the comparison and venv-setup scripts themselves, see [`tools/README_compare_cases.md`](../../tools/README_compare_cases.md).

---

## Directory layout

```
tests/system/
├── ud_test_sim.sh            # Output file comparison entry point — build, simulate, compare outputs
├── ud_test_preprocess.sh     # Input file comparison entry point — write inputs, compare against reference
├── experiments/              # Input files for each test case
│   ├── 001/                  # Case directories
|   ├── ...
│   └── 999/
├── ref_data/                 # Reference data to compare against
│   ├── 001/                  # Case directories
|   ├── ...
│   └── 999/
├── outputs/                  # Simulation outputs written here at runtime
│   ├── 001/                  # Case directories
|   ├── ...
│   └── 999/

tools/                        # (repo root)
├── python/
    ├── setup_venv.sh         # One-time setup of the Python tooling environment
├── ud_compare_outputs.py     # Output comparison script (called by ud_test_sim.sh)
├── ud_compare_inputs.py      # Input file comparison script (called by ud_test_preprocess.sh)
├── ud_compare_multiple_outputs.py
└── ud_compare_multiple_inputs.py
```

---

## Quick start

### Step 1 — Set up the Python virtual environment (once)

`ud_compare_outputs.py` and `ud_compare_inputs.py` require `netCDF4` and `numpy`, and the preprocessing tests use the same Python tooling environment. Run `tools/python/setup_venv.sh common` once to create an isolated virtual environment with all required dependencies.
Use `icl` instead of `common` on the Imperial HPC cluster and `ecmwf` on the ECMWF Atos HPC2020.

#### Requirements

- Python 3.10 or newer (the script will check and tell you if your version is too old)

#### Setup (run once from the repo root)

```bash
# For a local machine
bash tools/python/setup_venv.sh common

# For the Imperial HPC machine
bash tools/python/setup_venv.sh icl

# For the ECMWF Atos HPC2020
bash tools/python/setup_venv.sh ecmwf
```

The script will:
1. Run from any location when called with its full path.
2. Create a virtual environment at `tools/python/.venv/` (relative to the repo root).
3. Install the Python tooling dependencies, including `numpy` and `netCDF4`.
4. Verify the runtime imports work.
5. Build and verify the preprocessing tools requested by the setup target.

If `tools/python/.venv/` already exists, you will be asked whether to recreate it.
Answering **N** skips installation/rebuild and runs validation checks on the existing environment instead.

#### Recreating the Python environment (run from the repo root)

```bash
rm -rf tools/python/.venv

# For a local machine
bash tools/python/setup_venv.sh common

# For the Imperial HPC machine
bash tools/python/setup_venv.sh icl

# For the ECMWF Atos HPC2020
bash tools/python/setup_venv.sh ecmwf
```

#### Notes

- Virtual environment is reused across multiple runs.
- Dependencies only installed once (first run).
- Environment persists until manually deleted.

---

### Step 2a — Run the simulation test suite

Run from the `tests/system/` directory:

```bash
# CPU build — default cases (100, 218, 224, 242, 807)
./ud_test_sim.sh ref_data/

# CPU build — explicit
./ud_test_sim.sh ref_data/ --system common

# GPU build — default cases (402, 502, 452, 410, 411)
./ud_test_sim.sh ref_data/ --system gpu

# Override which cases to test
./ud_test_sim.sh 100 201 ref_data/
./ud_test_sim.sh 402 452 ref_data/ --system gpu

# Custom tolerance for all variables
./ud_test_sim.sh 224 ref_data/ --tolerance 1e-8

# Separate tolerances for general and temperature variables
./ud_test_sim.sh 224 ref_data/ --tolerance 1e-8 --tol-thl 1e-7
```

When specific cases are supplied, place them before `ref_data/`. The last positional argument is always the reference data path.

Progress and results are written to `logdir/test_common.log` or `logdir/test_gpu.log`.

`ud_test_sim.sh` detects `tools/python/.venv/` automatically and uses it for comparisons.

---

### Step 2b — Run the input file test suite

`ud_test_preprocess.sh` regenerates input files from scratch and compares them against a reference directory.

```bash
# Default cases against ref_data/ with the MATLAB preprocessing route
./ud_test_preprocess.sh -m ref_data/

# Generate inputs with the Python preprocessing route
./ud_test_preprocess.sh -p ref_data/

# Specific cases
./ud_test_preprocess.sh -m 100 224 ref_data/

# Custom tolerance for numeric files
./ud_test_preprocess.sh -m 224 ref_data/ --tolerance 1e-8

# On the Imperial HPC cluster (loads the icl modules for the build)
./ud_test_preprocess.sh -p 990 991 ref_data/ --system icl

# On the ECMWF Atos HPC2020: the whole suite is submitted as one Slurm job
./ud_test_preprocess.sh -p 990 991 /path/to/ref_data ecmwf
```

When specific cases are supplied, place them before `ref_data/`. The last positional argument is always the reference data path. The system can be given as `--system <common|icl|ecmwf>` or as the bare word, as for `build_executable.sh`; the default is `common`.

How the suite runs depends on the system:

- `common` and `icl`: the suite detaches from the terminal, so it keeps running when the terminal is closed. Everything that would have appeared on the terminal goes to `logdir/test_inputs.console.log` (follow it with `tail -f`); the detailed log is `logdir/test_inputs.log`. Pass `--foreground` to run in the current shell instead. When the output is already redirected (not a terminal) the suite runs in place.
- `ecmwf`: the suite is submitted as one Slurm job on `qos nf` that builds the preprocessing tools, writes the inputs of every case and compares them. The job output is `logdir/test_inputs.console.log`. The job uses `PREPROC_ACCOUNT` (default `spnlthee`), `PREPROC_QOS` (`nf`), `PREPROC_WALLTIME` (`04:00:00`) and `PREPROC_MEM` (`32gb`) from the environment, and asks for as many CPUs as the largest `nompthreads` among the selected cases (at least 8). Only the Python route is available on ECMWF.

Within the suite each case is preprocessed with `write_inputs.sh ... f <system>` (foreground), so a failure is detected from the exit code as well as from the log.

---

## What `ud_test_sim.sh` does

`ud_test_sim.sh` runs three phases per invocation, with an optional diagnostic phase:

| Phase | What happens |
|-------|-------------|
| **1 — Build** | Deletes `build/` and rebuilds u-dales (`tools/build_executable.sh common release` for `--system common`; loads `nvhpc/24.11` and uses `tools/build_executable.sh gpu release` for `--system gpu`) |
| **2 — Simulate** | For each case, deletes any previous `outputs/<case>/` directory then runs the simulation via `tools/local_execute.sh tests/system/experiments/<case>` |
| **3 — Compare outputs** | Calls `tools/ud_compare_outputs.py <case> outputs/ <ref_data_path> <tolerance> <tol_thl>`, comparing freshly produced NetCDF files in `outputs/<case>/` against `<ref_data_path>/<case>/`. `--tolerance` applies to non-temperature variables and defaults to `1e-6`; `--tol-thl` applies to temperature variables and defaults to `--tolerance` if not set |
| **4 — Compare inputs** *(optional)* | Runs automatically if Phase 3 fails. Calls `tools/ud_compare_inputs.py <case> experiments/ <ref_data_path>` to check whether input files differ from the reference, helping diagnose the root cause |

The script exits with code `0` (all passed) or `1` (any failure).

---

## What `ud_test_preprocess.sh` does

`ud_test_preprocess.sh` runs two phases per case:

| Phase | What happens |
|-------|-------------|
| **1 — Write inputs** | Cleans the experiment directory (keeping geometry, `namoptions`, `config`, `trees.inp`, `heatpump.inp`, and `scalarsource*`), then regenerates all other input files by calling `tools/write_inputs.sh -m experiments/<case> f <system>` when `ud_test_preprocess.sh -m ...` is used or `tools/write_inputs.sh -p experiments/<case> f <system>` when `ud_test_preprocess.sh -p ...` is used; `f` runs the preprocessing in the foreground so its exit code is checked |
| **2 — Compare inputs** | Calls `tools/ud_compare_inputs.py <case> experiments/ <ref_data_path>`, comparing all generated input files against the reference |

The script exits with code `0` (all passed) or `1` (any failure).

---

## Adding a new test case

1. Create `experiments/<NNN>/` with the simulation input files and a `config.sh`
   that exports `DA_EXPDIR`, `DA_BUILD`, `DA_WORKDIR`, and `NCPU`.
2. Run the simulation once manually to produce output, then copy the NetCDF
   output files (`stats_xy.<NNN>.nc`, `stats_t.<rank>.<NNN>.nc`,
   `ins_field.<rank>.<NNN>.nc`, and `stats_tree.<rank>.<NNN>.nc` if
   applicable) to `ref_data/<NNN>/`.
3. Optionally also copy the generated input files to `ref_data/<NNN>/` to
   enable input comparison.
4. Pass the case number on the command line:
   ```bash
   ./ud_test_sim.sh <NNN> ref_data/
   ./ud_test_preprocess.sh -m <NNN> ref_data/
   ```
