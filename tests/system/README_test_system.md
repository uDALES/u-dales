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
├── ud_set_nc_venv.sh         # One-time setup of the Python NetCDF environment
├── ud_compare_outputs.py     # Output comparison script (called by ud_test_sim.sh)
├── ud_compare_inputs.py      # Input file comparison script (called by ud_test_preprocess.sh)
├── ud_compare_multiple_outputs.py
└── ud_compare_multiple_inputs.py
```

---

## Quick start

### Step 1 — Set up the Python virtual environment (once)

`ud_compare_outputs.py` and `ud_compare_inputs.py` require `netCDF4` and `numpy`. Run `ud_set_nc_venv.sh` once to create an isolated Python virtual environment with all required dependencies.

#### Requirements

- Python 3.9 or newer (the script will check and tell you if your version is too old)

#### Setup (run once)

```bash
../../tools/ud_set_nc_venv.sh
```

The script will:
1. Can be run from any location with full path.
2. Create a virtual environment at `tools/.venv_netcdf/` (relative to the repo root).
3. Install `numpy` and `netCDF4`.
4. Verify all imports work.
5. Deactivate the environment.

If `.venv_netcdf/` already exists, you will be asked whether to recreate it.
Answering **N** skips installation and runs the import checks on the existing environment instead.

#### Recreating the NetCDF environment

```bash
rm -rf ../../tools/.venv_netcdf
../../tools/ud_set_nc_venv.sh
```

#### Notes

- Virtual environment is reused across multiple runs.
- Dependencies only installed once (first run).
- Environment persists until manually deleted.

---

### Step 2a — Run the simulation test suite

Run from the `tests/system/` directory:

```bash
# CPU build — default cases (224)
./ud_test_sim.sh ref_data/

# CPU build — explicit
./ud_test_sim.sh ref_data/ --system common

# GPU build — default cases (402, 502, 452, 410, 411)
./ud_test_sim.sh ref_data/ --system gpu

# Override which cases to test
./ud_test_sim.sh ref_data/ 100 201
./ud_test_sim.sh ref_data/ 402 452 --system gpu

# Custom tolerance
./ud_test_sim.sh ref_data/ 224 --tolerance 1e-8

# Separate tolerances for general and temperature variables
./ud_test_sim.sh ref_data/ 224 --tolerance 1e-8 --tol-thl 1e-7
```

Progress and results are written to `logdir/test_common.log` or `logdir/test_gpu.log`.

`ud_test_sim.sh` detects `tools/.venv_netcdf/` automatically and uses it for comparisons.

---

### Step 2b — Run the input file test suite

`ud_test_preprocess.sh` regenerates input files from scratch and compares them against a reference directory.

```bash
# Default cases against ref_data/
./ud_test_preprocess.sh ref_data/

# Specific cases
./ud_test_preprocess.sh ref_data/ 100 224

# Custom tolerance for numeric files
./ud_test_preprocess.sh ref_data/ 224 --tolerance 1e-8
```

Progress and results are written to `test_inputs.log`.

---

## What `ud_test_sim.sh` does

`ud_test_sim.sh` runs three phases per invocation, with an optional diagnostic phase:

| Phase | What happens |
|-------|-------------|
| **1 — Build** | Deletes `build/` and rebuilds u-dales (`tools/build_executable.sh common release` for `--system common`; loads `nvhpc/24.11` and uses `tools/build_executable.sh gpu release` for `--system gpu`) |
| **2 — Simulate** | For each case, deletes any previous `outputs/<case>/` directory then runs the simulation via `tools/local_execute.sh tests/system/experiments/<case>` |
| **3 — Compare outputs** | Calls `tools/ud_compare_outputs.py <case> outputs/ <ref_data_path> <tolerance> <tol_thl>`, comparing freshly produced NetCDF files in `outputs/<case>/` against `<ref_data_path>/<case>/`. `--tolerance` applies to all variables; `--tol-thl` applies to temperature variables (defaults to `--tolerance` if not set) |
| **4 — Compare inputs** *(optional)* | Runs automatically if Phase 3 fails. Calls `tools/ud_compare_inputs.py <case> experiments/ <ref_data_path>` to check whether input files differ from the reference, helping diagnose the root cause |

The script exits with code `0` (all passed) or `1` (any failure).

---

## What `ud_test_preprocess.sh` does

`ud_test_preprocess.sh` runs two phases per case:

| Phase | What happens |
|-------|-------------|
| **1 — Write inputs** | Cleans the experiment directory (keeping geometry, `namoptions`, `config`, `trees.inp`, `heatpump.inp`, and `scalarsource*`), then regenerates all other input files by calling `tools/write_inputs.sh experiments/<case>` |
| **2 — Compare inputs** | Calls `tools/ud_compare_inputs.py <case> experiments/ <ref_data_path>`, comparing all generated input files against the reference |

The script exits with code `0` (all passed) or `1` (any failure).

---

## Adding a new test case

1. Create `experiments/<NNN>/` with the simulation input files and a `config.sh`
   that exports `DA_EXPDIR`, `DA_BUILD`, `DA_WORKDIR`, and `NCPU`.
2. Run the simulation once manually to produce output, then copy the NetCDF
   output files (`xytdump.<NNN>.nc`, `tdump.<NNN>.nc`, `fielddump.<NNN>.nc`,
   and `treedump.<NNN>.nc` if applicable) to `ref_data/<NNN>/`.
3. Optionally also copy the generated input files to `ref_data/<NNN>/` to
   enable input comparison.
4. Pass the case number on the command line:
   ```bash
   ./ud_test_sim.sh ref_data/ <NNN>
   ./ud_test_preprocess.sh ref_data/ <NNN>
   ```
