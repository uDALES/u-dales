# Comparison Scripts

This directory contains standalone Python scripts for comparing u-DALES input and output files, and a shell script for setting up the required Python environment.

For documentation on running the automated test suite, see [`tests/system/README_test_system.md`](../tests/system/README_test_system.md).

---

## Directory layout

```
tools/
├── ud_set_nc_venv.sh              # One-time setup of the Python NetCDF environment
├── ud_compare_inputs.py           # Compare input files for a single case pair
├── ud_compare_outputs.py          # Compare NetCDF output files for a single case pair
├── ud_compare_multiple_inputs.py  # Compare input files across N cases (all pairs)
└── ud_compare_multiple_outputs.py # Compare output files across N cases (all pairs)
```

---

## Setting up the Python environment

`ud_compare_inputs.py` and `ud_compare_outputs.py` require `netCDF4` and `numpy`.
Run `ud_set_nc_venv.sh` once to create an isolated virtual environment with all required dependencies.

### Requirements

- Python 3.9 or newer (the script checks and reports if your version is too old)

### Setup (run once)

```bash
./ud_set_nc_venv.sh
```

The script will:
1. Can be run from any location with its full path.
2. Create a virtual environment at `tools/.venv_netcdf/` (always next to the script).
3. Install `numpy` and `netCDF4`.
4. Verify all imports work.
5. Deactivate the environment.

If `.venv_netcdf/` already exists, you will be asked whether to recreate it.
Answering **N** skips installation and runs the import checks on the existing environment instead.

### Recreating the environment

```bash
rm -rf .venv_netcdf
./ud_set_nc_venv.sh
```

---

## `ud_compare_outputs.py`

Compares four NetCDF output files between two experiment directories.

### Signature

```
ud_compare_outputs.py <exp_num> <exppath> <ref_data_path> [tolerance] [tol_thl]
```

| Argument | Description |
|----------|-------------|
| `exp_num` | Experiment number, 1–999 (e.g. `224`) |
| `exppath` | Parent outputs directory; the `<exp_num>` subdirectory is appended automatically |
| `ref_data_path` | Parent directory containing reference output cases |
| `tolerance` | Max absolute error (default: `1e-6`) |
| `tol_thl` | Tolerance for temperature variables (default: `1e-6`) |

### What it checks

| File | Variables checked |
|------|------------------|
| `xytdump.<case>.nc` | Horizontal-plane averages: velocities, TKE, Reynolds stresses, temperature, SGS fluxes |
| `tdump.<case>.nc` | Vertical profiles: velocities, variances, scalars, SGS variables |
| `fielddump.<case>.nc` | Full 3-D instantaneous fields: u, v, w, thl, qt, sca1 |
| `treedump.<case>.nc` | Tree-resolved fields: tr_u, tr_v, tr_w, tr_thl, tr_qt, tr_qtR, tr_qtA, tr_sv1, tr_sv2, tr_omega |

For each variable the maximum absolute difference between the two files is computed. A test **passes** if that value is ≤ the tolerance. Temperature-related variables (names containing `thl`) use `tol_thl`; all other variables use `tolerance`.

- If a dump file is absent from **both** directories it is skipped with `[SKIP]`.
- If a dump file is present in one directory only it is reported as `[FAIL]`.
- The same three-way logic applies at the variable level within each file.

### Examples

```bash
VENV_PYTHON=.venv_netcdf/bin/python3

# Default tolerance (1e-6)
$VENV_PYTHON ud_compare_outputs.py 224 tests/system/outputs/ tests/system/ref_data/

# Custom tolerance
$VENV_PYTHON ud_compare_outputs.py 224 tests/system/outputs/ tests/system/ref_data/ 1e-8

# Separate tolerances for general and temperature variables
$VENV_PYTHON ud_compare_outputs.py 224 tests/system/outputs/ tests/system/ref_data/ 1e-8 1e-7
```

---

## `ud_compare_inputs.py`

Compares input files between two experiment directories using text, numeric, or NetCDF methods depending on the file type.

### Signature

```
ud_compare_inputs.py <exp_num> <exppath> <ref_path> [tolerance]
```

| Argument | Description |
|----------|-------------|
| `exp_num` | Experiment number, 1–999 |
| `exppath` | Parent experiments directory; the `<exp_num>` subdirectory is appended automatically |
| `ref_path` | Parent reference directory; the `<exp_num>` subdirectory is appended automatically |
| `tolerance` | Max absolute error for numeric/NetCDF files (default: `1e-10`) |

### Comparison methods

| Method | Files |
|--------|-------|
| **Exact text** | `namoptions.<case>` |
| **Numeric** (max absolute error, `#` comments skipped) | `prof.inp`, `lscale.inp`, `scalar.inp`, `scalarsourcel.inp.N`, `facetarea.inp`, `facets.inp`, `facets_unused`, `factypes.inp`, `svf.inp`, `netsw.inp`, `Tfacinit.inp`, `heatpump.inp`, `trees.inp`, `Sdir.txt`, `facet_sections_*.txt`, `fluid_boundary_*.txt`, `solid_*.txt` |
| **NetCDF** (all numeric variables, max absolute error) | `vf.nc.inp.<case>` |

- If a file is absent from **both** directories it is skipped with `[SKIP]`.
- If a file is present in one directory but not the other it is reported as `[FAIL]`.

### Examples

```bash
VENV_PYTHON=.venv_netcdf/bin/python3

# Default tolerance (1e-10)
$VENV_PYTHON ud_compare_inputs.py 224 tests/system/experiments/ tests/system/ref_data/

# Custom tolerance
$VENV_PYTHON ud_compare_inputs.py 224 tests/system/experiments/ tests/system/ref_data/ 1e-8
```

---

## `ud_compare_multiple_outputs.py`

Compares NetCDF output files across N experiment cases (all pairwise combinations). Writes full variable-by-variable details to a timestamped log file; the terminal shows only a one-line result per pair plus an overall summary.

### Signature

```
ud_compare_multiple_outputs.py <exppath1> <exp1> <exppath2> <exp2> [exppath3 exp3 ...] [tolerance]
```

| Argument | Description |
|----------|-------------|
| `exppathN` | Parent outputs directory for case N |
| `expN` | Experiment number for case N (integer 1–999) |
| `tolerance` | Max absolute error (default: `1e-6`); detected as the last argument if it parses as a float |

At least two cases (4 arguments) are required.

### Examples

```bash
VENV_PYTHON=.venv_netcdf/bin/python3

# Compare two output directories for exp 224
$VENV_PYTHON ud_compare_multiple_outputs.py tests/system/outputs/ 224 tests/system/ref_data/ 224

# Compare three directories (three pairs)
$VENV_PYTHON ud_compare_multiple_outputs.py path_a/outputs/ 100 path_b/outputs/ 100 path_c/outputs/ 100

# Custom tolerance
$VENV_PYTHON ud_compare_multiple_outputs.py tests/system/outputs/ 224 tests/system/ref_data/ 224 1e-8
```

The log file `ud_compare_multiple_outputs_YYYYMMDD_HHMMSS.log` is written to the current working directory.

---

## `ud_compare_multiple_inputs.py`

Compares input files across N experiment cases (all pairwise combinations). Writes full variable-by-variable details to a timestamped log file; the terminal shows only a one-line result per pair plus an overall summary.

### Signature

```
ud_compare_multiple_inputs.py <exppath1> <exp1> <exppath2> <exp2> [exppath3 exp3 ...] [tolerance]
```

| Argument | Description |
|----------|-------------|
| `exppathN` | Parent experiments directory for case N |
| `expN` | Experiment number for case N (integer 1–999) |
| `tolerance` | Max absolute error (default: `1e-10`); detected as the last argument if it parses as a float |

At least two cases (4 arguments) are required.

### Examples

```bash
VENV_PYTHON=.venv_netcdf/bin/python3

# Compare two input directories for exp 224
$VENV_PYTHON ud_compare_multiple_inputs.py tests/system/experiments/ 224 tests/system/ref_data/ 224

# Compare three directories (three pairs)
$VENV_PYTHON ud_compare_multiple_inputs.py path_a/experiments/ 100 path_b/experiments/ 100 path_c/ 100

# Custom tolerance
$VENV_PYTHON ud_compare_multiple_inputs.py tests/system/experiments/ 224 tests/system/ref_data/ 224 1e-8
```

The log file `ud_compare_multiple_inputs_YYYYMMDD_HHMMSS.log` is written to the current working directory.
