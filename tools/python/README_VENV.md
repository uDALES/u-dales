# Python Virtual Environment Setup

The Linux setup script creates a repo-local virtual environment at `.venv/`,
installs the Python dependencies, and builds the standalone preprocessing
toolchain via `tools/build_preprocessing.sh`.

Python 3.9 or newer is required. Older interpreters such as Python 3.6 are not
supported by the current Python tooling.

## One-Time Setup

From the `u-dales` root:

```bash
./tools/python/setup_venv.sh
```

If your default `python3` does not ship the development headers required by
f2py, point the script at another interpreter:

```bash
PYTHON_BIN=/opt/pbs/python/bin/python3 ./tools/python/setup_venv.sh
```

## Daily Usage

Activate the environment from the repo root:

```bash
source .venv/bin/activate
```

Then run tooling or tests, for example:

```bash
python -m unittest discover -s tools/python/tests
python tools/write_inputs.py /path/to/case
```

When finished:

```bash
deactivate
```

## What The Setup Script Does

- creates `.venv/`
- installs `tools/python/requirements.txt`
- installs `tools/python/requirements-build.txt`
- builds `tools/preprocessing/build/bin/view3d`
- builds `tools/preprocessing/build/bin/IBM_preproc`
- builds `tools/python/udprep/directshortwave_f2py*.so`
- builds `tools/python/udprep/ibm_preproc_f2py*.so`

## Manual Rebuilds

Rebuild the full preprocessing toolchain:

```bash
./tools/build_preprocessing.sh common
```

Rebuild a single target:

```bash
./tools/build_preprocessing.sh common directshortwave_f2py
./tools/build_preprocessing.sh common ibm_preproc_f2py
./tools/build_preprocessing.sh common view3d
```

## Notes

- `directshortwave_f2py`, `ibm_preproc_f2py`, and View3D are required parts of
  the Python preprocessing workflow.
- `requirements-build.txt` covers Python-side build tools only. A system
  Fortran compiler such as `gfortran` must also be available.
