# Python Virtual Environment Setup

The Linux setup script creates a repo-local virtual environment at `.venv/`,
installs the Python dependencies, builds the mandatory
`directshortwave_f2py` wrapper, and builds `tools/View3D`.

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
- builds `tools/python/udprep/directshortwave_f2py*.so`
- builds `tools/View3D/build/src/view3d`

## Manual Rebuilds

Rebuild the f2py wrapper:

```bash
source .venv/bin/activate
./tools/python/fortran/build_f2py.sh
```

Rebuild View3D:

```bash
./tools/build_preprocessing.sh common
```

## Notes

- `directshortwave_f2py` and View3D are required by the Python radiation
  workflow; they are not optional extras.
- `requirements-build.txt` covers Python-side build tools only. A system
  Fortran compiler such as `gfortran` must also be available.
