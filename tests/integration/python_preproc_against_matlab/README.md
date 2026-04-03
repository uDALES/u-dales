# Python Preprocessing Against MATLAB

This directory contains a preprocessing parity test between committed
MATLAB-generated reference outputs and the Python preprocessing entry point.

Reusable workflow logic now lives in:

- [compare_preprocessing.py](/rds/general/user/mvr/home/udales/u-dales/tools/compare_preprocessing.py)

The file in this directory is the submitted integration test. MATLAB is not
required when running it.

Current scope:

- cases:
  - `tests/cases/100/`
  - `tests/cases/064/`
- includes:
  - a no-tree core preprocessing reference case (`100`)
  - a small SEB/View3D preprocessing reference case (`064`)
- committed MATLAB-generated preprocessing outputs live directly in the case
  directory
- the test copies the case to a temp directory, cleans generated preprocessing
  outputs, reruns Python `tools/write_inputs.py`, and compares the regenerated
  files against the committed references

Why it exists:

- it is a repo-level integration test, not a Python unit test
- it compares Python preprocessing against a committed MATLAB reference contract
- it still avoids vegetation for now
- it now covers both the core preprocessing path and a small SEB/View3D path

Reference refresh:

- committed references can be refreshed with:
  - [harvest_preprocessing_reference.py](/rds/general/user/mvr/home/udales/u-dales/tools/harvest_preprocessing_reference.py)
- this developer tool runs MATLAB in a temp copy of the case and copies the
  generated outputs back into the case directory

Exploratory real-case note:

- the reusable tool supports real-case exploratory runs and uses slightly
  looser tolerances for `svf.inp.*`, `vfsparse.inp.*`, and `netsw.inp.*`
- this is not because the Python preprocessing physics differs
- on cluster case `065`, the Python path matches a direct deterministic
  `View3D` run exactly, and the Python and MATLAB `.vs3` inputs are
  byte-identical
- the remaining drift comes from the MATLAB-driven `View3D` system-call path:
  72 sparse view-factor entries differ by `1e-6`, which then produces
  `O(1e-4)` differences in `netsw`
- the strict submitted integration test for case `100` is unaffected; this
  note applies only to the temporary real-case sweep

To run the submitted test:

```bash
source ../.venv/bin/activate
python tests/integration/python_preproc_against_matlab/test_python_preproc_against_matlab.py
```

To refresh the committed MATLAB references for a case:

```bash
module load tools/prod MATLAB/2024b Python/3.9.6-GCCcore-11.2.0
source /rds/general/user/mvr/home/udales/.venv/bin/activate
python tools/harvest_preprocessing_reference.py --case-dir tests/cases/100
python tools/harvest_preprocessing_reference.py --case-dir tests/cases/064
```
