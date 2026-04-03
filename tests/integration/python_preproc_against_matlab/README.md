# Python Preprocessing Against MATLAB

This directory contains a preprocessing parity test between the MATLAB and
Python preprocessing entry points.

Reusable workflow logic now lives in:

- [compare_preprocessing.py](/rds/general/user/mvr/home/udales/u-dales/tools/compare_preprocessing.py)

The files in this directory are thin test and exploratory wrappers around that
tool.

Current scope:

- case `tests/cases/100/`
- no-tree preprocessing only
- intended to compare the generated preprocessing text files from:
  - MATLAB `write_inputs.m`
  - Python `tools/write_inputs.py`

Why it exists:

- it is a repo-level integration test, not a Python unit test
- it compares two end-to-end preprocessing implementations on the same case
- it avoids vegetation for now so the test isolates the core preprocessing path

Current status:

- this suite captures the intended parity contract
- the Python preprocessing chain is still being completed, so this test may
  fail until the Python path can regenerate the same IBM/profile outputs
- if MATLAB is not available, the test skips cleanly

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

To run it:

```bash
source ../.venv/bin/activate
python tests/integration/python_preproc_against_matlab/test_python_preproc_against_matlab.py
```

To run the reusable tool directly:

```bash
source ../.venv/bin/activate
python tools/compare_preprocessing.py --case-dir tests/cases/100
python tools/compare_preprocessing.py --list-cases
python tools/compare_preprocessing.py 064 --keep-temp
```
