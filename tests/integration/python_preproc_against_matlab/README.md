# Python Preprocessing Against MATLAB

This directory contains a preprocessing parity test between the MATLAB and
Python preprocessing entry points.

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

To run it:

```bash
source ../.venv/bin/activate
python tests/integration/python_preproc_against_matlab/test_python_preproc_against_matlab.py
```
