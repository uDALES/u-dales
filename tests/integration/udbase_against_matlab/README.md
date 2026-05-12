# `UDBase` Against MATLAB

This integration suite compares selected Python `UDBase` behaviors against
committed reference outputs harvested once from the MATLAB `udbase` class.

It belongs under `tests/integration/` rather than `tools/python/tests/`
because it is anchored to real committed simulation directories under
`tests/cases/` and validates parity across two implementations.

Current coverage:
- `load_facsec('c')`
- `calculate_frontal_properties()`

Reference cases:
- `tests/cases/064`
- `tests/cases/101`

Normal test runs consume the committed JSON in `data/` and do not require
MATLAB. To regenerate the references:

```bash
module load tools/prod MATLAB/2024b Python/3.9.6-GCCcore-11.2.0
source /rds/general/user/mvr/home/udales/.venv/bin/activate
python tools/harvest_udbase_matlab_references.py --clean
```
