# `udgeom` MATLAB Reference Data

This directory contains committed MATLAB reference outputs for the offline
Python `UDGeom` unit tests.

Contents:
- small deterministic STL fixtures
- one real geometry fixture: `xie_castro_2008`
- JSON files harvested from the MATLAB `udgeom` class
- test-owned MATLAB harvest helper: `harvest_udgeom_reference.m`

Regeneration:
```bash
module load tools/prod MATLAB/2024b Python/3.9.6-GCCcore-11.2.0
source /rds/general/user/mvr/home/udales/.venv/bin/activate
python tools/harvest_udgeom_matlab_references.py --clean
```

Notes:
- normal Python unit tests consume the committed JSON and do not require MATLAB
- `flat_ground` uses a fallback in the MATLAB harvest helper because MATLAB
  `splitBuildings` errors when ground removal leaves no triangles; the harvested
  building-derived outputs for that fixture therefore represent the effective
  no-buildings behavior rather than a clean MATLAB success path
