# uDPrep Integration Tests

This directory contains preprocessing integration tests that exercise the real
`udprep` stack against committed repo cases or external preprocessing tools.

Current contents:

- `test_view3d.py`: validates View3D execution and cached path handling through
  `UDPrep.radiation.calc_view_factors`
- `test_scanline_contract.py`: validates scanline input construction and legacy
  serialization against a committed case and the compiled `directshortwave_f2py`
  wrapper

These tests are heavier than pure Python unit tests under `tools/python/tests`
because they depend on real case data, compiled wrappers, or external
preprocessing executables.
