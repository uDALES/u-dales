# Fortran preprocessing wrappers

This folder contains the Python-facing Fortran sources, wrapper glue, and
`.pyf` signatures used to build the preprocessing extension modules.

The build entry point is the standalone preprocessing CMake project:

```bash
./tools/build_preprocessing.sh common
```

Or directly:

```bash
cmake -S tools/preprocessing -B tools/preprocessing/build \
  -DPREPROCESSING_PYTHON_EXECUTABLE=$(command -v python)
cmake --build tools/preprocessing/build --target preprocessing_tools
```

Current wrapped targets include:

- `directshortwave_f2py`
- `ibm_preproc_f2py`

The built extension modules are written to `tools/python/udprep/`.
