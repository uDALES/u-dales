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

## Platform-specific binaries — rebuild when they don't match (S2)

The compiled `*_f2py` extensions are **platform- and interpreter-specific**
(e.g. a `cp312-win_amd64` build only loads under CPython 3.12 on 64-bit
Windows). On any other OS / Python version the import raises `RuntimeError`
with build instructions by design — rebuild for your platform:

```bash
./tools/build_preprocessing.sh common preprocessing_tools
```

A committed binary can also drift from the caller signature (e.g. a stale
zero-argument extension vs a ~20-argument call) — `UDPrep.run_all` then fails
until you rebuild. If IBM / direct-shortwave preprocessing errors on the f2py
call, rebuild the extensions before anything else.
