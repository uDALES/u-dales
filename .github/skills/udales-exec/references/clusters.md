# Cluster Profiles

This file records cluster-specific build/run/test instructions discovered over time.
Append new entries rather than rewriting existing ones unless explicitly asked.

## Template

```
### <cluster-name>

Host patterns:
- <hostname or regex>

Module stack:
```bash
module load <modules>
```

## ICL (RDS/EPHEMERAL)

Host patterns:
- `*.imperial.ac.uk`
- `/rds/general/user/*/home`
- `/rds/general/user/*/ephemeral`

Module stack:
```bash
module load intel/2025a netCDF/4.9.2-iimpi-2023a netCDF-Fortran/4.6.1-iimpi-2023a FFTW/3.3.9-intel-2021a CMake/3.29.3-GCCcore-13.3.0 git/2.45.1-GCCcore-13.3.0
```

Solver build:
```bash
./tools/build_executable.sh icl debug
./tools/build_executable.sh icl release
```

Preprocessing build:
```bash
./tools/build_preprocessing.sh icl
```

Python environment:
```bash
module load Python/3.9.6-GCCcore-11.2.0
source /rds/general/user/mvr/home/udales/.venv/bin/activate
```

MPI launcher notes:
- prefer the module-provided `mpiexec`
- if Open MPI is used, set `TMPDIR=/tmp` and `--oversubscribe` for login-node runs
- interactive runs can be slow; prefer repo wrappers

Tests:
```bash
python tests/run_tests.py supported --branch-a <branch_a> --branch-b <branch_b> --build-type <Debug|Release>
bash tests/integration/mpi_operators/run_test.sh
python tests/integration/processor_boundaries/test_processor_boundaries.py
```

Solver build:
```bash
./tools/build_executable.sh <system> <debug|release>
```

Preprocessing build:
```bash
./tools/build_preprocessing.sh <system>
```

Python environment:
```bash
module load Python/<version>
source /path/to/.venv/bin/activate
```

MPI launcher notes:
- <mpiexec path or required extra args>

Tests:
```bash
python tests/run_tests.py supported --branch-a <branch_a> --branch-b <branch_b> --build-type <Debug|Release>
bash tests/integration/mpi_operators/run_test.sh
```
```
