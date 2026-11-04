# Troubleshooting

This page collects common problems and their fixes, grouped by the stage of the [workflow](udales-workflow.md) where they tend to show up. If your problem isn't listed here, please open a [GitHub issue](https://github.com/uDALES/u-dales/issues).

## Build problems

**`Could NOT find NetCDF (missing: NETCDF_LIBRARIES NETCDF_INCLUDE_DIRS NETCDF_HAS_INTERFACES)`** — CMake could not locate the NetCDF-C and/or NetCDF-Fortran libraries and headers. On HPC systems in particular, CMake often cannot find NetCDF automatically; hint the install locations explicitly with `-DNETCDF_DIR=<path-to-netcdf-c> -DNETCDF_FORTRAN_DIR=<path-to-netcdf-fortran>` as shown in the [HPC build instructions](udales-installation.md#build-on-hpcs). On local systems, install the `libnetcdf-dev`/`libnetcdff-dev` (Linux) or `netcdf`/`netcdf-fortran` (macOS Homebrew) packages listed in [Installation](udales-installation.md#prerequisites).

**`Only GNU, Intel, and Cray Fortran compilers are supported`** — CMake found a Fortran compiler that isn't GNU (`gfortran`), Intel, or Cray, so it doesn't know which compile flags to apply (see the `elseif` chain in `CMakeLists.txt`). Set `FC`/`CC` (or `-DCMAKE_Fortran_COMPILER=...`) to point at a supported compiler before configuring, e.g. `FC=mpiifort cmake ...` as in the [HPC build example](udales-installation.md#build-on-hpcs).

**`Build type 'X' not supported.`** — `CMAKE_BUILD_TYPE` was set to something other than `Debug` or `Release`. Reconfigure with `-DCMAKE_BUILD_TYPE=Release` or `Debug` (see the [build options table](udales-installation.md#build-defaultsoptions)).

**CMake reports NMake generator errors, or `numpy.f2py`/build tools fail near a truncated OneDrive-style path** — both are native-Windows preprocessing build issues (wrong default generator, and paths containing spaces breaking f2py). These are already documented, with fixes, in [Windows setup for preprocessing libraries](udales-preprocessing-windows.md#troubleshooting) — use `-G Ninja` explicitly and build from a no-space `subst` alias path.

## Pre-processing problems

**`The preprocessing route and path to case/experiment folder must be set.`** / **`config.sh must be set inside <dir>`** — `write_inputs.sh` requires both a route flag (`-m` for MATLAB or `-p` for Python) and a case directory, and the case directory must contain a `config.sh` defining `DA_TOOLSDIR` and `DA_EXPDIR`. See the usage message and checks in `tools/write_inputs.sh` and the [pre-processing guide](udales-pre-processing.md) for a template `config.sh`.

**`Python virtual environment not found or not executable: <path>/.venv/bin/python`** — the Python pre-processing route (`write_inputs.sh -p`) expects a virtual environment at `$DA_TOOLSDIR/python/.venv`. Create it first with `bash $DA_TOOLSDIR/python/setup_venv.sh <common|icl>` as the script itself suggests, then rerun `write_inputs.sh -p`.

**`Unrecognised preprocessing route: X`** — `write_inputs.sh` only accepts `-m` (MATLAB) or `-p` (Python) as its first argument; any other value is rejected. See [Pre-processing](udales-pre-processing.md) for the two supported routes.

## Runtime problems

**`ERROR: Namoptions does not exist`** — the executable was started without a valid namoptions path, or the file isn't in the working directory the run script `cd`'d into. The run scripts (`local_execute.sh`, `hpc_execute.sh`, `archer_execute.sh`, `ecmwf_execute.sh`) start the executable from inside the copied output directory `$DA_WORKDIR/<exp>` with `namoptions.<exp>`, so the file must be named `namoptions.<exp>` and match the three-digit experiment number in the case path.

**`ERROR: Problem in namoptions <BLOCK>` with `iostat error: <n>`** — a namelist block (`RUN`, `DOMAIN`, `PHYSICS`, `BC`, `WALLS`, `TREES`, ...) failed to parse, usually because of a typo'd key, a missing value, or a stray/misplaced entry. Each namelist has its own check in `src/modstartup.f90` (and similar per-module checks, e.g. `NAMSUBGRID` in `src/modsubgrid.f90`, `NAMCHECKSIM` in `src/modchecksim.f90`); compare your `namoptions` block against the [input parameters overview](udales-namoptions-overview.md).

**`STOP ERROR IN NUMBER OF PROCESSORS` / `nprocx must divide itot!!!` / `nprocy must divide jtot!!!`** — `nprocx` and `nprocy` in `namoptions` must exactly divide `itot`/`jtot` (and `ktot` for `nprocy`) so the domain decomposes into equal pencils; see the [`nprocx`/`nprocy` constraints](udales-namoptions-overview.md#namelist-run) and the [parallelisation overview](udales-2decomp.md). This is distinct from the related 2DECOMP&FFT startup error **`Invalid 2D processor grid - nproc /= p_row*p_col`**, which fires when the number of MPI ranks you launched with doesn't equal `nprocx * nprocy`; you must set `nprocx * nprocy = NCPU` (or `NNODE * NCPU` on HPC), see [Running uDALES](udales-simulation-setup.md). The run scripts do not read `namoptions`, so they cannot catch a mismatch: it surfaces here, at runtime.

**`Number of CPU cores NCPU must be set inside <dir>/config.sh`** (the ECMWF script says `Number of MPI ranks on each node NCPU must be set…`) — the run scripts check for `NCPU` (and, on HPC, `NNODE`) before launching `mpiexec`/`mpirun`/`srun`; set these in `config.sh` so their product matches `nprocx * nprocy`, per the note above.

**`Executable DA_BUILD=<path> does not exist or is not executable`** — `DA_BUILD` in `config.sh` points at a file that is not there (a typo'd path, or a build that has not finished). The run scripts check this before copying anything or submitting a job. Build first, e.g. `tools/build_executable.sh <system> release`, and check `DA_BUILD` matches the executable produced (`u-dales/bin/u-dales`, or `u-dales/build/<release|debug>/u-dales`).

**`usage: ud_run <machine> sim|gather <case_directory>` / `Unknown machine '<name>'`** — `ud_run` needs exactly three arguments, and `<machine>` must be one of `icl`, `hx1`, `archer`, `ecmwf` or `common` (the names used by `tools/build_executable.sh`). See [Running uDALES](udales-simulation-setup.md).

**`DA_WORKDIR=<dir> is not under $SCRATCH`** (ECMWF) — `ecmwf_execute.sh` only submits runs whose `DA_WORKDIR` is under `$SCRATCH`, the only ECMWF filesystem sized for large outputs. Point `DA_WORKDIR` in `config.sh` at a directory under `$SCRATCH`, see [Run on ECMWF HPC2020](udales-simulation-setup.md#run-on-ecmwf-hpc2020).

**`Quality of service QOS must be set inside <dir>/config.sh`** (ECMWF) — `ecmwf_execute.sh` has no default queue. Set `QOS="np"` for real runs (exclusive nodes, all 240 GB) or `QOS="nf"` for small test runs (shared nodes, 1 node, at most 128 cores / 128 GB). The gather job does not use `QOS`: it always runs on `nf`.

**`Memory requirement MEM must be set inside <dir>/config.sh for QOS=<qos>`** (ECMWF) — any queue other than `np` shares nodes, where a job that requests no memory gets the shared pool default of 8 GB. Set `MEM` in Slurm units, e.g. `MEM="16G"`. Simulations on `QOS=np` need no `MEM`, because an exclusive node comes with all of its 240 GB.

**`Paths must not contain spaces: ECMWF's sbatch cannot submit from such a directory`** (ECMWF) — ECMWF's `sbatch` is a site wrapper that splits its arguments on whitespace, so it cannot submit a job from a directory whose path contains a space. Use paths without spaces for both the case directory and `DA_WORKDIR`.

**`NCPU=<n> is more than the 128 physical cores of an ECMWF compute node`** (ECMWF) — `NCPU` is the number of MPI ranks per node; use more nodes (`NNODE`) instead of more than 128 ranks per node.

**`MEM=<value> is not a memory size`** (ECMWF) — `MEM` needs a number and a unit, e.g. `MEM="16G"` or the PBS spelling `MEM="16gb"`; both are accepted in any case and normalised for Slurm. A bare number such as `16` is rejected as ambiguous, as are forms like `16gib`.

**`ERROR: no restartfile set`** — `lwarmstart` or `lstratstart` was set to `.true.` in the `RUN` namelist without a corresponding `startfile`. Set `startfile` to the name of an existing restart file (`initdNNNNNNNN_xxx_xxx.000`), see the [`RUN` namelist](udales-namoptions-overview.md#namelist-run).

**`ERROR: invalid itree_mode. Supported values are 1 (drag only), 2 (sveg), 99 (legacy SEB).`** / **`ERROR: legacy tree SEB (itree_mode=99) cannot be combined with lEB=.true.`** — `ltrees=.true.` requires `itree_mode` to be one of the supported drag/SEB modes, and the legacy tree SEB mode (`99`) cannot be combined with the building energy balance (`lEB=.true.`). Check the `TREES` namelist against `src/modstartup.f90`.

**`ERROR: zgrid.inf does not exist`** — a driven simulation (`idriver = 2`) couldn't find the `zgrid.inf`/`zgrid.inl` inflow files written by the corresponding precursor run. Make sure the precursor (`idriver = 1`) has completed and `driverjobnr` points at the right experiment number, per [Driver simulations](udales-driver-simulations.md).

**Simulation stops partway through a driven run** — for `idriver = 2`, `runtime` must be `(driverstore-1)*dtdriver` seconds or less, where `driverstore`/`dtdriver` come from the precursor; exceeding this runs out of saved inflow data, per the note in [Driver simulations](udales-driver-simulations.md#running-driven-simulations).

**Simulation crashes with a floating-point exception (SIGFPE), or diverges/produces `NaN` output** — GNU builds are compiled with `-ffpe-trap=invalid,zero,overflow` (`CMakeLists.txt`), so an invalid operation, division by zero, or overflow anywhere aborts the run immediately rather than silently propagating `NaN`s. This is usually a sign the timestep is too large or the flow has become numerically unstable. Watch the Courant number, diffusion number, and divergence printed periodically by `modchecksim` (controlled by `tcheck` in `NAMCHECKSIM`); if the Courant number is climbing, reduce `dtmax` or enable `ladaptive` so the solver keeps `courant` within its recommended `1 <= courant <= 2` range — see the [`dtmax`/`ladaptive`/`courant` entries](udales-namoptions-overview.md#namelist-run) and the [`tcheck` entry](udales-namoptions-overview.md#namelist-namchecksim).

**`Error: y-averaged statistics not currently implemented for nprocx > 1.`** / **`Error: constant x outflow only possible for nprocx = 1.`** / **`Error: constant y outflow only possible for nprocy = 1.`** — `lydump`/`lytdump`, `luoutflowr`, and `lvoutflowr` each have a decomposition restriction: run with a single pencil in the relevant direction (`nprocx = 1` or `nprocy = 1`), or disable the option. See `src/modstartup.f90`.
