# Running uDALES

Simulations are run with `u-dales/bin/ud_run`, a single command that works on every supported machine:

``` sh
# General syntax: ud_run <machine> sim|gather <case_directory>
./u-dales/bin/ud_run ecmwf sim experiments/009      # run (or submit) the simulation
./u-dales/bin/ud_run ecmwf gather experiments/009   # merge the per-CPU outputs afterwards
```

Running `ud_run` with no arguments prints its usage. `<machine>` uses the same names as `tools/build_executable.sh` for the machines listed below (that script also has a legacy `cca` target, which `ud_run` does not support), and `<case_directory>` is always the experiment directory containing `config.sh`. `ud_run` only picks the machine-specific script in `u-dales/tools` and runs it:

| `<machine>` | System | `sim` runs | `gather` runs |
|---|---|---|---|
| `common` | Local workstation | `local_execute.sh` (also gathers at the end) | `gather_outputs.sh` on `$DA_WORKDIR/<exp>` |
| `icl` | Imperial CX3 (PBS) | `hpc_execute.sh` | `hpc_gather.sh` |
| `hx1` | Imperial HX1 (PBS) | `hpc_execute.sh` | `hpc_gather.sh` |
| `archer` | ARCHER2 (Slurm) | `archer_execute.sh` | `archer_gather.sh` on `$DA_WORKDIR/<exp>` |
| `ecmwf` | ECMWF Atos HPC2020 (Slurm) | `ecmwf_execute.sh` | `ecmwf_gather.sh` |

On clusters, `sim` and `gather` submit batch jobs; submit `gather` once the simulation has finished (see [Post-processing](./udales-post-processing.md)). The machine-specific scripts can still be called directly, e.g. `./u-dales/tools/ecmwf_execute.sh experiments/009`. `ud_run` lives in `u-dales/bin`, so if you add that directory to your `PATH` (see below) you can call it as plain `ud_run` from anywhere.

The scripts require several variables to be set up. Below is an example setup for copying and pasting. You can also specify these parameters in a `config.sh` file within the example directory, which is then read by the scripts. We recommend keeping a `config.sh` in each example case directory with the appropriate variable setting.
The simulation workflow consists of three stages:

1. Pre-processing: create or update input files with `write_inputs.sh` (calls `write_inputs.m` or `write_inputs.py`) or other pre-processing tools.
2. Execution: launch the solver using `ud_run <machine> sim`.
3. Post-processing: merge output files using `ud_run <machine> gather`, then analyse them.

Example cases shipped with uDALES are located under `u-dales/examples/` and are suitable for testing an installation.

Note that you need to choose the number of CPUs you are using to run the simulation such that the product of `nprocx` and `nprocy` (in the `namoptions` input file) is equal to the total number of CPU asked, i.e., `nprocx * nprocy = NCPU` for local machines, and `nprocx * nprocy = NNODE * NCPU` for ICL HPC, ARCHER2 or ECMWF clusters.

## Run on common systems

``` sh
# Contents of the config.sh file
# We assume you are running the following commands from your
# top-level project directory.

export DA_EXPDIR=$(pwd)/experiments                     # Experiments top-level directory
export DA_TOOLSDIR=$(pwd)/u-dales/tools                 # Directory of scripts
export DA_BUILD=$(pwd)/u-dales/build/release/u-dales    # Build file
export DA_WORKDIR=$(pwd)/outputs                        # Output top-level directory
export NCPU=8                                           # Number of CPUs to use for a simulation

# It is recommended to write the full path instead of using $(pwd) in config.sh file
```

`tools/write_inputs.sh` scans `namoptions.###` for `nompthreads` and uses it for the preprocessing CPU request. If `nompthreads` is omitted, the preprocessing default is `8`; if it appears more than once anywhere in the file, the wrapper stops and asks for a single value. The wrapper exports the derived value internally as `PREPROC_NCPU` for PBS submission and for the default View3D thread count, unless `VIEW3D_NUM_THREADS` is explicitly overridden. For preprocessing, `DA_TOOLSDIR` defaults to the wrapper script directory and `DA_EXPDIR` defaults to the parent directory of the case directory unless these are set in `config.sh` or the calling environment. On Imperial HPC, compute-node preprocessing with `write_inputs.sh ... c` defaults to `PREPROC_WALLTIME="24:00:00"` and `PREPROC_MEM="128gb"`. `PREPROC_MEM` must be written as a number followed by lowercase `gb`, such as `128gb`; unitless values are rejected. Unless `VIEW3D_MAX_DENSE_MATRIX_GIB` is explicitly set, the default View3D dense-matrix guard leaves 16 GiB for overhead when `PREPROC_MEM` is above `16gb`, while smaller requests use the requested GiB value.

Then, to start the simulation, run:

``` sh
# We assume you are running the following commands from your
# top-level project directory.

# Runs local_execute.sh, which also gathers the outputs at the end
./u-dales/bin/ud_run common sim experiments/009
```

## Run on ICL clusters (CX3 and HX1)

``` sh
# Contents of the config.sh file
# We assume you are running the following commands from your
# top-level project directory.

export DA_EXPDIR=$(pwd)/experiments                     # Experiments top-level directory
export DA_TOOLSDIR=$(pwd)/u-dales/tools                 # Directory of scripts
export DA_BUILD=$(pwd)/u-dales/build/release/u-dales    # Build file
export DA_WORKDIR=$EPHEMERAL                            # Output top-level directory
export NCPU=128                                         # Number of CPUs to use for a simulation
export PREPROC_WALLTIME="24:00:00"                      # Optional preprocessing override; defaults to 24:00:00
export PREPROC_MEM="128gb"                              # Optional preprocessing override; defaults to 128gb
export NNODE=1                                          # Number of nodes to use for a simulation
export WALLTIME="00:30:00"                              # Maximum runtime for simulation in hours:minutes:seconds
export MEM="128gb"                                      # Memory request per node

# It is recommended to write the full path instead of using $(pwd) in config.sh file
```

For guidance on how to set the parameters on HPC, have a look at [Job sizing guidance](https://icl-rcs-user-guide.readthedocs.io/en/latest/hpc/queues/job-sizing-guidance/).
Then, to start the simulation, run:

``` sh
# We assume you are running the following commands from your
# top-level project directory.

# Runs hpc_execute.sh; use hx1 instead of icl on HX1
./u-dales/bin/ud_run icl sim experiments/009
```

## Run on ARCHER2

``` sh
# Contents of the config.sh file
# We assume you are running the following commands from your
# top-level project directory.

export DA_EXPDIR=/work/account/account/username/top_level_project_directory/experiments                     # Experiments top-level directory
export DA_TOOLSDIR=/work/account/account/username/top_level_project_directory/u-dales/tools                 # Directory of scripts
export DA_BUILD=/work/account/account/username/top_level_project_directory/u-dales/build/release/u-dales    # Build file
export DA_WORKDIR=/work/account/account/username/top_level_project_directory/outputs                        # Output top-level directory
export NCPU=128                                                                                             # Number of CPUs to use for a simulation
export NNODE=1                                                                                              # Number of nodes to use for a simulation
export WALLTIME="24:00:00"                                                                  # Maximum runtime for simulation in hours:minutes:seconds
export MEM="256gb"                                                                          # Memory request per node
export QOS="standard"                                                                       # Queue
```

For guidance on how to set the parameters on ARCHER2, have a look at the [ARCHER2 documentation](https://docs.archer2.ac.uk/user-guide/). In particular, make sure to edit the `archer_execute.sh` script (the line `#SBATCH --account=n02-ASSURE`) and set the account corresponds to one you use.
Then, to start the simulation, run:

``` sh
# We assume you are running the following commands from your
# top-level project directory.

# Runs archer_execute.sh
./u-dales/bin/ud_run archer sim experiments/009
```

## Run on ECMWF HPC2020

The ECMWF Atos HPC2020 uses Slurm. Build first with `tools/build_executable.sh ecmwf release` (see [Installation](./udales-installation.md#build-on-hpcs)).

``` sh
# Contents of the config.sh file
# We assume you are running the following commands from your
# top-level project directory.

export DA_EXPDIR=/home/username/top_level_project_directory/experiments           # Experiments top-level directory
export DA_TOOLSDIR=/home/username/top_level_project_directory/u-dales/tools       # Directory of scripts
export DA_BUILD=/home/username/top_level_project_directory/u-dales/bin/u-dales    # Executable
export DA_WORKDIR=$SCRATCH/outputs                                                # Output top-level directory; must be under $SCRATCH
export NCPU=128                                                                   # MPI ranks per node (at most 128)
export NNODE=2                                                                    # Number of nodes to use for a simulation
export WALLTIME="24:00:00"                                                        # Maximum runtime in hours:minutes:seconds (np and nf allow up to 48:00:00)
export QOS="np"                                                                   # Queue for the simulation; required (the gather always uses nf)
export MEM="16G"                                                                  # Memory, 16G or 16gb; required for the gather and for simulations not on np
```

Things to know on ECMWF:

- **Where outputs go.** `DA_WORKDIR` must be under `$SCRATCH`, otherwise `ecmwf_execute.sh` refuses to submit. `$SCRATCH` is the only filesystem sized for multi-TB output (50 TB quota); `$HOME`, `$PERM` and `$HPCPERM` are small and not meant for parallel I/O. Files on `$SCRATCH` are **deleted automatically 30 days after last access**, so move results you want to keep (e.g. to ECFS). The inputs are copied to `$DA_WORKDIR/<exp>` at submission, so editing the case directory while the job is queued does not affect it. The executable is not copied: the job runs it from `DA_BUILD`, so rebuilding while a job is queued *does* change what runs.
- **No spaces in paths.** ECMWF's `sbatch` wrapper splits its arguments on whitespace, so neither the case directory nor `DA_WORKDIR` may contain a space; the scripts refuse before submitting anything.
- **Nodes.** Each compute node has 128 physical cores and 240 GB of memory. The job runs `NCPU` MPI ranks per node on physical cores only (no hyperthreading), so `NCPU` must be at most 128, and `nprocx * nprocy` must equal `NNODE * NCPU`. For example, about 10,000 ranks needs `NNODE=79` with `NCPU=128`.
- **Queue (`QOS`).** `QOS` sets the queue for the **simulation** and must be present in `config.sh`; there is no default, and `ecmwf_execute.sh` stops with a message listing the sensible values if it is missing. Use `np` for real runs: nodes are exclusive and each comes with all 240 GB, so no memory request is needed. Use `nf` for small test runs: nodes are shared, limited to 1 node with at most 128 cores and 128 GB. Other queues are listed in the comments of `ecmwf_execute.sh` and in the [ECMWF batch system guide](https://confluence.ecmwf.int/display/UDOC/HPC2020%3A+Batch+system).
- **The gather job always runs on `nf`.** `ecmwf_gather.sh` fixes its own queue to `nf` and ignores `QOS`, because gathering is a single NCO process and belongs on a shared node. This is deliberate: on `np` that one process would be given a whole exclusive node, and all 128 cores would be billed for the duration (SBUs are charged on allocated cores × elapsed time).
- **Memory (`MEM`).** The gather job always requests `MEM`. The simulation requests it only when `QOS` is not `np`, because every other queue shares nodes, where a job that asks for nothing is given the shared pool default of **8 GB** — enough to get the solver killed mid-run with `oom-kill`. Both the Slurm spelling (`16G`) and the PBS one used on Imperial HPC (`16gb`) are accepted, in any case, and converted to what Slurm expects; the unit itself is required, so a bare `16` is rejected as ambiguous.
- **Account.** Jobs are charged to the Slurm account set by `ACCOUNT=` near the top of `ecmwf_execute.sh` and `ecmwf_gather.sh`. Change it to your own project; `account` lists the accounts you can use, and budgets are shown at [hpc-usage.ecmwf.int](https://hpc-usage.ecmwf.int).
- **Logs.** The job scripts (`job.<exp>.slurm`, `post-job.<exp>.slurm`), Slurm logs (`slurm-<jobid>.out`) and `output.<exp>.log` are written to `$DA_WORKDIR/<exp>`.

See the [ECMWF HPC2020 user guide](https://confluence.ecmwf.int/spaces/UDOC/pages/240851027/HPC2020+User+Guide) for more on queues, filesystems and accounting.
Then, to start the simulation, run:

``` sh
# We assume you are running the following commands from your
# top-level project directory.

# Runs ecmwf_execute.sh
./u-dales/bin/ud_run ecmwf sim experiments/009
```
