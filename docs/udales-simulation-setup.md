## Running uDALES

The scripts `local_execute.sh` (for local machines), `hpc_execute.sh` (for ICL cluster) and `archer_execute.sh` (for ARCHER2) in `u-dales/tools` are used as wrappers to run simulations. These scripts contain several helpers to run the simulations and merge outputs (see [Post-processing](./udales-post-processing.md) for more info) from several CPUs into a single file.

The scripts require several variables to be set up. Below is an example setup for copying and pasting. You can also specify these parameters in a `config.sh` file within the example directory, which is then read by the scripts. We recommend keeping a `config.sh` in each example case directory with the appropriate variable setting.
The simulation workflow consists of three stages:

1. Pre-processing: create or update input files with `write_input.sh` (calls `write_input.m`) or other pre-processing tools.
2. Execution: launch the solver using `local_execute.sh` (for desktop) or `hpc_execute.sh` (for clusters).
3. Post-processing: merge and analyse output files using `gather_outputs.sh` or other scripts.

Example cases shipped with uDALES are located under `u-dales/examples/` and are suitable for testing an installation.

Note that you need to choose the number of CPUs you are using to run the simulation such that the product of `nprocx` and `nprocy` (in the `namoptions` input file) is equal to the total number of CPU asked, i.e., `nprocx * nprocy = NCPU` for local machines, and `nprocx * nprocy = NNODE * NCPU` for ICL HPC or ARCHER2 clusters.

### Run on common systems

``` sh
# Contents of the config.sh file
# We assume you are running the following commands from your
# top-level project directory.

export DA_EXPDIR=$(pwd)/experiments                     # Experiments top-level directory
export DA_TOOLSDIR=$(pwd)/u-dales/tools                 # Directory of scripts
export DA_BUILD=$(pwd)/u-dales/build/release/u-dales    # Build file
export DA_WORKDIR=$(pwd)/outputs                        # Output top-level directory
export NCPU=8                                           # Number of CPUs to use for a simulation
export PREPROC_NCPU=8                                   # Number of CPUs to use for preprocessing; match &INPS nompthreads

# It is recommended to write the full path instead of using $(pwd) in config.sh file
```

`PREPROC_NCPU` is required when running `tools/write_inputs.sh`. It should be the same integer as `nompthreads` under the `&INPS` section of `namoptions.###`. The default value of `nompthreads` is `8`, so use `export PREPROC_NCPU=8` unless you explicitly set a different `nompthreads` value in `namoptions.###`. On Imperial HPC, compute-node preprocessing with `write_inputs.sh ... c` also requires `PREPROC_WALLTIME` and `PREPROC_MEM` in `config.sh`. The default View3D dense-matrix guard is derived as `PREPROC_MEM - 16 GiB`, unless `VIEW3D_MAX_DENSE_MATRIX_GIB` is explicitly set.

Then, to start the simulation, run:

``` sh
# We assume you are running the following commands from your
# top-level project directory.

# General syntax: local_execute.sh exp_directory
./u-dales/tools/local_execute.sh experiments/009
```

### Run on ICL cluster

``` sh
# Contents of the config.sh file
# We assume you are running the following commands from your
# top-level project directory.

export DA_EXPDIR=$(pwd)/experiments                     # Experiments top-level directory
export DA_TOOLSDIR=$(pwd)/u-dales/tools                 # Directory of scripts
export DA_BUILD=$(pwd)/u-dales/build/release/u-dales    # Build file
export DA_WORKDIR=$EPHEMERAL                            # Output top-level directory
export NCPU=128                                         # Number of CPUs to use for a simulation
export PREPROC_NCPU=8                                   # Number of CPUs to use for preprocessing; match &INPS nompthreads
export PREPROC_WALLTIME="48:00:00"                      # Walltime for Imperial HPC preprocessing jobs
export PREPROC_MEM="128gb"                              # Memory request for Imperial HPC preprocessing jobs
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

# General syntax: hpc_execute.sh exp_directory
./u-dales/tools/hpc_execute.sh experiments/009
```

### Run on ARCHER2

``` sh
# Contents of the config.sh file
# We assume you are running the following commands from your
# top-level project directory.

export DA_EXPDIR=/work/account/account/username/top_level_project_directory/experiments                     # Experiments top-level directory
export DA_TOOLSDIR=/work/account/account/username/top_level_project_directory/u-dales/tools                 # Directory of scripts
export DA_BUILD=/work/account/account/username/top_level_project_directory/u-dales/build/release/u-dales    # Build file
export DA_WORKDIR=/work/account/account/username/top_level_project_directory/outputs                        # Output top-level directory
export NCPU=128                                                                                             # Number of CPUs to use for a simulation
export PREPROC_NCPU=8                                                                                       # Number of CPUs to use for preprocessing; match &INPS nompthreads
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

# General syntax: hpc_execute.sh exp_directory
bash ./u-dales/tools/archer_execute.sh experiments/009
```
