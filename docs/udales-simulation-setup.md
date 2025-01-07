## Running uDALES

The scripts `local_execute.sh` (for local machines), `hpc_execute.sh` (for ICL cluster) and `archer_execute.sh` (for ARCHER2) in `u-dales/tools` are used as wrappers to run simulations. These scripts contain several helpers to run the simulations and merge outputs from several CPUs into a single file (see [Post-processing](./udales-post-processing.md) for more info about the individual scripts).

The scripts require several variables to be set up. Below is an example setup for copying and pasting. You can also specify these parameters in a `config.sh` file within the experiment directory, which is then read by the scripts.

Note that you need to choose the number of CPUs you are using to run the simulation such that the number of grid cells in the y-direction (`jtot` parameter in the `namoptions` input file) is a multiple of the number of CPUs.

### Run on common systems

``` sh
# We assume you are running the following commands from your
# top-level project directory.

export DA_TOOLSDIR=$(pwd)/u-dales/tools # Directory of scripts
export DA_BUILD=$(pwd)/u-dales/build/release/u-dales # Build file
export NCPU=2 # Number of CPUs to use for a simulation
export DA_WORKDIR=$(pwd)/outputs # Output top-level directory
```

Then, to start the simulation, run:

``` sh
# We assume you are running the following commands from your
# top-level project directory.

# General syntax: local_execute.sh exp_directory
./u-dales/tools/local_execute.sh experiments/009
```

### Run on ICL cluster

``` sh
export DA_TOOLSDIR=$(pwd)/u-dales/tools # Directory of scripts
export DA_BUILD=$(pwd)/u-dales/build/release/u-dales # Build file
export NCPU=24 # Number of CPUs to use for a simulation
export NNODE=1 # Number of nodes to use for a simulation
export WALLTIME="00:30:00" # Maximum runtime for simulation in hours:minutes:seconds
export MEM="128gb" # Memory request per node
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
export DA_TOOLSDIR=$(pwd)/u-dales/tools # Directory of scripts
export DA_BUILD=$(pwd)/u-dales/build/release/u-dales # Build file
export NCPU=128 # Number of CPUs to use for a simulation
export NNODE=1 # Number of nodes to use for a simulation
export WALLTIME="24:00:00" # Maximum runtime for simulation in hours:minutes:seconds
export MEM="256gb" # Memory request per node
export QOS="standard" # Queue
```

For guidance on how to set the parameters on ARCHER2, have a look at the [ARCHER2 documentation](https://docs.archer2.ac.uk/user-guide/). In particular, take care to edit the `archer_execute.sh` script so that the account corresponds to one you can use.
Then, to start the simulation, run:

``` sh
# We assume you are running the following commands from your
# top-level project directory.

# General syntax: hpc_execute.sh exp_directory
bash ./u-dales/tools/archer_execute.sh experiments/009
```
