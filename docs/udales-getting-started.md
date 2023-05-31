# Getting Started

Getting started with uDALES to set up your own experiments is straightforward. This guide goes through the steps required to [install](#installation) uDALES, and [set-up](#set-up) and [run](#run) a simple example. Results are outputted in netCDF format, for a quick inspection you can use GUI tools such as [Panoply](https://www.giss.nasa.gov/tools/panoply/) or [ncview](http://meteora.ucsd.edu/~pierce/ncview_home_page.html). To learn more about pre- and post-processing steps see the [what's next section](#whats-next).


## Singularity

If you have [Singularity](https://sylabs.io/) available on your system, you can use the provided scripts under `tools/singularity` to build and run uDALES cases locally or on HPC environments, for other options, see the sections below. If you are looking for information on how to install or use Singularity on your system, please refer to the [Singularity documentation ](https://sylabs.io/docs). The use of Singularity is undoubtedly the easiest way to build and run cases in uDALES as all dependencies are provided and uDALES will compile out of the box. Furthermore, users wishing to achieve a reasonable level of scientific reproducibility may archive software, tools, and data with their Singularity image containing OS and external libraries to an open access repository (e.g. [Meyer et al., 2020](https://doi.org/10.1029/2019MS001961)).

First clone the uDALES repository with:

```sh
https://github.com/uDALES/u-dales.git
```

Then, to build and download the Singularity image use:

```sh
singularity build --remote tools/singularity/image.sif tools/singularity/image.def
```

then, to install uDALES use:

```sh
# udales_build.sh <NPROC> [Debug, Release]
./tools/singularity/udales_build.sh 2 Release
```

Finally, to run an example case use:

```sh
# udales_run.sh <NPROC> <BUILD_TYPE> <PATH_TO_CASE> <NAMELIST>
./tools/singularity/udales_run.sh 2 Release examples/001 namoptions.001
```

If you are looking to run the build and run commands on HPC, we have provided a sample script under `tools/singularity/udales_pbs_submit.sh`, you can modify and run it with `qsub tools/singularity/udales_pbs_submit.sh`.


## Prerequisites when not using Singularity

### uDALES

uDALES is supported to run on Linux, macOS and Windows Subsystem for Linux (WSL). Please ensure that the latest version of the following libraries and software are available on your system:

- [CMake](https://cmake.org/) >= 3.9.
- [NetCDF-Fortran](https://www.unidata.ucar.edu/downloads/netcdf/index.jsp) >= 4.
- [GNU](https://gcc.gnu.org/wiki/GFortran) <= 9, [Intel](https://software.intel.com/content/www/us/en/develop/documentation/fortran-compiler-developer-guide-and-reference/top.html), or [Cray](https://pubs.cray.com/) Fortran compiler.
- A recent version of [MPICH](https://www.mpich.org/) or [Open-MPI](https://www.open-mpi.org/).
- [FFTW](http://www.fftw.org/) 

### Project setup

This guide helps you set up a project template for uDALES with a generic folder structure set-up that you can later use to set up your own experiments. For this you also need:

- [Git](https://git-scm.com/) >= 2.
- A [GitHub](https://github.com) account. (optional)
- [Python](https://www.python.org/) >= 3.6.

### Pre-processing

When you create your own experiments, you will need to set up specific input files. We have a system in place that does that for you, written in MATLAB. Information can be found under [pre-processing](./udales-pre-processing.md) and is not discussed in the getting-started set-up.

- [MATLAB](https://www.mathworks.com/products/matlab.html)

### Post-processing

For better organised netcdf output files, you will need:

- [netCDF Operators](https://github.com/nco/nco) (NCO).

On local systems, these software and libraries (except MATLAB) should be available from your system's package manager (e.g. APT, yum, Homebrew, etc.) and examples on how to install all the required libraries for Linux/WSL and macOS are given below.

On high performance computing (HPC) clusters, these software and libraries should have already been installed. Please refer to the specific documentation to load the above software and libraries. Alternatively, you can install all the required packages easily after installing [Linuxbrew](https://docs.brew.sh/Homebrew-on-Linux) and using the instructions for [macOS](#macos).

### Linux/WSL (Ubuntu)

```sh
sudo apt-get update && sudo apt-get upgrade -y
sudo apt-get install -y git cmake gfortran libopenmpi-dev openmpi-bin libnetcdf-dev libnetcdff-dev nco python3 python3-pip libfftw3-dev
```

### macOS

On macOS, use [Homebrew](https://docs.brew.sh) to install the required libraries. If you do not have Homebrew installed on your system, install it from the [Homebrew installation page](https://docs.brew.sh/Installation) then, to install all the required dependencies, including support for MPI, run the following commands from your terminal prompt:

```sh
brew update
brew install git cmake gcc@9 netcdf open-mpi nco python3 fftw
```


## Installation

The installation and set-up of uDALES is straightforward thanks to the use of a [Cookiecutter repository](https://github.com/uDALES/cookiecutter-u-dales) to create a project template for uDALES with a generic folder structure set-up that you can later use to set up your own experiments.

Please make sure that you have created a new (private) repository on your GitHub account and made note of the URL as we need this later to configure Cookiecutter.

### Repository set-up

First, install [Cookiecutter](https://github.com/cookiecutter/cookiecutter) from your command prompt:

``` sh
python3 -m pip install --user cookiecutter
```

Then, to create a new uDALES project within the current working directory:

``` sh
cookiecutter https://github.com/uDALES/cookiecutter-u-dales
```

and fill in the required fields when prompted. `<PROJECT_NAME>` is the name of the generated project directory and `<GITHUB_PROJECT_REMOTE>` is the URL to your remote GitHub account (this is optional, you can just press the return key to leave this empty). E.g.:

``` sh
directory_name [<PROJECT_NAME>]: neutral_experiments
github_project_remote [<GITHUB_PROJECT_REMOTE>]: https://github.com/<MY_GITHUB_USERNAME>/<MY_NEW_EMPTY_REPO>.git
```

This creates a Git repository for your own projects named `<PROJECT_NAME>` with the [uDALES model development repository](https://github.com/uDALES/u-dales) as submodule, and a generic tree that you can use to set up your own experiments:

``` sh
.
├── data        # Contains or links to any external data used by the experiments.
├── docs        # Relevant documentation or papers used for the experiment.
├── experiments # Configuration files grouped by experiment number.
│   └── <N>     # Any configurations files needed by uDALES to run experiment <N>.
│   └── ...
├── tools       # Additional or specialized tools other then the ones included with uDALES.
└── u-dales     # uDALES model development repository (submodule).
```

In the next steps we will assume your current working directory is the top-level project directory you just created with Cookiecutter.

### Build on common systems

On standard systems and configurations, you can build uDALES with the following commands:

```sh
# We assume you are running the following commands from your
# top-level project directory.

mkdir -p u-dales/build/release # in case you want to later create a build/debug
pushd u-dales/build/release
cmake -LA ../..
make
popd
```

You can compile in parallel mode by passing Make the `j` flag followed by the number of CPU cores to use. For example, to compile with 2 cores do `make -j2`.

### Build on HPCs

If you are a High Performance Cluster (HPC) user you are likely using the [Environment Modules package](http://modules.sourceforge.net/) for the dynamic modification of the user's environment via modulefiles and therefore you may need to hint CMake the PATH to NetCDF (see below how).

Here we show how to compile uDALES using the [HPC at ICL](https://www.imperial.ac.uk/admin-services/ict/self-service/research-support/rcs/) as an example, therefore please note that the specific names/versions installed on your system may be different.

``` sh
module list # list currently enabled modules -- should be empty!
module avail # list available modules
```

``` sh
# This is an example, please check with the previous command for the exact name of the
# modules available on your system. This will load NetCDF compiled with Intel Suite
# 2019.4 and add the correct version of icc and ifort to the PATH.
module load intel-suite/2017.6 mpi/intel-2018 cmake/3.14.0 git/2.14.3
```

Then, to build the uDALES executable, run the following commands:

``` sh
# We assume you are running the following commands from your
# top-level project directory.

mkdir -p u-dales/build/release
pushd u-dales/build/release
FC=mpiifort cmake -DNETCDF_DIR=/apps/netcdf/4.4.1-c -DNETCDF_FORTRAN_DIR=/apps/netcdf/4.4.4-fortran -LA ../..
make
popd
```

where `NETCDF_DIR` and `NETCDF_FORTRAN_DIR` indicates the absolute path to your NetCDF-C and NetCDF-Fortran installation directories. Here, we use the utilities `nc-config` and `nf-config` to hint CMake the location of NetCDF, but you can simply pass the absolute path to the NetCDF-C and NetCDF-Fortran manually instead. You can compile in parallel mode by passing Make the `j` flag followed by the number of CPU cores to use. For example, to compile with 2 cores do `make -j2`.

### Build defaults/options

By default uDALES will compile in `Release` mode. You can change this by specifying the option (or flag) at configure time. The general syntax for specifying an option in CMake is `-D<flag_name>=<flag_value>` where `<flag_name>` is the option/flag name and `<flag_value>` is the option/flag value. The following options can be specified when configuring uDALES:

| Name                            | Options            | Default   | Description                                   |
| ------------------------------- | ------------------ | --------- | --------------------------------------------- |
| `CMAKE_BUILD_TYPE`              | `Release`, `Debug` | `Release` | Whether to optimise/build with debug flags    |
| `NETCDF4_DIR`                   | `<path>`           | -         | Path to NetCDF-C installation directory       |
| `NETCDF_FORTRAN_DIR`            | `<path>`           | -         | Path to NetCDF-Fortran installation directory |
| `SKIP_UPDATE_EXTERNAL_PROJECTS` | `ON`, `OFF`        | `OFF`     | Whether to skip updating external projects    |

## Set-up

To set up a new simulation, `copy_inputs.sh` in `u-dales/tools/` is used to create a new simulation setup `new_exp_id` based on another simulation `old_exp_id`. All `exp_ids` are three digit numbers, e.g. 001, and are stored in directories of that name. Scripts requires several variables to be set up. You can do this by copying and pasting the snippet below or by including it in a bash script (or bash profile if you are unlikely to change them).

``` sh
# We assume you are running the following commands from your
# top-level project directory.

export DA_EXPDIR=$(pwd)/experiments #  The top-level directory of the simulation setups.
export DA_WORKDIR=$(pwd)/outputs # Output top-level directory

# If source directories (DA_EXPDIR_SRC, DA_WORKDIR_SRC) are not set,
# the experiment set-up folder will be copied from the same target directory.
# I.e. DA_EXPDIR_SRC==DA_EXPDIR and DA_WORKDIR_SRC==DA_WORKDIR.
export DA_EXPDIR_SRC=$(pwd)/u-dales/examples
export DA_WORKDIR_SRC=$(pwd)/u-dales/examples
```

If you set up a new experiment on HPC, also use:

``` sh
export DA_WORKDIR=$EPHEMERAL # Output top-level directory on HPC
export DA_WORKDIR_SRC=$EPHEMERAL
```

Now to set-up a new experiment (here we use case `009`) based on a previous example (here we use case `001`), run:

``` sh
# We assume you are running the following commands from your
# top-level project directory.

# General syntax: copy_inputs.sh new_exp_id old_exp_id
# To set up a new simulation starting from the restart files of another simulation
# ("warmstart"), use the 'w' flag. E.g.: copy_inputs.sh new_exp_id old_exp_id w
./u-dales/tools/copy_inputs.sh 009 001
```

## Run

The scripts `local_execute.sh` and `hpc_execute.sh` in `u-dales/tools` are used as wrappers to run simulations on your local machine and HPC at ICL respectively. These scripts contain several helpers to run the simulations and merge outputs from several CPUs into a single file (see [Post-processing](./udales-post-processing.md) for more info about the individual scripts).

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

### Run on HPCs

``` sh
export DA_TOOLSDIR=$(pwd)/u-dales/tools # Directory of scripts
export DA_BUILD=$(pwd)/u-dales/build/release/u-dales # Build file
export NCPU=2 # Number of CPUs to use for a simulation

export NNODE=1 # Number of nodes to use for a simulation
export WALLTIME="00:30:00" # Maximum runtime for simulation in hours:minutes:seconds
export MEM="2gb" # Memory request per node
```

For guidance on how to set the parameters on HPC, have a look at [Job sizing guidance](https://www.imperial.ac.uk/admin-services/ict/self-service/research-support/rcs/computing/job-sizing-guidance/).
Then, to start the simulation, run:

``` sh
# We assume you are running the following commands from your
# top-level project directory.

# General syntax: hpc_execute.sh exp_directory
./u-dales/tools/hpc_execute.sh experiments/009
```

## What's next?

This simple guide is meant to allow you to get started based on an existing example. We set up several [example simulations](./udales-example-simulations.md) to get you started. To learn more about pre- and post-processing steps in uDALES, please see the [pre-processing](./udales-pre-processing.md) and [post-processing](./udales-post-processing.md) pages.
