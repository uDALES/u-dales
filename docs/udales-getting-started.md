# Getting Started

Getting started with uDALES to set up your own experiments is straightforward. This guide goes through the steps required to [install](#installation) uDALES, and [set-up](#set-up) and [run](#run) a simple exmaple.

## Prerequisites

uDALES is supported to run on Linux, macOS and Windows Subsystem for Linux (WSL). Please ensure that the latest version of the following libraries and software are avalable on your system:

- [Git](https://git-scm.com/).
- [CMake](https://cmake.org/).
- [NetCDF-Fortran](https://www.unidata.ucar.edu/downloads/netcdf/index.jsp).
- A Fortran compiler (e.g. [GNU Fortran compiler](https://gcc.gnu.org/wiki/GFortran)).
- An MPI library implementation (e.g. [MPICH](https://www.mpich.org/)).
- [netCDF Operators](https://github.com/nco/nco) (NCO).
- [Python](https://www.python.org/) 3.5 or above.
- A [GitHub](https://github.com) account.

On local systems, these software and libraries should be available from your system's package manager (e.g. APT, yum, Homebrew, etc.) and examples on how to install all the required libraries for Linux/WSL and macOS are given below.

On high performance computing (HPC) clusters, these software and libraries should have already been installed. Please refer to the specific documentation to load the above software and libraries. Alternatively, you can install all the required packages easily after installing [Linuxbrew](https://docs.brew.sh/Homebrew-on-Linux) and using the instructions for [macOS](#macos).


### Linux/WSL (Ubuntu)

```sh
sudo apt-get update && sudo apt-get upgrade -y
sudo apt-get install -y git cmake gfortran libopenmpi-dev openmpi-bin libnetcdf-dev libnetcdff-dev nco python3 python3-pip
```


### macOS

On macOS, use [Homebrew](https://docs.brew.sh) to to install the required libraries. If you do not have Homebrew installed on your system, install it from the [Homebrew installation page](https://docs.brew.sh/Installation) then, to install all the required dependencies, including support for MPI, run the following commands from your terminal prompt:

```sh
brew update
brew install git cmake gcc netcdf open-mpi nco python3
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

and fill in the required fields when prompted. `<PROJECT_NAME>` is the name of the generated project directory and `<GITHUB_PROJECT_REMOTE>` is the URL to your remote GitHub account. E.g.:

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
├── tools       # Additional or specialized tools other then the ones already included with uDALES.
└── u-dales     # uDALES model development repository (submodule).
```

### Build on common systems

On standard systems and configurations, you can build uDALES with the following commands:

```sh
# We assume you are running the following commands from the uDALES project directory.
# I.e. the one you set up earlier with Cookiecutter.
mkdir u-dales/build
pushd u-dales/build
cmake ..
make -j$(nproc)
popd
```

Where `$(nproc)` will use all the number of CPU cores/threads available on your system. Note that using the maximum number of CPU cores/threads available may not necessarily be the fastest way to build the software, therefore you may want to manually specify the number of CPU cores/threads to use manually.


### Build on HPC

If you are an HPC user, you are, most likely, using the [Environment Modules package](http://modules.sourceforge.net/) for the dynamic modification of the user's environment via modulefiles. In that case, you will need to specify the path to the NetCDF manually _after_ loading all the libraries required to compile uDALES. For example:

``` sh
# This is an example, module names/versions may be different on your system
module list # list currently enabled modules
module avail # list available modules
module load cmake netcdf4 openmpi gnu # This is an example, please check with the previous command for the exact name of the modules available on your system.
```

If you do not know the location of NetCDF, you can locate it with the `nc-config --prefix`. Then, to build the uDALES executable, from the `u-dales` repository, run the following commands:

``` sh
# We assume you are running the following commands from the uDALES project directory.
# I.e. the one you set up earlier with Cookiecutter.
mkdir u-dales/build
pushd u-dales/build
cmake -DNETCDF_DIR=$NETCDF4_DIR -DNETCDF_FORTRAN_DIR=$DNETCDF_FORTRAN_DIR ..
make -j$(nproc)
popd
```

where `$NETCDF4_DIR` and `$DNETCDF_FORTRAN_DIR` are the absolute path to your NetCDF-C and NetCDF-Fortran installation directories.

#### Note for ICL users
If you are an HPC user at ICL, please note that while there is an issue with how the libraries are loaded and added to the `PATH` on This issue has been reported and should be resolved 'soon'. As a workaround, you can use the following instead:


``` sh
# The following assumes that you have no other modules loaded on your environment.
module load intel-suite/2017.6 mpi/intel-2018 cmake/3.14.0 git/2.14.3

# We assume you are running the following commands from the uDALES project directory.
# I.e. the one you set up earlier with Cookiecutter.
mkdir u-dales/build
pushd u-dales/build
cmake ..
FC=mpiifort cmake -DNETCDF_DIR=/apps/netcdf/4.4.1-c -DNETCDF_FORTRAN_DIR=/apps/netcdf/4.4.4-fortran ..
make -j$(nproc)
popd
```

### Build defaults/options

By default uDALES will compile in `Release` mode. You can change this by specifying the option (or flag) at configure time. The general syntax for specifying an option in CMake is `-D<flag_name>=<flag_value>` where `<flag_name>` is the option/flag name and `<flag_value>` is the option/flag value. The following options can be specified when configuring uDALES:

| Name                 | Options            | Default   | Description                                   |
| -------------------- | ------------------ | --------- | --------------------------------------------- |
| `CMAKE_BUILD_TYPE`   | `Release`, `Debug` | `Release` | Whether to optimise/build with debug flags    |
| `NETCDF4_DIR`        | `<path>`           | -         | Path to NetCDF-C installation directory       |
| `NETCDF_FORTRAN_DIR` | `<path>`           | -         | Path to NetCDF-Fortran installation directory |


## Set-up

To set up a new simulation, `da_prep.sh` in `u-dales/tools/utils` is used to create a new simulation setup `new_exp_id` based on another simulation `old_exp_id`. All `exp_ids` are three digit numbers, e.g. 001, and are stored in directories of that name. Scripts requires several variables to be set up. You can do this by copying and pasting the snippet below or by including it in a bash script (or bash profile if you are unlikely to change them).

``` sh
# We assume you are running the following commands from the uDALES project directory.
# I.e. the one you set up earlier with Cookiecutter.
export DA_TOPDIR=$(pwd) # This is your top level directory (i.e. Cookiecutter).
export DA_EXPDIR=$(pwd)/experiments #  The top-level directory of the simulation setups.
export DA_WORKDIR=$(pwd)/outputs # Output directory

export LOCAL_EXECUTE=1 # Do not set when executing on ICL HPC (used by `mergehelper.sh`).
export NCPU=2 # Change this to the number of CPUs you want to use.

# If source directories (DA_EXPDIR_SRC, DA_WORKDIR_SRC) are not set,
# the experiment set-up folder will be copied from the same target directory.
# I.e. DA_EXPDIR_SRC==DA_EXPDIR and DA_WORKDIR_SRC==DA_WORKDIR.
export DA_EXPDIR_SRC=$(pwd)/u-dales/examples
export DA_WORKDIR_SRC=$(pwd)/u-dales/outputs
```


Now to set-up a new experiment based on a previous example (here we use case `999`), run:

``` sh
# General syntax: da_prep.sh new_exp_id old_exp_id
# To set up a new simulation starting from the restart files of another simulation
# ("warmstart"), use the 'w' flag. E.g.: da_prep.sh new_exp_id old_exp_id w
./u-dales/tools/utils/da_prep.sh 999 999
```

## Run

The scripts `local_execute.sh` and `hpc_execute.sh` in `u-dales/tools/utils` are used as wrappers to run simulations on your local machine and HPC at ICL respectively. These scripts contain several helpers to run the simulations and merge outputs from several CPUs into a single file (see [Post-processing](./udales-post-processing.md) for more info about the individual scripts).

``` sh
# General syntax: local_execute.sh exp_id
# To run on HPC at ICL, run `hpc_execute.sh` instead.
./u-dales/tools/utils/local_execute.sh 999
```

## What's next?

This simple guide is meant to allow you to get started based on an existing example. To learn more about pre- and post-processing steps in uDALES, please see the [Pre-processing](./udales-pre-processing.md) and [Post-processing](./udales-post-processing.md) pages.