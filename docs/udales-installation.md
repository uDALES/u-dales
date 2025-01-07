# Installation

Getting started with uDALES to set up your own experiments is straightforward. This guide goes through the steps required to [install](#installation) uDALES, and [set-up](#set-up) and [run](#run) a simple example. Results are output in netCDF format, for a quick inspection you can use GUI tools such as [Panoply](https://www.giss.nasa.gov/tools/panoply/) or [ncview](http://meteora.ucsd.edu/~pierce/ncview_home_page.html). To learn more about pre- and post-processing steps see the [what's next section](#whats-next).

If you have [Singularity](https://sylabs.io/) available on your system, you can use the provided scripts under `tools/singularity` to build and run uDALES cases locally or on HPC environments. See [Singularity](#singularity) for instructions; otherwise, see the next section.

## Prerequisites

### uDALES

uDALES is supported to run on Linux, macOS and Windows Subsystem for Linux (WSL). Please ensure that the latest version of the following libraries and software are available on your system:

- [CMake](https://cmake.org/) >= 3.9.
- [NetCDF-Fortran](https://www.unidata.ucar.edu/downloads/netcdf/index.jsp) >= 4.
- [GNU](https://gcc.gnu.org/wiki/GFortran) >= 9, [Intel](https://software.intel.com/content/www/us/en/develop/documentation/fortran-compiler-developer-guide-and-reference/top.html), or [Cray](https://pubs.cray.com/) Fortran compiler.
- A recent version of [MPICH](https://www.mpich.org/) or [Open-MPI](https://www.open-mpi.org/). 
- [FFTW](http://www.fftw.org/) 

### Project setup

To set up a project template for uDALES with a generic folder structure that you can later use to set up your own experiments, you will need:

- [Git](https://git-scm.com/) >= 2.
- A [GitHub](https://github.com) account. (optional)
- [Python](https://www.python.org/) >= 3.6.

### Pre-processing

When you create your own experiments, the input data needs to be processed before the simulation can be performed. For this reason, uDALES has a suite of pre-processing routines written in several languages including MATLAB. For information on how perform the pre-processing, see [pre-processing](./udales-pre-processing.md).

- [MATLAB](https://www.mathworks.com/products/matlab.html) >= R2017b

### Post-processing

The uDALES output data are stored in netCDF files. However, as the code is parallel, these need to be merged before they can be conveniently visualised an processed using Python, MATLAB or other languages. The merging is performed with the NCO package:

- [netCDF Operators](https://github.com/nco/nco) (NCO).

On local systems, these software and libraries (except MATLAB) should be available from your system's package manager (e.g. APT, yum, Homebrew, etc.) and examples on how to install all the required libraries for Linux/WSL and macOS are given below.

On high performance computing (HPC) clusters, these software and libraries should have already been installed. Please refer to the specific documentation to load the above software and libraries. Alternatively, you can install all the required packages easily after installing [Linuxbrew](https://docs.brew.sh/Homebrew-on-Linux) and using the instructions for [macOS](#macos).

## Linux/WSL (Ubuntu)

```sh
sudo apt-get update && sudo apt-get upgrade -y
sudo apt-get install -y git cmake gfortran libomp-dev libopenmpi-dev openmpi-bin libnetcdf-dev libnetcdff-dev nco python3 python3-pip libfftw3-dev
```

## macOS

On macOS, use [Homebrew](https://docs.brew.sh) to install the required libraries. If you do not have Homebrew installed on your system, install it from the [Homebrew installation page](https://docs.brew.sh/Installation) then, to install all the required dependencies, including support for MPI, run the following commands from your terminal prompt:

```sh
brew update
brew install git cmake gcc netcdf netcdf-fortran mpich nco python3 fftw
```

## Repository set-up

Create a top-level directory, for example called "uDALES":

```sh
mkdir uDALES
```

Clone the u-dales repository into the top-level directory:

```sh
cd uDALES
git clone --recurse-submodules https://github.com/uDALES/u-dales.git
```

Create directories for experiment set-ups and output data:

```sh
mkdir experiments outputs
```

such that your directory tree resembles the following:

``` sh
.
myproject
│
├── experiments # Configuration files grouped by experiment number.
│   └── <N>     # Any configurations files needed by uDALES to run experiment <N> a three digit integer number.
│   └── ...
│
├── outputs     # Additional or specialized tools other then the ones included with uDALES.
│   └── <N>     # Output from experiment <N>.
│   └── ...
│
└── u-dales     # uDALES model development repository (submodule).
│   └── 2decomp-fft
│   └── ...
│   └── src
│   └── ...
│   └── tools
│   └── ...
```

In the next steps we will assume your current working directory is the top-level project directory.


## Build on common systems

To compile uDALES (in release mode) on common/local uvuntu or mac systems using helper script, run:

```sh
# We assume you are running the following commands from the u-dales directory
tools/build_executable.sh common release
```

OR, 
you can do it manually. On standard systems and configurations, you can build uDALES with the following commands:

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

## Build on HPCs

To compile uDALES (in release mode) on the ICL HPC cluster run:

```sh
# We assume you are running the following commands from the u-dales directory
tools/build_executable.sh icl release
```

To compile uDALES (in release mode) on ARCHER2, use:

```sh
# We assume you are running the following commands from the u-dales directory
tools/build_executable.sh archer release
```

Information for developers: if you are a High Performance Cluster (HPC) user you are likely using the [Environment Modules package](http://modules.sourceforge.net/) for the dynamic modification of the user's environment via modulefiles and therefore you may need to hint CMake the PATH to netCDF (see below how).

Here we show how to compile uDALES using the [HPC at ICL](https://www.imperial.ac.uk/admin-services/ict/self-service/research-support/rcs/) as an example, therefore please note that the specific names/versions installed on your system may be different.

``` sh
module list # list currently enabled modules -- should be empty!
module avail # list available modules
```

``` sh
# This is an example, please check with the previous command for the exact name of the
# modules available on your system. This will load netCDF compiled with Intel Suite
# 2020.2 and add the correct version of icc and ifort to the PATH.
module load intel-suite/2020.2 mpi/intel-2019.8.254 cmake/3.18.2 git/2.14.3
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

where `NETCDF_DIR` and `NETCDF_FORTRAN_DIR` indicates the absolute path to your netCDF-C and netCDF-Fortran installation directories. Here, we use the utilities `nc-config` and `nf-config` to hint CMake the location of netCDF, but you can simply pass the absolute path to the netCDF-C and netCDF-Fortran manually instead. You can compile in parallel mode by passing Make the `j` flag followed by the number of CPU cores to use. For example, to compile with 2 cores do `make -j2`.

## Build defaults/options

By default uDALES will compile in `Release` mode. You can change this by specifying the option (or flag) at configure time. The general syntax for specifying an option in CMake is `-D<flag_name>=<flag_value>` where `<flag_name>` is the option/flag name and `<flag_value>` is the option/flag value. The following options can be specified when configuring uDALES:

| Name                            | Options            | Default   | Description                                   |
| ------------------------------- | ------------------ | --------- | --------------------------------------------- |
| `CMAKE_BUILD_TYPE`              | `Release`, `Debug` | `Release` | Whether to optimise/build with debug flags    |
| `NETCDF4_DIR`                   | `<path>`           | -         | Path to netCDF-C installation directory       |
| `NETCDF_FORTRAN_DIR`            | `<path>`           | -         | Path to netCDF-Fortran installation directory |
| `SKIP_UPDATE_EXTERNAL_PROJECTS` | `ON`, `OFF`        | `OFF`     | Whether to skip updating external projects    |

# Singularity

If you are looking for information on how to install or use Singularity on your system, please refer to the [Singularity documentation ](https://sylabs.io/docs). The use of Singularity is undoubtedly the easiest way to build and run cases in uDALES as all dependencies are provided and uDALES will compile out of the box. Furthermore, users wishing to achieve a reasonable level of scientific reproducibility may archive software, tools, and data with their Singularity image containing OS and external libraries to an open access repository (e.g. [Meyer et al., 2020](https://doi.org/10.1029/2019MS001961)).

First clone the uDALES repository with:

```sh
git clone https://github.com/uDALES/u-dales.git
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


## What's next?

This simple guide is meant to allow you to get started based on an existing example. We set up several [example simulations](./udales-example-simulations.md) to get you started. To learn more about pre- and post-processing steps in uDALES, please see the [pre-processing](./udales-pre-processing.md) and [post-processing](./udales-post-processing.md) pages.
