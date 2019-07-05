# Required libraries for uDALES

This is a helper guide on how to install all the required libraries to build uDALES from source on Ubuntu or macOS. We assume that you have administrative privileges on your computer. If you do not have administrative privileges you should request these libraries to be installed by your system administrator so that they can tare of this for you and manage updates on your behalf.

## Ubuntu

To install all the required dependencies, including support for MPI, run the following commands from your terminal prompt:

```sh
sudo apt-get update && sudo apt-get upgrade -y
sudo apt-get install -y git cmake gfortran libopenmpi-dev openmpi-bin libnetcdf-dev libnetcdff-dev
```

After the installation is complete, continue with installing uDALES.

## macOS

On macOS, we can use [Homebrew](https://docs.brew.sh) to to install the required libraries. If you do not have Homebrew installed on your system, install it from the [Homebrew installation page](https://docs.brew.sh/Installation) then, to install all the required dependencies, including support for MPI, run the following commands from your terminal prompt:

```sh
brew update
brew install git cmake gcc netcdf open-mpi
```

After the installation is complete, continue with installing uDALES.
