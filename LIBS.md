# Required libraries for uDALES
The following libraries are required on your system to install uDALES from source: [Git](https://git-scm.com/),  [HDF5](https://support.hdfgroup.org/HDF5/), [NetCDF-C](https://www.unidata.ucar.edu/downloads/netcdf/index.jsp), [NetCDF-Fortran](https://www.unidata.ucar.edu/downloads/netcdf/index.jsp), and MPI.

For the installation of libraries, we rely on your system's package managers (APT, YUM, Homebrew). We assume that you have administrative privileges on your computer. If you do not have administrative privileges you should request these libraries to be installed by your system administrator so that they can tare of this for you and manage updates on your behalf.


## Ubuntu

To install all the required dependencies, including support for MPI, run the following commands from your terminal prompt:

```sh
sudo apt-get update && sudo apt-get upgrade -y
sudo apt-get install -y git cmake gfortran libnetcdf-dev libnetcdff-dev libmpich-dev
```

After the installation is complete, you can go back to installing uDALES.

## macOS

On macOS, we can use [Homebrew](https://docs.brew.sh) to to install the required libraries. If you do not have Homebrew installed on your system, install it from the [Homebrew installation page](https://docs.brew.sh/Installation) then, to install all the required dependencies, including support for MPI, run the following commands from your terminal prompt:

```sh
brew update
brew install git cmake gcc netcdf jasper open-mpi
```

After the installation is complete, you can go back to installing uDALES.
