## Installing

First, clone or download from the GitHub repository at https://github.com/uDALES/u-dales.

```
git clone https://github.com/uDALES/u-dales.git --branch dmey/cmake
cd u-dales
```

All the steps below assume that you are in the `u-dales` directory. If you are building on the HPC at ICL, please see [Note for HPC users at ICL](#Note-for-HPC-users-at-ICL)

### Build instructions for Linux and macOS

The following libraries are required on your system to install uDALES from source: [Git](https://git-scm.com/), [HDF5](https://support.hdfgroup.org/HDF5/), [NetCDF-C](https://www.unidata.ucar.edu/downloads/netcdf/index.jsp), [NetCDF-Fortran](https://www.unidata.ucar.edu/downloads/netcdf/index.jsp), and an MPI library implementation. The above libraries should be avalable from your system's package manager (e.g. APT, yum, Homebrew, etc.).  If you do not have the latest version of these libraries installed on your system, please see the [Libraries page](LIBS.md).

To build the uDALES executable, from the `u-dales` repository, run the following commands:

```sh
$ mkdir build && cd build
$ cmake -DCMAKE_BUILD_TYPE=Release ..
$ make -j$(nproc)
```

Where `$(nproc)` will use all the number of CPU cores/threads available on your system. Note that using the maximum number of CPU cores/threads available may not necessarily be the fastest way to build the software.


### Build instructions for HPC users

If you are an HPC user, you are most likely using the [Environment Modules package](http://modules.sourceforge.net/) for the dynamic modification of the user's environment via modulefiles. In that case, you will need to specify the path to the NetCDF manually _after_ loading all the libraries required to compile uDALES. For example:

``` sh
# This is an example, module names/versions may be different on your system
module list # list currently enabled modules
module avail # list available modules
module load cmake netcdf4 openmpi gnu # This is an example, please check with the previous command for the exact name of the modules available on your system.
```

If you do not know the location of NetCDF, you can locate it with the `nc-config --prefix`.

``` sh
cmake -DCMAKE_INSTALL_PREFIX=install -DNETCDF_DIR=$NETCDF4_DIR -DNETCDF_FORTRAN_DIR=$DNETCDF_FORTRAN_DIR ..
```

where `$NETCDF4_DIR` and `$DNETCDF_FORTRAN_DIR` are the absolute path to your NetCDF-C and NetCDF-Fortran installation directories.


#### Note for HPC users at ICL

There is an issue with how the libraries are loaded and added to the `PATH`. This issue has been reported and should be resolved soon.
For the moment, if you want to compile uDALES on the HPC at ICL, please load use the following workaround instead:


``` sh
module load intel-suite/2017.6 mpi/intel-2018 cmake/3.14.0 git/2.14.3
mkdir build && cd build
FC=mpiifort cmake -DNETCDF_DIR=/apps/netcdf/4.4.1-c -DNETCDF_FORTRAN_DIR=/apps/netcdf/4.4.4-fortran ..
make -j$(nproc)
```
