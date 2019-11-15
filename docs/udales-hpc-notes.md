# Notes on using uDALES with HPCs

## Build on HPCs

If you are a High Performance Cluster (HPC) user you are likely using the [Environment Modules package](http://modules.sourceforge.net/) for the dynamic modification of the user's environment via modulefiles and therefore you may need to hint CMake the PATH to NetCDF (see below how).

Here we show how to compile uDALES using the [HPC at ICL](https://www.imperial.ac.uk/admin-services/ict/self-service/research-support/rcs/) as an example, therefore please note that the specific names/versions installed on your system may be different.

``` sh
# This is an example, module names/versions may be different on your system
module list # list currently enabled modules -- should be empty!
module avail # list available modules
# This is an example, please check with the previous command for the exact name of the modules available on your system. This will load NetCDF compiled with Intel Suite 2019.4 and add the correct version of icc and ifort to the PATH.
module load intel-suite/2017.6 mpi/intel-2018 cmake/3.14.0 git/2.14.3
```

Then, to build the uDALES executable, run the following commands:

``` sh
# We assume you are running the following commands from your
# top-level project directory.
mkdir -p u-dales/build/release
pushd u-dales/build/release
FC=mpiifort cmake -DNETCDF_DIR=/apps/netcdf/4.4.1-c -DNETCDF_FORTRAN_DIR=/apps/netcdf/4.4.4-fortran ../..
make
popd
```

where `NETCDF_DIR` and `NETCDF_FORTRAN_DIR` indicates the absolute path to your NetCDF-C and NetCDF-Fortran installation directories. Here, we use the utilities `nc-config` and `nf-config` to hint CMake the location of NetCDF, but you can simply pass the absolute path to the NetCDF-C and NetCDF-Fortran manually instead. You can compile in parallel mode by passing Make the `j` flag followed by the number of CPU cores to use. For exmaple, to compile with 2 cores do `make -j2`.


## Run on HPC

Here we show how to execute uDALES using the [HPC at ICL](https://www.imperial.ac.uk/admin-services/ict/self-service/research-support/rcs/) as an example, therefore please note that the queuing system on your system may be set up differently.
The PBS queuing system at ICL requires you to submit a job file for your simulations. It specifies the number of processors, the amount of time the simulation needs to run, etc.
The script `hpc_execute.sh` in `u-dales/tools/utils` is used as wrapper to write a job file and run simulations with the queuing system on the HPC at ICL.
First you need to go to your experiment folder, e.g.:

``` sh
# We assume you are running the following commands from your
# top-level project directory.

cd experiments/009
```

Copy the wrapper script to the current directory using:

``` sh
cp ../../u-dales/tools/utils/hpc_execute.sh .
```

Next you can modify parameters set within the script. You may want to adapt the numbers of processors and nodes used for the experiments, or modify the walltime. For guidance on how to set the parameters, see [Job sizing guidance](https://www.imperial.ac.uk/admin-services/ict/self-service/research-support/rcs/computing/job-sizing-guidance/).

To queue the simulation for execution, run:

``` sh
./hpc_execute.sh
```

You can check the status of the simulation by the command:

``` sh
qstat
```

## Project setup on several machines

If you want to use the same project setup on several machines, e.g. both locally on your computer and on a High Performance Cluster (HPC), you need to clone the project git repository including its submodules with:

``` sh
git clone --recurse-submodules https://github.com/<MY_GITHUB_USERNAME>/<PROJECT_REPO>.git
```

The submodules will be in a detached HEAD state, you will need to checkout a working branch, e.g.

``` sh
cd u-dales
git checkout master
```
