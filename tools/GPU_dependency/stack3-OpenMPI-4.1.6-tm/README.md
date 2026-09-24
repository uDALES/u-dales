# Stack 3: OpenMPI 4.1.6 built for NVHPC 23.7 (multi-node capable)

The stack for multi-node GPU runs on Imperial HX1. It has two parts, built in
this order:

1. `openmpi/`: OpenMPI 4.1.6 compiled with `nvc`/`nvfortran` 23.7, configured
   with `--enable-ipv6` (HX1's node networks are IPv6-only), `--with-tm=/opt/pbs`
   (PBS launches the daemons, so no hostfile or ssh is needed), the site
   `UCX/1.14.1-GCCcore-12.3.0` with `UCX-CUDA`, PMIx, hwloc and libevent. This
   follows the site's own `OpenMPI/4.1.4-NVHPC-22.7-CUDA-11.7.0` recipe.
2. the dependency libraries, same scripts as stacks 1 and 2 but with
   `NVHPC_MPI` pointing at the OpenMPI installed in step 1.

## Step 1: OpenMPI

```bash
cd openmpi
curl -fL -o openmpi-4.1.6.tar.bz2 https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.6.tar.bz2
tar -xjf openmpi-4.1.6.tar.bz2
mkdir -p logs
./s2_configure.sh          # sources s1_set_env.sh; log in logs/configure.log
./s3_build_install.sh      # make, make install into ./openmpi, prints the components
```

`s1_set_env.sh` sets the prefix to `<this directory>/openmpi`. Two things the
scripts already handle and that will bite anyone re-deriving them: every flag
set needs `-fPIC` (nvfortran otherwise leaves one `mpi_f08` object
non-relocatable and the shared bindings fail to link), and `libmpi` needs
`libpbs.so.0` from `/opt/pbs/lib`, which is not on the system loader path,
hence `LDFLAGS=-Wl,-rpath,/opt/pbs/lib`. `s3` ends by printing the `plm:tm`,
`ras:tm`, `pml:ucx`, `btl:smcuda` and `coll:cuda` components, which must all
be present.

## Step 2: the libraries

Edit the `NVHPC_MPI` line of the four `*_s1_set_env_for_install.sh` scripts to
`<prefix from step 1>/bin` (they carry the path of the original build), then
run the same sequence as stacks 1 and 2:


Copy this directory to where the libraries should live (the scripts place
sources and installs beside themselves), then:

```bash
./hdf5_s0_download.sh
source ./hdf5_s1_set_env_for_install.sh
./hdf5_s2_configure_install.sh
./hdf5_s3_install.sh

./netcdf_c_s0_download.sh
./netcdf_c_s2_configure_install.sh
./netcdf_c_s3_build_install.sh

./netcdf_fortran_s0_download.sh
./netcdf_fortran_s2_configure_install.sh
./netcdf_fortran_s3_build_install.sh

./fftw_s0_download.sh
./fftw_s2_configure_install.sh
./fftw_s3_build_install.sh
```

`BUILD_JOBS=16` raises the parallelism of the `s3` steps. If `wget` is
missing, fetch the archives with `curl -fL -o <archive> <url>` and untar them
by hand. Installed prefixes: `hdf5-2.1.1/hdf5`, `netcdf-c-4.10.1/netcdfc`,
`netcdf-fortran-4.6.4/netcdff`, `fftw-3.3.11/fftw3`; verify with
`nc-config --has-parallel4` (yes) and `grep "Parallel HDF5" hdf5-2.1.1/hdf5/lib/libhdf5.settings`.

## Using the result

Point the `gpuhx1` block of `tools/build_executable.sh` at the installed
prefixes (`NETCDF_DIR`, `NETCDF_FORTRAN_DIR`, `FFTWDIR`) and the matching MPI
wrapper directory, and the `hx1:gpu` block of `tools/hpc_execute.sh` at the
matching `mpirun`. The three must belong to the same stack: the executable,
2DECOMP&FFT, these libraries and the launcher all share one `libmpi`.

## Runtime

Besides the `mpirun` of step 1, jobs need UCX's CUDA transports, which the
site keeps in a separate UCX-CUDA install and which load only with:

```bash
export UCX_MODULE_DIR=/gpfs/easybuild/prod/software/UCX-CUDA/1.14.1-GCCcore-12.3.0-CUDA-12.1.1/ucx
export EB_UCX_uct_MODULES=":ib:rdmacm:cma:cuda"
export EB_UCX_ucm_MODULES=":cuda"
export EB_UCX_uct_cuda_MODULES=":gdrcopy"
export LD_LIBRARY_PATH="$EBROOTNVHPC/Linux_x86_64/23.7/cuda/12.2/lib64:/gpfs/easybuild/prod/software/GDRCopy/2.3.1-GCCcore-12.3.0/lib:$LD_LIBRARY_PATH"
```

Without them UCX treats GPU buffers as host memory and multi-rank runs
segfault on small messages. `tools/hpc_execute.sh` writes this block into
every GPU job, and launches with `mpirun --map-by ppr:<NGPU>:node`. Limits and
measured performance are in
`.github/skills/udales-exec/references/clusters.md` (HX1 GPU section): a user
may run at most 12 GPUs at once, multi-node does not scale at 512^3, and the
1024^3 case runs on 8 GPUs over 2 nodes at about 5 s per step.
