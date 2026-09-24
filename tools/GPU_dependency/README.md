# GPU dependency stacks for uDALES on Imperial HX1

The GPU build (`tools/build_executable.sh gpuhx1`) needs HDF5, netCDF-C,
netCDF-Fortran and FFTW compiled with the NVIDIA HPC SDK against one specific
MPI, and the job must launch with that MPI's `mpirun`. Three such stacks have
been built; each subdirectory holds the install scripts and a README with the
run order. Only the scripts are kept here, not the sources or installs.

| Directory | MPI | Use |
|---|---|---|
| `stack1-NVHPC-OpenMPI-3.1.5` | OpenMPI 3.1.5 bundled with NVHPC 23.7 | single node; the original stack |
| `stack2-NVHPC-OpenMPI-4.1.5` | OpenMPI 4.1.5 bundled with NVHPC 23.7 (UCX, CUDA) | single node; ~20 % faster than stack 1 |
| `stack3-OpenMPI-4.1.6-tm` | OpenMPI 4.1.6 built here with IPv6, PBS launch and the site UCX-CUDA | single and multi-node; the current default |

The bundled OpenMPIs cannot run across HX1 nodes: the node networks are
IPv6-only and those builds lack IPv6 support, and their bundled UCX cannot use
the InfiniBand cards. Stack 3 exists for that reason.

Where the stack is selected: the `gpuhx1` block of `tools/build_executable.sh`
(MPI wrapper and library prefixes), the `hx1:gpu` block of
`tools/hpc_execute.sh` (runtime `mpirun` and UCX environment) and the GPU
launcher lines of `tests/hpc_run_all_tests.sh`. All three keep the lines of
the other stacks as comments so a stack can be switched by moving the
comment markers, always in all three files together.
