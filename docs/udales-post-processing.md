# Post-processing

uDALES saves the outputs as NetCDF files. If a simulation is run on several processors, each processor writes independent output files. The scripts `nco_concatenate_field_x.sh` and `nco_concatenate_field_y.sh` in the `tools` directory can be used together to gather these output files into a single file.
The wrapper script `gather_outputs.sh` does this automatically for all output fields of the simulation. The script is automatically called after a simulation run when using `local_execute.sh` for executing your simulation. However on Imperial HPC or ARCHER2, the file gathering process should be executed seperately using `hpc_gather.sh` or `archer_gather.sh` respectively, once the main simulation completes.

If you have separate output files of a continuous simulation, e.g. because one simulation is the warmstart of the other simulation, you can append these output files into a single file using the script `append_outputs.sh`.

## Gather output fields

To gather the output files of serveral processors from your simulation to a single file, use:

``` sh
# We assume you are running the following commands from your
# top-level project directory.

# General syntax: gather_outputs.sh <path-to-exp-outputs>
./u-dales/tools/gather_outputs.sh outputs/009
```

In the above command, replace 009 with the number of your simulation. 

When on Imperial HPC, one should carry out the gather operation on a compute node using the wrapper script `hpc_gather.sh` as below,

``` sh
# We assume you are running the following commands from your
# top-level project directory.

./u-dales/tools/hpc_gather.sh experiments/009
```
Note in the above command `path-to-exp` is passed as an input argument instead of `path-to-exp-outputs`. The hpc_gather.sh script will identify the output work directory from the `config.sh` file available in the experiment directory.

When on the ARCHER2 cluster, one should use the wrapper script `archer_gather.sh`. In this case, one should pass `path-to-exp-outputs` as the input argument as below,

``` sh
# We assume you are running the following commands from your
# top-level project directory.

./u-dales/tools/archer_gather.sh outputs/009
```

## Append two output files

We assume that simulation 1 was run before simulation 2, i.e. the time steps of simulation 1 are all before simulation 2. To append the output files of simulation 1 (009) to simulation 2 (010), use:

``` sh
# We assume you are running the following commands from your
# top-level project directory.

# General syntax: append_outputs.sh <path-to-simulation-1-outputs> <path-to-simulation-2-outputs>
./u-dales/tools/append_outputs.sh outputs/009 outputs/010
```

Replace 009 and 010 with the numbers of your simulations.

## Different output files explained

The output files generated depend on the parameters specified under `&OUPUT` in the `namoptions` file of your simulation (see [Configuration](udales-namoptions-overview.md) for details), and the name of the output file(s) matches the name of that switch, e.g. if `lxytdump` is selected for experiment `009` then there will be an output file called `xytdump.009.nc`. If `lfielddump` is selected, note that there will be a `fielddump.xxx.009.nc` file for each cpu.

## Reading output files

These output files are in netcdf format, and so it is possible to obtain a description of any particular file using the command `ncdisp('<top-level-directory>/outputs/009/xytdump.009.nc')` in Matlab. To read a variable, one can use e.g. `u = ncread(<top-level-directory>/outputs/009/xytdump.009.nc', 'uxyt')`
