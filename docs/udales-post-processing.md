# Post-processing

uDALES saves the outputs as NetCDF files. If a simulation is run on several processors, each processor writes independent output files. The script `concatenate_field.sh` in the `tools` directory can be used to gather these output files into a single file.
The wrapper script `da_concatenate.sh` does this automatically for all output fields of the simulation. The script is automatically called after a simulation run when using `local_execute.sh` or `hpc_execute.sh` for executing your simulation.

If you have separate output files of a continuous simulation, e.g. because one simulation is the warmstart of the other simulation, you can append these output files into a single file using the script `da_append.sh`.

## Gather output fields

To gather the output files of serveral processors from your simulation to a single file, use:

``` sh
# We assume you are running the following commands from your
# top-level project directory.

# General syntax: da_concatenate.sh <path-to-exp-outputs>
./u-dales/tools/da_concatenate.sh outputs/009
```

Replace 009 with the number of your simulation.

## Append two output files

We assume that simulation 1 was run before simulation 2, i.e. the time steps of simulation 1 are all before simulation 2. To append the output files of simulation 1 (009) to simulation 2 (010), use:

``` sh
# We assume you are running the following commands from your
# top-level project directory.

# General syntax: da_append.sh <path-to-simulation-1-outputs> <path-to-simulation-2-outputs>
./u-dales/tools/da_append.sh outputs/009 outputs/010
```

Replace 009 and 010 with the numbers of your simulations.
