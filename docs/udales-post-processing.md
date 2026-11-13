# Post-processing

uDALES saves the outputs as NetCDF files. If a simulation is run on several processors, each processor writes independent output files. The scripts `nco_concatenate_field_x.sh` and `nco_concatenate_field_y.sh` in the `tools` directory can be used together to gather these output files into a single file.
The wrapper script `gather_outputs.sh` does this automatically for all output fields of the simulation. It is called automatically at the end of a local run (`ud_run common sim`). On clusters (Imperial HPC, ARCHER2, ECMWF), the gathering is submitted as a separate job with `ud_run <machine> gather` once the main simulation completes.

If you have separate output files of a continuous simulation, e.g. because one simulation is the warmstart of the other simulation, you can append these output files into a single file using the script `append_outputs.sh`.

## Gather output fields

To gather the output files of several processors from your simulation into single files, use `ud_run` with the same `<machine>` and experiment directory you used for the simulation:

``` sh
# We assume you are running the following commands from your
# top-level project directory.

# General syntax: ud_run <machine> gather <path-to-exp>
./u-dales/bin/ud_run icl gather experiments/009
```

In the above command, replace `icl` with your machine (`common`, `icl`, `hx1`, `archer` or `ecmwf`) and 009 with the number of your simulation. The experiment directory (`path-to-exp`) is passed, not the output directory: the output directory `$DA_WORKDIR/009` is found from the `config.sh` file in the experiment directory. Depending on the machine, `ud_run` does the following:

| `<machine>` | What `gather` does |
|---|---|
| `common` | Runs `gather_outputs.sh` directly on `$DA_WORKDIR/<exp>`, appending to `output.<exp>.log`. Only needed to redo a gather, since `ud_run common sim` already gathers. |
| `icl`, `hx1` | Submits a PBS job with `hpc_gather.sh`. |
| `archer` | Submits a Slurm job with `archer_gather.sh`, passing it `$DA_WORKDIR/<exp>`. |
| `ecmwf` | Submits a Slurm job with `ecmwf_gather.sh`, always on the shared `nf` queue (it ignores `QOS`). `MEM` (`16G` or `16gb`) and `WALLTIME` must be set in `config.sh`; see [Run on ECMWF HPC2020](./udales-simulation-setup.md#run-on-ecmwf-hpc2020). |

The gather step can also be run without `ud_run`, by calling `gather_outputs.sh` on the output directory (on clusters, do this on a compute node, not a login node):

``` sh
# General syntax: gather_outputs.sh <path-to-exp-outputs>
./u-dales/tools/gather_outputs.sh outputs/009
```

Note that `archer_gather.sh`, when called directly, takes the output directory (`path-to-exp-outputs`), while `hpc_gather.sh` and `ecmwf_gather.sh` take the experiment directory.

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

The output files generated depend on the parameters specified under `&OUTPUT` in the `namoptions` file of your simulation (see [Configuration](udales-namoptions-overview.md) for details), and the name of the output file(s) matches the name of that switch, e.g. if `lxytdump` is selected for experiment `009` then there will be an output file called `xytdump.009.nc`. If `lfielddump` is selected, note that there will be a `fielddump.xxx.009.nc` file for each cpu.

## Reading output files

These output files are in netcdf format, and so it is possible to obtain a description of any particular file using the command `ncdisp('<top-level-directory>/outputs/009/xytdump.009.nc')` in Matlab. To read a variable, one can use e.g. `u = ncread(<top-level-directory>/outputs/009/xytdump.009.nc', 'uxyt')`
