## Post-processing

The output files generated depend on the parameters specified under `&OUPUT` in [namoptions](udales-namoptions-overview.md), and the name of the output file(s) matches the name of the switch, e.g. if `lxytdump` is selected for experiment `009` then there will be an output file called `xytdump.009.nc`. If `lfielddump` is selected, note that there will be a `fielddump` file for each cpu.

These output files are in netcdf format, and so it is possible to obtain a description of any particular file using the command `ncdisp('<top-level-directory>/outputs/009/xytdump.009.nc')` in Matlab. To read a variable, one can use e.g. `u = ncread(<top-level-directory>/outputs/009/xytdump.009.nc', 'uxyt')`

The scripts `da_append.sh`, `da_concatenate.sh` in the directory `tools/utils` can be used to merge the output of netcdf files into a single netcdf file. The file `concatenate_field.sh` is also needed because it is used by these scripts.

* To concatenate the output of serveral cpus from one simulation to a single file, use:
``` sh
da_concatenate.sh <path-to-exp-outputs>
```

* To append the output files of a simulation to the output files of another simulation, use:
``` sh
da_append.sh <path-to-exp-outputs1> <path-to-exp-outputs2>
```


