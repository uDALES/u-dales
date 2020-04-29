## Post-processing

The scripts `da_append.sh`, `da_concatenate.sh` in the directory `tools` can be used to merge the output of netcdf files into a single netcdf file. The file `da_merge.sh` is also needed because it is used by these scripts.

* To concatenate the output of serveral cpus from one simulation to a single file, use:
``` sh
da_concatenate.sh <path-to-exp-outputs>
```
* To append the output files of a simulation to the output files of another simulation, use:
``` sh
da_append.sh <path-to-exp-outputs1> <path-to-exp-outputs2>
```