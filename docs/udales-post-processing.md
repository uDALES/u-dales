## Post-processing

The scripts `da_merge.sh` and `mergehelper.sh` in the directory `tools/utils` can be used to merge the output of netcdf files into a single netcdf file.
This can be used to concatenate the output of serveral cpus from one simulation to a single file, or to append the output files of a simulation to the output files of another simulation.

* Concatenate output files of a single simulation:
``` sh
mergehelper.sh exp_id
```
* Appending output files of two simulations:
``` sh
mergehelper.sh exp1_id exp2_id
```