# uDALES [![Build Status](https://travis-ci.com/uDALES/u-dales.svg?token=3tqUbxqJuLtozjxqDymC&branch=master)](https://travis-ci.com/uDALES/u-dales)

TODO: uDALES is...

## Installing

See the [installation instructions](INSTALL.md)

## Running

A uDALES simulation can be manually started from a directory that contains all input files. In this example this is `exp/006`. Use

``` sh
cd exp/006
mpiexec -n NPROC ../../build/u-dales namoptions.006 > output.006 2>&1
```

Where `NPROC` is the number of CPUs to use in the simulation (e.g. 2).

Two scripts are available in the directory `tools/utils` that will set this up automatically. The script `local_execute.sh` is for execution on your local machine and `hpc_execute.sh` for simulations on the ICL HPC. Copy the script into the simulation directory `exp/exp_id` and execute it from there:

``` sh
cd exp/exp_id
cp ../../tools/utils/local_execute.sh .
./local_execute.sh
```

Changes to the execution setup, e.g. changing the number of cpus being used, can be made directly to the local copy of these files.
After the simulation the scripts are set up to automatically call a script which merges the output of several cpus into one file (see [Post-processing](#Post-processing)).

## Pre-processing

1) Setting up a new simulation: The script `da_prep.sh` in `tools/utils` can be used to create a new simulation setup `new_exp_id` based on another simulation `old_exp_id`. All exp_ids are three digit numbers, e.g. 001, and are stored in directories of that name. The script requires you to have the variables `${DA_EXPDIR}` and `${DA_WORKDIR}` set up, e.g. via your bash_profile. `${DA_EXPDIR}` points to the top-level directory of the simulation setups (here this is the folder `exp`) and `${DA_WORKDIR}` is your output directory.

Usage: 

``` sh
da_prep.sh new_exp_id old_exp_id
```

To set up a new simulation starting from the restart files of another simulation ("warmstart"), use 
``` sh
da_prep.sh new_exp_id old_exp_id w
```

2) Updating input files: For some changes in the simulation setup (e.g. domain or resolution, blocks, initial wind speed) the input files need to be reconfigured.

TODO: Describe pre-processing scripts

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

## Documentation

WIP at: https://uDALES.github.io/u-dales/0YiO263pFxExSdkMvWfId3qkVUSF4dREFnwM1jQD9y1KvzeAVAWzGykQemUrkJCM/html/index.html

## Testing

TODO: add info about current testing and CI.


## Versioning

TODO: need create first release as soon as cleanup is complete.

## Copyright and License

FIXME: need authorship and year

This project is licensed under the GNU General Public License - see the [LICENSE.md](LICENSE.md) file for details