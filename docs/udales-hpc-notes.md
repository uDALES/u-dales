# Notes on using uDALES with HPCs

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

## Run on HPC

Here we show how to execute uDALES using the [HPC at ICL](https://www.imperial.ac.uk/admin-services/ict/self-service/research-support/rcs/) as an example, therefore please note that the queuing system on your system may be different.

The script `hpc_execute.sh` in `u-dales/tools/utils` is used as wrapper to run simulations with the queuing system on the HPC at ICL.

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

Next you can modify parameters set within the script. You may want to adapt the numbers of processors and nodes used for the experiments, or modify the walltime.
For guidance on how to set the parameters, see [Job sizing guidance](https://www.imperial.ac.uk/admin-services/ict/self-service/research-support/rcs/computing/job-sizing-guidance/).

To queue the simulation for execution, run:

``` sh
./hpc_execute.sh
```

You can check the status of the simulation by the command:

``` sh
qstat
```