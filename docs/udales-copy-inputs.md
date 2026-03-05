# Copying simulation input files

To set up a new simulation, `copy_inputs.sh` in `u-dales/tools/` is used to create a new simulation setup `new_exp_id` based on another simulation `old_exp_id`. All `exp_ids` are three digit integer numbers, e.g. 001, and are stored in directories of that name. Each experiment case directory must contain a config.sh file where appropriate paths for DA_EXPDIR (experiments directory), DA_WORKDIR (outputs directory), DA_TOOLSDIR (u-dales/tools directory) are set using export.

<!--
Scripts requires several variables to be set up. You can do this by copying and pasting the snippet below or by including it in a bash script (or bash profile if you are unlikely to change them).

``` sh
# We assume you are running the following commands from your
# top-level project directory.

export DA_EXPDIR=$(pwd)/experiments #  The top-level directory of the simulation setups.
export DA_WORKDIR=$(pwd)/outputs # Output top-level directory

# If source directories (DA_EXPDIR_SRC, DA_WORKDIR_SRC) are not set,
# the experiment set-up folder will be copied from the same target directory.
# I.e. DA_EXPDIR_SRC==DA_EXPDIR and DA_WORKDIR_SRC==DA_WORKDIR.
export DA_EXPDIR_SRC=$(pwd)/u-dales/examples
export DA_WORKDIR_SRC=$(pwd)/u-dales/examples
```

If you set up a new experiment on HPC, also use:

``` sh
export DA_WORKDIR=$EPHEMERAL # Output top-level directory on HPC
export DA_WORKDIR_SRC=$EPHEMERAL
```
-->

Now to set-up a new experiment (here we use case `009`) based on a previous example (here we use case `001`), run:

``` sh
# We assume you are running the following commands from your
# top-level project directory.

# General syntax: copy_inputs.sh old_exp_id new_exp_id
u-dales/tools/copy_inputs.sh experiments/001 009

# To set up a new simulation starting from the restart files of another simulation
# ("warmstart"), use the 'w' flag. E.g.: copy_inputs.sh old_exp_id new_exp_id w
u-dales/tools/copy_inputs.sh experiments/001 009 w
```
