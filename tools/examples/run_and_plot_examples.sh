#!/usr/bin/env bash

# Runs and plots cases in the examples folders.
# Data and plots are saved under outputs.
# Usage: ./tools/examples/run_and_plot_examples.sh

set -e

if [ ! -d src ]; then
  echo "Please run this script from the root folder"
  exit 1
fi

export DA_TOOLSDIR=$(pwd)/tools
export DA_BUILD=$(pwd)/build/release/u-dales
export NCPU=4
export DA_WORKDIR=$(pwd)/outputs

for example in 001 002 101 102 201 501 502
do
    if [[ $example == 102 ]]; then
        # Download required files for warmstart simulation
        mkdir -p $DA_WORKDIR/$example
        pushd $DA_WORKDIR/$example
        curl -o examples_warmstart_102.zip -L https://www.dropbox.com/sh/20rsgpt0gh09gr7/AABuoCFtn6_zFTxx4k8pKqvLa?dl=1
        set +e # Unzip may raise warnings as errors
        unzip -o examples_warmstart_102.zip
        set -e
        popd
    fi

    if [[ $example == 502 ]]; then
        # Download required files for driver simulation
        mkdir -p $DA_WORKDIR/$example
        pushd $DA_WORKDIR/$example
        curl -o examples_driver_501.zip -L https://www.dropbox.com/sh/spld3hqipqe17j1/AAA0cuzW3qc9ftY6dvHcSSL8a?dl=1
        set +e # Unzip may raise warnings as errors
        unzip -o examples_driver_501.zip
        set -e
        popd
    fi

    ./tools/local_execute.sh examples/$example

done

eval "$(conda shell.bash hook)"
conda activate udales
python tools/examples/plot_examples.py