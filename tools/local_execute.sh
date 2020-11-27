#!/usr/bin/env bash

# Usage: ./tools/local_execute.sh <PATH_TO_CASE>

set -e

if (( $# < 1 ))
then
    echo "The path to case folder must be set."
    exit
fi

## go to experiment directory
pushd $1
inputdir=$(pwd)

## set experiment number via path
exp="${inputdir: -3}"

echo "Setting up uDALES for case $exp..."

## read in additional variables
if [ -f config.sh ]; then
    source config.sh
fi

## check if required variables are set
## or set default if not
if [ -z $NCPU ]; then
    NCPU=1
fi;
if [ -z $DA_WORKDIR ]; then
    echo "Output top-level directory DA_WORKDIR must be set"
    exit
fi;
if [ -z $DA_BUILD ]; then
    echo "Executable DA_BUILD must be set"
    exit
fi;
if [ -z $DA_TOOLSDIR ]; then
    echo "Script directory DA_TOOLSDIR must be set"
    exit
fi;

## set the experiment output directory
outdir=$DA_WORKDIR/$exp

echo "Starting job for case $exp..."

## always start afresh
rm -rf $outdir/*.log $outdir/*.nc 

## copy files to output directory
mkdir -p $outdir
cp ./* $outdir

## go to execution and output directory
pushd $outdir

## execute program with mpi
mpiexec -n $NCPU $DA_BUILD namoptions.$exp 2>&1 | tee -a run.$exp.log

## Merge output files across outputs.
if (($NCPU > 1 )); then
    echo "Merging outputs across cores into one..."
    export LOCAL_EXECUTE=1
    $DA_TOOLSDIR/gather_outputs.sh $outdir
fi

popd

echo "Simulation for case $exp ran sucesfully!"
