#!/usr/bin/env bash

# uDALES (https://github.com/uDALES/u-dales).

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Copyright (C) 2016-2019 the uDALES Team.

set -e

if (( $# < 1 ))
then
 echo "The experiment directory must be set."
 echo "usage: FROM THE TOP LEVEL DIRECTORY run: u-dales/tools/hpc_execute.sh <PATH_TO_CASE>"
 exit 1
fi

## go to experiment directory
pushd $1
inputdir=$(pwd)

## set experiment number via path
exp="${inputdir: -3}"

echo "experiment number: $exp"

## read in additional variables
if [ -f config.sh ]; then
    source config.sh
else
    echo "config.sh must be there inside $inputdir"
    exit 1
fi

## check if required variables are set
if [ -z $DA_WORKDIR ]; then
    echo "Script directory DA_WORKDIR must be set inside $inputdir/config.sh"
    exit 1
fi;
if [ -z $DA_BUILD ]; then
    echo "Script directory DA_BUILD must be set inside $inputdir/config.sh"
    exit 1
fi;
if [ -z $DA_TOOLSDIR ]; then
    echo "Script directory DA_TOOLSDIR must be set inside $inputdir/config.sh"
    exit 1
fi;
if [ -z $NNODE ]; then
    echo "Script directory NNODE must be set inside $inputdir/config.sh"
    exit 1
fi;
if [ -z $NCPU ]; then
    echo "Script directory NCPU must be set inside $inputdir/config.sh"
    exit 1
fi;
if [ -z $WALLTIME ]; then
    echo "Script directory WALLTIME must be set inside $inputdir/config.sh"
    exit 1
fi;
if [ -z $MEM ]; then
    echo "Script directory MEM must be set inside $inputdir/config.sh"
    exit 1
fi;

## set the output directory
outdir=$DA_WORKDIR/$exp

echo "writing job.$exp."

## write new job.exp file for HPC
cat <<EOF > job.$exp
#!/bin/bash
#PBS -l walltime=${WALLTIME}
#PBS -l select=${NNODE}:ncpus=${NCPU}:mpiprocs=$(( $NCPU * $NNODE )):mem=${MEM}
module load intel/2025a netCDF/4.9.2-iimpi-2023a netCDF-Fortran/4.6.1-iimpi-2023a FFTW/3.3.9-intel-2021a CMake/3.29.3-GCCcore-13.3.0 git/2.45.1-GCCcore-13.3.0
mkdir -p $outdir
cp -r $inputdir/* $outdir
pushd $outdir
mpiexec -n $(( $NCPU * $NNODE )) $DA_BUILD $outdir/namoptions.$exp > $outdir/output.$exp 2>&1
module load NCO/5.2.9-foss-2024a
$DA_TOOLSDIR/gather_outputs.sh $outdir
EOF

## submit job.exp file to queue
qsub job.$exp

echo "job.$exp submitted."
