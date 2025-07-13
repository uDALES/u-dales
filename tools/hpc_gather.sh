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
 echo "The output case directory must be set."
 echo "usage: FROM THE TOP LEVEL DIRECTORY run: bash ./u-dales/tools/hpc_gather.sh <PATH_TO_OUTPUT_CASE>"
 exit 1
fi

## go to output case directory
pushd $1
outdir=$(pwd)

## set experiment number via path
exp="${outdir: -3}"

echo "experiment number: $exp"

## read in additional variables
if [ -f config.sh ]; then
    source config.sh
else
    echo "config.sh must be there inside $outdir"
    exit 1
fi

## check if required variables are set
if [ -z $DA_TOOLSDIR ]; then
    echo "Script directory DA_TOOLSDIR must be set inside $outdir/config.sh"
    exit 1
fi;
if [ -z $NNODE ]; then
    echo "Script directory NNODE must be set inside $outdir/config.sh"
    exit 1
fi;
if [ -z $NCPU ]; then
    echo "Script directory NCPU must be set inside $outdir/config.sh"
    exit 1
fi;
if [ -z $WALLTIME ]; then
    echo "Script directory WALLTIME must be set inside $outdir/config.sh"
    exit 1
fi;
if [ -z $MEM ]; then
    echo "Script directory MEM must be set inside $outdir/config.sh"
    exit 1
fi;

## write post-job.exp file for HPC
cat <<EOF > post-job.$exp
#!/bin/bash
#PBS -l walltime=${WALLTIME}
#PBS -l select=${NNODE}:ncpus=${NCPU}:mem=${MEM}
module load intel/2024a netCDF/4.9.2-iimpi-2023a netCDF-Fortran/4.6.1-iimpi-2023a NCO/5.2.9-foss-2024a GSL/2.8-GCC-13.3.0
$DA_TOOLSDIR/gather_outputs.sh $outdir
EOF

## submit post-job.exp file to queue
qsub post-job.$exp

echo "post-job.$exp submitted."

popd
