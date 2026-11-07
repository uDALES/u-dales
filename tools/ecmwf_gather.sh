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

## Submit the gather_outputs.sh post-processing of a finished uDALES run on the
## ECMWF Atos HPC2020 (Slurm).

set -e

if (( $# < 1 ))
then
 echo "The experiment directory must be set."
 echo "usage: FROM THE TOP LEVEL DIRECTORY run: u-dales/tools/ecmwf_gather.sh <PATH_TO_CASE>"
 exit 1
fi

## Slurm account (project) charged for the job. Changing it only changes whose
## SBU budget pays; see `account` and https://hpc-usage.ecmwf.int.
ACCOUNT=spnlthee

## go to experiment directory
pushd "$1"
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
if [ -z "$DA_WORKDIR" ]; then
    echo "Output top-level directory DA_WORKDIR must be set inside $inputdir/config.sh"
    exit 1
fi;
if [ -z "$DA_TOOLSDIR" ]; then
    echo "Script directory DA_TOOLSDIR must be set inside $inputdir/config.sh"
    exit 1
fi;
if [ -z "$WALLTIME" ]; then
    echo "Wall clock time WALLTIME (HH:MM:SS, at most 48:00:00) must be set inside $inputdir/config.sh"
    exit 1
fi;
if [ -z "$MEM" ]; then
    echo "Memory requirement MEM must be set inside $inputdir/config.sh"
    echo "  export MEM=\"16G\"   # or the PBS spelling MEM=\"16gb\"; nf allows up to 128G"
    exit 1
fi;
## Accept both the Slurm spelling (16G) and the PBS one used on Imperial HPC
## (16gb), in any case, and normalise to what Slurm wants: digits + K/M/G/T.
if [[ "$MEM" =~ ^([0-9]+)([KkMmGgTt])[Bb]?$ ]]; then
    mem_slurm="${BASH_REMATCH[1]}${BASH_REMATCH[2]^^}"
else
    echo "MEM=$MEM is not a memory size; use e.g. MEM=16G or MEM=16gb inside $inputdir/config.sh"
    echo "A unit is required (K, M, G or T), so that the size is unambiguous."
    exit 1
fi
echo "MEM: $mem_slurm"

## Quality of service is fixed to nf and is deliberately not configurable:
## gathering is a single NCO process on shared nodes, which is exactly what nf
## is for (1 node, at most 128 cores / 128 GB, max 2 days). The QOS set in
## config.sh for the simulation is ignored here, so an np simulation does not
## drag the gather onto an exclusive 128-core node that would be billed in full
## (SBU = allocated cores x elapsed time) to run one process.
GATHER_QOS=nf
echo "QOS: $GATHER_QOS (fixed)"

## ECMWF's sbatch is a site wrapper (ecsbatch) that splits its arguments on
## whitespace, so it cannot submit from a directory whose path contains a space
## ("Unable to open file /lus/.../<first word>"). Catch that here rather than
## after the inputs have been copied.
if [[ "$DA_WORKDIR" =~ [[:space:]] ]] || [[ "$inputdir" =~ [[:space:]] ]]; then
    echo "Paths must not contain spaces: ECMWF's sbatch cannot submit from such a directory"
    echo "  DA_WORKDIR = $DA_WORKDIR"
    echo "  case dir   = $inputdir"
    exit 1
fi

## set the output directory
outdir="$DA_WORKDIR/$exp"

if [ ! -d "$outdir" ]; then
    echo "Output directory $outdir does not exist; run ecmwf_execute.sh first"
    exit 1
fi

echo "writing post-job.$exp.slurm"

## write post-job.exp.slurm file from inside the output directory: Slurm splits
## an #SBATCH directive on whitespace and does not honour quotes, so paths are
## kept out of the directives and the job inherits $outdir as its working
## directory (which is where slurm-%j.out then lands).
## gather_outputs.sh and nco_concatenate_field*.sh need ncks, ncpdq and ncrcat
## from NCO, plus ncdump from netCDF.
pushd "$outdir"
cat <<EOF > "post-job.$exp.slurm"
#!/bin/bash
#SBATCH --job-name=${exp}_gather
#SBATCH --account=${ACCOUNT}
#SBATCH --qos=${GATHER_QOS}
#SBATCH --time=${WALLTIME}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=${mem_slurm}
#SBATCH --output=slurm-%j.out
module load prgenv/intel intel/2021.4.0 netcdf4/4.10.0 nco/5.3.7
"$DA_TOOLSDIR/gather_outputs.sh" "$outdir" >> "$outdir/output.$exp.log" 2>&1
EOF

## submit post-job.exp.slurm file to queue
sbatch "post-job.$exp.slurm"

echo "post-job.$exp.slurm submitted."

popd
popd
