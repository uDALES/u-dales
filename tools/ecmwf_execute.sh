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

## Submit a uDALES run on the ECMWF Atos HPC2020 (Slurm).
## Build first with: tools/build_executable.sh ecmwf release

set -e

if (( $# < 1 ))
then
 echo "The experiment directory must be set."
 echo "usage: FROM THE TOP LEVEL DIRECTORY run: u-dales/tools/ecmwf_execute.sh <PATH_TO_CASE>"
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
if [ -z "$DA_BUILD" ]; then
    echo "Executable DA_BUILD must be set inside $inputdir/config.sh"
    exit 1
fi;
if [ ! -x "$DA_BUILD" ]; then
    echo "Executable DA_BUILD=$DA_BUILD does not exist or is not executable"
    echo "Build it first with: tools/build_executable.sh ecmwf release"
    exit 1
fi
if [ -z "$DA_TOOLSDIR" ]; then
    echo "Script directory DA_TOOLSDIR must be set inside $inputdir/config.sh"
    exit 1
fi;
if [ -z $NNODE ]; then
    echo "Number of nodes NNODE must be set inside $inputdir/config.sh"
    echo "Product of NNODE and NCPU set in $inputdir/config.sh must be equal to the product of nprocx and nprocy set in $inputdir/namoptions.$exp"
    exit 1
fi;
if [ -z $NCPU ]; then
    echo "Number of MPI ranks on each node NCPU must be set inside $inputdir/config.sh"
    echo "Product of NNODE and NCPU set in $inputdir/config.sh must be equal to the product of nprocx and nprocy set in $inputdir/namoptions.$exp"
    exit 1
fi;
if [ -z $WALLTIME ]; then
    echo "Wall clock time WALLTIME (HH:MM:SS, at most 48:00:00) must be set inside $inputdir/config.sh"
    exit 1
fi;

## Compute nodes have 128 physical cores (256 hyperthreads, which are not used
## here) and 240 GB of memory.
if (( NCPU > 128 )); then
    echo "NCPU=$NCPU is more than the 128 physical cores of an ECMWF compute node"
    exit 1
fi

## Quality of service. Must be set as QOS in config.sh: there is no default,
## because the queues differ in how nodes and memory are allocated and the
## wrong one silently changes what the job gets.
##   np  parallel jobs using more than half a node; nodes are exclusive and each
##       comes with all its 240 GB, so no memory request is needed. Max 2 days.
##   nf  shared nodes for serial and small parallel jobs: at most 128 cores on
##       1 node and 128 GB (sacctmgr: MaxTRES cpu=128,mem=128G,node=1). A job
##       that does not ask for memory gets the shared pool default of 8 GB, so
##       MEM is required below. Max 2 days.
##   ni  interactive-style shared jobs: 1 job, at most 32 cores / 32 GB. Max 7 days.
##   ef  ECS users only; ng / dg GPU and GPU-debug queues (not used by uDALES).
if [ -z "${QOS:-}" ]; then
    echo "Quality of service QOS must be set inside $inputdir/config.sh"
    echo "  export QOS=\"np\"   # parallel runs on exclusive nodes, 240 GB per node (usual choice)"
    echo "  export QOS=\"nf\"   # shared nodes for small test runs: 1 node, max 128 cores / 128 GB"
    echo "See https://confluence.ecmwf.int/display/UDOC/HPC2020%3A+Batch+system for the full list."
    exit 1
fi
echo "QOS: $QOS"

## Memory. Only np gives a node exclusively with all of its 240 GB; every other
## queue shares nodes, where a job that does not ask for memory is given the
## shared pool default of 8 GB (scontrol: DefMemPerNode=8000). So --mem is
## emitted for those queues only, and MEM must be set for them.
mem_directive=""
if [ "$QOS" != "np" ]; then
    if [ -z "${MEM:-}" ]; then
        echo "Memory requirement MEM must be set inside $inputdir/config.sh for QOS=$QOS"
        echo "  export MEM=\"16G\"   # per node, or the PBS spelling MEM=\"16gb\"; nf allows up to 128G"
        echo "QOS=np needs no MEM: an exclusive node comes with all its 240 GB."
        exit 1
    fi
    ## Accept both the Slurm spelling (16G) and the PBS one used on Imperial HPC
    ## (16gb), in any case, and normalise to what Slurm wants: digits + K/M/G/T.
    if [[ "$MEM" =~ ^([0-9]+)([KkMmGgTt])[Bb]?$ ]]; then
        mem_slurm="${BASH_REMATCH[1]}${BASH_REMATCH[2]^^}"
    else
        echo "MEM=$MEM is not a memory size; use e.g. MEM=16G or MEM=16gb inside $inputdir/config.sh"
        echo "A unit is required (K, M, G or T), so that the size is unambiguous."
        exit 1
    fi
    mem_directive="#SBATCH --mem=${mem_slurm}"
    echo "MEM: $mem_slurm (QOS=$QOS shares nodes)"
fi

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

## Output can reach several TB: only SCRATCH (50 TB quota) is sized for it.
## HOME, PERM and HPCPERM are small and not meant for parallel I/O.
## NOTE: SCRATCH files are deleted 30 days after last access; move results you
## want to keep (e.g. to ECFS).
## Check $SCRATCH before calling realpath on it: realpath -m "" fails, which
## under set -e would abort with a bare "realpath: '': No such file or
## directory" instead of the message below.
if [ -z "${SCRATCH:-}" ]; then
    echo "\$SCRATCH is not set, so DA_WORKDIR cannot be checked against it."
    echo "On ECMWF it is set for you; if config.sh overrides it, remove that."
    exit 1
fi
scratch_real=$(realpath -m "$SCRATCH")
workdir_real=$(realpath -m "$DA_WORKDIR")
if [[ "$workdir_real/" != "$scratch_real/"* ]]; then
    echo "DA_WORKDIR=$DA_WORKDIR is not under \$SCRATCH ($SCRATCH)."
    echo "Set DA_WORKDIR inside $inputdir/config.sh to a directory under \$SCRATCH."
    exit 1
fi

## set the output directory
outdir="$DA_WORKDIR/$exp"

## copy files to execution and output directory; the executable is copied too,
## so rebuilding while the job is queued does not change the run.
mkdir -p "$outdir"
cp -r -P "$inputdir"/* "$outdir"
cp "$DA_BUILD" "$outdir/u-dales"
pushd "$outdir"

echo "writing job.$exp.slurm"

## write new job.exp.slurm file
## Modules are the runtime the executable was built against; keep in step with
## the ecmwf block of tools/build_executable.sh.
## Paths are kept out of the #SBATCH directives: Slurm splits a directive on
## whitespace and does not honour quotes, so an absolute --chdir/--output would
## break for a DA_WORKDIR containing spaces. The job is submitted from $outdir,
## which is Slurm's default working directory, so slurm-%j.out lands there too.
sbatch_directives="#SBATCH --job-name=${exp}
#SBATCH --account=${ACCOUNT}
#SBATCH --qos=${QOS}
#SBATCH --time=${WALLTIME}
#SBATCH --nodes=${NNODE}
#SBATCH --ntasks-per-node=${NCPU}
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --output=slurm-%j.out"
if [ -n "$mem_directive" ]; then
    sbatch_directives="${sbatch_directives}
${mem_directive}"
fi

cat <<EOF > "job.$exp.slurm"
#!/bin/bash
${sbatch_directives}
module load prgenv/intel intel/2021.4.0 intel-mpi/2021.4.0 netcdf4/4.10.0 fftw/3.3.10
export OMP_NUM_THREADS=1
srun ./u-dales "$outdir/namoptions.$exp" >> "$outdir/output.$exp.log" 2>&1
EOF

## submit job.exp.slurm file to queue
sbatch "job.$exp.slurm"

echo "job.$exp.slurm submitted."

popd
