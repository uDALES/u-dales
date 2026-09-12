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

# Usage: FROM THE TOP LEVEL DIRECTORY run:
#   u-dales/tools/hpc_execute.sh <PATH_TO_CASE>
#
# One script for CX3 and HX1, CPU and GPU. The cluster is detected; the rest
# comes from config.sh in the case directory:
#
#   NNODE, NCPU, WALLTIME, MEM   required for every run
#   NGPU                         GPUs per node; setting it makes this a GPU run
#   QUEUE                        optional, adds "#PBS -q <QUEUE>"
#   PLACE                        optional, adds "#PBS -l place=<PLACE>", e.g. excl
#   GPU_TYPE                     GPU card for the select line, default A100;
#                                set to "" to omit the constraint
#
# Environment overrides:
#   UDALES_SYSTEM   cx3 | hx1    force the cluster instead of detecting it
#   UDALES_TARGET   cpu | gpu    force the target instead of deriving it from NGPU

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
    echo "Output top-level directory DA_WORKDIR must be set inside $inputdir/config.sh"
    exit 1
fi;
if [ -z $DA_BUILD ]; then
    echo "Executable DA_BUILD must be set inside $inputdir/config.sh"
    exit 1
fi;
if [ -z $DA_TOOLSDIR ]; then
    echo "Script directory DA_TOOLSDIR must be set inside $inputdir/config.sh"
    exit 1
fi;
if [ -z $NNODE ]; then
    echo "Number of nodes NNODE must be set inside $inputdir/config.sh"
    echo "Product of NNODE and NCPU set in $inputdir/config.sh must be equal to the product of nprocx and nprocy set in $inputdir/namoptions.$exp"
    exit 1
fi;
if [ -z $NCPU ]; then
    echo "Number of CPU cores on each node NCPU must be set inside $inputdir/config.sh"
    echo "Product of NNODE and NCPU set in $inputdir/config.sh must be equal to the product of nprocx and nprocy set in $inputdir/namoptions.$exp"
    exit 1
fi;
if [ -z $WALLTIME ]; then
    echo "Wall clock time WALLTIME must be set inside $inputdir/config.sh"
    exit 1
fi;
if [ -z $MEM ]; then
    echo "Memory requirement MEM must be set inside $inputdir/config.sh"
    exit 1
fi;

## ---------------------------------------------------------------------------
## Which cluster, and CPU or GPU
##
## The cluster is worked out from the login node so the same script serves CX3
## and HX1 with no editing; set UDALES_SYSTEM=<cx3|hx1> to override it. The run
## target follows NGPU - set NGPU in config.sh for a GPU run, leave it out for a
## CPU one - and UDALES_TARGET=<cpu|gpu> forces it either way.
## ---------------------------------------------------------------------------
if [ -z "${UDALES_SYSTEM:-}" ]; then
    # Three tells, tried in turn.
    # 1. The hostname: HX1 nodes are hx1-... (login: hx1-c12-login-1); CX3
    #    compute nodes are cx3-... but its login nodes are login-a, login-b,
    #    login-ai, login-bi.
    # 2. The module tree: HX1 keeps EasyBuild under /gpfs/easybuild/prod, CX3
    #    under /sw-eb.
    # 3. The PBS server from /etc/pbs.conf: pbs-6 serves HX1; Imperial's other
    #    server is CX3's.
    pbs_server=$(sed -n 's/^PBS_SERVER=//p' /etc/pbs.conf 2>/dev/null)
    case "$(hostname -s)" in
        hx1-*)         UDALES_SYSTEM=hx1; cluster_tell="hostname $(hostname -s)" ;;
        cx3-*|login-*) UDALES_SYSTEM=cx3; cluster_tell="hostname $(hostname -s)" ;;
    esac
    if [ -z "${UDALES_SYSTEM:-}" ]; then
        if [ -d /gpfs/easybuild/prod ]; then
            UDALES_SYSTEM=hx1; cluster_tell="module tree /gpfs/easybuild/prod"
        elif [ -d /sw-eb ]; then
            UDALES_SYSTEM=cx3; cluster_tell="module tree /sw-eb"
        fi
    fi
    if [ -z "${UDALES_SYSTEM:-}" ]; then
        case "$pbs_server" in
            pbs-6.*)        UDALES_SYSTEM=hx1; cluster_tell="PBS server $pbs_server" ;;
            *.hpc.ic.ac.uk) UDALES_SYSTEM=cx3; cluster_tell="PBS server $pbs_server" ;;
        esac
    fi
    if [ -z "${UDALES_SYSTEM:-}" ]; then
        echo "Could not tell which cluster this is: hostname $(hostname -s)," \
             "no /gpfs/easybuild/prod or /sw-eb, PBS server ${pbs_server:-none}."
        echo "Set UDALES_SYSTEM=cx3 or UDALES_SYSTEM=hx1 and run again."
        exit 1
    fi
else
    cluster_tell="UDALES_SYSTEM set by hand"
fi

if [ -z "${UDALES_TARGET:-}" ]; then
    if [ -n "${NGPU:-}" ] && [ "${NGPU:-0}" -gt 0 ]; then
        UDALES_TARGET=gpu
    else
        UDALES_TARGET=cpu
    fi
fi

if [ "$UDALES_TARGET" = "gpu" ]; then
    if [ -z "${NGPU:-}" ] || [ "${NGPU:-0}" -lt 1 ]; then
        echo "A GPU run needs the number of GPUs per node, NGPU, set inside $inputdir/config.sh"
        exit 1
    fi
    if [ ! -f "$DA_TOOLSDIR/bind.sh" ]; then
        echo "bind.sh not found in $DA_TOOLSDIR; it is needed to give each rank its own GPU."
        exit 1
    fi
fi

echo "cluster: $UDALES_SYSTEM ($cluster_tell)"
echo "target:  $UDALES_TARGET"

## The executable has to match the target. build_executable.sh writes CPU builds
## under build/cpu and GPU builds under build/gpu, so the target is readable from
## the path - and catching a mismatch here saves a job that would queue, start,
## and only then misbehave on the node. A build directory outside that layout
## (UDALES_BUILD_DIR, or a copied executable) says nothing either way, so it
## warns rather than refusing.
case "$DA_BUILD" in
    */build/cpu/*) build_kind=cpu ;;
    */build/gpu/*) build_kind=gpu ;;
    *)             build_kind=unknown ;;
esac

if [ "$build_kind" = "unknown" ]; then
    echo "warning: cannot tell from DA_BUILD whether this is a CPU or a GPU build:"
    echo "  $DA_BUILD"
    echo "warning: expected it under u-dales/build/cpu or u-dales/build/gpu. Continuing."
elif [ "$build_kind" != "$UDALES_TARGET" ]; then
    echo "This is a $UDALES_TARGET run, but DA_BUILD is a $build_kind build:"
    echo "  $DA_BUILD"
    if [ "$UDALES_TARGET" = "gpu" ]; then
        echo "Set DA_BUILD to u-dales/build/gpu/<release|debug>/u-dales in $inputdir/config.sh,"
        echo "or unset NGPU there to run on CPUs."
    else
        echo "Set DA_BUILD to u-dales/build/cpu/<release|debug>/u-dales in $inputdir/config.sh,"
        echo "or set NGPU there to run on GPUs."
    fi
    exit 1
fi

## The module sets below have to stay in step with the matching block of
## tools/build_executable.sh: the solver must run against the libraries it was
## built against. Single quotes keep $EBROOTNVHPC and $PATH for the compute
## node rather than expanding them here.
case "$UDALES_SYSTEM:$UDALES_TARGET" in
    cx3:cpu)
        # The self-consistent 2021a stack from the "icl" block of
        # build_executable.sh; tools/prod makes it resolve from a clean shell.
        job_modules='module load tools/prod
module load intel/2021a netCDF/4.8.0-iimpi-2021a netCDF-Fortran/4.5.3-iimpi-2021a FFTW/3.3.9-intel-2021a CMake/3.20.1-GCCcore-10.3.0 git/2.32.0-GCCcore-10.3.0-nodocs'
        ;;
    hx1:cpu)
        # The runtime part of the "hx1" block of build_executable.sh.
        job_modules='module load intel/2023a netCDF/4.9.2-iimpi-2023a netCDF-Fortran/4.6.1-iimpi-2023a FFTW/3.3.10-intel-compilers-2023.1.0'
        ;;
    hx1:gpu)
        # Mirrors "gpuhx1". netCDF and FFTW are reached through the executable's
        # RPATH, so NVHPC is the only module needed - but its mpirun lives under
        # comm_libs, which the module does not put on PATH.
        job_modules='module load NVHPC/23.7-CUDA-12.2.0
export PATH="$EBROOTNVHPC/Linux_x86_64/23.7/comm_libs/mpi/bin:$PATH"'
        ;;
    cx3:gpu)
        echo "There is no CX3 GPU target in tools/build_executable.sh."
        echo "Build with 'gpuhx1' and run on HX1, or add a CX3 GPU block there first."
        exit 1
        ;;
    *)
        echo "Unsupported combination: $UDALES_SYSTEM / $UDALES_TARGET"
        exit 1
        ;;
esac

## PBS resources and launcher.
##
## mpiprocs is a per-chunk count, so it is the per-node figure and not the
## total: PBS multiplies it by the number of chunks itself.
if [ "$UDALES_TARGET" = "gpu" ]; then
    NP=$(( NNODE * NGPU ))
    pbs_select="select=${NNODE}:ncpus=${NCPU}:mpiprocs=${NGPU}:ngpus=${NGPU}:mem=${MEM}"
    # HX1's a100 queue is chosen through the resource request rather than -q, and
    # every GPU example in the RCS user guide names the card. This branch is only
    # reachable on HX1, where A100 is the only GPU. Set GPU_TYPE to ask for
    # something else, or to the empty string to leave the constraint off.
    gpu_type="${GPU_TYPE-A100}"
    if [ -n "$gpu_type" ]; then
        pbs_select="${pbs_select}:gpu_type=${gpu_type}"
    fi
    # One rank per GPU; bind.sh pins each to its own device.
    launch="mpirun -n ${NP} ${DA_TOOLSDIR}/bind.sh ${DA_BUILD}"
else
    NP=$(( NNODE * NCPU ))
    pbs_select="select=${NNODE}:ncpus=${NCPU}:mpiprocs=${NCPU}:mem=${MEM}"
    launch="mpirun -v6 -n ${NP} ${DA_BUILD}"
fi

pbs_directives="#PBS -l walltime=${WALLTIME}
#PBS -l ${pbs_select}"
if [ -n "${QUEUE:-}" ]; then
    pbs_directives="${pbs_directives}
#PBS -q ${QUEUE}"
fi
# Optional placement, e.g. PLACE=excl for a node of your own. On the a100
# queue that holds a whole 4-GPU node whatever NGPU is, so use it deliberately.
if [ -n "${PLACE:-}" ]; then
    pbs_directives="${pbs_directives}
#PBS -l place=${PLACE}"
fi

## set the output directory
outdir=$DA_WORKDIR/$exp

echo "writing job.$exp with $NP MPI ranks."

## write new job.exp file for HPC
cat <<EOF > job.$exp
#!/bin/bash
${pbs_directives}
${job_modules}
EOF

## Report how long this job waited in the queue, from PBS's own timestamps.
## Quoted heredoc: nothing below is expanded at submit time.
cat <<'EOF' >> job.$exp
queue_wait_line() {
    local info qt st w
    if [ -n "${PBS_JOBID:-}" ] && command -v qstat >/dev/null 2>&1; then
        info=$(qstat -f "$PBS_JOBID" 2>/dev/null)
        qt=$(printf '%s\n' "$info" | awk -F' = ' '/^[[:space:]]*qtime = /{print $2}')
        st=$(printf '%s\n' "$info" | awk -F' = ' '/^[[:space:]]*stime = /{print $2}')
        if [ -n "$qt" ] && [ -n "$st" ]; then
            w=$(( $(date -d "$st" +%s) - $(date -d "$qt" +%s) ))
            printf 'PBS job %s on %s: queued %s, started %s, queue wait %dm %02ds\n' \
                "$PBS_JOBID" "$(hostname -s)" "$qt" "$st" $((w / 60)) $((w % 60))
            return
        fi
    fi
    echo "PBS queue wait: unknown (not under PBS, or qstat unavailable)"
}
EOF

## The queue-wait line goes into output.exp.log ahead of this run's solver
## output (master's convention, shared with local_execute.sh and hpc_gather.sh).
cat <<EOF >> job.$exp
mkdir -p $outdir
cp -r $inputdir/* $outdir
pushd $outdir
queue_wait_line >> $outdir/output.$exp.log
echo "cluster: $UDALES_SYSTEM ($cluster_tell), target: $UDALES_TARGET" >> $outdir/output.$exp.log
${launch} $outdir/namoptions.$exp >> $outdir/output.$exp.log 2>&1
EOF

## submit job.exp file to queue
qsub job.$exp

echo "job.$exp submitted."
