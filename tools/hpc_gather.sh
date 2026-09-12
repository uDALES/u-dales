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
 echo "usage: FROM THE TOP LEVEL DIRECTORY run: u-dales/tools/hpc_gather.sh <PATH_TO_CASE>"
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
if [ -z $DA_TOOLSDIR ]; then
    echo "Script directory DA_TOOLSDIR must be set inside $inputdir/config.sh"
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
## Which cluster. Detected from the login node so the same script serves CX3
## and HX1; set UDALES_SYSTEM=<cx3|hx1> to override it.
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
echo "cluster: $UDALES_SYSTEM ($cluster_tell)"

## gather_outputs.sh and nco_concatenate_field*.sh need ncks, ncpdq and ncrcat
## from NCO, plus ncdump from netCDF.
case "$UDALES_SYSTEM" in
    cx3)
        gather_modules='module load tools/prod
module load NCO/5.2.9-foss-2024a'
        ;;
    hx1)
        # HX1 ships no NCO at all, so it comes from a local build (see
        # $HOME/nco-5.2.9). The netCDF module is still needed at run time: its
        # libnetcdf.so finds HDF5 through LD_LIBRARY_PATH, and it supplies ncdump.
        nco_bin="${NCO_BIN:-$HOME/nco-5.2.9/nco/bin}"
        if [ ! -x "$nco_bin/ncks" ]; then
            echo "No NCO found at $nco_bin"
            echo "HX1 has no NCO module, so it has to be built locally:"
            echo "  bash \$HOME/nco-5.2.9/s2_configure_build_install.sh"
            echo "Or set NCO_BIN to a directory containing ncks, ncpdq and ncrcat."
            exit 1
        fi
        gather_modules="module load netCDF/4.9.2-gompi-2023a UDUNITS/2.2.28-GCCcore-12.3.0 GSL/2.7-GCC-12.3.0
export PATH=$nco_bin:\$PATH"
        ;;
esac

## set the output directory
outdir=$DA_WORKDIR/$exp

echo "writing post-job.$exp."

## write post-job.exp file for HPC
cat <<EOF > post-job.$exp
#!/bin/bash
#PBS -l walltime=${WALLTIME}
#PBS -l select=1:ncpus=1:mem=${MEM}
${gather_modules}
EOF

## Report how long this job waited in the queue, from PBS's own timestamps.
## Quoted heredoc: nothing below is expanded at submit time.
cat <<'EOF' >> post-job.$exp
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

cat <<EOF >> post-job.$exp
queue_wait_line >> $outdir/output.$exp.log
echo "cluster: $UDALES_SYSTEM ($cluster_tell)" >> $outdir/output.$exp.log
$DA_TOOLSDIR/gather_outputs.sh $outdir >> $outdir/output.$exp.log 2>&1
EOF

## submit post-job.exp file to queue
qsub post-job.$exp

echo "post-job.$exp submitted."
