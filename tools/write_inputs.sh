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

# script for setting key namoptions and writing *.inp.* files
# variables not specified in namoptions can be set here. It updates these in
# write_inputs.m or write_inputs.py and then runs the selected preprocessing
# route to produce the required text files.

# tg3315 20/07/2017, modified by SO 06/02/20, modified by DM 27/08/24

# NOTES
# (1) if no forcing found in namoptions then applies the initial velocities and uses the pressure terms

# Assuming running from top-level project directory.
# u-dales/tools/write_inputs.sh <-m|-p> <PATH_TO_CASE> [start] [system]
#
#   start   l  (default) run on the current node in the background: the script
#              returns at once and the output goes to write_inputs.<expnr>.log
#           f  run on the current node in the foreground: the script blocks,
#              the output goes to the terminal and to the log, and the exit
#              code of the preprocessing is returned
#           c  submit to a compute node through the scheduler of <system>
#              (icl: PBS, ecmwf: Slurm)
#   system  common (default), icl or ecmwf
#
# Scheduler settings (environment or the case's config.sh):
#   PREPROC_WALLTIME  default 24:00:00
#   PREPROC_MEM       default 128gb (number followed by lowercase gb)
#   PREPROC_ACCOUNT   ecmwf only: Slurm account, default spnlthee
#   PREPROC_QOS       ecmwf only: Slurm QOS, default nf (shared node, at most
#                     128 cores / 128 GB)

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
	echo "usage: FROM THE TOP LEVEL DIRECTORY run: u-dales/tools/write_inputs.sh <-m|-p> <PATH_TO_CASE> [start] [system]"
	echo "   -m: run MATLAB preprocessing route (common and icl only)"
	echo "   -p: run Python preprocessing route"
	echo "   start (optional): default l"
	echo "      l: current node, in the background; returns at once, output in write_inputs.<expnr>.log"
	echo "      f: current node, in the foreground; blocks, output on the terminal and in the log, exit code returned"
	echo "      c: compute node through the scheduler (icl: PBS, ecmwf: Slurm)"
	echo "   system (optional): common (default), icl or ecmwf"
}

default_nompthreads=8

extract_nompthreads() {
	local namoptions_file=$1
	local occurrence_count
	local nompthreads_line
	local nompthreads_value

	# grep -c returns 1 when nothing matches; that is a count of 0, not an error
	occurrence_count=$(grep -ic "nompthreads" "$namoptions_file" || true)

	if [ "$occurrence_count" -eq 0 ]; then
		return 0
	fi
	if [ "$occurrence_count" -gt 1 ]; then
		echo "MULTIPLE"
		return 0
	fi

	nompthreads_line=$(grep -i "nompthreads" "$namoptions_file")
	if [[ "$nompthreads_line" != *"="* ]]; then
		echo "INVALID"
		return 0
	fi

	nompthreads_value=${nompthreads_line#*=}
	nompthreads_value=${nompthreads_value%%!*}
	nompthreads_value=${nompthreads_value%%,*}
	nompthreads_value=${nompthreads_value%%/*}
	nompthreads_value="${nompthreads_value#"${nompthreads_value%%[![:space:]]*}"}"
	nompthreads_value="${nompthreads_value%"${nompthreads_value##*[![:space:]]}"}"

	if [[ "$nompthreads_value" =~ ^[0-9]+$ ]] && [ "$nompthreads_value" -gt 0 ]; then
		printf "%d\n" "$nompthreads_value"
	else
		echo "INVALID"
	fi
}

source_view3d_config() {
	local config_file=$1
	local status

	set -a
	if source "$config_file"; then
		set +a
		return 0
	else
		status=$?
		set +a
		echo "Failed to source View3D config: $config_file" >&2
		echo "Check the file syntax and required environment variables." >&2
		return "$status"
	fi
}

if (( $# < 2 ))
then
	echo "The preprocessing route and path to case/experiment folder must be set."
	usage
	echo "... execution terminated"
    exit 1
fi

route_arg=$1
case_path=$2
start=${3:-"l"}     # l: background on this node, f: foreground on this node, c: compute node through the scheduler
system=${4:-"common"}

case "$route_arg" in
	-m)
		route="matlab"
		;;
	-p)
		route="python"
		;;
	*)
		echo "Unrecognised preprocessing route: $route_arg"
		usage
		echo "... execution terminated"
		exit 1
		;;
esac

case "$start" in
	l|f|c) ;;
	*)
		echo "Unrecognised start mode: $start (expected l, f or c)"
		usage
		echo "... execution terminated"
		exit 1
		;;
esac

case "$system" in
	common|icl|ecmwf) ;;
	*)
		echo "Unrecognised system: $system (expected common, icl or ecmwf)"
		usage
		echo "... execution terminated"
		exit 1
		;;
esac

if [ "$start" == "c" ] && [ "$system" == "common" ]; then
	echo "start mode c (compute node) needs a scheduler: pass the system as well, e.g."
	echo "  u-dales/tools/write_inputs.sh $route_arg $case_path c icl"
	echo "  u-dales/tools/write_inputs.sh $route_arg $case_path c ecmwf"
	echo "... execution terminated"
	exit 1
fi

if [ "$route" == "matlab" ] && [ "$system" == "ecmwf" ]; then
	echo "The MATLAB preprocessing route is not available on ecmwf; use the Python route (-p)."
	echo "... execution terminated"
	exit 1
fi

# go to experiment directory
pushd "$case_path"
inputdir=$(pwd)

## set experiment number via path
iexpnr="${inputdir: -3}"

## read in additional variables
if [ -f config.sh ]; then
  source config.sh
fi

export DA_TOOLSDIR="${DA_TOOLSDIR:-$script_dir}"
export DA_EXPDIR="${DA_EXPDIR:-$(cd "$inputdir/.." && pwd)}"

namoptions_file="$inputdir/namoptions.$iexpnr"
if [ ! -f "$namoptions_file" ]; then
	echo "Namelist file not found: $namoptions_file"
	exit 1
fi

nompthreads=$(extract_nompthreads "$namoptions_file")
case "$nompthreads" in
	MULTIPLE)
		echo "Multiple nompthreads occurrences found in $namoptions_file"
		echo "Set nompthreads only once."
		exit 1
		;;
	INVALID)
		echo "Invalid nompthreads value in $namoptions_file"
		echo "nompthreads must be a positive integer."
		exit 1
		;;
esac
export PREPROC_NCPU="${nompthreads:-$default_nompthreads}"
export NUMBA_NUM_THREADS="${NUMBA_NUM_THREADS:-$PREPROC_NCPU}"
export PREPROC_WALLTIME="${PREPROC_WALLTIME:-24:00:00}"
export PREPROC_MEM="${PREPROC_MEM:-128gb}"
## ecmwf scheduler settings; see tools/ecmwf_execute.sh for the queues.
PREPROC_ACCOUNT="${PREPROC_ACCOUNT:-spnlthee}"
PREPROC_QOS="${PREPROC_QOS:-nf}"

if [[ ! "$PREPROC_MEM" =~ ^[0-9]+gb$ ]]; then
	echo "Memory requirement PREPROC_MEM must be set like 128gb"
	exit 1
fi

popd

export VIEW3D_CONFIG="${VIEW3D_CONFIG:-$DA_TOOLSDIR/view3d_config.sh}"
if [ ! -f "$VIEW3D_CONFIG" ]; then
	echo "View3D config file not found: $VIEW3D_CONFIG" >&2
	exit 1
fi
source_view3d_config "$VIEW3D_CONFIG" || exit 1

python_exe="$DA_TOOLSDIR/python/.venv/bin/python"
if [ "$route" == "python" ]; then
	if [ ! -x "$python_exe" ]; then
		echo "Python virtual environment not found or not executable: $python_exe"
		echo "Set it up with: bash $DA_TOOLSDIR/python/setup_venv.sh <common|icl|ecmwf>"
		echo "... execution terminated"
		exit 1
	fi
	# Keep f2py/Fortran status messages ordered with Python output in redirected logs.
	export PYTHONUNBUFFERED=1
	export GFORTRAN_UNBUFFERED_PRECONNECTED=1
fi

log_file="$inputdir/write_inputs.$iexpnr.log"

if command -v stdbuf >/dev/null 2>&1; then
	stdbuf_cmd="stdbuf -oL -eL"
else
	stdbuf_cmd=""
fi

## Body shared by the PBS and the Slurm job scripts: environment, View3D
## configuration and the preprocessing command. Written after the scheduler
## header and the module loads, which differ per system.
write_job_body() {
cat <<EOF
cd "$DA_TOOLSDIR"

export DA_TOOLSDIR="$DA_TOOLSDIR"
export DA_EXPDIR="$DA_EXPDIR"
export MATLAB_USE_USERWORK=0
export OMP_NUM_THREADS=$PREPROC_NCPU
export PYTHONUNBUFFERED=1
export GFORTRAN_UNBUFFERED_PRECONNECTED=1
export PREPROC_NCPU="$PREPROC_NCPU"
export NUMBA_NUM_THREADS="$NUMBA_NUM_THREADS"
export PREPROC_MEM="$PREPROC_MEM"
export VIEW3D_CONFIG="$VIEW3D_CONFIG"
source_view3d_config() {
	local config_file=\$1
	local status

	set -a
	if source "\$config_file"; then
		set +a
		return 0
	else
		status=\$?
		set +a
		echo "Failed to source View3D config: \$config_file" >&2
		echo "Check the file syntax and required environment variables." >&2
		return "\$status"
	fi
}
if [ ! -f "$VIEW3D_CONFIG" ]; then
	echo "View3D config file not found: $VIEW3D_CONFIG" >&2
	exit 1
fi
source_view3d_config "$VIEW3D_CONFIG" || exit 1

EOF

	if [ "$route" == "matlab" ]; then
		cat <<EOF
module load MATLAB/2024b
nohup matlab -nodesktop -noFigureWindows -nosplash -nodisplay -r "expnr=$iexpnr; write_inputs; quit" > "$log_file" 2>&1 < /dev/null

EOF
	else
		cat <<EOF
if command -v stdbuf >/dev/null 2>&1; then
	stdbuf_cmd="stdbuf -oL -eL"
else
	stdbuf_cmd=""
fi
\$stdbuf_cmd "$python_exe" -u "$DA_TOOLSDIR/write_inputs.py" "$inputdir" > "$log_file" 2>&1 < /dev/null

EOF
	fi
}

if [ "$start" == "c" ]; then

	cd "$inputdir"

	if [ "$system" == "icl" ]; then
		###### Imperial HPC: PBS job script
		job_file="pre-job.$iexpnr"
		{
cat <<EOF
#!/bin/bash
#PBS -l walltime=$PREPROC_WALLTIME
#PBS -l select=1:ncpus=$PREPROC_NCPU:mem=${PREPROC_MEM}

module load tools/prod
module load GCC/14.2.0
module load Python/3.13.1-GCCcore-14.2.0

EOF
		write_job_body
		} > "$job_file"

		## submit job file to queue
		if ! qsub "$job_file"; then
			echo "[ERROR] qsub failed; job script kept at $job_file" >&2
			exit 1
		fi
		echo "$job_file submitted."

	else
		###### ECMWF Atos HPC2020: Slurm job script
		## ECMWF's sbatch is a site wrapper (ecsbatch) that splits its arguments
		## on whitespace, so it cannot submit from a directory whose path
		## contains a space. Paths are also kept out of the #SBATCH directives.
		if [[ "$inputdir" =~ [[:space:]] ]] || [[ "$DA_TOOLSDIR" =~ [[:space:]] ]]; then
			echo "Paths must not contain spaces: ECMWF's sbatch cannot submit from such a directory"
			echo "  case dir    = $inputdir"
			echo "  DA_TOOLSDIR = $DA_TOOLSDIR"
			exit 1
		fi
		## Only the np queue gives a node exclusively (no --mem needed); every
		## other queue shares nodes and a job that does not ask for memory gets
		## the shared pool default of 8 GB. Preprocessing is one process with
		## PREPROC_NCPU threads, so nf is the usual choice.
		mem_directive=""
		if [ "$PREPROC_QOS" != "np" ]; then
			mem_directive="#SBATCH --mem=${PREPROC_MEM%gb}G"
		fi
		job_file="pre-job.$iexpnr.slurm"
		{
cat <<EOF
#!/bin/bash
#SBATCH --job-name=pre-$iexpnr
#SBATCH --account=$PREPROC_ACCOUNT
#SBATCH --qos=$PREPROC_QOS
#SBATCH --time=$PREPROC_WALLTIME
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$PREPROC_NCPU
#SBATCH --hint=nomultithread
#SBATCH --output=pre-job.$iexpnr.out
${mem_directive}
set -eo pipefail

## Runtime of the preprocessing tools; keep in step with the ecmwf block of
## tools/build_preprocessing.sh (the f2py modules and View3D link libgfortran
## and libgomp from this gcc).
module load prgenv/gnu gcc/14.2.0 python3/3.13.16-01

EOF
		write_job_body
		} > "$job_file"

		## submit job file to queue
		if ! sbatch "$job_file"; then
			echo "[ERROR] sbatch failed; job script kept at $job_file" >&2
			exit 1
		fi
		echo "$job_file submitted (account $PREPROC_ACCOUNT, qos $PREPROC_QOS, $PREPROC_NCPU cpus, ${PREPROC_MEM}, $PREPROC_WALLTIME)."
		echo "Scheduler output: $inputdir/pre-job.$iexpnr.out; preprocessing output: $log_file"
	fi

elif [ "$start" == "f" ]; then
	###### RUN SCRIPT in the foreground on this node; exit code is the preprocessing's
	cd "$DA_TOOLSDIR"
	set +e
	if [ "$route" == "matlab" ]; then
		matlab -nodesktop -noFigureWindows -nosplash -nodisplay -r "expnr=$iexpnr; write_inputs; quit" < /dev/null 2>&1 | tee "$log_file"
	else
		$stdbuf_cmd "$python_exe" -u "$DA_TOOLSDIR/write_inputs.py" "$inputdir" < /dev/null 2>&1 | tee "$log_file"
	fi
	status=${PIPESTATUS[0]}
	set -e
	if [ "$status" -ne 0 ]; then
		echo "Preprocessing of $inputdir failed with exit code $status (log: $log_file)" >&2
	fi
	exit "$status"

else
	###### RUN SCRIPT in the background on this node
	if [ "$system" == "ecmwf" ]; then
		echo "Note: this runs on the current node. ECMWF login nodes give a session very few cores and"
		echo "      end long processes; for anything but a small case use start mode c (Slurm job)."
	fi
	cd "$DA_TOOLSDIR"
	if [ "$route" == "matlab" ]; then
		nohup matlab -nodesktop -noFigureWindows -nosplash -nodisplay -r "expnr=$iexpnr; write_inputs; quit" > "$log_file" 2>&1 < /dev/null &
	else
		nohup $stdbuf_cmd "$python_exe" -u "$DA_TOOLSDIR/write_inputs.py" "$inputdir" > "$log_file" 2>&1 < /dev/null &
	fi
	cd "$DA_EXPDIR"
	cd ..
fi
