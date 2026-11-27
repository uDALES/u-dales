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

# Copyright (C) 2016-2026 the uDALES Team.

# Thermal-comfort indices from a gathered uDALES case in one command. Runs
# tools/python/udcomf/ud_run_thermal_comfort.py so that it survives a closed
# terminal: in the background on this node, in the foreground, or as a
# scheduler job (PBS on icl, Slurm on ecmwf), like tools/write_inputs.sh.
#
# Assuming running from top-level project directory.
# u-dales/tools/ud_run_thermal_comfort.sh <PATH_TO_CASE> [start] [system] -- --heights H [H ...] [wrapper options]
#
#   start   l  (default) run on the current node in the background: the script
#              returns at once; console output goes to
#              ud_run_thermal_comfort.<expnr>.out in the case directory
#           f  run on the current node in the foreground: the script blocks,
#              output goes to the terminal, the exit code is returned
#           c  submit to a compute node through the scheduler of <system>
#              (icl: PBS, ecmwf: Slurm)
#   system  common (default), icl or ecmwf
#
#   Everything after "--" is passed unchanged to ud_run_thermal_comfort.py;
#   --heights is required there. Run
#     tools/python/.venv/bin/python tools/python/udcomf/ud_run_thermal_comfort.py --help
#   for the full list. The Python chain writes its own log,
#   thermal_comfort.<expnr>.log in the case directory, in every start mode.
#
# Scheduler settings (environment or the case's config.sh):
#   COMFORT_NCPU      default: nompthreads of namoptions, else 8
#   COMFORT_WALLTIME  default 24:00:00
#   COMFORT_MEM       default 64gb (number followed by lowercase gb)
#   COMFORT_ACCOUNT   ecmwf only: Slurm account, default spnlthee
#   COMFORT_QOS       ecmwf only: Slurm QOS, default nf

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
	echo "usage: FROM THE TOP LEVEL DIRECTORY run:"
	echo "  u-dales/tools/ud_run_thermal_comfort.sh <PATH_TO_CASE> [start] [system] -- --heights H [H ...] [options]"
	echo "   start (optional): default l"
	echo "      l: current node, in the background; returns at once, console output in ud_run_thermal_comfort.<expnr>.out"
	echo "      f: current node, in the foreground; blocks, exit code returned"
	echo "      c: compute node through the scheduler (icl: PBS, ecmwf: Slurm)"
	echo "   system (optional): common (default), icl or ecmwf"
	echo "   options after --: passed to tools/python/udcomf/ud_run_thermal_comfort.py (--heights is required)"
}

default_ncpu=8

extract_nompthreads() {
	local namoptions_file=$1
	local line value
	if [ "$(grep -ic "nompthreads" "$namoptions_file" || true)" -ne 1 ]; then
		return 0
	fi
	line=$(grep -i "nompthreads" "$namoptions_file")
	value=${line#*=}
	value=${value%%!*}
	value=${value%%,*}
	value=${value%%/*}
	value="${value#"${value%%[![:space:]]*}"}"
	value="${value%"${value##*[![:space:]]}"}"
	if [[ "$value" =~ ^[0-9]+$ ]] && [ "$value" -gt 0 ]; then
		printf "%d\n" "$value"
	fi
}

if (( $# < 1 )); then
	echo "The path to the case/experiment folder must be set."
	usage
	echo "... execution terminated"
	exit 1
fi

case_path=$1
shift
start="l"
system="common"
if (( $# > 0 )) && [ "$1" != "--" ]; then
	start=$1
	shift
fi
if (( $# > 0 )) && [ "$1" != "--" ]; then
	system=$1
	shift
fi
if (( $# > 0 )) && [ "$1" == "--" ]; then
	shift
fi
wrapper_args=("$@")

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
	echo "  u-dales/tools/ud_run_thermal_comfort.sh $case_path c icl -- --heights 1.1"
	echo "  u-dales/tools/ud_run_thermal_comfort.sh $case_path c ecmwf -- --heights 1.1"
	echo "... execution terminated"
	exit 1
fi
if (( ${#wrapper_args[@]} == 0 )); then
	echo "No options for ud_run_thermal_comfort.py were given after --; at least --heights is required."
	usage
	echo "... execution terminated"
	exit 1
fi

# go to experiment directory
pushd "$case_path" > /dev/null
inputdir=$(pwd)
## set experiment number via path
iexpnr="${inputdir: -3}"
## read in additional variables
if [ -f config.sh ]; then
	# shellcheck disable=SC1091
	source config.sh
fi
export DA_TOOLSDIR="${DA_TOOLSDIR:-$script_dir}"
export DA_EXPDIR="${DA_EXPDIR:-$(cd "$inputdir/.." && pwd)}"
namoptions_file="$inputdir/namoptions.$iexpnr"
if [ ! -f "$namoptions_file" ]; then
	echo "Namelist file not found: $namoptions_file"
	exit 1
fi
popd > /dev/null

nompthreads=$(extract_nompthreads "$namoptions_file")
export COMFORT_NCPU="${COMFORT_NCPU:-${nompthreads:-$default_ncpu}}"
export COMFORT_WALLTIME="${COMFORT_WALLTIME:-24:00:00}"
export COMFORT_MEM="${COMFORT_MEM:-64gb}"
COMFORT_ACCOUNT="${COMFORT_ACCOUNT:-spnlthee}"
COMFORT_QOS="${COMFORT_QOS:-nf}"
if [[ ! "$COMFORT_MEM" =~ ^[0-9]+gb$ ]]; then
	echo "Memory requirement COMFORT_MEM must be set like 64gb"
	exit 1
fi

python_exe="$DA_TOOLSDIR/python/.venv/bin/python"
if [ ! -x "$python_exe" ]; then
	echo "Python virtual environment not found or not executable: $python_exe"
	echo "Set it up with: bash $DA_TOOLSDIR/python/setup_venv.sh <common|icl|ecmwf>"
	echo "and install the comfort dependency: $python_exe -m pip install -e '$DA_TOOLSDIR/python[comfort]'"
	echo "... execution terminated"
	exit 1
fi
wrapper="$DA_TOOLSDIR/python/udcomf/ud_run_thermal_comfort.py"
out_file="$inputdir/ud_run_thermal_comfort.$iexpnr.out"
chain_log="$inputdir/thermal_comfort.$iexpnr.log"
export PYTHONUNBUFFERED=1

## The command line, quoted once for the job scripts.
printf -v wrapper_cmd '%q ' "$python_exe" -u "$wrapper" --case-dir "$inputdir" "${wrapper_args[@]}"

write_job_body() {
cat <<JOBBODY
cd "$DA_TOOLSDIR/python"
export DA_TOOLSDIR="$DA_TOOLSDIR"
export DA_EXPDIR="$DA_EXPDIR"
export OMP_NUM_THREADS=$COMFORT_NCPU
export NUMBA_NUM_THREADS=$COMFORT_NCPU
export PYTHONUNBUFFERED=1
$wrapper_cmd > "$out_file" 2>&1 < /dev/null
JOBBODY
}

if [ "$start" == "c" ]; then
	cd "$inputdir"
	if [[ "$inputdir" =~ [[:space:]] ]] || [[ "$DA_TOOLSDIR" =~ [[:space:]] ]]; then
		echo "Paths must not contain spaces when submitting to a scheduler"
		echo "  case dir    = $inputdir"
		echo "  DA_TOOLSDIR = $DA_TOOLSDIR"
		exit 1
	fi
	if [ "$system" == "icl" ]; then
		###### Imperial HPC: PBS job script
		job_file="comfort-job.$iexpnr"
		{
cat <<PBSHEAD
#!/bin/bash
#PBS -l walltime=$COMFORT_WALLTIME
#PBS -l select=1:ncpus=$COMFORT_NCPU:mem=${COMFORT_MEM}
module load tools/prod
module load GCC/14.2.0
module load Python/3.13.1-GCCcore-14.2.0
PBSHEAD
		write_job_body
		} > "$job_file"
		if ! qsub "$job_file"; then
			echo "[ERROR] qsub failed; job script kept at $job_file" >&2
			exit 1
		fi
		echo "$job_file submitted. Console output: $out_file; chain log: $chain_log"
	else
		###### ECMWF Atos HPC2020: Slurm job script
		mem_directive=""
		if [ "$COMFORT_QOS" != "np" ]; then
			mem_directive="#SBATCH --mem=${COMFORT_MEM%gb}G"
		fi
		job_file="comfort-job.$iexpnr.slurm"
		{
cat <<SLURMHEAD
#!/bin/bash
#SBATCH --job-name=comfort-$iexpnr
#SBATCH --account=$COMFORT_ACCOUNT
#SBATCH --qos=$COMFORT_QOS
#SBATCH --time=$COMFORT_WALLTIME
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$COMFORT_NCPU
#SBATCH --hint=nomultithread
#SBATCH --output=comfort-job.$iexpnr.out
${mem_directive}
set -eo pipefail
## Runtime of the Python tools; keep in step with the ecmwf block of
## tools/build_preprocessing.sh.
module load prgenv/gnu gcc/14.2.0 python3/3.13.16-01
SLURMHEAD
		write_job_body
		} > "$job_file"
		if ! sbatch "$job_file"; then
			echo "[ERROR] sbatch failed; job script kept at $job_file" >&2
			exit 1
		fi
		echo "$job_file submitted (account $COMFORT_ACCOUNT, qos $COMFORT_QOS, $COMFORT_NCPU cpus, ${COMFORT_MEM}, $COMFORT_WALLTIME)."
		echo "Console output: $out_file; chain log: $chain_log"
	fi
elif [ "$start" == "f" ]; then
	###### RUN in the foreground on this node; exit code is the chain's
	cd "$DA_TOOLSDIR/python"
	export OMP_NUM_THREADS="$COMFORT_NCPU"
	export NUMBA_NUM_THREADS="$COMFORT_NCPU"
	set +e
	"$python_exe" -u "$wrapper" --case-dir "$inputdir" "${wrapper_args[@]}" < /dev/null
	status=$?
	set -e
	if [ "$status" -ne 0 ]; then
		echo "Thermal-comfort chain of $inputdir failed with exit code $status (log: $chain_log)" >&2
	fi
	exit "$status"
else
	###### RUN in the background on this node
	if [ "$system" == "ecmwf" ]; then
		echo "Note: this runs on the current node. ECMWF login nodes give a session very few cores and"
		echo "      end long processes; for anything but a small case use start mode c (Slurm job)."
	fi
	cd "$DA_TOOLSDIR/python"
	export OMP_NUM_THREADS="$COMFORT_NCPU"
	export NUMBA_NUM_THREADS="$COMFORT_NCPU"
	nohup "$python_exe" -u "$wrapper" --case-dir "$inputdir" "${wrapper_args[@]}" > "$out_file" 2>&1 < /dev/null &
	echo "Started in the background (PID $!). Console output: $out_file; chain log: $chain_log"
	cd "$DA_EXPDIR"
	cd ..
fi
