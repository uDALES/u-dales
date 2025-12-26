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
# variables not specified in namoptions can be set here. It updates these in write_inputs.m and then runs the matlab code to produce the required text files.

# tg3315 20/07/2017, modified by SO 06/02/20, modified by DM 27/08/24

# NOTES 
# (1) if no forcing found in namoptions then applies the initial velocities and uses the pressure terms

# Assuming running from top-level project directory.
# u-dales/tools/write_inputs.sh <PATH_TO_CASE>

set -e

if (( $# < 1 )) 
then
    echo "The path to case/experiment folder must be set."
	echo "usage: FROM THE TOP LEVEL DIRECTORY run: u-dales/tools/write_inputs.sh <PATH_TO_CASE> (start)"
	echo "   start (optional): (c)ompute node or (l)ogin node"
	echo "... execution terminated"
    exit 1
fi

start=${2:-"x"}     # pass 'c' if needs to be run on hpc compute node, or 'l' if to be run on login node

# go to experiment directory
pushd $1
	inputdir=$(pwd)

	## set experiment number via path
	iexpnr="${inputdir: -3}"

	## read in additional variables
	if [ -f config.sh ]; then
   	 source config.sh
	else
	 echo "config.sh must be set inside $inputdir"
     exit 1
	fi

	## check if required variables are set
	if [ -z $DA_TOOLSDIR ]; then
	    echo "Script directory DA_TOOLSDIR must be set inside $inputdir/config.sh"
	    exit 1
	fi;
	if [ -z $DA_EXPDIR ]; then
		echo "Experiment directory DA_EXPDIR must be set $inputdir/config.sh"
		exit 1
	fi;

popd

if [ $start == "c" ]; then

	cd $inputdir

###### RUN MATLAB SCRIPT through HPC job script
cat <<EOF > pre-job.$iexpnr
#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=64gb

module load tools/prod
module load MATLAB/2024b
module load GCC/14.2.0

cd $DA_TOOLSDIR

export DA_TOOLSDIR=$DA_TOOLSDIR
export DA_EXPDIR=$DA_EXPDIR
export MATLAB_USE_USERWORK=0

nohup matlab -nodesktop -noFigureWindows -nosplash -nodisplay -r "expnr=$iexpnr; write_inputs; quit" > $inputdir/write_inputs.$iexpnr.log 2>&1 < /dev/null

EOF

## submit job.exp file to queue
	qsub pre-job.$iexpnr
	echo "pre-job.$iexpnr submitted."
else
	###### RUN MATLAB SCRIPT
	cd $DA_TOOLSDIR
	nohup matlab -nodesktop -noFigureWindows -nosplash -nodisplay -r "expnr=$iexpnr; write_inputs; quit" > $inputdir/write_inputs.$iexpnr.log 2>&1 < /dev/null &
	cd $DA_EXPDIR
	cd ..
fi



###### alter files in namoptions
# cd $DA_EXPDIR/$iexpnr

# ####### modified function for sed -i
# function sedi { if [[ "$OSTYPE" == "darwin"* ]]; then
#         		sed -i '' "$1" "$2"
# 		elif [[ "$OSTYPE" == "linux"* ]]; then
#         		sed -i "$1" "$2"
# 		fi;}

# function update_namoptions() {
#     local filename=$1
#     local varname=$2
#     local sub_offset=$3  # number of garbage lines in the specific file

#     if [ -f "$DA_EXPDIR/$iexpnr/${filename}" ]; then
#         local count=$(wc -l < "$DA_EXPDIR/$iexpnr/${filename}")
# 		count=$(($count-$sub_offset))
#     else
#         local count=0
#     fi

#     if grep -w -q "$varname" "$DA_EXPDIR/$iexpnr/namoptions.$iexpnr"; then
#         sedi "/^$varname =/s/.*/$varname = $count/g" "$DA_EXPDIR/$iexpnr/namoptions.$iexpnr"
#     else
#         sedi '/&WALLS/a\'$'\n'"$varname = $count"$'\n' "$DA_EXPDIR/$iexpnr/namoptions.$iexpnr"
#     fi
# }

# # Call the function for blocks
# update_namoptions "facet_sections_c.txt" "nfctsecs_c" 1
# update_namoptions "facet_sections_w.txt" "nfctsecs_w" 1
# update_namoptions "facet_sections_v.txt" "nfctsecs_v" 1
# update_namoptions "facet_sections_u.txt" "nfctsecs_u" 1

# update_namoptions "fluid_boundary_c.txt" "nbndpts_c" 2
# update_namoptions "fluid_boundary_w.txt" "nbndpts_w" 2
# update_namoptions "fluid_boundary_v.txt" "nbndpts_v" 2
# update_namoptions "fluid_boundary_u.txt" "nbndpts_u" 2

# update_namoptions "solid_c.txt" "nsolpts_c" 2
# update_namoptions "solid_w.txt" "nsolpts_w" 2
# update_namoptions "solid_v.txt" "nsolpts_v" 2
# update_namoptions "solid_u.txt" "nsolpts_u" 2

# update_namoptions "facets.inp.$iexpnr" "nfcts" 1
