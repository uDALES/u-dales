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

# Function to check which variables exist in NetCDF files and return only those that are present
check_existing_variables() {
    local proposed_vars="$1"
    local sample_file="$2"
    
    if [ ! -f "$sample_file" ]; then
        echo "Error: Sample file $sample_file not found"
        return 1
    fi
    if ! command -v ncdump >/dev/null 2>&1; then
        echo "Error: ncdump not found in PATH; it is needed to detect the staggered variables (load a NetCDF module)" >&2
        exit 1
    fi
    
    # Get list of variables in the file
    local available_vars=$(ncdump -h "$sample_file" 2>/dev/null | grep -E "^\s+(int|float|double|byte|char|short|long)" | awk '{print $2}' | sed 's/(.*//' | tr '\n' ',' | sed 's/,$//')
    
    # Check which proposed variables actually exist
    local existing_vars=""
    IFS=',' read -ra VAR_ARRAY <<< "$proposed_vars"
    
    for var in "${VAR_ARRAY[@]}"; do
        var=$(echo "$var" | tr -d ' ')  # remove whitespace
        if echo ",$available_vars," | grep -q ",$var,"; then
            if [ -z "$existing_vars" ]; then
                existing_vars="$var"
            else
                existing_vars="$existing_vars,$var"
            fi
        fi
    done
    
    echo "$existing_vars"
}

# List the distinct dump types matching a glob (e.g. "stats_*", "ins_?slice"), whatever the
# block number of their files: slice files exist only for the blocks that contain a plane,
# so a ".000." file cannot be assumed.
list_dump_types() {
    local pattern="$1"
    ls ${pattern}.???.${expnr}.nc 2>/dev/null | sed "s/\.[0-9][0-9][0-9]\.${expnr}\.nc$//" | sort -u
}

# First existing file of a dump type (used as the sample for check_existing_variables)
first_file_of() {
    ls ${1}.???.${expnr}.nc 2>/dev/null | head -n 1
}

if (( $# == 1 )) ; then
    datapath=$1
    expnr="${datapath: -3}"  ## set experiment number via path
else
	echo "error: call script as `basename $0` experiment-directory."
	exit 0
fi

# Absolute path to this script directory
pushd $(dirname "${0}") > /dev/null
scriptdir=$(pwd -L)
popd > /dev/null
toolsdir=${scriptdir}  # assume same directory for nco_concatenate_field.sh

## go to files directory
cd ${datapath}


######### ## Gathering fields along spatial axis y. ## #########

## call loop for *dumps
for file in *dump.*.000.${expnr}.nc ; do
    if [ -f $file ]; then
        dumps=${file%.000.${expnr}.nc}

        if [ ${dumps:0:8} == "mintdump" ]; then
            echo "Merging mintdump along y-direction."
	        proposed_ymparam="vt,ym"
	        ymparam=$(check_existing_variables "$proposed_ymparam" "$file")
        else
            ymparam="ym"
        fi

        outfile="${dumps}.${expnr}.nc"

        echo "We are in ${datapath}."
        echo "Gathering ${dumps} files with ym-dependent variables ${ymparam}."
        echo "Saving output to ${outfile}."

        ${toolsdir}/nco_concatenate_field_y.sh $dumps $ymparam $outfile
        echo "Merging done."
    fi
done

## call loop for ins_jslice* and stats_jslice* (one file per y-block, gathered along y)
for dumps in $(list_dump_types "ins_jslice") $(list_dump_types "stats_jslice") ; do
    file=$(first_file_of $dumps)
    if [ -f "$file" ]; then

        if [ $dumps == "ins_jslice" ]; then
            echo "Merging $dumps along y-direction."
            proposed_ymparam="ym"
            ymparam=$(check_existing_variables "$proposed_ymparam" "$file")
        elif [ $dumps == "stats_jslice" ]; then
            # time-averaged statistics on xz-planes: same variables (and staggering) as stats_t
            echo "Merging $dumps along y-direction."
            proposed_ymparam="v,vpwp,upvp,vsgs,vpthlp,vps1p,vps2p,vps3p,vps4p,ym"
            ymparam=$(check_existing_variables "$proposed_ymparam" "$file")
        else
            ymparam="ym"
        fi

        outfile="${dumps}.${expnr}.nc"

        echo "We are in ${datapath}."
        echo "Gathering ${dumps} files with ym-dependent variables ${ymparam}."
        echo "Saving output to ${outfile}."

        ${toolsdir}/nco_concatenate_field_y.sh $dumps $ymparam $outfile
        echo "Merging done."
    fi
done


######### ## Gathering fields along spatial axis x. ## #########

## call loop for *dumps
for file in *dump.000.${expnr}.nc ; do
    if [ -f $file ]; then
        dumps=${file%.000.${expnr}.nc}

        if [ $dumps == "mintdump" ]; then
            echo "Merging mintdump along x-direction."
            proposed_xmparam="ut,xm"
            xmparam=$(check_existing_variables "$proposed_xmparam" "$file")
        else
            xmparam="xm"
        fi

        outfile="${dumps}.${expnr}.nc"

        echo "We are in ${datapath}."
        echo "Gathering ${dumps} files with xm-dependent variables ${xmparam}."
        echo "Saving output to ${outfile}."

        ${toolsdir}/nco_concatenate_field_x.sh $dumps $xmparam $outfile
        echo "Merging done."
    fi
done

## call loop for stats_* (stats_jslice is gathered along y above)
for dumps in $(list_dump_types "stats_*") ; do
    file=$(first_file_of $dumps)
    if [ -f "$file" ]; then

        if [ $dumps == "stats_jslice" ]; then
            continue
        fi

        if [ $dumps == "stats_t" ] || [ $dumps == "stats_islice" ] || [ $dumps == "stats_kslice" ]; then
            # stats_islice/stats_kslice hold the same variables and staggering as stats_t
	        echo "Merging $dumps along x-direction."	
            proposed_xmparam="u,upwp,upvp,usgs,upthlp,ups1p,ups2p,ups3p,ups4p,xm"
            xmparam=$(check_existing_variables "$proposed_xmparam" "$file")
        elif [ $dumps == "stats_yt" ] || [ $dumps == "stats_y" ]; then
	        echo "Merging $dumps along x-direction."	
            proposed_xmparam="u,upwp,uw,usgs,xm"
            xmparam=$(check_existing_variables "$proposed_xmparam" "$file")
        elif [ $dumps == "stats_tree" ]; then
	        echo "Merging $dumps along x-direction."	
            proposed_xmparam="tr_u,xm"
            xmparam=$(check_existing_variables "$proposed_xmparam" "$file")
        else
            xmparam="xm"
        fi

        outfile="${dumps}.${expnr}.nc"

        echo "We are in ${datapath}."
        echo "Gathering ${dumps} files with xm-dependent variables ${xmparam}."
        echo "Saving output to ${outfile}."

        ${toolsdir}/nco_concatenate_field_x.sh $dumps $xmparam $outfile
        echo "Merging done."
    fi
done

## call loop for ins_islice* and ins_kslice* (but not ins_jslice*)
for dumps in $(list_dump_types "ins_islice") $(list_dump_types "ins_kslice") ; do
    file=$(first_file_of $dumps)
    if [ -f "$file" ]; then

        if [ $dumps == "ins_kslice" ]; then
            echo "Merging $dumps along x-direction."
            proposed_xmparam="u,xm"
            xmparam=$(check_existing_variables "$proposed_xmparam" "$file")
        elif [ $dumps == "ins_islice" ]; then
            echo "Merging $dumps along x-direction."
            proposed_xmparam="xm"
            xmparam=$(check_existing_variables "$proposed_xmparam" "$file")
        else
            xmparam="xm"
        fi

        outfile="${dumps}.${expnr}.nc"

        echo "We are in ${datapath}."
        echo "Gathering ${dumps} files with xm-dependent variables ${xmparam}."
        echo "Saving output to ${outfile}."

        ${toolsdir}/nco_concatenate_field_x.sh $dumps $xmparam $outfile
        echo "Merging done."
    fi
done

## call loop for ins_field*
for file in ins_field.000.${expnr}.nc ; do
    if [ -f $file ]; then
        dumps=${file%.000.${expnr}.nc}

        echo "Merging $dumps along x-direction."
        proposed_xmparam="u,tau_x,mask_u,pup,xm"
        xmparam=$(check_existing_variables "$proposed_xmparam" "$file")
        outfile="${dumps}.${expnr}.nc"

        echo "We are in ${datapath}."
        echo "Gathering ${dumps} files with xm-dependent variables ${xmparam}."
        echo "Saving output to ${outfile}."

        ${toolsdir}/nco_concatenate_field_x.sh $dumps $xmparam $outfile
        echo "Merging done."
    fi
done
