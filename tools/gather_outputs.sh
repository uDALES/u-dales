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

        if [ ${dumps:0:9} == "fielddump" ]; then
            echo "Merging fielddump along y-direction."
	        proposed_ymparam="v,tau_y,ym"
	        ymparam=$(check_existing_variables "$proposed_ymparam" "$file")
        elif [ ${dumps:0:8} == "mintdump" ]; then
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

## call loop for ins_jslice*
for file in ins_?slice.000.${expnr}.nc ; do
    if [ -f $file ]; then
        dumps=${file%.000.${expnr}.nc}

        # Skip processing for ins_islice and ins_kslice
        if [ $dumps == "ins_islice" ] || [ $dumps == "ins_kslice" ]; then
            continue
        fi

        if [ $dumps == "ins_jslice" ]; then
            echo "Merging $dumps along y-direction."
            proposed_ymparam="ym"
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

        if [ $dumps == "fielddump" ]; then
	        echo "Merging fielddump along x-direction."
	        proposed_xmparam="u,tau_x,xm"
	        xmparam=$(check_existing_variables "$proposed_xmparam" "$file")
        elif [ $dumps == "mintdump" ]; then
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

## call loop for stats_*
for file in stats_*.000.${expnr}.nc ; do
    if [ -f $file ]; then
        dumps=${file%.000.${expnr}.nc}

        if [ $dumps == "stats_t" ]; then
	        echo "Merging $dumps along x-direction."	
            proposed_xmparam="u,upwp,upvp,usgs,ups1p,ups2p,ups3p,xm"
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
for file in ins_?slice.000.${expnr}.nc ; do
    if [ -f $file ]; then
        dumps=${file%.000.${expnr}.nc}

        # Skip processing for ins_jslice
        if [ $dumps == "ins_jslice" ]; then
            continue
        fi

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
