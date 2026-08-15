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

## call loop for *DUMPS
for file in *dump.*.000.${expnr}.nc ; do
    if [ -f $file ]; then
        ## echo "Gathering fields along spatial axis y."
        dumps=${file%.000.${expnr}.nc}

        if [ ${dumps:0:9} == "fielddump" ]; then
            echo "Merging fielddump along y-direction."
	        #ymparam="v,tau_y,ym"
	        ymparam="v,ym"
        elif [ ${dumps:0:8} == "mintdump" ]; then
            echo "Merging mintdump along y-direction."
	        ymparam="vt,ym"
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

## call loop for stats_*
for file in stats_*.*.000.${expnr}.nc ; do
    if [ -f $file ]; then
        ## Gathering fields along spatial axis y.
        dumps=${file%.000.${expnr}.nc}
        base=${dumps%%.*}

        if [ "$base" == "stats_t" ]; then
            echo "Merging $base along y-direction."
	        ymparam="v,vpwp,upvp,vsgs,ym"
        elif [ "$base" == "stats_tree" ]; then
            echo "Merging $base along y-direction."
	        ymparam="tr_v,ym"
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

## call loop for islice* and kslice*
for file in *slice.*.000.${expnr}.nc ; do
    if [ -f $file ]; then
        dumps=${file%.000.${expnr}.nc}

        if [ ${dumps:0:6} == "kslice" ] || [ ${dumps:0:6} == "islice" ]; then
            echo "Merging $dumps along y-direction."
	        ymparam="v,ym"
        fi

        outfile="${dumps}.${expnr}.nc"

        echo "We are in ${datapath}."
        echo "Gathering ${dumps} files with ym-dependent variables ${ymparam}."
        echo "Saving output to ${outfile}."

        ${toolsdir}/nco_concatenate_field_y.sh $dumps $ymparam $outfile
        echo "Merging done."

        if [ ${dumps:0:6} == "islice" ]; then
            # remove procx from name
            mv $outfile "islice.${expnr}.nc"
        fi
    fi
done

## call loop for jslice*
for file in jslice.???.???.${expnr}.nc ; do
	if [ -f $file ]; then
		procx=${file:7:3}
        # remove procy from name
        cp $file "jslice.${procx}.${expnr}.nc"
	fi
done


for file in *dump.000.${expnr}.nc ; do
    if [ -f $file ]; then
        ## echo "Gathering fields along spatial axis x."

        dumps=${file%.000.${expnr}.nc}

        if [ $dumps == "fielddump" ]; then
	        echo "Merging fielddump along x-direction."
            #xmparam="u,tau_x,xm"
	        xmparam="u,xm"
        elif [ $dumps == "mintdump" ]; then
            echo "Merging mintdump along x-direction."
            xmparam="ut,xm"
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

## Gathering fields along spatial axis x.
for file in stats_*.000.${expnr}.nc ; do
    if [ -f $file ]; then
        dumps=${file%.000.${expnr}.nc}

        # Skip raw files accidentally matched (e.g. stats_t.000.000.expnr.nc gives dumps=stats_t.000)
        [[ "$dumps" == *.* ]] && continue  # if dumps contains a dot anywhere (e.g. stats_t.000), skip to next file

        if [ $dumps == "stats_t" ]; then
	        echo "Merging $dumps along x-direction."	
            xmparam="u,upwp,upvp,usgs,xm"
        elif [ $dumps == "stats_yt" ] || [ $dumps == "stats_y" ]; then
	        echo "Merging $dumps along x-direction."	
            xmparam="u,upwp,uw,usgs,xm"
        elif [ $dumps == "stats_tree" ]; then
	        echo "Merging $dumps along x-direction."	
            xmparam="tr_u,xm"
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

for file in *slice.000.${expnr}.nc ; do
    if [ -f $file ]; then
        dumps=${file%.000.${expnr}.nc}

        if [ $dumps == "kslice" ] || [ $dumps == "jslice" ]; then
            echo "Merging $dumps along x-direction."
            xmparam="u,xm"
        fi

        outfile="${dumps}.${expnr}.nc"

        echo "We are in ${datapath}."
        echo "Gathering ${dumps} files with xm-dependent variables ${xmparam}."
        echo "Saving output to ${outfile}."

        ${toolsdir}/nco_concatenate_field_x.sh $dumps $xmparam $outfile
        echo "Merging done."

	    if [ $dumps == "jslice" ]; then
            rm ${dumps}.???.${expnr}.nc
    	fi
    fi
done
