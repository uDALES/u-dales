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

if [ -z $LOCAL_EXECUTE ]; then
    module load nco gsl
fi;

## go to files directory
cd ${datapath}

## call loop for *DUMPS
for file in *dump.*.000.${expnr}.nc ; do
    if [ -f $file ]; then

        ## Gathering fields along spatial axis.
        #echo "Gathering fields along spatial axis."

        dumps=${file%.000.${expnr}.nc}

        if [ ${dumps:0:9} == "fielddump" ]; then
            echo "Merging fielddump along y-direction."
	    #ymparam="v,tau_y,ym"
	    ymparam="v,ym"

        elif [ ${dumps:0:5} == "tdump" ]; then
            echo "Merging tdump along y-direction."
	    ymparam="vt,vpwpt,upvpt,ym"

        elif [ ${dumps:0:8} == "mintdump" ]; then
            echo "Merging mintdump along y-direction."
	    ymparam="vt,ym"

        elif [ ${dumps:0:10} == "kslicedump" ]; then
            echo "Merging kslicedump along y-direction."
	    ymparam="v_kslice,ym"

        elif [ ${dumps:0:10} == "islicedump" ]; then
            echo "Merging islicedump along y-direction."
	    ymparam="v_islice,ym"

        else
            ymparam="ym"
        fi

        outfile="${dumps}.${expnr}.nc"

        echo "We are in ${datapath}."
        echo "Gathering ${dumps} files with ym-dependent variables ${ymparam}."
        echo "Saving output to ${outfile}."

        ${toolsdir}/nco_concatenate_field_y.sh $dumps $ymparam $outfile
        echo "Merging done."

	if [ ${dumps:0:10} == "islicedump" ]; then
	    # remove procx from name
            mv $outfile "islicedump.${expnr}.nc"
    	fi

    fi

done


#for file in islicedump.???.${expnr}.nc ; do
#        if [ -f $file ]; then
#                procy=${file:15:3}
#                # remove procx from name
#                mv $file "islicedump.${expnr}.nc"
#
#        fi
#done


for file in jslicedump.???.???.${expnr}.nc ; do
	if [ -f $file ]; then
		procx=${file:11:3}
                # remove procy from name
                cp $file "jslicedump.${procx}.${expnr}.nc"
	fi
done


for file in *dump.000.${expnr}.nc ; do

    if [ -f $file ]; then

        ## Gathering fields along spatial axis.
        #echo "Gathering fields along spatial axis."

        dumps=${file%.000.${expnr}.nc}

        if [ $dumps == "fielddump" ]; then
	    echo "Merging fielddump along x-direction."
            #xmparam="u,tau_x,xm"
	    xmparam="u,xm"

        elif [ $dumps == "tdump" ]; then
	    echo "Merging tdump along x-direction."	
            xmparam="ut,upwpt,upvpt,xm"

	elif [ $dumps == "mintdump" ]; then
            echo "Merging mintdump along x-direction."
	    xmparam="ut,xm"

	elif [ $dumps == "kslicedump" ]; then
            echo "Merging kslicedump along x-direction."
            xmparam="u_kslice,xm"

        elif [ $dumps == "jslicedump" ]; then
            echo "Merging jslicedump along x-direction."
            xmparam="u_jslice,xm"

        else
            xmparam="xm"
        fi

        outfile="${dumps}.${expnr}.nc"

        echo "We are in ${datapath}."
        echo "Gathering ${dumps} files with xm-dependent variables ${xmparam}."
        echo "Saving output to ${outfile}."

        ${toolsdir}/nco_concatenate_field_x.sh $dumps $xmparam $outfile
        echo "Merging done."

	if [ $dumps == "jslicedump" ]; then
            rm ${dumps}.???.${expnr}.nc
    	fi

    fi

done


if [ -z $LOCAL_EXECUTE ]; then
    module unload nco
fi;

