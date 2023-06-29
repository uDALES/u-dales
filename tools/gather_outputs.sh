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
    module load nco
fi;

## go to files directory
cd ${datapath}

## call loop for *DUMPS
for file in *dump.*.000.${expnr}.nc ; do
    if [ -f $file ]; then

        ## Gathering fields along spatial axis.
        echo "Gathering fields along spatial axis."

        dumps=${file%.000.${expnr}.nc}

        if [ ${dumps:0:9} == "fielddump" ]; then
            ymparam="v,tau_y,ym"

        elif [ ${dumps:0:5} == "tdump" ]; then
            ymparam="vt,vpwpt,upvpt,ym"

        elif [ ${dumps:0:9} == "slicedump" ]; then
            ymparam="v_2,v_20,ym"
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

for file in *dump.000.${expnr}.nc ; do

    if [ -f $file ]; then

        ## Gathering fields along spatial axis.
        echo "Gathering fields along spatial axis."

        dumps=${file%.000.${expnr}.nc}

        if [ $dumps == "fielddump" ]; then
	          #xmparam="xm"
            xmparam="u,tau_x,xm"

        elif [ $dumps == "tdump" ]; then
            xmparam="ut,upwpt,upvpt,xm"

        elif [ $dumps == "slicedump" ]; then
            xmparam="u_2,u_20,xm"
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


if [ -z $LOCAL_EXECUTE ]; then
    module unload nco
fi;

