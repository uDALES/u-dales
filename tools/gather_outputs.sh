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
	echo "error: call scipt as `basename $0` experiment-directory."
	exit 0
fi

# Absolute path to this script directory
pushd $(dirname "${0}") > /dev/null
scriptdir=$(pwd -L)
popd > /dev/null
toolsdir=${scriptdir}  # assume same directory for nco_concatenate_field.sh

if [ -z $LOCAL_EXECUTE ]; then
    echo "cluster"
    module load intel-suite udunits nco/4.6.2
fi;

## go to files directory
cd ${datapath}
echo ${datapath}

## call loop for *DUMPS

for file in *dump.000.000.${expnr}.nc ; do
#for file in *dump* ; do
    echo ${file}
    if [ -f $file ]; then

        ## Gathering fields along spatial axis.
        echo "Gathering fields along spatial axis."

        dumps=${file%.000.000.${expnr}.nc}
        echo ${dumps}

        if [ $dumps == "fielddump" ]; then
            # ymparam="ym"
            ymparam="v,ym"

        elif [ $dumps == "tdump" ]; then
            ymparam="vt,vpwpt,upvpt,ym"

        elif [ $dumps == "slicedump" ]; then
            ymparam="v_2,v_20,ym"
        else
            ymparam="ym"
        fi

        outfile="${dumps}.${expnr}.nc"

        echo "We are in ${datapath}."
        echo "Gathering ${dumps} files with ym-dependent variables ${ymparam}."
        echo "Saving output to ${outfile}."

        ${toolsdir}/nco_concatenate_field.sh $dumps $ymparam $outfile
        echo "Merging done."

    fi

done

if [ -z $LOCAL_EXECUTE ]; then
    module unload intel-suite udunits nco/4.6.2
fi;
