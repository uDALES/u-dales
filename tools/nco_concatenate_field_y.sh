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

## call: dumptype, ym-variables(list separated by comma, no space, include ym), out.nc
## e.g. nco_concatenate_field.sh fielddump v,ym out.nc

if (( $# == 3 )) ; then
	dumps=$1
	ymparam=$2
	outfile=$3
else
	echo "Wrong function call. Call nco_concatenate_field.sh dumptype ym-variables out.nc"
	exit 1
fi

n=$(ls ${dumps}.*.*.nc | wc -l)
## start if more than one *dumps
if (( n > 1 )) ; then

    if [ $ymparam = "ym" ] ; then
        echo "No variable dependent on ym other than ym itself. omitting ym."
        for file in ${dumps}.*.*.nc; do
            ncpdq -64 -a yt,time,zt,zm,xt,xm -C -x -v ${ymparam} ${file} "tmp${file#${dumps}}"
        done
        ncrcat -64 tmp.*.nc fieldtmp2.nc
        rm tmp.*.nc

        ## revert to time as record dimension
        ncpdq -64 -a time,zt,zm,yt,ym,xt,xm fieldtmp2.nc ${outfile}
        rm fieldtmp2.nc

    else
        ## start with the ym-variables and gather them in a new file
        for file in ${dumps}.*.*.nc; do
                ncpdq -64 -a ym,time,zt,zm,xt,xm -v ${ymparam} ${file} "tmp1${file#${dumps}}"
        done
        ncrcat -64 tmp1.*.nc fieldtmp1.nc
        rm tmp1.*.nc

        ## take only remaining variables and gather them in a new file
        for file in ${dumps}.*.*.nc; do
            ncpdq -64 -a yt,time,zt,zm,xt,xm -C -x -v ${ymparam} ${file} "tmp${file#${dumps}}"
        done
        ncrcat -64 tmp.*.nc fieldtmp2.nc
        rm tmp.*.nc

        ## revert to time as record dimension. needs to be done before merging
        ncpdq -64 -a time,zt,zm,yt,ym,xt,xm fieldtmp1.nc fieldtmp3.nc
        ncpdq -64 -a time,zt,zm,yt,ym,xt,xm fieldtmp2.nc ${outfile}
        rm fieldtmp1.nc fieldtmp2.nc

        ## merges the two files
        ncks -64 -A fieldtmp3.nc ${outfile}
        rm fieldtmp3.nc
    fi
fi

