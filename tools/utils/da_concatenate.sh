#!/usr/bin/env bash

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
utilsdir=${scriptdir}  # assume same utils dir for da_merge.sh

if [ -z $LOCAL_EXECUTE ]; then
    echo "cluster"
    module load intel-suite udunits nco/4.6.2
fi;

## Merging fields along spatial axis.

echo "Merging fields along spatial axis."

## go to files directory
cd ${datapath}

## call loop for *DUMPS
dumpslist=$(ls *dump.000.${expnr}.nc)

for file in $dumpslist ; do

    dumps=${file%.000.${expnr}.nc}

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
    echo "Merging ${dumps} files with ym-dependent variables ${ymparam}."
    echo "Saving output to ${outfile}."

    ${utilsdir}/da_merge.sh $dumps $ymparam $outfile
    echo "Merging done."

done

if [ -z $LOCAL_EXECUTE ]; then
    module unload intel-suite udunits nco/4.6.2
fi;
