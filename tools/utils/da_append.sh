#!/usr/bin/env bash

set -e

if (( $# == 2 )) ; then
    datapath1=$1
    expnr1="${datapath1: -3}"  ## set experiment number via path
    datapath2=$2
    expnr2="${datapath2: -3}"  ## set experiment number via path
else
	echo "error: call scipt as `basename $0` experiment-1-directory experiment-2-directory for appending two experiments."
	exit 0
fi

if [ -z $LOCAL_EXECUTE ]; then
    echo "cluster"
    module load intel-suite udunits nco/4.6.2
fi;

echo "Merging two fields along time."

## get absolute path for first data directory
pushd ${datapath1}  > /dev/null
datapath1_update=$(pwd -L)
popd  > /dev/null

## go to files directory
cd "${datapath2}"

## call loop for *DUMPS
# dumpslist="ztdump.${expnr2}.nc"
# dumpslist="fielddump.${expnr2}.nc"
dumpslist=$(ls *dump.${expnr2}.nc)
# dumpslist="ztdump.${expnr2}.nc fielddump.${expnr2}.nc"

for file in $dumpslist ; do

    
    dump=${file%.${expnr2}.nc}
    
    mv $file "short-${file}"

    # outfile="${dump}.${iexpnr2}-merged.nc"
    outfile=$file
    
    echo "Merging ${dump} files. Saving output to ${datapath2}/${outfile}."
    
    ncrcat "${datapath1_update}/${dump}.${expnr1}.nc" "short-${file}" ${outfile}
    
    echo "Done merging ${dump}."    
    
    # rename so that this is not picked up as nc file
    mv "short-${file}" "short-${file%.nc}"
    
    echo "Testing now for time overlap."

    times=$(ncks -s "%.0f\n" -H -C -v time $outfile)

    told=0
    counter=0
tcut=""  # setting tcut to empty string, testing later whether it stays like that
    overlap=0  # setting overlap to false
for t in $times; do
        if [[ $told -lt $t ]] ; then  # as long as monotone increase
    # loop for when an overlap occured, check where the overlap finishes
    if [[ $overlap -eq 1 ]] ; then  # if overlap has occured in previous times
        if [[ $t -ge $tdouble ]] ; then  # current time step greater or equal to overlap time
            tcont=$((counter+1))
            overlap=0  # overlap time fixed, stop checking
        fi
            fi
        else  # not monotone!
            tdouble=$told
            tcut=$((counter-1))
    overlap=1
        fi
        # increase counters
        told=$t
        ((counter++))
    done
    
    if [[ ! -z $tcut ]] ; then  # testing if tcut is not empty
        echo "Time overlap occurred."
        echo "Take field until time index $tcut and continue with index $tcont. Cutting file now." 
        
        outfile2="long-${outfile}"
        mv $outfile $outfile2
        
        ncks -d time,,$tcut -d time,$tcont, $outfile2 $outfile
        echo "Done fixing time overlap."
    else
        echo "No time overlap."
    fi

done
       
if [ -z $LOCAL_EXECUTE ]; then
    module unload intel-suite udunits nco/4.6.2
fi;
