#!/usr/bin/env bash

if (( $# == 1 )) ; then
    iexpnr=$1
    type="field"
elif (( $# == 2 )) ; then
    iexpnr1=$1
    iexpnr2=$2
    type="time"
else
	echo "error: call scipt as `basename $0` expnr for spatial merging or as `basename $0` expnr1 expnr2 for merging two experiments."
	exit 0
fi

## Merging fields along spatial axis.
if [ $type == "field" ]; then
    
    echo "Merging fields along spatial axis."
    
    ## go to files directory
    cd ${DA_WORKDIR}/${iexpnr}

    ## call loop for FIELDDUMPS
    dumps="fielddump"
    # ymparam="v,ym"
    ymparam="ym"
    outfile="${dumps}.${iexpnr}.nc"

    echo "We are in ${DA_WORKDIR}/${iexpnr}."
    echo "Merging ${dumps} files with ym-dependent variables ${ymparam}."
    echo "Saving output to ${outfile}."

    . da_merge.sh $dumps $ymparam $outfile
    echo "Merging done."

    ### SAME LOOP FOR SLICEDUMPS
    # dumps=slicedump
    # ymparam="v_2,v_20,ym"
    # outfile="${dumps}.${iexpnr}.nc"
    #
    # echo "We are in ${DA_WORKDIR}/${iexpnr}.
    # Merging ${dumps} files with ym-dependent variables ${ymparam}.
    # Saving output to ${outfile}."
    #
    # . da_merge.sh $dumps $ymparam $outfile
    #
    # echo "Merging done."
    
elif  [ $type == "time" ]; then

    echo "Merging two fields along time."
    
    ## go to files directory
    cd "${DA_WORKDIR}/${iexpnr2}"
    
    ## call loop for *DUMPS
    dumpslist=$(ls *dump.${iexpnr2}.nc)
    # dumpslist="ztdump.${iexpnr2}.nc fielddump.${iexpnr2}.nc"
    
    for file in $dumpslist ; do
    
        
        dump=${file%.${iexpnr2}.nc}
        
        mv $file "short-${file}"
    
        # outfile="${dump}.${iexpnr2}-merged.nc"
        outfile=$file
        
        echo "Merging ${dump} files. Saving output to ${DA_WORKDIR}/${iexpnr2}/${outfile}."
        
        ncrcat "../${iexpnr1}/${dump}.${iexpnr1}.nc" "short-${file}" ${outfile}
        
        echo "Done merging ${dump}."    
        
        echo "Testing now for time overlap."

        times=$(ncks -s "%.0f\n" -H -C -v time $outfile)

        told=0
        counter=0

        for t in $times; do 
            if [[ $told -lt $t ]] ; then 
                # loop for when there is an overlap where the overlap finishes
                if [[ $t -eq $tdouble ]] ; then 
                    tcont=$((counter+1))
                fi
            else
                tdouble=$told
                tcut=$((counter-1))
            fi
            # increase counters
            told=$t
            ((counter++))
        done
        
        if [[ ! -z $tcut ]] ; then
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
       
else
    echo "Something went wrong."

fi