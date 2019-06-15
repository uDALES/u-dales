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

datapath=${EPHEMERAL}
# datapath = ${DA_WORKDIR}

module load intel-suite udunits nco/4.6.2

## Merging fields along spatial axis.
if [ $type == "field" ]; then
    
    echo "Merging fields along spatial axis."
    
    ## go to files directory
    cd ${datapath}/${iexpnr}
    
    ## call loop for *DUMPS
    dumpslist=$(ls *dump.000.${iexpnr}.nc)
    
    for file in $dumpslist ; do
    
        dumps=${file%.000.${iexpnr}.nc}

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
        
        outfile="${dumps}.${iexpnr}.nc"

        echo "We are in ${datapath}/${iexpnr}."
        echo "Merging ${dumps} files with ym-dependent variables ${ymparam}."
        echo "Saving output to ${outfile}."

        . da_merge.sh $dumps $ymparam $outfile
        echo "Merging done."

    done
    
    
elif  [ $type == "time" ]; then

    echo "Merging two fields along time."
    
    ## go to files directory
    cd "${datapath}/${iexpnr2}"
    
    ## call loop for *DUMPS
    # dumpslist="ztdump.${iexpnr2}.nc"
    # dumpslist="fielddump.${iexpnr2}.nc"
    dumpslist=$(ls *dump.${iexpnr2}.nc)
    # dumpslist="ztdump.${iexpnr2}.nc fielddump.${iexpnr2}.nc"
    
    for file in $dumpslist ; do
    
        
        dump=${file%.${iexpnr2}.nc}
        
        mv $file "short-${file}"
    
        # outfile="${dump}.${iexpnr2}-merged.nc"
        outfile=$file
        
        echo "Merging ${dump} files. Saving output to ${datapath}/${iexpnr2}/${outfile}."
        
        ncrcat "../${iexpnr1}/${dump}.${iexpnr1}.nc" "short-${file}" ${outfile}
        
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
       
else
    echo "Something went wrong."

fi

module unload intel-suite udunits nco/4.6.2
