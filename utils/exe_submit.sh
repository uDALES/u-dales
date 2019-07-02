#!/bin/bash

## go to inps directory
cd ${expdir}/${iexpnr}

## check if exe is already in folder, then use it, otherwise copy from source
    if [ -x ${exe} ]; then
       echo "using existing executable in directory."
    else
       cp "${srcdir}/dales-urban/${exe}" "${expdir}/${iexpnr}"
       echo "copy new executable from source directory."
    fi
	
    echo "running $exe on $location"
    echo "running with $nnode x $ncpu cpus"

## test experiment number ##
if [ "$(ls -A *.${iexpnr})" ]; then
     echo "experiment numbers and files match"
else
    echo "files do not match experiment number"
    exit 1
fi

if grep iexpnr namoptions.${iexpnr} | grep -q ${iexpnr} ; then
    echo "experiment number matches with namoptions"
else
    echo "experiment numbers does not match with namoptions"
    exit 1
fi

## run on cx1 or cx2
if [ "$location" == "hpc" ] || [ "$location" == "cx1" ] || [ "$location" == "cx2" ] ; then

    ## note: necessary libraries have to be defined in the jobfile or bash profile
    ## copy generic job file to exp directory
    cp "${utildir}/job.${location}" "job.${iexpnr}"

    ## substituting strings in exp job file to current values
    sed -i "s/iexpnr=.../iexpnr=$iexpnr/g" "job.${iexpnr}"
    sed -i "s/exe=.*/exe=${exe}/g" "job.${iexpnr}"
    expdir=$(echo "$expdir" | sed 's/\//\\\//g')  # escaping all slashes
    sed -i "s/expdir=.*/expdir=${expdir}/g" "job.${iexpnr}"
    workdir=$(echo "$workdir" | sed 's/\//\\\//g')  # escaping all slashes
    sed -i "s/workdir=.*/workdir=${workdir}/g" "job.${iexpnr}"
    sed -i "s/select=.*n/select=${nnode}:n/g" "job.${iexpnr}"
    sed -i "s/ncpus=.*:/ncpus=${ncpu}:/g" "job.${iexpnr}"
    sed -i "s/walltime=.*/walltime=${walltime}/g" "job.${iexpnr}"
    sed -i "s/mem=.*/mem=${mem}gb/g" "job.${iexpnr}"
    sed -i "s/.*#PBS -q.*/${queue}/g" "job.${iexpnr}"
    modules=$(echo "$modules" | sed 's/\//\\\//g')  # escaping all slashes
    sed -i "s/.*module load.*/module load ${modules}/g" "job.${iexpnr}"

    echo "submitting job.${iexpnr} to ${location}"
    qsub "job.${iexpnr}"

## run locally (ubuntu or macos)
elif [ "$location" == "local" ] || [ "$location" == "macos" ] ; then

    outdir=${workdir}/${iexpnr}
    
    ## create output directory
    echo "creating outdir"
    echo "----------------"
    [[ -d ${outdir} ]] || mkdir ${outdir}

    echo "copying files"
    cp *inp.$iexpnr $outdir
    cp namoptions.$iexpnr $outdir
    cp $exe $outdir

    #checking if warmstart  #I'm sure there's better code, but it works.. ils13
    if grep warmstart namoptions.$iexpnr | grep -q .true. ; then  #look in namoptions for the string "warmstart" combined with ".true."
        echo "warmstart"
        #look in namoptions for the string "startfile", extract the substring in singlequotes "'" and replace the "x"s with "?"s
        startfiles=$(grep "startfile" namoptions.$iexpnr | sed "s/[^']*'\([^']*\).*/\1/" | sed "s/x/?/g")
        echo ${startfiles}
        cp ${startfiles} ${outdir}
    fi

    echo "going to $outdir"
    cd ${outdir}

    # execute the program
    echo "execute program"
    echo "mpiexec ./$exe namoptions.$iexpnr > output.$iexpnr 2>&1"
    mpiexec -n $ncpu ./$exe namoptions.$iexpnr > output.$iexpnr

    echo "... job.$iexpnr done"

else
    echo "execution failed"
    exit 1
fi
