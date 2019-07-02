#!/bin/bash

## set environment variables
topdir=${DA_TOPDIR}
#utildir=${DA_UTILSDIR}
srcdir=${DA_SRCDIR}
expdir=${DA_EXPDIR}
workdir=${DA_WORKDIR}
location=${PLATFORM}

## set platform dependent variables

if [ "$location" == "hpc" ] ||  [ "$location" == "cx1" ] || [ "$location" == "cx2" ] ; then
    exe=dalesurban_hpc
    ncpu=12
    nnode=5
    walltime="24:00:00"
    mem=20
    queue=""  # do not select a queue: line empty
    read modules < $utildir/dalesmodules # read in list of modules


elif [ "$location" == "local" ] ; then
    exe=dalesurban_local
    ncpu=4
    
elif [ "$location" == "macos" ] ; then
    exe=dalesurban_macos
    ncpu=1
else
    echo "location $location not specified"
    exit 1
fi

#################################################


