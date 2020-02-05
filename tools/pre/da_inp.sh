#!/bin/bash                                                                                        
# script for setting key namoptions and writing *.inp.* files
# variables not specified in namoptions can be set here. It updates these in da_inp.m and then runs the matlab code to produce the required text files.

# tg3315 20/07/2017

# edit needs to be done so is irrespective of white spaces

# NOTES 
# (1) if no forcing found in namoptions then applies the initial velocities and uses the pressure terms
# (2) if libm is false it should overwrite and produce no blocks

iexpnr=$1

# Assuming running from top-level project directory.
if [ -f config.sh ]; then
    source config.sh
fi

#DA_GITDIR=$(git rev-parse --show-toplevel)
#DA_EXPDIR=$DA_GITDIR/experiments
#DA_PREDIR=$DA_GITDIR/u-dales/tools/pre
export PATH=$DA_PREDIR:$PATH

cd $DA_EXPDIR/$iexpnr

#if [[ "$OSTYPE" == "darwin"* ]]; then
#	sedi="sed -i ''"
#elif [[ "$OSTYPE" == "linux-gnu" ]]; then
#	sedi="sed -i"
#fi

function sedi { if [[ "$OSTYPE" == "darwin"* ]]; then
        		sed -i '' "$1" "$2"
		elif [[ "$OSTYPE" == "linux-gnu" ]]; then
        		sed -i "$1" "$2"
		fi;}

####### set iexpnir in matlab file
sedi "/expnr = '/s/.*/expnr = '$iexpnr';/g" $DA_PREDIR"/da_inp.m"
###### set # CPUS from execute to test domain size !edit : should maybe multiply by nnode
sedi  "/CPUS = /s/.*/CPUS = $(grep -m 1 'ncpu=' ../../u-dales/tools/utils/local_execute.sh | cut -d "=" -f 2 | cut -d " " -f 1 | tr -d ' ');       % # cpus/g" $DA_PREDIR"/da_inp.m" 

###### RUN MATLAB SCRIPT FOR .inp. files

#octave $DA_TOPDIR/bin/da_inp.m
cd $DA_PREDIR/
matlab -nodesktop -nosplash -r "da_inp; quit"
cd $DA_EXPDIR/$iexpnr

##### alter files in namoptions (thl_top, nblocks etc.)
nblocks=$(wc -l < $DA_EXPDIR/$iexpnr/blocks.inp.$iexpnr)
nblocks=$(($nblocks-2))

if grep -q nblocks $DA_EXPDIR/$iexpnr"/namoptions."$iexpnr; then
	sedi "/nblocks/s/.*/nblocks    =  $nblocks/g" $DA_EXPDIR/$iexpnr"/namoptions."$iexpnr
else
	sedi '/&DOMAIN/a\'$'\n''nblocks    =  '$nblocks''$'\n' $DA_EXPDIR/$iexpnr"/namoptions."$iexpnr
fi

nfcts=$(wc -l < $DA_EXPDIR/$iexpnr/facets.inp.$iexpnr)
nfcts=$(($nfcts-1))

if grep -q nfcts $DA_EXPDIR/$iexpnr"/namoptions."$iexpnr; then
	sedi "/nfcts/s/.*/nfcts      = $nfcts/g" $DA_EXPDIR/$iexpnr"/namoptions."$iexpnr
else
	sedi '/&ENERGYBALANCE/a\'$'\n''nfcts      = '$nfcts''$'\n' $DA_EXPDIR/$iexpnr"/namoptions."$iexpnr
fi
