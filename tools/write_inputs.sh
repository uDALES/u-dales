#!/bin/bash                                                                                        
# script for setting key namoptions and writing *.inp.* files
# variables not specified in namoptions can be set here. It updates these in write_inputs.m and then runs the matlab code to produce the required text files.

# tg3315 20/07/2017, modified by SO 06/02/20

# NOTES 
# (1) if no forcing found in namoptions then applies the initial velocities and uses the pressure terms
# (2) if libm is false it should overwrite and produce no blocks

iexpnr=$1

# Assuming running from top-level project directory.
if [ -f ./experiments/$iexpnr/config.sh ]; then
    	echo "found config.sh"
	source ./experiments/$iexpnr/config.sh
fi

export PATH=$DA_TOOLSDIR:$PATH
cd $DA_EXPDIR/$iexpnr

function sedi { if [[ "$OSTYPE" == "darwin"* ]]; then
        		sed -i '' "$1" "$2"
		elif [[ "$OSTYPE" == "linux-gnu" ]]; then
        		sed -i "$1" "$2"
		fi;}

####### set iexpnir in matlab file
sedi "/expnr = '/s/.*/expnr = '$iexpnr';/g" $DA_TOOLSDIR"/write_inputs.m"
###### set # CPUS from execute to test domain size !edit : should maybe multiply by nnode (bss116: yes! this needs to change)
sedi  "/CPUS = /s/.*/CPUS = $(grep -m 1 'ncpu=' ../../u-dales/tools/local_execute.sh | cut -d "=" -f 2 | cut -d " " -f 1 | tr -d ' ');       % # cpus/g" $DA_TOOLSDIR"/write_inputs.m"

###### RUN MATLAB SCRIPT FOR .inp. files

cd $DA_TOOLSDIR/
matlab -nodesktop -nosplash -r "write_inputs; quit"
cd $DA_EXPDIR/$iexpnr

##### alter files in namoptions (thl_top, nblocks etc.)
### Blocks

if [ -f $DA_EXPDIR/$iexpnr/blocks.inp.$iexpnr ]; then
	nblocks=$(wc -l < $DA_EXPDIR/$iexpnr/blocks.inp.$iexpnr)
	nblocks=$(($nblocks-2))
else
	nblocks=0
fi

if grep -q nblocks $DA_EXPDIR/$iexpnr"/namoptions."$iexpnr; then
	sedi "/nblocks/s/.*/nblocks      = $nblocks/g" $DA_EXPDIR/$iexpnr"/namoptions."$iexpnr
else
	sedi '/&WALLS/a\'$'\n''nblocks      = '$nblocks''$'\n' $DA_EXPDIR/$iexpnr"/namoptions."$iexpnr
fi

### Facets

if [ -f $DA_EXPDIR/$iexpnr/blocks.inp.$iexpnr ]; then
	nfcts=$(wc -l < $DA_EXPDIR/$iexpnr/facets.inp.$iexpnr)
	nfcts=$(($nfcts-1))
else
	nfcts=0
fi

if grep -q nfcts $DA_EXPDIR/$iexpnr"/namoptions."$iexpnr; then
	sedi "/nfcts/s/.*/nfcts        = $nfcts/g" $DA_EXPDIR/$iexpnr"/namoptions."$iexpnr
else
	sedi '/&WALLS/a\'$'\n''nfcts        = '$nfcts''$'\n' $DA_EXPDIR/$iexpnr"/namoptions."$iexpnr
fi
