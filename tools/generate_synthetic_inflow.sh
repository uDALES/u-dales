## To convert a text file filename.sh to an exicutable bash file filename.sh, use the following command
##	~$ chmod 777 filename.sh


## Usage: ./u-dales/tools/generate_synthetic_inflow.sh <PATH_TO_CASE>
## FROM THE PARENT DIRECTORY

#!/usr/bin/env bash

set -e

if (( $# < 1 )) 
then
    echo "The path to case/experiment folder must be set."
	echo "usage: FROM THE PARENT DIRECTORY ./u-dales/tools/generate_synthetic_inflow.sh <PATH_TO_CASE> "
	echo "... execution terminated"
    exit 0
fi

## go to experiment directory
pushd $1
inputdir=$(pwd)

## set experiment number via path
exp="${inputdir: -3}"

echo "Creating driver files for uDALES target simulation case $exp..."

## read in additional variables
if [ -f config.sh ]; then
    source config.sh
fi;

popd

## check if required variables are set
if [ -z $DA_TOOLSDIR ]; then
    echo "Script directory DA_TOOLSDIR must be set"
    exit
fi;
if [ -z $DA_EXPDIR ]; then
    echo "Experiment directory DA_EXPDIR must be set"
    exit
fi;
if [ -z $DA_WORKDIR ]; then
    echo "Work directory DA_WORKDIR must be set"
    exit
fi;

if [ -d $DA_WORKDIR/$exp ]; then
    echo "target directory $exp already exists in $DA_WORKDIR"
else
    mkdir -p $DA_WORKDIR/$exp
    echo "Creating target work directory, $exp, in $DA_WORKDIR"
fi;

pushd $DA_WORKDIR/$exp


# check if driver files already exists. If so, ask how to proceed.
if [ -f 'tdriver_000.'$exp ]; then
  echo "Few or all driver files *driver_*.$exp already exists inside $DA_WORKDIR/$exp/"
  echo "continue with overwrite/abort? (o/a)"
  read answer
  if [ "$answer" == "o" ]; then
    case=1
  else
    echo "abort"
    exit 1
  fi
else
  case=1
fi

if [ case==1 ]; then

	## compiling synthetic inflow generator source code
	gfortran -O2 -fopenmp $DA_TOOLSDIR/syntheticInflow/modSyntheticInflow.f90 -o synInflow_executable

	## copying syntheticInfow inputs to the work directory   (may be modified later)
	cp $DA_EXPDIR/$exp/syntheticInflow_inputs/*.txt .
	cp $DA_EXPDIR/$exp/*'.'$exp .

	## setting up OpenMP environmental variables
	## export OMP_NUM_THREADS=8    ## please increase this for large job based on availability of cores
	export OMP_PLACES=cores

	## run
	./synInflow_executable $exp 2>&1 | tee -a log_synthetic_inflow_computation.$exp

	## cleaning work dircetory
	rm -f synInflow_executable
	rm -f Reynolds_stress_profiles_velocity.txt
	rm -f Reynolds_stress_profiles_temp.txt
	rm -f Reynolds_stress_profiles_moist.txt
	rm -f length_time_scales_u.txt length_time_scales_v.txt length_time_scales_w.txt
	rm -f length_time_scales_temp.txt length_time_scales_qt.txt
	
fi;

popd
