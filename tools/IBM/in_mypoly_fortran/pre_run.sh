#To convert a text file filename.sh to an exicutable bash file filename.sh, use the following command
#	~$ chmod 777 filename.sh


# Usage: ./u-dales/tools/pre_run.sh <PATH_TO_CASE>
# FROM THE PARENT DIRECTORY


set -e

if (( $# < 1 ))
then
    echo "The path to case/experiment folder must be set."
    exit
fi

## go to experiment directory
pushd $1
	inputdir=$(pwd)

	## set experiment number via path
	exp="${inputdir: -3}"

#	echo "Setting up uDALES input files for case $exp..."

	## read in additional variables
	if [ -f config.sh ]; then
   	 source config.sh
	fi

	## check if required variables are set
	if [ -z $DA_TOOLSDIR ]; then
	    echo "Script directory DA_TOOLSDIR must be set"
	    exit
	fi;

popd

pushd $DA_TOOLSDIR/IBM/in_mypoly_fortran/

#	echo "Compiling the Fortran pre-processing code..."
#	gfortran -O2 in_mypoly_functions.f90 IBM_flagging.f90 -o pre.out
	gfortran -O2 -fopenmp in_mypoly_functions.f90 IBM_flagging.f90 -o pre.out
	cp pre.out $inputdir
	rm pre.out in_mypoly_functions.mod
popd

pushd $inputdir

	./pre.out
	rm pre.out inmypoly_inp_info.txt Stl_data.txt vertices.txt zfgrid.txt zhgrid.txt

popd

#cd $DA_TOOLSDIR
#matlab -r "edit my_write_inputs.m"
