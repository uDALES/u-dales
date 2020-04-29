#!/bin/bash
#examples:
#   new sim: "da_prep 2 1" or "da_prep 2 1 c"
#   continue with old sim: "da_prep 1 1 w"
#   continue with new sim  "da_prep 2 1 w"

set -xe

if [ -z $DA_EXPDIR_SRC ]; then
  DA_EXPDIR_SRC=$DA_EXPDIR
fi;
if [ -z $DA_WORKDIR_SRC ]; then
  DA_WORKDIR_SRC=$DA_WORKDIR
fi;

if (( $# < 2 ))
then
 echo "usage: `basename $0` sim#1 sim#2 (start)"
 echo "prepares a new simulation with DALES"
 echo "   sim#1: number of new simulation"
 echo "   sim#2: number of simulation upon which the newone is based"
 echo "   start (optional): (c)old- or (w)arm-start, c is default"
 echo "... execution terminated"
 exit 0
fi

tar=$(printf "%03.0f" $1)    #pad argument 1 (target simulation) with zeros
src=$(printf "%03.0f" $2)    #pad argument 2 (origin) with zeros
start=${3:-c}                #default argument 3 to c (coldstart) if not provided
case=1                       #default value 1 (target folder does not exist or overwrite, 2=only copy files that don't exist, 3=interactive)


if [ -d $DA_EXPDIR/$tar ]; then  #check if target simulation already exists. If so, ask how to proceed.
  echo "target directory $tar already exists in $DA_EXPDIR"
  echo "continue with overwrite/copy nonexistent files/interactive/abort? (o/c/i/a)"
  read answer
  if [ "$answer" == "o" ]; then
  case=1
  elif [ "$answer" == "c" ]; then
  case=2
  elif [ "$answer" == "i" ]; then
  case=3
  else
  echo "abort"
  exit 1
  fi
else
mkdir -p $DA_EXPDIR/$tar
fi
echo "continuing"

if [ -d $DA_WORKDIR/$tar ]; then
  echo "target directory $tar already exists in $DA_WORKDIR"
else
  mkdir -p $DA_WORKDIR/$tar
  echo "Creating target work directory, $tar, in $DA_WORKDIR"
fi

if [ ! -d $DA_EXPDIR_SRC/$src ]; then  #test if original simulation exists
  echo "original simulation $src does not exist in $DA_EXPDIR"
  echo "exit"
  exit 1
fi




if [ $start == "c" ]; then   #coldstart: copy and rename files. change simulation number in namoptions
#to copy
declare -a tocopy=("/namoptions." "/lscale.inp." "/prof.inp." "/scalar.inp." "/xgrid.inp." "/zgrid.inp." "/blocks.inp." "/purifs.inp." "/trees.inp." "/scals.inp." "/lad.inp." "/facetarea.inp." "/facetnumbers.inp." "/facets.inp." "/netsw.inp." "/pf1sf2.inp." "/svf.inp." "/Tfacinit.inp." "/Tfacinitnudged.inp." "/vf.inp." "/walltypes.inp.")

case $case in
     1)  #target folder didn't exist or just overwrite all files
            for i in "${tocopy[@]}"
            do
		cp $DA_EXPDIR_SRC/$src$i$src $DA_EXPDIR/$tar$i$tar 2>/dev/null || :
		#echo $i
            done
          ;; 
     2)  #don't overwrite files in target folder #cp -n does not overwrite an existing file
            for i in "${tocopy[@]}"
            do
		cp -n $DA_EXPDIR_SRC/$src$i$src $DA_EXPDIR/$tar$i$tar 2>/dev/null || :
		#echo $i
            done
          ;; 
     3)  #ask for every file if to overwrite or not
            for i in "${tocopy[@]}"
            do
		cp -i $DA_EXPDIR_SRC/$src$i$src $DA_EXPDIR/$tar$i$tar 2>/dev/null || :
		#echo $i
            done
          ;; 
esac
sed -i.bak -e '/lwarmstart/s/.*/lwarmstart =  .false./g' $DA_EXPDIR/$tar"/namoptions."$tar  #set warmstart to false in namoptions
sed -i.bak -e "/iexpnr/s/.*/iexpnr      = $tar/g" $DA_EXPDIR/$tar"/namoptions."$tar  #change experiment number in namoptions

rm $DA_EXPDIR/$tar"/namoptions."$tar".bak"	
	
	
	
	
elif [ $start == "w" ]; then  #copy and rename files,
declare -a tocopy=("/namoptions." "/lscale.inp." "/prof.inp." "/scalar.inp." "/xgrid.inp." "/zgrid.inp." "/blocks.inp." "/job." "/purifs.inp." "/trees.inp." "/scals.inp." "/lad.inp." "/facetarea.inp." "/facetnumbers.inp." "/facets.inp." "/netsw.inp." "/pf1sf2.inp." "/svf.inp." "/Tfacinit.inp." "/Tfacinitnudged.inp." "/vf.inp." "/walltypes.inp.")
case $case in
     1)  #target folder didn't exist or just overwrite all files
            for i in "${tocopy[@]}"
            do
		cp $DA_EXPDIR_SRC/$src$i$src $DA_EXPDIR/$tar$i$tar 2>/dev/null || :
		#echo $i
            done
          ;; 
     2)  #don't overwrite files in target folder #cp -n does not overwrite an existing file
            for i in "${tocopy[@]}"
            do
		cp -n $DA_EXPDIR_SRC/$src$i$src $DA_EXPDIR/$tar$i$tar 2>/dev/null || :
		#echo $i
            done
          ;; 
     3)  #ask for every file if to overwrite or not
            for i in "${tocopy[@]}"
            do
		cp -i $DA_EXPDIR_SRC/$src$i$src $DA_EXPDIR/$tar$i$tar 2>/dev/null || :
		#echo $i
            done
          ;; 
esac

### link newest warmstart files
startfilen=$(ls -t $DA_WORKDIR_SRC/$src"/initd"* | head -1)
if [ -z "$startfilen" ]; then
echo "no restart files found in $DA_WORKDIR_SRC/$src"
echo "exit"
exit 1
fi
startfilen=${startfilen##*/}  # retain the part after the last slash
startfilen=${startfilen%_*}   # retain the part before the underscore
ln -s $DA_WORKDIR_SRC/$src"/"*$startfilen* $DA_EXPDIR/$tar/
for f in $DA_EXPDIR/$tar/*.$src; do 
mv $f "${f%.$src}.$tar"
done
# we do not copy files because it promtps errors in mac. a link is enough.
# cp $DA_WORKDIR_SRC/$src"/"*$startfilen* $DA_WORKDIR/$tar
# for f in $DA_WORKDIR/$tar/*.$src; do 
# mv $f "${f%.$src}.$tar"
# done
# ln -s $DA_WORKDIR/$tar"/"*$startfilen* $DA_EXPDIR/$tar/

sed -i.bak -e '/lwarmstart/s/.*/lwarmstart =  .true./g' $DA_EXPDIR/$tar"/namoptions."$tar #set warmstart to true in namoptions
sed -i.bak -e "/iexpnr/s/.*/iexpnr      = $tar/g" $DA_EXPDIR/$tar"/namoptions."$tar  #change experiment number in namoptions
sed -i.bak -e "/startfile/s/.*/startfile  =  '$startfilen\_xxx.$tar'/g" $DA_EXPDIR/$tar"/namoptions."$tar # change startfile to newest restartfiles

rm $DA_EXPDIR/$tar"/namoptions."$tar".bak"

fi

## copy config script for execution
cp $DA_EXPDIR_SRC/$src/config.sh $DA_EXPDIR/$tar
