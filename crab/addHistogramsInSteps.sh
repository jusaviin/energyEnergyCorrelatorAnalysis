#!/bin/bash

if [ "$#" -ne 4 ]; then
  echo "Usage of the script:"
  echo "$0 [folderName] [baseName] [nFiles] [nSteps]"
  echo "folderName = The name of the folder where the files are located in EOS"
  echo "baseName = Output file name without .root extension"
  echo "nFiles = Number of files"
  echo "nSteps = Number of steps used to merge the files"
  exit
fi

FOLDERPATH=${1%/}
BASENAME=$2
NFILES=$3
NSTEPS=$4

# Determine the ranges of added files
declare -a FILELOWINDEX
declare -a FILEHIGHINDEX
RUNNINGNUMBER=0
FILESPERSTEP=$((${NFILES}/${NSTEPS}))
MODULO=$((${NFILES}%${NSTEPS}))

for i in `seq 1 $NSTEPS`;
do
    FILELOWINDEX[$i-1]=$((${RUNNINGNUMBER}+1))
    RUNNINGNUMBER=$((${RUNNINGNUMBER}+${FILESPERSTEP}))
    FILEHIGHINDEX[$i-1]=$RUNNINGNUMBER
done

# Add the reminder files to the last merge job
FILEHIGHINDEX[${NSTEPS}-1]=$((${FILEHIGHINDEX[${NSTEPS}-1]}+${MODULO}))

# Merge the histograms in the defined number of steps
declare -a TEMPFILES
for i in `seq 1 $NSTEPS`;
do
    TEMPFILES[$i-1]="${BASENAME}_part${i}.root"
    hadd -ff ${TEMPFILES[$i-1]} `xrdfs root://cmseos.fnal.gov ls -u $FOLDERPATH | grep '\.root' | sed -n ${FILELOWINDEX[$i-1]},${FILEHIGHINDEX[$i-1]}p`
done

hadd -ff ${BASENAME}.root ${TEMPFILES[*]}

# After the histograms are merged, run the JCard skimmer
root -l -b -q 'skimJCardEEC.C("'${BASENAME}'.root",'${NFILES}')'

# In the end, delete the temporary files that were created in the intermediate step
for filename in "${TEMPFILES[@]}";
do
    rm $filename
done
