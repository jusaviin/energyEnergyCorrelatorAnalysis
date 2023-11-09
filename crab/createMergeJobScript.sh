#!/bin/bash

if [ "$#" -ne 2 ]; then
  echo "Create a merge script that can merge all the jobs defined in input list"
  echo "Usage of the script:"
  echo "$0 inputList mergeScriptName"
  echo "inputList = List of job names that need to be merged"
  echo "mergeScriptName = Name given to the merging script"
  exit
fi

INPUT=$1
MERGESCRIPT=$2

# If the script exist already, delete it
if test -f ${MERGESCRIPT}; then
  rm ${MERGESCRIPT}
fi

# Expand the names for the input file to temporary file
TEMPFILE="tempExpandedFile.txt"
./expandFile.sh $INPUT $TEMPFILE

# From each line in the file, add a merge command to the merge script
exec 3< ${INPUT}    # Assign file descriptor 3 to the input file
exec 4< ${TEMPFILE} # Assign file descriptor 4 to the temporary file
while read -r LINE <&3 && read -r EXPANDEDLINE <&4; do
  EOSFOLDER=(`./findEosFolder.sh ${LINE}`)
  NFILES=(`eos root://cmseos.fnal.gov ls ${EOSFOLDER} | wc -l`)
  echo "./addHistogramsInSteps.sh \"${EOSFOLDER}\" \"${EXPANDEDLINE}\" ${NFILES} 2" >> ${MERGESCRIPT}
done < ${INPUT} < ${TEMPFILE}
exec 3<&-  # Close file descriptor 3
exec 4<&-  # Close file descriptor 4

# Delete the temporary expanded file
rm $TEMPFILE

# Make the created script executable
chmod +x ${MERGESCRIPT}
