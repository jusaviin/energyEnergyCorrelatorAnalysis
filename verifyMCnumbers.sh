#!/bin/bash

if [ "$#" -ne 1 ]; then
  echo "Script checking if MC files have the correct number of events"
  echo "Usage of the script:"
  echo "$0 inputList"
  echo "inputList = List of all the jobs that were submitted"
  exit
fi

INPUTLIST=$1

# Loop over the input files and run the checker macro for each
while read -r CURRENTFILE
do
  root -l -b -q 'checkNumberOfEventsMC.C("'${CURRENTFILE}'")'
done < <(cat $INPUTLIST)
