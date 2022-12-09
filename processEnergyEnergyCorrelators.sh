#!/bin/bash

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
  echo "Usage of the script:"
  echo "$0 [fileName] <outputFileName>"
  echo "fileName = Name of the file for which energy-energy correlator histograms are processed"
  echo "outputFileName = If given, instead of updating the file fileName, a new file called outputFileName is created with the processed histograms"
  exit
fi

FILENAME=$1    # Name of the files for processing

# if second argument is given, use that as output file. Otherwise update input file.
if [ "$#" -gt 1 ]; then
  OUTPUTFILE=$2
else
  OUTPUTFILE=$FILENAME
fi

# Find the git hash of the current commit
GITHASH=`git rev-parse HEAD`

# Replace the placeholder string in the projection code by git hash
sed -i '' 's/GITHASHHERE/'${GITHASH}'/' plotting/processEEChistograms.C

# Process the energy-energy correlator histograms
root -l -b -q 'plotting/processEEChistograms.C("'${FILENAME}'","'${OUTPUTFILE}'")'

# Put the placeholder string back to the histogram projection file
sed -i '' 's/'${GITHASH}'/GITHASHHERE/' plotting/processEEChistograms.C
