#!/bin/bash

if [ "$#" -ne 1 ]; then
  echo "Usage of the script:"
  echo "$0 [inputFile]"
  echo "fileName = Name of the file for which energy-energy correlator histograms are processed"
  exit
fi

FILENAME=$1    # Name of the files for processing

# Find the git hash of the current commit
GITHASH=`git rev-parse HEAD`

# Replace the placeholder string in the projection code by git hash
sed -i '' 's/GITHASHHERE/'${GITHASH}'/' plotting/processEEChistograms.C

# Process the energy-energy correlator histograms
root -l -b -q 'plotting/processEEChistograms.C("'${FILENAME}'")'

# Put the placeholder string back to the histogram projection file
sed -i '' 's/'${GITHASH}'/GITHASHHERE/' plotting/processEEChistograms.C
