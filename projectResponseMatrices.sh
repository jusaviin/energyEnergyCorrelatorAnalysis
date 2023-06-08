#!/bin/bash

if [ "$#" -ne 2 ]; then
  echo "Usage of the script:"
  echo "$0 [inputFile] [outputFile]"
  echo "inputFile = Name of the input file"
  echo "outputFile = Name of the output file"
  exit
fi

INPUT=$1    # Name of the input file
OUTPUT=$2   # Name of the output file

# Find the git hash of the current commit
GITHASH=`git rev-parse HEAD`

# Replace the placeholder string in the projection code by git hash
sed -i '' 's/GITHASHHERE/'${GITHASH}'/' plotting/projectEEChistograms.C

# Project event information and jet histograms
#root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",3)'

# Project track histograms
#root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",4)'

# Project uncorrected track histograms
#root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",8)'

# Project multiplicity histograms within the jet cone
#root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",16)'

# Project the track density around the jet axis histograms
#root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",32)'

# Project the track pT density around the jet axis histograms
#root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",64)'

# Project the maximum particle pT within the jet cone histograms
#root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",128)'

# Project regular energy-energy correlator histograms
#root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",256)'

# Project energy-energy correlator histograms with positive track efficiency variation
#root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",512)'

# Project energy-energy correlator histograms with negative track efficiency variation
#root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",1024)'

# Project jet pT closure histograms
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",2048)'

# Project jet pT response matrices
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",4096)'

# Put the placeholder string back to the histogram projection file
sed -i '' 's/'${GITHASH}'/GITHASHHERE/' plotting/projectEEChistograms.C
