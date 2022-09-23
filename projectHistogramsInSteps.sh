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

# Project event information, track, and jet histograms
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",15)'

# Project regular energy-energy correlator histograms
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",16)'

# Project jet pT weighted energy-energy correlator histograms
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",32)'

# Project uncorrected energy-energy correlator histograms
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",64)'

# Project uncorrected jet pT weighted energy-energy correlator histograms
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",128)'

# Project jet pT closure histograms
# root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",256)'
