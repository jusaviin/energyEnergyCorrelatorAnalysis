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

# Project event information and jet histograms
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",3)'

# Project track histograms
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",4)'

# Project uncorrected track histograms
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",8)'

# Project multiplicity histograms within the jet cone
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",16)'

# Preject the track density around the jet axis histograms
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",32)'

# Preject the track pT density around the jet axis histograms
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",64)'

# Project regular energy-energy correlator histograms
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",128)'

# Project jet pT weighted energy-energy correlator histograms
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",256)'

# Project uncorrected energy-energy correlator histograms
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",512)'

# Project uncorrected jet pT weighted energy-energy correlator histograms
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",1024)'

# Project jet pT closure histograms
# root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",128)'
