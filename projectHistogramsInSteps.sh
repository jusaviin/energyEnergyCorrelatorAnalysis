#!/bin/bash

if [ "$#" -lt 2 ]; then
  echo "Usage of the script:"
  echo "$0 inputFile outputFile [-e]"
  echo "inputFile = Name of the input file"
  echo "outputFile = Name of the output file"
  echo "-e = Project also energy-energy correlators used for tracking systematics"
  exit
fi

INPUT=$1    # Name of the input file
OUTPUT=$2   # Name of the output file
shift 2     # Shift the positional parameters to read the optional ones

# Read the optional arguments. (Semicolon after letter: expects argument)
while getopts ":e" opt; do
case $opt in
e) ERROR=true
;;
\?) echo "Invalid option -$OPTARG" >&2
exit 1
;;
esac
done

# Set default values to optional arguments if they are not given
ERROR=${ERROR:-false}

# Find the git hash of the current commit
GITHASH=`git rev-parse HEAD`

# Replace the placeholder string in the projection code by git hash
sed -i '' 's/GITHASHHERE/'${GITHASH}'/' plotting/projectEEChistograms.C

# Project event information and jet histograms
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",3)'

# Project track histograms
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",4)'

# Project uncorrected track histograms
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",8)'

# Project multiplicity histograms within the jet cone
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",16)'

# Project the track density around the jet axis histograms
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",32)'

# Project the track pT density around the jet axis histograms
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",64)'

# Project the maximum particle pT within the jet cone histograms
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",128)'

# Project regular energy-energy correlator histograms
root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",256)'

if $ERROR; then

  # Project energy-energy correlator histograms with track efficiency variations
  root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",512)'

  # Project energy-energy correlator histograms with track pair efficiency variations
  root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",1024)'

fi

# Project jet pT closure histograms
#root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",2048)'

# Project jet pT response matrices
#root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",4096)'

# Put the placeholder string back to the histogram projection file
sed -i '' 's/'${GITHASH}'/GITHASHHERE/' plotting/projectEEChistograms.C
