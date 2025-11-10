#!/bin/bash

if [ "$#" -lt 1 ]; then
  echo "Usage of the script:"
  echo "$0 fileName [-o outputFileName] [-b backgroundSubtraction] [-u systematicUncertainty]"
  echo "fileName = Name of the file for which energy-energy correlator histograms are processed"
  echo "-o outputFileName = If given, instead of updating the file fileName, a new file called outputFileName is created with the processed histograms"
  echo "-b backgroundMethod = Index for used background estimation method. Default = 0."
  echo "  0 = Mixed cone background"
  echo "  1 = Perpendicular cone background"
  echo "  2 = Reflected cone background"
  echo "  3 = No background subtraction"
  echo "  4 = Truth level separation for signal and background in Monte Carlo"
  echo "-u systematicUncertainty = Index for systematic uncertainty. Default = 0."
  echo "  0 = Nominal result"
  echo "  1 = Background scale from 2% centrality shift"
  echo "  2 = Background scale from 6% centrality shift"
  exit
fi

FILENAME=$1    # Name of the files for processing
shift # Shift the positional argument such that optional ones are properly seen afterwards

# Read the optional arguments
while getopts ":o:b:u:" opt; do
case $opt in
o) OUTPUTFILE="$OPTARG"
;;
b) BACKGROUND="$OPTARG"
;;
u) SYSTEMATIC="$OPTARG"
;;
\?) echo "Invalid option -$OPTARG" >&2
exit 1
;;
esac
done

# Set default values to optional arguments if they are not given
OUTPUTFILE=${OUTPUTFILE:-$FILENAME}
BACKGROUND=${BACKGROUND:-0}
SYSTEMATIC=${SYSTEMATIC:-0}

# Find the git hash of the current commit
GITHASH=`git rev-parse HEAD`

# Check which operating system we are using
# The sed command takes different arguments for Mac and Linux, so we need to adjust that accordingly
OS=$(uname)

# Replace the placeholder string in the projection code by git hash
if [ "$OS" == "Darwin" ]; then
  # For Mac, we specify that no backup file is needed with the argument ''
  sed -i '' 's/GITHASHHERE/'${GITHASH}'/' plotting/processEEChistograms.C
else
  # For Linux, '' is not a valid argument, so it needs to be removed from command
  sed -i 's/GITHASHHERE/'${GITHASH}'/' plotting/processEEChistograms.C
fi

# Process the energy-energy correlator histograms
root -l -b -q 'plotting/processEEChistograms.C("'${FILENAME}'","'${OUTPUTFILE}'",'${BACKGROUND}','${SYSTEMATIC}')'

# Put the placeholder string back to the histogram projection file
if [ "$OS" == "Darwin" ]; then
  # For Mac, we specify that no backup file is needed with the argument ''
  sed -i '' 's/'${GITHASH}'/GITHASHHERE/' plotting/processEEChistograms.C
else
  # For Linux, '' is not a valid argument, so it needs to be removed from command
  sed -i 's/'${GITHASH}'/GITHASHHERE/' plotting/processEEChistograms.C
fi
