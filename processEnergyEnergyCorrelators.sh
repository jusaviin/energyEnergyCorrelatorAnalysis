#!/bin/bash

if [ "$#" -lt 1 ]; then
  echo "Usage of the script:"
  echo "$0 fileName [-o outputFileName] [-u systematicUncertainty]"
  echo "fileName = Name of the file for which energy-energy correlator histograms are processed"
  echo "-o outputFileName = If given, instead of updating the file fileName, a new file called outputFileName is created with the processed histograms"
  echo "-u systematicUncertainty = Index for systematic uncertainty. 0 = Nominal result, 1 = Background scale from 2% centrality shift, 2 = Background scale from 6% centrality shift"
  exit
fi

FILENAME=$1    # Name of the files for processing
shift # Shift the positional argument such that optional ones are properly seen afterwards

# Read the optional arguments
while getopts ":o:u:" opt; do
case $opt in
o) OUTPUTFILE="$OPTARG"
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
SYSTEMATIC=${SYSTEMATIC:-0}

# Find the git hash of the current commit
GITHASH=`git rev-parse HEAD`

# Replace the placeholder string in the projection code by git hash
sed -i '' 's/GITHASHHERE/'${GITHASH}'/' plotting/processEEChistograms.C

# Process the energy-energy correlator histograms
root -l -b -q 'plotting/processEEChistograms.C("'${FILENAME}'","'${OUTPUTFILE}'",'${SYSTEMATIC}')'

# Put the placeholder string back to the histogram projection file
sed -i '' 's/'${GITHASH}'/GITHASHHERE/' plotting/processEEChistograms.C
