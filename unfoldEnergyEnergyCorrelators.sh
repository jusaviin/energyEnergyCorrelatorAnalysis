#!/bin/bash

if [ "$#" -lt 1 ]; then
  echo "Usage of the script:"
  echo "$0 fileName [-o outputFileName] [-s split] [-u systematicUncertainty] [-e energyEnergyCorrelatorType]"
  echo "fileName = Name of the file containing the energy-energy correlators that needs unfolding"
  echo "-o outputFileName = If given, instead of updating the file fileName, a new file called outputFileName is created with the unfolded histograms"
  echo "-s split = Split index for response matrix. Default: 0."
  echo "   0 = Full MC statistics"
  echo "   1 = First half of statistics"
  echo "   2 = Second half of statistics"
  echo "-u systematicUncertainty = Index for systematic uncertainty. Default: 0."
  echo "   0 = Nominal result"
  echo "   1 = Jet energy resolution uncertainty"
  echo "   2 = Jet energy scale uncertainty"
  echo "-e energyEnergyCorrelatorType = Index for energy-energy correlators. Default: 0."
  echo "   0 = Regular energy-energy correlator"
  echo "   1 = Energy-energy correlator with positive track efficiency variation"
  echo "   2 = Energy-energy correlator with negative track efficiency variation"
  echo "   3 = Energy-energy correlator with positive track pair efficiency variation"
  echo "   4 = Energy-energy correlator with negative track pair efficiency variation"
  exit
fi

# The first argument is the mandatory input file
FILENAME=$1
shift # Shift the positional argument such that optional ones are properly seen afterwards

# Read the optional arguments
while getopts ":o:s:u:e:" opt; do
case $opt in
o) OUTPUTFILE="$OPTARG"
;;
s) SPLIT="$OPTARG"
;;
u) SYSTEMATIC="$OPTARG"
;;
e) ENERGYENERGYCORRELATOR="$OPTARG"
;;
\?) echo "Invalid option -$OPTARG" >&2
exit 1
;;
esac
done

# Set default values to optional arguments if they are not given
OUTPUTFILE=${OUTPUTFILE:-$FILENAME}
SPLIT=${SPLIT:-0}
SYSTEMATIC=${SYSTEMATIC:-0}
ENERGYENERGYCORRELATOR=${ENERGYENERGYCORRELATOR:-0}

# Find the git hash of the current commit
GITHASH=`git rev-parse HEAD`

# Replace the placeholder string in the unfolding code by git hash
sed -i '' 's/GITHASHHERE/'${GITHASH}'/' plotting/unfoldEEChistograms.C

# Unfold the energy-energy correlator histograms
root -l -b -q 'plotting/unfoldEEChistograms.C("'${FILENAME}'","'${OUTPUTFILE}'",'${SPLIT}','${SYSTEMATIC}','${ENERGYENERGYCORRELATOR}')'

# Put the placeholder string back to the histogram unfolding file
sed -i '' 's/'${GITHASH}'/GITHASHHERE/' plotting/unfoldEEChistograms.C
