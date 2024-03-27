#!/bin/bash

if [ "$#" -lt 1 ]; then
  echo "Usage of the script:"
  echo "$0 fileName [-o outputFileName] [-b backgroundMethod] [-u systematicUncertainty] [-e energyEnergyCorrelatorType]"
  echo "fileName = Name of the file for which energy-energy correlator histograms are processed"
  echo "-o outputFileName = If given, instead of updating the file fileName, a new file called outputFileName is created with the processed histograms"
  echo "-b backgroundMethod = Index for used background estimation method. Default = 0."
  echo "  0 = Mixed cone background"
  echo "  1 = Reflected cone background" 
  echo "-u systematicUncertainty = Index for systematic uncertainty. Default = 0."
  echo "  0 = Nominal result"
  echo "  1 = Background scale from 2% centrality shift" 
  echo "  2 = Background scale from 6% centrality shift"
  echo "  3 = Lower scaling estimate for signal-to-background ratio after unfolding"
  echo "  4 = Higher scaling estimate for signal-to-background ratio after unfolding"
  echo "-e energyEnergyCorrelatorType = Index for energy-energy correlators. Default: 0."
  echo "   0 = Regular energy-energy correlator"
  echo "   1 = Energy-energy correlator with positive track efficiency variation"
  echo "   2 = Energy-energy correlator with negative track efficiency variation"
  echo "   3 = Energy-energy correlator with positive track pair efficiency variation"
  echo "   4 = Energy-energy correlator with negative track pair efficiency variation"
  exit
fi

FILENAME=$1    # Name of the files for processing
shift # Shift the positional argument such that optional ones are properly seen afterwards

# Read the optional arguments
while getopts ":o:b:u:e:" opt; do
case $opt in
o) OUTPUTFILE="$OPTARG"
;;
b) BACKGROUND="$OPTARG"
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
BACKGROUND=${BACKGROUND:-0}
SYSTEMATIC=${SYSTEMATIC:-0}
ENERGYENERGYCORRELATOR=${ENERGYENERGYCORRELATOR:-0}

# Find the git hash of the current commit
GITHASH=`git rev-parse HEAD`

# Replace the placeholder string in the processing code by git hash
sed -i '' 's/GITHASHHERE/'${GITHASH}'/' plotting/processUnfoldedEEChistograms.C

# Process the unfolded energy-energy correlator histograms
root -l -b -q 'plotting/processUnfoldedEEChistograms.C("'${FILENAME}'","'${OUTPUTFILE}'",'${BACKGROUND}','${SYSTEMATIC}','${ENERGYENERGYCORRELATOR}')'

# Put the placeholder string back to the histogram processing file
sed -i '' 's/'${GITHASH}'/GITHASHHERE/' plotting/processUnfoldedEEChistograms.C
