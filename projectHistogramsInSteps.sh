#!/bin/bash

if [ "$#" -lt 2 ]; then
  echo "Usage of the script:"
  echo "$0 inputFile outputFile [-n] [-e] [-c] [-r] [-o] [-m] [-a] [-u]"
  echo "inputFile = Name of the input file"
  echo "outputFile = Name of the output file"
  echo "-n = Do not project nominal histograms"
  echo "-e = Project energy-energy correlators used for tracking systematics"
  echo "-c = Project jet pT closure histograms"
  echo "-r = Project jet pT response matrices"
  echo "-o = Only project the jet histograms from nominal set of histograms"
  echo "-m = Project track/particle matching study histograms"
  echo "-a = Project histograms for one-dimensional jet pT unfolding"
  echo "-u = Project covariance matrices used in jet pT unfolding"
  exit
fi

INPUT=$1    # Name of the input file
OUTPUT=$2   # Name of the output file
shift 2     # Shift the positional parameters to read the optional ones

# Read the optional arguments. (Semicolon after letter: expects argument)
while getopts ":necromau" opt; do
case $opt in
n) NOMINAL=false
;;
e) ERROR=true
;;
c) CLOSURE=true
;;
r) RESPONSE=true
;;
o) ONLYJETS=true
;;
m) MATCHINGSTUDY=true
;;
a) ONEDIMENSIONALUNFOLD=true
;;
u) COVARIANCE=true
;;
\?) echo "Invalid option -$OPTARG" >&2
exit 1
;;
esac
done

# Set default values to optional arguments if they are not given
NOMINAL=${NOMINAL:-true}
ERROR=${ERROR:-false}
CLOSURE=${CLOSURE:-false}
RESPONSE=${RESPONSE:-false}
MATCHINGSTUDY=${MATCHINGSTUDY:-false}
ONLYJETS=${ONLYJETS:-false}
ONEDIMENSIONALUNFOLD=${ONEDIMENSIONALUNFOLD:-false}
COVARIANCE=${COVARIANCE:-false}

# Find the git hash of the current commit
GITHASH=`git rev-parse HEAD`

# Replace the placeholder string in the projection code by git hash
sed -i '' 's/GITHASHHERE/'${GITHASH}'/' plotting/projectEEChistograms.C

if $ONLYJETS; then

  # Project event information and jet histograms
  root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",3)'

elif $NOMINAL; then

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

fi

if $ERROR; then

  # Project energy-energy correlator histograms with track efficiency variations
  root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",512)'

  # Project energy-energy correlator histograms with track pair efficiency variations
  root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",1024)'

fi

if $CLOSURE; then

  # Project jet pT closure histograms
  root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",2048)'

  # Project jet pT response matrices from closure histograms
  root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",4096)'

fi

if $RESPONSE; then

  # Project jet pT response matrices for unfolding
  root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",8192)'

fi

if $MATCHINGSTUDY; then

  # Project histograms related to track/particle matching study
  root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",16384)'

fi

if $ONEDIMENSIONALUNFOLD; then

  # Project histograms related to one dimensional jet pT unfolding study
  root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",32768)'

fi

if $COVARIANCE; then

  # Project covariiance matrices needed in jet pT unfolding
  root -l -b -q 'plotting/projectEEChistograms.C("'${INPUT}'","'${OUTPUT}'",65536)'

fi

# Put the placeholder string back to the histogram projection file
sed -i '' 's/'${GITHASH}'/GITHASHHERE/' plotting/projectEEChistograms.C
