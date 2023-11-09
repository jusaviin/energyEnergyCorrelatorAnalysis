#!/bin/bash

if [ "$#" -ne 2 ]; then
  echo "Expand the file names such that they can easily be copied to histogram added script"
  echo "Usage of the script:"
  echo "$0 inputFile outputFile"
  echo "inputFile = File containing all the file names that will be expanded"
  echo "outputFile = New file with all the file names in expanded form"
  exit
fi

INPUT=$1
OUTPUT=$2
cp $INPUT $OUTPUT

sed -i -e "s#eWSquared#energyWeightSquared#" ${OUTPUT}
sed -i -e "s#resMat#responseMatrix#" ${OUTPUT}
sed -i -e "s#JECunc_#JECuncertainty_#" ${OUTPUT}
sed -i -e "s#nomSmear#nominalSmear#" ${OUTPUT}
sed -i -e "s#uncSm#uncertaintySm#" ${OUTPUT}
sed -i -e "s#truthRef_#truthReference_#" ${OUTPUT}
sed -i -e "s#recoRef#reconstructedReference#" ${OUTPUT}
sed -i -e "s#jetPtW_#jetPtWeight_#" ${OUTPUT}
sed -i -e "s#4pC_#4pCentShift_#" ${OUTPUT}
sed -i -e "s#2pC_#2pCentShift_#" ${OUTPUT}
sed -i -e "s#6pC_#6pCentShift_#" ${OUTPUT}
