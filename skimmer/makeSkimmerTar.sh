#!/bin/bash

# Define the name for the output tar ball
OUTPUTTAR="skimmer.tar.gz"

# If the tar ball already exists, remove it
if [ -f $OUTPUTTAR ]; then 
  rm $OUTPUTTAR 
fi

# Make sure there are no object files going to the tar
make clean

# Create the new tar ball
tar -cvzf megaSkimmerPPb.tar.gz megaSkimmerPPb.cxx Makefile TrackingEfficiencyInterface.cxx TrackingEfficiencyInterface.h trackingEfficiency2016pPb.cxx trackingEfficiency2016pPb.h trackCorrectionTables
