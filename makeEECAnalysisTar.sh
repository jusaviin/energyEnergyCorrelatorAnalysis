#!/bin/bash

# Read the optional arguments
while getopts ":s" opt; do
case $opt in
s) SMEARING=true
;;
\?) echo "Invalid option -$OPTARG" >&2
exit 1
;;
esac
done

# Set the default value if optional argument is not given
SMEARING=${SMEARING:-false}

# Find the git hash of the current commit
GITHASH=`git rev-parse HEAD`

# Define the name for the output tar ball
OUTPUTTAR="eec5TeV.tar.gz"

# If the tar ball already exists, remove it
if [ -f $OUTPUTTAR ]; then 
rm $OUTPUTTAR 
fi

# Replace the placeholder string in the main analysis file by git hash
sed -i '' 's/GITHASHHERE/'${GITHASH}'/' eecAnalysis.cxx

# Make sure there are no object files going to the tar
make clean

# Create the new tar ball
if $SMEARING; then
  tar -cvzf $OUTPUTTAR Makefile eecAnalysis.cxx jetEnergyCorrections src trackCorrectionTables smearingFiles mixingFileList
else
  tar -cvzf $OUTPUTTAR Makefile eecAnalysis.cxx jetEnergyCorrections src trackCorrectionTables mixingFileList
fi

# Put placeholder string back to the main analysis file
sed -i '' 's/'${GITHASH}'/GITHASHHERE/' eecAnalysis.cxx
