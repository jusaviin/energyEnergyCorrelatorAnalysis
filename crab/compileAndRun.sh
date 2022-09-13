#!/bin/bash
echo "This is job number $1"

# Read the script arguments in format name=value
# The line below only reads the value from the above format
CARD=${2#*=}
OUTPUT=${3#*=}
LOCATION=${4#*=}

# Untar the input file list
tar xf input_files.tar.gz

# Unzip tar ball
tar -xvzf eec5TeV.tar.gz

# Compile the code
make

# Run the code
./eecAnalysis $1 $CARD $OUTPUT $LOCATION
