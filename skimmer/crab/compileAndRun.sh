#!/bin/bash
echo "This is job number $1"

# Read the script arguments in format name=value
# The line below only reads the value from the above format
ISMC=${2#*=}
OUTPUT=${3#*=}
LOCATION=${4#*=}

# Untar the input file list
tar xf input_files.tar.gz

# Unzip tar ball
tar -xvzf skimmer.tar.gz

# Compile the code
make

# Run the code
./skimmer $1 $ISMC $OUTPUT $LOCATION
