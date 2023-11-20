#!/bin/bash

if [ "$#" -lt 2 ]; then
  echo "Script to generate another script to projects histograms from THnSparses"
  echo "Usage of the script:"
  echo "$0 inputList outputScript [-e]"
  echo "inputList = List of all the files for which the projection is done"
  echo "outputScript = Name given to the histogram projection script"
  echo "-e = Execute the script after it is created"
  exit
fi

# Read the used input
INPUT=$1
OUTPUT=$2
shift 2

# Read the optional arguments. (Semicolon after letter: expects argument)
while getopts ":e" opt; do
case $opt in
e) EXECUTE=true
;;
\?) echo "Invalid option -$OPTARG" >&2
exit 1
;;
esac
done

# Set default values to optional arguments if they are not given
EXECUTE=${EXECUTE:-false}

# Add shebang to the output script
echo "#!/bin/bash" > ${OUTPUT}
echo "" >> ${OUTPUT}

# Loop over all the files defined in the input file list
while read -r LINE; do

  # Define a string for options to give to histogram projecting 
  OPTIONS=""

  # First we need to determine from the file name which histograms to project
  # If we are dealing with response matrices, project those
  if [[ $LINE == *responseMatrix* ]]; then
    OPTIONS+=" -r -o"
  fi

  # Extract the date part from the file name
  DATE=$(echo "$LINE" | grep -o -E '\d{4}-\d{2}-\d{2}')

  # Add the text "processed_" before the date
  PROCESSEDLINE=$(echo "$LINE" | sed 's'/${DATE}'/processed_'${DATE}'/g')

  # Print the command to the output script
  echo "./projectHistogramsInSteps.sh ${LINE} ${PROCESSEDLINE}${OPTIONS}" >> ${OUTPUT}
done < $INPUT

# Make the output script file executable
chmod +x $OUTPUT

# If selected, execute the created script
if $EXECUTE; then
  ./${OUTPUT}
  rm $OUTPUT
fi
