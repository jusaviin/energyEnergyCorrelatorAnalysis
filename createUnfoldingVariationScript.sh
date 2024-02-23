#!/bin/bash

if [ "$#" -lt 4 ]; then
  echo "Script to generate another script that does all the unfolding variations for systematic uncertainty analysis."
  echo "Tracking related systematics are not included in this script"
  echo "Usage of the script:"
  echo "$0 inputFile templateString dataType outputScript [-n] [-b] [-e]"
  echo "  inputFile = Name of the input file for which the unfolding is done"
  echo "  templateString = String in the file name that will be replaced to tell which unfolding configuration is used"
  echo "  dataType = PbPb or pp"
  echo "  outputScript = Name given to the unfolding variation script"
  echo "  -n = Include unfolding for nominal results"
  echo "  -b = Add background subtraction for results to the script"
  echo "  -e = Execute the created script"
  exit
fi

# Read the user input
INPUT=$1
TEMPLATE=$2
DATATYPE=$3
OUTPUT=$4
shift 4

# Read the optional arguments. (Semicolon after letter: expects argument)
while getopts ":nbe" opt; do
case $opt in
n) NOMINAL=true
;;
b) BACKGROUND=true
;;
e) EXECUTE=true
;;
\?) echo "Invalid option -$OPTARG" >&2
exit 1
;;
esac
done

# Set default values to optional arguments if they are not given
NOMINAL=${NOMINAL:-false}
BACKGROUND=${BACKGROUND:-false}
EXECUTE=${EXECUTE:-false}

# Enable case insensitive matching
shopt -s nocasematch

# Check that the data type is either PbPb or pp
if [[ $DATATYPE == "PbPb" ]]; then
  ISPBPBDATA=true
elif [[ $DATATYPE == "pp" ]]; then
  ISPBPBDATA=false
else
  echo "Unknown data type: ${DATATYPE}"
  echo "Please give either PbPb or pp"
  exit
fi

# Disable case insensitive matching
shopt -u nocasematch

# Add shebang to the output script
echo "#!/bin/bash" > ${OUTPUT}
echo "" >> ${OUTPUT}

declare -a VARIATIONS
declare -a OPTIONS
declare -a BACKGROUNDOPTIONS

if $NOMINAL; then
  VARIATIONS+=("unfoldingWithNominalSmear")
  OPTIONS+=("-u 0")
  BACKGROUNDOPTIONS+=("-u 0")
fi

# Make arrays of all the configurations that will be done for unfolding
VARIATIONS+=("unfoldingWithUncertaintySmearDown" "unfoldingWithUncertaintySmearUp" "unfoldingWithMinusJetEnergyScale" "unfoldingWithPlusJetEnergyScale" "unfoldingWithModifiedPrior" "unfoldingWith3Iterations" "unfoldingWith5Iterations")
OPTIONS+=("-u 1" "-u 2" "-u 3" "-u 4" "-u 5" "-u 8" "-u 9")
BACKGROUNDOPTIONS+=("-u 0" "-u 0" "-u 0" "-u 0" "-u 0" "-u 0" "-u 0")

if $ISPBPBDATA; then
  VARIATIONS+=("unfoldingWith2pCentShift" "unfoldingWith6pCentShift" "unfoldingWithNominalSmear_backgroundScaleUncertainty2pShift" "unfoldingWithNominalSmear_backgroundScaleUncertainty6pShift" "unfoldingWithNominalSmear_lowSignalToBackgroundScaleEstimateAfterUnfolding" "unfoldingWithNominalSmear_highSignalToBackgroundScaleEstimateAfterUnfolding")
  OPTIONS+=("-u 6" "-u 7" "-u 0" "-u 0" "-u 0" "-u 0")
  BACKGROUNDOPTIONS+=("-u 1" "-u 2" "-u 1" "-u 2" "-u 3" "-u 4")
else
  VARIATIONS+=("unfoldingWithNominalSmear_backgroundSubtractionSystematics")
  OPTIONS+=("-u 0")
  BACKGROUNDOPTIONS+=("-u 1")
fi

# Loop over all the unfolding variations
for i in "${!VARIATIONS[@]}"; do

  # First, create a new file, where the string template part is replaced by a name from variations array
  NEWNAME=$(echo ${INPUT} | sed s/${TEMPLATE}/${VARIATIONS[i]}/g)
  cp $INPUT $NEWNAME

  # Once the new file is created, write the unfolding script with correct options for the new file
  echo "./unfoldEnergyEnergyCorrelators.sh ${NEWNAME} ${OPTIONS[i]}" >> $OUTPUT

  # If we also want to add the background subtraction commands, write them
  if $BACKGROUND; then
    echo "./processUnfoldedEnergyEnergyCorrelators.sh ${NEWNAME} ${BACKGROUNDOPTIONS[i]}" >> $OUTPUT
  fi
done

# Make the output script file executable
chmod +x $OUTPUT

# If selected, execute the created script
if $EXECUTE; then
  ./${OUTPUT}
  rm $OUTPUT
fi
