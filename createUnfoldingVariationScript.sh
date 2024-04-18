#!/bin/bash

if [ "$#" -lt 4 ]; then
  echo "Script to generate another script that does all the unfolding variations for systematic uncertainty analysis."
  echo "Tracking related systematics are not included in this script"
  echo "Usage of the script:"
  echo "$0 inputFile templateString dataType outputScript [-n] [-b backgroundTemplate] [-e]"
  echo "  inputFile = Name of the input file for which the unfolding is done"
  echo "  templateString = String in the file name that will be replaced to tell which unfolding configuration is used"
  echo "  dataType = PbPb or pp"
  echo "  outputScript = Name given to the unfolding variation script"
  echo "  -n = Include unfolding for nominal results"
  echo "  -b backgroundTemplate = Background template string that is replaced with the actual background selection. Default: allBackgrounds"
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
while getopts ":nb:e" opt; do
case $opt in
n) NOMINAL=true
;;
b) BACKGROUNDTEMPLATE="$OPTARG"
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
BACKGROUNDTEMPLATE=${BACKGROUNDTEMPLATE:-"allBackgrounds"}
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
declare -a BACKGROUNDSTRING

if $NOMINAL; then
  VARIATIONS+=("unfoldingWithNominalSmear")
  OPTIONS+=("-u 0")
  BACKGROUNDOPTIONS+=("-u 0")
  if $ISPBPBDATA; then
    BACKGROUNDSTRING+=("mixedConeBackground")
  else
    BACKGROUNDSTRING+=("reflectedConeBackground")
  fi
fi

# Make arrays of all the configurations that will be done for unfolding
VARIATIONS+=("unfoldingWithUncertaintySmearDown" "unfoldingWithUncertaintySmearUp" "unfoldingWithMinusJetEnergyScale" "unfoldingWithPlusJetEnergyScale" "unfoldingWithModifiedPrior" "unfoldingWith3Iterations" "unfoldingWith5Iterations")
OPTIONS+=("-u 1" "-u 2" "-u 3" "-u 4" "-u 5" "-u 8" "-u 9")
BACKGROUNDOPTIONS+=("-u 0" "-u 0" "-u 0" "-u 0" "-u 0" "-u 0" "-u 0")
BACKGROUNDSTRING+=("mixedConeBackground" "mixedConeBackground" "mixedConeBackground" "mixedConeBackground" "mixedConeBackground" "mixedConeBackground" "mixedConeBackground")

if $ISPBPBDATA; then
  VARIATIONS+=("unfoldingWith2pCentShift" "unfoldingWith6pCentShift" "unfoldingWithNominalSmear" "unfoldingWithNominalSmear_lowSignalToBackgroundScaleEstimateAfterUnfolding" "unfoldingWithNominalSmear_highSignalToBackgroundScaleEstimateAfterUnfolding")
  OPTIONS+=("-u 6" "-u 7" "-u 0" "-u 0" "-u 0")
  BACKGROUNDOPTIONS+=("-u 1" "-u 2" "-b 1" "-u 3" "-u 4")
  BACKGROUNDSTRING+=("mixedConeBackground" "mixedConeBackground" "reflectedConeBackground" "mixedConeBackground" "mixedConeBackground")
else
  VARIATIONS+=("unfoldingWithNominalSmear_backgroundSubtractionSystematics")
  OPTIONS+=("-u 0")
  BACKGROUNDOPTIONS+=("-u 1")
  BACKGROUNDSTRING+=("reflectedConeBackground")
fi

# Loop over all the unfolding variations
for i in "${!VARIATIONS[@]}"; do

  # First, create a new file, where the string template part is replaced by a name from variations array, and background template by the background variation.
  NEWNAME=$(echo ${INPUT} | sed s/${TEMPLATE}/${VARIATIONS[i]}/g | sed s/${BACKGROUNDTEMPLATE}/${BACKGROUNDSTRING[i]}/g)
  cp $INPUT $NEWNAME

  # Once the new file is created, write the unfolding script with correct options for the new file
  echo "./unfoldEnergyEnergyCorrelators.sh ${NEWNAME} ${OPTIONS[i]}" >> $OUTPUT

  # Write also the commands for background subtraction
  echo "./processUnfoldedEnergyEnergyCorrelators.sh ${NEWNAME} ${BACKGROUNDOPTIONS[i]}" >> $OUTPUT
done

# Make the output script file executable
chmod +x $OUTPUT

# If selected, execute the created script
if $EXECUTE; then
  ./${OUTPUT}
  rm $OUTPUT
fi
