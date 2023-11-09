#!/bin/bash

if [ "$#" -lt 1 ]; then
  echo "Usage of the script:"
  echo "$0 systemString [-c comment] [-n] [-s]"
  echo "systemString = PbPb or pp"
  echo "-c comment = Comment added to all sent jobs. Default: \"\""
  echo "-n = Send jobs for nominal unfolding parameters"
  echo "-s = Send jobs for all systematic variations"
  exit
fi

SYSTEM=$1  # Collision system for the sent jobs
shift      # Shift the positional parameters to read the optional ones

# Enable case insensitive matching
# shopt -s nocasematch

# Check that the system string is sensible
if [[ $SYSTEM == "PbPb" ]]; then
  CARD="responseMatrixVariationCardPbPbMC.input"
  CRAB="responseMatrixVariationCrabPbPbMC.py"
elif [[ $SYSTEM == "pp" ]]; then
  CARD="responseMatrixVariationCardPpMC.input"
  CRAB="responseMatrixVariationCrabPpMC.py"
else
  echo "Unknown system: ${SYSTEM}"
  echo "Please give either PbPb or pp"
  exit
fi

# Disable case insensitive matching
# shopt -u nocasematch

# Read the optional arguments. (Semicolon after letter: expects argument)
while getopts ":c:ns" opt; do
case $opt in
c) COMMENT="$OPTARG"
# If a comment is added, make sure it ends with "_"
if [[ $COMMENT != *_ ]]; then
  COMMENT+="_"
fi
;;
n) NOMINAL=true
;;
s) VARIATIONS=true
;;
\?) echo "Invalid option -$OPTARG" >&2
exit 1
;;
esac
done

# Set default values to optional arguments if they are not given
COMMENT=${COMMENT:-""}
NOMINAL=${NOMINAL:-false}
VARIATIONS=${VARIATIONS:-false}

#################################
##       Nominal results       ##
#################################

if $NOMINAL; then

  ##### Update configuration in EECCard #####

  # Select response matrix histograms to be filled in the input card
  LINE=`grep "FilledHistograms" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "131")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define centrality bins only if we are looking at PbPb system
  if [[ $SYSTEM == "PbPb" ]]; then
    # Define 4% centrality shift
    LINE=`grep "CentralityBinEdges" ${CARD}`
    NEWLINE=$(./updateString.sh "${LINE}" "1 2 3 4 5" "4 14 34 54 94")
    sed -i "s/${LINE}/${NEWLINE}/" $CARD
  fi

  # Define nominal jet pT smearing
  LINE=`grep "JetUncertainty" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "3")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define no jet pT weighting
  LINE=`grep "JetPtWeight" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "0")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Select GenGen as correlation type
  LINE=`grep "McCorrelationType" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "3")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD


  ##### Update CRAB configuration #####

  # Select first half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'1\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Define GenGen as the system name
  LINE=`grep "system =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'GenGen\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Define centrality bins only if we are looking at PbPb system
  if [[ $SYSTEM == "PbPb" ]]; then
    # Define 4% centrality shift
    LINE=`grep "centralityShift =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'4pC\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB
  fi

  # Define comment to be given to this specific job
  LINE=`grep "comment =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}nomSmear_resMat\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  # Select second half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'2\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  ##### Truth reference for nominal results #####

  # Select regular energy-energy correlator histograms to be filled in the input card
  LINE=`grep "FilledHistograms" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "31")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define comment to be given to this specific job
  LINE=`grep "comment =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}nomSmear_truthRef\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  # Select first half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'1\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}


  ##### Reconstructed reference for nominal results #####

  # Select RecoGen as correlation type
  LINE=`grep "McCorrelationType" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "1")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define RecoGen as the system name
  LINE=`grep "system =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'RecoGen\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Define comment to be given to this specific job
  LINE=`grep "comment =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}nomSmear_recoRef\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  # Select second half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'2\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

fi

if $VARIATIONS; then

  # Centrality variations are done only for PbPb system
  if [[ $SYSTEM == "PbPb" ]]; then

    #########################################
    ##  Results for 2% shifted centrality  ##
    #########################################

    ##### Update configuration in EECCard #####

    # Select response matrix histograms to be filled in the input card
    LINE=`grep "FilledHistograms" ${CARD}`
    NEWLINE=$(./updateString.sh "${LINE}" "1" "131")
    sed -i "s/${LINE}/${NEWLINE}/" $CARD

    # Define 2% centrality shift
    LINE=`grep "CentralityBinEdges" ${CARD}`
    NEWLINE=$(./updateString.sh "${LINE}" "1 2 3 4 5" "2 12 32 52 92")
    sed -i "s/${LINE}/${NEWLINE}/" $CARD

    # Define nominal jet pT smearing
    LINE=`grep "JetUncertainty" ${CARD}`
    NEWLINE=$(./updateString.sh "${LINE}" "1" "3")
    sed -i "s/${LINE}/${NEWLINE}/" $CARD

    # Define no jet pT weighting
    LINE=`grep "JetPtWeight" ${CARD}`
    NEWLINE=$(./updateString.sh "${LINE}" "1" "0")
    sed -i "s/${LINE}/${NEWLINE}/" $CARD

    # Select GenGen as correlation type
    LINE=`grep "McCorrelationType" ${CARD}`
    NEWLINE=$(./updateString.sh "${LINE}" "1" "3")
    sed -i "s/${LINE}/${NEWLINE}/" $CARD


    ##### Update CRAB configuration #####

    # Select first half of MC statistics
    LINE=`grep "iPart =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'1\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB

    # Define GenGen as the system name
    LINE=`grep "system =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'GenGen\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB

    # Define 2% centrality shift
    LINE=`grep "centralityShift =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'2pC\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB

    # Define comment to be given to this specific job
    LINE=`grep "comment =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}nomSmear_resMat\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB

    # Send jobs
    crab submit -c ${CRAB}

    # Select second half of MC statistics
    LINE=`grep "iPart =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'2\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB

    # Send jobs
    crab submit -c ${CRAB}

    ##### Truth reference for 2% shifted results #####

    # Select regular energy-energy correlator histograms to be filled in the input card
    LINE=`grep "FilledHistograms" ${CARD}`
    NEWLINE=$(./updateString.sh "${LINE}" "1" "31")
    sed -i "s/${LINE}/${NEWLINE}/" $CARD

    # Define comment to be given to this specific job
    LINE=`grep "comment =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}nomSmear_truthRef\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB

    # Send jobs
    crab submit -c ${CRAB}

    # Select first half of MC statistics
    LINE=`grep "iPart =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'1\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB

    # Send jobs
    crab submit -c ${CRAB}


    ##### Reconstructed reference for 2% shifted results #####

    # Select RecoGen as correlation type
    LINE=`grep "McCorrelationType" ${CARD}`
    NEWLINE=$(./updateString.sh "${LINE}" "1" "1")
    sed -i "s/${LINE}/${NEWLINE}/" $CARD

    # Define RecoGen as the system name
    LINE=`grep "system =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'RecoGen\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB

    # Define comment to be given to this specific job
    LINE=`grep "comment =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}nomSmear_recoRef\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB

    # Send jobs
    crab submit -c ${CRAB}

    # Select second half of MC statistics
    LINE=`grep "iPart =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'2\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB

    # Send jobs
    crab submit -c ${CRAB}

    #########################################
    ##  Results for 6% shifted centrality  ##
    #########################################

    ##### Update configuration in EECCard #####

    # Select response matrix histograms to be filled in the input card
    LINE=`grep "FilledHistograms" ${CARD}`
    NEWLINE=$(./updateString.sh "${LINE}" "1" "131")
    sed -i "s/${LINE}/${NEWLINE}/" $CARD

    # Define 6% centrality shift
    LINE=`grep "CentralityBinEdges" ${CARD}`
    NEWLINE=$(./updateString.sh "${LINE}" "1 2 3 4 5" "6 16 36 56 96")
    sed -i "s/${LINE}/${NEWLINE}/" $CARD

    # Define nominal jet pT smearing
    LINE=`grep "JetUncertainty" ${CARD}`
    NEWLINE=$(./updateString.sh "${LINE}" "1" "3")
    sed -i "s/${LINE}/${NEWLINE}/" $CARD

    # Define no jet pT weighting
    LINE=`grep "JetPtWeight" ${CARD}`
    NEWLINE=$(./updateString.sh "${LINE}" "1" "0")
    sed -i "s/${LINE}/${NEWLINE}/" $CARD

    # Select GenGen as correlation type
    LINE=`grep "McCorrelationType" ${CARD}`
    NEWLINE=$(./updateString.sh "${LINE}" "1" "3")
    sed -i "s/${LINE}/${NEWLINE}/" $CARD


    ##### Update CRAB configuration #####

    # Select first half of MC statistics
    LINE=`grep "iPart =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'1\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB

    # Define GenGen as the system name
    LINE=`grep "system =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'GenGen\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB

    # Define 6% centrality shift
    LINE=`grep "centralityShift =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'6pC\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB

    # Define comment to be given to this specific job
    LINE=`grep "comment =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}nomSmear_resMat\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB

    # Send jobs
    crab submit -c ${CRAB}

    # Select second half of MC statistics
    LINE=`grep "iPart =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'2\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB

    # Send jobs
    crab submit -c ${CRAB}

    ##### Truth reference for 6% shifted results #####

    # Select regular energy-energy correlator histograms to be filled in the input card
    LINE=`grep "FilledHistograms" ${CARD}`
    NEWLINE=$(./updateString.sh "${LINE}" "1" "31")
    sed -i "s/${LINE}/${NEWLINE}/" $CARD

    # Define comment to be given to this specific job
    LINE=`grep "comment =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}nomSmear_truthRef\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB

    # Send jobs
    crab submit -c ${CRAB}

    # Select first half of MC statistics
    LINE=`grep "iPart =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'1\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB

    # Send jobs
    crab submit -c ${CRAB}


    ##### Reconstructed reference for 6% shifted results #####

    # Select RecoGen as correlation type
    LINE=`grep "McCorrelationType" ${CARD}`
    NEWLINE=$(./updateString.sh "${LINE}" "1" "1")
    sed -i "s/${LINE}/${NEWLINE}/" $CARD

    # Define RecoGen as the system name
    LINE=`grep "system =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'RecoGen\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB

    # Define comment to be given to this specific job
    LINE=`grep "comment =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}nomSmear_recoRef\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB

    # Send jobs
    crab submit -c ${CRAB}

    # Select second half of MC statistics
    LINE=`grep "iPart =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'2\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB

    # Send jobs
    crab submit -c ${CRAB}

  fi

  ####################################################
  ##  Results with negative jet energy scale shift  ##
  ####################################################

  ##### Update configuration in EECCard #####

  # Select response matrix histograms to be filled in the input card
  LINE=`grep "FilledHistograms" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "131")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define centrality bins only if we are looking at PbPb system
  if [[ $SYSTEM == "PbPb" ]]; then
    # Define 4% centrality shift
    LINE=`grep "CentralityBinEdges" ${CARD}`
    NEWLINE=$(./updateString.sh "${LINE}" "1 2 3 4 5" "4 14 34 54 94")
    sed -i "s/${LINE}/${NEWLINE}/" $CARD
  fi

  # Define negative jet energy scale shift
  LINE=`grep "JetUncertainty" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "1")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define no jet pT weighting
  LINE=`grep "JetPtWeight" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "0")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Select GenGen as correlation type
  LINE=`grep "McCorrelationType" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "3")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD


  ##### Update CRAB configuration #####

  # Select first half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'1\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Define GenGen as the system name
  LINE=`grep "system =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'GenGen\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Define centrality bins only if we are looking at PbPb system
  if [[ $SYSTEM == "PbPb" ]]; then
    # Define 4% centrality shift
    LINE=`grep "centralityShift =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'4pC\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB
  fi

  # Define comment to be given to this specific job
  LINE=`grep "comment =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}minusJECunc_resMat\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  # Select second half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'2\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  ##### Truth reference for results with negative jet energy scale shift #####

  # Select regular energy-energy correlator histograms to be filled in the input card
  LINE=`grep "FilledHistograms" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "31")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define comment to be given to this specific job
  LINE=`grep "comment =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}minusJECunc_truthRef\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  # Select first half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'1\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}


  ##### Reconstructed reference for results with negative jet energy scale shift #####

  # Select RecoGen as correlation type
  LINE=`grep "McCorrelationType" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "1")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define RecoGen as the system name
  LINE=`grep "system =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'RecoGen\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Define comment to be given to this specific job
  LINE=`grep "comment =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}minusJECunc_recoRef\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  # Select second half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'2\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}


  ####################################################
  ##  Results with positive jet energy scale shift  ##
  ####################################################

  ##### Update configuration in EECCard #####

  # Select response matrix histograms to be filled in the input card
  LINE=`grep "FilledHistograms" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "131")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define centrality bins only if we are looking at PbPb system
  if [[ $SYSTEM == "PbPb" ]]; then
    # Define 4% centrality shift
    LINE=`grep "CentralityBinEdges" ${CARD}`
    NEWLINE=$(./updateString.sh "${LINE}" "1 2 3 4 5" "4 14 34 54 94")
    sed -i "s/${LINE}/${NEWLINE}/" $CARD
  fi

  # Define negative jet energy scale shift
  LINE=`grep "JetUncertainty" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "2")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define no jet pT weighting
  LINE=`grep "JetPtWeight" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "0")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Select GenGen as correlation type
  LINE=`grep "McCorrelationType" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "3")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD


  ##### Update CRAB configuration #####

  # Select first half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'1\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Define GenGen as the system name
  LINE=`grep "system =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'GenGen\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Define centrality bins only if we are looking at PbPb system
  if [[ $SYSTEM == "PbPb" ]]; then
    # Define 4% centrality shift
    LINE=`grep "centralityShift =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'4pC\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB
  fi

  # Define comment to be given to this specific job
  LINE=`grep "comment =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}plusJECunc_resMat\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  # Select second half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'2\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  ##### Truth reference for results with positive jet energy scale shift #####

  # Select regular energy-energy correlator histograms to be filled in the input card
  LINE=`grep "FilledHistograms" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "31")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define comment to be given to this specific job
  LINE=`grep "comment =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}plusJECunc_truthRef\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  # Select first half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'1\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}


  ##### Reconstructed reference for result with positive jet energy scale shift #####

  # Select RecoGen as correlation type
  LINE=`grep "McCorrelationType" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "1")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define RecoGen as the system name
  LINE=`grep "system =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'RecoGen\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Define comment to be given to this specific job
  LINE=`grep "comment =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}plusJECunc_recoRef\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  # Select second half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'2\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}


  ########################################
  ##  Results with less smearing in MC  ##
  ########################################

  ##### Update configuration in EECCard #####

  # Select response matrix histograms to be filled in the input card
  LINE=`grep "FilledHistograms" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "131")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define centrality bins only if we are looking at PbPb system
  if [[ $SYSTEM == "PbPb" ]]; then
    # Define 4% centrality shift
    LINE=`grep "CentralityBinEdges" ${CARD}`
    NEWLINE=$(./updateString.sh "${LINE}" "1 2 3 4 5" "4 14 34 54 94")
    sed -i "s/${LINE}/${NEWLINE}/" $CARD
  fi

  # Define less jet pT smearing
  LINE=`grep "JetUncertainty" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "4")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define no jet pT weighting
  LINE=`grep "JetPtWeight" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "0")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Select GenGen as correlation type
  LINE=`grep "McCorrelationType" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "3")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD


  ##### Update CRAB configuration #####

  # Select first half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'1\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Define GenGen as the system name
  LINE=`grep "system =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'GenGen\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Define centrality bins only if we are looking at PbPb system
  if [[ $SYSTEM == "PbPb" ]]; then
    # Define 4% centrality shift
    LINE=`grep "centralityShift =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'4pC\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB
  fi

  # Define comment to be given to this specific job
  LINE=`grep "comment =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}uncSmearDown_resMat\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  # Select second half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'2\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  ##### Truth reference for results with less smearing in MC #####

  # Select regular energy-energy correlator histograms to be filled in the input card
  LINE=`grep "FilledHistograms" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "31")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define comment to be given to this specific job
  LINE=`grep "comment =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}uncSmearDown_truthRef\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  # Select first half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'1\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}


  ##### Reconstructed reference for results with less smearing in MC #####

  # Select RecoGen as correlation type
  LINE=`grep "McCorrelationType" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "1")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define RecoGen as the system name
  LINE=`grep "system =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'RecoGen\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Define comment to be given to this specific job
  LINE=`grep "comment =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}uncSmearDown_recoRef\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  # Select second half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'2\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}


  ########################################
  ##  Results with more smearing in MC  ##
  ########################################

  ##### Update configuration in EECCard #####

  # Select response matrix histograms to be filled in the input card
  LINE=`grep "FilledHistograms" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "131")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define centrality bins only if we are looking at PbPb system
  if [[ $SYSTEM == "PbPb" ]]; then
    # Define 4% centrality shift
    LINE=`grep "CentralityBinEdges" ${CARD}`
    NEWLINE=$(./updateString.sh "${LINE}" "1 2 3 4 5" "4 14 34 54 94")
    sed -i "s/${LINE}/${NEWLINE}/" $CARD
  fi

  # Define more jet pT smearing
  LINE=`grep "JetUncertainty" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "5")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define no jet pT weighting
  LINE=`grep "JetPtWeight" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "0")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Select GenGen as correlation type
  LINE=`grep "McCorrelationType" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "3")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD


  ##### Update CRAB configuration #####

  # Select first half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'1\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Define GenGen as the system name
  LINE=`grep "system =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'GenGen\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Define centrality bins only if we are looking at PbPb system
  if [[ $SYSTEM == "PbPb" ]]; then
    # Define 4% centrality shift
    LINE=`grep "centralityShift =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'4pC\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB
  fi

  # Define comment to be given to this specific job
  LINE=`grep "comment =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}uncSmearUp_resMat\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  # Select second half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'2\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  ##### Truth reference for results with more smearing in MC #####

  # Select regular energy-energy correlator histograms to be filled in the input card
  LINE=`grep "FilledHistograms" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "31")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define comment to be given to this specific job
  LINE=`grep "comment =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}uncSmearUp_truthRef\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  # Select first half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'1\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}


  ##### Reconstructed reference for results with more smearing in MC #####

  # Select RecoGen as correlation type
  LINE=`grep "McCorrelationType" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "1")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define RecoGen as the system name
  LINE=`grep "system =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'RecoGen\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Define comment to be given to this specific job
  LINE=`grep "comment =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}uncSmearUp_recoRef\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  # Select second half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'2\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}


  #############################################
  ##  Results with weighted jet pT spectrum  ##
  #############################################

  ##### Update configuration in EECCard #####

  # Select response matrix histograms to be filled in the input card
  LINE=`grep "FilledHistograms" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "131")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define centrality bins only if we are looking at PbPb system
  if [[ $SYSTEM == "PbPb" ]]; then
    # Define 4% centrality shift
    LINE=`grep "CentralityBinEdges" ${CARD}`
    NEWLINE=$(./updateString.sh "${LINE}" "1 2 3 4 5" "4 14 34 54 94")
    sed -i "s/${LINE}/${NEWLINE}/" $CARD
  fi

  # Define nominal jet pT smearing
  LINE=`grep "JetUncertainty" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "3")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define jet pT weighting
  LINE=`grep "JetPtWeight" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "1")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Select GenGen as correlation type
  LINE=`grep "McCorrelationType" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "3")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD


  ##### Update CRAB configuration #####

  # Select first half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'1\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Define GenGen as the system name
  LINE=`grep "system =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'GenGen\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Define centrality bins only if we are looking at PbPb system
  if [[ $SYSTEM == "PbPb" ]]; then
    # Define 4% centrality shift
    LINE=`grep "centralityShift =" ${CRAB}`
    NEWLINE=$(./updateString.sh "${LINE}" "2" "\'4pC\'")
    sed -i "s/${LINE}/${NEWLINE}/" $CRAB
  fi

  # Define comment to be given to this specific job
  LINE=`grep "comment =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}nomSmear_jetPtW_resMat\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  # Select second half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'2\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  ##### Truth reference for jet pT weighted results #####

  # Select regular energy-energy correlator histograms to be filled in the input card
  LINE=`grep "FilledHistograms" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "31")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define comment to be given to this specific job
  LINE=`grep "comment =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}nomSmear_jetPtW_truthRef\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  # Select first half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'1\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}


  ##### Reconstructed reference for jet pT weighted results #####

  # Select RecoGen as correlation type
  LINE=`grep "McCorrelationType" ${CARD}`
  NEWLINE=$(./updateString.sh "${LINE}" "1" "1")
  sed -i "s/${LINE}/${NEWLINE}/" $CARD

  # Define RecoGen as the system name
  LINE=`grep "system =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'RecoGen\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Define comment to be given to this specific job
  LINE=`grep "comment =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'${COMMENT}nomSmear_jetPtW_recoRef\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

  # Select second half of MC statistics
  LINE=`grep "iPart =" ${CRAB}`
  NEWLINE=$(./updateString.sh "${LINE}" "2" "\'2\'")
  sed -i "s/${LINE}/${NEWLINE}/" $CRAB

  # Send jobs
  crab submit -c ${CRAB}

fi
