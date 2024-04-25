#include "EECCard.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECHistogramManager.h"

/*
 * Macro for processing the unfolded energy-energy correlator histograms. Currently the following processing steps are done:
 *  1) Normalize the reflected cone background to match the distribution level after unfolding
 *  2) Subtract the normalized reflected cone background from the unfolded distribution to extract signal
 *
 * Notice that the input file must contain both raw and unfolded energy-energy correlators for this to work.
 *
 *  Arguments:
 *   TString baseFileName = File to which the mixed cone background histograms are added
 *   TString addedFileName = File from which the added mixed cone backgroung histograms can be found
 *   const int iEnergyEnergyCorrelator = Energy-energy correlator index for the unfolded correlator. Indices are explained in EECHistogramManager.h
 */
void combineMixedConeBackgrounds(TString baseFileName, TString addedFileName, const int iEnergyEnergyCorrelator = EECHistogramManager::kEnergyEnergyCorrelator){

  // Print the file name to console
  cout << "Adding mixed cone histograms from " << addedFileName.Data() << " to " << baseFileName.Data() << endl;
  
  // We want to update more information to the file
  const char* fileWriteMode = "UPDATE";
  
  // Open the input file
  TFile* inputFile[2];
  inputFile[0] = TFile::Open(baseFileName);
  
  if(inputFile[0] == NULL){
    cout << "Error! The file " << baseFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }

  // Open the file from which the mixed cone histograms are added
  inputFile[1] = TFile::Open(addedFileName);

  if(inputFile[1] == NULL){
    cout << "Error! The file " << addedFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }
  
  // Load the cards from the files
  EECCard* card[2];
  for(int iFile = 0; iFile < 2; iFile++){
    card[iFile] = new EECCard(inputFile[iFile]);
  }
  
  // ============================ //
  //     EECHistogramManager      //
  // ============================ //
    
  // Create and setup a new histogram managers to handle the histograms
  EECHistogramManager* histograms[2];
  for(int iFile = 0; iFile < 2; iFile++){
    histograms[iFile] = new EECHistogramManager(inputFile[iFile],card[iFile]);
  
    // Load all energy-energy correlators
    histograms[iFile]->SetLoadEnergyEnergyCorrelators(iEnergyEnergyCorrelator == EECHistogramManager::kEnergyEnergyCorrelator);
    histograms[iFile]->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationPlus(iEnergyEnergyCorrelator == EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationPlus);
    histograms[iFile]->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationMinus(iEnergyEnergyCorrelator == EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationMinus);
    histograms[iFile]->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationPlus(iEnergyEnergyCorrelator == EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationPlus);
    histograms[iFile]->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationMinus(iEnergyEnergyCorrelator == EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationMinus);

    // Set the bin ranges to those for which the unfolded histograms are available
    histograms[iFile]->SetCentralityBinRange(0, card[iFile]->GetNCentralityBins()-1);
    histograms[iFile]->SetTrackPtBinRangeEEC(0, card[iFile]->GetNTrackPtBinsEEC()-1);
    histograms[iFile]->SetJetPtBinRangeEEC(0, card[iFile]->GetNJetPtBinsEEC()-1);

    // Load the histograms from the file
    histograms[iFile]->LoadProcessedHistograms();
  }
  
  // Subtract the background from the unfolded energy-energy correlator histograms
  histograms[0]->CombineMixedConeBackgrounds(histograms[1]);
  
  // Save the processed histograms to the file
  histograms[0]->WriteCombinedMixedConeHistograms(baseFileName,fileWriteMode);
  
}
