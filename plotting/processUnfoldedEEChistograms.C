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
 *   TString fileName = File from which the histograms are read and to which the processed histograms are written
 *   TString outputFileName = If given, the processed histograms are written to this file. Otherwise the fileName file is updated.
 *   const int iSystematic = Index for systematic uncertainty estimation for background subtraction.
 *   const int iBackgroundMethod = Index defining the used background subtraction method
 *                                 0: Mixed cone background subtraction
 *                                 1: Reflected cone background subtraction
 *   const int iSystematic = Index for systematic uncertainty estimation for background subtraction.
 *                           0: Nominal results, no systematic uncertainty estimation
 *                           1: Systematic uncertainty derived from 2% centrality shifted simulation
 *                           2: Systematic uncertainty derived from 6% centrality shifted simulation
 *                           3: Lower scaling estimate for signal-to-background ratio after unfolding
 *                           4: Higher scaling estimate for signal-to-background ratio after unfolding
 *   const int iEnergyEnergyCorrelator = Energy-energy correlator index for the unfolded correlator. Indices are explained in EECHistogramManager.h
 */
void processUnfoldedEEChistograms(TString fileName, TString outputFileName, const int iBackgroundMethod, const int iSystematic, const int iEnergyEnergyCorrelator = EECHistogramManager::kEnergyEnergyCorrelator){

  // Print the file name to console
  cout << "Processing unfolded histograms from " << fileName.Data() << endl;
  
  // We want to update more information to the file
  const char* fileWriteMode = "UPDATE";
  
  // If output file name is not given, just add the processed histograms to the input file
  if(outputFileName == "") outputFileName = fileName;
  
  // Open the input file
  TFile *inputFile = TFile::Open(fileName);
  
  if(inputFile == NULL){
    cout << "Error! The file " << fileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }
  
  // Load the card from the file and read the collision system
  EECCard* card = new EECCard(inputFile);

  // The git hash here will be replaced by the latest commit hash by processUnfoldedEnergyEnergyCorrelators.sh script
  const char* gitHash = "GITHASHHERE";
  card->AddProcessGitHash(gitHash);
  card->AddOneDimensionalVector(EECCard::kBackgroundMethod, iBackgroundMethod);
  card->AddOneDimensionalVector(EECCard::kBackgroundSystematic, iSystematic);
  
  // ============================ //
  //     EECHistogramManager    //
  // ============================ //
    
  // Create and setup a new histogram manager to project and handle the histograms
  EECHistogramManager* histograms = new EECHistogramManager(inputFile,card);
  
  // Load all energy-energy correlators
  histograms->SetLoadEnergyEnergyCorrelators(iEnergyEnergyCorrelator == EECHistogramManager::kEnergyEnergyCorrelator);
  histograms->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationPlus(iEnergyEnergyCorrelator == EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationPlus);
  histograms->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationMinus(iEnergyEnergyCorrelator == EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationMinus);
  histograms->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationPlus(iEnergyEnergyCorrelator == EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationPlus);
  histograms->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationMinus(iEnergyEnergyCorrelator == EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationMinus);

  // Set the bin ranges to those for which the unfolded histograms are available
  histograms->SetCentralityBinRange(card->GetFirstUnfoldedCentralityBin(), card->GetLastUnfoldedCentralityBin());
  histograms->SetTrackPtBinRangeEEC(card->GetFirstUnfoldedTrackPtBin(), card->GetLastUnfoldedTrackPtBin());
  histograms->SetJetPtBinRangeEEC(card->GetFirstUnfoldedJetPtBin(), card->GetLastUnfoldedJetPtBin());

  // Load the histograms from the file
  histograms->LoadProcessedHistograms();
  
  // Subtract the background from the unfolded energy-energy correlator histograms
  histograms->SubtractBackgroundFromUnfolded(iBackgroundMethod, iSystematic);
  
  // Save the processed histograms to the file
  histograms->WriteProcessedAfterUnfolding(outputFileName,fileWriteMode);
  
}
