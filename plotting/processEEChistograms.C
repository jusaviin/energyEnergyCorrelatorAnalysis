#include "EECCard.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECHistogramManager.h"

/*
 * Macro for processing energy-energy correlator histograms. Currently the following processing steps are done:
 *  1) The energy-energy correlator is normalized
 *  2) The reflected cone background estimate is normalized
 *  3) The reflected cone background estimate is subtracted from the normalized energy-energy correlator to extract the signal
 *
 *  Arguments:
 *   TString fileName = File from which the histograms are read and to which the processed histograms are written
 *   TString outputFileName = If given, the processed histograms are written to this file. Otherwise the fileName file is updated.
 *   const int iBackgroundMethod = Index defining the used background subtraction method
 *                                 0: Mixed cone background subtraction
 *                                 1: Reflected cone background subtraction
 *   const int iSystematic = Index for systematic uncertainty estimation for background subtraction.
 *                           0: Nominal results, no systematic uncertainty estimation
 *                           1: Systematic uncertainty derived from 2% centrality shifted simulation
 *                           2: Systematic uncertainty derived from 6% centrality shifted simulation
 *                           3: Lower scaling estimate for signal-to-background ratio after unfolding
 *                           4: Higher scaling estimate for signal-to-background ratio after unfolding
 *   
 */
void processEEChistograms(TString fileName, TString outputFileName, const int iBackgroundMethod, const int iSystematic = 0){

  // Print the file name to console
  cout << "Processing histograms from " << fileName.Data() << endl;
  
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
  EECCard *card = new EECCard(inputFile);
  
  // The git hash here will be replaced by the latest commit hash by processEnergyEnergyCorrelators.sh script
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
  histograms->SetLoadEnergyEnergyCorrelators(true);
  histograms->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationPlus(false);
  histograms->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationMinus(false);
  histograms->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationPlus(false);
  histograms->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationMinus(false);

  // Load the histograms from the file
  histograms->LoadProcessedHistograms();
  
  // Subtract the background from the energy-energy correlator histograms
  histograms->SubtractBackground(iBackgroundMethod, iSystematic);
  
  // Add the processed histograms to the file
  histograms->WriteProcessed(outputFileName,fileWriteMode);
  
}
