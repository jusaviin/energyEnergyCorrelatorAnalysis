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
 */
void processEEChistograms(TString fileName = "veryCoolData_processed.root"){

  // Print the file name to console
  cout << "Processing histograms from " << fileName.Data() << endl;
  
  // We want to update more information to the file
  const char* fileWriteMode = "UPDATE";
  
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
  
  // ============================ //
  //     EECHistogramManager    //
  // ============================ //
    
  // Create and setup a new histogram manager to project and handle the histograms
  EECHistogramManager *histograms = new EECHistogramManager(inputFile,card);
  
  // Load all energy-energy correlators
  histograms->SetLoadEnergyEnergyCorrelators(true);
  histograms->SetLoadEnergyEnergyCorrelatorsJetPt(true);
  histograms->SetLoadEnergyEnergyCorrelatorsUncorrected(true);
  histograms->SetLoadEnergyEnergyCorrelatorsJetPtUncorrected(true);

  // Load the histograms from the file
  histograms->LoadProcessedHistograms();
  
  // Subtract the background from the energy-energy correlator histograms
  histograms->SubtractBackground();
  
  // Add the processed histograms to the input file
  histograms->WriteProcessed(fileName,fileWriteMode);
  
}
