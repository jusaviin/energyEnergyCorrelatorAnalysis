#include "../src/EECHistograms.h"

/*
 * Macro for finding average effective prescale from the data. This is calculated using the formula: (Only Jet100) / (Jet80 && Jet100)
 */
void findEffectivePrescale(){

  // File from which the integrals are calculated
  TString inputFileName = "data/eecAnalysis_akFlowJet_findEffectivePrescale_processed_2023-02-28.root";
  
  // Open the input file
  TFile *inputFile = TFile::Open(inputFileName);
  
  if(inputFile == NULL){
    cout << "Error! The file " << inputFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }
  
  // Load the trigger selection histogram from the file
  TH1D *triggerSelectionHistogram = (TH1D*) inputFile->Get("triggers");
  
  // Calculate the effective prescale
  cout << "Average effective prescale: " << (triggerSelectionHistogram->GetBinContent(1+EECHistograms::kOnlyCaloJet100) + triggerSelectionHistogram->GetBinContent(1+EECHistograms::kCaloJet80And100)) / triggerSelectionHistogram->GetBinContent(1+EECHistograms::kCaloJet80And100) << endl;
  
}
