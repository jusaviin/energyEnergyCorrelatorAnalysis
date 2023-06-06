/*
 * Implementation of SystematicUncertaintyOrganizer
 */

#include "SystematicUncertaintyOrganizer.h"

/*
 * Default constructor
 */
SystematicUncertaintyOrganizer::SystematicUncertaintyOrganizer()
{
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBinsEEC; iTrackPt++){
        for(int iUncertainty = 0; iUncertainty < knUncertaintySources; iUncertainty++){
          fhEnergyEnergyCorrelatorUncertainty[iCentrality][iJetPt][iTrackPt][iUncertainty] = NULL;
        } // Uncertainty loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
}

/*
 * Constructor
 */
SystematicUncertaintyOrganizer::SystematicUncertaintyOrganizer(TFile* inputFile) :
  SystematicUncertaintyOrganizer()
{
  ReadInputFile(inputFile);
}

/*
 * Copy constructor
 */
SystematicUncertaintyOrganizer::SystematicUncertaintyOrganizer(const SystematicUncertaintyOrganizer& in)
{
  // Copy constructor
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBinsEEC; iTrackPt++){
        for(int iUncertainty = 0; iUncertainty < knUncertaintySources; iUncertainty++){
          fhEnergyEnergyCorrelatorUncertainty[iCentrality][iJetPt][iTrackPt][iUncertainty] = in.fhEnergyEnergyCorrelatorUncertainty[iCentrality][iJetPt][iTrackPt][iUncertainty];
        } // Uncertainty loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
}

// Read the input file containing the uncertainty histograms
void SystematicUncertaintyOrganizer::ReadInputFile(TFile* inputFile){
  
  // The systematic uncertainty estimation is done for all unfolded bins. Read the indices for those from the card
  EECCard* card = new EECCard(inputFile);
  
  int firstStudiedCentralityBin = card->GetFirstUnfoldedCentralityBin();
  int lastStudiedCentralityBin = card->GetLastUnfoldedCentralityBin();
  
  int firstStudiedJetPtBinEEC = card->GetFirstUnfoldedJetPtBin();
  int lastStudiedJetPtBinEEC = card->GetLastUnfoldedJetPtBin();
  
  int firstStudiedTrackPtBinEEC = card->GetFirstUnfoldedTrackPtBin();
  int lastStudiedTrackPtBinEEC = card->GetLastUnfoldedTrackPtBin();

  // Helper variables
  TString saveName;
  
  // Read the uncertainties from the file
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
        for(int iUncertainty = 0; iUncertainty < knUncertaintySources; iUncertainty++){
          saveName = Form("systematicUncertainty_%s_C%dJ%dT%d", fSystematicUncertaintyName[iUncertainty].Data(), iCentrality, iJetPt, iTrackPt);
          fhEnergyEnergyCorrelatorUncertainty[iCentrality][iJetPt][iTrackPt][iUncertainty] = (TH1D*) inputFile->Get(saveName);
        } // Uncertainty loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
  
  // ============================================================================================= //
  // For all the uncertainty sources that are not defined in the file, set the uncertainty to zero //
  // ============================================================================================= //

  // First, find at least one histogram that is is not NULL
  TH1D* exampleHistogram;
  bool exampleFound = false;
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin && !exampleFound; iCentrality++){
    for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC && !exampleFound; iJetPt++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC && !exampleFound; iTrackPt++){
        for(int iUncertainty = 0; iUncertainty < knUncertaintySources; iUncertainty++){
          if(fhEnergyEnergyCorrelatorUncertainty[iCentrality][iJetPt][iTrackPt][iUncertainty] != NULL){
            exampleHistogram = fhEnergyEnergyCorrelatorUncertainty[iCentrality][iJetPt][iTrackPt][iUncertainty];
            exampleFound = true;
            break;
          }
        } // Uncertainty loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

  // Find all histograms that are NULL, and set them as error zero histograms
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
        for(int iUncertainty = 0; iUncertainty < knUncertaintySources; iUncertainty++){
          if(fhEnergyEnergyCorrelatorUncertainty[iCentrality][iJetPt][iTrackPt][iUncertainty] == NULL){
            fhEnergyEnergyCorrelatorUncertainty[iCentrality][iJetPt][iTrackPt][iUncertainty] = (TH1D*) exampleHistogram->Clone(Form("zeroErrorHistogram%d%d%d%d", iCentrality, iJetPt, iTrackPt, iUncertainty));
            for(int iBin = 1; iBin <= exampleHistogram->GetNbinsX(); iBin++){
              fhEnergyEnergyCorrelatorUncertainty[iCentrality][iJetPt][iTrackPt][iUncertainty]->SetBinError(iBin, 0);
            } 
          }
        } // Uncertainty loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
  
}

// Getter for absolute systematic uncertainty for energy-energy correlators
TH1D* SystematicUncertaintyOrganizer::GetSystematicUncertainty(const int iCentrality, const int iJetPt, const int iTrackPt, const int iUncertainty) const{
  
  // Return the uncertainty in the selected bin
  return fhEnergyEnergyCorrelatorUncertainty[iCentrality][iJetPt][iTrackPt][iUncertainty];
}

// Getter for a name for the source of systematic uncertainty
TString SystematicUncertaintyOrganizer::GetSystematicUncertaintyName(const int iUncertainty) const{
  return fSystematicUncertaintyName[iUncertainty];
}

// Getter for an axis name for the source of systematic uncertainty
TString SystematicUncertaintyOrganizer::GetUncertaintyAxisName(const int iUncertainty) const{
  return fUncertaintyAxisName[iUncertainty];
}

// Getter for an axis name for the source of systematic uncertainty
int SystematicUncertaintyOrganizer::GetNUncertaintySources() const{
  return knUncertaintySources;
}