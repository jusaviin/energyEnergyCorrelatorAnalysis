/*
 * Implementation of the TrackPairEfficiencyCorrector class
 */

// Own includes
#include "TrackPairEfficiencyCorrector.h"

/*
 * Contructor
 */
TrackPairEfficiencyCorrector::TrackPairEfficiencyCorrector() :
  fCard(NULL),
  fInputFile(NULL),
  fnCentralityBins(0),
  fnTrackPtBins(0),
  fnJetPtBins(0),
  fDisableCorrection(false),
  fUseSmoothedCorrection(false)
{
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    for(int iTriggerPt = 0; iTriggerPt < kMaxTrackPtBins; iTriggerPt++){
      for(int iAssociatedPt = 0; iAssociatedPt < kMaxTrackPtBins; iAssociatedPt++){
        fCorrectionTable[iCentrality][iTriggerPt][iAssociatedPt] = NULL;
        for(int iJetPt = 0; iJetPt < kMaxJetPtBins; iJetPt++){
          fCorrectionTableCloseToJet[iCentrality][iTriggerPt][iAssociatedPt][iJetPt] = NULL;
        } // Jet pT loop
      } // Associated pT loop
    } // Trigger pT loop
  } // Centrality loop
}

/*
 * Custom constructor
 */
TrackPairEfficiencyCorrector::TrackPairEfficiencyCorrector(TString inputFileName, bool useSmoothedCorrection) :
  fDisableCorrection(false)
{
  fInputFile = TFile::Open(inputFileName, "r");
  fUseSmoothedCorrection = useSmoothedCorrection;
  ReadCorrectionTables();
}

/*
 * Read the track pair efficiency correction tables from the file
 */
void TrackPairEfficiencyCorrector::ReadCorrectionTables(){
  
  // Read card from inputfile
  fCard = new TrackPairEfficiencyCard(fInputFile);
  
  // Read the dimensions of the arrays from the card
  fnCentralityBins = fCard->GetNCentralityBins();
  fnTrackPtBins = fCard->GetNTrackPairPtBins();
  fnJetPtBins = fCard->GetNJetPtBins();
  
  // Read the centrality bin borders from the card
  for(int iCentrality = 0; iCentrality < fnCentralityBins; iCentrality++){
    fCentralityBinBorders[iCentrality] = fCard->GetLowBinBorderCentrality(iCentrality);
  }
  fCentralityBinBorders[fnCentralityBins] = fCard->GetHighBinBorderCentrality(fnCentralityBins-1);
  
  // Read the track pT bin borders from the card
  for(int iTrackPt = 0; iTrackPt < fnTrackPtBins; iTrackPt++){
    fTrackPtBinBorders[iTrackPt] = fCard->GetLowBinBorderTrackPairPt(iTrackPt);
  }
  fTrackPtBinBorders[fnTrackPtBins] = fCard->GetHighBinBorderTrackPairPt(fnTrackPtBins-1);

  // Read the jet pT bin borders from the card
  for(int iJetPt = 0; iJetPt < fnJetPtBins; iJetPt++){
    fJetPtBinBorders[iJetPt] = fCard->GetLowBinBorderJetPt(iJetPt);
  }
  fJetPtBinBorders[fnJetPtBins] = fCard->GetHighBinBorderJetPt(fnJetPtBins-1);
  
  // Select the correct name for the corrections
  const char* histogramName = fUseSmoothedCorrection ? "smoothedTrackPairEfficiencyCorrection" : "trackPairEfficiencyCorrection";  

  // Read all the tables from the file
  for(int iCentrality = 0; iCentrality < fnCentralityBins; iCentrality++){
    for(int iTriggerPt = 0; iTriggerPt < fnTrackPtBins; iTriggerPt++){
      for(int iAssociatedPt = 0; iAssociatedPt <= iTriggerPt; iAssociatedPt++){
        fCorrectionTable[iCentrality][iTriggerPt][iAssociatedPt] = (TH1D*) fInputFile->Get(Form("%s_C%dT%dA%d", histogramName, iCentrality, iTriggerPt, iAssociatedPt));
        for(int iJetPt = 0; iJetPt < fnJetPtBins; iJetPt++){
          fCorrectionTableCloseToJet[iCentrality][iTriggerPt][iAssociatedPt][iJetPt] = (TH1D*) fInputFile->Get(Form("trackPairEfficiencyCorrectionCloseToJet_C%dT%dA%dJ%d", iCentrality, iTriggerPt, iAssociatedPt, iJetPt));
        }
      }
    }
  }
  
}

/*
 * Find a bin index from an array
 */
int TrackPairEfficiencyCorrector::FindBinIndex(const double* array, const int nBins, const double value) const{
  
  // If the value is below the lowest index, just use the lowest index
  if(value < array[0]) return 0;
  
  // Find the value from the arrey and return the bin index
  for(int iBin = 0; iBin < nBins; iBin++){
    if(value < array[iBin+1]) return iBin;
  }
  
  // If the value is higher than the upper bin border, just use the highest index
  return nBins-1;
}


/*
 * Find centrality bin index
 */
int TrackPairEfficiencyCorrector::FindCentralityBin(const double centrality) const{
  return FindBinIndex(fCentralityBinBorders, fnCentralityBins, centrality);
}

/*
 * Find track pT bin index
 */
int TrackPairEfficiencyCorrector::FindTrackPtBin(const double trackPt) const{
  return FindBinIndex(fTrackPtBinBorders, fnTrackPtBins, trackPt);
}

/*
 * Find jet pT bin index
 */
int TrackPairEfficiencyCorrector::FindJetPtBin(const double jetPt) const{
  return FindBinIndex(fJetPtBinBorders, fnJetPtBins, jetPt);
}

/*
 * Get the weight to be given to the reflected cone track based on a data study
 *
 *  const double deltaR = DeltaR between the two tracks
 *  const double centrality = Centrality of the event
 *  const double triggerPt = Trigger particle pT
 *  const double associatedPt = Associated particle pT
 */
double TrackPairEfficiencyCorrector::GetTrackPairEfficiencyCorrection(const double deltaR, const double centrality, const double triggerPt, const double associatedPt, const double jetPt) const{
  
  if(fDisableCorrection) return 1;

  // Determine the correct bin from the centrality and track pT values
  int iCentrality = FindCentralityBin(centrality);
  int iTriggerPt = FindTrackPtBin(triggerPt);
  int iAssociatedPt = FindTrackPtBin(associatedPt);

  // If jet pT is given, use the correction with jet pT bins
  if(jetPt > 0) { // There is no statistics for very high pT jets with very high pT tracks, so relax the jet requirement there
    if(deltaR > 0.4) return 1;
    int iJetPt = FindJetPtBin(jetPt);
    // Ensure that the trigger pT is higher than associated pT
    double correction = 1;
    if(iTriggerPt >= iAssociatedPt){
      correction = 1.0/fCorrectionTableCloseToJet[iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->GetBinContent(fCorrectionTableCloseToJet[iCentrality][iTriggerPt][iAssociatedPt][iJetPt]->FindBin(deltaR));
    } else {
      correction = 1.0/fCorrectionTableCloseToJet[iCentrality][iAssociatedPt][iTriggerPt][iJetPt]->GetBinContent(fCorrectionTableCloseToJet[iCentrality][iAssociatedPt][iTriggerPt][iJetPt]->FindBin(deltaR));
    }

    // Check that the correction is sane
    if(1.0/correction > 0.01 && 1.0/correction < 2){
      return correction;
    }

    // If the correction does not make sense, use the correction without jet pT bins instead
    std::cout << "Insane track pair efficiency correction in bin: triggerPt: " << iTriggerPt << " assoc pT: " << iAssociatedPt << " jetPt: " << iJetPt << std::endl;
    std::cout << "Using the correction without jets instead" << std::endl;
    
  }
    
  // Ensure that the trigger pT is higher than associated pT
  if(iTriggerPt >= iAssociatedPt){
    return 1.0/fCorrectionTable[iCentrality][iTriggerPt][iAssociatedPt]->GetBinContent(fCorrectionTable[iCentrality][iTriggerPt][iAssociatedPt]->FindBin(deltaR));
  }
  return 1.0/fCorrectionTable[iCentrality][iAssociatedPt][iTriggerPt]->GetBinContent(fCorrectionTable[iCentrality][iAssociatedPt][iTriggerPt]->FindBin(deltaR));
}

// Setter for disabling the correction
void TrackPairEfficiencyCorrector::SetDisableCorrection(const bool disable){
  fDisableCorrection = disable;
}
