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
  fDisableCorrection(false)
{
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    for(int iTriggerPt = 0; iTriggerPt < kMaxTrackPtBins; iTriggerPt++){
      for(int iAssociatedPt = 0; iAssociatedPt < kMaxTrackPtBins; iAssociatedPt++){
        fCorrectionTable[iCentrality][iTriggerPt][iAssociatedPt] = NULL;
      } // Associated pT loop
    } // Trigger pT loop
  } // Centrality loop
}

/*
 * Custom constructor
 */
TrackPairEfficiencyCorrector::TrackPairEfficiencyCorrector(TString inputFileName) :
  fDisableCorrection(false)
{
  fInputFile = TFile::Open(inputFileName, "r");
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
  
  // Read all the tables from the file
  for(int iCentrality = 0; iCentrality < fnCentralityBins; iCentrality++){
    for(int iTriggerPt = 0; iTriggerPt < fnTrackPtBins; iTriggerPt++){
      for(int iAssociatedPt = 0; iAssociatedPt <= iTriggerPt; iAssociatedPt++){
        fCorrectionTable[iCentrality][iTriggerPt][iAssociatedPt] = (TH1D*) fInputFile->Get(Form("smoothedTrackPairEfficiencyCorrection_C%dT%dA%d", iCentrality, iTriggerPt, iAssociatedPt));
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
 * Get the weight to be given to the reflected cone track based on a data study
 *
 *  const double deltaR = DeltaR between the two tracks
 *  const double centrality = Centrality of the event
 *  const double triggerPt = Trigger particle pT
 *  const double associatedPt = Associated particle pT
 */
double TrackPairEfficiencyCorrector::GetTrackPairEfficiencyCorrection(const double deltaR, const double centrality, const double triggerPt, const double associatedPt) const{
  
  // Determine the correct bin from the centrality and track pT values
  int iCentrality = FindCentralityBin(centrality);
  int iTriggerPt = FindTrackPtBin(triggerPt);
  int iAssociatedPt = FindTrackPtBin(associatedPt);
    
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
