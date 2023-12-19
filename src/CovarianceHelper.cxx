/*
 * Implementation of the CovarianceHelper class
 */

// Own includes
#include "CovarianceHelper.h"

/*
 * Contructor
 */
CovarianceHelper::CovarianceHelper() :
  fCard(NULL),
  fInputFile(NULL),
  fnCentralityBins(0),
  fnTrackPtBins(0),
  fnJetPtBins(0),
  fIsPbPbData(false)
{
  // Initialize all histograms to NULL
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < kMaxJetPtBins; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins; iTrackPt++){
        fAverageEECvalue[iCentrality][iJetPt][iTrackPt] = NULL;
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
}

/*
 * Custom constructor
 */
CovarianceHelper::CovarianceHelper(TString inputFileName){

  fInputFile = TFile::Open(inputFileName);

  // Initialize all histograms to NULL
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < kMaxJetPtBins; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins; iTrackPt++){
        fAverageEECvalue[iCentrality][iJetPt][iTrackPt] = NULL;
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

  ReadAverageValueTables();
}

/*
 * Read the smearing tables from the file
 */
void CovarianceHelper::ReadAverageValueTables(){
  
  // Read card from inputfile
  fCard = new BinningCard(fInputFile);

  // Check if we are dealing with pp or PbPb data
  fIsPbPbData = fCard->GetDataType().Contains("PbPb");
  
  // Read the dimensions of the arrays from the card
  fnCentralityBins = fCard->GetNCentralityBins();
  fnTrackPtBins = fCard->GetNTrackPtBinsEEC();
  fnJetPtBins = fCard->GetNJetPtBins();
  
  // Read the centrality bin borders from the card
  for(int iCentrality = 0; iCentrality < fnCentralityBins; iCentrality++){
    fCentralityBinBorders[iCentrality] = fCard->GetLowBinBorderCentrality(iCentrality);
  }
  fCentralityBinBorders[fnCentralityBins] = fCard->GetHighBinBorderCentrality(fnCentralityBins-1);
  
  // Read the track pT bin borders from the card
  for(int iTrackPt = 0; iTrackPt < fnTrackPtBins; iTrackPt++){
    fTrackPtBinBorders[iTrackPt] = fCard->GetLowBinBorderTrackPtEEC(iTrackPt);
  }
  fTrackPtBinBorders[fnTrackPtBins] = fCard->GetHighBinBorderTrackPtEEC(fnTrackPtBins-1);

  // Read the jet pT bin borders from the card
  for(int iJetPt = 0; iJetPt < fnJetPtBins; iJetPt++){
    fJetPtBinBorders[iJetPt] = fCard->GetLowBinBorderJetPt(iJetPt);
  }
  fJetPtBinBorders[fnJetPtBins] = fCard->GetHighBinBorderJetPt(fnJetPtBins-1);

  // Read the average EEC values histograms from the file
  for(int iCentrality = 0; iCentrality < fnCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < fnJetPtBins; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < fnTrackPtBins; iTrackPt++){
        fAverageEECvalue[iCentrality][iJetPt][iTrackPt] = (TH1D*) fInputFile->Get(Form("%s_C%dJ%dT%d", fHistogramName.Data(), iCentrality, iJetPt, iTrackPt));
      }
    }
  }
  
}

/*
 * Find a bin index from an array
 */
int CovarianceHelper::FindBinIndex(const double* array, const int nBins, const double value) const{
  
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
int CovarianceHelper::FindCentralityBin(const double centrality) const{
  return FindBinIndex(fCentralityBinBorders, fnCentralityBins, centrality);
}

/*
 * Find track pT bin index
 */
int CovarianceHelper::FindTrackPtBin(const double trackPt) const{
  return FindBinIndex(fTrackPtBinBorders, fnTrackPtBins, trackPt);
}

/*
 * Find jet pT bin index
 */
int CovarianceHelper::FindJetPtBin(const double jetPt) const{
  return FindBinIndex(fJetPtBinBorders, fnJetPtBins, jetPt);
}

/*
 * Find the average EEC value in the bin corresponding to input values
 *
 *  const double deltaR = DeltaR value for which we need average value
 *  const double centrality = Centrality of the event
 *  const double jetPt = Jet pT
 *  const double trackPt = Track pT
 *
 *  return: Average EEC value in the bin corresponding to input values
 */
double CovarianceHelper::GetAverageValue(const double deltaR, const double centrality, const double jetPt, const double trackPt) const{

  // Determine the correct bin indices from the input values
  const int iCentrality = fIsPbPbData ? FindCentralityBin(centrality) : 0;
  const int iJetPt = FindJetPtBin(jetPt);
  const int iTrackPt = FindTrackPtBin(trackPt);
  const int iDeltaR = fAverageEECvalue[iCentrality][iJetPt][iTrackPt]->FindBin(deltaR);

  return fAverageEECvalue[iCentrality][iJetPt][iTrackPt]->GetBinContent(iDeltaR);

}

/*
 * Find the average EEC value in the bin corresponding to input values
 *
 *  const double deltaR = Value that we want to smear
 *  const double centrality = Centrality of the event
 *  const double jetPt = Jet pT
 *  const double trackPt = Track pT
 *
 *  return: Average EEC value in the bin corresponding to input values
 */
double CovarianceHelper::GetAverageValue(const int iDeltaR, const double centrality, const double jetPt, const double trackPt) const{
  
  // Determine the correct bin indices from the input values
  const int iCentrality = fIsPbPbData ? FindCentralityBin(centrality) : 0;
  const int iJetPt = FindJetPtBin(jetPt);
  const int iTrackPt = FindTrackPtBin(trackPt);

  return fAverageEECvalue[iCentrality][iJetPt][iTrackPt]->GetBinContent(iDeltaR);
 
}
