/*
 * Implementation of the SmearingProvider class
 */

// Own includes
#include "SmearingProvider.h"

/*
 * Contructor
 */
SmearingProvider::SmearingProvider() :
  fCard(NULL),
  fInputFile(NULL),
  fRng(0),
  fnCentralityBins(0),
  fnTrackPtBins(0),
  fnJetPtBins(0),
  fDisableSmearing(false),
  fHistogramName("particleDeltaRResponseMatrix"),
  fUseSmearingFactors(false),
  fIsPbPbData(false)
{
  // Initialize all histograms to NULL
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < kMaxJetPtBins; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins; iTrackPt++){
        fSmearingFactors[iCentrality][iJetPt][iTrackPt] = NULL;
        fResponseMatrix[iCentrality][iJetPt][iTrackPt] = NULL;
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
}

/*
 * Custom constructor
 */
SmearingProvider::SmearingProvider(TString inputFileName, TString histogramName, bool useSmearingFactors, bool isPbPbData) :
  fDisableSmearing(false),
  fHistogramName(histogramName),
  fUseSmearingFactors(useSmearingFactors),
  fIsPbPbData(isPbPbData)
{
  fInputFile = TFile::Open(inputFileName);

  // Initialize the random number generator with a random seed
  fRng = new TRandom3();
  fRng->SetSeed(0);

  // Initialize all histograms to NULL
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < kMaxJetPtBins; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins; iTrackPt++){
        fSmearingFactors[iCentrality][iJetPt][iTrackPt] = NULL;
        fResponseMatrix[iCentrality][iJetPt][iTrackPt] = NULL;
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

  ReadSmearingTables();
}

/*
 * Read the smearing tables from the file
 */
void SmearingProvider::ReadSmearingTables(){
  
  // Read card from inputfile
  fCard = new BinningCard(fInputFile);
  
  // Read the dimensions of the arrays from the card
  fnCentralityBins = fCard->GetNCentralityBins();
  fnTrackPtBins = fCard->GetNTrackPtBins();
  fnJetPtBins = fCard->GetNJetPtBins();
  
  // Read the centrality bin borders from the card
  for(int iCentrality = 0; iCentrality < fnCentralityBins; iCentrality++){
    fCentralityBinBorders[iCentrality] = fCard->GetLowBinBorderCentrality(iCentrality);
  }
  fCentralityBinBorders[fnCentralityBins] = fCard->GetHighBinBorderCentrality(fnCentralityBins-1);
  
  // Read the track pT bin borders from the card
  for(int iTrackPt = 0; iTrackPt < fnTrackPtBins; iTrackPt++){
    fTrackPtBinBorders[iTrackPt] = fCard->GetLowBinBorderTrackPt(iTrackPt);
  }
  fTrackPtBinBorders[fnTrackPtBins] = fCard->GetHighBinBorderTrackPt(fnTrackPtBins-1);

  // Read the jet pT bin borders from the card
  for(int iJetPt = 0; iJetPt < fnJetPtBins; iJetPt++){
    fJetPtBinBorders[iJetPt] = fCard->GetLowBinBorderJetPt(iJetPt);
  }
  fJetPtBinBorders[fnJetPtBins] = fCard->GetHighBinBorderJetPt(fnJetPtBins-1);

  // Read all the smearing factors or response matrices from the file
  if(fUseSmearingFactors){
    for(int iCentrality = 0; iCentrality < fnCentralityBins; iCentrality++){
      for(int iJetPt = 0; iJetPt < fnJetPtBins; iJetPt++){
        for(int iTrackPt = 0; iTrackPt < fnTrackPtBins; iTrackPt++){
          fSmearingFactors[iCentrality][iJetPt][iTrackPt] = (TH1D*) fInputFile->Get(Form("%s/%s_C%dJ%dT%d", fHistogramName.Data(), fHistogramName.Data(), iCentrality, iJetPt, iTrackPt));
        }
      }
    }
  } else {
    for(int iCentrality = 0; iCentrality < fnCentralityBins; iCentrality++){
      for(int iJetPt = 0; iJetPt < fnJetPtBins; iJetPt++){
        for(int iTrackPt = 0; iTrackPt < fnTrackPtBins; iTrackPt++){
          fResponseMatrix[iCentrality][iJetPt][iTrackPt] = (TH2D*) fInputFile->Get(Form("%s/%s_C%dJ%dT%d", fHistogramName.Data(), fHistogramName.Data(), iCentrality, iJetPt, iTrackPt));
        }
      }
    }
  }
  
}

/*
 * Find a bin index from an array
 */
int SmearingProvider::FindBinIndex(const double* array, const int nBins, const double value) const{
  
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
int SmearingProvider::FindCentralityBin(const double centrality) const{
  return FindBinIndex(fCentralityBinBorders, fnCentralityBins, centrality);
}

/*
 * Find track pT bin index
 */
int SmearingProvider::FindTrackPtBin(const double trackPt) const{
  return FindBinIndex(fTrackPtBinBorders, fnTrackPtBins, trackPt);
}

/*
 * Find jet pT bin index
 */
int SmearingProvider::FindJetPtBin(const double jetPt) const{
  return FindBinIndex(fJetPtBinBorders, fnJetPtBins, jetPt);
}

/*
 * Get a smeared value based on input value and kinematic region
 *
 *  const double valueInNeedOfSmearing = Value that we want to smear
 *  const double centrality = Centrality of the event
 *  const double triggerPt = Trigger particle pT
 *  const double associatedPt = Associated particle pT
 *
 *  return: Value that is smeared according to response matrix
 */
double SmearingProvider::GetSmearedValue(const double valueInNeedOfSmearing, const double centrality, const double jetPt, const double trackPt) const{
  
  // 
  if(fDisableSmearing) return valueInNeedOfSmearing;

  // Determine the correct bin from the centrality and track pT values
  int iCentrality = fIsPbPbData ? FindCentralityBin(centrality) : 0;
  int iJetPt = FindJetPtBin(jetPt);
  int iTrackPt = FindTrackPtBin(trackPt);

  // Smearing using smearing factors
  if(fUseSmearingFactors){
    double smearingFactor = fSmearingFactors[iCentrality][iJetPt][iTrackPt]->GetRandom(fRng);
    return smearingFactor*valueInNeedOfSmearing;
  }

  // Smearing using response matrix
  int yBin = fResponseMatrix[iCentrality][iJetPt][iTrackPt]->GetYaxis()->FindBin(valueInNeedOfSmearing);
  return fResponseMatrix[iCentrality][iJetPt][iTrackPt]->ProjectionX("dummy",yBin,yBin)->GetRandom(fRng);

 
}

// Setter for disabling the smearing
void SmearingProvider::SetDisableSmearing(const bool disable){
  fDisableSmearing = disable;
}
