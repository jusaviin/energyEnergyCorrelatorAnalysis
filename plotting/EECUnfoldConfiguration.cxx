/*
 * Implementation of the EECUnfoldConfiguration class
 */

// Own includes
#include "EECUnfoldConfiguration.h"


/*
 * Contructor
 */
EECUnfoldConfiguration::EECUnfoldConfiguration():
  fIsPbPbData(true),
  fResponseMatrixFileName("")

{
  // Just add four to all best iteration numbers
  for(int iCentrality = 0; iCentrality < kNCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < kNTrackPtBins; iTrackPt++){
      fBestNumberOfIterations[iCentrality][iTrackPt] = 4;
    }
  }

  // Find the actual numbers based on the study 
  InitializeArrays();
}

/*
 * Custom constructor
 */
EECUnfoldConfiguration::EECUnfoldConfiguration(EECCard* card){
  TString collisionSystem = card->GetDataType();
  fIsPbPbData = collisionSystem.Contains("PbPb");
  
  // Just add four to all best iteration numbers
  for(int iCentrality = 0; iCentrality < kNCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < kNTrackPtBins; iTrackPt++){
      fBestNumberOfIterations[iCentrality][iTrackPt] = 4;
    }
  }

  // Find the actual numbers based on the study 
  InitializeArrays();
}

/*
 * Destructor
 */
EECUnfoldConfiguration::~EECUnfoldConfiguration(){

}

/*
 * Initialization for the number of iterations with Bayesian iteration
 */
void EECUnfoldConfiguration::InitializeArrays(){
  
  // ======================================================== //
  //                                                          //
  //  M     M   EEEEEE   TTTTTTT   H    H    OOOOO    DDDD    //
  //  MMM MMM   E           T      H    H   O     O   D   D   //
  //  M  M  M   EEEEEE      T      HHHHHH   O     O   D    D  //
  //  M     M   E           T      H    H   O     O   D   D   //
  //  M     M   EEEEEE      T      H    H    OOOOO    DDDD    //
  //                                                          //
  // ======================================================== //
  
  // The best number of iterations is determined from doing a chi2 fit to the unfolded to truth ratio.
  // The fit is done separately for each jet pT bin, and then the chi2 values are summed together for each response matrix.
  // For these summer chi2 values, we find a minimum of the distribution.
  // The arrays below show the minimum of these chi2 distributions.
  // There is also the response matrix names for which these numbers are determined included in the configurations.

  // The latest configuration for PbPb
  if(fIsPbPbData){

    fResponseMatrixFileName = "data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_moreLowPtBins_unfoldingHistograms_part1_processed_2023-05-20.root";

    fBestNumberOfIterations[0][5] = 5;  // Centrality = 0-10, track pT > 3 GeV
    fBestNumberOfIterations[1][5] = 4;  // Centrality = 10-30, track pT > 3 GeV
    fBestNumberOfIterations[2][5] = 4;  // Centrality = 30-50, track pT > 3 GeV
    fBestNumberOfIterations[3][5] = 3;  // Centrality = 50-90, track pT > 3 GeV

  // The latest configuration for pp
  } else {

    fResponseMatrixFileName = "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_unfoldingTestPart1_processed_2023-05-09.root";

    fBestNumberOfIterations[0][5] = 5;  // track pT > 3 GeV

  }
  
}

/*
 * Get the background scaling factor for the given bin
 *
 *  const std::pair<double,double> centralityBinBorders = Centrality bin borders for the scaling factor
 *  double trackPtBorderLow = Lower track pT border for the studied track pT bin
 *
 *  return: Scaling factor corresponding to the bin with the given bin borders
 */
double EECUnfoldConfiguration::GetNumberOfIterations(const std::pair<double,double> centralityBinBorders, double trackPtBorderLow) const{

  // ******************************************************************** //
  // First, find the bin indices that correspond to the input bin borders //
  // ******************************************************************** //

  // Small number
  const double epsilon = 0.0001;

  // Search if the given centrality bin borders are included in the background scale table
  int centralityIndex = -1;
  if(fIsPbPbData){
    for(int iCentrality = 0; iCentrality < kNCentralityBins; iCentrality++){
      if(TMath::Abs(fCentralityBinBorderLow[iCentrality] - centralityBinBorders.first) < epsilon){
        if(TMath::Abs(fCentralityBinBorderHigh[iCentrality] - centralityBinBorders.second) < epsilon){
        centralityIndex = iCentrality;
          break;
        }
      }
      if(TMath::Abs(fCentralityBinBorderLowShifted[iCentrality] - centralityBinBorders.first) < epsilon){
        if(TMath::Abs(fCentralityBinBorderHighShifted[iCentrality] - centralityBinBorders.second) < epsilon){
        centralityIndex = iCentrality;
          break;
        }
      }
    } // Loop over centrality bins included in the scaling table

    // If centrality bin is not found, print an error message and return -1 to show we did not find a proper scaling factor.
    if(centralityIndex == -1){
      std::cout << "EECUnfoldConfiguration::Error! Centrality bin " << centralityBinBorders.first << "-" << centralityBinBorders.second << " % was not found from the table." << std::endl;
      std::cout << "Cannot provide a background scaling factor. Please check your code for errors." << std::endl;
      return -1;
    }
  } else {
    centralityIndex = 0;
  }

  // Search if the given track pT bin borders are included in the background scale table
  int trackPtIndex = -1;
  for(int iTrackPt = 0; iTrackPt < kNTrackPtBins; iTrackPt++){
    if(TMath::Abs(fTrackPtBinBorderLow[iTrackPt] - trackPtBorderLow) < epsilon){
      trackPtIndex = iTrackPt;
      break;
    }
  } // Loop over track pT bins included in the scaling table

  // If track pT bin is not found, print an error message and return -1 to show we did not find a proper scaling factor.
  if(trackPtIndex == -1){
    std::cout << "EECUnfoldConfiguration::Error! Track pT lower border " << trackPtBorderLow << " GeV was not found from the table." << std::endl;
    std::cout << "Cannot provide a background scaling factor. Please check your code for errors." << std::endl;
    return -1;
  }

  // *********************************************************************************** //
  // After the correct bin has been found, return the corresponding number of iterations //
  // *********************************************************************************** //

  return fBestNumberOfIterations[centralityIndex][trackPtIndex];

}

// Getter for the name of the file containing the response matrices
TString EECUnfoldConfiguration::GetResponseFileName() const{
  return fResponseMatrixFileName;
}
