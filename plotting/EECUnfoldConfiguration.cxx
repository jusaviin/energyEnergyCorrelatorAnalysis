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
  fSystematicIndex(0),
  fSplitIndex(0),
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
EECUnfoldConfiguration::EECUnfoldConfiguration(EECCard* card, const int iSplit, const int iSystematic){
  TString collisionSystem = card->GetDataType();
  fIsPbPbData = collisionSystem.Contains("PbPb");

  fSplitIndex = iSplit;
  fSystematicIndex = iSystematic;
  
  // Just add four to all best iteration numbers
  for(int iCentrality = 0; iCentrality < kNCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < kNTrackPtBins; iTrackPt++){
      fBestNumberOfIterations[iCentrality][iTrackPt] = 4;
    }
  }

    // Check that the split number is sensible
  if(fSplitIndex < 0 || fSplitIndex >= kNDatasetSplits){
    std::cout << "EECUnfoldConfiguration::ERROR! Selected split index for the response matrix is out or range!" << std::endl;
    std::cout << "You selected " << fSplitIndex << " while the allowed range is 0-" << kNDatasetSplits-1 << std::endl;
    std::cout << "Next the code will crash, or you will see undefined behavior. Please check your code." << std::endl;
    return;
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

  // Helper to select response matrix based on the selected Monte Carlo split
  const char* splitName[kNDatasetSplits] = {"", "_part1", "_part2"};

  // The latest configuration for PbPb
  if(fIsPbPbData){

    if(fSystematicIndex == kDefault){
      // Configuration for default response matrix
      // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
      // Input: chi2Histograms_PbPb_regular_threeTrackPtBins_2023-06-02.root
      //        chi2Histograms_PbPb_swapped_threeTrackPtBins_2023-06-02.root 

      fResponseMatrixFileName = Form("data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_moreLowPtBins_unfoldingHistograms%s_processed_2023-05-20.root", splitName[fSplitIndex]);

      fBestNumberOfIterations[0][3] = 4;  // Centrality = 0-10, track pT > 2 GeV
      fBestNumberOfIterations[1][3] = 4;  // Centrality = 10-30, track pT > 2 GeV
      fBestNumberOfIterations[2][3] = 3;  // Centrality = 30-50, track pT > 2 GeV
      fBestNumberOfIterations[3][3] = 3;  // Centrality = 50-90, track pT > 2 GeV

      fBestNumberOfIterations[0][4] = 4;  // Centrality = 0-10, track pT > 2.5 GeV
      fBestNumberOfIterations[1][4] = 4;  // Centrality = 10-30, track pT > 2.5 GeV
      fBestNumberOfIterations[2][4] = 3;  // Centrality = 30-50, track pT > 2.5 GeV
      fBestNumberOfIterations[3][4] = 3;  // Centrality = 50-90, track pT > 2.5 GeV

      fBestNumberOfIterations[0][5] = 5;  // Centrality = 0-10, track pT > 3 GeV
      fBestNumberOfIterations[1][5] = 4;  // Centrality = 10-30, track pT > 3 GeV
      fBestNumberOfIterations[2][5] = 4;  // Centrality = 30-50, track pT > 3 GeV
      fBestNumberOfIterations[3][5] = 3;  // Centrality = 50-90, track pT > 3 GeV

    } else if(fSystematicIndex == kJetPtResolutionUncertainty){
      // Configuration for jet pT resolution uncertainty evaluation
      // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
      // Input: chi2Histograms_PbPb_split1_systematicForJetPtResolution_2023-06-05.root
      //        chi2Histograms_PbPb_split2_systematicForJetPtResolution_2023-06-05.root

      fResponseMatrixFileName = Form("data/PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_noJetPtWeight_smearJetPtResolution_responseMatrix%s_processed_2023-06-02.root", splitName[fSplitIndex]);

      fBestNumberOfIterations[0][3] = 1;  // Centrality = 0-10, track pT > 2 GeV
      fBestNumberOfIterations[1][3] = 2;  // Centrality = 10-30, track pT > 2 GeV
      fBestNumberOfIterations[2][3] = 2;  // Centrality = 30-50, track pT > 2 GeV
      fBestNumberOfIterations[3][3] = 2;  // Centrality = 50-90, track pT > 2 GeV

      fBestNumberOfIterations[0][4] = 1;  // Centrality = 0-10, track pT > 2.5 GeV
      fBestNumberOfIterations[1][4] = 2;  // Centrality = 10-30, track pT > 2.5 GeV
      fBestNumberOfIterations[2][4] = 2;  // Centrality = 30-50, track pT > 2.5 GeV
      fBestNumberOfIterations[3][4] = 2;  // Centrality = 50-90, track pT > 2.5 GeV

      fBestNumberOfIterations[0][5] = 1;  // Centrality = 0-10, track pT > 3 GeV
      fBestNumberOfIterations[1][5] = 2;  // Centrality = 10-30, track pT > 3 GeV
      fBestNumberOfIterations[2][5] = 2;  // Centrality = 30-50, track pT > 3 GeV
      fBestNumberOfIterations[3][5] = 2;  // Centrality = 50-90, track pT > 3 GeV

    } else if(fSystematicIndex == kJetEnergyScaleUncertainty){
      // Configuration for jet energy scale uncertainty evaluation
      // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
      // Input: chi2Histograms_PbPb_split1_systematicForJetPtUncertainty_2023-06-05.root
      //        chi2Histograms_PbPb_split2_systematicForJetPtUncertainty_2023-06-05.root

      fResponseMatrixFileName = Form("data/PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_noJetPtWeight_smearJetPtUncertainty_responseMatrix%s_processed_2023-06-02.root", splitName[fSplitIndex]);

      fBestNumberOfIterations[0][3] = 3;  // Centrality = 0-10, track pT > 2 GeV
      fBestNumberOfIterations[1][3] = 3;  // Centrality = 10-30, track pT > 2 GeV
      fBestNumberOfIterations[2][3] = 3;  // Centrality = 30-50, track pT > 2 GeV
      fBestNumberOfIterations[3][3] = 3;  // Centrality = 50-90, track pT > 2 GeV

      fBestNumberOfIterations[0][4] = 3;  // Centrality = 0-10, track pT > 2.5 GeV
      fBestNumberOfIterations[1][4] = 3;  // Centrality = 10-30, track pT > 2.5 GeV
      fBestNumberOfIterations[2][4] = 3;  // Centrality = 30-50, track pT > 2.5 GeV
      fBestNumberOfIterations[3][4] = 2;  // Centrality = 50-90, track pT > 2.5 GeV

      fBestNumberOfIterations[0][5] = 4;  // Centrality = 0-10, track pT > 3 GeV
      fBestNumberOfIterations[1][5] = 3;  // Centrality = 10-30, track pT > 3 GeV
      fBestNumberOfIterations[2][5] = 3;  // Centrality = 30-50, track pT > 3 GeV
      fBestNumberOfIterations[3][5] = 3;  // Centrality = 50-90, track pT > 3 GeV

    } else if(fSystematicIndex == kJetPtPriorUncertainty){
      // Configuration for jet pT prior uncertainty evaluation
      // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
      // Input: chi2Histograms_PbPb_split1_jetPtWeight_2023-06-16.root
      //        chi2Histograms_PbPb_split2_jetPtWeight_2023-06-16.root

      fResponseMatrixFileName = Form("data/PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_jetPtWeight_responseMatrix%s_processed_2023-06-15.root", splitName[fSplitIndex]);

      fBestNumberOfIterations[0][3] = 4;  // Centrality = 0-10, track pT > 2 GeV
      fBestNumberOfIterations[1][3] = 3;  // Centrality = 10-30, track pT > 2 GeV
      fBestNumberOfIterations[2][3] = 3;  // Centrality = 30-50, track pT > 2 GeV
      fBestNumberOfIterations[3][3] = 3;  // Centrality = 50-90, track pT > 2 GeV

      fBestNumberOfIterations[0][4] = 4;  // Centrality = 0-10, track pT > 2.5 GeV
      fBestNumberOfIterations[1][4] = 4;  // Centrality = 10-30, track pT > 2.5 GeV
      fBestNumberOfIterations[2][4] = 3;  // Centrality = 30-50, track pT > 2.5 GeV
      fBestNumberOfIterations[3][4] = 3;  // Centrality = 50-90, track pT > 2.5 GeV

      fBestNumberOfIterations[0][5] = 4;  // Centrality = 0-10, track pT > 3 GeV
      fBestNumberOfIterations[1][5] = 4;  // Centrality = 10-30, track pT > 3 GeV
      fBestNumberOfIterations[2][5] = 3;  // Centrality = 30-50, track pT > 3 GeV
      fBestNumberOfIterations[3][5] = 3;  // Centrality = 50-90, track pT > 3 GeV

    } else {
      std::cout << "EECUnfoldConfiguration::ERROR! " << fSystematicIndex  << " is undefined index for unfolding parameter set!" << std::endl;
      std::cout << "Yout code will crash soon if you are reading the response matrix file name from EECUnfoldConfiguration." << std::endl;
      std::cout << "If not, default number 4 is used in all cases. Please check your code and define a number in range 0-" << kNParameterSets-1 << std::endl;
    }

  // The latest configuration for pp
  } else {

    if(fSystematicIndex == kDefault){
      // Configuration for default response matrix
      // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
      // Input: chi2Histograms_pp_split1_2023-06-05.root
      //        chi2Histograms_pp_split2_2023-06-05.root

      fResponseMatrixFileName = Form("data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_32deltaRBins_responseMatrix%s_processed_2023-06-02.root", splitName[fSplitIndex]);

      fBestNumberOfIterations[0][3] = 6;  // track pT > 2 GeV
      fBestNumberOfIterations[0][4] = 6;  // track pT > 2.5 GeV
      fBestNumberOfIterations[0][5] = 6;  // track pT > 3 GeV

    } else if(fSystematicIndex == kJetPtResolutionUncertainty){
      // Configuration for jet pT resolution uncertainty evaluation
      // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
      // Input: chi2Histograms_pp_split1_systematicForJetPtResolution_2023-06-05.root
      //        chi2Histograms_pp_split2_systematicForJetPtResolution_2023-06-05.root

      fResponseMatrixFileName = Form("data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_32deltaRBins_smearJetPtResolution_responseMatrix%s_processed_2023-06-02.root", splitName[fSplitIndex]);

      fBestNumberOfIterations[0][3] = 10;  // track pT > 2 GeV
      fBestNumberOfIterations[0][4] = 10;  // track pT > 2.5 GeV
      fBestNumberOfIterations[0][5] = 10;  // track pT > 3 GeV

    } else if(fSystematicIndex == kJetEnergyScaleUncertainty){
      // Configuration for jet energy scale uncertainty evaluation
      // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
      // Input: chi2Histograms_pp_split1_systematicForJetPtUncertainty_2023-06-05.root
      //        chi2Histograms_pp_split2_systematicForJetPtUncertainty_2023-06-05.root

      fResponseMatrixFileName = Form("data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_32deltaRBins_smearJetPtUncertainty_responseMatrix%s_processed_2023-06-02.root", splitName[fSplitIndex]);

      fBestNumberOfIterations[0][3] = 6;  // track pT > 2 GeV
      fBestNumberOfIterations[0][4] = 6;  // track pT > 2.5 GeV
      fBestNumberOfIterations[0][5] = 6;  // track pT > 3 GeV

    } else if(fSystematicIndex == kJetPtPriorUncertainty){
      // Configuration for jet energy scale uncertainty evaluation
      // Macro from which the numbers are determined: drawUnfoldingChi2Test.C
      // Input: chi2Histograms_pp_split1_jetPtWeight_2023-06-16.root
      //        chi2Histograms_pp_split2_jetPtWeight_2023-06-16.root

      fResponseMatrixFileName = Form("data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_32deltaRBins_jetPtWeight_responseMatrix%s_processed_2023-06-15.root", splitName[fSplitIndex]);

      fBestNumberOfIterations[0][3] = 6;  // track pT > 2 GeV
      fBestNumberOfIterations[0][4] = 6;  // track pT > 2.5 GeV
      fBestNumberOfIterations[0][5] = 6;  // track pT > 3 GeV

    } else {
      std::cout << "EECUnfoldConfiguration::ERROR! " << fSystematicIndex << " is undefined index for unfolding parameter set!" << std::endl;
      std::cout << "Yout code will crash soon if you are reading the response matrix file name from EECUnfoldConfiguration." << std::endl;
      std::cout << "If not, default number 4 is used in all cases. Please check your code and define a number in range 0-" << kNParameterSets-1 << std::endl;
    }

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
