#include "EECCard.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"
#include "AlgorithmLibrary.h"
#include "EECUnfoldConfiguration.h"

/*
 * Draw the chi2 test results to determine a good number of iterations for unfolding
 */
void drawUnfoldingChi2Test(){

  // **********************************
  //       Open the input files
  // **********************************

  const int nInputFiles = 2;
  TString inputFileName[nInputFiles];
  inputFileName[0] = "chi2Files/chi2Histograms_PbPb_nominalEnergyWeight_optimizedUnfoldingBinning_split1_nominalSmear_4pCentShift_2024-01-17.root";
  inputFileName[1] = "chi2Files/chi2Histograms_PbPb_nominalEnergyWeight_optimizedUnfoldingBinning_split2_nominalSmear_4pCentShift_2024-01-17.root";
  // chi2Histograms_pp_split1_2023-06-05.root
  // chi2Histograms_pp_split2_2023-06-05.root
  // chi2Histograms_PbPb_energyWeightSquared_split1_nominalSmear_4pCentShift_2023-10-26.root
  // chi2Histograms_pp_energyWeightSquared_split1_nominalSmear_2023-10-31.root
  // chi2Histograms_PbPb_energyWeightSquared_optimizedUnfoldingBinning_split1_nominalSmear_4pCentShift_2024-01-11.root
  // chi2Histograms_PbPb_nominalEnergyWeight_optimizedUnfoldingBinning_split1_nominalSmear_4pCentShift_2024-01-17.root

  // Open the input files and read the card containing binning information
  TFile* inputFile[nInputFiles];
  EECCard* unfoldingCard[nInputFiles];
  for(int iFile = 0; iFile < nInputFiles; iFile++){
    inputFile[iFile] = TFile::Open(inputFileName[iFile]);

    if(inputFile[iFile] == NULL) {
      cout << "Error! The file " << inputFileName[iFile].Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }

    unfoldingCard[iFile] = new EECCard(inputFile[iFile]);
  }


  // Determine if we are dealing with pp or PbPb data
  TString collisionSystem = unfoldingCard[0]->GetDataType();
  bool isPbPbData = collisionSystem.Contains("PbPb");
  
  // **********************************************************************************
  //       Check that the centrality, jet pT and track pT bins match with the files
  // **********************************************************************************

  for(int iFile = 1; iFile < nInputFiles; iFile++) {

    if(isPbPbData != (unfoldingCard[iFile]->GetDataType().Contains("PbPb"))) {
      cout << "You are trying to compare chi2 histograms between pp with PbPb or vice versa!" << endl;
      cout << "Please ensure that the chi2 histograms are from the same system!" << endl;
    }

    // Only check compatible centrality bins for PbPb
    if(isPbPbData) {
      if(unfoldingCard[0]->GetNCentralityBins() != unfoldingCard[iFile]->GetNCentralityBins()) {
        cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
        return;
      }
      for(int iCentrality = 0; iCentrality < unfoldingCard[0]->GetNCentralityBins(); iCentrality++) {
        if(TMath::Abs(unfoldingCard[0]->GetLowBinBorderCentrality(iCentrality) - unfoldingCard[iFile]->GetLowBinBorderCentrality(iCentrality)) > 5) {
          cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
          return;
        }
        if(TMath::Abs(unfoldingCard[0]->GetHighBinBorderCentrality(iCentrality) - unfoldingCard[iFile]->GetHighBinBorderCentrality(iCentrality)) > 5) {
          cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
          return;
        }
      }
    }  // Checking compatible centrality bins

    // Check that track pT bins are compatible
    if(unfoldingCard[0]->GetNTrackPtBinsEEC() != unfoldingCard[iFile]->GetNTrackPtBinsEEC()) {
      cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
    for(int iTrackPt = 0; iTrackPt < unfoldingCard[0]->GetNTrackPtBinsEEC(); iTrackPt++) {
      if(TMath::Abs(unfoldingCard[0]->GetLowBinBorderTrackPtEEC(iTrackPt) - unfoldingCard[iFile]->GetLowBinBorderTrackPtEEC(iTrackPt)) > 0.01) {
        cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
        return;
      }
      if(TMath::Abs(unfoldingCard[0]->GetHighBinBorderTrackPtEEC(iTrackPt) - unfoldingCard[iFile]->GetHighBinBorderTrackPtEEC(iTrackPt)) > 0.01) {
        cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
        return;
      }
    }  // Checking the compatibility of track pT bins

    // Check that jet pT bins are compatible
    if(unfoldingCard[0]->GetNJetPtBinsEEC() != unfoldingCard[iFile]->GetNJetPtBinsEEC()) {
      cout << "Error! Jet pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
    for(int iTrackPt = 0; iTrackPt < unfoldingCard[0]->GetNJetPtBinsEEC(); iTrackPt++) {
      if(TMath::Abs(unfoldingCard[0]->GetLowBinBorderJetPtEEC(iTrackPt) - unfoldingCard[iFile]->GetLowBinBorderJetPtEEC(iTrackPt)) > 0.01) {
        cout << "Error! Jet pT bins do not match between the two files! Cannot execute the code!" << endl;
        return;
      }
      if(TMath::Abs(unfoldingCard[0]->GetHighBinBorderJetPtEEC(iTrackPt) - unfoldingCard[iFile]->GetHighBinBorderJetPtEEC(iTrackPt)) > 0.01) {
        cout << "Error! Jet pT bins do not match between the two files! Cannot execute the code!" << endl;
        return;
      }
    }  // Checking the compatibility of track pT bins
  } // Loop to check all the studied files have the same binning information

  // ********************************************************
  //       Binning configuration for the unfolding study
  // ********************************************************
  
  // Find the number of bins from the card
  const int nCentralityBins = (isPbPbData) ? unfoldingCard[0]->GetNCentralityBins() : 1;
  const int nTrackPtBinsEEC = unfoldingCard[0]->GetNTrackPtBinsEEC();
  const int nJetPtBins = unfoldingCard[0]->GetNJetPtBinsEEC();

  // Bin range to be studied
  int firstStudiedCentralityBin = 0;
  int lastStudiedCentralityBin = nCentralityBins-1;
  
  int firstStudiedTrackPtBinEEC = 3;
  int lastStudiedTrackPtBinEEC = 5;

  int firstStudiedJetPtBin = 0;
  int lastStudiedJetPtBin = 3;

  const bool drawChi2map = true;                      // Draw the chi2 values for individual jet pT bins
  const bool drawChi2combined = false;                 // Draw single good chi2 value for each response matrix determined from relevent region
  const bool drawChi2mapForwardFolded = false;          // Draw the chi2 values for individual jet pT bins from forward folded distributions
  const bool drawChi2combinedForwardFolded = false;    // Draw single good chi2 value for each response matrix determined from relevent region from forward folded distributions
  const bool drawUnfoldedToTruthComparison = false;    // Compare unfolded distributions to truth
  const bool drawForwardFoldedToMeasuredComparison = false;  // Compare measured distribution to forward folded distribution
  const bool drawBestIterationRatioComparison = false; // Draw unfolded to truth ratios for the selected number of iterations
  const bool oneIterationPerMatrix = false;            // If drawing best iteration ratio, use single iteration number for each matrix 

  int bestNumberOfIterations[nCentralityBins][nJetPtBins][nTrackPtBinsEEC];
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){
        bestNumberOfIterations[iCentrality][iJetPt][iTrackPt] = 4;
      }
    }
  }

  // String for number of iterations
  TString iterationString = "NIT iter";

  // Create unfolding configuration object to read the best number of iterations when one iteration per matrix is used
  EECUnfoldConfiguration* grandIterationOracle = new EECUnfoldConfiguration(unfoldingCard[0], 0, 0, 2);

  // Set manually the best number of iteration determined totally objectively
  if(isPbPbData){
    bestNumberOfIterations[0][0][3] = 7;  // Centrality = 0-10, track pT > 2 GeV, 120 < jet pT < 140
    bestNumberOfIterations[0][1][3] = 3;  // Centrality = 0-10, track pT > 2 GeV, 140 < jet pT < 160
    bestNumberOfIterations[0][2][3] = 2;  // Centrality = 0-10, track pT > 2 GeV, 160 < jet pT < 180
    bestNumberOfIterations[0][3][3] = 2;  // Centrality = 0-10, track pT > 2 GeV, 180 < jet pT < 200
    bestNumberOfIterations[1][0][3] = 6;  // Centrality = 10-30, track pT > 2 GeV, 120 < jet pT < 140
    bestNumberOfIterations[1][1][3] = 3;  // Centrality = 10-30, track pT > 2 GeV, 140 < jet pT < 160
    bestNumberOfIterations[1][2][3] = 3;  // Centrality = 10-30, track pT > 2 GeV, 160 < jet pT < 180
    bestNumberOfIterations[1][3][3] = 3;  // Centrality = 10-30, track pT > 2 GeV, 180 < jet pT < 200
    bestNumberOfIterations[2][0][3] = 3;  // Centrality = 30-50, track pT > 2 GeV, 120 < jet pT < 140
    bestNumberOfIterations[2][1][3] = 3;  // Centrality = 30-50, track pT > 2 GeV, 140 < jet pT < 160
    bestNumberOfIterations[2][2][3] = 3;  // Centrality = 30-50, track pT > 2 GeV, 160 < jet pT < 180
    bestNumberOfIterations[2][3][3] = 5;  // Centrality = 30-50, track pT > 2 GeV, 180 < jet pT < 200
    bestNumberOfIterations[3][0][3] = 4;  // Centrality = 50-90, track pT > 2 GeV, 120 < jet pT < 140
    bestNumberOfIterations[3][1][3] = 3;  // Centrality = 50-90, track pT > 2 GeV, 140 < jet pT < 160
    bestNumberOfIterations[3][2][3] = 3;  // Centrality = 50-90, track pT > 2 GeV, 160 < jet pT < 180
    bestNumberOfIterations[3][3][3] = 4;  // Centrality = 50-90, track pT > 2 GeV, 180 < jet pT < 200

    bestNumberOfIterations[0][0][4] = 8;  // Centrality = 0-10, track pT > 2.5 GeV, 120 < jet pT < 140
    bestNumberOfIterations[0][1][4] = 4;  // Centrality = 0-10, track pT > 2.5 GeV, 140 < jet pT < 160
    bestNumberOfIterations[0][2][4] = 3;  // Centrality = 0-10, track pT > 2.5 GeV, 160 < jet pT < 180
    bestNumberOfIterations[0][3][4] = 3;  // Centrality = 0-10, track pT > 2.5 GeV, 180 < jet pT < 200
    bestNumberOfIterations[1][0][4] = 7;  // Centrality = 10-30, track pT > 2.5 GeV, 120 < jet pT < 140
    bestNumberOfIterations[1][1][4] = 5;  // Centrality = 10-30, track pT > 2.5 GeV, 140 < jet pT < 160
    bestNumberOfIterations[1][2][4] = 3;  // Centrality = 10-30, track pT > 2.5 GeV, 160 < jet pT < 180
    bestNumberOfIterations[1][3][4] = 3;  // Centrality = 10-30, track pT > 2.5 GeV, 180 < jet pT < 200
    bestNumberOfIterations[2][0][4] = 5;  // Centrality = 30-50, track pT > 2.5 GeV, 120 < jet pT < 140
    bestNumberOfIterations[2][1][4] = 3;  // Centrality = 30-50, track pT > 2.5 GeV, 140 < jet pT < 160
    bestNumberOfIterations[2][2][4] = 3;  // Centrality = 30-50, track pT > 2.5 GeV, 160 < jet pT < 180
    bestNumberOfIterations[2][3][4] = 5;  // Centrality = 30-50, track pT > 2.5 GeV, 180 < jet pT < 200
    bestNumberOfIterations[3][0][4] = 4;  // Centrality = 50-90, track pT > 2.5 GeV, 120 < jet pT < 140
    bestNumberOfIterations[3][1][4] = 3;  // Centrality = 50-90, track pT > 2.5 GeV, 140 < jet pT < 160
    bestNumberOfIterations[3][2][4] = 3;  // Centrality = 50-90, track pT > 2.5 GeV, 160 < jet pT < 180
    bestNumberOfIterations[3][3][4] = 3;  // Centrality = 50-90, track pT > 2.5 GeV, 180 < jet pT < 200

    bestNumberOfIterations[0][0][5] = 10;  // Centrality = 0-10, track pT > 3 GeV, 120 < jet pT < 140
    bestNumberOfIterations[0][1][5] = 5;  // Centrality = 0-10, track pT > 3 GeV, 140 < jet pT < 160
    bestNumberOfIterations[0][2][5] = 3;  // Centrality = 0-10, track pT > 3 GeV, 160 < jet pT < 180
    bestNumberOfIterations[0][3][5] = 3;  // Centrality = 0-10, track pT > 3 GeV, 180 < jet pT < 200
    bestNumberOfIterations[1][0][5] = 7;  // Centrality = 10-30, track pT > 3 GeV, 120 < jet pT < 140
    bestNumberOfIterations[1][1][5] = 4;  // Centrality = 10-30, track pT > 3 GeV, 140 < jet pT < 160
    bestNumberOfIterations[1][2][5] = 3;  // Centrality = 10-30, track pT > 3 GeV, 160 < jet pT < 180
    bestNumberOfIterations[1][3][5] = 3;  // Centrality = 10-30, track pT > 3 GeV, 180 < jet pT < 200
    bestNumberOfIterations[2][0][5] = 5;  // Centrality = 30-50, track pT > 3 GeV, 120 < jet pT < 140
    bestNumberOfIterations[2][1][5] = 3;  // Centrality = 30-50, track pT > 3 GeV, 140 < jet pT < 160
    bestNumberOfIterations[2][2][5] = 3;  // Centrality = 30-50, track pT > 3 GeV, 160 < jet pT < 180
    bestNumberOfIterations[2][3][5] = 5;  // Centrality = 30-50, track pT > 3 GeV, 180 < jet pT < 200
    bestNumberOfIterations[3][0][5] = 4;  // Centrality = 50-90, track pT > 3 GeV, 120 < jet pT < 140
    bestNumberOfIterations[3][1][5] = 3;  // Centrality = 50-90, track pT > 3 GeV, 140 < jet pT < 160
    bestNumberOfIterations[3][2][5] = 3;  // Centrality = 50-90, track pT > 3 GeV, 160 < jet pT < 180
    bestNumberOfIterations[3][3][5] = 3;  // Centrality = 50-90, track pT > 3 GeV, 180 < jet pT < 200

  } else {    
    bestNumberOfIterations[0][0][5] = 11;  // pp, track pT > 3 GeV, 120 < jet pT < 140
    bestNumberOfIterations[0][1][5] = 4;  // pp, track pT > 3 GeV, 140 < jet pT < 160
    bestNumberOfIterations[0][2][5] = 10;  // pp, track pT > 3 GeV, 160 < jet pT < 180
    bestNumberOfIterations[0][3][5] = 6;  // pp, track pT > 3 GeV, 180 < jet pT < 200

  }

  // Best number of iteratios determined for combining all jet pT bins
  if(oneIterationPerMatrix) {
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){
          bestNumberOfIterations[iCentrality][iJetPt][iTrackPt] = grandIterationOracle->GetNumberOfIterations(unfoldingCard[0]->GetBinBordersCentrality(iCentrality), unfoldingCard[0]->GetLowBinBorderTrackPtEEC(iTrackPt));
        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
  }

  bool saveFigures = true;
  TString saveComment = "";
  TString figureFormat = "pdf";

  if(unfoldingCard[0]->GetWeightExponent() == 1){
    saveComment.Prepend("_nominalEnergyWeight");
  } else {
    saveComment.Prepend("_energyWeightSquared");
  }

  // Histograms for chi2 and error2 map
  TH1D* hChi2map[nCentralityBins][nJetPtBins][nTrackPtBinsEEC][nInputFiles+1];
  TH1D* hChi2combined[nCentralityBins][nTrackPtBinsEEC][nInputFiles+1];
  TH1D* hChi2mapForwardFolded[nCentralityBins][nJetPtBins][nTrackPtBinsEEC][nInputFiles+1];
  TH1D* hChi2combinedForwardFolded[nCentralityBins][nTrackPtBinsEEC][nInputFiles+1];
  TH1D* hErrorMap[nCentralityBins][nJetPtBins][nTrackPtBinsEEC][nInputFiles+1];

  // Initialize all the histograms to null
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
      for(int iFile = 0; iFile < nInputFiles+1; iFile++){
        for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){
          hChi2map[iCentrality][iJetPt][iTrackPt][iFile] = NULL;
          hChi2mapForwardFolded[iCentrality][iJetPt][iTrackPt][iFile] = NULL;
          hErrorMap[iCentrality][iJetPt][iTrackPt][iFile] = NULL;
        } // Jet pT loop
        hChi2combined[iCentrality][iTrackPt][iFile] = NULL;
        hChi2combinedForwardFolded[iCentrality][iTrackPt][iFile] = NULL;
      } // File loop
    } // Track pT loop
  }  // Centrality loop

  // Read the selected histograms from the input files
  for(int iFile = 0; iFile < nInputFiles; iFile++){
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
        for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){
          hChi2map[iCentrality][iJetPt][iTrackPt][iFile] = (TH1D*) inputFile[iFile]->Get(Form("chi2map/chi2map%d%d%d", iCentrality, iTrackPt, iJetPt));
          hChi2mapForwardFolded[iCentrality][iJetPt][iTrackPt][iFile] = (TH1D*) inputFile[iFile]->Get(Form("chi2mapForwardFold/chi2mapForwardFold%d%d%d", iCentrality, iTrackPt, iJetPt));
          hErrorMap[iCentrality][iJetPt][iTrackPt][iFile] = (TH1D*) inputFile[iFile]->Get(Form("errorMap/errorMap%d%d%d", iCentrality, iTrackPt, iJetPt));
        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
  } // File loop

  // Add all the files together for a combined result
  if(drawChi2combined || drawChi2combinedForwardFolded || drawChi2map || drawChi2mapForwardFolded){
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
        for(int iJetPt = firstStudiedJetPtBin; iJetPt <= lastStudiedJetPtBin; iJetPt++){
          hChi2map[iCentrality][iJetPt][iTrackPt][nInputFiles] = (TH1D*) hChi2map[iCentrality][iJetPt][iTrackPt][0]->Clone(Form("averageChi2map%d%d%d", iCentrality, iTrackPt, iJetPt));
          hChi2mapForwardFolded[iCentrality][iJetPt][iTrackPt][nInputFiles] = (TH1D*) hChi2mapForwardFolded[iCentrality][iJetPt][iTrackPt][0]->Clone(Form("averageChi2mapForwardFolded%d%d%d", iCentrality, iTrackPt, iJetPt));
          for(int iFile = 1; iFile < nInputFiles; iFile++){
            hChi2map[iCentrality][iJetPt][iTrackPt][nInputFiles]->Add(hChi2map[iCentrality][iJetPt][iTrackPt][iFile]);
            hChi2mapForwardFolded[iCentrality][iJetPt][iTrackPt][nInputFiles]->Add(hChi2mapForwardFolded[iCentrality][iJetPt][iTrackPt][iFile]);
          }
          for(int iBin = 1; iBin <= hChi2map[iCentrality][iJetPt][iTrackPt][nInputFiles]->GetNbinsX(); iBin++){
            hChi2map[iCentrality][iJetPt][iTrackPt][nInputFiles]->SetBinContent(iBin, hChi2map[iCentrality][iJetPt][iTrackPt][nInputFiles]->GetBinContent(iBin) / nInputFiles);
            hChi2mapForwardFolded[iCentrality][iJetPt][iTrackPt][nInputFiles]->SetBinContent(iBin, hChi2mapForwardFolded[iCentrality][iJetPt][iTrackPt][nInputFiles]->GetBinContent(iBin) / nInputFiles);
          }
        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop

    // Combine the chi2 values from each response matrix to a single grand chi2 value
    for(int iFile = 0; iFile < nInputFiles+1; iFile++){
      for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
          hChi2combined[iCentrality][iTrackPt][iFile] = (TH1D*) hChi2map[iCentrality][firstStudiedJetPtBin][iTrackPt][iFile]->Clone(Form("combinedChi2%d%d%d", iCentrality, iTrackPt, iFile));
          hChi2combinedForwardFolded[iCentrality][iTrackPt][iFile] = (TH1D*) hChi2mapForwardFolded[iCentrality][firstStudiedJetPtBin][iTrackPt][iFile]->Clone(Form("combinedChi2ForwardFolded%d%d%d", iCentrality, iTrackPt, iFile));
          for(int iJetPt = firstStudiedJetPtBin+1; iJetPt <= lastStudiedJetPtBin; iJetPt++){
            hChi2combined[iCentrality][iTrackPt][iFile]->Add(hChi2map[iCentrality][iJetPt][iTrackPt][iFile]);
            hChi2combinedForwardFolded[iCentrality][iTrackPt][iFile]->Add(hChi2mapForwardFolded[iCentrality][iJetPt][iTrackPt][iFile]);
          }
        } // Track pT loop
      } // Centrality loop
    } // File loop
  }

  // Determine the number of iterations from the chi2 histograms
  const int nIterations = hChi2map[firstStudiedCentralityBin][0][firstStudiedTrackPtBinEEC][0]->GetNbinsX();

  // Histograms for energy-energy correlator distributions and their ratios
  TH1D* hMeasured[nCentralityBins][nJetPtBins][nTrackPtBinsEEC][nInputFiles];
  TH1D* hTruth[nCentralityBins][nJetPtBins][nTrackPtBinsEEC][nInputFiles];
  TH1D* hUnfolded[nCentralityBins][nJetPtBins][nTrackPtBinsEEC][nIterations][nInputFiles];
  TH1D* hMeasuredToTruthRatio[nCentralityBins][nJetPtBins][nTrackPtBinsEEC][nInputFiles];
  TH1D* hUnfoldedToTruthRatio[nCentralityBins][nJetPtBins][nTrackPtBinsEEC][nIterations][nInputFiles];
  TH1D* hForwardFoldedToMeasuredRatio[nCentralityBins][nJetPtBins][nTrackPtBinsEEC][nIterations][nInputFiles];

  // Initialize the energy-energy correlator distributions and their ratios to NULL
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){
        for(int iFile = 0; iFile < nInputFiles; iFile++){
          hMeasured[iCentrality][iJetPt][iTrackPt][iFile] = NULL;
          hTruth[iCentrality][iJetPt][iTrackPt][iFile] = NULL;
          hMeasuredToTruthRatio[iCentrality][iJetPt][iTrackPt][iFile] = NULL;
          for(int iIteration = 0; iIteration < nIterations; iIteration++){
            hUnfolded[iCentrality][iJetPt][iTrackPt][iIteration][iFile] = NULL;
            hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][iIteration][iFile] = NULL;
            hForwardFoldedToMeasuredRatio[iCentrality][iJetPt][iTrackPt][iIteration][iFile] = NULL;
          }
        } // File loop
      } // Jet pT loop
    } // Track pT loop
  }  // Centrality loop

  // Read the energy-energy correlator histograms and ratios from the input file
  for(int iFile = 0; iFile < nInputFiles; iFile++) {
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++) {
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++) {
        for(int iJetPt = firstStudiedJetPtBin; iJetPt <= lastStudiedJetPtBin; iJetPt++) {
          hMeasured[iCentrality][iJetPt][iTrackPt][iFile] = (TH1D*) inputFile[iFile]->Get(Form("measured/hMeasured%d%d%d", iCentrality, iTrackPt, iJetPt));
          hTruth[iCentrality][iJetPt][iTrackPt][iFile] = (TH1D*) inputFile[iFile]->Get(Form("truth/hTruth%d%d%d", iCentrality, iTrackPt, iJetPt));
          hMeasuredToTruthRatio[iCentrality][iJetPt][iTrackPt][iFile] = (TH1D*) inputFile[iFile]->Get(Form("measuredToTruthRatio/measuredToTruthRatio%d%d%d", iCentrality, iTrackPt, iJetPt));
          for(int iIteration = 0; iIteration < nIterations; iIteration++) {
            hUnfolded[iCentrality][iJetPt][iTrackPt][iIteration][iFile] = (TH1D*) inputFile[iFile]->Get(Form("unfolded/hUnfolded%d%d%d%d", iCentrality, iTrackPt, iJetPt, iIteration));
            hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][iIteration][iFile] = (TH1D*) inputFile[iFile]->Get(Form("unfoldedToTruthRatio/unfoldedToTruthRatio%d%d%d%d", iCentrality, iTrackPt, iJetPt, iIteration));
            hForwardFoldedToMeasuredRatio[iCentrality][iJetPt][iTrackPt][iIteration][iFile] = (TH1D*) inputFile[iFile]->Get(Form("forwardFoldedToMeasuredRatio/unfoldedToTruthRatio%d%d%d%d", iCentrality, iTrackPt, iJetPt, iIteration));
          } // Iteration loop
        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
  } // File loop

  // ******************************************************************
  //         Draw different comparisons to check the performance
  // ******************************************************************

  // Prepare a JDrawer for drawing purposes
  JDrawer* drawer = new JDrawer();
  drawer->SetRelativeCanvasSize(1.1,1.2);
  drawer->SetLeftMargin(0.15);
  drawer->SetTopMargin(0.07);
  drawer->SetRightMargin(0.05);
  drawer->SetTitleOffsetY(1.2);

  // Helper variables
  TLegend* legend;
  TLegend* lowLegend;
  const int nIterationLegends = 4;
  TLegend* iterationLegend[nIterationLegends];
  TString centralityString;
  TString compactCentralityString;
  TString jetPtString;
  TString compactJetPtString;
  TString trackPtString;
  TString compactTrackPtString;

  int fileColor[] = {kBlack, kRed, kBlue, kGreen+3, kMagenta, kCyan, kOrange, kViolet+3, kPink-7, kSpring+3, kAzure-7};
  int fileMarker[] = {kFullSquare, kFullCircle, kFullDiamond, kFullCross};
  std::pair<double,double> histogramYrange;
  AlgorithmLibrary* methods = new AlgorithmLibrary();
  double normalizationRegionLow = 0.008;
  double normalizationRegionHigh = 0.39;

  // Automatically extract colors for iterations:
  gStyle->SetPalette(kBird);
  auto niceColors =  TColor::GetPalette();
  double step = 254.0/nIterations;
  int iterationColor[nIterations];
  for(int iIteration = 0; iIteration < nIterations; iIteration++){
    iterationColor[iIteration] = niceColors.At(TMath::Ceil(iIteration*step));
  }

  TLine* oneLine = new TLine(normalizationRegionLow, 1, normalizationRegionHigh, 1);
  oneLine->SetLineStyle(2);

  if(drawChi2map){

    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){

      // Set the centrality information for legends and figure saving
      if(isPbPbData) {
        centralityString = Form("Pythia+Hydjet: %.0f-%.0f", unfoldingCard[0]->GetLowBinBorderCentrality(iCentrality), unfoldingCard[0]->GetHighBinBorderCentrality(iCentrality));
        compactCentralityString = Form("_C%.0f-%.0f", unfoldingCard[0]->GetLowBinBorderCentrality(iCentrality), unfoldingCard[0]->GetHighBinBorderCentrality(iCentrality));
      } else {
        centralityString = "Pythia8";
        compactCentralityString = "_pythia8";
      }

      for(int iJetPt = firstStudiedJetPtBin; iJetPt <= lastStudiedJetPtBin; iJetPt++){

        // Set the jet pT information for legends and figure saving
        jetPtString = Form("%.0f < jet p_{T} < %.0f", unfoldingCard[0]->GetLowBinBorderJetPtEEC(iJetPt), unfoldingCard[0]->GetHighBinBorderJetPtEEC(iJetPt));
        compactJetPtString = Form("_J=%.0f-%.0f", unfoldingCard[0]->GetLowBinBorderJetPtEEC(iJetPt), unfoldingCard[0]->GetHighBinBorderJetPtEEC(iJetPt));

        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

          // Set the track pT information for legends and figure saving
          trackPtString = Form("%.1f < track p_{T}", unfoldingCard[0]->GetLowBinBorderTrackPtEEC(iTrackPt));
          compactTrackPtString = Form("_T%.1f", unfoldingCard[0]->GetLowBinBorderTrackPtEEC(iTrackPt));
          compactTrackPtString.ReplaceAll(".","v");

          // Find good y-ranges for plotting
          histogramYrange = std::make_pair(10e10, 0);
          for(int iFile = 0; iFile < nInputFiles+1; iFile++){
            histogramYrange = methods->FindHistogramMinMax(hChi2map[iCentrality][iJetPt][iTrackPt][iFile], histogramYrange);
          }
          // Add some empty space to the top and the bottom of the histogram
          histogramYrange.first = TMath::Max(0.0, histogramYrange.first - histogramYrange.second*0.05);
          histogramYrange.second = histogramYrange.second + histogramYrange.second*0.06;

          // Draw the chi2 map to the canvas
          hChi2map[iCentrality][iJetPt][iTrackPt][0]->GetYaxis()->SetRangeUser(histogramYrange.first, histogramYrange.second);
          hChi2map[iCentrality][iJetPt][iTrackPt][0]->SetMarkerStyle(fileMarker[0]);
          hChi2map[iCentrality][iJetPt][iTrackPt][0]->SetMarkerColor(fileColor[0]);
          drawer->DrawHistogram(hChi2map[iCentrality][iJetPt][iTrackPt][0], "Number of iterations", "#chi^{2}", " ", "p");

          for(int iFile = 1; iFile < nInputFiles+1; iFile++){
            hChi2map[iCentrality][iJetPt][iTrackPt][iFile]->SetMarkerStyle(fileMarker[iFile]);
            hChi2map[iCentrality][iJetPt][iTrackPt][iFile]->SetMarkerColor(fileColor[iFile]);
            hChi2map[iCentrality][iJetPt][iTrackPt][iFile]->Draw("p,same");
          }

          // Add a legend to the figure
          legend = new TLegend(0.18, 0.54, 0.43, 0.87);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);

          legend->AddEntry((TObject*)0, centralityString.Data(), "");
          legend->AddEntry((TObject*)0, jetPtString.Data(), "");
          legend->AddEntry((TObject*)0, trackPtString.Data(), "");

          for(int iFile = 0; iFile < nInputFiles; iFile++){
            legend->AddEntry(hChi2map[iCentrality][iJetPt][iTrackPt][iFile], Form("Split %d", iFile+1), "p");
          }
          legend->AddEntry(hChi2map[iCentrality][iJetPt][iTrackPt][nInputFiles], Form("Split average"), "p");

          legend->Draw();

          // Save the figures to a file
          if(saveFigures) {
            gPad->GetCanvas()->SaveAs(Form("figures/chi2map%s%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
          }

        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Drawing chi2 map

  // Draw the chi2 map from forward folded distributions
  if(drawChi2mapForwardFolded){

    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){

      // Set the centrality information for legends and figure saving
      if(isPbPbData) {
        centralityString = Form("Pythia+Hydjet: %.0f-%.0f", unfoldingCard[0]->GetLowBinBorderCentrality(iCentrality), unfoldingCard[0]->GetHighBinBorderCentrality(iCentrality));
        compactCentralityString = Form("_C%.0f-%.0f", unfoldingCard[0]->GetLowBinBorderCentrality(iCentrality), unfoldingCard[0]->GetHighBinBorderCentrality(iCentrality));
      } else {
        centralityString = "Pythia8";
        compactCentralityString = "_pythia8";
      }

      for(int iJetPt = firstStudiedJetPtBin; iJetPt <= lastStudiedJetPtBin; iJetPt++){

        // Set the jet pT information for legends and figure saving
        jetPtString = Form("%.0f < jet p_{T} < %.0f", unfoldingCard[0]->GetLowBinBorderJetPtEEC(iJetPt), unfoldingCard[0]->GetHighBinBorderJetPtEEC(iJetPt));
        compactJetPtString = Form("_J=%.0f-%.0f", unfoldingCard[0]->GetLowBinBorderJetPtEEC(iJetPt), unfoldingCard[0]->GetHighBinBorderJetPtEEC(iJetPt));

        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

          // Set the track pT information for legends and figure saving
          trackPtString = Form("%.1f < track p_{T}", unfoldingCard[0]->GetLowBinBorderTrackPtEEC(iTrackPt));
          compactTrackPtString = Form("_T%.0f", unfoldingCard[0]->GetLowBinBorderTrackPtEEC(iTrackPt));

          // Find good y-ranges for plotting
          histogramYrange = std::make_pair(10e10, 0);
          for(int iFile = 0; iFile < nInputFiles+1; iFile++){
            histogramYrange = methods->FindHistogramMinMax(hChi2map[iCentrality][iJetPt][iTrackPt][iFile], histogramYrange);
          }
          // Add some empty space to the top and the bottom of the histogram
          histogramYrange.first = TMath::Max(0.0, histogramYrange.first - histogramYrange.second*0.05);
          histogramYrange.second = histogramYrange.second + histogramYrange.second*0.06;

          // Draw the chi2 map to the canvas
          hChi2mapForwardFolded[iCentrality][iJetPt][iTrackPt][0]->GetYaxis()->SetRangeUser(histogramYrange.first, histogramYrange.second);
          hChi2mapForwardFolded[iCentrality][iJetPt][iTrackPt][0]->SetMarkerStyle(fileMarker[0]);
          hChi2mapForwardFolded[iCentrality][iJetPt][iTrackPt][0]->SetMarkerColor(fileColor[0]);
          drawer->DrawHistogram(hChi2mapForwardFolded[iCentrality][iJetPt][iTrackPt][0], "Number of iterations", "#chi^{2}", " ", "p");

          for(int iFile = 1; iFile < nInputFiles+1; iFile++){
            hChi2mapForwardFolded[iCentrality][iJetPt][iTrackPt][iFile]->SetMarkerStyle(fileMarker[iFile]);
            hChi2mapForwardFolded[iCentrality][iJetPt][iTrackPt][iFile]->SetMarkerColor(fileColor[iFile]);
            hChi2mapForwardFolded[iCentrality][iJetPt][iTrackPt][iFile]->Draw("p,same");
          }

          // Add a legend to the figure
          legend = new TLegend(0.18, 0.54, 0.43, 0.87);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);

          legend->AddEntry((TObject*)0, centralityString.Data(), "");
          legend->AddEntry((TObject*)0, jetPtString.Data(), "");
          legend->AddEntry((TObject*)0, trackPtString.Data(), "");

          for(int iFile = 0; iFile < nInputFiles; iFile++){
            legend->AddEntry(hChi2mapForwardFolded[iCentrality][iJetPt][iTrackPt][iFile], Form("Split %d", iFile+1), "p");
          }
          legend->AddEntry(hChi2mapForwardFolded[iCentrality][iJetPt][iTrackPt][nInputFiles], Form("Split average"), "p");

          legend->Draw();

          // Save the figures to a file
          if(saveFigures) {
            gPad->GetCanvas()->SaveAs(Form("figures/chi2mapForwardFolded%s%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
          }

        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Drawing chi2 map

  if(drawChi2combined){

    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){

      // Set the centrality information for legends and figure saving
      if(isPbPbData) {
        centralityString = Form("Pythia+Hydjet: %.0f-%.0f", unfoldingCard[0]->GetLowBinBorderCentrality(iCentrality), unfoldingCard[0]->GetHighBinBorderCentrality(iCentrality));
        compactCentralityString = Form("_C%.0f-%.0f", unfoldingCard[0]->GetLowBinBorderCentrality(iCentrality), unfoldingCard[0]->GetHighBinBorderCentrality(iCentrality));
      } else {
        centralityString = "Pythia8";
        compactCentralityString = "_pythia8";
      }

      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

        // Set the track pT information for legends and figure saving
        trackPtString = Form("%.1f < track p_{T}", unfoldingCard[0]->GetLowBinBorderTrackPtEEC(iTrackPt));
        compactTrackPtString = Form("_T%.1f", unfoldingCard[0]->GetLowBinBorderTrackPtEEC(iTrackPt));
        compactTrackPtString.ReplaceAll(".","v");

        // Find good y-ranges for plotting
        histogramYrange = std::make_pair(10e10, 0);
        for(int iFile = 0; iFile < nInputFiles+1; iFile++){
          histogramYrange = methods->FindHistogramMinMax(hChi2combined[iCentrality][iTrackPt][iFile], histogramYrange);
        }
        // Add some empty space to the top and the bottom of the histogram
        histogramYrange.first = TMath::Max(0.0, histogramYrange.first - histogramYrange.second*0.05);
        histogramYrange.second = histogramYrange.second + histogramYrange.second*0.06;

        // Draw the chi2 map to the canvas
        hChi2combined[iCentrality][iTrackPt][0]->GetYaxis()->SetRangeUser(histogramYrange.first, histogramYrange.second);
        hChi2combined[iCentrality][iTrackPt][0]->SetMarkerStyle(fileMarker[0]);
        hChi2combined[iCentrality][iTrackPt][0]->SetMarkerColor(fileColor[0]);
        drawer->DrawHistogram(hChi2combined[iCentrality][iTrackPt][0], "Number of iterations", "#chi^{2}", " ", "p");

        for(int iFile = 1; iFile < nInputFiles+1; iFile++){
          hChi2combined[iCentrality][iTrackPt][iFile]->SetMarkerStyle(fileMarker[iFile]);
          hChi2combined[iCentrality][iTrackPt][iFile]->SetMarkerColor(fileColor[iFile]);
          hChi2combined[iCentrality][iTrackPt][iFile]->Draw("p,same");
        }

        // Add a legend to the figure
        if(isPbPbData){
          legend = new TLegend(0.18, 0.54, 0.43, 0.87);
        } else {
          legend = new TLegend(0.55, 0.24, 0.8, 0.57);
        }
        legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);

        legend->AddEntry((TObject*)0, centralityString.Data(), "");
        legend->AddEntry((TObject*)0, trackPtString.Data(), "");

        for(int iFile = 0; iFile < nInputFiles; iFile++){
          legend->AddEntry(hChi2combined[iCentrality][iTrackPt][iFile], Form("Split %d", iFile+1), "p");
        }
        legend->AddEntry(hChi2combined[iCentrality][iTrackPt][nInputFiles], Form("Split average"), "p");

        legend->Draw();

        // Save the figures to a file
        if(saveFigures) {
          gPad->GetCanvas()->SaveAs(Form("figures/chi2combined%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
         }

      } // Track pT loop
    } // Centrality loop
  } // Drawing chi2 map

  if(drawChi2combinedForwardFolded){

    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){

      // Set the centrality information for legends and figure saving
      if(isPbPbData) {
        centralityString = Form("Pythia+Hydjet: %.0f-%.0f", unfoldingCard[0]->GetLowBinBorderCentrality(iCentrality), unfoldingCard[0]->GetHighBinBorderCentrality(iCentrality));
        compactCentralityString = Form("_C%.0f-%.0f", unfoldingCard[0]->GetLowBinBorderCentrality(iCentrality), unfoldingCard[0]->GetHighBinBorderCentrality(iCentrality));
      } else {
        centralityString = "Pythia8";
        compactCentralityString = "_pythia8";
      }

      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

        // Set the track pT information for legends and figure saving
        trackPtString = Form("%.1f < track p_{T}", unfoldingCard[0]->GetLowBinBorderTrackPtEEC(iTrackPt));
        compactTrackPtString = Form("_T%.0f", unfoldingCard[0]->GetLowBinBorderTrackPtEEC(iTrackPt));

        // Find good y-ranges for plotting
        histogramYrange = std::make_pair(10e10, 0);
        for(int iFile = 0; iFile < nInputFiles+1; iFile++){
          histogramYrange = methods->FindHistogramMinMax(hChi2combined[iCentrality][iTrackPt][iFile], histogramYrange);
        }
        // Add some empty space to the top and the bottom of the histogram
        histogramYrange.first = TMath::Max(0.0, histogramYrange.first - histogramYrange.second*0.05);
        histogramYrange.second = histogramYrange.second + histogramYrange.second*0.06;

        // Draw the chi2 map to the canvas
        hChi2combinedForwardFolded[iCentrality][iTrackPt][0]->GetYaxis()->SetRangeUser(histogramYrange.first, histogramYrange.second);
        hChi2combinedForwardFolded[iCentrality][iTrackPt][0]->SetMarkerStyle(fileMarker[0]);
        hChi2combinedForwardFolded[iCentrality][iTrackPt][0]->SetMarkerColor(fileColor[0]);
        drawer->DrawHistogram(hChi2combinedForwardFolded[iCentrality][iTrackPt][0], "Number of iterations", "#chi^{2}", " ", "p");

        for(int iFile = 1; iFile < nInputFiles+1; iFile++){
          hChi2combinedForwardFolded[iCentrality][iTrackPt][iFile]->SetMarkerStyle(fileMarker[iFile]);
          hChi2combinedForwardFolded[iCentrality][iTrackPt][iFile]->SetMarkerColor(fileColor[iFile]);
          hChi2combinedForwardFolded[iCentrality][iTrackPt][iFile]->Draw("p,same");
        }

        // Add a legend to the figure
        legend = new TLegend(0.18, 0.54, 0.43, 0.87);
        legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);

        legend->AddEntry((TObject*)0, centralityString.Data(), "");
        legend->AddEntry((TObject*)0, trackPtString.Data(), "");

        for(int iFile = 0; iFile < nInputFiles; iFile++){
          legend->AddEntry(hChi2combinedForwardFolded[iCentrality][iTrackPt][iFile], Form("Split %d", iFile+1), "p");
        }
        legend->AddEntry(hChi2combinedForwardFolded[iCentrality][iTrackPt][nInputFiles], Form("Split average"), "p");

        legend->Draw();

        // Save the figures to a file
        if(saveFigures) {
          gPad->GetCanvas()->SaveAs(Form("figures/chi2combinedForwardFolded%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
         }

      } // Track pT loop
    } // Centrality loop
  } // Drawing chi2 map

  drawer->SetLogX(true);

  if(drawBestIterationRatioComparison){

    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){

      // Set the centrality information for legends and figure saving
      if(isPbPbData) {
        centralityString = Form("Pythia+Hydjet: %.0f-%.0f", unfoldingCard[0]->GetLowBinBorderCentrality(iCentrality), unfoldingCard[0]->GetHighBinBorderCentrality(iCentrality));
        compactCentralityString = Form("_C%.0f-%.0f", unfoldingCard[0]->GetLowBinBorderCentrality(iCentrality), unfoldingCard[0]->GetHighBinBorderCentrality(iCentrality));
      } else {
        centralityString = "Pythia8";
        compactCentralityString = "_pythia8";
      }

      for(int iJetPt = firstStudiedJetPtBin; iJetPt <= lastStudiedJetPtBin; iJetPt++){

        // Set the jet pT information for legends and figure saving
        jetPtString = Form("%.0f < jet p_{T} < %.0f", unfoldingCard[0]->GetLowBinBorderJetPtEEC(iJetPt), unfoldingCard[0]->GetHighBinBorderJetPtEEC(iJetPt));
        compactJetPtString = Form("_J=%.0f-%.0f", unfoldingCard[0]->GetLowBinBorderJetPtEEC(iJetPt), unfoldingCard[0]->GetHighBinBorderJetPtEEC(iJetPt));

        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

          // Set the track pT information for legends and figure saving
          trackPtString = Form("%.1f < track p_{T}", unfoldingCard[0]->GetLowBinBorderTrackPtEEC(iTrackPt));
          compactTrackPtString = Form("_T%.1f", unfoldingCard[0]->GetLowBinBorderTrackPtEEC(iTrackPt));
          compactTrackPtString.ReplaceAll(".","v");

          // Draw the first ratio to canvas
          hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][bestNumberOfIterations[iCentrality][iJetPt][iTrackPt]][0]->GetYaxis()->SetRangeUser(0.8, 1.2);
          hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][bestNumberOfIterations[iCentrality][iJetPt][iTrackPt]][0]->GetXaxis()->SetRangeUser(normalizationRegionLow, normalizationRegionHigh);
          hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][bestNumberOfIterations[iCentrality][iJetPt][iTrackPt]][0]->SetLineColor(fileColor[0]);
          drawer->DrawHistogram(hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][bestNumberOfIterations[iCentrality][iJetPt][iTrackPt]][0], "#Deltar", "Unfolded to truth ratio", " ", "");

          // Add the ratios from other splits to the same figure
          for(int iFile = 1; iFile < nInputFiles; iFile++){
            hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][bestNumberOfIterations[iCentrality][iJetPt][iTrackPt]][iFile]->SetLineColor(fileColor[iFile]);
            hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][bestNumberOfIterations[iCentrality][iJetPt][iTrackPt]][iFile]->Draw("same");
          }

          // Add a legend to the figure
          legend = new TLegend(0.25, 0.73, 0.5, 0.91);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);

          lowLegend = new TLegend(0.25, 0.18, 0.5, 0.36);
          lowLegend->SetFillStyle(0); lowLegend->SetBorderSize(0); lowLegend->SetTextSize(0.05); lowLegend->SetTextFont(62);

          legend->AddEntry((TObject*)0, centralityString.Data(), "");
          legend->AddEntry((TObject*)0, Form("%s, %s", jetPtString.Data(), trackPtString.Data()), "");

          for(int iFile = 0; iFile < nInputFiles; iFile++){
            lowLegend->AddEntry(hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][bestNumberOfIterations[iCentrality][iJetPt][iTrackPt]][iFile], Form("Split %d, %d iterations", iFile+1, bestNumberOfIterations[iCentrality][iJetPt][iTrackPt]), "l");
          }

          legend->Draw();
          lowLegend->Draw();
          
          // Draw a helpful line to one
          oneLine->Draw();

          // Save the figures to a file
          if(saveFigures) {
            gPad->GetCanvas()->SaveAs(Form("figures/unfoldingNumberOfIterationsCheck%s%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
          }
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Best iteration comparison

  if(drawUnfoldedToTruthComparison){

    // In order not to draw too much curves to one plot, treat the two splits separately
    for(int iSplit = 0; iSplit < nInputFiles; iSplit++){

      for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){

         // Set the centrality information for legends and figure saving
        if(isPbPbData) {
          centralityString = Form("Pythia+Hydjet: %.0f-%.0f", unfoldingCard[iSplit]->GetLowBinBorderCentrality(iCentrality), unfoldingCard[iSplit]->GetHighBinBorderCentrality(iCentrality));
          compactCentralityString = Form("_C%.0f-%.0f", unfoldingCard[iSplit]->GetLowBinBorderCentrality(iCentrality), unfoldingCard[iSplit]->GetHighBinBorderCentrality(iCentrality));
        } else {
          centralityString = "Pythia8";
          compactCentralityString = "_pythia8";
        }

        for(int iJetPt = firstStudiedJetPtBin; iJetPt <= lastStudiedJetPtBin; iJetPt++){

          // Set the jet pT information for legends and figure saving
          jetPtString = Form("%.0f < jet p_{T} < %.0f", unfoldingCard[iSplit]->GetLowBinBorderJetPtEEC(iJetPt), unfoldingCard[iSplit]->GetHighBinBorderJetPtEEC(iJetPt));
          compactJetPtString = Form("_J=%.0f-%.0f", unfoldingCard[iSplit]->GetLowBinBorderJetPtEEC(iJetPt), unfoldingCard[iSplit]->GetHighBinBorderJetPtEEC(iJetPt));

          for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

            // Set the track pT information for legends and figure saving
            trackPtString = Form("%.1f < track p_{T}", unfoldingCard[iSplit]->GetLowBinBorderTrackPtEEC(iTrackPt));
            compactTrackPtString = Form("_T%.1f", unfoldingCard[iSplit]->GetLowBinBorderTrackPtEEC(iTrackPt));
            compactTrackPtString.ReplaceAll(".","v");
              
            // Draw the ratios to the canvas
            for(int iIteration = 0; iIteration < nIterations; iIteration++){
              hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][iIteration][iSplit]->SetLineColor(iterationColor[iIteration]);
              if(iIteration == 0){
                hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][iIteration][iSplit]->GetYaxis()->SetRangeUser(0.8, 1.2);
                hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][iIteration][iSplit]->GetXaxis()->SetRangeUser(normalizationRegionLow, normalizationRegionHigh);
                drawer->DrawHistogram(hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][iIteration][iSplit], "#Deltar", "Ratio to truth", " ", "");
              } else {
                hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][iIteration][iSplit]->Draw("same");
              }
            }

            // Add a legend to the figure
            legend = new TLegend(0.2, 0.75, 0.5, 0.9);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);

            legend->AddEntry((TObject*)0, Form("Split %d   %s", iSplit+1, centralityString.Data()), "");
            legend->AddEntry((TObject*)0, Form("%s   %s", jetPtString.Data(), trackPtString.Data()), "");
              

            // Put the iteration numbers to different legend
            for(int iIterationLegend = 0; iIterationLegend < nIterationLegends; iIterationLegend++){
              iterationLegend[iIterationLegend] = new TLegend(0.17+0.19*iIterationLegend, 0.18, 0.4+0.19*iIterationLegend, 0.38);
              iterationLegend[iIterationLegend]->SetFillStyle(0); 
              iterationLegend[iIterationLegend]->SetBorderSize(0); 
              iterationLegend[iIterationLegend]->SetTextSize(0.05); 
              iterationLegend[iIterationLegend]->SetTextFont(62);
            }

            for(int iIteration = 0; iIteration < nIterations; iIteration++){
              iterationLegend[iIteration/(nIterations/nIterationLegends)]->AddEntry(hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][iIteration][iSplit], Form("%s", iterationString.ReplaceAll("NIT",Form("%d",iIteration+1)).Data()), "l");
              iterationString.ReplaceAll(Form("%d",iIteration+1),"NIT");
            }

            oneLine->Draw();
            legend->Draw();
            for(int iIterationLegend = 0; iIterationLegend < nIterationLegends; iIterationLegend++){
              iterationLegend[iIterationLegend]->Draw();
            }

            // Save the figures to a file
            if(saveFigures) {
              gPad->GetCanvas()->SaveAs(Form("figures/unfoldedToTruthRatioCheck%s_split%d%s%s%s.%s", saveComment.Data(), iSplit+1, compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
            }
          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // MC split loop
  } // Unfolded to truth comparison

  if(drawForwardFoldedToMeasuredComparison){

    // In order not to draw too much curves to one plot, treat the two splits separately
    for(int iSplit = 0; iSplit < nInputFiles; iSplit++){

      for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){

         // Set the centrality information for legends and figure saving
        if(isPbPbData) {
          centralityString = Form("PbPb: %.0f-%.0f", unfoldingCard[iSplit]->GetLowBinBorderCentrality(iCentrality), unfoldingCard[iSplit]->GetHighBinBorderCentrality(iCentrality));
          compactCentralityString = Form("_C%.0f-%.0f", unfoldingCard[iSplit]->GetLowBinBorderCentrality(iCentrality), unfoldingCard[iSplit]->GetHighBinBorderCentrality(iCentrality));
        } else {
          centralityString = "pp";
          compactCentralityString = "_pp";
        }

        for(int iJetPt = firstStudiedJetPtBin; iJetPt <= lastStudiedJetPtBin; iJetPt++){

          // Set the jet pT information for legends and figure saving
          jetPtString = Form("%.0f < jet p_{T} < %.0f", unfoldingCard[iSplit]->GetLowBinBorderJetPtEEC(iJetPt), unfoldingCard[iSplit]->GetHighBinBorderJetPtEEC(iJetPt));
          compactJetPtString = Form("_J=%.0f-%.0f", unfoldingCard[iSplit]->GetLowBinBorderJetPtEEC(iJetPt), unfoldingCard[iSplit]->GetHighBinBorderJetPtEEC(iJetPt));

          for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

            // Set the track pT information for legends and figure saving
            trackPtString = Form("%.1f < track p_{T}", unfoldingCard[iSplit]->GetLowBinBorderTrackPtEEC(iTrackPt));
            compactTrackPtString = Form("_T%.1f", unfoldingCard[iSplit]->GetLowBinBorderTrackPtEEC(iTrackPt));
            compactTrackPtString.ReplaceAll(".","v");
              
            // Draw the ratios to the canvas
            for(int iIteration = 0; iIteration < nIterations; iIteration++){
              hForwardFoldedToMeasuredRatio[iCentrality][iJetPt][iTrackPt][iIteration][iSplit]->SetLineColor(iterationColor[iIteration]);
              if(iIteration == 0){
                hForwardFoldedToMeasuredRatio[iCentrality][iJetPt][iTrackPt][iIteration][iSplit]->GetYaxis()->SetRangeUser(0.8, 1.2);
                hForwardFoldedToMeasuredRatio[iCentrality][iJetPt][iTrackPt][iIteration][iSplit]->GetXaxis()->SetRangeUser(normalizationRegionLow, normalizationRegionHigh);
                drawer->DrawHistogram(hForwardFoldedToMeasuredRatio[iCentrality][iJetPt][iTrackPt][iIteration][iSplit], "#Deltar", "Ratio to measured", " ", "");
              } else {
                hForwardFoldedToMeasuredRatio[iCentrality][iJetPt][iTrackPt][iIteration][iSplit]->Draw("same");
              }
            }

            // Add a legend to the figure
            legend = new TLegend(0.2, 0.75, 0.5, 0.9);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);

            legend->AddEntry((TObject*)0, Form("Split %d   %s", iSplit+1, centralityString.Data()), "");
            legend->AddEntry((TObject*)0, Form("%s   %s", jetPtString.Data(), trackPtString.Data()), "");
              

            // Put the iteration numbers to different legend
            for(int iIterationLegend = 0; iIterationLegend < nIterationLegends; iIterationLegend++){
              iterationLegend[iIterationLegend] = new TLegend(0.17+0.19*iIterationLegend, 0.18, 0.4+0.19*iIterationLegend, 0.38);
              iterationLegend[iIterationLegend]->SetFillStyle(0); 
              iterationLegend[iIterationLegend]->SetBorderSize(0); 
              iterationLegend[iIterationLegend]->SetTextSize(0.05); 
              iterationLegend[iIterationLegend]->SetTextFont(62);
            }

            for(int iIteration = 0; iIteration < nIterations; iIteration++){
              iterationLegend[iIteration/(nIterations/nIterationLegends)]->AddEntry(hForwardFoldedToMeasuredRatio[iCentrality][iJetPt][iTrackPt][iIteration][iSplit], Form("%s", iterationString.ReplaceAll("NIT",Form("%d",iIteration+1)).Data()), "l");
              iterationString.ReplaceAll(Form("%d",iIteration+1),"NIT");
            }

            oneLine->Draw();
            legend->Draw();
            for(int iIterationLegend = 0; iIterationLegend < nIterationLegends; iIterationLegend++){
              iterationLegend[iIterationLegend]->Draw();
            }

            // Save the figures to a file
            if(saveFigures) {
              gPad->GetCanvas()->SaveAs(Form("figures/forwardFoldedToMeasuredRatioCheck%s_split%d%s%s%s.%s", saveComment.Data(), iSplit+1, compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
            }
          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // MC split loop
  } // Unfolded to truth comparison

}
