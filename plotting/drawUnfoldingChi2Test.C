#include "EECCard.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"
#include "AlgorithmLibrary.h"

/*
 * Draw the chi2 test results to determine a good number of iterations for unfolding
 */
void drawUnfoldingChi2Test(){

  // **********************************
  //       Open the input files
  // **********************************

  const int nInputFiles = 2;
  TString inputFileName[] = {"chi2Files/chi2Histograms_PbPb_regular_threeTrackPtBins_2023-06-02.root", "chi2Files/chi2Histograms_PbPb_swapped_threeTrackPtBins_2023-06-02.root"};
  // chi2Histograms_pp_regular_2023-05-26.root
  // chi2Histograms_pp_swapped_2023-05-26.root
  // chi2Histograms_PbPb_regular_threeTrackPtBins_2023-06-02.root
  // chi2Histograms_PbPb_swapped_threeTrackPtBins_2023-06-02.root

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
  
  int firstStudiedTrackPtBinEEC = 4;
  int lastStudiedTrackPtBinEEC = 4;

  int firstStudiedJetPtBin = 0;
  int lastStudiedJetPtBin = nJetPtBins-1;

  const bool drawChi2map = true;                      // Draw the chi2 values for individual jet pT bins
  const bool drawChi2combined = false;                 // Draw single good chi2 value for each response matrix determined from relevent region
  const bool drawUnfoldedToTruthComparison = false;    // Compare unfolded distributions to truth
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

  // Set manually the best number of iteration determined totally objectively
  if(isPbPbData){
    bestNumberOfIterations[0][0][3] = 5;  // Centrality = 0-10, track pT > 2 GeV, 120 < jet pT < 140
    bestNumberOfIterations[0][1][3] = 3;  // Centrality = 0-10, track pT > 2 GeV, 140 < jet pT < 160
    bestNumberOfIterations[0][2][3] = 3;  // Centrality = 0-10, track pT > 2 GeV, 160 < jet pT < 180
    bestNumberOfIterations[0][3][3] = 3;  // Centrality = 0-10, track pT > 2 GeV, 180 < jet pT < 200
    bestNumberOfIterations[1][0][3] = 5;  // Centrality = 10-30, track pT > 2 GeV, 120 < jet pT < 140
    bestNumberOfIterations[1][1][3] = 3;  // Centrality = 10-30, track pT > 2 GeV, 140 < jet pT < 160
    bestNumberOfIterations[1][2][3] = 3;  // Centrality = 10-30, track pT > 2 GeV, 160 < jet pT < 180
    bestNumberOfIterations[1][3][3] = 5;  // Centrality = 10-30, track pT > 2 GeV, 180 < jet pT < 200
    bestNumberOfIterations[2][0][3] = 4;  // Centrality = 30-50, track pT > 2 GeV, 120 < jet pT < 140
    bestNumberOfIterations[2][1][3] = 3;  // Centrality = 30-50, track pT > 2 GeV, 140 < jet pT < 160
    bestNumberOfIterations[2][2][3] = 5;  // Centrality = 30-50, track pT > 2 GeV, 160 < jet pT < 180
    bestNumberOfIterations[2][3][3] = 4;  // Centrality = 30-50, track pT > 2 GeV, 180 < jet pT < 200
    bestNumberOfIterations[3][0][3] = 3;  // Centrality = 50-90, track pT > 2 GeV, 120 < jet pT < 140
    bestNumberOfIterations[3][1][3] = 3;  // Centrality = 50-90, track pT > 2 GeV, 140 < jet pT < 160
    bestNumberOfIterations[3][2][3] = 3;  // Centrality = 50-90, track pT > 2 GeV, 160 < jet pT < 180
    bestNumberOfIterations[3][3][3] = 8;  // Centrality = 50-90, track pT > 2 GeV, 180 < jet pT < 200

    bestNumberOfIterations[0][0][4] = 6;  // Centrality = 0-10, track pT > 2.5 GeV, 120 < jet pT < 140
    bestNumberOfIterations[0][1][4] = 3;  // Centrality = 0-10, track pT > 2.5 GeV, 140 < jet pT < 160
    bestNumberOfIterations[0][2][4] = 3;  // Centrality = 0-10, track pT > 2.5 GeV, 160 < jet pT < 180
    bestNumberOfIterations[0][3][4] = 3;  // Centrality = 0-10, track pT > 2.5 GeV, 180 < jet pT < 200
    bestNumberOfIterations[1][0][4] = 5;  // Centrality = 10-30, track pT > 2.5 GeV, 120 < jet pT < 140
    bestNumberOfIterations[1][1][4] = 3;  // Centrality = 10-30, track pT > 2.5 GeV, 140 < jet pT < 160
    bestNumberOfIterations[1][2][4] = 3;  // Centrality = 10-30, track pT > 2.5 GeV, 160 < jet pT < 180
    bestNumberOfIterations[1][3][4] = 3;  // Centrality = 10-30, track pT > 2.5 GeV, 180 < jet pT < 200
    bestNumberOfIterations[2][0][4] = 4;  // Centrality = 30-50, track pT > 2.5 GeV, 120 < jet pT < 140
    bestNumberOfIterations[2][1][4] = 3;  // Centrality = 30-50, track pT > 2.5 GeV, 140 < jet pT < 160
    bestNumberOfIterations[2][2][4] = 5;  // Centrality = 30-50, track pT > 2.5 GeV, 160 < jet pT < 180
    bestNumberOfIterations[2][3][4] = 4;  // Centrality = 30-50, track pT > 2.5 GeV, 180 < jet pT < 200
    bestNumberOfIterations[3][0][4] = 3;  // Centrality = 50-90, track pT > 2.5 GeV, 120 < jet pT < 140
    bestNumberOfIterations[3][1][4] = 3;  // Centrality = 50-90, track pT > 2.5 GeV, 140 < jet pT < 160
    bestNumberOfIterations[3][2][4] = 3;  // Centrality = 50-90, track pT > 2.5 GeV, 160 < jet pT < 180
    bestNumberOfIterations[3][3][4] = 8;  // Centrality = 50-90, track pT > 2.5 GeV, 180 < jet pT < 200

    bestNumberOfIterations[0][0][5] = 8;  // Centrality = 0-10, track pT > 3 GeV, 120 < jet pT < 140
    bestNumberOfIterations[0][1][5] = 4;  // Centrality = 0-10, track pT > 3 GeV, 140 < jet pT < 160
    bestNumberOfIterations[0][2][5] = 3;  // Centrality = 0-10, track pT > 3 GeV, 160 < jet pT < 180
    bestNumberOfIterations[0][3][5] = 3;  // Centrality = 0-10, track pT > 3 GeV, 180 < jet pT < 200
    bestNumberOfIterations[1][0][5] = 5;  // Centrality = 10-30, track pT > 3 GeV, 120 < jet pT < 140
    bestNumberOfIterations[1][1][5] = 3;  // Centrality = 10-30, track pT > 3 GeV, 140 < jet pT < 160
    bestNumberOfIterations[1][2][5] = 3;  // Centrality = 10-30, track pT > 3 GeV, 160 < jet pT < 180
    bestNumberOfIterations[1][3][5] = 4;  // Centrality = 10-30, track pT > 3 GeV, 180 < jet pT < 200
    bestNumberOfIterations[2][0][5] = 4;  // Centrality = 30-50, track pT > 3 GeV, 120 < jet pT < 140
    bestNumberOfIterations[2][1][5] = 8;  // Centrality = 30-50, track pT > 3 GeV, 140 < jet pT < 160
    bestNumberOfIterations[2][2][5] = 5;  // Centrality = 30-50, track pT > 3 GeV, 160 < jet pT < 180
    bestNumberOfIterations[2][3][5] = 4;  // Centrality = 30-50, track pT > 3 GeV, 180 < jet pT < 200
    bestNumberOfIterations[3][0][5] = 3;  // Centrality = 50-90, track pT > 3 GeV, 120 < jet pT < 140
    bestNumberOfIterations[3][1][5] = 3;  // Centrality = 50-90, track pT > 3 GeV, 140 < jet pT < 160
    bestNumberOfIterations[3][2][5] = 3;  // Centrality = 50-90, track pT > 3 GeV, 160 < jet pT < 180
    bestNumberOfIterations[3][3][5] = 8;  // Centrality = 50-90, track pT > 3 GeV, 180 < jet pT < 200

    // Best number of iteratios determined for combining all jet pT bins
    if(oneIterationPerMatrix) {
      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++) bestNumberOfIterations[0][iJetPt][3] = 4;  // Centrality = 0-10, track pT > 2 GeV
      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++) bestNumberOfIterations[1][iJetPt][3] = 4;  // Centrality = 10-30, track pT > 2 GeV
      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++) bestNumberOfIterations[2][iJetPt][3] = 3;  // Centrality = 30-50, track pT > 2 GeV
      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++) bestNumberOfIterations[3][iJetPt][3] = 3;  // Centrality = 50-90, track pT > 2 GeV

      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++) bestNumberOfIterations[0][iJetPt][4] = 4;  // Centrality = 0-10, track pT > 2.5 GeV
      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++) bestNumberOfIterations[1][iJetPt][4] = 4;  // Centrality = 10-30, track pT > 2.5 GeV
      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++) bestNumberOfIterations[2][iJetPt][4] = 3;  // Centrality = 30-50, track pT > 2.5 GeV
      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++) bestNumberOfIterations[3][iJetPt][4] = 3;  // Centrality = 50-90, track pT > 2.5 GeV

      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++) bestNumberOfIterations[0][iJetPt][5] = 5;  // Centrality = 0-10, track pT > 3 GeV
      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++) bestNumberOfIterations[1][iJetPt][5] = 4;  // Centrality = 10-30, track pT > 3 GeV
      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++) bestNumberOfIterations[2][iJetPt][5] = 4;  // Centrality = 30-50, track pT > 3 GeV
      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++) bestNumberOfIterations[3][iJetPt][5] = 3;  // Centrality = 50-90, track pT > 3 GeV
    }
  } else {
    bestNumberOfIterations[0][0][5] = 7;  // pp, track pT > 3 GeV, 120 < jet pT < 140
    bestNumberOfIterations[0][1][5] = 3;  // pp, track pT > 3 GeV, 140 < jet pT < 160
    bestNumberOfIterations[0][2][5] = 7;  // pp, track pT > 3 GeV, 160 < jet pT < 180
    bestNumberOfIterations[0][3][5] = 2;  // pp, track pT > 3 GeV, 180 < jet pT < 200

    // Best number of iteratios determined for combining all jet pT bins
    if(oneIterationPerMatrix) {
      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++) bestNumberOfIterations[0][iJetPt][5] = 2;  // pp, track pT > 3 GeV
    }
  }

  bool saveFigures = false;
  TString saveComment = "_splitComparisonSingleIteration";
  TString figureFormat = "pdf";

  // Histograms for chi2 and error2 map
  TH1D* hChi2map[nCentralityBins][nJetPtBins][nTrackPtBinsEEC][nInputFiles+1];
  TH1D* hChi2combined[nCentralityBins][nTrackPtBinsEEC][nInputFiles+1];
  TH1D* hErrorMap[nCentralityBins][nJetPtBins][nTrackPtBinsEEC][nInputFiles+1];

  // Initialize all the histograms to null
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
      for(int iFile = 0; iFile < nInputFiles+1; iFile++){
        for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){
          hChi2map[iCentrality][iJetPt][iTrackPt][iFile] = NULL;
          hErrorMap[iCentrality][iJetPt][iTrackPt][iFile] = NULL;
        } // Jet pT loop
        hChi2combined[iCentrality][iTrackPt][iFile] = NULL;
      } // File loop
    } // Track pT loop
  }  // Centrality loop

  // Read the selected histograms from the input files
  for(int iFile = 0; iFile < nInputFiles; iFile++){
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
        for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){
          hChi2map[iCentrality][iJetPt][iTrackPt][iFile] = (TH1D*) inputFile[iFile]->Get(Form("chi2map/chi2map%d%d%d", iCentrality, iTrackPt, iJetPt));
          hErrorMap[iCentrality][iJetPt][iTrackPt][iFile] = (TH1D*) inputFile[iFile]->Get(Form("errorMap/errorMap%d%d%d", iCentrality, iTrackPt, iJetPt));
        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
  } // File loop

  // Add all the files together for a combined result
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){
        hChi2map[iCentrality][iJetPt][iTrackPt][nInputFiles] = (TH1D*) hChi2map[iCentrality][iJetPt][iTrackPt][0]->Clone(Form("averageChi2map%d%d%d", iCentrality, iTrackPt, iJetPt));
        for(int iFile = 1; iFile < nInputFiles; iFile++){
          hChi2map[iCentrality][iJetPt][iTrackPt][nInputFiles]->Add(hChi2map[iCentrality][iJetPt][iTrackPt][iFile]);
        }
        for(int iBin = 1; iBin <= hChi2map[iCentrality][iJetPt][iTrackPt][nInputFiles]->GetNbinsX(); iBin++){
          hChi2map[iCentrality][iJetPt][iTrackPt][nInputFiles]->SetBinContent(iBin, hChi2map[iCentrality][iJetPt][iTrackPt][nInputFiles]->GetBinContent(iBin) / nInputFiles);
        }
      } // Jet pT loop
    } // Track pT loop
  } // Centrality loop

  // Combine the chi2 values from each response matrix to a single grand chi2 value
  for(int iFile = 0; iFile < nInputFiles+1; iFile++){
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
        hChi2combined[iCentrality][iTrackPt][iFile] = (TH1D*) hChi2map[iCentrality][0][iTrackPt][iFile]->Clone(Form("combinedChi2%d%d%d", iCentrality, iTrackPt, iFile));
        for(int iJetPt = 1; iJetPt < nJetPtBins; iJetPt++){
          hChi2combined[iCentrality][iTrackPt][iFile]->Add(hChi2map[iCentrality][iJetPt][iTrackPt][iFile]);
        }
      } // Track pT loop
    } // Centrality loop
  } // File loop

  // Determine the number of iterations from the chi2 histograms
  const int nIterations = hChi2map[firstStudiedCentralityBin][0][firstStudiedTrackPtBinEEC][0]->GetNbinsX();

  // Histograms for energy-energy correlator distributions and their ratios
  TH1D* hMeasured[nCentralityBins][nJetPtBins][nTrackPtBinsEEC][nInputFiles];
  TH1D* hTruth[nCentralityBins][nJetPtBins][nTrackPtBinsEEC][nInputFiles];
  TH1D* hUnfolded[nCentralityBins][nJetPtBins][nTrackPtBinsEEC][nIterations][nInputFiles];
  TH1D* hMeasuredToTruthRatio[nCentralityBins][nJetPtBins][nTrackPtBinsEEC][nInputFiles];
  TH1D* hUnfoldedToTruthRatio[nCentralityBins][nJetPtBins][nTrackPtBinsEEC][nIterations][nInputFiles];

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
          }
        } // File loop
      } // Jet pT loop
    } // Track pT loop
  }  // Centrality loop

  // Read the energy-energy correlator histograms and ratios from the input file
  for(int iFile = 0; iFile < nInputFiles; iFile++) {
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++) {
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++) {
        for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++) {
          hMeasured[iCentrality][iJetPt][iTrackPt][iFile] = (TH1D*) inputFile[iFile]->Get(Form("measured/hMeasured%d%d%d", iCentrality, iTrackPt, iJetPt));
          hTruth[iCentrality][iJetPt][iTrackPt][iFile] = (TH1D*) inputFile[iFile]->Get(Form("truth/hTruth%d%d%d", iCentrality, iTrackPt, iJetPt));
          hMeasuredToTruthRatio[iCentrality][iJetPt][iTrackPt][iFile] = (TH1D*) inputFile[iFile]->Get(Form("measuredToTruthRatio/measuredToTruthRatio%d%d%d", iCentrality, iTrackPt, iJetPt));
          for(int iIteration = 0; iIteration < nIterations; iIteration++) {
            hUnfolded[iCentrality][iJetPt][iTrackPt][iIteration][iFile] = (TH1D*) inputFile[iFile]->Get(Form("unfolded/hUnfolded%d%d%d%d", iCentrality, iTrackPt, iJetPt, iIteration));
            hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][iIteration][iFile] = (TH1D*) inputFile[iFile]->Get(Form("unfoldedToTruthRatio/unfoldedToTruthRatio%d%d%d%d", iCentrality, iTrackPt, iJetPt, iIteration));
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
  double normalizationRegionLow = 0.006;
  double normalizationRegionHigh = 0.39;

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

      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){

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
        legend = new TLegend(0.18, 0.54, 0.43, 0.87);
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

      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){

        // Set the jet pT information for legends and figure saving
        jetPtString = Form("%.0f < jet p_{T} < %.0f", unfoldingCard[0]->GetLowBinBorderJetPtEEC(iJetPt), unfoldingCard[0]->GetHighBinBorderJetPtEEC(iJetPt));
        compactJetPtString = Form("_J=%.0f-%.0f", unfoldingCard[0]->GetLowBinBorderJetPtEEC(iJetPt), unfoldingCard[0]->GetHighBinBorderJetPtEEC(iJetPt));

        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

          // Set the track pT information for legends and figure saving
          trackPtString = Form("%.1f < track p_{T}", unfoldingCard[0]->GetLowBinBorderTrackPtEEC(iTrackPt));
          compactTrackPtString = Form("_T%.0f", unfoldingCard[0]->GetLowBinBorderTrackPtEEC(iTrackPt));

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
  } // Unfolded to truth comparison

  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);

  if(drawUnfoldedToTruthComparison){

    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){

      // Set the centrality information for legends and figure saving
      if(isPbPbData) {
        centralityString = Form("Pythia+Hydjet: %.0f-%.0f", unfoldingCard[0]->GetLowBinBorderCentrality(iCentrality), unfoldingCard[0]->GetHighBinBorderCentrality(iCentrality));
        compactCentralityString = Form("_C%.0f-%.0f", unfoldingCard[0]->GetLowBinBorderCentrality(iCentrality), unfoldingCard[0]->GetHighBinBorderCentrality(iCentrality));
      } else {
        centralityString = "Pythia8";
        compactCentralityString = "_pythia8";
      }

      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){

        // Set the jet pT information for legends and figure saving
        jetPtString = Form("%.0f < jet p_{T} < %.0f", unfoldingCard[0]->GetLowBinBorderJetPtEEC(iJetPt), unfoldingCard[0]->GetHighBinBorderJetPtEEC(iJetPt));
        compactJetPtString = Form("_J=%.0f-%.0f", unfoldingCard[0]->GetLowBinBorderJetPtEEC(iJetPt), unfoldingCard[0]->GetHighBinBorderJetPtEEC(iJetPt));

        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

          // Create a new canvas for the plot
          drawer->CreateSplitCanvas();

          // Set the track pT information for legends and figure saving
          trackPtString = Form("%.1f < track p_{T}", unfoldingCard[0]->GetLowBinBorderTrackPtEEC(iTrackPt));
          compactTrackPtString = Form("_T%.0f", unfoldingCard[0]->GetLowBinBorderTrackPtEEC(iTrackPt));

          // Draw first the generator level distribution
          hTruth[iCentrality][iJetPt][iTrackPt][0]->SetLineColor(kBlack);
          hTruth[iCentrality][iJetPt][iTrackPt][0]->GetXaxis()->SetRangeUser(normalizationRegionLow, normalizationRegionHigh);
          drawer->SetLogY(true);
          drawer->DrawHistogramToUpperPad(hTruth[iCentrality][iJetPt][iTrackPt][0], "#Deltar", "EEC", " ", "");

          // Add the reconstructed and unfolded distributions to the same plot
          hMeasured[iCentrality][iJetPt][iTrackPt][0]->SetLineColor(kRed);
          hMeasured[iCentrality][iJetPt][iTrackPt][0]->Draw("same");

          //for(int iIteration = 0; iIteration < nIterations; iIteration++){
          //  hUnfolded[iCentrality][iJetPt][iTrackPt][iIteration]->SetLineColor(iterationColor[iIteration]);
          //  hUnfolded[iCentrality][iJetPt][iTrackPt][iIteration]->Draw("same");
          //} 

          // Add a legend to the figure
          legend = new TLegend(0.25, 0.18, 0.5, 0.48);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);

          legend->AddEntry((TObject*)0, centralityString.Data(), "");
          legend->AddEntry((TObject*)0, jetPtString.Data(), "");
          legend->AddEntry((TObject*)0, trackPtString.Data(), "");

          legend->AddEntry(hTruth[iCentrality][iJetPt][iTrackPt][0], "Generator level reference", "l");
          legend->AddEntry(hMeasured[iCentrality][iJetPt][iTrackPt][0], "Reconstructed jets + gen particles", "l");
          //for(int iIteration = 0; iIteration < nIterations; iIteration++){
          //  legend->AddEntry(hUnfolded[iCentrality][iJetPt][iTrackPt][iIteration], Form("Unfolded result (%s)", unfoldName.ReplaceAll("NIT",Form("%d",iIteration+1)).Data()), "l");
          //  unfoldName.ReplaceAll(Form("%d",iIteration+1),"NIT");
          //}

          legend->Draw();

          // Draw the ratios to lower pad
          drawer->SetLogY(false);
          hMeasuredToTruthRatio[iCentrality][iJetPt][iTrackPt][0]->SetLineColor(kRed);
          hMeasuredToTruthRatio[iCentrality][iJetPt][iTrackPt][0]->GetYaxis()->SetRangeUser(0.8, 1.2);
          hMeasuredToTruthRatio[iCentrality][iJetPt][iTrackPt][0]->GetXaxis()->SetRangeUser(normalizationRegionLow, normalizationRegionHigh);
          drawer->DrawHistogramToLowerPad(hMeasuredToTruthRatio[iCentrality][iJetPt][iTrackPt][0], "#Deltar", "Ratio to truth", " ", "");

          //for(int iIteration = 0; iIteration < nIterations; iIteration++){
          //  hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][iIteration]->SetLineColor(iterationColor[iIteration]);
          //  hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][iIteration]->Draw("same");
          //}

          oneLine->Draw();

          // Save the figures to a file
          if(saveFigures) {
            gPad->GetCanvas()->SaveAs(Form("figures/jetPtUnfoldingWithRooUnfoldTruthCheck%s%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
          }
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Unfolded to truth comparison

}
