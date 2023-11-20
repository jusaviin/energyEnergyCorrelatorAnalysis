#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Reser underflow and overflow bins for a one-dimensional histogram
 */
void removeOutOfRange(TH1D* histogramInNeedOfTrimming){
  histogramInNeedOfTrimming->SetBinContent(0, 0);
  histogramInNeedOfTrimming->SetBinContent(histogramInNeedOfTrimming->GetNbinsX()+1, 0);
}

/*
 * Reser underflow and overflow bins for a two-dimensional histogram
 */
void removeOutOfRange(TH2D* histogramInNeedOfTrimming)
{
  int NX = histogramInNeedOfTrimming->GetNbinsX();
  int NY = histogramInNeedOfTrimming->GetNbinsY();

  for(int iX = 0; iX <= NX + 1; iX++){
    histogramInNeedOfTrimming->SetBinContent(iX, 0, 0);
    histogramInNeedOfTrimming->SetBinContent(iX, NY+1, 0);
  }
  for(int iY = 0; iY <= NY + 1; iY++){
    histogramInNeedOfTrimming->SetBinContent(0, iY, 0);
    histogramInNeedOfTrimming->SetBinContent(NX+1, iY, 0);
  }
}

/*
 * Forward fold the generator level distribution using the response matrix to get the detector level distribution
 *
 * Arguments:
 *  TH1D* unfoldedHistogram = Histogram with unfolded results. This will be forward folded
 *  TH2D* responseMatrix = Response matrix from which the forward folding is calculated
 *
 *  return: Forward folded histogram
 */
TH1D* forwardFold(TH1D* unfoldedHistogram, TH2D* responseMatrix){

  // Cannot perform forward folding is the histograms do not exist
  if(unfoldedHistogram == NULL || responseMatrix == NULL) return NULL;

  // Find the number of bins in generator and reconstructed levels
  int nBinsGen = responseMatrix->GetNbinsY();
  int nBinsReco = responseMatrix->GetNbinsX();

  // Create a histogram to hold the forward folded distribution
  TH1D* foldedHistogram = (TH1D*)unfoldedHistogram->Clone(Form("%sFolded", unfoldedHistogram->GetName()));
  foldedHistogram->Reset();

  // Do the actual forward folding
  double recoSum;
  double recoFractionOfTotal, recoFractionOfTotalError;
  double truthContent, truthError;
  double valueUpdate, errorUpdate;

  // Loop over all generator level bins
  for(int iGenBin = 1; iGenBin <= nBinsGen; iGenBin++){

    // Sum all the reconstructed values corresponding to given generator level value
    recoSum = 0;
    for(int iRecoBin = 1; iRecoBin <= nBinsReco; iRecoBin++){
      recoSum = recoSum + responseMatrix->GetBinContent(iRecoBin, iGenBin);
    }

    // If there are no reconstructed entries, also the entry in the folded distribution will be 0
    if(recoSum == 0) continue;

    // Loop over all the reconstructed bins corresponding to the generator level bin again and update their values
    for(int iRecoBin = 1; iRecoBin <= nBinsReco; iRecoBin++){

      // Fraction of the total reconstructed value that is found from this generator level bin
      recoFractionOfTotal = responseMatrix->GetBinContent(iRecoBin, iGenBin) / recoSum;

      // Content of the generator level bin
      truthContent = unfoldedHistogram->GetBinContent(iGenBin);

      // Relative uncertianties of the two values above
      recoFractionOfTotalError = responseMatrix->GetBinError(iRecoBin, iGenBin) / responseMatrix->GetBinContent(iRecoBin, iGenBin);
      truthError = unfoldedHistogram->GetBinError(iGenBin) / unfoldedHistogram->GetBinContent(iGenBin);

      // If the bin content in matrix is 0, ensure that also error for this is 0
      if(responseMatrix->GetBinContent(iRecoBin, iGenBin) == 0) recoFractionOfTotalError = 0;

      // Take the fraction of the truth as an input for the folded histogram
      valueUpdate = recoFractionOfTotal * truthContent;
      errorUpdate = TMath::Sqrt(recoFractionOfTotalError * recoFractionOfTotalError + truthError * truthError) * valueUpdate;

      foldedHistogram->SetBinContent(iRecoBin, foldedHistogram->GetBinContent(iRecoBin) + valueUpdate);
      foldedHistogram->SetBinError(iRecoBin, TMath::Sqrt(foldedHistogram->GetBinError(iRecoBin)*foldedHistogram->GetBinError(iRecoBin) + errorUpdate*errorUpdate));
    } // Reconstructed level bin loop
  } // Generator level bin loop

  return foldedHistogram;
}

/*
 * Do the one dimensional jet pT unfolding using RooUnfold
 */
void unfoldJetPtWithRooUnfoldOneDimensional(){

  // Enumeration variables
  enum enumEnergyEnergyCorrelatorFileType{kUnfoldingFile, kReferenceFile, kNFileTypes};
  enum enumUnfoldingMethods{kMatrixInversionUnfold, kBayesianUnfold, kNUnfoldMethods};

  TString unfoldName[kNUnfoldMethods] = {"martix inversion", "Bayes, NIT iterations"};

  // **********************************
  //       Open the input files
  // **********************************

  // Define the name for the file containing histograms needed for unfolding and reference distributions for checking the quality of unfolding
  TString inputFileName[kNFileTypes];
  inputFileName[kUnfoldingFile] = "data/PbPbMC2018_GenGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_nominalSmear_1Dunfolding_part1_processed_2023-09-14.root";
  // PbPbMC2018_GenGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_nominalSmear_1Dunfolding_part1_processed_2023-09-14.root
  // PbPbMC2018_GenGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_nominalSmear_1Dunfolding_part1_noBinWidthNormalization_processed_2023-09-14.root
  // PbPbMC2018_GenGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_nominalSmear_1DunfoldMoreBins_part1_processed_2023-09-19.root
  inputFileName[kReferenceFile] = "data/PbPbMC2018_GenGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_nominalSmear_1Dunfolding_part2_processed_2023-09-14.root";
  // PbPbMC2018_GenGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_nominalSmear_1Dunfolding_part2_processed_2023-09-14.root
  // PbPbMC2018_GenGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_nominalSmear_1Dunfolding_part2_noBinWidthNormalization_processed_2023-09-14.root
  // PbPbMC2018_GenGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_nominalSmear_1DunfoldMoreBins_part2_processed_2023-09-19.root

  // Open the input files
  TFile* inputFile[kNFileTypes];
  EECCard* unfoldingCard[kNFileTypes];

  for(int iFileType = 0; iFileType < kNFileTypes; iFileType++){
  
    inputFile[iFileType] = TFile::Open(inputFileName[iFileType]);

    if(inputFile[iFileType] == NULL) {
      cout << "Error! The file " << inputFileName[iFileType].Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }

    // For the binning information, read the card
    unfoldingCard[iFileType] = new EECCard(inputFile[iFileType]);

  }

  // Determine if we are dealing with pp or PbPb data
  TString collisionSystem = unfoldingCard[kUnfoldingFile]->GetDataType();
  bool isPbPbData = collisionSystem.Contains("PbPb");
  
  // **********************************************************************************
  //       Check that the centrality, jet pT and track pT bins match with the files
  // **********************************************************************************

  for(int iFileType = kReferenceFile; iFileType < kNFileTypes; iFileType++) {

    if(isPbPbData != (unfoldingCard[iFileType]->GetDataType().Contains("PbPb"))) {
      cout << "You are trying to unfold pp with PbPb unfolding histograms or vice versa!" << endl;
      cout << "Please ensure that the data and unfolding histograms are comparible!" << endl;
    }

    // Only check compatible centrality bins for PbPb
    if(isPbPbData) {
      if(unfoldingCard[kUnfoldingFile]->GetNCentralityBins() != unfoldingCard[iFileType]->GetNCentralityBins()) {
        cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
        return;
      }
      for(int iCentrality = 0; iCentrality < unfoldingCard[kUnfoldingFile]->GetNCentralityBins(); iCentrality++) {
        if(TMath::Abs(unfoldingCard[kUnfoldingFile]->GetLowBinBorderCentrality(iCentrality) - unfoldingCard[iFileType]->GetLowBinBorderCentrality(iCentrality)) > 5) {
          cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
          return;
        }
        if(TMath::Abs(unfoldingCard[kUnfoldingFile]->GetHighBinBorderCentrality(iCentrality) - unfoldingCard[iFileType]->GetHighBinBorderCentrality(iCentrality)) > 5) {
          cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
          return;
        }
      }
    }  // Checking compatible centrality bins

    if(unfoldingCard[kUnfoldingFile]->GetNJetPtBinsEEC() != unfoldingCard[iFileType]->GetNJetPtBinsEEC()) {
      cout << "Error! Jet pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
    for(int iTrackPt = 0; iTrackPt < unfoldingCard[kUnfoldingFile]->GetNJetPtBinsEEC(); iTrackPt++) {
      if(TMath::Abs(unfoldingCard[kUnfoldingFile]->GetLowBinBorderJetPtEEC(iTrackPt) - unfoldingCard[iFileType]->GetLowBinBorderJetPtEEC(iTrackPt)) > 0.01) {
        cout << "Error! Jet pT bins do not match between the two files! Cannot execute the code!" << endl;
        return;
      }
      if(TMath::Abs(unfoldingCard[kUnfoldingFile]->GetHighBinBorderJetPtEEC(iTrackPt) - unfoldingCard[iFileType]->GetHighBinBorderJetPtEEC(iTrackPt)) > 0.01) {
        cout << "Error! Jet pT bins do not match between the two files! Cannot execute the code!" << endl;
        return;
      }
    }  // Checking the compatibility of jet pT bins
  }

  // ********************************************************
  //       Binning configuration for the unfolding study
  // ********************************************************
  
  // Find the number of bins from the card
  const int nCentralityBins = (isPbPbData) ? unfoldingCard[kUnfoldingFile]->GetNCentralityBins() : 1;
  const int nJetPtBinsEEC = unfoldingCard[kUnfoldingFile]->GetNJetPtBinsEEC();

  // Bin range to be studied
  int firstStudiedCentralityBin = 0;
  int lastStudiedCentralityBin = 3;

  // Default binning ranges for reference
  // centrality = {0,10,30,50,90}

  const int iUnfoldMethod = kBayesianUnfold; // Select the unfolding method: kMatrixInversionUnfold or kBayesianUnfold
  const int nIterations = (iUnfoldMethod == kMatrixInversionUnfold) ? 1 : 5; // Number of iterations used with the regularized unfolding
  const int iterationKey[] = {1,2,3,4,5};
  const bool calculateConditionNumber = false;

  // Flag to unfold the same dataset from which the response matrix in constructed
  const bool unfoldResponseMatrixDataset = true;

  const bool drawUnfoldedToTruthComparison = true;    // Compare unfolded distribution to truth reference
  const bool drawResponseMatrix = true;               // Draw the used response matrices
  const bool drawRefoldingTest = false;                // Compare refolded distribution to the original measured distribution

  std::pair<double,double> ratioZoom = std::make_pair(0,2);

  bool saveFigures = false;
  TString saveComment = "_bayesUnfold";
  TString figureFormat = "pdf";
    
  // ***************************************************************
  //    Create histogram managers and load the needed histograms
  // ***************************************************************

  // Load the histograms needed for unfolding from the unfolding histogram manager
  EECHistogramManager* histogramManager[kNFileTypes];
  for(int iFileType = 0; iFileType < kNFileTypes; iFileType++){
    histogramManager[iFileType] = new EECHistogramManager(inputFile[iFileType], unfoldingCard[iFileType]);
    histogramManager[iFileType]->SetLoadJetPtOneDimensionalUnfoldingHistograms(true);
    histogramManager[iFileType]->SetCentralityBinRange(firstStudiedCentralityBin,lastStudiedCentralityBin);
    histogramManager[iFileType]->LoadProcessedHistograms();
  }

  // Histograms that are needed to create the unfolding response
  TH1D* hUnfoldingMeasured[nCentralityBins];
  TH1D* hUnfoldingTruth[nCentralityBins];
  TH2D* hUnfoldingResponse[nCentralityBins];
  TH1D* hUnfoldedDistribution[nCentralityBins][nIterations];
  TH1D* hForwardFoldedDistribution[nCentralityBins][nIterations];
  TH1D* hReferenceDistributionMeasured[nCentralityBins];
  TH1D* hReferenceDistributionTruth[nCentralityBins];

  // RooUnfold response object
  RooUnfoldResponse* rooResponse[nCentralityBins];

  // Initialize all the histograms to null
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      hUnfoldingMeasured[iCentrality] = NULL;
      hUnfoldingTruth[iCentrality] = NULL;
      hUnfoldingResponse[iCentrality] = NULL;
      hReferenceDistributionMeasured[iCentrality] = NULL;
      hReferenceDistributionTruth[iCentrality] = NULL;
      rooResponse[iCentrality] = NULL;
      for(int iIteration = 0; iIteration < nIterations; iIteration++){
        hUnfoldedDistribution[iCentrality][iIteration] = NULL;
        hForwardFoldedDistribution[iCentrality][iIteration] = NULL;
      }
  }  // Centrality loop

  // Read the histograms needed for the unfolding response and create the RooUnfold response objects
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      hUnfoldingMeasured[iCentrality] = histogramManager[kUnfoldingFile]->GetHistogramJetPtOneDimensionalUnfoldingMeasured(iCentrality);
      hUnfoldingTruth[iCentrality] = histogramManager[kUnfoldingFile]->GetHistogramJetPtOneDimensionalUnfoldingTruth(iCentrality);
      hUnfoldingResponse[iCentrality] = histogramManager[kUnfoldingFile]->GetHistogramJetPtOneDimensionalUnfoldingResponse(iCentrality);

      // Remove the overflow and underflow bins from the histograms
      removeOutOfRange(hUnfoldingMeasured[iCentrality]);
      removeOutOfRange(hUnfoldingTruth[iCentrality]);
      removeOutOfRange(hUnfoldingResponse[iCentrality]);

      // Create the response object with trimmed histograms
      rooResponse[iCentrality] = new RooUnfoldResponse(hUnfoldingMeasured[iCentrality], hUnfoldingTruth[iCentrality], hUnfoldingResponse[iCentrality]);

      // Read also the reference histograms
      hReferenceDistributionMeasured[iCentrality] = histogramManager[kReferenceFile]->GetHistogramJetPtOneDimensionalUnfoldingMeasured(iCentrality);
      hReferenceDistributionTruth[iCentrality] = histogramManager[kReferenceFile]->GetHistogramJetPtOneDimensionalUnfoldingTruth(iCentrality);

      // Remove the overflow and underflow bins from the histograms
      removeOutOfRange(hReferenceDistributionMeasured[iCentrality]);
      removeOutOfRange(hReferenceDistributionTruth[iCentrality]);

  } // Centrality loop

  // ******************************************************************
  //       Calculating condition number for the response matrix
  // ******************************************************************

  if(calculateConditionNumber){
    TDecompSVD* svd[nCentralityBins];  // TDecompSVD is the singular value decomposition (SVD) matrix.
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      svd[iCentrality] = new TDecompSVD(rooResponse[iCentrality]->Mresponse()); // The response->Mresponse() returns the normalized migration matrix
      auto singular_values = svd[iCentrality]->GetSig(); // this is a vector with the singular values, i.e., the diagonal elements of S. They are ordered from largest to smallest.
      cout << "Condition number for centrality bin " << iCentrality << " is: " << endl;
      cout << singular_values.Max() / singular_values.Min() << endl;
    } // Centrality loop
  }

  // ********************************************************
  //                     The actual unfolding
  // ********************************************************

  if(iUnfoldMethod == kMatrixInversionUnfold){

    // Unfolding using matrix inversion
    RooUnfoldInvert* invertUnfold[nCentralityBins];
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      invertUnfold[iCentrality] = new RooUnfoldInvert(rooResponse[iCentrality], unfoldResponseMatrixDataset ? hUnfoldingMeasured[iCentrality] : hReferenceDistributionMeasured[iCentrality]);
      hUnfoldedDistribution[iCentrality][0] = (TH1D*)invertUnfold[iCentrality]->Hunfold();

      // Forward fold the unfolded distribution
      hForwardFoldedDistribution[iCentrality][0] = forwardFold(hUnfoldedDistribution[iCentrality][0], hUnfoldingResponse[iCentrality]);

    } // Centrality loop

  } else if(iUnfoldMethod == kBayesianUnfold){

    // Unfolding using Bayesian unfolding
    RooUnfoldBayes* bayesUnfold[nCentralityBins][nIterations];
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iIteration = 0; iIteration < nIterations; iIteration++){
        bayesUnfold[iCentrality][iIteration] = new RooUnfoldBayes(rooResponse[iCentrality], unfoldResponseMatrixDataset ? hUnfoldingMeasured[iCentrality] : hReferenceDistributionMeasured[iCentrality], iterationKey[iIteration]);
        hUnfoldedDistribution[iCentrality][iIteration] = (TH1D*)bayesUnfold[iCentrality][iIteration]->Hunfold();

        // Forward fold the unfolded distribution
        hForwardFoldedDistribution[iCentrality][iIteration] = forwardFold(hUnfoldedDistribution[iCentrality][iIteration], hUnfoldingResponse[iCentrality]);
        
      } // Iteration loop
    } // Centrality loop
  }


  // **************************************************************************************
  //     Do bin width and total normalization for the distributions and calculate ratios
  // **************************************************************************************

  // Define needed ratio histograms
  TH1D* hMeasuredToTruthRatio[nCentralityBins];
  TH1D* hUnfoldedToTruthRatio[nCentralityBins][nIterations];
  TH1D* hForwardFoldedToMeasuredRatio[nCentralityBins][nIterations];

  // Initialize the dissected histograms to null
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    hMeasuredToTruthRatio[iCentrality] = NULL;
    for(int iIteration = 0; iIteration < nIterations; iIteration++){
      hUnfoldedToTruthRatio[iCentrality][iIteration] = NULL;
      hForwardFoldedToMeasuredRatio[iCentrality][iIteration] = NULL;
    } // Iteration loop
  } // Centrality loop

  // Also normalize the bins to bin width and histograms to one in the defined range
  double normalizationRegionLow = 80+0.1;
  double normalizationRegionHigh = 500-0.1;
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){

    // First do the bin width normalization for all distributions
    /*hUnfoldingMeasured[iCentrality]->Scale(1.0,"width");
    hUnfoldingTruth[iCentrality]->Scale(1.0,"width");
    hReferenceDistributionMeasured[iCentrality]->Scale(1.0,"width");
    hReferenceDistributionTruth[iCentrality]->Scale(1.0,"width");
    for(int iIteration = 0; iIteration < nIterations; iIteration++){
      hUnfoldedDistribution[iCentrality][iIteration]->Scale(1.0,"width");
      hForwardFoldedDistribution[iCentrality][iIteration]->Scale(1.0,"width");
    }*/

    // Next, normalize all the distributions to one in the range 80-500
    hUnfoldingMeasured[iCentrality]->Scale(1.0 / hUnfoldingMeasured[iCentrality]->Integral(hUnfoldingMeasured[iCentrality]->FindBin(normalizationRegionLow), hUnfoldingMeasured[iCentrality]->FindBin(normalizationRegionHigh), "width"));
    hUnfoldingTruth[iCentrality]->Scale(1.0 / hUnfoldingTruth[iCentrality]->Integral(hUnfoldingTruth[iCentrality]->FindBin(normalizationRegionLow), hUnfoldingTruth[iCentrality]->FindBin(normalizationRegionHigh), "width"));
    hReferenceDistributionMeasured[iCentrality]->Scale(1.0 / hReferenceDistributionMeasured[iCentrality]->Integral(hReferenceDistributionMeasured[iCentrality]->FindBin(normalizationRegionLow), hReferenceDistributionMeasured[iCentrality]->FindBin(normalizationRegionHigh), "width"));
    hReferenceDistributionTruth[iCentrality]->Scale(1.0 / hReferenceDistributionTruth[iCentrality]->Integral(hReferenceDistributionTruth[iCentrality]->FindBin(normalizationRegionLow), hReferenceDistributionTruth[iCentrality]->FindBin(normalizationRegionHigh), "width"));

    for(int iIteration = 0; iIteration < nIterations; iIteration++){
      hUnfoldedDistribution[iCentrality][iIteration]->Scale(1.0 / hUnfoldedDistribution[iCentrality][iIteration]->Integral(hUnfoldedDistribution[iCentrality][iIteration]->FindBin(normalizationRegionLow), hUnfoldedDistribution[iCentrality][iIteration]->FindBin(normalizationRegionHigh), "width"));
      hForwardFoldedDistribution[iCentrality][iIteration]->Scale(1.0 / hForwardFoldedDistribution[iCentrality][iIteration]->Integral(hForwardFoldedDistribution[iCentrality][iIteration]->FindBin(normalizationRegionLow), hForwardFoldedDistribution[iCentrality][iIteration]->FindBin(normalizationRegionHigh), "width"));
    }
  } // Centrality loop

  // Calculate the ratios of the other histograms to truth
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    hMeasuredToTruthRatio[iCentrality] = unfoldResponseMatrixDataset ? (TH1D*) hUnfoldingMeasured[iCentrality]->Clone(Form("measuredToTruthRatio%d", iCentrality)) : (TH1D*) hReferenceDistributionMeasured[iCentrality]->Clone(Form("measuredToTruthRatio%d", iCentrality));
    hMeasuredToTruthRatio[iCentrality]->Divide(unfoldResponseMatrixDataset ? hUnfoldingTruth[iCentrality] : hReferenceDistributionTruth[iCentrality]);

    for(int iIteration = 0; iIteration < nIterations; iIteration++){
      hUnfoldedToTruthRatio[iCentrality][iIteration] = (TH1D*) hUnfoldedDistribution[iCentrality][iIteration]->Clone(Form("unfoldedToTruthRatio%d%d", iCentrality, iIteration));
      hUnfoldedToTruthRatio[iCentrality][iIteration]->Divide(unfoldResponseMatrixDataset ? hUnfoldingTruth[iCentrality] : hReferenceDistributionTruth[iCentrality]);

      hForwardFoldedToMeasuredRatio[iCentrality][iIteration] = (TH1D*) hForwardFoldedDistribution[iCentrality][iIteration]->Clone(Form("forwardFoldedToMeasuredRatio%d%d", iCentrality, iIteration));
      hForwardFoldedToMeasuredRatio[iCentrality][iIteration]->Divide(unfoldResponseMatrixDataset ? hUnfoldingMeasured[iCentrality] : hReferenceDistributionMeasured[iCentrality]);
    } // Iteration loop
  } // Centrality loop
  
  // ******************************************************************
  //         Draw different comparisons to check the perfoemance
  // ******************************************************************

  // Prepare a JDrawer for drawing purposes
  JDrawer* drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  
  TLine* oneLine = new TLine(normalizationRegionLow, 1, normalizationRegionHigh, 1);
  oneLine->SetLineStyle(2);

  // Helper variables
  TLegend* legend;
  TString centralityString;
  TString compactCentralityString;

  int iterationColor[] = {kBlue, kGreen+3, kMagenta, kCyan, kOrange, kViolet+3, kPink-7, kSpring+3, kAzure-7};

  if(drawUnfoldedToTruthComparison){

    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){

      // Set the centrality information for legends and figure saving
      if(isPbPbData) {
        centralityString = Form("Pythia+Hydjet: %.0f-%.0f", unfoldingCard[kUnfoldingFile]->GetLowBinBorderCentrality(iCentrality), unfoldingCard[kUnfoldingFile]->GetHighBinBorderCentrality(iCentrality));
        compactCentralityString = Form("_C%.0f-%.0f", unfoldingCard[kUnfoldingFile]->GetLowBinBorderCentrality(iCentrality), unfoldingCard[kUnfoldingFile]->GetHighBinBorderCentrality(iCentrality));
      } else {
        centralityString = "Pythia8";
        compactCentralityString = "_pythia8";
      }

      // Create a new canvas for the plot
      drawer->CreateSplitCanvas();

      // Draw first the generator level distribution
      hUnfoldingTruth[iCentrality]->SetLineColor(kBlack);
      hUnfoldingTruth[iCentrality]->GetXaxis()->SetRangeUser(normalizationRegionLow, normalizationRegionHigh);
      hReferenceDistributionTruth[iCentrality]->SetLineColor(kBlack);
      hReferenceDistributionTruth[iCentrality]->GetXaxis()->SetRangeUser(normalizationRegionLow, normalizationRegionHigh);
      drawer->SetLogY(true);
      drawer->DrawHistogramToUpperPad(unfoldResponseMatrixDataset ? hUnfoldingTruth[iCentrality] : hReferenceDistributionTruth[iCentrality], "Jet p_{T}", "dN/dp_{T}", " ", "");

      // Add the reconstructed and unfolded distributions to the same plot
      hUnfoldingMeasured[iCentrality]->SetLineColor(kRed);
      hReferenceDistributionMeasured[iCentrality]->SetLineColor(kRed);
      unfoldResponseMatrixDataset ? hUnfoldingMeasured[iCentrality]->Draw("same") : hReferenceDistributionMeasured[iCentrality]->Draw("same");

      for(int iIteration = 0; iIteration < nIterations; iIteration++){
        hUnfoldedDistribution[iCentrality][iIteration]->SetLineColor(iterationColor[iIteration]);
        hUnfoldedDistribution[iCentrality][iIteration]->Draw("same");
      } 

      // Add a legend to the figure
      legend = new TLegend(0.25, 0.18- 0.03*nIterations, 0.5, 0.48 + 0.02*nIterations);
      legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);

      legend->AddEntry((TObject*)0, centralityString.Data(), "");
      legend->AddEntry(unfoldResponseMatrixDataset ? hUnfoldingTruth[iCentrality] : hReferenceDistributionTruth[iCentrality], "Generator level jets", "l");
      legend->AddEntry(unfoldResponseMatrixDataset ? hUnfoldingMeasured[iCentrality] : hReferenceDistributionMeasured[iCentrality], "Reconstructed jets", "l");
      for(int iIteration = 0; iIteration < nIterations; iIteration++){
        legend->AddEntry(hUnfoldedDistribution[iCentrality][iIteration], Form("Unfolded result (%s)", unfoldName[iUnfoldMethod].ReplaceAll("NIT",Form("%d",iterationKey[iIteration])).Data()), "l");
        unfoldName[iUnfoldMethod].ReplaceAll(Form("%d",iterationKey[iIteration]),"NIT");
      }

      legend->Draw();

      // Draw the ratios to lower pad
      drawer->SetLogY(false);
      hMeasuredToTruthRatio[iCentrality]->SetLineColor(kRed);
      hMeasuredToTruthRatio[iCentrality]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
      hMeasuredToTruthRatio[iCentrality]->GetXaxis()->SetRangeUser(normalizationRegionLow, normalizationRegionHigh);
      drawer->DrawHistogramToLowerPad(hMeasuredToTruthRatio[iCentrality], "Jet p_{T}", "Ratio to truth", " ", "");

      for(int iIteration = 0; iIteration < nIterations; iIteration++){
        hUnfoldedToTruthRatio[iCentrality][iIteration]->SetLineColor(iterationColor[iIteration]);
        hUnfoldedToTruthRatio[iCentrality][iIteration]->Draw("same");
      }

      oneLine->Draw();

      // Save the figures to a file
      if(saveFigures) {
        gPad->GetCanvas()->SaveAs(Form("figures/jetPtOneDimensionalUnfoldingTruthCheck%s%s.%s", saveComment.Data(), compactCentralityString.Data(), figureFormat.Data()));
      }
    } // Centrality loop
  } // Unfolded to truth comparison

  if(drawResponseMatrix){

    drawer->Reset();
    drawer->SetRelativeCanvasSize(1.1,1.2);
    drawer->SetLogZ(true);
    drawer->SetLeftMargin(0.12);
    drawer->SetTopMargin(0.07);
    drawer->SetRightMargin(0.12);
    drawer->SetTitleOffsetY(0.8);

    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){

      // Set the centrality information for legends and figure saving
      if(isPbPbData) {
        centralityString = Form("Pythia+Hydjet: %.0f-%.0f", unfoldingCard[kUnfoldingFile]->GetLowBinBorderCentrality(iCentrality), unfoldingCard[kUnfoldingFile]->GetHighBinBorderCentrality(iCentrality));
        compactCentralityString = Form("_C%.0f-%.0f", unfoldingCard[kUnfoldingFile]->GetLowBinBorderCentrality(iCentrality), unfoldingCard[kUnfoldingFile]->GetHighBinBorderCentrality(iCentrality));
      } else {
        centralityString = "Pythia8";
        compactCentralityString = "_pythia8";
      }

      hUnfoldingResponse[iCentrality]->GetXaxis()->SetRangeUser(60,300);
      hUnfoldingResponse[iCentrality]->GetYaxis()->SetRangeUser(60,300);
      drawer->DrawHistogram(hUnfoldingResponse[iCentrality], "Reco jet p_{T}", "Gen jet p_{T}", Form("Response for %s", centralityString.Data()), "colz");

      // Save the figures to a file
      if(saveFigures) {
        gPad->GetCanvas()->SaveAs(Form("figures/jetPtOneDimensionalUnfoldingResponse%s%s.%s", saveComment.Data(), compactCentralityString.Data(), figureFormat.Data()));
      }

    } // Centrality loop
  } // Drawing response matrices

  // Refold the unfolded distribution and compare it with the measured distribution
  if(drawRefoldingTest){

    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){

      // Set the centrality information for legends and figure saving
      if(isPbPbData) {
        centralityString = Form("PbPb data: %.0f-%.0f", unfoldingCard[kUnfoldingFile]->GetLowBinBorderCentrality(iCentrality), unfoldingCard[kUnfoldingFile]->GetHighBinBorderCentrality(iCentrality));
        compactCentralityString = Form("_C%.0f-%.0f", unfoldingCard[kUnfoldingFile]->GetLowBinBorderCentrality(iCentrality), unfoldingCard[kUnfoldingFile]->GetHighBinBorderCentrality(iCentrality));
      } else {
        centralityString = "Pythia8";
        compactCentralityString = "_pythia8";
      }

      // Create a new canvas for the plot
      drawer->CreateSplitCanvas();

      // Draw first the reconstructed distribution
      hReferenceDistributionMeasured[iCentrality]->SetLineColor(kBlack);
      hReferenceDistributionMeasured[iCentrality]->GetXaxis()->SetRangeUser(normalizationRegionLow, normalizationRegionHigh);
      drawer->SetLogY(true);
      drawer->DrawHistogramToUpperPad(hReferenceDistributionMeasured[iCentrality], "Jet p_{T}", "dN/dp_{T}", " ", "");

      // Add the forward folded distributions to the same plot
      for(int iIteration = 0; iIteration < nIterations; iIteration++){
        hForwardFoldedDistribution[iCentrality][iIteration]->SetLineColor(iterationColor[iIteration]);
        hForwardFoldedDistribution[iCentrality][iIteration]->Draw("same");
      } 

      // Add a legend to the figure
      legend = new TLegend(0.25, 0.18- 0.03*nIterations, 0.5, 0.48 + 0.02*nIterations);
      legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);

      legend->AddEntry((TObject*)0, centralityString.Data(), "");

      legend->AddEntry(hReferenceDistributionMeasured[iCentrality], "Measured distribution", "l");
      for(int iIteration = 0; iIteration < nIterations; iIteration++){
        legend->AddEntry(hForwardFoldedDistribution[iCentrality][iIteration], Form("Forward folded (%s)", unfoldName[iUnfoldMethod].ReplaceAll("NIT",Form("%d",iterationKey[iIteration])).Data()), "l");
        unfoldName[iUnfoldMethod].ReplaceAll(Form("%d",iterationKey[iIteration]),"NIT");
      }

      legend->Draw();

      // Draw the ratios to lower pad
      drawer->SetLogY(false);

      for(int iIteration = 0; iIteration < nIterations; iIteration++){
        hForwardFoldedToMeasuredRatio[iCentrality][iIteration]->SetLineColor(iterationColor[iIteration]);
        hForwardFoldedToMeasuredRatio[iCentrality][iIteration]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
        hForwardFoldedToMeasuredRatio[iCentrality][iIteration]->GetXaxis()->SetRangeUser(normalizationRegionLow, normalizationRegionHigh);
        if(iIteration == 0){
          drawer->DrawHistogramToLowerPad(hForwardFoldedToMeasuredRatio[iCentrality][iIteration], "Jet p_{T}", "Ratio to measured", " ", "");
        } else {
          hForwardFoldedToMeasuredRatio[iCentrality][iIteration]->Draw("same");
        }
      }

      oneLine->Draw();

      // Save the figures to a file
      if(saveFigures) {
        gPad->GetCanvas()->SaveAs(Form("figures/jetPtOneDimensionalUnfoldingRefoldingCheck%s%s.%s", saveComment.Data(), compactCentralityString.Data(), figureFormat.Data()));
      }
    } // Centrality loop
  } // Unfolded to truth comparison
}
