#include "SystematicUncertaintyOrganizer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECHistogramManager.h"
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Logarithmic function used to fit the ratio plot
 */
double logFit(double* x, double* par){
  return par[1] * TMath::Log10(x[0]) + par[0];
}

/*
 * Logarithmic third order polynomial
 */
double logPol3(double* x, double* par){
  return par[0] + par[1] * TMath::Log10(x[0]) + par[2] * TMath::Power(TMath::Log10(x[0]),2) + par[3] * TMath::Power(TMath::Log10(x[0]),3);
}


/*
 * Shift the distributions up and down based on the third distribution uncertainties
 */
void shiftByErrors(TH1D* upShiftedHistogram, TH1D* downShiftedHistogram, TH1D* shiftHistogram){
  double binContent, binError;
  int criticalBin = shiftHistogram->GetXaxis()->FindBin(0.055);
  for(int iBin = 1; iBin <= criticalBin; iBin++){
    binError = shiftHistogram->GetBinError(iBin);
    binContent = upShiftedHistogram->GetBinContent(iBin);
    upShiftedHistogram->SetBinContent(iBin,binContent+binError);
    downShiftedHistogram->SetBinContent(iBin,binContent-binError);
  }
  for(int iBin = criticalBin+1; iBin <= shiftHistogram->GetNbinsX(); iBin++){
    binError = shiftHistogram->GetBinError(iBin);
    binContent = upShiftedHistogram->GetBinContent(iBin);
    upShiftedHistogram->SetBinContent(iBin,binContent-binError);
    downShiftedHistogram->SetBinContent(iBin,binContent+binError);
  }
}

/*
 * Find the bigger difference of the two numbers compared to the first one
 */
double getBiggerDifference(double nominal, double shiftUp, double shiftDown){
  double differenceUp = TMath::Abs(nominal - shiftUp);
  double differenceDown = TMath::Abs(nominal - shiftDown);
  if(differenceUp > differenceDown) return differenceUp;
  return differenceDown;
}

/*
 * Macro for making final result plots comparing energy-energy correlators between pp and PbPb
 */
void fitPbPbToPpRatio(){

  enum enumDataType{kPbPb, kPp, kNDataTypes};
  enum enumUncertaintyType{kCorrelatedUncertainty, kUncorrelatedUncertainty, kShiftedUp, kShiftedDown, kNUncertaintyTypes};

  TString uncertaintyString[kNUncertaintyTypes] = {"Correlated", "Regular", "Upshifted", "Downshifted"};
  TString compactUncertaintyString[kNUncertaintyTypes] = {"_correlatedUncertainty", "_uncorrelatedUncertainty", "_upshiftedDistribution", "_downshiftedDistribution"};

  // ============= //
  // Configuration //
  // ============= //
  
  // Input files
  TString inputFileName[kNDataTypes];
  inputFileName[kPbPb] = "data/eecAnalysis_akFlowJet_wtaAxis_newTrackPairEfficiencySmoothed_unfoldingWithNominalSmear_processed_2023-07-13.root";
  inputFileName[kPp] = "data/ppData_pfJets_wtaAxis_newTrackPairEfficiency_unfoldingWithNominalSmear_processed_2023-07-13.root";
  TString uncertaintyFileName[kNDataTypes];
  uncertaintyFileName[kPbPb] = "systematicUncertainties/systematicUncertainties_jetMetUpdate_includeMCnonClosure_2023-07-16.root";
  uncertaintyFileName[kPp] = "systematicUncertainties/systematicUncertaintiesForPp_jetMetUpdate_includeMCnonClosure_2023-07-16.root";
  
  TFile* inputFile[kNDataTypes];
  TFile* uncertaintyFile[kNDataTypes];
  EECCard* card[kNDataTypes];
  EECCard* uncertaintyCard[kNDataTypes];
  for(int iFile = 0; iFile < kNDataTypes; iFile++){

    // Load the input file
    inputFile[iFile] = TFile::Open(inputFileName[iFile]);

    // Check that the input file exists
    if(inputFile[iFile] == NULL){
      cout << "Error! The file " << inputFileName[iFile].Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }
    
    // Load the card from the file and read the collision system
    card[iFile] = new EECCard(inputFile[iFile]);

    // Load the uncertainty file
    uncertaintyFile[iFile] = TFile::Open(uncertaintyFileName[iFile]);

    // Check that the uncertianty file exists
    if(uncertaintyFile[iFile] == NULL){
      cout << "Error! The file " << uncertaintyFileName[iFile].Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the systematicUncertainties/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }
    
    // Load the card from the file and read the collision system
    uncertaintyCard[iFile] = new EECCard(uncertaintyFile[iFile]);
  }

  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Find the number of bins from the PbPb card
  const int nCentralityBins = card[kPbPb]->GetNCentralityBins();
  const int nJetPtBinsEEC = card[kPbPb]->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card[kPbPb]->GetNTrackPtBinsEEC();
  
  // The final results are available for all the bins that are unfolded
  int firstDrawnCentralityBin = card[kPbPb]->GetFirstUnfoldedCentralityBin();
  int lastDrawnCentralityBin = 0;
  
  int firstDrawnJetPtBinEEC = card[kPbPb]->GetFirstUnfoldedJetPtBin();
  int lastDrawnJetPtBinEEC = card[kPbPb]->GetLastUnfoldedJetPtBin();
  
  int firstDrawnTrackPtBinEEC = card[kPbPb]->GetFirstUnfoldedTrackPtBin();
  int lastDrawnTrackPtBinEEC = card[kPbPb]->GetLastUnfoldedTrackPtBin();

  // Choose which plots to draw
  bool drawRatiosWithFits = true;
  bool drawRatiosWithAlternativeFits = false;
  bool drawWiggleGraph = false;
  bool printWiggleTable = false;
  bool drawRatioOfFitsAndRatios = false;
  
  // Save the final plots
  const bool saveFigures = false;
  TString saveComment = "_errorIllustration";

  // Illustrate error
  const bool drawErrorIllustration = true;

  // Ratio zoom settings
  std::pair<double, double> analysisDeltaR = std::make_pair(0.006, 0.39); // DeltaR span in which the analysis is done
  std::pair<double, double> ratioZoom = std::make_pair(0.4, 1.6);
  
  // Marker colors and styles
  int markerStylePbPb[] = {kFullSquare, kFullCircle, kFullCross, kFullFourTrianglesPlus};
  int markerStylePp = kFullDiamond;
  int markerColorPbPb[] = {kRed, kBlue, kMagenta, kGreen+3};
  int markerColorPp = kBlack;
  int correlatedColorPbPb[] = {kMagenta, kCyan, kMagenta, kCyan};

  // =============================================== //
  // Read the histograms from the histogram managers //
  // =============================================== //

  // Create histogram managers for result files and systematic uncertainty organizers for systematic uncertainty files
  EECHistogramManager* histograms[kNDataTypes];
  SystematicUncertaintyOrganizer* uncertainties[kNDataTypes];
  
  for(int iDataType = 0; iDataType < kNDataTypes; iDataType++){

    // Create a new histogram manager
    histograms[iDataType] = new EECHistogramManager(inputFile[iDataType], card[iDataType]);
  
    // Load all unfolded energy-energy correlators
    histograms[iDataType]->SetLoadEnergyEnergyCorrelators(true);
    histograms[iDataType]->SetCentralityBinRange(card[iDataType]->GetFirstUnfoldedCentralityBin(), card[iDataType]->GetLastUnfoldedCentralityBin());
    histograms[iDataType]->SetTrackPtBinRangeEEC(card[iDataType]->GetFirstUnfoldedTrackPtBin(), card[iDataType]->GetLastUnfoldedTrackPtBin());
    histograms[iDataType]->SetJetPtBinRangeEEC(card[iDataType]->GetFirstUnfoldedJetPtBin(), card[iDataType]->GetLastUnfoldedJetPtBin());

    // Load the histograms from the file
    histograms[iDataType]->LoadProcessedHistograms();

    // Create a new systematic uncertainty organizer
    uncertainties[iDataType] = new SystematicUncertaintyOrganizer(uncertaintyFile[iDataType]);

  }
 
  const int nFitParts = 3;
  TH1D* energyEnergyCorrelatorSignalPbPb[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorSignalPp[nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorPbPbToPpRatio[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyForPbPb[kNUncertaintyTypes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyForPp[kNUncertaintyTypes][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyPbPbToPpRatio[kNUncertaintyTypes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* distributionRatioToFitRatio[kNUncertaintyTypes][nFitParts][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];

  TF1* fitToRatio[kNUncertaintyTypes][nFitParts][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TF1* alternativeFitToRatio[kNUncertaintyTypes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];

  // Ranges for fits, should be optimized bin by bin
  int fitColor[4][nFitParts] = {{kBlack, kBlue, kBlack},  // 0-10% centrality
                                {kRed, kMagenta, kRed},  // 10-30% centrality
                                {kBlack, kBlue, kBlack},  // 30-50% centrality
                                {kRed, kMagenta, kRed}}; // 50-90% centrality

  double fitLowRange[4][3][4][nFitParts] =
  // Centrality 0-10 %
  {
    // Track pT > 2 GeV
    {{{0.04,0.1,0.24},    // 120 < jet pT < 140 GeV
     {0.04,0.1,0.24},     // 140 < jet pT < 160 GeV
     {0.04,0.085,0.24},     // 160 < jet pT < 180 GeV
     {0.04,0.085,0.24}},    // 180 < jet pT < 200 GeV
    // Track pT > 2.5 GeV
    {{0.04,0.1,0.24},    // 120 < jet pT < 140 GeV
     {0.04,0.1,0.24},     // 140 < jet pT < 160 GeV
     {0.04,0.085,0.24},     // 160 < jet pT < 180 GeV
     {0.04,0.085,0.24}},    // 180 < jet pT < 200 GeV
    // Track pT > 3 GeV
    {{0.04,0.1,0.24},     // 120 < jet pT < 140 GeV
     {0.04,0.1,0.24},     // 140 < jet pT < 160 GeV
     {0.04,0.085,0.24},     // 160 < jet pT < 180 GeV
     {0.04,0.085,0.24}}},   // 180 < jet pT < 200 GeV

  // Centrality 10-30 %
   // Track pT > 2 GeV
    {{{0.04,0.12,0.24},    // 120 < jet pT < 140 GeV
     {0.04,0.1,0.24},     // 140 < jet pT < 160 GeV
     {0.04,0.085,0.24},     // 160 < jet pT < 180 GeV
     {0.04,0.085,0.24}},    // 180 < jet pT < 200 GeV
    // Track pT > 2.5 GeV
    {{0.04,0.12,0.24},    // 120 < jet pT < 140 GeV
     {0.04,0.1,0.24},     // 140 < jet pT < 160 GeV
     {0.04,0.085,0.24},     // 160 < jet pT < 180 GeV
     {0.02,0.085,0.24}},    // 180 < jet pT < 200 GeV
    // Track pT > 3 GeV
    {{0.04,0.12,0.24},     // 120 < jet pT < 140 GeV
     {0.04,0.1,0.24},     // 140 < jet pT < 160 GeV
     {0.04,0.085,0.24},     // 160 < jet pT < 180 GeV
     {0.04,0.085,0.24}}},   // 180 < jet pT < 200 GeV

   // Centrality 30-50 %
   // Track pT > 2 GeV
    {{{0.02,0.1,0.21},    // 120 < jet pT < 140 GeV
     {0.02,0.1,0.21},     // 140 < jet pT < 160 GeV
     {0.02,0.085,0.21},     // 160 < jet pT < 180 GeV
     {0.02,0.085,0.21}},    // 180 < jet pT < 200 GeV
    // Track pT > 2.5 GeV
    {{0.02,0.1,0.21},    // 120 < jet pT < 140 GeV
     {0.02,0.1,0.21},     // 140 < jet pT < 160 GeV
     {0.02,0.085,0.21},     // 160 < jet pT < 180 GeV
     {0.02,0.085,0.21}},    // 180 < jet pT < 200 GeV
    // Track pT > 3 GeV
    {{0.02,0.1,0.21},     // 120 < jet pT < 140 GeV
     {0.02,0.1,0.21},     // 140 < jet pT < 160 GeV
     {0.02,0.085,0.21},     // 160 < jet pT < 180 GeV
     {0.02,0.085,0.21}}},   // 180 < jet pT < 200 GeV

   // Centrality 50-90%
   // Track pT > 2 GeV
    {{{0.02,0.1,0.21},    // 120 < jet pT < 140 GeV
     {0.02,0.1,0.21},     // 140 < jet pT < 160 GeV
     {0.02,0.085,0.21},     // 160 < jet pT < 180 GeV
     {0.02,0.085,0.21}},    // 180 < jet pT < 200 GeV
    // Track pT > 2.5 GeV
    {{0.02,0.1,0.21},    // 120 < jet pT < 140 GeV
     {0.02,0.1,0.21},     // 140 < jet pT < 160 GeV
     {0.02,0.085,0.21},     // 160 < jet pT < 180 GeV
     {0.02,0.085,0.21}},    // 180 < jet pT < 200 GeV
    // Track pT > 3 GeV
    {{0.02,0.1,0.21},     // 120 < jet pT < 140 GeV
     {0.02,0.1,0.21},     // 140 < jet pT < 160 GeV
     {0.02,0.085,0.21},     // 160 < jet pT < 180 GeV
     {0.02,0.085,0.21}}}   // 180 < jet pT < 200 GeV
  };

    double fitHighRange[4][3][4][nFitParts] =
  
  // Centrality 0-10 %
  {
    // Track pT > 2 GeV
    {{{0.06,0.19,0.39},    // 120 < jet pT < 140 GeV
     {0.06,0.16,0.39},     // 140 < jet pT < 160 GeV
     {0.06,0.16,0.39},     // 160 < jet pT < 180 GeV
     {0.06,0.16,0.39}},    // 180 < jet pT < 200 GeV
    // Track pT > 2.5 GeV
    {{0.06,0.19,0.39},    // 120 < jet pT < 140 GeV
     {0.06,0.16,0.39},     // 140 < jet pT < 160 GeV
     {0.06,0.16,0.39},     // 160 < jet pT < 180 GeV
     {0.06,0.16,0.39}},    // 180 < jet pT < 200 GeV
    // Track pT > 3 GeV
    {{0.06,0.19,0.39},     // 120 < jet pT < 140 GeV
     {0.06,0.16,0.39},     // 140 < jet pT < 160 GeV
     {0.06,0.16,0.39},     // 160 < jet pT < 180 GeV
     {0.06,0.16,0.39}}},   // 180 < jet pT < 200 GeV

  // Centrality 10-30 %
   // Track pT > 2 GeV
    {{{0.06,0.19,0.39},    // 120 < jet pT < 140 GeV
     {0.06,0.19,0.39},     // 140 < jet pT < 160 GeV
     {0.06,0.16,0.39},     // 160 < jet pT < 180 GeV
     {0.06,0.14,0.32}},    // 180 < jet pT < 200 GeV
    // Track pT > 2.5 GeV
    {{0.06,0.19,0.39},    // 120 < jet pT < 140 GeV
     {0.06,0.19,0.39},     // 140 < jet pT < 160 GeV
     {0.06,0.16,0.39},     // 160 < jet pT < 180 GeV
     {0.06,0.14,0.32}},    // 180 < jet pT < 200 GeV
    // Track pT > 3 GeV
    {{0.06,0.19,0.39},     // 120 < jet pT < 140 GeV
     {0.06,0.19,0.39},     // 140 < jet pT < 160 GeV
     {0.06,0.16,0.39},     // 160 < jet pT < 180 GeV
     {0.06,0.14,0.32}}},   // 180 < jet pT < 200 GeV

   // Centrality 30-50 %
   // Track pT > 2 GeV
    {{{0.075,0.19,0.39},    // 120 < jet pT < 140 GeV
     {0.075,0.19,0.39},     // 140 < jet pT < 160 GeV
     {0.06,0.16,0.39},     // 160 < jet pT < 180 GeV
     {0.06,0.16,0.39}},    // 180 < jet pT < 200 GeV
    // Track pT > 2.5 GeV
    {{0.075,0.19,0.39},    // 120 < jet pT < 140 GeV
     {0.075,0.19,0.39},     // 140 < jet pT < 160 GeV
     {0.06,0.16,0.39},     // 160 < jet pT < 180 GeV
     {0.06,0.16,0.39}},    // 180 < jet pT < 200 GeV
    // Track pT > 3 GeV
    {{0.075,0.19,0.39},     // 120 < jet pT < 140 GeV
     {0.075,0.19,0.39},     // 140 < jet pT < 160 GeV
     {0.06,0.16,0.39},     // 160 < jet pT < 180 GeV
     {0.06,0.16,0.39}}},   // 180 < jet pT < 200 GeV

   // Centrality 50-90%
   // Track pT > 2 GeV
    {{{0.075,0.19,0.39},    // 120 < jet pT < 140 GeV
     {0.075,0.19,0.39},     // 140 < jet pT < 160 GeV
     {0.06,0.16,0.39},     // 160 < jet pT < 180 GeV
     {0.06,0.16,0.39}},    // 180 < jet pT < 200 GeV
    // Track pT > 2.5 GeV
    {{0.075,0.19,0.39},    // 120 < jet pT < 140 GeV
     {0.075,0.19,0.39},     // 140 < jet pT < 160 GeV
     {0.06,0.16,0.39},     // 160 < jet pT < 180 GeV
     {0.06,0.16,0.39}},    // 180 < jet pT < 200 GeV
    // Track pT > 3 GeV
    {{0.075,0.19,0.39},     // 120 < jet pT < 140 GeV
     {0.075,0.19,0.39},     // 140 < jet pT < 160 GeV
     {0.06,0.16,0.39},     // 160 < jet pT < 180 GeV
     {0.06,0.16,0.39}}},   // 180 < jet pT < 200 GeV
  };

  // Initialize histograms to NULL
  for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
      energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt] = NULL;

      for(int iUncertainty = 0; iUncertainty < kNUncertaintyTypes; iUncertainty++){
        systematicUncertaintyForPp[iUncertainty][iJetPt][iTrackPt] = NULL;
      }
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt] = NULL;
        energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt] = NULL;

        for(int iUncertainty = 0; iUncertainty < kNUncertaintyTypes; iUncertainty++){
          systematicUncertaintyForPbPb[iUncertainty][iCentrality][iJetPt][iTrackPt] = NULL;
          systematicUncertaintyPbPbToPpRatio[iUncertainty][iCentrality][iJetPt][iTrackPt] = NULL;

          // Initialize the fit functions to NULL
          for(int iFit = 0; iFit < nFitParts; iFit++){
            fitToRatio[iUncertainty][iFit][iCentrality][iJetPt][iTrackPt] = NULL;
            distributionRatioToFitRatio[iUncertainty][iFit][iCentrality][iJetPt][iTrackPt] = NULL;
          }
          alternativeFitToRatio[iUncertainty][iCentrality][iJetPt][iTrackPt] = NULL;

        } // Uncertainty loop
      } // Centrality loop
    } // Track pT loop
  } // Jet pT loop

  // Read the histograms from managers
  int iTrackPtMatchedPp, iTrackPtMatchedPbPbUncertainty, iTrackPtMatchedPpUncertainty;
  int iJetPtMatchedPp, iJetPtMatchedPbPbUncertainty, iJetPtMatchedPpUncertainty;
  int iCentralityMatched;

  for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
    iJetPtMatchedPp = card[kPp]->FindBinIndexJetPtEEC(card[kPbPb]->GetBinBordersJetPtEEC(iJetPt));
    iJetPtMatchedPbPbUncertainty = uncertaintyCard[kPbPb]->FindBinIndexJetPtEEC(card[kPbPb]->GetBinBordersJetPtEEC(iJetPt));
    iJetPtMatchedPpUncertainty = uncertaintyCard[kPp]->FindBinIndexJetPtEEC(card[kPbPb]->GetBinBordersJetPtEEC(iJetPt));
    for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
      iTrackPtMatchedPp = card[kPp]->FindBinIndexTrackPtEEC(card[kPbPb]->GetBinBordersTrackPtEEC(iTrackPt));
      iTrackPtMatchedPbPbUncertainty = uncertaintyCard[kPbPb]->FindBinIndexTrackPtEEC(card[kPbPb]->GetBinBordersTrackPtEEC(iTrackPt));
      iTrackPtMatchedPpUncertainty = uncertaintyCard[kPp]->FindBinIndexTrackPtEEC(card[kPbPb]->GetBinBordersTrackPtEEC(iTrackPt));

      // Read the pp histograms that do not have centrality binning
      energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt] = histograms[kPp]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatchedPp, iTrackPtMatchedPp, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);
      systematicUncertaintyForPp[kCorrelatedUncertainty][iJetPt][iTrackPt] = uncertainties[kPp]->GetCorrelatedSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty);
      systematicUncertaintyForPp[kUncorrelatedUncertainty][iJetPt][iTrackPt] = uncertainties[kPp]->GetUncorrelatedSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty);

      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        iCentralityMatched = uncertaintyCard[kPbPb]->FindBinIndexCentrality(card[kPbPb]->GetBinBordersCentrality(iCentrality));

        // Read the PbPb histograms
        energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt] = histograms[kPbPb]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);
        systematicUncertaintyForPbPb[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt] = uncertainties[kPbPb]->GetCorrelatedSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty);
        systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt] = uncertainties[kPbPb]->GetUncorrelatedSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty);

        // Initialize the fit functions with the defined function
        for(int iFit = 0; iFit < nFitParts; iFit++){
          for(int iUncertainty = 0; iUncertainty <kNUncertaintyTypes; iUncertainty++){
            fitToRatio[iUncertainty][iFit][iCentrality][iJetPt][iTrackPt] = new TF1(Form("fit%d%d%d%d%d", iUncertainty, iCentrality, iJetPt, iTrackPt, iFit), logFit, fitLowRange[iCentrality - firstDrawnCentralityBin][iTrackPt - firstDrawnTrackPtBinEEC][iJetPt - firstDrawnJetPtBinEEC][iFit], fitHighRange[iCentrality - firstDrawnCentralityBin][iTrackPt - firstDrawnTrackPtBinEEC][iJetPt - firstDrawnJetPtBinEEC][iFit], 2);
            fitToRatio[iUncertainty][iFit][iCentrality][iJetPt][iTrackPt]->SetParameters(-1, 1);
            fitToRatio[iUncertainty][iFit][iCentrality][iJetPt][iTrackPt]->SetLineColor(fitColor[iCentrality][iFit]);
          }
        }

        // Initialize the alternative fit functions
        for(int iUncertainty = 0; iUncertainty <kNUncertaintyTypes; iUncertainty++){
          alternativeFitToRatio[iUncertainty][iCentrality][iJetPt][iTrackPt] = new TF1(Form("alternativeFit%d%d%d%d", iUncertainty, iCentrality, iJetPt, iTrackPt), logPol3, 0.02, 0.39, 4);
          alternativeFitToRatio[iUncertainty][iCentrality][iJetPt][iTrackPt]->SetParameters(0.3, -1.5, -1.5, -0.5);
          alternativeFitToRatio[iUncertainty][iCentrality][iJetPt][iTrackPt]->SetLineColor(fitColor[iCentrality][0]);
        }
      } // Centrality loop
    } // Track pT loop
  } // Jet pT loop

  // ================================================ //
  //   Normalize the histograms and calculate ratios  //
  // ================================================ //

  double epsilon = 0.0001;
  int lowAnalysisBin = energyEnergyCorrelatorSignalPbPb[firstDrawnCentralityBin][firstDrawnJetPtBinEEC][firstDrawnTrackPtBinEEC]->GetXaxis()->FindBin(analysisDeltaR.first + epsilon);
  int highAnalysisBin = energyEnergyCorrelatorSignalPbPb[firstDrawnCentralityBin][firstDrawnJetPtBinEEC][firstDrawnTrackPtBinEEC]->GetXaxis()->FindBin(analysisDeltaR.second - epsilon);
  for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
    for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
      energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->Scale(1.0 / energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
      systematicUncertaintyForPp[kUncorrelatedUncertainty][iJetPt][iTrackPt]->Scale(1.0 / systematicUncertaintyForPp[kUncorrelatedUncertainty][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
      systematicUncertaintyForPp[kCorrelatedUncertainty][iJetPt][iTrackPt]->Scale(1.0 / systematicUncertaintyForPp[kCorrelatedUncertainty][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

      // Shift the points in the uncorrelated ratio by the correlated uncertainties
      systematicUncertaintyForPp[kShiftedUp][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPp[kUncorrelatedUncertainty][iJetPt][iTrackPt]->Clone(Form("ppShiftedUp%d%d", iJetPt, iTrackPt));
      systematicUncertaintyForPp[kShiftedDown][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPp[kUncorrelatedUncertainty][iJetPt][iTrackPt]->Clone(Form("ppShiftedDown%d%d", iJetPt, iTrackPt));
      shiftByErrors(systematicUncertaintyForPp[kShiftedUp][iJetPt][iTrackPt], systematicUncertaintyForPp[kShiftedDown][iJetPt][iTrackPt], systematicUncertaintyForPp[kCorrelatedUncertainty][iJetPt][iTrackPt]);

      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->Scale(1.0 / energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
        systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Scale(1.0 / systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
        systematicUncertaintyForPbPb[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Scale(1.0 / systematicUncertaintyForPbPb[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

        // Shift the points in the uncorrelated ratio by the correlated uncertainties
        systematicUncertaintyForPbPb[kShiftedUp][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("PbPbShiftedUp%d%d%d", iCentrality, iJetPt, iTrackPt));
        systematicUncertaintyForPbPb[kShiftedDown][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("PbPbShiftedDown%d%d%d", iCentrality, iJetPt, iTrackPt));
        shiftByErrors(systematicUncertaintyForPbPb[kShiftedUp][iCentrality][iJetPt][iTrackPt], systematicUncertaintyForPbPb[kShiftedDown][iCentrality][iJetPt][iTrackPt], systematicUncertaintyForPbPb[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]);

        energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt] = (TH1D*) energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->Clone(Form("energyEnergyCorrelatorRatio%d%d%d", iCentrality, iJetPt, iTrackPt));
        energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->Divide(energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]);

        systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("systematicUncertaintyUncorrelatedRatio%d%d%d", iCentrality, iJetPt, iTrackPt));
        systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Divide(systematicUncertaintyForPp[kUncorrelatedUncertainty][iJetPt][iTrackPt]);
        systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPbPb[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("systematicUncertaintyUncorrelatedRatio%d%d%d", iCentrality, iJetPt, iTrackPt));
        systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Divide(systematicUncertaintyForPp[kCorrelatedUncertainty][iJetPt][iTrackPt]);

        // Shift the points in the uncorrelated ratio by the correlated uncertainties
        systematicUncertaintyPbPbToPpRatio[kShiftedUp][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("ratioShiftedUp%d%d%d", iCentrality, iJetPt, iTrackPt));
        systematicUncertaintyPbPbToPpRatio[kShiftedDown][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("ratioShiftedDown%d%d%d", iCentrality, iJetPt, iTrackPt));
        shiftByErrors(systematicUncertaintyPbPbToPpRatio[kShiftedUp][iCentrality][iJetPt][iTrackPt], systematicUncertaintyPbPbToPpRatio[kShiftedDown][iCentrality][iJetPt][iTrackPt], systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]);

        // Fit the thingy
        for(int iUncertainty = 0; iUncertainty < kNUncertaintyTypes; iUncertainty++){
          for(int iFit = 0; iFit < nFitParts; iFit++){
            systematicUncertaintyPbPbToPpRatio[iUncertainty][iCentrality][iJetPt][iTrackPt]->Fit(fitToRatio[iUncertainty][iFit][iCentrality][iJetPt][iTrackPt], "NQ", "", fitLowRange[iCentrality - firstDrawnCentralityBin][iTrackPt - firstDrawnTrackPtBinEEC][iJetPt - firstDrawnJetPtBinEEC][iFit], fitHighRange[iCentrality - firstDrawnCentralityBin][iTrackPt - firstDrawnTrackPtBinEEC][iJetPt - firstDrawnJetPtBinEEC][iFit]);
          }

          systematicUncertaintyPbPbToPpRatio[iUncertainty][iCentrality][iJetPt][iTrackPt]->Fit(alternativeFitToRatio[iUncertainty][iCentrality][iJetPt][iTrackPt], "NQ", "", 0.02, 0.39);
        }

      } // Centrality loop
    } // Track pT loop
  } // Jet pT loop

  // Find the deltaR value corresponding to the points where the different fit region overlap
  double firstTurningPoint[kNUncertaintyTypes+1][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  double secondTurningPoint[kNUncertaintyTypes+1][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  double firstTurningPointError[kNUncertaintyTypes+1][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  double secondTurningPointError[kNUncertaintyTypes+1][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  for(int iUncertainty = 0; iUncertainty < kNUncertaintyTypes+1; iUncertainty++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
          firstTurningPoint[iUncertainty][iCentrality][iJetPt][iTrackPt] = 0;
          secondTurningPoint[iUncertainty][iCentrality][iJetPt][iTrackPt] = 0;
          firstTurningPointError[iUncertainty][iCentrality][iJetPt][iTrackPt] = 0;
          secondTurningPointError[iUncertainty][iCentrality][iJetPt][iTrackPt] = 0;
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Uncertainty loop

  double a,b,c,d;
  double ea, eb, ec, ed;
  for(int iUncertainty = 0; iUncertainty < kNUncertaintyTypes; iUncertainty++){
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
        for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
        
          // Calculate the x-value for which the two functions overlap
          // a*log(x)+b = c*log(x)+d => log(x) = (d - b) / (a - c) => x = 10^[(d - b) / (a - c)]
          a = fitToRatio[iUncertainty][0][iCentrality][iJetPt][iTrackPt]->GetParameter(1);
          b = fitToRatio[iUncertainty][0][iCentrality][iJetPt][iTrackPt]->GetParameter(0);
          c = fitToRatio[iUncertainty][1][iCentrality][iJetPt][iTrackPt]->GetParameter(1);
          d = fitToRatio[iUncertainty][1][iCentrality][iJetPt][iTrackPt]->GetParameter(0);
          firstTurningPoint[iUncertainty][iCentrality][iJetPt][iTrackPt] = TMath::Power(10.0,((d-b)/(a-c)));

          // Calculate the error using the law of error propagation
          // Chain rule: d/dx f(g(x)) = f'(g(x)) * g'(x)
          // d/da (d-b)/(a-c) E(a) = (d-b) d/da (a-c)^{-1} E(a) = (b-d)/(a-c)^2 E(a)
          // d/db (d-b)/(a-c) E(b) = -1/(a-c) E(b)
          // d/dc (d-b)/(a-c) E(c) = (d-b)/(a-c)^2 E(c)
          // d/dd (d-b)/(a-c) E(d) = 1/(a-c) E(d)
          // d/dx a^x = a^x * ln(a)
          ea = (b-d)/(TMath::Power(a-c,2)) * fitToRatio[iUncertainty][0][iCentrality][iJetPt][iTrackPt]->GetParError(1);
          eb = -1/(a-c) * fitToRatio[iUncertainty][0][iCentrality][iJetPt][iTrackPt]->GetParError(0);
          ec = (d-b)/(TMath::Power(a-c,2)) * fitToRatio[iUncertainty][1][iCentrality][iJetPt][iTrackPt]->GetParError(1);
          ed = 1/(a-c) * fitToRatio[iUncertainty][1][iCentrality][iJetPt][iTrackPt]->GetParError(0);
          firstTurningPointError[iUncertainty][iCentrality][iJetPt][iTrackPt] = TMath::Sqrt(TMath::Power(firstTurningPoint[iUncertainty][iCentrality][iJetPt][iTrackPt] * TMath::Log(10) * ea,2) + TMath::Power(firstTurningPoint[iUncertainty][iCentrality][iJetPt][iTrackPt] * TMath::Log(10) * eb,2) + TMath::Power(firstTurningPoint[iUncertainty][iCentrality][iJetPt][iTrackPt] * TMath::Log(10) * ec,2) + TMath::Power(firstTurningPoint[iUncertainty][iCentrality][iJetPt][iTrackPt] * TMath::Log(10) * ed,2));

          // Do the same for the second turning point
          a = fitToRatio[iUncertainty][2][iCentrality][iJetPt][iTrackPt]->GetParameter(1);
          b = fitToRatio[iUncertainty][2][iCentrality][iJetPt][iTrackPt]->GetParameter(0);
          secondTurningPoint[iUncertainty][iCentrality][iJetPt][iTrackPt] = TMath::Power(10.0,((d-b)/(a-c)));

          ea = (b-d)/(TMath::Power(a-c,2)) * fitToRatio[iUncertainty][2][iCentrality][iJetPt][iTrackPt]->GetParError(1);
          eb = -1/(a-c) * fitToRatio[iUncertainty][0][iCentrality][iJetPt][iTrackPt]->GetParError(0);
          ec = (d-b)/(TMath::Power(a-c,2)) * fitToRatio[iUncertainty][2][iCentrality][iJetPt][iTrackPt]->GetParError(1);
          ed = 1/(a-c) * fitToRatio[iUncertainty][1][iCentrality][iJetPt][iTrackPt]->GetParError(0);
          secondTurningPointError[iUncertainty][iCentrality][iJetPt][iTrackPt] = TMath::Sqrt(TMath::Power(secondTurningPoint[iUncertainty][iCentrality][iJetPt][iTrackPt] * TMath::Log(10) * ea,2) + TMath::Power(secondTurningPoint[iUncertainty][iCentrality][iJetPt][iTrackPt] * TMath::Log(10) * eb,2) + TMath::Power(secondTurningPoint[iUncertainty][iCentrality][iJetPt][iTrackPt] * TMath::Log(10) * ec,2) + TMath::Power(secondTurningPoint[iUncertainty][iCentrality][iJetPt][iTrackPt] * TMath::Log(10) * ed,2));

        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Uncertainty loop

  // Calculate the uncertainty manually using the difference in values between nominal and shifted results
  for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
    for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
      for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){

        firstTurningPoint[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt] = firstTurningPoint[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt];
        secondTurningPoint[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt] = secondTurningPoint[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt];
        firstTurningPointError[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt] = getBiggerDifference(firstTurningPoint[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], firstTurningPoint[kShiftedUp][iCentrality][iJetPt][iTrackPt], firstTurningPoint[kShiftedDown][iCentrality][iJetPt][iTrackPt]);
        secondTurningPointError[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt] = getBiggerDifference(secondTurningPoint[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], secondTurningPoint[kShiftedUp][iCentrality][iJetPt][iTrackPt], secondTurningPoint[kShiftedDown][iCentrality][iJetPt][iTrackPt]);

      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

  // Make ratio between the linear fits and the distribution to which those fits are made
  for(int iUncertainty = 0; iUncertainty < kNUncertaintyTypes; iUncertainty++){
    for(int iFit = 0; iFit < nFitParts; iFit++){
      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
          for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
            distributionRatioToFitRatio[iUncertainty][iFit][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyPbPbToPpRatio[iUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("ratioBetweenRatioAndTheFitToTheSaidRatio%d%d%d%d%d",iUncertainty,iFit,iCentrality,iJetPt,iTrackPt));
            distributionRatioToFitRatio[iUncertainty][iFit][iCentrality][iJetPt][iTrackPt]->Divide(fitToRatio[iUncertainty][iFit][iCentrality][iJetPt][iTrackPt]);
          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop 
    }  // Fit part loop
  } // Uncertainty type loop


  // Create graphs with the locations of first and second wiggling angles
  TGraph *wiggleGraph[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  double yPointsWiggleGraph[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][2];
  double xPointsWiggleGraph[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][2];
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        wiggleGraph[iCentrality][iJetPt][iTrackPt] = NULL;
        for(int iWiggle = 0; iWiggle < 2; iWiggle++){
          xPointsWiggleGraph[iCentrality][iJetPt][iTrackPt][iWiggle] = 0;
          yPointsWiggleGraph[iCentrality][iJetPt][iTrackPt][iWiggle] = 0;
        } // Wiggle loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

  for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
    for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
      for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
        yPointsWiggleGraph[iCentrality][iJetPt][iTrackPt][0] = 0.2 + 0.2 * (iJetPt - firstDrawnJetPtBinEEC);
        yPointsWiggleGraph[iCentrality][iJetPt][iTrackPt][1] = 0.2 + 0.2 * (iJetPt - firstDrawnJetPtBinEEC);
        xPointsWiggleGraph[iCentrality][iJetPt][iTrackPt][0] = firstTurningPoint[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt];
        xPointsWiggleGraph[iCentrality][iJetPt][iTrackPt][1] = secondTurningPoint[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt];
        wiggleGraph[iCentrality][iJetPt][iTrackPt] = new TGraph(2, xPointsWiggleGraph[iCentrality][iJetPt][iTrackPt], yPointsWiggleGraph[iCentrality][iJetPt][iTrackPt]);
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

  // ========================================= //
  //   Draw the histograms in separate plots   //
  // ========================================= //

  JDrawer* drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  drawer->SetLogX(true); // Logarithmic deltaR axis

  TString centralityString;
  TString compactCentralityString;
  TString jetPtString;
  TString compactJetPtString;
  TString trackPtString;
  TString compactTrackPtString;

  // Common variables for different plots
  TLegend* legend;
  TLatex* mainTitle;
  TLine* errorDrawer = new TLine(0,1,1,1); errorDrawer->SetLineStyle(2);

  // Draw individual plots with all centralities mixed together in a single figure
  if(drawRatiosWithFits){

    for(int iUncertainty = kUncorrelatedUncertainty; iUncertainty < kNUncertaintyTypes; iUncertainty++){
      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        centralityString = Form("Cent: %.0f-%.0f%%", card[kPbPb]->GetLowBinBorderCentrality(iCentrality), card[kPbPb]->GetHighBinBorderCentrality(iCentrality));
        compactCentralityString = Form("_C=%.0f-%.0f", card[kPbPb]->GetLowBinBorderCentrality(iCentrality), card[kPbPb]->GetHighBinBorderCentrality(iCentrality));

        for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++) {
          jetPtString = Form("%.0f < jet p_{T} < %.0f", card[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));
          compactJetPtString = Form("_J=%.0f-%.0f", card[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));

          for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++) {
            trackPtString = Form("%.1f < track p_{T}", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));
            compactTrackPtString = Form("_T>%.1f", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));
            compactTrackPtString.ReplaceAll(".", "v");

            // Setup the legend for plots
            legend = new TLegend(0.53, 0.6, 0.83, 0.85);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);
            legend->AddEntry((TObject*)0, uncertaintyString[iUncertainty].Data(), "");
            legend->AddEntry((TObject*)0, centralityString.Data(), "");
            legend->AddEntry((TObject*)0, jetPtString.Data(), "");
            legend->AddEntry((TObject*)0, trackPtString.Data(), "");

            // Linear scale for the ratio
            drawer->SetLogY(false);

            // Set the axis drawing ranges
            systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
            systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

            // Set the style for histograms
            for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++) {
              systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(correlatedColorPbPb[iCentrality - firstDrawnCentralityBin], 0.4);
              systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);
              systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);

              systematicUncertaintyPbPbToPpRatio[iUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin]);
              systematicUncertaintyPbPbToPpRatio[iUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);

              systematicUncertaintyPbPbToPpRatio[iUncertainty][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin], 0.4);
              systematicUncertaintyPbPbToPpRatio[iUncertainty][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
              systematicUncertaintyPbPbToPpRatio[iUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
            }

            drawer->SetGridY(true);
            drawer->DrawHistogram(systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], "#Deltar", "#frac{PbPb}{pp}", " ", "e3");
            systematicUncertaintyPbPbToPpRatio[iUncertainty][iCentrality][iJetPt][iTrackPt]->Draw("same,p,e2");

            // Draw the fits
            for(int iFit = 0; iFit < nFitParts; iFit++) {
              fitToRatio[iUncertainty][iFit][iCentrality][iJetPt][iTrackPt]->SetLineWidth(5);
              fitToRatio[iUncertainty][iFit][iCentrality][iJetPt][iTrackPt]->Draw("same");
            }
            drawer->SetGridY(false);

            // Draw the legend
            legend->Draw();

            // Illustrate the error estimate in the plot
            if(drawErrorIllustration){

              if(iUncertainty == kUncorrelatedUncertainty){
                errorDrawer->SetLineColor(kMagenta);
                errorDrawer->DrawLine(firstTurningPoint[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt], 0.7, firstTurningPoint[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt], 1.1);
                errorDrawer->SetLineColor(kRed);
                errorDrawer->DrawLine(firstTurningPoint[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt]-firstTurningPointError[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt], 0.7, firstTurningPoint[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt]-firstTurningPointError[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt], 1.1);
                errorDrawer->DrawLine(firstTurningPoint[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt]+firstTurningPointError[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt], 0.7, firstTurningPoint[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt]+firstTurningPointError[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt], 1.1);
                errorDrawer->SetLineColor(kCyan);
                errorDrawer->DrawLine(secondTurningPoint[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt], 0.7, secondTurningPoint[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt], 1.1);
                errorDrawer->SetLineColor(kBlue);
                errorDrawer->DrawLine(secondTurningPoint[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt]-secondTurningPointError[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt], 0.7, secondTurningPoint[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt]-secondTurningPointError[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt], 1.1);
                errorDrawer->DrawLine(secondTurningPoint[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt]+secondTurningPointError[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt], 0.7, secondTurningPoint[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt]+secondTurningPointError[kNUncertaintyTypes][iCentrality][iJetPt][iTrackPt], 1.1);
              } else {
                errorDrawer->SetLineColor(kMagenta);
                errorDrawer->DrawLine(firstTurningPoint[iUncertainty][iCentrality][iJetPt][iTrackPt], 0.7, firstTurningPoint[iUncertainty][iCentrality][iJetPt][iTrackPt], 1.1);
                errorDrawer->SetLineColor(kCyan);
                errorDrawer->DrawLine(secondTurningPoint[iUncertainty][iCentrality][iJetPt][iTrackPt], 0.7, secondTurningPoint[iUncertainty][iCentrality][iJetPt][iTrackPt], 1.1);
              }
            }


            // If a plot name is given, save the plot in a file
            if(saveFigures) {
              gPad->GetCanvas()->SaveAs(Form("figures/energyEnergyCorrelator_ratioFit%s%s%s%s%s.pdf", saveComment.Data(), compactUncertaintyString[iUncertainty].Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data()));
            }

          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Uncertainty loop
  } // Drawing individual canvases with all centralities in one canvas

  // Draw individual plots with all centralities mixed together in a single figure
  if(drawRatiosWithAlternativeFits){

    for(int iUncertainty = kUncorrelatedUncertainty; iUncertainty < kNUncertaintyTypes; iUncertainty++){
      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        centralityString = Form("Cent: %.0f-%.0f%%", card[kPbPb]->GetLowBinBorderCentrality(iCentrality), card[kPbPb]->GetHighBinBorderCentrality(iCentrality));
        compactCentralityString = Form("_C=%.0f-%.0f", card[kPbPb]->GetLowBinBorderCentrality(iCentrality), card[kPbPb]->GetHighBinBorderCentrality(iCentrality));

        for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++) {
          jetPtString = Form("%.0f < jet p_{T} < %.0f", card[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));
          compactJetPtString = Form("_J=%.0f-%.0f", card[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));

          for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++) {
            trackPtString = Form("%.1f < track p_{T}", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));
            compactTrackPtString = Form("_T>%.1f", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));
            compactTrackPtString.ReplaceAll(".", "v");

            // Setup the legend for plots
            legend = new TLegend(0.53, 0.6, 0.83, 0.85);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);
            legend->AddEntry((TObject*)0, uncertaintyString[iUncertainty].Data(), "");
            legend->AddEntry((TObject*)0, centralityString.Data(), "");
            legend->AddEntry((TObject*)0, jetPtString.Data(), "");
            legend->AddEntry((TObject*)0, trackPtString.Data(), "");

            // Linear scale for the ratio
            drawer->SetLogY(false);

            // Set the axis drawing ranges
            systematicUncertaintyPbPbToPpRatio[iUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
            systematicUncertaintyPbPbToPpRatio[iUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

            // Set the style for histograms
            for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++) {
              systematicUncertaintyPbPbToPpRatio[iUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin]);
              systematicUncertaintyPbPbToPpRatio[iUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);

              systematicUncertaintyPbPbToPpRatio[iUncertainty][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin], 0.4);
              systematicUncertaintyPbPbToPpRatio[iUncertainty][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
              systematicUncertaintyPbPbToPpRatio[iUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
            }

            drawer->SetGridY(true);
            drawer->DrawHistogram(systematicUncertaintyPbPbToPpRatio[iUncertainty][iCentrality][iJetPt][iTrackPt], "#Deltar", "#frac{PbPb}{pp}", " ", "p,e2");

            // Draw the alternative fit
            alternativeFitToRatio[iUncertainty][iCentrality][iJetPt][iTrackPt]->SetLineWidth(5);
            alternativeFitToRatio[iUncertainty][iCentrality][iJetPt][iTrackPt]->Draw("same");
            drawer->SetGridY(false);

            // Draw the legend
            legend->Draw();

            // If a plot name is given, save the plot in a file
            if(saveFigures) {
              gPad->GetCanvas()->SaveAs(Form("figures/energyEnergyCorrelator_ratioAlternativeFit%s%s%s%s%s.pdf", saveComment.Data(), compactUncertaintyString[iUncertainty].Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data()));
            }

          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Uncertainty loop
  } // Drawing individual canvases with all centralities in one canvas

  // Draw ratios of ratios and fits to ratios
  if(drawRatioOfFitsAndRatios){

    for(int iFit = 0; iFit < nFitParts; iFit++){
      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        centralityString = Form("Cent: %.0f-%.0f%%", card[kPbPb]->GetLowBinBorderCentrality(iCentrality), card[kPbPb]->GetHighBinBorderCentrality(iCentrality));
        compactCentralityString = Form("_C=%.0f-%.0f", card[kPbPb]->GetLowBinBorderCentrality(iCentrality), card[kPbPb]->GetHighBinBorderCentrality(iCentrality));

        for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++) {
          jetPtString = Form("%.0f < jet p_{T} < %.0f", card[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));
          compactJetPtString = Form("_J=%.0f-%.0f", card[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));

          for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++) {
            trackPtString = Form("%.1f < track p_{T}", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));
            compactTrackPtString = Form("_T>%.1f", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));
            compactTrackPtString.ReplaceAll(".", "v");

            // Setup the legend for plots
            legend = new TLegend(0.53, 0.6, 0.83, 0.85);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);
            legend->AddEntry((TObject*)0, Form("Fit region %d", iFit), "");
            legend->AddEntry((TObject*)0, centralityString.Data(), "");
            legend->AddEntry((TObject*)0, jetPtString.Data(), "");
            legend->AddEntry((TObject*)0, trackPtString.Data(), "");

            // Linear scale for the ratio
            drawer->SetLogY(false);

            // Set the axis drawing ranges
            distributionRatioToFitRatio[kUncorrelatedUncertainty][iFit][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
            distributionRatioToFitRatio[kUncorrelatedUncertainty][iFit][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
            distributionRatioToFitRatio[kUncorrelatedUncertainty][iFit][iCentrality][iJetPt][iTrackPt]->SetLineColor(kBlack);

            drawer->SetGridY(true);
            drawer->DrawHistogram(distributionRatioToFitRatio[kUncorrelatedUncertainty][iFit][iCentrality][iJetPt][iTrackPt], "#Deltar", "#frac{Ratio}{Fit}", " ");

            // Draw the alternative fit
            distributionRatioToFitRatio[kShiftedUp][iFit][iCentrality][iJetPt][iTrackPt]->SetLineColor(kBlue);
            distributionRatioToFitRatio[kShiftedUp][iFit][iCentrality][iJetPt][iTrackPt]->Draw("same");
            distributionRatioToFitRatio[kShiftedDown][iFit][iCentrality][iJetPt][iTrackPt]->SetLineColor(kRed);
            distributionRatioToFitRatio[kShiftedDown][iFit][iCentrality][iJetPt][iTrackPt]->Draw("same");
            drawer->SetGridY(false);

            for(int iUncertainty = kUncorrelatedUncertainty; iUncertainty <= kShiftedDown; iUncertainty++){
              legend->AddEntry(distributionRatioToFitRatio[iUncertainty][iFit][iCentrality][iJetPt][iTrackPt], uncertaintyString[iUncertainty].Data(), "l");
            }

            // Draw the legend
            legend->Draw();

            // If a plot name is given, save the plot in a file
            if(saveFigures) {
              gPad->GetCanvas()->SaveAs(Form("figures/energyEnergyCorrelator_ratioOfFitAndRatio%s_fit%d%s%s%s.pdf", saveComment.Data(), iFit, compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data()));
            }

          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Uncertainty loop
  } // Drawing individual canvases with all centralities in one canvas

  // Draw graphs illustrating wiggly behavior
  if(drawWiggleGraph){

    drawer->SetDefaultAppearanceGraph();
    drawer->SetNDivisionsX(510);
    drawer->SetNDivisionsY(510);
    drawer->SetBottomMargin(0.16);
    drawer->SetTitleOffsetX(1.3);
    drawer->SetLabelOffsetX(0.02);
    drawer->SetTitleOffsetY(1.77);
    drawer->SetLabelOffsetY(0.01);

    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      centralityString = Form("Cent: %.0f-%.0f%%", card[kPbPb]->GetLowBinBorderCentrality(iCentrality), card[kPbPb]->GetHighBinBorderCentrality(iCentrality));
      compactCentralityString = Form("_C=%.0f-%.0f", card[kPbPb]->GetLowBinBorderCentrality(iCentrality), card[kPbPb]->GetHighBinBorderCentrality(iCentrality));
      for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
        trackPtString = Form("%.1f < track p_{T}", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));
          compactTrackPtString = Form("_T>%.1f", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));
          compactTrackPtString.ReplaceAll(".", "v");

          // Setup the legend for plots
          legend = new TLegend(0.48, 0.2, 0.78, 0.85);
          legend->SetFillStyle(0);
          legend->SetBorderSize(0);
          legend->SetTextSize(0.05);
          legend->SetTextFont(62);
          legend->AddEntry((TObject*)0, centralityString.Data(), "");
          legend->AddEntry((TObject*)0, trackPtString.Data(), "");

          // Set the style for graphs
          for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
            wiggleGraph[iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iJetPt - firstDrawnJetPtBinEEC]);
            wiggleGraph[iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);
            wiggleGraph[iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iJetPt - firstDrawnJetPtBinEEC]);
            wiggleGraph[iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iJetPt - firstDrawnJetPtBinEEC]);

            jetPtString = Form("%.0f < jet p_{T} < %.0f", card[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));
            legend->AddEntry(wiggleGraph[iCentrality][iJetPt][iTrackPt], jetPtString.Data(), "p");
          }

          // Add ghost points to the first graph to properly adjust drawing range in x-direction
          wiggleGraph[iCentrality][firstDrawnJetPtBinEEC][iTrackPt]->AddPoint(-1,0);
          wiggleGraph[iCentrality][firstDrawnJetPtBinEEC][iTrackPt]->AddPoint(1,0);

          // Draw the graphs to canvas
          drawer->DrawGraphCustomAxes(wiggleGraph[iCentrality][firstDrawnJetPtBinEEC][iTrackPt], 0, 0.7, 0, 1, "#Deltar", "Jet p_{T} bin", " ", "a,p");

          for(int iJetPt = firstDrawnJetPtBinEEC+1; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
            wiggleGraph[iCentrality][iJetPt][iTrackPt]->Draw("p,same");
          }

          legend->Draw();

          // If a plot name is given, save the plot in a file
          if(saveFigures) {
            gPad->GetCanvas()->SaveAs(Form("figures/wiggleGrahps%s%s%s.pdf", saveComment.Data(), compactCentralityString.Data(), compactTrackPtString.Data()));
          }


      }
    }

  }

  // Print a table of the location of the wiggle
  if(printWiggleTable){

    // Find the indices of the bins included in the table
    int iCentrality0to10 = card[kPbPb]->FindBinIndexCentrality(0.0,10.0);
    int iCentrality10to30 = card[kPbPb]->FindBinIndexCentrality(10.0,30.0);
    int trackPtBinFor2GeV = card[kPbPb]->GetBinIndexTrackPtEEC(2.0);
    int trackPtBinFor3GeV = card[kPbPb]->GetBinIndexTrackPtEEC(3.0);
    int iJetPt120to140 = card[kPbPb]->FindBinIndexJetPtEEC(120.0,140.0);
    int iJetPt140to160 = card[kPbPb]->FindBinIndexJetPtEEC(140.0,160.0);
    int iJetPt160to180 = card[kPbPb]->FindBinIndexJetPtEEC(160.0,180.0);
    int iJetPt180to200 = card[kPbPb]->FindBinIndexJetPtEEC(180.0,200.0);

    for(int iUncertainty = kUncorrelatedUncertainty; iUncertainty <= kNUncertaintyTypes; iUncertainty++){

      // Write the beginning of the table
      cout << endl;
      cout << "Wiggle type: " << iUncertainty << endl;
      cout << endl;
      cout << "\\begin{table}[htbp]" << endl;
      cout << "  \\centering{" << endl;
      cout << "    \\topcaption{Location of the wiggle triangulated using three lines.}" << endl;
      cout << "    \\label{tab:wiggleTable}" << endl;
      cout << "    \\begin{tabular}{cccc}" << endl;
      cout << "      \\ptCh (GeV) & Jet \\pt (GeV) & Wiggle location 0--10\\%  & Wiggle location 10--30\\% \\\\" << endl;
      cout << "      \\hline" << endl;
      cout << Form("      \\multirow{4}*{$\\ptCh>2$} & $120<\\pt<140$ & $%.3f \\pm %.3f < \\Delta r < %.3f \\pm %.3f$ & $%.3f \\pm %.3f < \\Delta r < %.3f \\pm %.3f$ \\\\", firstTurningPoint[iUncertainty][iCentrality0to10][iJetPt120to140][trackPtBinFor2GeV], firstTurningPointError[iUncertainty][iCentrality0to10][iJetPt120to140][trackPtBinFor2GeV], secondTurningPoint[iUncertainty][iCentrality0to10][iJetPt120to140][trackPtBinFor2GeV], secondTurningPointError[iUncertainty][iCentrality0to10][iJetPt120to140][trackPtBinFor2GeV], firstTurningPoint[iUncertainty][iCentrality10to30][iJetPt120to140][trackPtBinFor2GeV], firstTurningPointError[iUncertainty][iCentrality10to30][iJetPt120to140][trackPtBinFor2GeV], secondTurningPoint[iUncertainty][iCentrality10to30][iJetPt120to140][trackPtBinFor2GeV], secondTurningPointError[iUncertainty][iCentrality10to30][iJetPt120to140][trackPtBinFor2GeV]) << endl;
      cout << Form("                               & $140<\\pt<160$ & $%.3f \\pm %.3f < \\Delta r < %.3f \\pm %.3f$ & $%.3f \\pm %.3f < \\Delta r < %.3f \\pm %.3f$ \\\\", firstTurningPoint[iUncertainty][iCentrality0to10][iJetPt140to160][trackPtBinFor2GeV], firstTurningPointError[iUncertainty][iCentrality0to10][iJetPt140to160][trackPtBinFor2GeV], secondTurningPoint[iUncertainty][iCentrality0to10][iJetPt140to160][trackPtBinFor2GeV], secondTurningPointError[iUncertainty][iCentrality0to10][iJetPt140to160][trackPtBinFor2GeV], firstTurningPoint[iUncertainty][iCentrality10to30][iJetPt140to160][trackPtBinFor2GeV], firstTurningPointError[iUncertainty][iCentrality10to30][iJetPt140to160][trackPtBinFor2GeV], secondTurningPoint[iUncertainty][iCentrality10to30][iJetPt140to160][trackPtBinFor2GeV], secondTurningPointError[iUncertainty][iCentrality10to30][iJetPt140to160][trackPtBinFor2GeV]) << endl;
      cout << Form("                               & $160<\\pt<180$ & $%.3f \\pm %.3f < \\Delta r < %.3f \\pm %.3f$ & $%.3f \\pm %.3f < \\Delta r < %.3f \\pm %.3f$ \\\\", firstTurningPoint[iUncertainty][iCentrality0to10][iJetPt160to180][trackPtBinFor2GeV], firstTurningPointError[iUncertainty][iCentrality0to10][iJetPt160to180][trackPtBinFor2GeV], secondTurningPoint[iUncertainty][iCentrality0to10][iJetPt160to180][trackPtBinFor2GeV], secondTurningPointError[iUncertainty][iCentrality0to10][iJetPt160to180][trackPtBinFor2GeV], firstTurningPoint[iUncertainty][iCentrality10to30][iJetPt160to180][trackPtBinFor2GeV], firstTurningPointError[iUncertainty][iCentrality10to30][iJetPt160to180][trackPtBinFor2GeV], secondTurningPoint[iUncertainty][iCentrality10to30][iJetPt160to180][trackPtBinFor2GeV], secondTurningPointError[iUncertainty][iCentrality10to30][iJetPt160to180][trackPtBinFor2GeV]) << endl;
      cout << Form("                               & $180<\\pt<200$ & $%.3f \\pm %.3f < \\Delta r < %.3f \\pm %.3f$ & $%.3f \\pm %.3f < \\Delta r < %.3f \\pm %.3f$ \\\\[\\cmsTabSkip]", firstTurningPoint[iUncertainty][iCentrality0to10][iJetPt180to200][trackPtBinFor2GeV], firstTurningPointError[iUncertainty][iCentrality0to10][iJetPt180to200][trackPtBinFor2GeV], secondTurningPoint[iUncertainty][iCentrality0to10][iJetPt180to200][trackPtBinFor2GeV], secondTurningPointError[iUncertainty][iCentrality0to10][iJetPt180to200][trackPtBinFor2GeV], firstTurningPoint[iUncertainty][iCentrality10to30][iJetPt180to200][trackPtBinFor2GeV], firstTurningPointError[iUncertainty][iCentrality10to30][iJetPt180to200][trackPtBinFor2GeV], secondTurningPoint[iUncertainty][iCentrality10to30][iJetPt180to200][trackPtBinFor2GeV], secondTurningPointError[iUncertainty][iCentrality10to30][iJetPt180to200][trackPtBinFor2GeV]) << endl;
      cout << Form("      \\multirow{4}*{$\\ptCh>3$} & $120<\\pt<140$ & $%.3f \\pm %.3f < \\Delta r < %.3f \\pm %.3f$ & $%.3f \\pm %.3f < \\Delta r < %.3f \\pm %.3f$ \\\\", firstTurningPoint[iUncertainty][iCentrality0to10][iJetPt120to140][trackPtBinFor3GeV], firstTurningPointError[iUncertainty][iCentrality0to10][iJetPt120to140][trackPtBinFor3GeV], secondTurningPoint[iUncertainty][iCentrality0to10][iJetPt120to140][trackPtBinFor3GeV], secondTurningPointError[iUncertainty][iCentrality0to10][iJetPt120to140][trackPtBinFor3GeV], firstTurningPoint[iUncertainty][iCentrality10to30][iJetPt120to140][trackPtBinFor3GeV], firstTurningPointError[iUncertainty][iCentrality10to30][iJetPt120to140][trackPtBinFor3GeV], secondTurningPoint[iUncertainty][iCentrality10to30][iJetPt120to140][trackPtBinFor3GeV], secondTurningPointError[iUncertainty][iCentrality10to30][iJetPt120to140][trackPtBinFor3GeV]) << endl;
      cout << Form("                               & $140<\\pt<160$ & $%.3f \\pm %.3f < \\Delta r < %.3f \\pm %.3f$ & $%.3f \\pm %.3f < \\Delta r < %.3f \\pm %.3f$ \\\\", firstTurningPoint[iUncertainty][iCentrality0to10][iJetPt140to160][trackPtBinFor3GeV], firstTurningPointError[iUncertainty][iCentrality0to10][iJetPt140to160][trackPtBinFor3GeV], secondTurningPoint[iUncertainty][iCentrality0to10][iJetPt140to160][trackPtBinFor3GeV], secondTurningPointError[iUncertainty][iCentrality0to10][iJetPt140to160][trackPtBinFor3GeV], firstTurningPoint[iUncertainty][iCentrality10to30][iJetPt140to160][trackPtBinFor3GeV], firstTurningPointError[iUncertainty][iCentrality10to30][iJetPt140to160][trackPtBinFor3GeV], secondTurningPoint[iUncertainty][iCentrality10to30][iJetPt140to160][trackPtBinFor3GeV], secondTurningPointError[iUncertainty][iCentrality10to30][iJetPt140to160][trackPtBinFor3GeV]) << endl;
      cout << Form("                               & $160<\\pt<180$ & $%.3f \\pm %.3f < \\Delta r < %.3f \\pm %.3f$ & $%.3f \\pm %.3f < \\Delta r < %.3f \\pm %.3f$ \\\\", firstTurningPoint[iUncertainty][iCentrality0to10][iJetPt160to180][trackPtBinFor3GeV], firstTurningPointError[iUncertainty][iCentrality0to10][iJetPt160to180][trackPtBinFor3GeV], secondTurningPoint[iUncertainty][iCentrality0to10][iJetPt160to180][trackPtBinFor3GeV], secondTurningPointError[iUncertainty][iCentrality0to10][iJetPt160to180][trackPtBinFor3GeV], firstTurningPoint[iUncertainty][iCentrality10to30][iJetPt160to180][trackPtBinFor3GeV], firstTurningPointError[iUncertainty][iCentrality10to30][iJetPt160to180][trackPtBinFor3GeV], secondTurningPoint[iUncertainty][iCentrality10to30][iJetPt160to180][trackPtBinFor3GeV], secondTurningPointError[iUncertainty][iCentrality10to30][iJetPt160to180][trackPtBinFor3GeV]) << endl;
      cout << Form("                               & $180<\\pt<200$ & $%.3f \\pm %.3f < \\Delta r < %.3f \\pm %.3f$ & $%.3f \\pm %.3f < \\Delta r < %.3f \\pm %.3f$ \\\\", firstTurningPoint[iUncertainty][iCentrality0to10][iJetPt180to200][trackPtBinFor3GeV], firstTurningPointError[iUncertainty][iCentrality0to10][iJetPt180to200][trackPtBinFor3GeV], secondTurningPoint[iUncertainty][iCentrality0to10][iJetPt180to200][trackPtBinFor3GeV], secondTurningPointError[iUncertainty][iCentrality0to10][iJetPt180to200][trackPtBinFor3GeV], firstTurningPoint[iUncertainty][iCentrality10to30][iJetPt180to200][trackPtBinFor3GeV], firstTurningPointError[iUncertainty][iCentrality10to30][iJetPt180to200][trackPtBinFor3GeV], secondTurningPoint[iUncertainty][iCentrality10to30][iJetPt180to200][trackPtBinFor3GeV], secondTurningPointError[iUncertainty][iCentrality10to30][iJetPt180to200][trackPtBinFor3GeV]) << endl;
      cout << "    \\end{tabular}" << endl;
      cout << "  }" << endl;
      cout << "\\end{table}" << endl;

    } // Uncertainty loop

  } // Print a table with infomation about the wiggle
}
