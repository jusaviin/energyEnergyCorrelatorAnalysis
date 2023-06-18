#include "SystematicUncertaintyOrganizer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECHistogramManager.h"
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Logarithmic function used to fit the ratio plot
 */
double logFit(double *x, double *par){
  return par[1] * TMath::Log10(x[0]) + par[0];
}

/*
 * Macro for making final result plots comparing energy-energy correlators between pp and PbPb
 */
void fitPbPbToPpRatio(){

  enum enumDataType{kPbPb, kPp, kNDataTypes};

  // ============= //
  // Configuration //
  // ============= //
  
  // Input files
  TString inputFileName[kNDataTypes];
  inputFileName[kPbPb] = "data/eecAnalysis_akFlowJet_wtaAxis_binningForUnfolding_processed_nominalResult_2023-05-23.root";
  inputFileName[kPp] = "data/ppData_pfJets_wtaAxis_nominalResults_processed_2023-06-12.root";
  TString uncertaintyFileName[kNDataTypes];
  uncertaintyFileName[kPbPb] = "systematicUncertainties/systematicUncertainties_firstLook_2023-06-13.root";
  uncertaintyFileName[kPp] = "systematicUncertainties/systematicUncertaintiesForPp_firstLook_2023-06-13.root";
  
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
  int lastDrawnCentralityBin = card[kPbPb]->GetLastUnfoldedCentralityBin();
  
  int firstDrawnJetPtBinEEC = card[kPbPb]->GetFirstUnfoldedJetPtBin();
  int lastDrawnJetPtBinEEC = card[kPbPb]->GetLastUnfoldedJetPtBin();
  
  int firstDrawnTrackPtBinEEC = card[kPbPb]->GetFirstUnfoldedTrackPtBin();
  int lastDrawnTrackPtBinEEC = card[kPbPb]->GetLastUnfoldedTrackPtBin();

  // Choose which plots to draw
  bool drawRatiosWithFits = true;
  bool drawWiggleGraph = true;
  
  // Save the final plots
  const bool saveFigures = true;
  TString saveComment = "_firstLook";

  // Ratio zoom settings
  std::pair<double, double> analysisDeltaR = std::make_pair(0.006, 0.39); // DeltaR span in which the analysis is done
  std::pair<double, double> ratioZoom = std::make_pair(0.4, 1.6);
  
  // Marker colors and styles
  int markerStylePbPb[] = {kFullSquare, kFullCircle, kFullCross, kFullFourTrianglesPlus};
  int markerStylePp = kFullDiamond;
  int markerColorPbPb[] = {kRed, kBlue, kMagenta, kGreen+3};
  int markerColorPp = kBlack;

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
 
  TH1D* energyEnergyCorrelatorSignalPbPb[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorSignalPp[nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorPbPbToPpRatio[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyForPbPb[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyForPp[nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyPbPbToPpRatio[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];

  const int nFitParts = 3;
  TF1* fitToRatio[nFitParts][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];

  // Ranges for fits, should be optimized bin by bin
  int fitColor[4][nFitParts] = {{kBlack, kBlue, kBlack},  // 0-10% centrality
                                {kRed, kMagenta, kRed},  // 10-30% centrality
                                {kBlack, kBlue, kBlack},  // 30-50% centrality
                                {kRed, kMagenta, kRed}}; // 50-90% centrality

  double fitLowRange[4][3][4][nFitParts] =
  // Centrality 0-10 %
  {
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

  // Centrality 10-30 %
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

  // Centrality 10-30 %
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
      systematicUncertaintyForPp[iJetPt][iTrackPt] = NULL;
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt] = NULL;
        energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt] = NULL;
        systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt] = NULL;
        systematicUncertaintyPbPbToPpRatio[iCentrality][iJetPt][iTrackPt] = NULL;

        // Initialize the fit functions to NULL
        for(int iFit = 0; iFit < nFitParts; iFit++){
          fitToRatio[iFit][iCentrality][iJetPt][iTrackPt] = NULL;
        }
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
      systematicUncertaintyForPp[iJetPt][iTrackPt] = uncertainties[kPp]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty);

      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        iCentralityMatched = uncertaintyCard[kPbPb]->FindBinIndexCentrality(card[kPbPb]->GetBinBordersCentrality(iCentrality));

        // Read the PbPb histograms
        energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt] = histograms[kPbPb]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);
        systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt] = uncertainties[kPbPb]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty);

        // Initialize the fit functions with the defined function
        for(int iFit = 0; iFit < nFitParts; iFit++){
          fitToRatio[iFit][iCentrality][iJetPt][iTrackPt] = new TF1(Form("fit%d%d%d%d", iCentrality, iJetPt, iTrackPt, iFit), logFit, fitLowRange[iCentrality - firstDrawnCentralityBin][iTrackPt - firstDrawnTrackPtBinEEC][iJetPt - firstDrawnJetPtBinEEC][iFit], fitHighRange[iCentrality - firstDrawnCentralityBin][iTrackPt - firstDrawnTrackPtBinEEC][iJetPt - firstDrawnJetPtBinEEC][iFit], 2);
          fitToRatio[iFit][iCentrality][iJetPt][iTrackPt]->SetParameters(-1, 1);
          fitToRatio[iFit][iCentrality][iJetPt][iTrackPt]->SetLineColor(fitColor[iCentrality][iFit]);
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
      systematicUncertaintyForPp[iJetPt][iTrackPt]->Scale(1.0 / systematicUncertaintyForPp[iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->Scale(1.0 / energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
        systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->Scale(1.0 / systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

        energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt] = (TH1D*) energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->Clone(Form("energyEnergyCorrelatorRatio%d%d%d", iCentrality, iJetPt, iTrackPt));
        energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->Divide(energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]);

        systematicUncertaintyPbPbToPpRatio[iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->Clone(Form("systematicUncertaintyRatio%d%d%d", iCentrality, iJetPt, iTrackPt));
        systematicUncertaintyPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->Divide(systematicUncertaintyForPp[iJetPt][iTrackPt]);

        // TODO: For fitting, combine the systematic and statistical uncertainties

        // Fit the thingy
        for(int iFit = 0; iFit < nFitParts; iFit++){
          systematicUncertaintyPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->Fit(fitToRatio[iFit][iCentrality][iJetPt][iTrackPt], "NQ", "", fitLowRange[iCentrality - firstDrawnCentralityBin][iTrackPt - firstDrawnTrackPtBinEEC][iJetPt - firstDrawnJetPtBinEEC][iFit], fitHighRange[iCentrality - firstDrawnCentralityBin][iTrackPt - firstDrawnTrackPtBinEEC][iJetPt - firstDrawnJetPtBinEEC][iFit]);
        }

      } // Centrality loop
    } // Track pT loop
  } // Jet pT loop

  // Find the deltaR value corresponding to the points where the different fit region overlap
  double firstTurningPoint[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  double secondTurningPoint[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        firstTurningPoint[iCentrality][iJetPt][iTrackPt] = 0;
        secondTurningPoint[iCentrality][iJetPt][iTrackPt] = 0;
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

  double a,b,c,d;
  for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
    for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
      for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
        
        // Calculate the x-value for which the two functions overlap
        // a*log(x)+b = c*log(x)+d => log(x) = (d - b) / (a - c) => x = 10^[(d - b) / (a - c)]
        a = fitToRatio[0][iCentrality][iJetPt][iTrackPt]->GetParameter(1);
        b = fitToRatio[0][iCentrality][iJetPt][iTrackPt]->GetParameter(0);
        c = fitToRatio[1][iCentrality][iJetPt][iTrackPt]->GetParameter(1);
        d = fitToRatio[1][iCentrality][iJetPt][iTrackPt]->GetParameter(0);
        firstTurningPoint[iCentrality][iJetPt][iTrackPt] = TMath::Power(10.0,((d-b)/(a-c)));

        a = fitToRatio[2][iCentrality][iJetPt][iTrackPt]->GetParameter(1);
        b = fitToRatio[2][iCentrality][iJetPt][iTrackPt]->GetParameter(0);
        secondTurningPoint[iCentrality][iJetPt][iTrackPt] = TMath::Power(10.0,((d-b)/(a-c)));
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

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
        xPointsWiggleGraph[iCentrality][iJetPt][iTrackPt][0] = firstTurningPoint[iCentrality][iJetPt][iTrackPt];
        xPointsWiggleGraph[iCentrality][iJetPt][iTrackPt][1] = secondTurningPoint[iCentrality][iJetPt][iTrackPt];
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

  // Draw individual plots with all centralities mixed together in a single figure
  if(drawRatiosWithFits){

    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      centralityString = Form("Cent: %.0f-%.0f%%", card[kPbPb]->GetLowBinBorderCentrality(iCentrality), card[kPbPb]->GetHighBinBorderCentrality(iCentrality));
      compactCentralityString = Form("_C=%.0f-%.0f", card[kPbPb]->GetLowBinBorderCentrality(iCentrality), card[kPbPb]->GetHighBinBorderCentrality(iCentrality));

      for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
        jetPtString = Form("%.0f < jet p_{T} < %.0f", card[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));
        compactJetPtString = Form("_J=%.0f-%.0f", card[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));

        for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
          trackPtString = Form("%.1f < track p_{T}", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));
          compactTrackPtString = Form("_T>%.1f", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));
          compactTrackPtString.ReplaceAll(".", "v");

          // Setup the legend for plots
          legend = new TLegend(0.53, 0.6, 0.83, 0.85);
          legend->SetFillStyle(0);
          legend->SetBorderSize(0);
          legend->SetTextSize(0.05);
          legend->SetTextFont(62);
          legend->AddEntry((TObject*)0, centralityString.Data(), "");
          legend->AddEntry((TObject*)0, jetPtString.Data(), "");
          legend->AddEntry((TObject*)0, trackPtString.Data(), "");

          // Linear scale for the ratio
          drawer->SetLogY(false);

          // Set the axis drawing ranges
          systematicUncertaintyPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

          // Set the style for histograms
          for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++) {
            systematicUncertaintyPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin]);
            systematicUncertaintyPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);
            energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin]);
            energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);

            systematicUncertaintyPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin], 0.4);
            systematicUncertaintyPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
            energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
            systematicUncertaintyPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
            energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
          }

          drawer->SetGridY(true);
          drawer->DrawHistogram(systematicUncertaintyPbPbToPpRatio[iCentrality][iJetPt][iTrackPt], "#Deltar", "#frac{PbPb}{pp}", " ", "e2");
          energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->Draw("same,p");

          // Draw the fits
          for(int iFit = 0; iFit < nFitParts; iFit++){
            fitToRatio[iFit][iCentrality][iJetPt][iTrackPt]->Draw("same");
          }
          drawer->SetGridY(false);

          // Draw the legend
          legend->Draw();

          // If a plot name is given, save the plot in a file
          if(saveFigures) {
            gPad->GetCanvas()->SaveAs(Form("figures/energyEnergyCorrelator_ratioFit%s%s%s%s.pdf", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data()));
          }

        }  // Track pT loop
      }    // Jet pT loop
    }      // Centrality loop
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

}
