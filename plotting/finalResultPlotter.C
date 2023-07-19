#include "SystematicUncertaintyOrganizer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECHistogramManager.h"
#include "EECCard.h"
#include "JDrawer.h"
#include "SplitCanvas.h"

/*
 * Determine up and down uncertainty bands for shape correlated uncertainties
 */
void calculateCorrelatedBands(TH1D* correlatedHistogram, TH1D* correlatedBandUpHistogram, TH1D* correlatedBandDownHistogram){
  double binContent, binError;
  int criticalBin = correlatedHistogram->GetXaxis()->FindBin(0.055);
  for(int iBin = 1; iBin <= criticalBin; iBin++){
    binContent = correlatedHistogram->GetBinContent(iBin);
    binError = correlatedHistogram->GetBinError(iBin);
    correlatedBandUpHistogram->SetBinContent(iBin, binContent+binError/2.0);
    correlatedBandUpHistogram->SetBinError(iBin, binError/2.0);
    correlatedBandDownHistogram->SetBinContent(iBin, binContent-binError/2.0);
    correlatedBandDownHistogram->SetBinError(iBin, binError/2.0);
  }
  for(int iBin = criticalBin+1; iBin <= correlatedHistogram->GetNbinsX(); iBin++){
    binContent = correlatedHistogram->GetBinContent(iBin);
    binError = correlatedHistogram->GetBinError(iBin);
    correlatedBandUpHistogram->SetBinContent(iBin, binContent-binError/2.0);
    correlatedBandUpHistogram->SetBinError(iBin, binError/2.0);
    correlatedBandDownHistogram->SetBinContent(iBin, binContent+binError/2.0);
    correlatedBandDownHistogram->SetBinError(iBin, binError/2.0);
  }
}

/*
 * Macro for making final result plots comparing energy-energy correlators between pp and PbPb
 */
void finalResultPlotter(){

  enum enumDataType{kPbPb, kPp, kNDataTypes};
  enum enumSystematicUncertaintyType{kUncorrelatedUncertainty, kCorrelatedUncertainty, kCorrelatedUncertaintyShapeUp, kCorrelatedUncertaintyShapeDown, knSystematicUncertaintyTypes};

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
  int lastDrawnCentralityBin = card[kPbPb]->GetLastUnfoldedCentralityBin();
  
  int firstDrawnJetPtBinEEC = card[kPbPb]->GetFirstUnfoldedJetPtBin();
  int lastDrawnJetPtBinEEC = card[kPbPb]->GetLastUnfoldedJetPtBin();
  
  int firstDrawnTrackPtBinEEC = card[kPbPb]->GetFirstUnfoldedTrackPtBin();
  int lastDrawnTrackPtBinEEC = card[kPbPb]->GetLastUnfoldedTrackPtBin();

  // Choose which plots to draw
  bool drawIndividualPlotsAllCentralities = false;
  bool drawBigCanvasDistributions = false;
  bool drawBigCanvasRatios = false;
  bool drawDoubleRatios = false;
  bool drawDoubleRatioToSingleCanvas = true;

  // Select the bins to be drawn for double ratio plots
  std::pair<double, double> doubleRatioCentralityBin1 = std::make_pair(0.0,10.0);
  std::pair<double, double> doubleRatioCentralityBin2 = std::make_pair(10.0,30.0);
  std::pair<double, double> doubleRatioJetPtBin = std::make_pair(160,180);
  int doubleRatioCentralityBinIndex1;
  int doubleRatioCentralityBinIndex2;
  int doubleRatioJetPtBinIndex;
  
  // Save the final plots
  const bool saveFigures = true;
  TString saveComment = "_uncertaintyUpdate";

  // Ratio zoom settings
  std::pair<double, double> ratioZoom = std::make_pair(0.4, 1.6);
  std::pair<double, double> analysisDeltaR = std::make_pair(0.006, 0.39); // DeltaR span in which the analysis is done
  
  // Marker colors and styles
  int markerStylePbPb[] = {kFullSquare, kFullCircle, kFullCross, kFullFourTrianglesPlus};
  int markerStylePp = kFullDiamond;
  int markerColorPbPb[] = {kRed, kBlue, kMagenta, kGreen+3};
  int markerColorPp = kBlack;
  int bandColorUpPbPb[] = {kOrange+7, kViolet-3, kPink-3, kOrange-3};
  int bandColorDownPbPb[] = {kPink+9, kAzure+8, kViolet+6, kSpring};

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
 
  // Energy-energy correlators and PbPb to pp ratios
  TH1D* energyEnergyCorrelatorSignalPbPb[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorSignalPp[nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorPbPbToPpRatio[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyForPbPb[knSystematicUncertaintyTypes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyForPp[knSystematicUncertaintyTypes][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyPbPbToPpRatio[knSystematicUncertaintyTypes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];

  // Double ratios from energy-energy correlator histograms
  TH1D* energyEnergyCorrelatorForDoubleRatioFromPbPb[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorForDoubleRatioFromPp[nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorDoubleRatio[nCentralityBins][nJetPtBinsEEC];
  TH1D* systematicUncertaintyForDoubleRatioFromPbPb[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyForDoubleRatioFromPp[nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyDoubleRatio[nCentralityBins][nJetPtBinsEEC];

  // Initialize histograms to NULL
  for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      energyEnergyCorrelatorDoubleRatio[iCentrality][iJetPt] = NULL;
      systematicUncertaintyDoubleRatio[iCentrality][iJetPt] = NULL;
    }
    for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
      energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt] = NULL;

      energyEnergyCorrelatorForDoubleRatioFromPp[iJetPt][iTrackPt] = NULL;
      systematicUncertaintyForDoubleRatioFromPp[iJetPt][iTrackPt] = NULL;

      for(int iUncertainty = 0; iUncertainty < knSystematicUncertaintyTypes; iUncertainty++){
        systematicUncertaintyForPp[iUncertainty][iJetPt][iTrackPt] = NULL;
      } // Uncertainty type loop

      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt] = NULL;
        energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt] = NULL;

        energyEnergyCorrelatorForDoubleRatioFromPbPb[iCentrality][iJetPt][iTrackPt] = NULL;
        systematicUncertaintyForDoubleRatioFromPbPb[iCentrality][iJetPt][iTrackPt] = NULL;

        for(int iUncertainty = 0; iUncertainty < knSystematicUncertaintyTypes; iUncertainty++){
          systematicUncertaintyForPbPb[iUncertainty][iCentrality][iJetPt][iTrackPt] = NULL;
          systematicUncertaintyPbPbToPpRatio[iUncertainty][iCentrality][iJetPt][iTrackPt] = NULL;
        }
      } // Centrality loop
    } // Track pT loop
  } // Jet pT loop

  // Read the histograms from managers
  int iTrackPtMatchedPp, iTrackPtMatchedPbPbUncertainty, iTrackPtMatchedPpUncertainty;
  int iJetPtMatchedPp, iJetPtMatchedPbPbUncertainty, iJetPtMatchedPpUncertainty;
  int iCentralityMatched;

  // Define helper histograms to determine the uncertainties relevant for the double ratio
  TH1D* uncertaintyTrackSelection;
  TH1D* uncertaintyBackgroundSubtraction;
  TH1D* uncertaintyTrackPairEfficiency;
  TH1D* uncertaintyMCnonclosure;
  double trackCorrelation = 0.2;  // TODO: Check what is the actual number of overlapping tracks

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
      systematicUncertaintyForPp[kUncorrelatedUncertainty][iJetPt][iTrackPt] = uncertainties[kPp]->GetUncorrelatedSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty);
      systematicUncertaintyForPp[kCorrelatedUncertainty][iJetPt][iTrackPt] = uncertainties[kPp]->GetCorrelatedSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty);

      // Read the relevant systematic uncertainties for double ratio from pp
      uncertaintyBackgroundSubtraction = uncertainties[kPp]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kBackgroundSubtraction);
      uncertaintyTrackPairEfficiency = uncertainties[kPp]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kTrackPairEfficiency);
      uncertaintyMCnonclosure = uncertainties[kPp]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kMonteCarloNonClosure);

      systematicUncertaintyForDoubleRatioFromPp[iJetPt][iTrackPt] = (TH1D*) uncertaintyBackgroundSubtraction->Clone(Form("doubleRatioUncertaintyFromPp%d%d", iJetPt, iTrackPt));
      for(int iBin = 1; iBin <= uncertaintyBackgroundSubtraction->GetNbinsX(); iBin++){
        systematicUncertaintyForDoubleRatioFromPp[iJetPt][iTrackPt]->SetBinError(iBin, TMath::Sqrt(TMath::Power(uncertaintyBackgroundSubtraction->GetBinError(iBin),2) + TMath::Power(uncertaintyTrackPairEfficiency->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintyMCnonclosure->GetBinError(iBin),2)));
      }

      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        iCentralityMatched = uncertaintyCard[kPbPb]->FindBinIndexCentrality(card[kPbPb]->GetBinBordersCentrality(iCentrality));

        // Read the PbPb histograms
        energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt] = histograms[kPbPb]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);
        systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt] = uncertainties[kPbPb]->GetUncorrelatedSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty);
        systematicUncertaintyForPbPb[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt] = uncertainties[kPbPb]->GetCorrelatedSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty);

        // Read the relevant systematic uncertainties for double ratio from PbPb
        uncertaintyBackgroundSubtraction = uncertainties[kPbPb]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kBackgroundSubtraction);
        uncertaintyTrackPairEfficiency = uncertainties[kPbPb]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kTrackPairEfficiency);
        uncertaintyTrackSelection = uncertainties[kPbPb]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kTrackSelection);
        uncertaintyMCnonclosure = uncertainties[kPbPb]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kMonteCarloNonClosure);

        systematicUncertaintyForDoubleRatioFromPbPb[iCentrality][iJetPt][iTrackPt] = (TH1D*) uncertaintyBackgroundSubtraction->Clone(Form("doubleRatioUncertaintyFromPbPb%d%d%d", iCentrality, iJetPt, iTrackPt));
        for(int iBin = 1; iBin <= uncertaintyBackgroundSubtraction->GetNbinsX(); iBin++){
          systematicUncertaintyForDoubleRatioFromPbPb[iCentrality][iJetPt][iTrackPt]->SetBinError(iBin, TMath::Sqrt(TMath::Power(uncertaintyBackgroundSubtraction->GetBinError(iBin),2) + TMath::Power(uncertaintyTrackPairEfficiency->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintyTrackSelection->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintyMCnonclosure->GetBinError(iBin),2)));
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
      systematicUncertaintyForDoubleRatioFromPp[iJetPt][iTrackPt]->Scale(1.0 / systematicUncertaintyForDoubleRatioFromPp[iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->Scale(1.0 / energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
        systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Scale(1.0 / systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
        systematicUncertaintyForPbPb[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Scale(1.0 / systematicUncertaintyForPbPb[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
        systematicUncertaintyForDoubleRatioFromPbPb[iCentrality][iJetPt][iTrackPt]->Scale(1.0 / systematicUncertaintyForDoubleRatioFromPbPb[iCentrality][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

        energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt] = (TH1D*) energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->Clone(Form("energyEnergyCorrelatorRatio%d%d%d", iCentrality, iJetPt, iTrackPt));
        energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->Divide(energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]);

        systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("systematicUncertaintyRatioUncorrelated%d%d%d", iCentrality, iJetPt, iTrackPt));
        systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Divide(systematicUncertaintyForPp[kUncorrelatedUncertainty][iJetPt][iTrackPt]);

        systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPbPb[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("systematicUncertaintyRatioCorrelated%d%d%d", iCentrality, iJetPt, iTrackPt));
        systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Divide(systematicUncertaintyForPp[kCorrelatedUncertainty][iJetPt][iTrackPt]);
      } // Centrality loop
    } // Track pT loop
  } // Jet pT loop

  // Reduce the statistical uncertainties by the overlapping statistics fraction for the double ratios
  for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
    for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
      energyEnergyCorrelatorForDoubleRatioFromPp[iJetPt][iTrackPt] = (TH1D*) energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->Clone(Form("energyEnergyCorrelatorForDoubleRatioFromPp%d%d", iJetPt, iTrackPt));

      for(int iBin = 1; iBin <= energyEnergyCorrelatorForDoubleRatioFromPp[iJetPt][iTrackPt]->GetNbinsX(); iBin++){
        energyEnergyCorrelatorForDoubleRatioFromPp[iJetPt][iTrackPt]->SetBinError(iBin, energyEnergyCorrelatorForDoubleRatioFromPp[iJetPt][iTrackPt]->GetBinError(iBin) * trackCorrelation);
      }

      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        energyEnergyCorrelatorForDoubleRatioFromPbPb[iCentrality][iJetPt][iTrackPt] = (TH1D*) energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->Clone(Form("energyEnergyCorrelatorForDoubleRatioFromPbPb%d%d%d", iCentrality, iJetPt, iTrackPt));

      for(int iBin = 1; iBin <= energyEnergyCorrelatorForDoubleRatioFromPbPb[iCentrality][iJetPt][iTrackPt]->GetNbinsX(); iBin++){
        energyEnergyCorrelatorForDoubleRatioFromPbPb[iCentrality][iJetPt][iTrackPt]->SetBinError(iBin, energyEnergyCorrelatorForDoubleRatioFromPbPb[iCentrality][iJetPt][iTrackPt]->GetBinError(iBin) * trackCorrelation);
      }

      } // Centrality loop
    } // Track pT loop
  } // Jet pT loop

  // After regular ratios have been calculated, proceed to calculating double ratio. Here we need to use different histograms as above to properly take into account systematic uncertainty cancellation due to correlated datasets.
  int trackPtBinFor2GeV = card[kPbPb]->GetBinIndexTrackPtEEC(2.0);
  int trackPtBinFor3GeV = card[kPbPb]->GetBinIndexTrackPtEEC(3.0);
  for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
    for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
      for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
        // Calculate the single ratios with properly handled double ratio uncertainties
        energyEnergyCorrelatorForDoubleRatioFromPbPb[iCentrality][iJetPt][iTrackPt]->Divide(energyEnergyCorrelatorForDoubleRatioFromPp[iJetPt][iTrackPt]);
        systematicUncertaintyForDoubleRatioFromPbPb[iCentrality][iJetPt][iTrackPt]->Divide(systematicUncertaintyForDoubleRatioFromPp[iJetPt][iTrackPt]);
      }
      // Calculate the double ratios from the single ratios with properly handled uncertainties
      energyEnergyCorrelatorDoubleRatio[iCentrality][iJetPt] = (TH1D*) energyEnergyCorrelatorForDoubleRatioFromPbPb[iCentrality][iJetPt][trackPtBinFor3GeV]->Clone(Form("energyEnergyCorrelatorDoubleRatio%d%d", iCentrality, iJetPt));
      energyEnergyCorrelatorDoubleRatio[iCentrality][iJetPt]->Divide(energyEnergyCorrelatorForDoubleRatioFromPbPb[iCentrality][iJetPt][trackPtBinFor2GeV]);
      systematicUncertaintyDoubleRatio[iCentrality][iJetPt] = (TH1D*) systematicUncertaintyForDoubleRatioFromPbPb[iCentrality][iJetPt][trackPtBinFor3GeV]->Clone(Form("systematicUncertaintyDoubleRatio%d%d", iCentrality, iJetPt));
      systematicUncertaintyDoubleRatio[iCentrality][iJetPt]->Divide(systematicUncertaintyForDoubleRatioFromPbPb[iCentrality][iJetPt][trackPtBinFor2GeV]);
    }
  }

  // For illustration purposes, create up and down shifted uncertainty bands
  for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
    for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){

      systematicUncertaintyForPp[kCorrelatedUncertaintyShapeUp][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPp[kCorrelatedUncertainty][iJetPt][iTrackPt]->Clone(Form("correlatedUpBandPp%d%d", iJetPt, iTrackPt));
      systematicUncertaintyForPp[kCorrelatedUncertaintyShapeDown][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPp[kCorrelatedUncertainty][iJetPt][iTrackPt]->Clone(Form("correlatedDownBandPp%d%d", iJetPt, iTrackPt));
      calculateCorrelatedBands(systematicUncertaintyForPp[kCorrelatedUncertainty][iJetPt][iTrackPt], systematicUncertaintyForPp[kCorrelatedUncertaintyShapeUp][iJetPt][iTrackPt], systematicUncertaintyForPp[kCorrelatedUncertaintyShapeDown][iJetPt][iTrackPt]);

      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){

        systematicUncertaintyForPbPb[kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPbPb[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("correlatedUpBand%d%d%d", iCentrality, iJetPt, iTrackPt));
        systematicUncertaintyForPbPb[kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPbPb[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("correlatedDownBand%d%d%d", iCentrality, iJetPt, iTrackPt));
        calculateCorrelatedBands(systematicUncertaintyForPbPb[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], systematicUncertaintyForPbPb[kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt], systematicUncertaintyForPbPb[kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]);

        systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("correlatedUpBandForRatio%d%d%d", iCentrality, iJetPt, iTrackPt));
        systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("correlatedDownBandForRatio%d%d%d", iCentrality, iJetPt, iTrackPt));
        calculateCorrelatedBands(systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt], systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]);

      } // Centrality loop
    } // Track pT loop
  } // Jet pT loop

  // After all the ratios are calculated, calculate the double ratios 

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
  TLine* oneLine = new TLine(analysisDeltaR.first, 1, analysisDeltaR.second, 1);
  TBox* box;
  oneLine->SetLineColor(kBlack);
  oneLine->SetLineStyle(2);
  int canvasIndex;
  double bottomRowScale, bottomPadMargin, leftPadMargin;
  double leftMarginAdder, bottomMarginAdder;
  double thisPadScale;

  // Draw individual plots with all centralities mixed together in a single figure
  if(drawIndividualPlotsAllCentralities){

    for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
      jetPtString = Form("%.0f < jet p_{T} < %.0f", card[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));
      compactJetPtString = Form("_J=%.0f-%.0f", card[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));

      for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
        trackPtString = Form("%.1f < track p_{T}", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));
        compactTrackPtString = Form("_T>%.1f", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));
        compactTrackPtString.ReplaceAll(".","v");

        // Create a new canvas for the plot
        drawer->CreateSplitCanvas();

        // Use logarithmic axis for EEC
        drawer->SetLogY(true);

        // Setup the legend for plots
        legend = new TLegend(0.23, 0.05, 0.53, 0.6);
        legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);
        legend->AddEntry((TObject*)0, jetPtString.Data(), "");
        legend->AddEntry((TObject*)0, trackPtString.Data(), "");

        // Set the drawing style for pp histogram
        energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->SetMarkerStyle(kFullDiamond);
        energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->SetMarkerColor(kBlack);

        // Set the x-axis drawing range
        energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);

        // Draw the pp correlator to upper canves
        drawer->DrawHistogramToUpperPad(energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt], "#Deltar", "EEC Signal", " ", "p");
        legend->AddEntry(energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt], "pp", "p");

        // Draw the different centrality bins to the same plot
        for(int iCentrality = lastDrawnCentralityBin; iCentrality >= firstDrawnCentralityBin; iCentrality--) {
          energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
          energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin]);
          energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->Draw("same,p");
          legend->AddEntry(energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt], Form("PbPb %.0f-%.0f%%", card[kPbPb]->GetLowBinBorderCentrality(iCentrality), card[kPbPb]->GetHighBinBorderCentrality(iCentrality)), "p");
        }

        // Draw the legend to the upper pad
        legend->Draw();

        // Linear scale for the ratio
        drawer->SetLogY(false);

        // Set the axis drawing ranges
        energyEnergyCorrelatorPbPbToPpRatio[lastDrawnCentralityBin][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
        energyEnergyCorrelatorPbPbToPpRatio[lastDrawnCentralityBin][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

        // Set the style for histograms
        for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++) {
          energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
          energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin]);
        }

        drawer->SetGridY(true);
        drawer->DrawHistogramToLowerPad(energyEnergyCorrelatorPbPbToPpRatio[lastDrawnCentralityBin][iJetPt][iTrackPt], "#Deltar", "#frac{PbPb}{pp}", " ");
        for(int iCentrality = lastDrawnCentralityBin - 1; iCentrality >= firstDrawnCentralityBin; iCentrality--) {
          energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->Draw("same,p");
        }
        drawer->SetGridY(false);

        // If a plot name is given, save the plot in a file
        if(saveFigures) {
          gPad->GetCanvas()->SaveAs(Form("figures/energyEnergyCorrelator_centralityComparison%s%s%s.pdf", saveComment.Data(), compactJetPtString.Data(), compactTrackPtString.Data()));
        }

      } // Track pT loop
    } // Jet pT loop
  } // Drawing individual canvases with all centralities in one canvas

  // ============================================= //
  // ===== Drawing style with one big canvas ===== //
  // ============================================= //
  
  if(drawBigCanvasDistributions){
    
    // Draw all the distributions to big canvases
    SplitCanvas* bigCanvas[nTrackPtBinsEEC];
    mainTitle = new TLatex();

    for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
      compactTrackPtString = Form("_T>%.1f", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));
      compactTrackPtString.ReplaceAll(".","v");

      // Draw a big canvas and put all the plots in it
      bigCanvas[iTrackPt] = new SplitCanvas(Form("bigCanvas%d", iTrackPt), "", 1800, 1500);
      bigCanvas[iTrackPt]->SetMargin(0.08, 0.01, 0.08, 0.11);
      bigCanvas[iTrackPt]->DivideNeatly(4, 4);

      bottomRowScale = bigCanvas[iTrackPt]->GetBottomRowScale();
      bottomPadMargin = bigCanvas[iTrackPt]->GetBottomPadMargin();
      leftPadMargin = bigCanvas[iTrackPt]->GetLeftPadMargin();

      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++) {
          jetPtString = Form("%.0f < jet p_{T} < %.0f", card[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));

          canvasIndex = (iCentrality - firstDrawnCentralityBin) * (lastDrawnJetPtBinEEC - firstDrawnJetPtBinEEC + 1) + (iJetPt - firstDrawnJetPtBinEEC);
          bigCanvas[iTrackPt]->CD(canvasIndex);

          gPad->SetLogy();
          gPad->SetLogx();

          // The titles in the bottow row need to be scaled, since it has more margin than other rows
          thisPadScale = (iCentrality == lastDrawnCentralityBin) ? bottomRowScale : 1;

          // Set the axis titles and labels
          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleOffset(0.85);
          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleSize(0.2 * thisPadScale);
          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelOffset(0.01);
          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelSize(0.12 * thisPadScale);

          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleOffset(0.5 / thisPadScale);
          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleSize(0.2 * thisPadScale);
          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelOffset(0.01 / thisPadScale);
          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelSize(0.12 * thisPadScale);

          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitle("EEC");
          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetStats(0);
          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetNdivisions(505);
          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->CenterTitle();
          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitle("#Deltar");
          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetNdivisions(505);

          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetTitle("");

          // Set the drawing ranges
          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.008, 0.39);
          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.15, 30);
          systematicUncertaintyForPbPb[kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);
          systematicUncertaintyForPbPb[kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);
          systematicUncertaintyForPp[kCorrelatedUncertainty][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);

          // Set the drawing style for histograms
          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin]);
          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);
          energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin]);
          energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);

          systematicUncertaintyForPp[kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePp);
          systematicUncertaintyForPp[kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerSize(1.2);
          energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->SetMarkerStyle(markerStylePp);
          energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->SetMarkerSize(1.2);

          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin], 0.4);
          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
          energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
          energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);

          systematicUncertaintyForPp[kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPp, 0.4);
          systematicUncertaintyForPp[kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetLineColor(markerColorPp);
          energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->SetLineColor(markerColorPp);
          systematicUncertaintyForPp[kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerColor(markerColorPp);
          energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->SetMarkerColor(markerColorPp);

          // Different drawing styles for correlated uncertainties
          systematicUncertaintyForPp[kCorrelatedUncertainty][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPp, 0.4);
          systematicUncertaintyForPbPb[kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorUpPbPb[iCentrality - firstDrawnCentralityBin], 0.4);
          systematicUncertaintyForPbPb[kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorDownPbPb[iCentrality - firstDrawnCentralityBin], 0.4);
          systematicUncertaintyForPbPb[kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
          systematicUncertaintyForPbPb[kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);
          systematicUncertaintyForPbPb[kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
          systematicUncertaintyForPbPb[kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);

          // Draw the histograms
          systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Draw("e2");
          systematicUncertaintyForPp[kUncorrelatedUncertainty][iJetPt][iTrackPt]->Draw("same,e2");
          systematicUncertaintyForPp[kCorrelatedUncertainty][iJetPt][iTrackPt]->Draw("same,e3");
          systematicUncertaintyForPbPb[kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          systematicUncertaintyForPbPb[kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->Draw("same,p");
          energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->Draw("same,p");

          // Create a legend for each pad
          leftMarginAdder = (iJetPt == firstDrawnJetPtBinEEC) ? leftPadMargin : 0;
          bottomMarginAdder = (iCentrality == lastDrawnCentralityBin) ? bottomPadMargin : 0;
          legend = new TLegend(0.03 + leftMarginAdder, 0.05 + bottomMarginAdder, 0.5 / (1 - leftMarginAdder), 0.45 / thisPadScale);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
          legend->AddEntry((TObject*)0, jetPtString.Data(), "");
          legend->AddEntry(systematicUncertaintyForPbPb[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], Form("PbPb %.0f-%.0f%%", card[kPbPb]->GetLowBinBorderCentrality(iCentrality), card[kPbPb]->GetHighBinBorderCentrality(iCentrality)), "p");
          legend->AddEntry(systematicUncertaintyForPp[kUncorrelatedUncertainty][iJetPt][iTrackPt], "pp", "p");

          legend->Draw();

        }  // Jet pT loop
      }    // Centrality loop

      // Draw all necessary CMS text to the plot
      bigCanvas[iTrackPt]->cd(0);

      mainTitle->SetTextFont(62);
      mainTitle->SetTextSize(0.065);
      mainTitle->DrawLatexNDC(0.08, 0.93, "CMS");

      mainTitle->SetTextFont(42);
      mainTitle->SetTextSize(0.025);
      mainTitle->DrawLatexNDC(0.23, 0.96, "PbPb #sqrt{s_{NN}} = 5.02 TeV, 1.70 nb^{-1}");
      mainTitle->DrawLatexNDC(0.23, 0.92, "pp #sqrt{s} = 5.02 TeV, 302 pb^{-1}");
      mainTitle->DrawLatexNDC(0.53, 0.96, "anti-k_{T} R = 0.4");
      mainTitle->DrawLatexNDC(0.53, 0.92, "|#eta_{jet}| < 1.6");

      mainTitle->SetTextSize(0.035);
      if(iTrackPt == 4){
        mainTitle->DrawLatexNDC(0.75, 0.94, Form("p_{T}^{ch} > %.1f GeV", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt)));
      } else {
        mainTitle->DrawLatexNDC(0.75, 0.94, Form("p_{T}^{ch} > %.0f GeV", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt)));
      }

      // Save the figures to file
      if(saveFigures) {
        gPad->GetCanvas()->SaveAs(Form("figures/finalBigCanvas%s%s.pdf", saveComment.Data(), compactTrackPtString.Data()));
      }

    }  // Track pT loop
  } // If for drawing big canvases

  // ======================================================= //
  // ===== Drawing style with ratios in one big canvas ===== //
  // ======================================================= //
  
  if(drawBigCanvasRatios){
    
    // Draw all the distributions to big canvases
    SplitCanvas* bigRatioCanvas[nTrackPtBinsEEC];
    mainTitle = new TLatex();

    for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
      compactTrackPtString = Form("_T>%.1f", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));
      compactTrackPtString.ReplaceAll(".","v");

      // Draw a big canvas and put all the plots in it
      bigRatioCanvas[iTrackPt] = new SplitCanvas(Form("bigRatioCanvas%d", iTrackPt), "", 1800, 1500);
      bigRatioCanvas[iTrackPt]->SetMargin(0.08, 0.01, 0.08, 0.11);
      bigRatioCanvas[iTrackPt]->DivideNeatly(4, 4);

      bottomRowScale = bigRatioCanvas[iTrackPt]->GetBottomRowScale();
      bottomPadMargin = bigRatioCanvas[iTrackPt]->GetBottomPadMargin();
      leftPadMargin = bigRatioCanvas[iTrackPt]->GetLeftPadMargin();

      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++) {
          jetPtString = Form("%.0f < jet p_{T} < %.0f", card[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));

          canvasIndex = (iCentrality - firstDrawnCentralityBin) * (lastDrawnJetPtBinEEC - firstDrawnJetPtBinEEC + 1) + (iJetPt - firstDrawnJetPtBinEEC);
          bigRatioCanvas[iTrackPt]->CD(canvasIndex);

          gPad->SetLogy(0); // Linear y-axis scale for ratios
          gPad->SetLogx();

          // The titles in the bottow row need to be scaled, since it has more margin than other rows
          thisPadScale = (iCentrality == lastDrawnCentralityBin) ? bottomRowScale : 1;

          // Set the axis titles and labels
          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleOffset(0.85);
          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleSize(0.2 * thisPadScale);
          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelOffset(0.01);
          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelSize(0.12 * thisPadScale);

          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleOffset(0.5 / thisPadScale);
          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleSize(0.2 * thisPadScale);
          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelOffset(0.01 / thisPadScale);
          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelSize(0.12 * thisPadScale);

          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitle("#frac{PbPb}{pp}");
          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetStats(0);
          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetNdivisions(505);
          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitle("#Deltar");
          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetNdivisions(505);

          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetTitle("");

          // Set the drawing ranges
          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.008, 0.39);
          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.2, 1.8);
          systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);
          systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);

          // Set the drawing style for histograms
          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin]);
          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);
          energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin]);
          energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);

          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin], 0.4);
          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
          energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
          energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);

          // Different drawing styles for correlated uncertainties
          systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorUpPbPb[iCentrality - firstDrawnCentralityBin], 0.4);
          systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorDownPbPb[iCentrality - firstDrawnCentralityBin], 0.4);
          systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
          systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);
          systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
          systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);

          // Draw the histograms
          systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Draw("e2");
          systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          systematicUncertaintyPbPbToPpRatio[kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->Draw("same,p");

          // Create a legend for each pad
          leftMarginAdder = (iJetPt == firstDrawnJetPtBinEEC) ? leftPadMargin : 0;
          bottomMarginAdder = (iCentrality == lastDrawnCentralityBin) ? bottomPadMargin : 0;
          legend = new TLegend(0.05 + leftMarginAdder, 0.05 + bottomMarginAdder, 0.5 / (1 - leftMarginAdder), 0.33 / (thisPadScale*0.81));
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
          legend->AddEntry((TObject*)0, jetPtString.Data(), "");
          legend->AddEntry(systematicUncertaintyPbPbToPpRatio[kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], Form("PbPb %.0f-%.0f%% / pp", card[kPbPb]->GetLowBinBorderCentrality(iCentrality), card[kPbPb]->GetHighBinBorderCentrality(iCentrality)), "p,e2");

          legend->Draw();

          // Draw a line to one
          oneLine->Draw();

        }  // Jet pT loop
      }    // Centrality loop

      // Draw all necessary CMS text to the plot
      bigRatioCanvas[iTrackPt]->cd(0);

      mainTitle->SetTextFont(62);
      mainTitle->SetTextSize(0.065);
      mainTitle->DrawLatexNDC(0.08, 0.93, "CMS");

      mainTitle->SetTextFont(42);
      mainTitle->SetTextSize(0.025);
      mainTitle->DrawLatexNDC(0.23, 0.96, "PbPb #sqrt{s_{NN}} = 5.02 TeV, 1.70 nb^{-1}");
      mainTitle->DrawLatexNDC(0.23, 0.92, "pp #sqrt{s} = 5.02 TeV, 302 pb^{-1}");
      mainTitle->DrawLatexNDC(0.53, 0.96, "anti-k_{T} R = 0.4");
      mainTitle->DrawLatexNDC(0.53, 0.92, "|#eta_{jet}| < 1.6");

      mainTitle->SetTextSize(0.035);
      if(iTrackPt == 4){
        mainTitle->DrawLatexNDC(0.75, 0.94, Form("p_{T}^{ch} > %.1f GeV", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt)));
      } else {
        mainTitle->DrawLatexNDC(0.75, 0.94, Form("p_{T}^{ch} > %.0f GeV", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt)));
      }

      // Save the figures to file
      if(saveFigures) {
        gPad->GetCanvas()->SaveAs(Form("figures/finalBigRatioCanvas%s%s.pdf", saveComment.Data(), compactTrackPtString.Data()));
      }

    }  // Track pT loop
  } // If for drawing big canvases

  if(drawDoubleRatios){

    // For the double ratios, only draw one example bin
    doubleRatioCentralityBinIndex1 = card[kPbPb]->FindBinIndexCentrality(doubleRatioCentralityBin1);
    doubleRatioCentralityBinIndex2 = card[kPbPb]->FindBinIndexCentrality(doubleRatioCentralityBin2);
    doubleRatioJetPtBinIndex = card[kPbPb]->FindBinIndexJetPtEEC(doubleRatioJetPtBin);

    // Create a SplitCanvas object for the two figures
    SplitCanvas* doubleRatioCanvas = new SplitCanvas("doubleRatioCanvas", "", 900, 400);
    doubleRatioCanvas->SetMargin(0.16, 0.01, 0.2, 0.05);
    doubleRatioCanvas->DivideNeatly(1, 2);

    bottomRowScale = doubleRatioCanvas->GetBottomRowScale();
    bottomPadMargin = doubleRatioCanvas->GetBottomPadMargin();
    leftPadMargin = doubleRatioCanvas->GetLeftPadMargin();

    mainTitle = new TLatex();
    canvasIndex = 0;

    // Go through the two centrality bin
    for(int iCentrality = doubleRatioCentralityBinIndex1; iCentrality <= doubleRatioCentralityBinIndex2; iCentrality += (doubleRatioCentralityBinIndex2-doubleRatioCentralityBinIndex1)){
      centralityString = Form("Cent: %.0f-%.0f%%", card[kPbPb]->GetLowBinBorderCentrality(iCentrality), card[kPbPb]->GetHighBinBorderCentrality(iCentrality));
      compactCentralityString = Form("_C=%.0f-%.0f", card[kPbPb]->GetLowBinBorderCentrality(iCentrality), card[kPbPb]->GetHighBinBorderCentrality(iCentrality));

      doubleRatioCanvas->CD(canvasIndex++);

      gPad->SetLogy(false);
      gPad->SetLogx();

      // Set the axis titles and labels
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetTitleOffset(1);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetTitleSize(0.09);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetLabelOffset(0.01);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetLabelSize(0.07);

      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetTitleOffset(1.4);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetTitleSize(0.08);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetLabelOffset(0.01);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetLabelSize(0.07);

      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetTitle("#frac{PbPb/pp (p_{T}^{ch} > 3 GeV)}{PbPb/pp (p_{T}^{ch} > 2 GeV)}");
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->CenterTitle();
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->SetStats(0);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetNdivisions(505);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->CenterTitle();
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetTitle("#Deltar");
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->CenterTitle();
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetNdivisions(505);

      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->SetTitle("");

      // Set the drawing ranges
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetRangeUser(0.008, 0.39);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetRangeUser(0.8, 1.2);

      // Set the drawing style for histograms
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin]);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->SetMarkerSize(1.2);
      energyEnergyCorrelatorDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin]);
      energyEnergyCorrelatorDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->SetMarkerSize(1.2);

      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin], 0.4);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
      energyEnergyCorrelatorDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
      energyEnergyCorrelatorDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);

      // Draw the histograms
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->Draw("e2");
      energyEnergyCorrelatorDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->Draw("same,p");

      // Show the centrality bin in the legend
      leftMarginAdder = (iCentrality == doubleRatioCentralityBinIndex1) ? leftPadMargin-0.03 : 0;
      legend = new TLegend(0.05+leftMarginAdder, 0.3, 0.45 / (1-leftMarginAdder), 0.5);
      legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.07); legend->SetTextFont(62);
      legend->AddEntry((TObject*)0, Form("%.f < jet p_{T} < %.0f GeV",  card[kPbPb]->GetLowBinBorderJetPtEEC(doubleRatioJetPtBinIndex), card[kPbPb]->GetHighBinBorderJetPtEEC(doubleRatioJetPtBinIndex)), "");
      legend->AddEntry(systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex], Form("Centrality: %.0f-%.0f%%", card[kPbPb]->GetLowBinBorderCentrality(iCentrality), card[kPbPb]->GetHighBinBorderCentrality(iCentrality)), "p");

      legend->Draw();

      // Draw the line to one
      oneLine->Draw();

    }

    // Draw all necessary CMS text to the plot
    doubleRatioCanvas->cd(0);

    mainTitle->SetTextFont(62);
    mainTitle->SetTextSize(0.09);
    mainTitle->DrawLatexNDC(0.19, 0.81, "CMS");

    mainTitle->SetTextFont(42);
    mainTitle->SetTextSize(0.07);
    mainTitle->DrawLatexNDC(0.595, 0.83, "PbPb #sqrt{s_{NN}} = 5.02 TeV, 1.70 nb^{-1}");
    mainTitle->DrawLatexNDC(0.595, 0.73, "pp #sqrt{s} = 5.02 TeV, 302 pb^{-1}");
    mainTitle->DrawLatexNDC(0.35, 0.83, "anti-k_{T} R = 0.4");
    mainTitle->DrawLatexNDC(0.35, 0.73, "|#eta_{jet}| < 1.6");

    box = new TBox();
    box->SetFillColor(kWhite);
    box->DrawBox(0.12,0.9, 0.155, 0.99);
    mainTitle->DrawLatexNDC(0.121, 0.93, "1.2");

    // Save the figures to file
    if(saveFigures) {
      gPad->GetCanvas()->SaveAs(Form("figures/finalDoubleRatio%s.pdf", saveComment.Data()));
    }

  } // If for drawing double ratios

  if(drawDoubleRatioToSingleCanvas){

    // For the double ratios, only draw one example bin
    doubleRatioJetPtBinIndex = card[kPbPb]->FindBinIndexJetPtEEC(doubleRatioJetPtBin);

    // Create a TCanvas for the double ratio
    TCanvas* theGreatCanvasOfDoubleRatio = new SplitCanvas("theGreatCanvasOfDoubleRatio", "", 1200, 800);
    theGreatCanvasOfDoubleRatio->SetMargin(0.24, 0.01, 0.2, 0.05);
    theGreatCanvasOfDoubleRatio->cd();

    gPad->SetLogy(false);
    gPad->SetLogx();

    mainTitle = new TLatex();
    canvasIndex = 0;

    // Set the drawing and label styles
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){

      // Set the axis titles and labels
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetTitleOffset(1);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetTitleSize(0.09);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetLabelOffset(0.01);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetLabelSize(0.07);

      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetTitleOffset(1.2);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetTitleSize(0.08);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetLabelOffset(0.01);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetLabelSize(0.07);

      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetTitle("#frac{PbPb/pp (p_{T}^{ch} > 3 GeV)}{PbPb/pp (p_{T}^{ch} > 2 GeV)}");
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->CenterTitle();
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->SetStats(0);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetNdivisions(505);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->CenterTitle();
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetTitle("#Deltar");
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->CenterTitle();
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetNdivisions(505);

      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->SetTitle("");

      // Set the drawing ranges
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetRangeUser(0.008, 0.39);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetRangeUser(0.8, 1.2);

      // Set the drawing styles 
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin]);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->SetMarkerSize(1.2);
      energyEnergyCorrelatorDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin]);
      energyEnergyCorrelatorDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->SetMarkerSize(1.2);

      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin], 0.4);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
      energyEnergyCorrelatorDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
      energyEnergyCorrelatorDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
    }

    // Draw the histograms
    systematicUncertaintyDoubleRatio[firstDrawnCentralityBin][doubleRatioJetPtBinIndex]->Draw("e2");

    for(int iCentrality = firstDrawnCentralityBin+1; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->Draw("same,e2");
    }
    for(int iCentrality = lastDrawnCentralityBin; iCentrality >= firstDrawnCentralityBin; iCentrality--){
      energyEnergyCorrelatorDoubleRatio[iCentrality][doubleRatioJetPtBinIndex]->Draw("same,p");
    }

    // Show the centrality bin in the legend
    legend = new TLegend(0.28, 0.22, 0.53, 0.54);
    legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.055); legend->SetTextFont(62);
    legend->AddEntry((TObject*)0, Form("%.f < jet p_{T} < %.0f GeV",  card[kPbPb]->GetLowBinBorderJetPtEEC(doubleRatioJetPtBinIndex), card[kPbPb]->GetHighBinBorderJetPtEEC(doubleRatioJetPtBinIndex)), "");
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      legend->AddEntry(systematicUncertaintyDoubleRatio[iCentrality][doubleRatioJetPtBinIndex], Form("Centrality: %.0f-%.0f%%", card[kPbPb]->GetLowBinBorderCentrality(iCentrality), card[kPbPb]->GetHighBinBorderCentrality(iCentrality)), "p");
    }

    legend->Draw();

    // Draw the line to one
    oneLine->Draw();

    mainTitle->SetTextFont(62);
    mainTitle->SetTextSize(0.07);
    mainTitle->DrawLatexNDC(0.542, 0.666, "CMS");

    mainTitle->SetTextFont(42);
    mainTitle->SetTextSize(0.055);
    mainTitle->DrawLatexNDC(0.5, 0.87, "PbPb #sqrt{s_{NN}} = 5.02 TeV, 1.70 nb^{-1}");
    mainTitle->DrawLatexNDC(0.5, 0.78, "pp #sqrt{s} = 5.02 TeV, 302 pb^{-1}");
    mainTitle->DrawLatexNDC(0.27, 0.87, "anti-k_{T} R = 0.4");
    mainTitle->DrawLatexNDC(0.27, 0.78, "|#eta_{jet}| < 1.6");

    // Save the figures to file
    if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/finalDoubleRatioSingleCanvas%s.pdf", saveComment.Data()));
    }

  } // If for drawing double ratios
}
