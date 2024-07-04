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

  // Select the weight exponent that is used for the figures
  int weightExponent = 1;

  // Check that the selected weight exponent is reasonable
  if(weightExponent < 1 || weightExponent > 2){
    cout << "ERROR! You selected " << weightExponent << " for the value for weight exponent." << endl;
    cout << "However, only 1 and 2 are implemented." << endl;
    cout << "Please select one of the implemented values." << endl;
    return;
  }
  
  const int nWeightExponents = 2;

  // Input files
  TString inputFileName[kNDataTypes][nWeightExponents];
  inputFileName[kPbPb][0] = "data/eecAnalysis_akFlowJet_nominalEnergyWeight_combinedMixedConeBackground_fixedCovarianceMatrix_processed_2024-05-28.root";
  inputFileName[kPbPb][1] = "data/eecAnalysis_akFlowJet_energyWeightSquared_combinedMixedConeBackground_fixedCovarianceMatrix_processed_2024-05-28.root";
  inputFileName[kPp][0] = "data/ppData_pfJets_wtaAxis_nominalEnergyWeight_optimizedUnfoldingBins_fixedCovarianceMatrix_jet60or80triggers_addLowPtBins_processed_2024-04-18.root";
  inputFileName[kPp][1] = "data/ppData_pfJets_wtaAxis_energyWeightSquared_optimizedUnfoldingBins_fixedCovarianceMatrix_jet60or80triggers_unfoldingWithCovariance_processed_2024-01-23.root";
  TString uncertaintyFileName[kNDataTypes][nWeightExponents];
  uncertaintyFileName[kPbPb][0] = "systematicUncertainties/systematicUncertainties_PbPb_nominalEnergyWeight_combinedMixedConeBackground_noMCnonClosure_2024-05-28.root";
  uncertaintyFileName[kPbPb][1] = "systematicUncertainties/systematicUncertainties_PbPb_energyWeightSquared_combinedMixedConeBackground_noMCnonClosure_2024-05-28.root";
  uncertaintyFileName[kPp][0] = "systematicUncertainties/systematicUncertainties_pp_nominalEnergyWeight_noMCnonClosure_2024-05-02.root";
  uncertaintyFileName[kPp][1] = "systematicUncertainties/systematicUncertainties_pp_energyWeightSquared_noMCnonClosure_2024-05-02.root";

  TString shiftedPtFileName[nWeightExponents];
  shiftedPtFileName[0] = "data/ppData_pfJets_wtaAxis_shiftedJetPtBins_nominalEnergyWeight_jet60or80trigger_unfoldWithShiftedResponse_processed_2024-02-19.root";
  shiftedPtFileName[1] = "data/ppData_pfJets_wtaAxis_shiftedJetPtBins_energyWeightSquared_jet60or80trigger_unfoldWithShiftedResponse_processed_2024-02-19.root";

  // eecAnalysis_akFlowJet_nominalEnergyWeight_optimizedUnfoldingBins_fixedCovarianceMatrix_updatedBackgroundSubtraction_processed_2024-02-23.root
  // eecAnalysis_akFlowJet_energyWeightSquared_optimizedUnfoldingBins_fixedCovarianceMatrix_updatedBackgroundSubtraction_processed_2024-02-23.root
  // ppData_pfJets_wtaAxis_nominalEnergyWeight_optimizedUnfoldingBins_fixedCovarianceMatrix_jet60or80triggers_unfoldingWithCovariance_processed_2024-01-23.root
  // ppData_pfJets_wtaAxis_nominalEnergyWeight_optimizedUnfoldingBins_jet60or80triggers_unfoldingWithNominalSmear_processed_2024-01-17.root
  // ppData_pfJets_wtaAxis_energyWeightSquared_optimizedUnfoldingBins_fixedCovarianceMatrix_jet60or80triggers_unfoldingWithCovariance_processed_2024-01-23.root
  // ppData_pfJets_wtaAxis_energyWeightSquared_optimizedUnfoldingBins_jet60or80triggers_unfoldingWithNominalSmear_processed_2024-01-17.root
  
  TFile* inputFile[kNDataTypes][nWeightExponents];
  TFile* uncertaintyFile[kNDataTypes][nWeightExponents];
  TFile* shiftedPtFile[nWeightExponents];
  EECCard* card[kNDataTypes][nWeightExponents];
  EECCard* uncertaintyCard[kNDataTypes][nWeightExponents];
  EECCard* shiftedPtCard[nWeightExponents];
  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    for(int iFile = 0; iFile < kNDataTypes; iFile++){

      // Load the input file
      inputFile[iFile][iWeightExponent] = TFile::Open(inputFileName[iFile][iWeightExponent]);

      // Check that the input file exists
      if(inputFile[iFile][iWeightExponent] == NULL){
        cout << "Error! The file " << inputFileName[iFile][iWeightExponent].Data() << " does not exist!" << endl;
        cout << "Maybe you forgot the data/ folder path?" << endl;
        cout << "Will not execute the code" << endl;
        return;
      }
    
      // Load the card from the file and read the collision system
      card[iFile][iWeightExponent] = new EECCard(inputFile[iFile][iWeightExponent]);

      // Load the uncertainty file
      uncertaintyFile[iFile][iWeightExponent] = TFile::Open(uncertaintyFileName[iFile][iWeightExponent]);

      // Check that the uncertianty file exists
      if(uncertaintyFile[iFile][iWeightExponent] == NULL){
        cout << "Error! The file " << uncertaintyFileName[iFile][iWeightExponent].Data() << " does not exist!" << endl;
        cout << "Maybe you forgot the systematicUncertainties/ folder path?" << endl;
        cout << "Will not execute the code" << endl;
        return;
      }
    
      // Load the card from the file and read the collision system
      uncertaintyCard[iFile][iWeightExponent] = new EECCard(uncertaintyFile[iFile][iWeightExponent]);
    }  // File loop

    // Load the shifted pT file
    shiftedPtFile[iWeightExponent] = TFile::Open(shiftedPtFileName[iWeightExponent]);

    // Check that the input file exists
    if(shiftedPtFile[iWeightExponent] == NULL){
      cout << "Error! The file " << shiftedPtFileName[iWeightExponent].Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }
    
    // Load the card from the file and read the collision system
    shiftedPtCard[iWeightExponent] = new EECCard(shiftedPtFile[iWeightExponent]);

  } // Weight exponent loop

  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Find the number of bins from the PbPb card
  const int nCentralityBins = card[kPbPb][0]->GetNCentralityBins();
  const int nJetPtBinsEEC = card[kPbPb][0]->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card[kPbPb][0]->GetNTrackPtBinsEEC();
  
  // The final results are available for all the bins that are unfolded
  int firstDrawnCentralityBin[nWeightExponents];
  int lastDrawnCentralityBin[nWeightExponents];
  
  int firstDrawnJetPtBinEEC[nWeightExponents];
  int lastDrawnJetPtBinEEC[nWeightExponents];
  
  int firstDrawnTrackPtBinEEC[nWeightExponents];
  int lastDrawnTrackPtBinEEC[nWeightExponents];

  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    firstDrawnCentralityBin[iWeightExponent] = card[kPbPb][iWeightExponent]->GetFirstUnfoldedCentralityBin();
    lastDrawnCentralityBin[iWeightExponent] = card[kPbPb][iWeightExponent]->GetLastUnfoldedCentralityBin();
    firstDrawnJetPtBinEEC[iWeightExponent] = card[kPbPb][iWeightExponent]->GetFirstUnfoldedJetPtBin();
    lastDrawnJetPtBinEEC[iWeightExponent] = card[kPbPb][iWeightExponent]->GetLastUnfoldedJetPtBin();
    firstDrawnTrackPtBinEEC[iWeightExponent] = card[kPbPb][iWeightExponent]->GetFirstUnfoldedTrackPtBin();
    lastDrawnTrackPtBinEEC[iWeightExponent] = card[kPbPb][iWeightExponent]->GetLastUnfoldedTrackPtBin();
  }

  // Choose which plots to draw
  bool drawIndividualPlotsAllCentralities = false;
  bool drawBigCanvasDistributions = false;
  bool drawBigCanvasRatios = false;
  bool drawDoubleRatios = false;
  bool drawDoubleRatioToSingleCanvas = false;
  bool drawBigCanvasAllDistributions = false; // Draw distributions with all defined weight exponents to the same figure
  bool drawBigCanvasAllRatios = false; // Draw ratios with all defined energy weight exponents to the same figure
  bool drawLetterPaperDistributions = false; // Draw the energy-energy correlator distribution for letter paper format
  bool drawLetterPaperRatios = false; // Draw the energy-energy correlators ratios for letter paper format

  // Supplementary plots
  bool drawSupplementaryEECJetPt = true; // Draw all jet pT bins for selected centrality, track pT, and energy weight to the same plot
  bool drawSupplementaryEECCentrality = false; // Draw all centrality bins for selected jet pT, track pT and energy weight to the smae plot
  bool drawSupplementaryEECRatioCentrality = false; // Draw PbPb/pp ratios from all centrality bins to the same figure
  bool drawSupplementaryEECRatioWeightExponent = false; // Draw PbPb/pp ratios from all weight exponents to the same figure

  // Configuration for supplementary plots
  double supplementaryLegendTextSize = 0.045;

  // Preliminary tag
  bool addPreliminaryTag = true;

  bool drawVerticalLines = false; // Draw illustrative vertical lines
  bool drawShiftedPtRatio = false;
  if(drawBigCanvasAllRatios || drawBigCanvasAllDistributions || drawLetterPaperDistributions || drawLetterPaperRatios || drawSupplementaryEECRatioWeightExponent) weightExponent = 1; // For technical reasons, use number 1 here when both ratios are drawn to same figure

  // Normalize all distributions to 2 GeV integral
  bool normalizeTo2GeV = false;
  int trackPtBinFor2GeV[nWeightExponents];
  std::pair<double, double> trackPtCutsForDoubleDatio = std::make_pair(1.0, 2.0);
  std::pair<int, int> trackPtBinsForDoubleRatio = std::make_pair(card[kPbPb][0]->GetBinIndexTrackPtEEC(trackPtCutsForDoubleDatio.first), card[kPbPb][0]->GetBinIndexTrackPtEEC(trackPtCutsForDoubleDatio.second));

  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    trackPtBinFor2GeV[iWeightExponent] = card[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(2.0);
  }

  // Binning for letter paper figures
  std::pair<double, double> letterPaperCentralityBin1 = std::make_pair(0.0,10.0);
  std::pair<double, double> letterPaperCentralityBin2 = std::make_pair(50.0,90.0);
  std::pair<double, double> letterPaperJetPtBin = std::make_pair(140.0, 160.0);
  double letterPaperTrackPtCut1 = 1;
  double letterPaperTrackPtCut2 = 2;
  std::vector<int> letterPaperCentralityBinIndex;
  std::vector<int> letterPaperJetPtBinIndex;
  std::vector<int> letterPaperTrackPtCutIndex;

  // Select the bins to be drawn for double ratio plots
  std::pair<double, double> doubleRatioCentralityBin1 = std::make_pair(0.0,10.0);
  std::pair<double, double> doubleRatioCentralityBin2 = std::make_pair(10.0,30.0);
  std::pair<double, double> doubleRatioJetPtBin = std::make_pair(180,200);
  int doubleRatioCentralityBinIndex1;
  int doubleRatioCentralityBinIndex2;
  int doubleRatioJetPtBinIndex;
  
  // Save the final plots
  const bool saveFigures = true;
  TString energyWeightString[nWeightExponents] = {"_nominalEnergyWeight", "_energyWeightSquared"};
  TString saveComment =  "_preliminaryTag";
  TString figureFormat = "pdf";
  if(!drawBigCanvasAllRatios && !drawBigCanvasAllDistributions && !drawLetterPaperDistributions && !drawLetterPaperRatios && !drawSupplementaryEECRatioWeightExponent){
    saveComment.Prepend(energyWeightString[weightExponent-1]);
  }

  // Ratio zoom settings
  std::pair<double, double> ratioZoom = std::make_pair(0.4, 1.6);
  std::pair<double, double> analysisDeltaR = std::make_pair(0.008, 0.39); // DeltaR span in which the analysis is done
  
  // Marker colors and styles
  int markerStylePbPb[] = {kFullSquare, kFullCircle, kFullCross, kFullFourTrianglesPlus};
  int markerStylePp = kFullDiamond;
  int markerColorPbPb[] = {kRed, kBlue, kMagenta, kGreen+3};
  int markerColorPp = kBlack;
  int bandColorUpPbPb[] = {kOrange+7, kViolet-3, kPink-3, kOrange-3};
  int bandColorDownPbPb[] = {kPink+9, kAzure+8, kViolet+6, kSpring};
  int markerStylePpEnergyWeight[] = {kFullDiamond, kFullDoubleDiamond};
  int markerColorPpEnergyWeight[] = {kViolet-2, kCyan+1};

  // Definition on how much the jet pT is shifted for shifted pT figures
  std::vector<std::pair<std::pair<double,double>,std::pair<double,double>>> jetPtShiftDefinition;
  jetPtShiftDefinition.push_back(std::make_pair(std::make_pair(120,140),std::make_pair(135,155)));
  jetPtShiftDefinition.push_back(std::make_pair(std::make_pair(140,160),std::make_pair(155,175)));
  jetPtShiftDefinition.push_back(std::make_pair(std::make_pair(160,180),std::make_pair(175,195)));
  jetPtShiftDefinition.push_back(std::make_pair(std::make_pair(180,200),std::make_pair(195,215)));

  TLine* lineDrawer = new TLine();
  lineDrawer->SetLineStyle(2);
  lineDrawer->SetLineColor(kBlack);

  // =============================================== //
  // Read the histograms from the histogram managers //
  // =============================================== //

  // Create histogram managers for result files and systematic uncertainty organizers for systematic uncertainty files
  EECHistogramManager* histograms[kNDataTypes][nWeightExponents];
  EECHistogramManager* shiftedPtHistograms[nWeightExponents];
  SystematicUncertaintyOrganizer* uncertainties[kNDataTypes][nWeightExponents];
  
  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    for(int iDataType = 0; iDataType < kNDataTypes; iDataType++){

      // Create a new histogram manager
      histograms[iDataType][iWeightExponent] = new EECHistogramManager(inputFile[iDataType][iWeightExponent], card[iDataType][iWeightExponent]);
  
      // Load all unfolded energy-energy correlators
      histograms[iDataType][iWeightExponent]->SetLoadEnergyEnergyCorrelators(true);
      histograms[iDataType][iWeightExponent]->SetCentralityBinRange(0, card[iDataType][iWeightExponent]->GetNCentralityBins());
      histograms[iDataType][iWeightExponent]->SetTrackPtBinRangeEEC(0, card[iDataType][iWeightExponent]->GetNTrackPtBinsEEC());
      histograms[iDataType][iWeightExponent]->SetJetPtBinRangeEEC(0, card[iDataType][iWeightExponent]->GetNJetPtBinsEEC());

      // Load the histograms from the file
      histograms[iDataType][iWeightExponent]->LoadProcessedHistograms();

      // Create a new systematic uncertainty organizer
      uncertainties[iDataType][iWeightExponent] = new SystematicUncertaintyOrganizer(uncertaintyFile[iDataType][iWeightExponent]);

    } // Data type loop

    // Create a new histogram manager for histograms with shifted pT
    if(drawShiftedPtRatio){
      shiftedPtHistograms[iWeightExponent] = new EECHistogramManager(shiftedPtFile[iWeightExponent], shiftedPtCard[iWeightExponent]);
  
      // Load all unfolded energy-energy correlators
      shiftedPtHistograms[iWeightExponent]->SetLoadEnergyEnergyCorrelators(true);
      shiftedPtHistograms[iWeightExponent]->SetCentralityBinRange(0, shiftedPtCard[iWeightExponent]->GetNCentralityBins());
      shiftedPtHistograms[iWeightExponent]->SetTrackPtBinRangeEEC(0, shiftedPtCard[iWeightExponent]->GetNTrackPtBinsEEC());
      shiftedPtHistograms[iWeightExponent]->SetJetPtBinRangeEEC(0, shiftedPtCard[iWeightExponent]->GetNJetPtBinsEEC());

      // Load the histograms from the file
      shiftedPtHistograms[iWeightExponent]->LoadProcessedHistograms();
    }

  } // Weight exponent loop
 
  // Energy-energy correlators and PbPb to pp ratios
  TH1D* energyEnergyCorrelatorRawPbPb[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorRawPp[nWeightExponents][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorSignalPbPb[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorSignalPp[nWeightExponents][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorSignalShiftedPt[nWeightExponents][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorPbPbToPpRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorShiftedPtRatio[nWeightExponents][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyForPbPb[nWeightExponents][knSystematicUncertaintyTypes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyForPp[nWeightExponents][knSystematicUncertaintyTypes][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyPbPbToPpRatio[nWeightExponents][knSystematicUncertaintyTypes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];

  // Double ratios from energy-energy correlator histograms
  TH1D* energyEnergyCorrelatorForDoubleRatioFromPbPb[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorForDoubleRatioFromPp[nWeightExponents][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorDoubleRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC];
  TH1D* systematicUncertaintyForDoubleRatioFromPbPb[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyForDoubleRatioFromPp[nWeightExponents][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyDoubleRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC];

  // Initialize histograms to NULL
  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        energyEnergyCorrelatorDoubleRatio[iWeightExponent][iCentrality][iJetPt] = NULL;
        systematicUncertaintyDoubleRatio[iWeightExponent][iCentrality][iJetPt] = NULL;
      }
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        energyEnergyCorrelatorRawPp[iWeightExponent][iJetPt][iTrackPt] = NULL;
        energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt] = NULL;
        energyEnergyCorrelatorSignalShiftedPt[iWeightExponent][iJetPt][iTrackPt] = NULL;
        energyEnergyCorrelatorShiftedPtRatio[iWeightExponent][iJetPt][iTrackPt] = NULL;

        energyEnergyCorrelatorForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt] = NULL;
        systematicUncertaintyForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt] = NULL;

        for(int iUncertainty = 0; iUncertainty < knSystematicUncertaintyTypes; iUncertainty++){
          systematicUncertaintyForPp[iWeightExponent][iUncertainty][iJetPt][iTrackPt] = NULL;
        } // Uncertainty type loop

        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          energyEnergyCorrelatorRawPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt] = NULL;
          energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt] = NULL;
          energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt] = NULL;

          energyEnergyCorrelatorForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt] = NULL;
          systematicUncertaintyForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt] = NULL;

          for(int iUncertainty = 0; iUncertainty < knSystematicUncertaintyTypes; iUncertainty++){
            systematicUncertaintyForPbPb[iWeightExponent][iUncertainty][iCentrality][iJetPt][iTrackPt] = NULL;
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][iUncertainty][iCentrality][iJetPt][iTrackPt] = NULL;
          }
        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // Weight exponent loop

  // Read the histograms from managers
  int iTrackPtMatchedPp, iTrackPtMatchedPbPbUncertainty, iTrackPtMatchedPpUncertainty, iTrackPtMatchedShifted;
  int iJetPtMatchedPp, iJetPtMatchedPbPbUncertainty, iJetPtMatchedPpUncertainty, iJetPtMatchedShifted;
  int iCentralityMatched;
  int referenceTrackPtBin, lowAnalysisBin, highAnalysisBin;
  double lowPtIntegral, highPtIntegral;
  double epsilon = 0.0001;

  // Define helper histograms to determine the uncertainties relevant for the double ratio
  TH1D* uncertaintyTrackSelection;
  TH1D* uncertaintyBackgroundSubtraction;
  TH1D* uncertaintyTrackPairEfficiency;
  TH1D* uncertaintyMCnonclosure;
  TH1D* uncertaintySignalToBackgroundRatio;
  double trackCorrelation;
  std::pair<double,double> currentJetPtBin;
  std::pair<double,double> shiftedJetPtBin;

  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    for(int iJetPt = firstDrawnJetPtBinEEC[iWeightExponent]; iJetPt <= lastDrawnJetPtBinEEC[iWeightExponent]; iJetPt++){

      // Find the current jet pT bin
      currentJetPtBin = card[kPbPb][iWeightExponent]->GetBinBordersJetPtEEC(iJetPt);

      // Find the shifted jet pT bin correcponding to the current pT bin
      for(auto ptMatchingRecipe : jetPtShiftDefinition){
        if(ptMatchingRecipe.first == currentJetPtBin){
          shiftedJetPtBin = ptMatchingRecipe.second;
          break;
        }
        shiftedJetPtBin = std::make_pair(-1,-1);
      }

      // Check that we found a shift definition
      if(shiftedJetPtBin.first == -1){
        cout << "ERROR! Could not find matching pT bin for jet pT bin " << currentJetPtBin.first << "-" << currentJetPtBin.second << endl;
      }

      iJetPtMatchedPp = card[kPp][iWeightExponent]->FindBinIndexJetPtEEC(currentJetPtBin);
      iJetPtMatchedPbPbUncertainty = uncertaintyCard[kPbPb][iWeightExponent]->FindBinIndexJetPtEEC(currentJetPtBin);
      iJetPtMatchedPpUncertainty = uncertaintyCard[kPp][iWeightExponent]->FindBinIndexJetPtEEC(currentJetPtBin);
      if(drawShiftedPtRatio) iJetPtMatchedShifted = shiftedPtCard[iWeightExponent]->FindBinIndexJetPtEEC(shiftedJetPtBin);

      // Before the main track pT loop, loop over all raw energy-energy correlators to have the information available to scale the tracking related uncertainties in double ratios.
      for(int iTrackPt = firstDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt <= lastDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt++){
        iTrackPtMatchedPp = card[kPp][iWeightExponent]->FindBinIndexTrackPtEEC(card[kPbPb][iWeightExponent]->GetBinBordersTrackPtEEC(iTrackPt));

        // Read the raw pp energy-energy correlator distributions
        energyEnergyCorrelatorRawPp[iWeightExponent][iJetPt][iTrackPt] = histograms[kPp][iWeightExponent]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatchedPp, iTrackPtMatchedPp);

        for(int iCentrality = firstDrawnCentralityBin[iWeightExponent]; iCentrality <= lastDrawnCentralityBin[iWeightExponent]; iCentrality++){

          // Read the raw PbPb energy-energy correlator distributions
          energyEnergyCorrelatorRawPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt] = histograms[kPbPb][iWeightExponent]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt);
        }
      }

      // Then go to the main track pT loop
      for(int iTrackPt = firstDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt <= lastDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt++){
        iTrackPtMatchedPp = card[kPp][iWeightExponent]->FindBinIndexTrackPtEEC(card[kPbPb][iWeightExponent]->GetBinBordersTrackPtEEC(iTrackPt));
        iTrackPtMatchedPbPbUncertainty = uncertaintyCard[kPbPb][iWeightExponent]->FindBinIndexTrackPtEEC(card[kPbPb][iWeightExponent]->GetBinBordersTrackPtEEC(iTrackPt));
        iTrackPtMatchedPpUncertainty = uncertaintyCard[kPp][iWeightExponent]->FindBinIndexTrackPtEEC(card[kPbPb][iWeightExponent]->GetBinBordersTrackPtEEC(iTrackPt));
        if(drawShiftedPtRatio) iTrackPtMatchedShifted = shiftedPtCard[iWeightExponent]->FindBinIndexTrackPtEEC(card[kPbPb][iWeightExponent]->GetBinBordersTrackPtEEC(iTrackPt));

        // Read the pp histograms that do not have centrality binning
        energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt] = histograms[kPp][iWeightExponent]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatchedPp, iTrackPtMatchedPp, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);
        systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt] = uncertainties[kPp][iWeightExponent]->GetUncorrelatedSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty);
        systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt] = uncertainties[kPp][iWeightExponent]->GetCorrelatedSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty);

        // Read the relevant systematic uncertainties for double ratio from pp
        uncertaintyBackgroundSubtraction = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kBackgroundSubtraction);
        uncertaintyTrackPairEfficiency = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kTrackPairEfficiency);
        uncertaintyMCnonclosure = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kMonteCarloNonClosure);
        uncertaintyTrackSelection = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kTrackSelection);

        // For double ratio, we can find the proper scaling factor for the tracking uncertainties by comparing integrals of raw
        // energy-energy correlators between different track pT bins. The higher track pT cuts are subsets of the lower cuts.
        // Thus integrals over the analysis region tell the degree of overlap of these regions.
        lowAnalysisBin = energyEnergyCorrelatorRawPp[iWeightExponent][iJetPt][iTrackPt]->GetXaxis()->FindBin(analysisDeltaR.first + epsilon);
        highAnalysisBin = energyEnergyCorrelatorRawPp[iWeightExponent][iJetPt][iTrackPt]->GetXaxis()->FindBin(analysisDeltaR.second - epsilon);
        lowPtIntegral = energyEnergyCorrelatorRawPp[iWeightExponent][iJetPt][trackPtBinsForDoubleRatio.first]->Integral(lowAnalysisBin, highAnalysisBin, "width");
        highPtIntegral = energyEnergyCorrelatorRawPp[iWeightExponent][iJetPt][trackPtBinsForDoubleRatio.second]->Integral(lowAnalysisBin, highAnalysisBin, "width");
        trackCorrelation = 1 - highPtIntegral / lowPtIntegral;

        // Now for the double ratio, we can combine the relevant uncertainties while scaling down the tracking related ones
        // by the expected overlap.
        systematicUncertaintyForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt] = (TH1D*) uncertaintyBackgroundSubtraction->Clone(Form("doubleRatioUncertaintyFromPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));
        for(int iBin = 1; iBin <= uncertaintyBackgroundSubtraction->GetNbinsX(); iBin++){
          systematicUncertaintyForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt]->SetBinError(iBin, TMath::Sqrt(TMath::Power(uncertaintyBackgroundSubtraction->GetBinError(iBin),2) + TMath::Power(uncertaintyTrackPairEfficiency->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintyTrackSelection->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintyMCnonclosure->GetBinError(iBin),2)));
        }

        // Read the histograms with shifted pT spectrum
        if(drawShiftedPtRatio){
          energyEnergyCorrelatorSignalShiftedPt[iWeightExponent][iJetPt][iTrackPt] = shiftedPtHistograms[iWeightExponent]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatchedShifted, iTrackPtMatchedShifted, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);
        }


        for(int iCentrality = firstDrawnCentralityBin[iWeightExponent]; iCentrality <= lastDrawnCentralityBin[iWeightExponent]; iCentrality++){
          iCentralityMatched = uncertaintyCard[kPbPb][iWeightExponent]->FindBinIndexCentrality(card[kPbPb][iWeightExponent]->GetBinBordersCentrality(iCentrality));

          // Read the PbPb histograms
          energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt] = histograms[kPbPb][iWeightExponent]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);
          systematicUncertaintyForPbPb[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt] = uncertainties[kPbPb][iWeightExponent]->GetUncorrelatedSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty);
          systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt] = uncertainties[kPbPb][iWeightExponent]->GetCorrelatedSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty);

          // Read the relevant systematic uncertainties for double ratio from PbPb
          uncertaintyBackgroundSubtraction = uncertainties[kPbPb][iWeightExponent]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kBackgroundSubtraction);
          uncertaintyTrackPairEfficiency = uncertainties[kPbPb][iWeightExponent]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kTrackPairEfficiency);
          uncertaintyTrackSelection = uncertainties[kPbPb][iWeightExponent]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kTrackSelection);
          uncertaintyMCnonclosure = uncertainties[kPbPb][iWeightExponent]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kMonteCarloNonClosure);
          uncertaintySignalToBackgroundRatio = uncertainties[kPbPb][iWeightExponent]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kSignalToBackgroundRatio);

          // For double ratio, we can find the proper scaling factor for the tracking uncertainties by comparing integrals of raw
          // energy-energy correlators between different track pT bins. The higher track pT cuts are subsets of the lower cuts.
          // Thus integrals over the analysis region tell the degree of overlap of these regions.
          lowAnalysisBin = energyEnergyCorrelatorRawPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(analysisDeltaR.first + epsilon);
          highAnalysisBin = energyEnergyCorrelatorRawPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(analysisDeltaR.second - epsilon);
          lowPtIntegral = energyEnergyCorrelatorRawPbPb[iWeightExponent][iCentrality][iJetPt][trackPtBinsForDoubleRatio.first]->Integral(lowAnalysisBin, highAnalysisBin, "width");
          highPtIntegral = energyEnergyCorrelatorRawPbPb[iWeightExponent][iCentrality][iJetPt][trackPtBinsForDoubleRatio.second]->Integral(lowAnalysisBin, highAnalysisBin, "width");
          trackCorrelation = 1 - highPtIntegral / lowPtIntegral;

          // Now for the double ratio, we can combine the relevant uncertainties while scaling down the tracking related ones
          // by the expected overlap.
          systematicUncertaintyForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt] = (TH1D*) uncertaintyBackgroundSubtraction->Clone(Form("doubleRatioUncertaintyFromPbPb%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          for(int iBin = 1; iBin <= uncertaintyBackgroundSubtraction->GetNbinsX(); iBin++){
            systematicUncertaintyForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetBinError(iBin, TMath::Sqrt(TMath::Power(uncertaintyBackgroundSubtraction->GetBinError(iBin),2) + TMath::Power(uncertaintyTrackPairEfficiency->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintyTrackSelection->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintySignalToBackgroundRatio->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintyMCnonclosure->GetBinError(iBin),2)));
          }
        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // Weight exponent loop

  // ================================================ //
  //   Normalize the histograms and calculate ratios  //
  // ================================================ //

  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    referenceTrackPtBin = trackPtBinFor2GeV[iWeightExponent];
    lowAnalysisBin = energyEnergyCorrelatorSignalPbPb[iWeightExponent][firstDrawnCentralityBin[iWeightExponent]][firstDrawnJetPtBinEEC[iWeightExponent]][firstDrawnTrackPtBinEEC[iWeightExponent]]->GetXaxis()->FindBin(analysisDeltaR.first + epsilon);
    highAnalysisBin = energyEnergyCorrelatorSignalPbPb[iWeightExponent][firstDrawnCentralityBin[iWeightExponent]][firstDrawnJetPtBinEEC[iWeightExponent]][firstDrawnTrackPtBinEEC[iWeightExponent]]->GetXaxis()->FindBin(analysisDeltaR.second - epsilon);
    for(int iJetPt = firstDrawnJetPtBinEEC[iWeightExponent]; iJetPt <= lastDrawnJetPtBinEEC[iWeightExponent]; iJetPt++){
      for(int iTrackPt = lastDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt >= firstDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt--){
        if(!normalizeTo2GeV) referenceTrackPtBin = iTrackPt;
        energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->Scale(1.0 / energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][referenceTrackPtBin]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
        systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt]->Scale(energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->GetBinContent(10) / systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetBinContent(10));
        systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt]->Scale(energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->GetBinContent(10) / systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt]->GetBinContent(10));
        systematicUncertaintyForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt]->Scale(energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->GetBinContent(10) / systematicUncertaintyForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt]->GetBinContent(10));

        // Shifted pT to regular pT ratios
        if(drawShiftedPtRatio){
          energyEnergyCorrelatorSignalShiftedPt[iWeightExponent][iJetPt][iTrackPt]->Scale(1.0 / energyEnergyCorrelatorSignalShiftedPt[iWeightExponent][iJetPt][referenceTrackPtBin]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

          energyEnergyCorrelatorShiftedPtRatio[iWeightExponent][iJetPt][iTrackPt] = (TH1D*) energyEnergyCorrelatorSignalShiftedPt[iWeightExponent][iJetPt][iTrackPt]->Clone(Form("shiftedPtRatio%d%d%d", iWeightExponent, iJetPt, iTrackPt));
          energyEnergyCorrelatorShiftedPtRatio[iWeightExponent][iJetPt][iTrackPt]->Divide(energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]);
        }

        for(int iCentrality = firstDrawnCentralityBin[iWeightExponent]; iCentrality <= lastDrawnCentralityBin[iWeightExponent]; iCentrality++){
          energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Scale(1.0 / energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][referenceTrackPtBin]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
          systematicUncertaintyForPbPb[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Scale(energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->GetBinContent(10) / systematicUncertaintyForPbPb[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetBinContent(10));
          systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Scale(energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->GetBinContent(10) / systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetBinContent(10));
          systematicUncertaintyForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Scale(energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->GetBinContent(10) / systematicUncertaintyForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->GetBinContent(10));

          energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt] = (TH1D*) energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Clone(Form("energyEnergyCorrelatorRatio%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Divide(energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]);

          systematicUncertaintyPbPbToPpRatio[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPbPb[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("systematicUncertaintyRatioUncorrelated%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          systematicUncertaintyPbPbToPpRatio[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Divide(systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt]);

          systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("systematicUncertaintyRatioCorrelated%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Divide(systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt]);
        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // Weight exponent loop

  // Reduce the statistical uncertainties by the overlapping statistics fraction for the double ratios
  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    for(int iJetPt = firstDrawnJetPtBinEEC[iWeightExponent]; iJetPt <= lastDrawnJetPtBinEEC[iWeightExponent]; iJetPt++){
      for(int iTrackPt = firstDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt <= lastDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt++){
        energyEnergyCorrelatorForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt] = (TH1D*) energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->Clone(Form("energyEnergyCorrelatorForDoubleRatioFromPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));

        for(int iBin = 1; iBin <= energyEnergyCorrelatorForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt]->GetNbinsX(); iBin++){
          energyEnergyCorrelatorForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt]->SetBinError(iBin, energyEnergyCorrelatorForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt]->GetBinError(iBin) * trackCorrelation);
        }

        for(int iCentrality = firstDrawnCentralityBin[iWeightExponent]; iCentrality <= lastDrawnCentralityBin[iWeightExponent]; iCentrality++){
          energyEnergyCorrelatorForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt] = (TH1D*) energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Clone(Form("energyEnergyCorrelatorForDoubleRatioFromPbPb%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));

          for(int iBin = 1; iBin <= energyEnergyCorrelatorForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->GetNbinsX(); iBin++){
            energyEnergyCorrelatorForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetBinError(iBin, energyEnergyCorrelatorForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->GetBinError(iBin) * trackCorrelation);
          }

        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // Weight exponent loop

  // After regular ratios have been calculated, proceed to calculating double ratio. Here we need to use different histograms as above to properly take into account systematic uncertainty cancellation due to correlated datasets.
  if(drawDoubleRatios || drawDoubleRatioToSingleCanvas){
    for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
      for(int iCentrality = firstDrawnCentralityBin[iWeightExponent]; iCentrality <= lastDrawnCentralityBin[iWeightExponent]; iCentrality++){
        for(int iJetPt = firstDrawnJetPtBinEEC[iWeightExponent]; iJetPt <= lastDrawnJetPtBinEEC[iWeightExponent]; iJetPt++){
          for(int iTrackPt = firstDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt <= lastDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt++){
            // Calculate the single ratios with properly handled double ratio uncertainties
            energyEnergyCorrelatorForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Divide(energyEnergyCorrelatorForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt]);
            systematicUncertaintyForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Divide(systematicUncertaintyForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt]);
          } // Track pT loop
          // Calculate the double ratios from the single ratios with properly handled uncertainties
          energyEnergyCorrelatorDoubleRatio[iWeightExponent][iCentrality][iJetPt] = (TH1D*) energyEnergyCorrelatorForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][trackPtBinsForDoubleRatio.second]->Clone(Form("energyEnergyCorrelatorDoubleRatio%d%d%d", iWeightExponent, iCentrality, iJetPt));
          energyEnergyCorrelatorDoubleRatio[iWeightExponent][iCentrality][iJetPt]->Divide(energyEnergyCorrelatorForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][trackPtBinsForDoubleRatio.first]);
          systematicUncertaintyDoubleRatio[iWeightExponent][iCentrality][iJetPt] = (TH1D*) systematicUncertaintyForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][trackPtBinsForDoubleRatio.second]->Clone(Form("systematicUncertaintyDoubleRatio%d%d%d", iWeightExponent, iCentrality, iJetPt));
          systematicUncertaintyDoubleRatio[iWeightExponent][iCentrality][iJetPt]->Divide(systematicUncertaintyForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][trackPtBinsForDoubleRatio.first]);
        } // Jet pT loop
      } // Centrality llop
    } // Weight exponent loop
  }

  // For illustration purposes, create up and down shifted uncertainty bands
  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    for(int iJetPt = firstDrawnJetPtBinEEC[iWeightExponent]; iJetPt <= lastDrawnJetPtBinEEC[iWeightExponent]; iJetPt++){
      for(int iTrackPt = firstDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt <= lastDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt++){

        systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertaintyShapeUp][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt]->Clone(Form("correlatedUpBandPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));
        systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertaintyShapeDown][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt]->Clone(Form("correlatedDownBandPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));
        calculateCorrelatedBands(systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt], systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertaintyShapeUp][iJetPt][iTrackPt], systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertaintyShapeDown][iJetPt][iTrackPt]);

        for(int iCentrality = firstDrawnCentralityBin[iWeightExponent]; iCentrality <= lastDrawnCentralityBin[iWeightExponent]; iCentrality++){

          systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("correlatedUpBand%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("correlatedDownBand%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          calculateCorrelatedBands(systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt], systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]);

          systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("correlatedUpBandForRatio%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("correlatedDownBandForRatio%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          calculateCorrelatedBands(systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt], systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]);

        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // Weight exponent loop

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
  TLegend* anotherLegend;
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

    for(int iJetPt = firstDrawnJetPtBinEEC[weightExponent-1]; iJetPt <= lastDrawnJetPtBinEEC[weightExponent-1]; iJetPt++){
      jetPtString = Form("%.0f < jet p_{T} < %.0f", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt));
      compactJetPtString = Form("_J=%.0f-%.0f", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt));

      for(int iTrackPt = firstDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt <= lastDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt++){
        trackPtString = Form("%.1f < track p_{T}", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt));
        compactTrackPtString = Form("_T>%.1f", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt));
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
        energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->SetMarkerStyle(kFullDiamond);
        energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->SetMarkerColor(kBlack);

        // Set the x-axis drawing range
        energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);

        // Draw the pp correlator to upper canves
        drawer->DrawHistogramToUpperPad(energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt], "#Deltar", "EEC Signal", " ", "p");
        legend->AddEntry(energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt], "pp", "p");

        // Draw the different centrality bins to the same plot
        for(int iCentrality = lastDrawnCentralityBin[weightExponent-1]; iCentrality >= firstDrawnCentralityBin[weightExponent-1]; iCentrality--) {
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->Draw("same,p");
          legend->AddEntry(energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt], Form("PbPb %.0f-%.0f%%", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality)), "p");
        }

        // Draw the legend to the upper pad
        legend->Draw();

        // Linear scale for the ratio
        drawer->SetLogY(false);

        // Set the axis drawing ranges
        energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][lastDrawnCentralityBin[weightExponent-1]][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
        energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][lastDrawnCentralityBin[weightExponent-1]][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

        // Set the style for histograms
        for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++) {
          energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
        }

        drawer->SetGridY(true);
        drawer->DrawHistogramToLowerPad(energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][lastDrawnCentralityBin[weightExponent-1]][iJetPt][iTrackPt], "#Deltar", "#frac{PbPb}{pp}", " ");
        for(int iCentrality = lastDrawnCentralityBin[weightExponent-1] - 1; iCentrality >= firstDrawnCentralityBin[weightExponent-1]; iCentrality--) {
          energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]->Draw("same,p");
        }
        drawer->SetGridY(false);

        // If a plot name is given, save the plot in a file
        if(saveFigures) {
          gPad->GetCanvas()->SaveAs(Form("figures/energyEnergyCorrelator_centralityComparison%s%s%s.%s", saveComment.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
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

    // Zoom for the y-axis
    double yZoom1 = 0.15;
    double yZoom2 = 30;
    if(weightExponent == 2){
      yZoom1 = 0.04;
      yZoom2 = 80;
    }

    for(int iTrackPt = firstDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt <= lastDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt++){
      compactTrackPtString = Form("_T>%.1f", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt));
      compactTrackPtString.ReplaceAll(".","v");

      // Draw a big canvas and put all the plots in it
      bigCanvas[iTrackPt] = new SplitCanvas(Form("bigCanvas%d", iTrackPt), "", 1800, 1500);
      bigCanvas[iTrackPt]->SetMargin(0.08, 0.01, 0.08, 0.11);
      bigCanvas[iTrackPt]->DivideNeatly(4, 4);

      bottomRowScale = bigCanvas[iTrackPt]->GetBottomRowScale();
      bottomPadMargin = bigCanvas[iTrackPt]->GetBottomPadMargin();
      leftPadMargin = bigCanvas[iTrackPt]->GetLeftPadMargin();

      for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
        for(int iJetPt = firstDrawnJetPtBinEEC[weightExponent-1]; iJetPt <= lastDrawnJetPtBinEEC[weightExponent-1]; iJetPt++) {
          jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt));

          canvasIndex = (iCentrality - firstDrawnCentralityBin[weightExponent-1]) * (lastDrawnJetPtBinEEC[weightExponent-1] - firstDrawnJetPtBinEEC[weightExponent-1] + 1) + (iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]);
          bigCanvas[iTrackPt]->CD(canvasIndex);

          gPad->SetLogy();
          gPad->SetLogx();

          // The titles in the bottow row need to be scaled, since it has more margin than other rows
          thisPadScale = (iCentrality == lastDrawnCentralityBin[weightExponent-1]) ? bottomRowScale : 1;

          // Set the axis titles and labels
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleOffset(0.85);
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleSize(0.2 * thisPadScale);
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelOffset(0.01);
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelSize(0.12 * thisPadScale);

          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleOffset(0.5 / thisPadScale);
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleSize(0.2 * thisPadScale);
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelOffset(0.01 / thisPadScale);
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelSize(0.12 * thisPadScale);

          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitle("EEC");
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetStats(0);
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetNdivisions(505);
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->CenterTitle();
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitle("#Deltar");
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetNdivisions(505);

          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetTitle("");

          // Set the drawing ranges
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(yZoom1, yZoom2);
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);
          systematicUncertaintyForPp[weightExponent-1][kCorrelatedUncertainty][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);

          // Set the drawing style for histograms
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);

          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePp);
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerSize(1.2);
          energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePp);
          energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->SetMarkerSize(1.2);

          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);

          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPp, 0.4);
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetLineColor(markerColorPp);
          energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->SetLineColor(markerColorPp);
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerColor(markerColorPp);
          energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->SetMarkerColor(markerColorPp);

          // Different drawing styles for correlated uncertainties
          systematicUncertaintyForPp[weightExponent-1][kCorrelatedUncertainty][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPp, 0.4);
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorUpPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorDownPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);

          // Draw the histograms
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Draw("e2");
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->Draw("same,e2");
          systematicUncertaintyForPp[weightExponent-1][kCorrelatedUncertainty][iJetPt][iTrackPt]->Draw("same,e3");
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->Draw("same,p");
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->Draw("same,p");

          // Create a legend for jet pT to each pad
          leftMarginAdder = (iJetPt == firstDrawnJetPtBinEEC[weightExponent-1]) ? leftPadMargin : 0;
          bottomMarginAdder = (iCentrality == lastDrawnCentralityBin[weightExponent-1]) ? bottomPadMargin : 0;
          //legend = new TLegend(0.05 + leftMarginAdder, 0.7 + bottomMarginAdder*0.5, 0.5 / (1 - leftMarginAdder), 0.8 / (thisPadScale*0.81) - bottomMarginAdder*1.4);
          //legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
          
          //legend->Draw();

          // Create a legend for collision systems to each pad
          legend = new TLegend(0.03 + leftMarginAdder, 0.05 + bottomMarginAdder, 0.5 / (1 - leftMarginAdder), 0.45 / thisPadScale);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
          legend->AddEntry(systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], Form("PbPb %.0f-%.0f%%", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality)), "p");
          legend->AddEntry(systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt], "pp", "p");
          legend->AddEntry((TObject*)0, jetPtString.Data(), "");
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
      if(iTrackPt == 0 || iTrackPt == 2 || iTrackPt == 4){
        mainTitle->DrawLatexNDC(0.7, 0.94, Form("p_{T}^{ch} > %.1f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt)));
        mainTitle->DrawLatexNDC(0.89, 0.94, Form("n=%d", weightExponent));
      } else {
        mainTitle->DrawLatexNDC(0.7, 0.94, Form("p_{T}^{ch} > %.0f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt)));
        mainTitle->DrawLatexNDC(0.87, 0.94, Form("n=%d", weightExponent));
      }

      // Save the figures to file
      if(saveFigures) {
        gPad->GetCanvas()->SaveAs(Form("figures/finalBigCanvas%s%s.%s", saveComment.Data(), compactTrackPtString.Data(), figureFormat.Data()));
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

    for(int iTrackPt = firstDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt <= lastDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt++){
      compactTrackPtString = Form("_T>%.1f", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt));
      compactTrackPtString.ReplaceAll(".","v");

      // Draw a big canvas and put all the plots in it
      bigRatioCanvas[iTrackPt] = new SplitCanvas(Form("bigRatioCanvas%d", iTrackPt), "", 1800, 1500);
      bigRatioCanvas[iTrackPt]->SetMargin(0.08, 0.01, 0.08, 0.11);
      bigRatioCanvas[iTrackPt]->DivideNeatly(4, 4);

      bottomRowScale = bigRatioCanvas[iTrackPt]->GetBottomRowScale();
      bottomPadMargin = bigRatioCanvas[iTrackPt]->GetBottomPadMargin();
      leftPadMargin = bigRatioCanvas[iTrackPt]->GetLeftPadMargin();

      for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
        for(int iJetPt = firstDrawnJetPtBinEEC[weightExponent-1]; iJetPt <= lastDrawnJetPtBinEEC[weightExponent-1]; iJetPt++) {
          jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt));

          canvasIndex = (iCentrality - firstDrawnCentralityBin[weightExponent-1]) * (lastDrawnJetPtBinEEC[weightExponent-1] - firstDrawnJetPtBinEEC[weightExponent-1] + 1) + (iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]);
          bigRatioCanvas[iTrackPt]->CD(canvasIndex);

          gPad->SetLogy(0); // Linear y-axis scale for ratios
          gPad->SetLogx();

          // The titles in the bottow row need to be scaled, since it has more margin than other rows
          thisPadScale = (iCentrality == lastDrawnCentralityBin[weightExponent-1]) ? bottomRowScale : 1;

          // Set the axis titles and labels
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleOffset(0.85);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleSize(0.2 * thisPadScale);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelOffset(0.01);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelSize(0.12 * thisPadScale);

          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleOffset(0.5 / thisPadScale);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleSize(0.2 * thisPadScale);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelOffset(0.01 / thisPadScale);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelSize(0.12 * thisPadScale);

          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitle("#frac{PbPb}{pp}");
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetStats(0);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetNdivisions(505);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitle("#Deltar");
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetNdivisions(505);

          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetTitle("");

          // Set the drawing ranges
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.2, 1.8);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);

          // Set the drawing style for histograms
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);

          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);

          // Different drawing styles for correlated uncertainties
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorUpPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorDownPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);

          // Draw the histograms
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Draw("e2");
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]->Draw("same,p");

          // Add the jet pT legend for each pad
          leftMarginAdder = (iJetPt == firstDrawnJetPtBinEEC[weightExponent-1]) ? leftPadMargin : 0;
          bottomMarginAdder = (iCentrality == lastDrawnCentralityBin[weightExponent-1]) ? bottomPadMargin : 0;
          legend = new TLegend(0.05 + leftMarginAdder, 0.7 + bottomMarginAdder*0.5, 0.5 / (1 - leftMarginAdder), 0.8 / (thisPadScale*0.81) - bottomMarginAdder*1.4);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
          legend->AddEntry((TObject*)0, jetPtString.Data(), "");
          legend->Draw();

          // Add the centrality legend for each pad
          legend = new TLegend(0.05 + leftMarginAdder, 0.05 + bottomMarginAdder, 0.5 / (1 - leftMarginAdder), 0.15 / (thisPadScale*0.6) + bottomMarginAdder*0.45);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
          legend->AddEntry(systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], Form("PbPb %.0f-%.0f%% / pp", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality)), "p,e2");
          legend->Draw();

          // Draw a line to one
          oneLine->Draw();

          // Draw the shifted pT to regular pT ratio
          if(drawShiftedPtRatio){
            energyEnergyCorrelatorShiftedPtRatio[weightExponent-1][iJetPt][iTrackPt]->SetLineStyle(2);
            energyEnergyCorrelatorShiftedPtRatio[weightExponent-1][iJetPt][iTrackPt]->SetLineColor(kBlack);
            energyEnergyCorrelatorShiftedPtRatio[weightExponent-1][iJetPt][iTrackPt]->Draw("same,C");
          }

          // Draw lines to 0.08 and 0.2 to
          if(drawVerticalLines){
            lineDrawer->DrawLine(0.08, 0.8-(weightExponent-1)*0.2, 0.08, 1.2);
            lineDrawer->DrawLine(0.2, 0.8-(weightExponent-1)*0.2, 0.2, 1.2);
          }

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
      if(iTrackPt == 0 || iTrackPt == 2 || iTrackPt == 4){
        mainTitle->DrawLatexNDC(0.7, 0.94, Form("p_{T}^{ch} > %.1f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt)));
        mainTitle->DrawLatexNDC(0.89, 0.94, Form("n=%d", weightExponent));
      } else {
        mainTitle->DrawLatexNDC(0.7, 0.94, Form("p_{T}^{ch} > %.0f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt)));
        mainTitle->DrawLatexNDC(0.87, 0.94, Form("n=%d", weightExponent));
      }

      // Save the figures to file
      if(saveFigures) {
        gPad->GetCanvas()->SaveAs(Form("figures/finalBigRatioCanvas%s%s.%s", saveComment.Data(), compactTrackPtString.Data(), figureFormat.Data()));
      }

    }  // Track pT loop
  } // If for drawing big canvases

  if(drawDoubleRatios){

    // For the double ratios, only draw one example bin
    doubleRatioCentralityBinIndex1 = card[kPbPb][weightExponent-1]->FindBinIndexCentrality(doubleRatioCentralityBin1);
    doubleRatioCentralityBinIndex2 = card[kPbPb][weightExponent-1]->FindBinIndexCentrality(doubleRatioCentralityBin2);
    doubleRatioJetPtBinIndex = card[kPbPb][weightExponent-1]->FindBinIndexJetPtEEC(doubleRatioJetPtBin);

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
      centralityString = Form("Cent: %.0f-%.0f%%", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality));
      compactCentralityString = Form("_C=%.0f-%.0f", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality));

      doubleRatioCanvas->CD(canvasIndex++);

      gPad->SetLogy(false);
      gPad->SetLogx();

      // Set the axis titles and labels
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetTitleOffset(1);
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetTitleSize(0.09);
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetLabelOffset(0.01);
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetLabelSize(0.07);

      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetTitleOffset(1.4);
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetTitleSize(0.08);
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetLabelOffset(0.01);
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetLabelSize(0.07);

      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetTitle(Form("#frac{PbPb/pp (p_{T}^{ch} > %.0f GeV)}{PbPb/pp (p_{T}^{ch} > %.0f GeV)}", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(trackPtBinsForDoubleRatio.second), card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(trackPtBinsForDoubleRatio.first)));
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->CenterTitle();
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->SetStats(0);
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetNdivisions(505);
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->CenterTitle();
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetTitle("#Deltar");
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->CenterTitle();
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetNdivisions(505);

      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->SetTitle("");

      // Set the drawing ranges
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetRangeUser(0.8, 1.2);

      // Set the drawing style for histograms
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->SetMarkerSize(1.2);
      energyEnergyCorrelatorDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
      energyEnergyCorrelatorDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->SetMarkerSize(1.2);

      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
      energyEnergyCorrelatorDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
      energyEnergyCorrelatorDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);

      // Draw the histograms
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->Draw("e2");
      energyEnergyCorrelatorDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->Draw("same,p");

      // Show the centrality bin in the legend
      leftMarginAdder = (iCentrality == doubleRatioCentralityBinIndex1) ? leftPadMargin-0.03 : 0;
      legend = new TLegend(0.05+leftMarginAdder, 0.3, 0.45 / (1-leftMarginAdder), 0.5);
      legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.07); legend->SetTextFont(62);
      legend->AddEntry((TObject*)0, Form("%.f < jet p_{T} < %.0f GeV",  card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(doubleRatioJetPtBinIndex), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(doubleRatioJetPtBinIndex)), "");
      legend->AddEntry(systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex], Form("Centrality: %.0f-%.0f%%", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality)), "p");

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
      gPad->GetCanvas()->SaveAs(Form("figures/finalDoubleRatio%s.%s", saveComment.Data(), figureFormat.Data()));
    }

  } // If for drawing double ratios

  if(drawDoubleRatioToSingleCanvas){

    TCanvas* theGreatCanvasOfDoubleRatio;

    double doubleRatioZoomMagnitude = 0.4;
    if(weightExponent == 1){
      if(trackPtBinsForDoubleRatio.first < 2){
        doubleRatioZoomMagnitude = 0.5;
      } else {
        doubleRatioZoomMagnitude = 0.3;
      }
    }

    TString doubleRatioJetPtString;
    TString doubleRatioTrackPtString;

    // For the double ratios, make separate canvases for each jet pT bin
    for(int iJetPt = firstDrawnJetPtBinEEC[weightExponent-1]; iJetPt <= lastDrawnJetPtBinEEC[weightExponent-1]; iJetPt++){

      // Create a TCanvas for the double ratio
      theGreatCanvasOfDoubleRatio = new SplitCanvas(Form("theGreatCanvasOfDoubleRatio%d", iJetPt), "", 1200, 800);
      theGreatCanvasOfDoubleRatio->SetMargin(0.24, 0.01, 0.2, 0.15);
      theGreatCanvasOfDoubleRatio->cd();

      gPad->SetLogy(false);
      gPad->SetLogx();

      mainTitle = new TLatex();
      canvasIndex = 0;

      // Set the drawing and label styles
      for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){

        // Set the axis titles and labels
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->GetXaxis()->SetTitleOffset(1);
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->GetXaxis()->SetTitleSize(0.09);
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->GetXaxis()->SetLabelOffset(0.01);
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->GetXaxis()->SetLabelSize(0.07);

        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->GetYaxis()->SetTitleOffset(1.2);
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->GetYaxis()->SetTitleSize(0.08);
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->GetYaxis()->SetLabelOffset(0.01);
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->GetYaxis()->SetLabelSize(0.07);

        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->GetYaxis()->SetTitle(Form("#frac{PbPb/pp (p_{T}^{ch} > %.0f GeV)}{PbPb/pp (p_{T}^{ch} > %.0f GeV)}", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(trackPtBinsForDoubleRatio.second), card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(trackPtBinsForDoubleRatio.first)));
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->GetXaxis()->CenterTitle();
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->SetStats(0);
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->GetYaxis()->SetNdivisions(505);
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->GetYaxis()->CenterTitle();
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->GetXaxis()->SetTitle("#Deltar");
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->GetXaxis()->CenterTitle();
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->GetXaxis()->SetNdivisions(505);

        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->SetTitle("");

        // Set the drawing ranges
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);

        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->GetYaxis()->SetRangeUser(1 - doubleRatioZoomMagnitude, 1 + doubleRatioZoomMagnitude);

        // Set the drawing styles 
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->SetMarkerSize(1.2);
        energyEnergyCorrelatorDoubleRatio[weightExponent-1][iCentrality][iJetPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
        energyEnergyCorrelatorDoubleRatio[weightExponent-1][iCentrality][iJetPt]->SetMarkerSize(1.2);

        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
        energyEnergyCorrelatorDoubleRatio[weightExponent-1][iCentrality][iJetPt]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
        energyEnergyCorrelatorDoubleRatio[weightExponent-1][iCentrality][iJetPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
      }

      // Draw the histograms
      systematicUncertaintyDoubleRatio[weightExponent-1][firstDrawnCentralityBin[weightExponent-1]][iJetPt]->Draw("e2");

      for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]+1; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->Draw("same,e2");
      }
      for(int iCentrality = lastDrawnCentralityBin[weightExponent-1]; iCentrality >= firstDrawnCentralityBin[weightExponent-1]; iCentrality--){
        energyEnergyCorrelatorDoubleRatio[weightExponent-1][iCentrality][iJetPt]->Draw("same,p");
      }

      // Show the centrality bin in the legend
      legend = new TLegend(0.28, 0.22, 0.53, 0.5);
      legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.055); legend->SetTextFont(62);
      for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
        legend->AddEntry(systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt], Form("Centrality: %.0f-%.0f%%", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality)), "p");
      }

      legend->Draw();

      // Draw the line to one
      oneLine->Draw();

      mainTitle->SetTextFont(62);
      mainTitle->SetTextSize(0.08);

      if(addPreliminaryTag){
        mainTitle->DrawLatexNDC(0.29, 0.93, "CMS");

        mainTitle->SetTextFont(52);
        mainTitle->SetTextSize(0.05);
        mainTitle->DrawLatexNDC(0.27, 0.88, "Preliminary");
      } else {
        mainTitle->DrawLatexNDC(0.3, 0.9, "CMS");
      }

      mainTitle->SetTextFont(42);
      mainTitle->SetTextSize(0.055);
      mainTitle->DrawLatexNDC(0.47, 0.945, "PbPb #sqrt{s_{NN}} = 5.02 TeV, 1.70 nb^{-1}");
      mainTitle->DrawLatexNDC(0.47, 0.87, "pp #sqrt{s} = 5.02 TeV, 302 pb^{-1}");
      mainTitle->DrawLatexNDC(0.3, 0.68, "anti-k_{T} R = 0.4,");
      mainTitle->DrawLatexNDC(0.54, 0.68, Form("n=%d", weightExponent));
      mainTitle->DrawLatexNDC(0.3, 0.77, Form("%.f < p_{T,jet} < %.0f GeV, |#eta_{jet}| < 1.6",  card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt)));

      // Save the figures to file
      if(saveFigures){
        doubleRatioJetPtBin = card[kPbPb][weightExponent-1]->GetBinBordersJetPtEEC(iJetPt);
        doubleRatioJetPtString = Form("_%.0f<jetpt<%.0f", doubleRatioJetPtBin.first, doubleRatioJetPtBin.second);
        doubleRatioTrackPtString = Form("_track%.0fover%.0f", trackPtCutsForDoubleDatio.second, trackPtCutsForDoubleDatio.first);
        gPad->GetCanvas()->SaveAs(Form("figures/finalDoubleRatioSingleCanvas%s%s%s.%s", saveComment.Data(), doubleRatioJetPtString.Data(), doubleRatioTrackPtString.Data(), figureFormat.Data()));
      }

    } // Jet pT loop

  } // If for drawing double ratios

  // ========================================================================================= //
  // ===== Drawing style with all energy weight exponent distributions in one big canvas ===== //
  // ========================================================================================= //
  
  if(drawBigCanvasAllDistributions){
    
    // Draw all the distributions to big canvases
    SplitCanvas* bigDualDistributionCanvas[nTrackPtBinsEEC];
    mainTitle = new TLatex();

    for(int iTrackPt = firstDrawnTrackPtBinEEC[0]; iTrackPt <= lastDrawnTrackPtBinEEC[0]; iTrackPt++){
      compactTrackPtString = Form("_T>%.1f", card[kPbPb][0]->GetLowBinBorderTrackPtEEC(iTrackPt));
      compactTrackPtString.ReplaceAll(".","v");

      // Draw a big canvas and put all the plots in it
      bigDualDistributionCanvas[iTrackPt] = new SplitCanvas(Form("bigDualDistributionCanvas%d", iTrackPt), "", 1800, 1500);
      bigDualDistributionCanvas[iTrackPt]->SetMargin(0.08, 0.01, 0.08, 0.11);
      bigDualDistributionCanvas[iTrackPt]->DivideNeatly(4, 4);

      bottomRowScale = bigDualDistributionCanvas[iTrackPt]->GetBottomRowScale();
      bottomPadMargin = bigDualDistributionCanvas[iTrackPt]->GetBottomPadMargin();
      leftPadMargin = bigDualDistributionCanvas[iTrackPt]->GetLeftPadMargin();

      for(int iCentrality = firstDrawnCentralityBin[0]; iCentrality <= lastDrawnCentralityBin[0]; iCentrality++){
        for(int iJetPt = firstDrawnJetPtBinEEC[0]; iJetPt <= lastDrawnJetPtBinEEC[0]; iJetPt++) {
          jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", card[kPbPb][0]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][0]->GetHighBinBorderJetPtEEC(iJetPt));

          canvasIndex = (iCentrality - firstDrawnCentralityBin[0]) * (lastDrawnJetPtBinEEC[0] - firstDrawnJetPtBinEEC[0] + 1) + (iJetPt - firstDrawnJetPtBinEEC[0]);
          bigDualDistributionCanvas[iTrackPt]->CD(canvasIndex);

          gPad->SetLogy();
          gPad->SetLogx();

          // The titles in the bottow row need to be scaled, since it has more margin than other rows
          thisPadScale = (iCentrality == lastDrawnCentralityBin[0]) ? bottomRowScale : 1;

          // Set the axis titles and labels
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleOffset(0.85);
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleSize(0.2 * thisPadScale);
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelOffset(0.01);
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelSize(0.12 * thisPadScale);

          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleOffset(0.5 / thisPadScale);
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleSize(0.2 * thisPadScale);
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelOffset(0.01 / thisPadScale);
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelSize(0.12 * thisPadScale);

          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitle("EEC");
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetStats(0);
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetNdivisions(505);
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->CenterTitle();
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitle("#Deltar");
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetNdivisions(505);

          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetTitle("");

          // Set the drawing ranges
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.04, 160);
          systematicUncertaintyForPbPb[0][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);
          systematicUncertaintyForPbPb[0][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);

          // Set the drawing style for histograms
          for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
            systematicUncertaintyForPbPb[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iWeightExponent]);
            systematicUncertaintyForPbPb[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);
            energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iWeightExponent]);
            energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);

            systematicUncertaintyForPbPb[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPbPb[iWeightExponent], 0.4);
            systematicUncertaintyForPbPb[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iWeightExponent]);
            energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iWeightExponent]);
            systematicUncertaintyForPbPb[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iWeightExponent]);
            energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iWeightExponent]);

            systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePpEnergyWeight[iWeightExponent]);
            systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerSize(1.2);
            energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePpEnergyWeight[iWeightExponent]);
            energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->SetMarkerSize(1.2);

            systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPpEnergyWeight[iWeightExponent], 0.4);
            systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetLineColor(markerColorPpEnergyWeight[iWeightExponent]);
            energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->SetLineColor(markerColorPpEnergyWeight[iWeightExponent]);
            systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerColor(markerColorPpEnergyWeight[iWeightExponent]);
            energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->SetMarkerColor(markerColorPpEnergyWeight[iWeightExponent]);          

            // Different drawing styles for correlated uncertainties
            systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorUpPbPb[iWeightExponent], 0.4);
            systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorDownPbPb[iWeightExponent], 0.4);
            systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
            systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);
            systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
            systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);

            systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPpEnergyWeight[iWeightExponent], 0.4);
          }

          // Draw the histograms
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Draw("e2");
          systematicUncertaintyForPbPb[0][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          systematicUncertaintyForPbPb[0][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          energyEnergyCorrelatorSignalPbPb[0][iCentrality][iJetPt][iTrackPt]->Draw("same,p");
          systematicUncertaintyForPbPb[1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Draw("same,e2");
          systematicUncertaintyForPbPb[1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          systematicUncertaintyForPbPb[1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          energyEnergyCorrelatorSignalPbPb[1][iCentrality][iJetPt][iTrackPt]->Draw("same,p");

          systematicUncertaintyForPp[0][kUncorrelatedUncertainty][iJetPt][iTrackPt]->Draw("same,e2");
          systematicUncertaintyForPp[0][kCorrelatedUncertainty][iJetPt][iTrackPt]->Draw("same,e3");
          energyEnergyCorrelatorSignalPp[0][iJetPt][iTrackPt]->Draw("same,p");
          systematicUncertaintyForPp[1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->Draw("same,e2");
          systematicUncertaintyForPp[1][kCorrelatedUncertainty][iJetPt][iTrackPt]->Draw("same,e3");
          energyEnergyCorrelatorSignalPp[1][iJetPt][iTrackPt]->Draw("same,p");

          // Legend adders to put the legends to correct positions
          leftMarginAdder = (iJetPt == firstDrawnJetPtBinEEC[0]) ? leftPadMargin : 0;
          bottomMarginAdder = (iCentrality == lastDrawnCentralityBin[0]) ? bottomPadMargin : 0;

          // Add the jet pT binning to the top row legend
          if(iCentrality == firstDrawnCentralityBin[0]){
            legend = new TLegend(0.05 + leftMarginAdder/1.1, 0.85, 0.5 / (1 - leftMarginAdder), 0.95);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry((TObject*)0, Form("%s", jetPtString.Data()), "");
            legend->Draw();

          }

          // Add the centrality bin information to the leftmost column
          if(iJetPt == firstDrawnJetPtBinEEC[0]){
            legend = new TLegend(0.19, 0.085 + bottomMarginAdder*1.05, 0.7, 0.2 / (thisPadScale*0.78) + bottomMarginAdder*0.1 + 0.05);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry((TObject*)0, Form("PbPb %.0f-%.0f%%", card[kPbPb][0]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][0]->GetHighBinBorderCentrality(iCentrality)), "");
            legend->Draw();
          }

          // Add information about the nominal weight PbPb histogram to selected pad
          if(iJetPt == firstDrawnJetPtBinEEC[0]+1 && iCentrality == firstDrawnCentralityBin[0]+1){
            legend = new TLegend(0.4, 0.83, 0.8, 0.93);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry(systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], "PbPb, n=1", "lpf");
            legend->Draw();
          }

          // Add information about the squared weight PbPb histogram to selected pad
          if(iJetPt == firstDrawnJetPtBinEEC[0]+2 && iCentrality == firstDrawnCentralityBin[0]+1){
            legend = new TLegend(0.4, 0.83, 0.8, 0.93);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry(systematicUncertaintyForPbPb[1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], "PbPb, n=2", "lpf");
            legend->Draw();
          }

          // Add information about the nominal weight PbPb histogram to selected pad
          if(iJetPt == firstDrawnJetPtBinEEC[0]+1 && iCentrality == firstDrawnCentralityBin[0]+2){
            legend = new TLegend(0.4, 0.83, 0.8, 0.93);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry(systematicUncertaintyForPp[0][kUncorrelatedUncertainty][iJetPt][iTrackPt], "pp, n=1", "lpf");
            legend->Draw();
          }

          // Add information about the squared weight PbPb histogram to selected pad
          if(iJetPt == firstDrawnJetPtBinEEC[0]+2 && iCentrality == firstDrawnCentralityBin[0]+2){
            legend = new TLegend(0.4, 0.83, 0.8, 0.93);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry(systematicUncertaintyForPp[1][kUncorrelatedUncertainty][iJetPt][iTrackPt], "pp, n=2", "lpf");
            legend->Draw();
          }

          // Create a legend for each pad
          //legend = new TLegend(0.05 + leftMarginAdder, 0.05 + bottomMarginAdder, 0.5 / (1 - leftMarginAdder), 0.2 / (thisPadScale*0.81));
          //legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
          //legend->AddEntry(systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], Form("#frac{PbPb %.0f-%.0f%%}{pp} #rightarrow", card[kPbPb][0]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][0]->GetHighBinBorderCentrality(iCentrality)), "p,e2");

          //legend->Draw();

          // Draw a line to one
          //oneLine->Draw();

          // Draw lines to 0.08 and 0.2 to
          if(drawVerticalLines){
            lineDrawer->DrawLine(0.08, 0.6, 0.08, 1.2);
            lineDrawer->DrawLine(0.2, 0.6, 0.2, 1.2);
          }

        }  // Jet pT loop
      }    // Centrality loop

      // Draw all necessary CMS text to the plot
      bigDualDistributionCanvas[iTrackPt]->cd(0);

      mainTitle->SetTextFont(62);

      if(addPreliminaryTag){
        mainTitle->SetTextSize(0.05);
        mainTitle->DrawLatexNDC(0.08, 0.945, "CMS");

        mainTitle->SetTextFont(52);
        mainTitle->SetTextSize(0.035);
        mainTitle->DrawLatexNDC(0.06, 0.91, "Preliminary");
      } else {
        mainTitle->SetTextSize(0.065);
        mainTitle->DrawLatexNDC(0.08, 0.93, "CMS");
      }

      mainTitle->SetTextFont(42);
      mainTitle->SetTextSize(0.025);
      mainTitle->DrawLatexNDC(0.23, 0.96, "PbPb #sqrt{s_{NN}} = 5.02 TeV, 1.70 nb^{-1}");
      mainTitle->DrawLatexNDC(0.23, 0.92, "pp #sqrt{s} = 5.02 TeV, 302 pb^{-1}");
      mainTitle->DrawLatexNDC(0.53, 0.96, "anti-k_{T} R = 0.4");
      mainTitle->DrawLatexNDC(0.53, 0.92, "|#eta_{jet}| < 1.6");

      mainTitle->SetTextSize(0.035);
      if(iTrackPt == 4){
        mainTitle->DrawLatexNDC(0.75, 0.94, Form("p_{T}^{ch} > %.1f GeV", card[kPbPb][0]->GetLowBinBorderTrackPtEEC(iTrackPt)));
      } else {
        mainTitle->DrawLatexNDC(0.75, 0.94, Form("p_{T}^{ch} > %.0f GeV", card[kPbPb][0]->GetLowBinBorderTrackPtEEC(iTrackPt)));
      }

      // Save the figures to file
      if(saveFigures) {
        gPad->GetCanvas()->SaveAs(Form("figures/finalBigDualDistributionCanvas%s%s.%s", saveComment.Data(), compactTrackPtString.Data(), figureFormat.Data()));
      }

    }  // Track pT loop
  } // If for drawing big canvases

  // ======================================================================================== //
  // ===== Drawing style with ratios from all energy weight exponents in one big canvas ===== //
  // ======================================================================================== //
  
  if(drawBigCanvasAllRatios){
    
    // Draw all the distributions to big canvases
    SplitCanvas* bigDualRatioCanvas[nTrackPtBinsEEC];
    mainTitle = new TLatex();

    for(int iTrackPt = firstDrawnTrackPtBinEEC[0]; iTrackPt <= lastDrawnTrackPtBinEEC[0]; iTrackPt++){
      compactTrackPtString = Form("_T>%.1f", card[kPbPb][0]->GetLowBinBorderTrackPtEEC(iTrackPt));
      compactTrackPtString.ReplaceAll(".","v");

      // Draw a big canvas and put all the plots in it
      bigDualRatioCanvas[iTrackPt] = new SplitCanvas(Form("bigDualRatioCanvas%d", iTrackPt), "", 1800, 1500);
      bigDualRatioCanvas[iTrackPt]->SetMargin(0.08, 0.01, 0.08, 0.11);
      bigDualRatioCanvas[iTrackPt]->DivideNeatly(4, 4);

      bottomRowScale = bigDualRatioCanvas[iTrackPt]->GetBottomRowScale();
      bottomPadMargin = bigDualRatioCanvas[iTrackPt]->GetBottomPadMargin();
      leftPadMargin = bigDualRatioCanvas[iTrackPt]->GetLeftPadMargin();

      for(int iCentrality = firstDrawnCentralityBin[0]; iCentrality <= lastDrawnCentralityBin[0]; iCentrality++){
        for(int iJetPt = firstDrawnJetPtBinEEC[0]; iJetPt <= lastDrawnJetPtBinEEC[0]; iJetPt++) {
          jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", card[kPbPb][0]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][0]->GetHighBinBorderJetPtEEC(iJetPt));

          canvasIndex = (iCentrality - firstDrawnCentralityBin[0]) * (lastDrawnJetPtBinEEC[0] - firstDrawnJetPtBinEEC[0] + 1) + (iJetPt - firstDrawnJetPtBinEEC[0]);
          bigDualRatioCanvas[iTrackPt]->CD(canvasIndex);

          gPad->SetLogy(0); // Linear y-axis scale for ratios
          gPad->SetLogx();

          // The titles in the bottow row need to be scaled, since it has more margin than other rows
          thisPadScale = (iCentrality == lastDrawnCentralityBin[0]) ? bottomRowScale : 1;

          // Set the axis titles and labels
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleOffset(0.85);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleSize(0.2 * thisPadScale);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelOffset(0.01);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelSize(0.12 * thisPadScale);

          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleOffset(0.5 / thisPadScale);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleSize(0.2 * thisPadScale);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelOffset(0.01 / thisPadScale);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelSize(0.12 * thisPadScale);

          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitle("#frac{PbPb}{pp}");
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetStats(0);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetNdivisions(505);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitle("#Deltar");
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetNdivisions(505);

          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetTitle("");

          // Set the drawing ranges
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.2, 1.8);
          systematicUncertaintyPbPbToPpRatio[0][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);
          systematicUncertaintyPbPbToPpRatio[0][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);

          // Set the drawing style for histograms
          for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iWeightExponent]);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);
            energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iWeightExponent]);
            energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);

            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPbPb[iWeightExponent], 0.4);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iWeightExponent]);
            energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iWeightExponent]);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iWeightExponent]);
            energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iWeightExponent]);

            // Different drawing styles for correlated uncertainties
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorUpPbPb[iWeightExponent], 0.4);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorDownPbPb[iWeightExponent], 0.4);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);
          }

          // Draw the histograms
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Draw("e2");
          systematicUncertaintyPbPbToPpRatio[0][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          systematicUncertaintyPbPbToPpRatio[0][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          energyEnergyCorrelatorPbPbToPpRatio[0][iCentrality][iJetPt][iTrackPt]->Draw("same,p");
          systematicUncertaintyPbPbToPpRatio[1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Draw("same,e2");
          systematicUncertaintyPbPbToPpRatio[1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          systematicUncertaintyPbPbToPpRatio[1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          energyEnergyCorrelatorPbPbToPpRatio[1][iCentrality][iJetPt][iTrackPt]->Draw("same,p");

          // Legend adders to put the legends to correct positions
          leftMarginAdder = (iJetPt == firstDrawnJetPtBinEEC[0]) ? leftPadMargin : 0;
          bottomMarginAdder = (iCentrality == lastDrawnCentralityBin[0]) ? bottomPadMargin : 0;

          // Add the jet pT binning to the top row legend
          if(iCentrality == firstDrawnCentralityBin[0]){
            legend = new TLegend(0.05 + leftMarginAdder/1.1, 0.83, 0.5 / (1 - leftMarginAdder), 0.93);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry((TObject*)0, Form("%s", jetPtString.Data()), "");
            legend->Draw();

          }

          // Add the centrality bin information to the leftmost column
          if(iJetPt == firstDrawnJetPtBinEEC[0]){
            legend = new TLegend(0.19, 0.05 + bottomMarginAdder*1.1, 0.7, 0.2 / (thisPadScale*0.78) + bottomMarginAdder*0.1);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry((TObject*)0, Form("#frac{PbPb %.0f-%.0f%%}{pp}", card[kPbPb][0]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][0]->GetHighBinBorderCentrality(iCentrality)), "");
            legend->Draw();
          }

          // Add information about the nominal weight histogram to selected oad
          if(iJetPt == firstDrawnJetPtBinEEC[0]+1 && iCentrality == firstDrawnCentralityBin[0]+1){
            legend = new TLegend(0.1, 0.83, 0.5, 0.93);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry(systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], "EEC ratio with n=1", "lpf");
            legend->Draw();
          }

          // Add information about the squared weight histogram to selected oad
          if(iJetPt == firstDrawnJetPtBinEEC[0]+2 && iCentrality == firstDrawnCentralityBin[0]+1){
            legend = new TLegend(0.1, 0.83, 0.5, 0.93);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry(systematicUncertaintyPbPbToPpRatio[1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], "EEC ratio with n=2", "lpf");
            legend->Draw();
          }

          // Create a legend for each pad
          //legend = new TLegend(0.05 + leftMarginAdder, 0.05 + bottomMarginAdder, 0.5 / (1 - leftMarginAdder), 0.2 / (thisPadScale*0.81));
          //legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
          //legend->AddEntry(systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], Form("#frac{PbPb %.0f-%.0f%%}{pp} #rightarrow", card[kPbPb][0]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][0]->GetHighBinBorderCentrality(iCentrality)), "p,e2");

          //legend->Draw();

          // Draw a line to one
          oneLine->Draw();

          // Draw lines to 0.08 and 0.2 to
          if(drawVerticalLines){
            lineDrawer->DrawLine(0.08, 0.6, 0.08, 1.2);
            lineDrawer->DrawLine(0.2, 0.6, 0.2, 1.2);
          }

        }  // Jet pT loop
      }    // Centrality loop

      // Draw all necessary CMS text to the plot
      bigDualRatioCanvas[iTrackPt]->cd(0);

      mainTitle->SetTextFont(62);

      if(addPreliminaryTag){
        mainTitle->SetTextSize(0.05);
        mainTitle->DrawLatexNDC(0.08, 0.945, "CMS");

        mainTitle->SetTextFont(52);
        mainTitle->SetTextSize(0.035);
        mainTitle->DrawLatexNDC(0.06, 0.91, "Preliminary");
      } else {
        mainTitle->SetTextSize(0.065);
        mainTitle->DrawLatexNDC(0.08, 0.93, "CMS");
      }

      mainTitle->SetTextFont(42);
      mainTitle->SetTextSize(0.025);
      mainTitle->DrawLatexNDC(0.23, 0.96, "PbPb #sqrt{s_{NN}} = 5.02 TeV, 1.70 nb^{-1}");
      mainTitle->DrawLatexNDC(0.23, 0.92, "pp #sqrt{s} = 5.02 TeV, 302 pb^{-1}");
      mainTitle->DrawLatexNDC(0.53, 0.96, "anti-k_{T} R = 0.4");
      mainTitle->DrawLatexNDC(0.53, 0.92, "|#eta_{jet}| < 1.6");

      mainTitle->SetTextSize(0.035);
      if(iTrackPt == 4){
        mainTitle->DrawLatexNDC(0.75, 0.94, Form("p_{T}^{ch} > %.1f GeV", card[kPbPb][0]->GetLowBinBorderTrackPtEEC(iTrackPt)));
      } else {
        mainTitle->DrawLatexNDC(0.75, 0.94, Form("p_{T}^{ch} > %.0f GeV", card[kPbPb][0]->GetLowBinBorderTrackPtEEC(iTrackPt)));
      }

      // Save the figures to file
      if(saveFigures) {
        gPad->GetCanvas()->SaveAs(Form("figures/finalBigDualRatioCanvas%s%s.%s", saveComment.Data(), compactTrackPtString.Data(), figureFormat.Data()));
      }

    }  // Track pT loop
  } // If for drawing big canvases

  // =================================================================================== //
  // ===== Drawing style for distributions to be shown in the letter format paper  ===== //
  // =================================================================================== //
  
  if(drawLetterPaperDistributions){

    // Only selected bins are shown in the letter paper format
    letterPaperCentralityBinIndex.clear();
    letterPaperJetPtBinIndex.clear();
    letterPaperTrackPtCutIndex.clear();

    letterPaperCentralityBinIndex.push_back(card[kPbPb][weightExponent-1]->FindBinIndexCentrality(letterPaperCentralityBin1));
    letterPaperCentralityBinIndex.push_back(card[kPbPb][weightExponent-1]->FindBinIndexCentrality(letterPaperCentralityBin2));
    letterPaperJetPtBinIndex.push_back(card[kPbPb][weightExponent-1]->FindBinIndexJetPtEEC(letterPaperJetPtBin));
    letterPaperTrackPtCutIndex.push_back(card[kPbPb][weightExponent-1]->GetBinIndexTrackPtEEC(letterPaperTrackPtCut1));
    
    // Define canvas and margins
    SplitCanvas* letterPaperDistributionCanvas;
    mainTitle = new TLatex();

    letterPaperDistributionCanvas = new SplitCanvas("letterPaperDistributionCanvas", "", 450, 800);
    letterPaperDistributionCanvas->SetMargin(0.22, 0.01, 0.11, 0.01);
    letterPaperDistributionCanvas->DivideNeatly(2, 1);

    bottomRowScale = letterPaperDistributionCanvas->GetBottomRowScale();
    bottomPadMargin = letterPaperDistributionCanvas->GetBottomPadMargin();
    leftPadMargin = letterPaperDistributionCanvas->GetLeftPadMargin();

    canvasIndex = 0;

    for(int iCentrality : letterPaperCentralityBinIndex){
      centralityString = Form("Cent: %.0f-%.0f%%", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality));
      for(int iJetPt : letterPaperJetPtBinIndex){
        jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][0]->GetHighBinBorderJetPtEEC(iJetPt));
        for(int iTrackPt : letterPaperTrackPtCutIndex){
          trackPtString = Form("p_{T}^{ch} > %.0f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt));
          
          letterPaperDistributionCanvas->CD(canvasIndex++);

          gPad->SetLogy();
          gPad->SetLogx();

          // The titles in the bottow row need to be scaled, since it has more margin than other rows
          thisPadScale = (canvasIndex > 1) ? bottomRowScale : 1;

          // Set the axis titles and labels
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleOffset(0.65);
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleSize(0.15 * thisPadScale);
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelOffset(0.0001);
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelSize(0.07 * thisPadScale);

          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleOffset(0.6 / thisPadScale);
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleSize(0.15 * thisPadScale);
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelOffset(0.01 / thisPadScale);
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelSize(0.07 * thisPadScale);

          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitle("EEC");
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetStats(0);
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetNdivisions(505);
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->CenterTitle();
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitle("#Deltar");
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetNdivisions(505);

          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetTitle("");

          // Set the drawing ranges
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.04, 160);
          systematicUncertaintyForPbPb[0][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);
          systematicUncertaintyForPbPb[0][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);

          // Set the drawing style for histograms
          for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
            systematicUncertaintyForPbPb[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iWeightExponent]);
            systematicUncertaintyForPbPb[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);
            energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iWeightExponent]);
            energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);

            systematicUncertaintyForPbPb[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPbPb[iWeightExponent], 0.4);
            systematicUncertaintyForPbPb[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iWeightExponent]);
            energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iWeightExponent]);
            systematicUncertaintyForPbPb[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iWeightExponent]);
            energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iWeightExponent]);

            systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePpEnergyWeight[iWeightExponent]);
            systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerSize(1.2);
            energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePpEnergyWeight[iWeightExponent]);
            energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->SetMarkerSize(1.2);

            systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPpEnergyWeight[iWeightExponent], 0.4);
            systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetLineColor(markerColorPpEnergyWeight[iWeightExponent]);
            energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->SetLineColor(markerColorPpEnergyWeight[iWeightExponent]);
            systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerColor(markerColorPpEnergyWeight[iWeightExponent]);
            energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->SetMarkerColor(markerColorPpEnergyWeight[iWeightExponent]);          

            // Different drawing styles for correlated uncertainties
            systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorUpPbPb[iWeightExponent], 0.4);
            systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorDownPbPb[iWeightExponent], 0.4);
            systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
            systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);
            systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
            systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);

            systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPpEnergyWeight[iWeightExponent], 0.4);
          }

          // Draw the histograms
          systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Draw("e2");
          systematicUncertaintyForPbPb[0][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          systematicUncertaintyForPbPb[0][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          energyEnergyCorrelatorSignalPbPb[0][iCentrality][iJetPt][iTrackPt]->Draw("same,p");
          systematicUncertaintyForPbPb[1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Draw("same,e2");
          systematicUncertaintyForPbPb[1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          systematicUncertaintyForPbPb[1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          energyEnergyCorrelatorSignalPbPb[1][iCentrality][iJetPt][iTrackPt]->Draw("same,p");

          systematicUncertaintyForPp[0][kUncorrelatedUncertainty][iJetPt][iTrackPt]->Draw("same,e2");
          systematicUncertaintyForPp[0][kCorrelatedUncertainty][iJetPt][iTrackPt]->Draw("same,e3");
          energyEnergyCorrelatorSignalPp[0][iJetPt][iTrackPt]->Draw("same,p");
          systematicUncertaintyForPp[1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->Draw("same,e2");
          systematicUncertaintyForPp[1][kCorrelatedUncertainty][iJetPt][iTrackPt]->Draw("same,e3");
          energyEnergyCorrelatorSignalPp[1][iJetPt][iTrackPt]->Draw("same,p");

          // Legend adders to put the legends to correct positions
          leftMarginAdder = (iJetPt == firstDrawnJetPtBinEEC[0]) ? leftPadMargin : 0;
          bottomMarginAdder = (iCentrality == lastDrawnCentralityBin[0]) ? bottomPadMargin : 0;

          // Upper of the two pads
          if(canvasIndex == 1){
            legend = new TLegend(0.5, 0.74, 0.95, 0.92);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.09 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry((TObject*)0, centralityString.Data(), "");
            legend->AddEntry((TObject*)0, trackPtString.Data(), "");
            legend->Draw();

          }

          // Lower of the two pads
          if(canvasIndex == 2){
            legend = new TLegend(0.47, 0.79, 0.92, 0.94);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.09 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry((TObject*)0, centralityString.Data(), "");
            legend->AddEntry((TObject*)0, trackPtString.Data(), "");
            legend->Draw();

            anotherLegend = new TLegend(0.27, 0.24, 0.57, 0.64);
            anotherLegend->SetFillStyle(0); anotherLegend->SetBorderSize(0); anotherLegend->SetTextSize(0.065 * thisPadScale); anotherLegend->SetTextFont(62);
            anotherLegend->AddEntry(systematicUncertaintyForPbPb[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], "PbPb, n=1", "lpf");
            anotherLegend->AddEntry(systematicUncertaintyForPbPb[1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], "PbPb, n=2", "lpf");
            anotherLegend->AddEntry(systematicUncertaintyForPp[0][kUncorrelatedUncertainty][iJetPt][iTrackPt], "pp, n=1", "lpf");
            anotherLegend->AddEntry(systematicUncertaintyForPp[1][kUncorrelatedUncertainty][iJetPt][iTrackPt], "pp, n=2", "lpf");
            anotherLegend->Draw();
          }

        } // Track pT loop
      }  // Jet pT loop
    }  // Centrality loop

    // Draw all necessary CMS text to the plot
    letterPaperDistributionCanvas->cd(0);

    mainTitle->SetTextFont(62);
    mainTitle->SetTextSize(0.08);
    mainTitle->DrawLatexNDC(0.26, 0.94, "CMS");

    mainTitle->SetTextFont(42);
    mainTitle->SetTextSize(0.05);
    mainTitle->DrawLatexNDC(0.265, 0.575, "PbPb #sqrt{s_{NN}} = 5.02 TeV, 1.70 nb^{-1}");
    mainTitle->DrawLatexNDC(0.265, 0.625, "pp #sqrt{s} = 5.02 TeV, 302 pb^{-1}");
    mainTitle->DrawLatexNDC(0.265, 0.725, "anti-k_{T} R = 0.4");
    mainTitle->DrawLatexNDC(0.265, 0.775, "|#eta_{jet}| < 1.6");
    //mainTitle->DrawLatexNDC(0.265, 0.775, Form("p_{T}^{ch} > %.0f GeV", letterPaperTrackPtCut1));
    mainTitle->DrawLatexNDC(0.265, 0.675, jetPtString.Data());

    // Save the figures to file
    if(saveFigures) {
      gPad->GetCanvas()->SaveAs(Form("figures/finalLetterPaperDistribution%s.%s", saveComment.Data(), figureFormat.Data()));
    }

  } // If for drawing distributions for the letter paper format

  // ================================================================ //
  // ===== Drawing style for ratio plots in letter paper format ===== //
  // ================================================================ //
  
  if(drawLetterPaperRatios){

    // Only selected bins are shown in the letter paper format
    letterPaperCentralityBinIndex.clear();
    letterPaperJetPtBinIndex.clear();
    letterPaperTrackPtCutIndex.clear();

    letterPaperCentralityBinIndex.push_back(card[kPbPb][weightExponent-1]->FindBinIndexCentrality(letterPaperCentralityBin1));
    letterPaperCentralityBinIndex.push_back(card[kPbPb][weightExponent-1]->FindBinIndexCentrality(letterPaperCentralityBin2));
    letterPaperJetPtBinIndex.push_back(card[kPbPb][weightExponent-1]->FindBinIndexJetPtEEC(letterPaperJetPtBin));
    letterPaperTrackPtCutIndex.push_back(card[kPbPb][weightExponent-1]->GetBinIndexTrackPtEEC(letterPaperTrackPtCut1));
    letterPaperTrackPtCutIndex.push_back(card[kPbPb][weightExponent-1]->GetBinIndexTrackPtEEC(letterPaperTrackPtCut2));
    
    // Define canvas and margins
    SplitCanvas* letterPaperRatioCanvas;
    mainTitle = new TLatex();

    letterPaperRatioCanvas = new SplitCanvas("letterPaperRatioCanvas", "", 1000, 1000);
    letterPaperRatioCanvas->SetMargin(0.13, 0.01, 0.13, 0.01);
    letterPaperRatioCanvas->DivideNeatly(2,2);

    bottomRowScale = letterPaperRatioCanvas->GetBottomRowScale();
    bottomPadMargin = letterPaperRatioCanvas->GetBottomPadMargin();
    leftPadMargin = letterPaperRatioCanvas->GetLeftPadMargin();

    canvasIndex = 0;

    for(int iCentrality : letterPaperCentralityBinIndex){
      centralityString = Form("Cent: %.0f-%.0f%%", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality));
      for(int iJetPt : letterPaperJetPtBinIndex){
        jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][0]->GetHighBinBorderJetPtEEC(iJetPt));
        for(int iTrackPt : letterPaperTrackPtCutIndex){
          trackPtString = Form("p_{T}^{ch} > %.0f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt));
          
          letterPaperRatioCanvas->CD(canvasIndex++);

          gPad->SetLogy(0); // Linear y-axis scale for ratios
          gPad->SetLogx();

          // The titles in the bottow row need to be scaled, since it has more margin than other rows
          thisPadScale = (canvasIndex == 3) ? bottomRowScale : 1;


          // Set the axis titles and labels
          if(canvasIndex == 4){
            systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleOffset(0.66);
            systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelOffset(-0.027);
          } else {
            systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleOffset(0.88);
            systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelOffset(0.001);
          }
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleSize(0.13 * thisPadScale);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelSize(0.1 * thisPadScale);

          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleOffset(0.63 / thisPadScale);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleSize(0.13 * thisPadScale);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelOffset(0.01 / thisPadScale);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelSize(0.1 * thisPadScale);

          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitle("#frac{PbPb}{pp}");
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetStats(0);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetNdivisions(505);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitle("#Deltar");
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetNdivisions(505);

          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetTitle("");

          // Set the drawing ranges
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.2, 1.8);
          systematicUncertaintyPbPbToPpRatio[0][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);
          systematicUncertaintyPbPbToPpRatio[0][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);

          // Set the drawing style for histograms
          for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iWeightExponent]);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);
            energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iWeightExponent]);
            energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);

            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPbPb[iWeightExponent], 0.4);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iWeightExponent]);
            energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iWeightExponent]);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iWeightExponent]);
            energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iWeightExponent]);

            // Different drawing styles for correlated uncertainties
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorUpPbPb[iWeightExponent], 0.4);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorDownPbPb[iWeightExponent], 0.4);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);
          }

          // Draw the histograms
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Draw("e2");
          systematicUncertaintyPbPbToPpRatio[0][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          systematicUncertaintyPbPbToPpRatio[0][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          energyEnergyCorrelatorPbPbToPpRatio[0][iCentrality][iJetPt][iTrackPt]->Draw("same,p");
          systematicUncertaintyPbPbToPpRatio[1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Draw("same,e2");
          systematicUncertaintyPbPbToPpRatio[1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          systematicUncertaintyPbPbToPpRatio[1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          energyEnergyCorrelatorPbPbToPpRatio[1][iCentrality][iJetPt][iTrackPt]->Draw("same,p");

          // Legend adders to put the legends to correct positions
          leftMarginAdder = (canvasIndex % 2 == 1) ? leftPadMargin : 0;
          bottomMarginAdder = (canvasIndex > 2) ? bottomPadMargin : 0;

          // Add centrality and track pT cut information to each panel
          if(canvasIndex == 1){
            legend = new TLegend(0.5, 0.75, 0.8, 0.95);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.08 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry((TObject*)0, centralityString.Data(), "");
            legend->AddEntry((TObject*)0, trackPtString.Data(), "");
            legend->Draw();

            anotherLegend = new TLegend(0.28, 0.05, 0.88, 0.15);
            anotherLegend->SetFillStyle(0); anotherLegend->SetBorderSize(0); anotherLegend->SetTextSize(0.08 * thisPadScale); anotherLegend->SetTextFont(62);
            anotherLegend->AddEntry(systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], "EEC ratio with n=1", "lpf");
            anotherLegend->Draw();
          }

          if(canvasIndex == 2){
            legend = new TLegend(0.4, 0.75, 0.7, 0.95);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.08 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry((TObject*)0, centralityString.Data(), "");
            legend->AddEntry((TObject*)0, trackPtString.Data(), "");
            legend->Draw();

            anotherLegend = new TLegend(0.08, 0.05, 0.83, 0.15);
            anotherLegend->SetFillStyle(0); anotherLegend->SetBorderSize(0); anotherLegend->SetTextSize(0.08 * thisPadScale); anotherLegend->SetTextFont(62);
            anotherLegend->AddEntry(systematicUncertaintyPbPbToPpRatio[1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], "EEC ratio with n=2", "lpf");
            anotherLegend->Draw();
          }

          if(canvasIndex == 3){
            legend = new TLegend(0.5, 0.8, 0.8, 0.95);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.08 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry((TObject*)0, centralityString.Data(), "");
            legend->AddEntry((TObject*)0, trackPtString.Data(), "");
            legend->Draw();
          }

          if(canvasIndex == 4){
            legend = new TLegend(0.4, 0.8, 0.7, 0.95);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.08 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry((TObject*)0, centralityString.Data(), "");
            legend->AddEntry((TObject*)0, trackPtString.Data(), "");
            legend->Draw();
          }

          // Draw a line to one
          oneLine->Draw();

        } // Track pT loop
      }  // Jet pT loop
    }  // Centrality loop

    // Draw all necessary CMS text to the plot
    letterPaperRatioCanvas->cd(0);

    mainTitle->SetTextFont(62);
    mainTitle->SetTextSize(0.062);
    mainTitle->DrawLatexNDC(0.16, 0.69, "CMS");

    mainTitle->SetTextFont(42);
    mainTitle->SetTextSize(0.03);
    mainTitle->DrawLatexNDC(0.16, 0.22, "PbPb #sqrt{s_{NN}} = 5.02 TeV, 1.70 nb^{-1}");
    mainTitle->DrawLatexNDC(0.16, 0.165, "pp #sqrt{s} = 5.02 TeV, 302 pb^{-1}");
    mainTitle->DrawLatexNDC(0.6, 0.165, "|#eta_{jet}| < 1.6, anti-k_{T} R = 0.4");
    mainTitle->DrawLatexNDC(0.6, 0.22, jetPtString.Data());

    // Save the figures to file
    if(saveFigures) {
      gPad->GetCanvas()->SaveAs(Form("figures/finalLetterPaperRatio%s.%s", saveComment.Data(), figureFormat.Data()));
    }

  } // If for ratio plots for letter paper

  // Axis drawing ranges for supplementary distributions
  double supplementaryDistributionZoomLow[2] = {0.4, 0.04};
  double supplementaryDistributionZoomHigh[2] = {25, 60};

  // Draw supplementary plots where all jet pT bins are drawn into a single canvas
  if(drawSupplementaryEECJetPt){

    TCanvas* theGreatCanvasOfJetPt;
    int firstJetPtIndex = firstDrawnJetPtBinEEC[weightExponent-1];

    // Track pT loop
    for(int iTrackPt = firstDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt <= lastDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt++){

      trackPtString = Form("p_{T}^{ch} > %.0f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt));
      compactTrackPtString = Form("_T>%.1f", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt));
      compactTrackPtString.ReplaceAll(".","v");

      // ===================================== //
      // Draw the pp energy-energy correlators //
      // ===================================== //

      // Create a TCanvas for each centrality+track pT bin
      theGreatCanvasOfJetPt = new SplitCanvas(Form("theGreatCanvasOfJetPtPp%d", iTrackPt), "", 1100, 1100);
      theGreatCanvasOfJetPt->SetMargin(0.15, 0.01, 0.15, 0.07);
      theGreatCanvasOfJetPt->cd();

      gPad->SetLogy();
      gPad->SetLogx();

      mainTitle = new TLatex();

      // Set the axis titles and labels
      systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][firstJetPtIndex][iTrackPt]->GetXaxis()->SetTitleOffset(0.85);
      systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][firstJetPtIndex][iTrackPt]->GetXaxis()->SetTitleSize(0.08);
      systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][firstJetPtIndex][iTrackPt]->GetXaxis()->SetLabelOffset(0.001);
      systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][firstJetPtIndex][iTrackPt]->GetXaxis()->SetLabelSize(0.07);

      systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][firstJetPtIndex][iTrackPt]->GetYaxis()->SetTitleOffset(0.9);
      systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][firstJetPtIndex][iTrackPt]->GetYaxis()->SetTitleSize(0.08);
      systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][firstJetPtIndex][iTrackPt]->GetYaxis()->SetLabelOffset(0.01);
      systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][firstJetPtIndex][iTrackPt]->GetYaxis()->SetLabelSize(0.07);

      systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][firstJetPtIndex][iTrackPt]->GetYaxis()->SetTitle("EEC");
      systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][firstJetPtIndex][iTrackPt]->GetXaxis()->CenterTitle();
      systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][firstJetPtIndex][iTrackPt]->SetStats(0);
      systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][firstJetPtIndex][iTrackPt]->GetYaxis()->SetNdivisions(505);
      systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][firstJetPtIndex][iTrackPt]->GetYaxis()->CenterTitle();
      systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][firstJetPtIndex][iTrackPt]->GetXaxis()->SetTitle("#Deltar");
      systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][firstJetPtIndex][iTrackPt]->GetXaxis()->CenterTitle();
      systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][firstJetPtIndex][iTrackPt]->GetXaxis()->SetNdivisions(505);

      systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][firstJetPtIndex][iTrackPt]->SetTitle("");

      // Set the drawing ranges
      for(int iJetPt = firstDrawnJetPtBinEEC[weightExponent-1]; iJetPt <= lastDrawnJetPtBinEEC[weightExponent-1]; iJetPt++){
        systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
        systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(supplementaryDistributionZoomLow[weightExponent-1], supplementaryDistributionZoomHigh[weightExponent-1]);
        systematicUncertaintyForPp[weightExponent-1][kCorrelatedUncertainty][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);

        // Set the drawing style for histograms
        systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);
        systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerSize(1.2);
        energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);
        energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->SetMarkerSize(1.2);

        systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]], 0.4);
        systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);
        energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);
        systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);
        energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);

        // Different drawing styles for correlated uncertainties
        if(iJetPt % 2 == 0){
          systematicUncertaintyForPp[weightExponent-1][kCorrelatedUncertainty][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorDownPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]], 0.4);
        } else {
          systematicUncertaintyForPp[weightExponent-1][kCorrelatedUncertainty][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorUpPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]], 0.4);
        }
        systematicUncertaintyForPp[weightExponent-1][kCorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerSize(0);
        systematicUncertaintyForPp[weightExponent-1][kCorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerStyle(9);

        // Draw the systematic uncertainties
        if(iJetPt == firstDrawnJetPtBinEEC[weightExponent-1]){
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->Draw("e2");
        } else {
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->Draw("same,e2");
        }

        systematicUncertaintyForPp[weightExponent-1][kCorrelatedUncertainty][iJetPt][iTrackPt]->Draw("same,e3");
      }

      // Make a legend for the jet pT bins
      legend = new TLegend(0.18, 0.18, 0.43, 0.46);
      legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(supplementaryLegendTextSize); legend->SetTextFont(62);

      // Draw the data points on top of the error bands
      for(int iJetPt = firstDrawnJetPtBinEEC[weightExponent-1]; iJetPt <= lastDrawnJetPtBinEEC[weightExponent-1]; iJetPt++){
        energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->Draw("same,p");
        legend->AddEntry(systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt], Form("%.0f < p_{T,jet} < %.0f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt)), "lpf");
      }
 
      legend->Draw();

      mainTitle->SetTextFont(62);
      mainTitle->SetTextSize(0.055);
      mainTitle->DrawLatexNDC(0.15, 0.95, "CMS");

      if(addPreliminaryTag){
        mainTitle->SetTextFont(52);
        mainTitle->SetTextSize(0.04);
        mainTitle->DrawLatexNDC(0.27, 0.95, "Preliminary");
      }

      mainTitle->SetTextFont(42);
      mainTitle->SetTextSize(0.04);
      mainTitle->DrawLatexNDC(0.626, 0.95, "302 pb^{-1} pp (5.02 TeV)");

      mainTitle->SetTextFont(62);
      mainTitle->SetTextSize(supplementaryLegendTextSize);
      mainTitle->DrawLatexNDC(0.8, 0.85, Form("pp, n=%d", weightExponent));
      mainTitle->DrawLatexNDC(0.726, 0.78, trackPtString.Data());

      mainTitle->DrawLatexNDC(0.19, 0.56, "anti-k_{T} R = 0.4");
      mainTitle->DrawLatexNDC(0.19, 0.49, "|#eta_{jet}| < 1.6");

      // Save the figures to file
      if(saveFigures){
         gPad->GetCanvas()->SaveAs(Form("figures/finalEECAllJetPt%s_pp%s.%s", saveComment.Data(), compactTrackPtString.Data(), figureFormat.Data()));
      }

      // ======================================= //
      // Draw the PbPb energy-energy correlators //
      // ======================================= //
    
      // Centrality loop
      for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
        
        centralityString = Form("PbPb %.0f-%.0f%%", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality));
        compactCentralityString = Form("_C=%.0f-%.0f", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality));
      

        // Create a TCanvas for each centrality+track pT bin
        theGreatCanvasOfJetPt = new SplitCanvas(Form("theGreatCanvasOfJetPt%d%d", iCentrality, iTrackPt), "", 1100, 1100);
        theGreatCanvasOfJetPt->SetMargin(0.15, 0.01, 0.15, 0.07);
        theGreatCanvasOfJetPt->cd();

        gPad->SetLogy();
        gPad->SetLogx();

        mainTitle = new TLatex();

        // Set the axis titles and labels
        systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][firstJetPtIndex][iTrackPt]->GetXaxis()->SetTitleOffset(0.85);
        systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][firstJetPtIndex][iTrackPt]->GetXaxis()->SetTitleSize(0.08);
        systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][firstJetPtIndex][iTrackPt]->GetXaxis()->SetLabelOffset(0.001);
        systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][firstJetPtIndex][iTrackPt]->GetXaxis()->SetLabelSize(0.07);

        systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][firstJetPtIndex][iTrackPt]->GetYaxis()->SetTitleOffset(0.9);
        systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][firstJetPtIndex][iTrackPt]->GetYaxis()->SetTitleSize(0.08);
        systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][firstJetPtIndex][iTrackPt]->GetYaxis()->SetLabelOffset(0.01);
        systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][firstJetPtIndex][iTrackPt]->GetYaxis()->SetLabelSize(0.07);

        systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][firstJetPtIndex][iTrackPt]->GetYaxis()->SetTitle("EEC");
        systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][firstJetPtIndex][iTrackPt]->GetXaxis()->CenterTitle();
        systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][firstJetPtIndex][iTrackPt]->SetStats(0);
        systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][firstJetPtIndex][iTrackPt]->GetYaxis()->SetNdivisions(505);
        systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][firstJetPtIndex][iTrackPt]->GetYaxis()->CenterTitle();
        systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][firstJetPtIndex][iTrackPt]->GetXaxis()->SetTitle("#Deltar");
        systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][firstJetPtIndex][iTrackPt]->GetXaxis()->CenterTitle();
        systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][firstJetPtIndex][iTrackPt]->GetXaxis()->SetNdivisions(505);

        systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][firstJetPtIndex][iTrackPt]->SetTitle("");

        // Set the drawing ranges
        for(int iJetPt = firstDrawnJetPtBinEEC[weightExponent-1]; iJetPt <= lastDrawnJetPtBinEEC[weightExponent-1]; iJetPt++){
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(supplementaryDistributionZoomLow[weightExponent-1], supplementaryDistributionZoomHigh[weightExponent-1]);
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);

          // Set the drawing style for histograms
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);


          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]], 0.4);
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);

          // Different drawing styles for correlated uncertainties
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorUpPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]], 0.4);
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorDownPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]], 0.4);
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);

          // Draw the systematic uncertainties
          if(iJetPt == firstDrawnJetPtBinEEC[weightExponent-1]){
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Draw("e2");
          } else {
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Draw("same,e2");
          }

          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
        }

        // Make a legend for the jet pT bins
        legend = new TLegend(0.18, 0.18, 0.43, 0.46);
        legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(supplementaryLegendTextSize); legend->SetTextFont(62);

        // Draw the data points on top of the error bands
        for(int iJetPt = firstDrawnJetPtBinEEC[weightExponent-1]; iJetPt <= lastDrawnJetPtBinEEC[weightExponent-1]; iJetPt++){
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->Draw("same,p");
          legend->AddEntry(systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], Form("%.0f < p_{T,jet} < %.0f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt)), "lpf");
        }
          

        legend->Draw();

        mainTitle->SetTextFont(62);
        mainTitle->SetTextSize(0.055);
        mainTitle->DrawLatexNDC(0.15, 0.95, "CMS");

        if(addPreliminaryTag){
          mainTitle->SetTextFont(52);
          mainTitle->SetTextSize(0.04);
          mainTitle->DrawLatexNDC(0.27, 0.95, "Preliminary");
        }

        mainTitle->SetTextFont(42);
        mainTitle->SetTextSize(0.04);
        mainTitle->DrawLatexNDC(0.568, 0.95, "1.70 nb^{-1} PbPb (5.02 TeV)");

        mainTitle->SetTextFont(62);
        mainTitle->SetTextSize(supplementaryLegendTextSize);
        if(iCentrality == 0){
          mainTitle->DrawLatexNDC(0.618, 0.85, Form("%s, n=%d", centralityString.Data(), weightExponent));
        } else {
          mainTitle->DrawLatexNDC(0.594, 0.85, Form("%s, n=%d", centralityString.Data(), weightExponent));
        }
        mainTitle->DrawLatexNDC(0.726, 0.78, trackPtString.Data());

        mainTitle->DrawLatexNDC(0.19, 0.56, "anti-k_{T} R = 0.4");
        mainTitle->DrawLatexNDC(0.19, 0.49, "|#eta_{jet}| < 1.6");


        // Save the figures to file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/finalEECAllJetPt%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
        }
      } // Track pT loop
    } // Centrality loop
  } // If for drawing all jet pT bins in one plot

  // Draw supplementary plots where all centrality bins are drawn into a single canvas
  if(drawSupplementaryEECCentrality){

    TCanvas* theGreatCanvasOfCentrality;
    const int nDrawingStyles = 7;
    TString drawingStyleString[nDrawingStyles] = {"_onlyC0", "_onlyC1", "_onlyC2", "_onlyC3", "_onlyPp", "_ppAndC0", ""};

    // Track pT loop
    for(int iTrackPt = firstDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt <= lastDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt++){

      trackPtString = Form("p_{T}^{ch} > %.0f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt));
      compactTrackPtString = Form("_T>%.1f", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt));
      compactTrackPtString.ReplaceAll(".","v");

      // Jet pT loop
      for(int iJetPt = firstDrawnJetPtBinEEC[weightExponent-1]; iJetPt <= lastDrawnJetPtBinEEC[weightExponent-1]; iJetPt++){

        jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt));
        compactJetPtString = Form("_J=%.0f-%.0f", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt));

        // Draw a canvas for each component separately for conference purposes
        for(int iDrawingStyle = 0; iDrawingStyle < nDrawingStyles; iDrawingStyle++){

          // Create a TCanvas for each centrality+track pT bin
          theGreatCanvasOfCentrality = new SplitCanvas(Form("theGreatCanvasOfCentrality%d%d%d", iTrackPt, iJetPt, iDrawingStyle), "", 1100, 1100);
          theGreatCanvasOfCentrality->SetMargin(0.15, 0.01, 0.15, 0.07);
          theGreatCanvasOfCentrality->cd();

          gPad->SetLogy();
          gPad->SetLogx();

          mainTitle = new TLatex();

          // Set the axis titles and labels
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetXaxis()->SetTitleOffset(0.85);
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetXaxis()->SetTitleSize(0.08);
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetXaxis()->SetLabelOffset(0.001);
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetXaxis()->SetLabelSize(0.07);

          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetYaxis()->SetTitleOffset(0.9);
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetYaxis()->SetTitleSize(0.08);
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetYaxis()->SetLabelOffset(0.01);
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetYaxis()->SetLabelSize(0.07);

          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetYaxis()->SetTitle("EEC");
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetStats(0);
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetYaxis()->SetNdivisions(505);
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetYaxis()->CenterTitle();
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetXaxis()->SetTitle("#Deltar");
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetXaxis()->SetNdivisions(505);

          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetTitle("");

          // Do the same for PbPb histograms
          for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleOffset(0.85);
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleSize(0.08);
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelOffset(0.001);
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelSize(0.07);

            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleOffset(0.9);
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleSize(0.08);
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelOffset(0.01);
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelSize(0.07);

            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitle("EEC");
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetStats(0);
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetNdivisions(505);
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->CenterTitle();
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitle("#Deltar");
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetNdivisions(505);

            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetTitle("");
          }

          // Set the drawing ranges
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(supplementaryDistributionZoomLow[weightExponent-1], supplementaryDistributionZoomHigh[weightExponent-1]);
          systematicUncertaintyForPp[weightExponent-1][kCorrelatedUncertainty][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);

          // Set the drawing style for pp
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePp);
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerSize(1.2);
          energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePp);
          energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->SetMarkerSize(1.2);

          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPp, 0.4);
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetLineColor(markerColorPp);
          energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->SetLineColor(markerColorPp);
          systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerColor(markerColorPp);
          energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->SetMarkerColor(markerColorPp);

          systematicUncertaintyForPp[weightExponent-1][kCorrelatedUncertainty][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPp, 0.4);
          systematicUncertaintyForPp[weightExponent-1][kCorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerSize(0);
          systematicUncertaintyForPp[weightExponent-1][kCorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerStyle(9);
      

          // Set the style also for PbPb plots
          for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){

            // Set drawing ranges
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(supplementaryDistributionZoomLow[weightExponent-1], supplementaryDistributionZoomHigh[weightExponent-1]);
            systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);
            systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);

            // Set the drawing style for the PbPb histograms
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);
            energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
            energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);

            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
            energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
            energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);

            // Different drawing styles for correlated uncertainties
            systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorUpPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
            systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorDownPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
            systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
            systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);
            systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
            systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);
          }

          // First, draw all uncorrelated systematic uncertainties
          if(iDrawingStyle == 6){
            systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->Draw("e2");
            for(int iCentrality = lastDrawnCentralityBin[weightExponent-1]; iCentrality >= firstDrawnCentralityBin[weightExponent-1]; iCentrality--){
              systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Draw("same,e2");
            }
          } else if (iDrawingStyle == 4){
            systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->Draw("e2");
          } else if (iDrawingStyle == 5){
            systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->Draw("e2");
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][0][iJetPt][iTrackPt]->Draw("same,e2");
          } else {
            systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iDrawingStyle][iJetPt][iTrackPt]->Draw("e2");
          }

          // Then, draw all correlated uncertainties
          if(iDrawingStyle == 6){
            systematicUncertaintyForPp[weightExponent-1][kCorrelatedUncertainty][iJetPt][iTrackPt]->Draw("same,e3");
            for(int iCentrality = lastDrawnCentralityBin[weightExponent-1]; iCentrality >= firstDrawnCentralityBin[weightExponent-1]; iCentrality--){
              systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
              systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
            }
          } else if(iDrawingStyle == 4){
            systematicUncertaintyForPp[weightExponent-1][kCorrelatedUncertainty][iJetPt][iTrackPt]->Draw("same,e3");
          } else if(iDrawingStyle == 5){
            systematicUncertaintyForPp[weightExponent-1][kCorrelatedUncertainty][iJetPt][iTrackPt]->Draw("same,e3");
            systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeUp][0][iJetPt][iTrackPt]->Draw("same,e3");
            systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeDown][0][iJetPt][iTrackPt]->Draw("same,e3");
          } else {
            systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeUp][iDrawingStyle][iJetPt][iTrackPt]->Draw("same,e3");
            systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeDown][iDrawingStyle][iJetPt][iTrackPt]->Draw("same,e3");
          }

          // Finally, draw the data points themselves
          if(iDrawingStyle == 6){
            for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
              energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->Draw("same,p");
            }
            energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->Draw("same,p");
          } else if(iDrawingStyle == 4){
            energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->Draw("same,p");
          } else if(iDrawingStyle == 5){
            energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->Draw("same,p");
            energyEnergyCorrelatorSignalPbPb[weightExponent-1][0][iJetPt][iTrackPt]->Draw("same,p");
          } else {
            energyEnergyCorrelatorSignalPbPb[weightExponent-1][iDrawingStyle][iJetPt][iTrackPt]->Draw("same,p");
          }

          // Make a legend for the jet pT bins
          legend = new TLegend(0.18, 0.18, 0.47, 0.52);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(supplementaryLegendTextSize); legend->SetTextFont(62);

          // Add legends from each distribution
          for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
            legend->AddEntry(systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], Form("PbPb %.0f-%.0f%%", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality)), "lpf");
          }

          legend->AddEntry(systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt], "pp", "lpf");

          // Draw the legend
          legend->Draw();

          // Semi-transparent boxes for some drawing styles
          box = new TBox();
          box->SetFillColorAlpha(kWhite, 0.95);

          // The pp text needs to be hidden for drawing styles 0-3
          if(iDrawingStyle < 4){
            box->DrawBox(0.008, 0.49, 0.015, 0.7); // Hide the pp text
          }

          // The 50-90% text needs to be drawn for styles 3 and 6
          if(iDrawingStyle != 3 && iDrawingStyle !=6){
            box->DrawBox(0.008, 0.7, 0.041, 0.92); // Hide the 50-90% text
          }
          
          // The 30-50% text needs to be drawn for styles 2 and 6
          if(iDrawingStyle != 2 && iDrawingStyle != 6){
            box->DrawBox(0.008, 0.92, 0.041, 1.33); // Hide the 30-50% text
          }
          
          // The 10-30% text needs to be drawn for styles 1 and 6
          if(iDrawingStyle != 1 && iDrawingStyle != 6){
            box->DrawBox(0.008, 1.33, 0.041, 1.88); // Hide the 10-30% text
          }

          // The 0-10% text needs to be drawn for styles 0, 5 and 6
          if(iDrawingStyle > 0 && iDrawingStyle < 5){
            box->DrawBox(0.008, 1.88, 0.04, 2.75); // Hide the 0-10% text
          }

          // Add some interesting information to the plot
          mainTitle->SetTextFont(62);
          mainTitle->SetTextSize(0.06);
          mainTitle->DrawLatexNDC(0.185, 0.55, "CMS");

          if(addPreliminaryTag){
            mainTitle->SetTextFont(52);
            mainTitle->SetTextSize(0.045);
            mainTitle->DrawLatexNDC(0.315, 0.55, "Preliminary");
          }

          mainTitle->SetTextFont(42);
          mainTitle->SetTextSize(0.04);
          mainTitle->DrawLatexNDC(0.165, 0.95, "1.70 nb^{-1} PbPb (5.02 TeV) + 302 pb^{-1} pp (5.02 TeV)");

          mainTitle->SetTextFont(62);
          mainTitle->SetTextSize(supplementaryLegendTextSize);

          mainTitle->DrawLatexNDC(0.522, 0.85, jetPtString.Data());
          mainTitle->DrawLatexNDC(0.675, 0.78, "anti-k_{T} R = 0.4");
          mainTitle->DrawLatexNDC(0.765, 0.71, "|#eta_{jet}| < 1.6");
          mainTitle->DrawLatexNDC(0.726, 0.64, trackPtString.Data());
          mainTitle->DrawLatexNDC(0.87, 0.57, Form("n=%d", weightExponent));

          // Save the figures to file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/finalEECAllCentrality%s%s%s%s.%s", drawingStyleString[iDrawingStyle].Data(), saveComment.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
          }
        } // Drawing style loop
      } // Track pT loop
    } // Centrality loop
  } // If for drawing all centrality bins in one plot

  // Draw supplementary plots for PbPb/pp ratios with all centrality bins in the same figure
  if(drawSupplementaryEECRatioCentrality){

    TCanvas* theGreatCanvasOfCentralityRatio;
    int firstJetPtIndex = firstDrawnJetPtBinEEC[weightExponent-1];
    const int nDrawingStylesRatio = 5;
    TString drawingStyleStringRatio[nDrawingStylesRatio] = {"_onlyC0", "_onlyC1", "_onlyC2", "_onlyC3", ""};

    // Track pT loop
    for(int iTrackPt = firstDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt <= lastDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt++){

      trackPtString = Form("p_{T}^{ch} > %.0f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt));
      compactTrackPtString = Form("_T>%.1f", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt));
      compactTrackPtString.ReplaceAll(".","v");

      // Jet pT loop
      for(int iJetPt = firstDrawnJetPtBinEEC[weightExponent-1]; iJetPt <= lastDrawnJetPtBinEEC[weightExponent-1]; iJetPt++){

        jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt));
        compactJetPtString = Form("_J=%.0f-%.0f", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt));

        // Draw a canvas for each component separately for conference purposes
        for(int iDrawingStyle = 0; iDrawingStyle < nDrawingStylesRatio; iDrawingStyle++){

          theGreatCanvasOfCentralityRatio = new SplitCanvas(Form("theGreatCanvasOfCentralityRatio%d%d%d", iTrackPt, iJetPt, iDrawingStyle), "", 1300, 1100);
          theGreatCanvasOfCentralityRatio->SetMargin(0.18, 0.01, 0.15, 0.07);
          theGreatCanvasOfCentralityRatio->cd();

          gPad->SetLogy(false);
          gPad->SetLogx();

          mainTitle = new TLatex();

          // Set the axis titles and labels
          for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleOffset(0.85);
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleSize(0.08);
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelOffset(0.001);
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelSize(0.07);

            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleOffset(0.9);
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleSize(0.08);
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelOffset(0.01);
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelSize(0.07);

            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitle("#frac{PbPb}{pp}");
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetStats(0);
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetNdivisions(505);
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->CenterTitle();
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitle("#Deltar");
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetNdivisions(505);

            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetTitle("");

            // Set drawing ranges
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.25, 1.75);
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);

            // Set the drawing style
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);
            energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
            energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);

            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
            energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
            energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);

            // Different drawing styles for correlated uncertainties
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorUpPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorDownPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);
          }

          // First, draw all uncorrelated systematic uncertainties
          if(iDrawingStyle == 4){
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][lastDrawnCentralityBin[weightExponent-1]][iJetPt][iTrackPt]->Draw("e2");
            for(int iCentrality = lastDrawnCentralityBin[weightExponent-1]-1; iCentrality >= firstDrawnCentralityBin[weightExponent-1]; iCentrality--){
              systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Draw("same,e2");
            }
          } else {
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iDrawingStyle][iJetPt][iTrackPt]->Draw("e2");
          }

          // Then, draw all correlated uncertainties
          if(iDrawingStyle == 4){
            for(int iCentrality = lastDrawnCentralityBin[weightExponent-1]; iCentrality >= firstDrawnCentralityBin[weightExponent-1]; iCentrality--){
              systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
              systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
            }
          } else {
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeUp][iDrawingStyle][iJetPt][iTrackPt]->Draw("same,e3");
            systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeDown][iDrawingStyle][iJetPt][iTrackPt]->Draw("same,e3");
          }

          // Finally, draw the data points themselves
          if(iDrawingStyle == 4){
            for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
              energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]->Draw("same,p");
            }
          }  else {
            energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iDrawingStyle][iJetPt][iTrackPt]->Draw("same,p");
          }

          // Make a legend for the jet pT bins
          legend = new TLegend(0.22, 0.18, 0.54, 0.45);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(supplementaryLegendTextSize); legend->SetTextFont(62);

          // Add legends from each distribution
          for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
            legend->AddEntry(systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], Form("PbPb %.0f-%.0f%%", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality)), "lpf");
          }


          // Draw the legend
          legend->Draw();

          // Draw a lone to one
          oneLine->Draw();

          // Semi-transparent boxes for some drawing styles
          box = new TBox();
          box->SetFillColorAlpha(kWhite, 0.95);

          // The 50-90% text needs to be drawn for styles 3 and 4
          if(iDrawingStyle < 3){
            box->DrawBox(0.008, 0.32, 0.043, 0.42); // Hide the 50-90% text
          }
          
          // The 30-50% text needs to be drawn for styles 2 and 4
          if(iDrawingStyle != 2 && iDrawingStyle != 4){
            box->DrawBox(0.008, 0.42, 0.043, 0.55); // Hide the 30-50% text
          }
          
          // The 10-30% text needs to be drawn for styles 1 and 4
          if(iDrawingStyle != 1 && iDrawingStyle != 4){
            box->DrawBox(0.008, 0.55, 0.043, 0.68); // Hide the 10-30% text
          }

          // The 0-10% text needs to be drawn for styles 0 and 4
          if(iDrawingStyle != 0 && iDrawingStyle != 4){
            box->DrawBox(0.008, 0.68, 0.04, 0.81); // Hide the 0-10% text
          }

          // Add some interesting information to the plot
          mainTitle->SetTextFont(62);
          mainTitle->SetTextSize(0.06);
          mainTitle->DrawLatexNDC(0.2, 0.86, "CMS");

          if(addPreliminaryTag){
            mainTitle->SetTextFont(52);
            mainTitle->SetTextSize(0.045);
            mainTitle->DrawLatexNDC(0.31, 0.86, "Preliminary");
          }

          mainTitle->SetTextFont(42);
          mainTitle->SetTextSize(0.04);
          mainTitle->DrawLatexNDC(0.29, 0.95, "1.70 nb^{-1} PbPb (5.02 TeV) + 302 pb^{-1} pp (5.02 TeV)");

          mainTitle->SetTextFont(62);
          mainTitle->SetTextSize(supplementaryLegendTextSize);

          mainTitle->DrawLatexNDC(0.53, 0.86, jetPtString.Data());
          mainTitle->DrawLatexNDC(0.657, 0.79, "anti-k_{T} R = 0.4");
          mainTitle->DrawLatexNDC(0.735, 0.72, "|#eta_{jet}| < 1.6");
          mainTitle->DrawLatexNDC(0.621, 0.65, Form("%s, n=%d", trackPtString.Data(), weightExponent));

          // Save the figures to file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/finalEECRatioCentrality%s%s%s%s.%s", drawingStyleStringRatio[iDrawingStyle].Data(), saveComment.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
          }
        } // Drawing style loop
      } // Track pT loop
    } // Centrality loop
  } // If for drawing PbPb/pp ratios from all centrality bins to the same figure

  // Draw supplementary plots for PbPb/pp ratios with all centrality bins in the same figure
  if(drawSupplementaryEECRatioWeightExponent){

    TCanvas* theGreatCanvasOfRatioWeightExponentEdition;

    // Track pT loop
    for(int iTrackPt = firstDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt <= lastDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt++){

      trackPtString = Form("p_{T}^{ch} > %.0f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt));
      compactTrackPtString = Form("_T>%.1f", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt));
      compactTrackPtString.ReplaceAll(".","v");

      // Jet pT loop
      for(int iJetPt = firstDrawnJetPtBinEEC[weightExponent-1]; iJetPt <= lastDrawnJetPtBinEEC[weightExponent-1]; iJetPt++){

        jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt));
        compactJetPtString = Form("_J=%.0f-%.0f", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt));

        // Centrality loop
        for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
        
          centralityString = Form("PbPb %.0f-%.0f%%", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality));
          compactCentralityString = Form("_C=%.0f-%.0f", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality));

          // Create a TCanvas for each centrality+jet pT+track pT bin
          theGreatCanvasOfRatioWeightExponentEdition = new SplitCanvas(Form("theGreatCanvasOfWeightExponentRatio%d%d%d", iTrackPt, iJetPt, iCentrality), "", 1300, 1100);
          theGreatCanvasOfRatioWeightExponentEdition->SetMargin(0.18, 0.01, 0.15, 0.07);
          theGreatCanvasOfRatioWeightExponentEdition->cd();

          gPad->SetLogy(false);
          gPad->SetLogx();

          mainTitle = new TLatex();

          // Set the axis titles and labels
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleOffset(0.85);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleSize(0.08);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelOffset(0.001);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelSize(0.07);

          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleOffset(0.9);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleSize(0.08);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelOffset(0.01);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelSize(0.07);

          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitle("#frac{PbPb}{pp}");
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetStats(0);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetNdivisions(505);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitle("#Deltar");
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetNdivisions(505);

          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetTitle("");

          // Set drawing ranges
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.25, 1.75);
          for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.4);


            // Set the drawing style
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);
            energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iWeightExponent]);
            energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);

            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPbPb[iWeightExponent], 0.4);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iWeightExponent]);
            energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iWeightExponent]);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iWeightExponent]);
            energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iWeightExponent]);

            // Different drawing styles for correlated uncertainties
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorUpPbPb[iWeightExponent], 0.4);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(bandColorDownPbPb[iWeightExponent], 0.4);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);
          }

          // First, draw all uncorrelated systematic uncertainties
          systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Draw("e2");
          for(int iWeightExponent = 1; iWeightExponent < nWeightExponents; iWeightExponent++){
              systematicUncertaintyPbPbToPpRatio[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Draw("same,e2");
          }
          
          // Then, draw all correlated uncertainties
          for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          }

          // Finally, draw the data points themselves
          for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
            energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Draw("same,p");
          }

          // Make a legend for the jet pT bins
          legend = new TLegend(0.22, 0.2, 0.48, 0.34);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(supplementaryLegendTextSize); legend->SetTextFont(62);

          // Add legends from each distribution
          for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
            legend->AddEntry(systematicUncertaintyPbPbToPpRatio[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], Form("n = %d", iWeightExponent+1), "lpf");
          }

          // Draw the legend
          legend->Draw();

          // Draw a lone to one
          oneLine->Draw();

          // Add some interesting information to the plot
          mainTitle->SetTextFont(62);
          mainTitle->SetTextSize(0.06);
          mainTitle->DrawLatexNDC(0.2, 0.86, "CMS");

          if(addPreliminaryTag){
            mainTitle->SetTextFont(52);
            mainTitle->SetTextSize(0.045);
            mainTitle->DrawLatexNDC(0.31, 0.86, "Preliminary");
          }

          mainTitle->SetTextFont(42);
          mainTitle->SetTextSize(0.04);
          mainTitle->DrawLatexNDC(0.29, 0.95, "1.70 nb^{-1} PbPb (5.02 TeV) + 302 pb^{-1} pp (5.02 TeV)");

          mainTitle->SetTextFont(62);
          mainTitle->SetTextSize(supplementaryLegendTextSize);

          mainTitle->DrawLatexNDC(0.53, 0.86, jetPtString.Data());
          mainTitle->DrawLatexNDC(0.657, 0.79, "anti-k_{T} R = 0.4");
          mainTitle->DrawLatexNDC(0.735, 0.72, "|#eta_{jet}| < 1.6");
          mainTitle->DrawLatexNDC(0.701, 0.65, trackPtString.Data());

          mainTitle->DrawLatexNDC(0.22, 0.38, Form("#frac{%s}{pp}", centralityString.Data()));

          // Add some interesting information to the plot
          /*mainTitle->SetTextFont(62);
          mainTitle->SetTextSize(0.08);

          if(addPreliminaryTag){
            mainTitle->DrawLatexNDC(0.065, 0.93, "CMS");

            mainTitle->SetTextFont(52);
            mainTitle->SetTextSize(0.05);
            mainTitle->DrawLatexNDC(0.035, 0.88, "Preliminary");
          } else {
            mainTitle->DrawLatexNDC(0.02, 0.89, "CMS");
          }

          mainTitle->SetTextFont(42);
          mainTitle->SetTextSize(0.055);

          if(addPreliminaryTag){
            mainTitle->DrawLatexNDC(0.33, 0.95, "1.70 nb^{-1} PbPb (5.02 TeV)");
            mainTitle->DrawLatexNDC(0.33, 0.88, "302 pb^{-1} pp (5.02 TeV)");
          } else {
            mainTitle->DrawLatexNDC(0.24, 0.95, "1.70 nb^{-1} PbPb (5.02 TeV)");
            mainTitle->DrawLatexNDC(0.24, 0.88, "302 pb^{-1} pp (5.02 TeV)");
          }

          mainTitle->DrawLatexNDC(0.34, 0.78, jetPtString.Data());
          mainTitle->DrawLatexNDC(0.57, 0.71, trackPtString.Data());
          mainTitle->DrawLatexNDC(0.51, 0.63, "anti-k_{T} R = 0.4");
          mainTitle->DrawLatexNDC(0.61, 0.56, "|#eta_{jet}| < 1.6");

          
          */

          // Save the figures to file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/finalEECRatioEnergyWeight%s%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
          }
        } // Drawing style loop
      } // Track pT loop
    } // Centrality loop
  } // If for drawing PbPb/pp ratios from all energy weights to the same figure
}
