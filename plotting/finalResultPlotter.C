#include "SystematicUncertaintyOrganizer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECHistogramManager.h"
#include "EECCard.h"
#include "JDrawer.h"
#include "SplitCanvas.h"
#include "AlgorithmLibrary.h"

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
void finalResultPlotter(bool drawOnlyPaperPlots = false){

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
  uncertaintyFileName[kPbPb][0] = "systematicUncertainties/systematicUncertainties_PbPb_nominalEnergyWeight_combinedMxedConeBackground_includeMCstats_2025-04-18.root";
  uncertaintyFileName[kPbPb][1] = "systematicUncertainties/systematicUncertainties_PbPb_energyWeightSquared_combinedMxedConeBackground_includeMCstats_2025-04-18.root";
  uncertaintyFileName[kPp][0] = "systematicUncertainties/systematicUncertainties_pp_nominalEnergyWeight_includeMCstats_2025-04-21.root";
  uncertaintyFileName[kPp][1] = "systematicUncertainties/systematicUncertainties_pp_energyWeightSquared_includeMCstats_2025-04-21.root";

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
  bool drawSupplementaryEECJetPt = false; // Draw all jet pT bins for selected centrality, track pT, and energy weight to the same plot
  bool drawSupplementaryEECCentrality = true; // Draw all centrality bins for selected jet pT, track pT and energy weight to the smae plot
  bool drawSupplementaryEECRatioCentrality = false; // Draw PbPb/pp ratios from all centrality bins to the same figure
  bool drawSupplementaryEECRatioWeightExponent = false; // Draw PbPb/pp ratios from all weight exponents to the same figure
  bool drawSupplementaryEECShiftIllustration = false; // Draw two pp distributions and one PbPb distribution to the same figure

  // Option to hide the shift histogram for shift illustration
  bool hideShiftInIllustration = false; // Hide the shifted pp distribution in the shift illustration

  // Save all the drawn histograms from HepData
  bool saveHepDataFiles = false;
  TString hepDataTrackPtString;

  // Configuration for supplementary plots
  double supplementaryLegendTextSize = 0.045;

  // Preliminary/Supplementary tag
  enum enumTagType{kNoTag, kPreliminaryTag, kSupplementaryTag, knTagTypes};
  int tagSelection = kNoTag;
  const char* tagName[knTagTypes];
  tagName[kNoTag] = "";
  tagName[kPreliminaryTag] = "Preliminary";
  tagName[kSupplementaryTag] = "Supplementary";

  double tagShift[knTagTypes] = {0,0,-0.07};

  // Simplified drawing style for physics briefing
  bool physicsBriefingStyle = false;

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
  TString saveComment =  "";
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
  int markerStylePpEnergyWeight[] = {kFullCrossX, kFullFourTrianglesPlus};
  int markerColorPpEnergyWeight[] = {kOrange, kCyan};

  // Definition on how much the jet pT is shifted for shifted pT figures
  std::vector<std::pair<std::pair<double,double>,std::pair<double,double>>> jetPtShiftDefinition;
  jetPtShiftDefinition.push_back(std::make_pair(std::make_pair(120,140),std::make_pair(135,155)));
  jetPtShiftDefinition.push_back(std::make_pair(std::make_pair(140,160),std::make_pair(155,175)));
  jetPtShiftDefinition.push_back(std::make_pair(std::make_pair(160,180),std::make_pair(175,195)));
  jetPtShiftDefinition.push_back(std::make_pair(std::make_pair(180,200),std::make_pair(195,215)));

  TLine* lineDrawer = new TLine();
  lineDrawer->SetLineStyle(2);
  lineDrawer->SetLineColor(kBlack);

  // If we are drawing only the paper plots, setup configuration for that
  if(drawOnlyPaperPlots){
    drawIndividualPlotsAllCentralities = false;
    drawBigCanvasDistributions = false;
    drawBigCanvasRatios = false;
    drawDoubleRatios = false;
    drawDoubleRatioToSingleCanvas = false;
    drawBigCanvasAllDistributions = true; // Draw distributions with all defined weight exponents to the same figure
    drawBigCanvasAllRatios = true; // Draw ratios with all defined energy weight exponents to the same figure
    drawLetterPaperDistributions = false; // Draw the energy-energy correlator distribution for letter paper format
    drawLetterPaperRatios = false; // Draw the energy-energy correlators ratios for letter paper format

    // Supplementary plots
    drawSupplementaryEECJetPt = false; // Draw all jet pT bins for selected centrality, track pT, and energy weight to the same plot
    drawSupplementaryEECCentrality = false; // Draw all centrality bins for selected jet pT, track pT and energy weight to the smae plot
    drawSupplementaryEECRatioCentrality = false; // Draw PbPb/pp ratios from all centrality bins to the same figure
    drawSupplementaryEECRatioWeightExponent = false; // Draw PbPb/pp ratios from all weight exponents to the same figure
    drawSupplementaryEECShiftIllustration = false; // Draw two pp distributions and one PbPb distribution to the same figure

    // Save all the drawn histograms from HepData
    saveHepDataFiles = false;

    // No tags for the figures
    tagSelection = kNoTag;

    // Default figure naming
    saveComment =  "";
    figureFormat = "pdf";
  }

  // Algorithm library for any useful manipulation
  AlgorithmLibrary* manipulator = new AlgorithmLibrary();

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

      // Create a new systematic uncertainty organizer
      uncertainties[iDataType][iWeightExponent] = new SystematicUncertaintyOrganizer(uncertaintyFile[iDataType][iWeightExponent]);

    } // Data type loop

    // Create a new histogram manager for histograms with shifted pT
    if(drawShiftedPtRatio){
      shiftedPtHistograms[iWeightExponent] = new EECHistogramManager(shiftedPtFile[iWeightExponent], shiftedPtCard[iWeightExponent]);
    }

  } // Weight exponent loop

  // Efficiency note: The element access operator [] for std::map works in logarithmic time, while for std::vector 
  // the time is constant. However, I think single map object for each type type of histogram is cleaner 
  // implementation compared to nested vectors. Since we are not pressed on processing time in a plotting macro,
  // std::map is a reasonable choice for a container to keep the histograms in.
 
  // Energy-energy correlators and PbPb to pp ratios
  std::map<std::tuple<int, int, int, int>, TH1D*> energyEnergyCorrelatorRawPbPb; // Dimensions: weight exponent, centrality, jet pT, track pT
  std::map<std::tuple<int, int, int>, TH1D*> energyEnergyCorrelatorRawPp; // Dimensions: weight exponent, jet pT, track pT
  std::map<std::tuple<int, int, int, int>, TH1D*> energyEnergyCorrelatorSignalPbPb; // Dimensions: weight exponent, centrality, jet pT, track pT
  std::map<std::tuple<int, int, int>, TH1D*> energyEnergyCorrelatorSignalPp; // Dimensions: weight exponent, jet pT, track pT
  std::map<std::tuple<int, int, int>, TH1D*> energyEnergyCorrelatorSignalShiftedPt; // Dimensions: weight exponent, jet pT, track pT
  std::map<std::tuple<int, int, int, int>, TH1D*> energyEnergyCorrelatorPbPbToPpRatio; // Dimensions: weight exponent, centrality, jet pT, track pT
  std::map<std::tuple<int, int, int>, TH1D*> energyEnergyCorrelatorShiftedPtRatio; // Dimensions: weight exponent, jet pT, track pT
  std::map<std::tuple<int, int, int, int, int>, TH1D*> systematicUncertaintyForPbPb; // Dimensions: weight exponent, systematic type, centrality, jet pT, track pT
  std::map<std::tuple<int, int, int, int>, TH1D*> systematicUncertaintyForPp; // Dimensions: weight exponent, systematic type, jet pT, track pT
  std::map<std::tuple<int, int, int, int, int>, TH1D*> systematicUncertaintyPbPbToPpRatio; // Dimensions: weight exponent, systematic type, centrality, jet pT, track pT

  // Double ratios from energy-energy correlator histograms
  std::map<std::tuple<int, int, int, int>, TH1D*> energyEnergyCorrelatorForDoubleRatioFromPbPb; // Dimensions: weight exponent, centrality, jet pT, track pT
  std::map<std::tuple<int, int, int>, TH1D*> energyEnergyCorrelatorForDoubleRatioFromPp; // Dimensions: weight exponent, jet pT, track pT
  std::map<std::tuple<int, int, int>, TH1D*> energyEnergyCorrelatorDoubleRatio; // Dimensions: weight exponent, centrality, jet pT
  std::map<std::tuple<int, int, int, int>, TH1D*> systematicUncertaintyForDoubleRatioFromPbPb; // Dimensions: weight exponent, centrality, jet pT, track pT
  std::map<std::tuple<int, int, int>, TH1D*> systematicUncertaintyForDoubleRatioFromPp; // Dimensions: weight exponent, jet pT, track pT
  std::map<std::tuple<int, int, int>, TH1D*> systematicUncertaintyDoubleRatio; // Dimensions: weight exponent, centrality, jet pT

  // Binnning variables for the histograms
  std::tuple<int, int, int, int> eecPbPbBin;
  std::tuple<int, int, int, int> eecPbPbReferenceBin;
  std::tuple<int, int, int, int, int> eecPbPbCorrelatedSystematicsBin;
  std::tuple<int, int, int, int, int> eecPbPbCorrelatedSystematicsBinShapeUp;
  std::tuple<int, int, int, int, int> eecPbPbCorrelatedSystematicsBinShapeDown;
  std::tuple<int, int, int, int, int> eecPbPbUncorrelatedSystematicsBin;

  std::tuple<int, int, int> eecPpBin;
  std::tuple<int, int, int> eecPpReferenceBin;
  std::tuple<int, int, int, int> eecPpCorrelatedSystematicsBin;
  std::tuple<int, int, int, int> eecPpCorrelatedSystematicsBinShapeUp;
  std::tuple<int, int, int, int> eecPpCorrelatedSystematicsBinShapeDown;
  std::tuple<int, int, int, int> eecPpUncorrelatedSystematicsBin;
  std::tuple<int, int, int, int> eecPpReferenceUncorrelatedSystematicsBin;
  std::tuple<int, int, int, int> eecPpReferenceCorrelatedSystematicsBin;

  std::tuple<int, int, int> eecDoubleRatioBin;

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
        eecPpBin = std::make_tuple(iWeightExponent, iJetPt, iTrackPt);

        // Read the raw pp energy-energy correlator distributions
        energyEnergyCorrelatorRawPp[eecPpBin] = histograms[kPp][iWeightExponent]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatchedPp, iTrackPtMatchedPp);

        for(int iCentrality = firstDrawnCentralityBin[iWeightExponent]; iCentrality <= lastDrawnCentralityBin[iWeightExponent]; iCentrality++){
          eecPbPbBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt);

          // Read the raw PbPb energy-energy correlator distributions
          energyEnergyCorrelatorRawPbPb[eecPbPbBin] = histograms[kPbPb][iWeightExponent]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt);
        }
      }

      // Then go to the main track pT loop
      for(int iTrackPt = firstDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt <= lastDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt++){
        iTrackPtMatchedPp = card[kPp][iWeightExponent]->FindBinIndexTrackPtEEC(card[kPbPb][iWeightExponent]->GetBinBordersTrackPtEEC(iTrackPt));
        iTrackPtMatchedPbPbUncertainty = uncertaintyCard[kPbPb][iWeightExponent]->FindBinIndexTrackPtEEC(card[kPbPb][iWeightExponent]->GetBinBordersTrackPtEEC(iTrackPt));
        iTrackPtMatchedPpUncertainty = uncertaintyCard[kPp][iWeightExponent]->FindBinIndexTrackPtEEC(card[kPbPb][iWeightExponent]->GetBinBordersTrackPtEEC(iTrackPt));
        if(drawShiftedPtRatio) iTrackPtMatchedShifted = shiftedPtCard[iWeightExponent]->FindBinIndexTrackPtEEC(card[kPbPb][iWeightExponent]->GetBinBordersTrackPtEEC(iTrackPt));

        eecPpBin = std::make_tuple(iWeightExponent, iJetPt, iTrackPt);
        eecPpUncorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kUncorrelatedUncertainty, iJetPt, iTrackPt);
        eecPpCorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kCorrelatedUncertainty, iJetPt, iTrackPt);

        // Read the pp histograms that do not have centrality binning
        energyEnergyCorrelatorSignalPp[eecPpBin] = histograms[kPp][iWeightExponent]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatchedPp, iTrackPtMatchedPp, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);
        systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin] = uncertainties[kPp][iWeightExponent]->GetUncorrelatedSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty);
        systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin] = uncertainties[kPp][iWeightExponent]->GetCorrelatedSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty);

        // Read the relevant systematic uncertainties for double ratio from pp
        uncertaintyBackgroundSubtraction = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kBackgroundSubtraction);
        uncertaintyTrackPairEfficiency = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kTrackPairEfficiency);
        uncertaintyMCnonclosure = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kMonteCarloNonClosure);
        uncertaintyTrackSelection = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kTrackSelection);

        // For double ratio, we can find the proper scaling factor for the tracking uncertainties by comparing integrals of raw
        // energy-energy correlators between different track pT bins. The higher track pT cuts are subsets of the lower cuts.
        // Thus integrals over the analysis region tell the degree of overlap of these regions.
        lowAnalysisBin = energyEnergyCorrelatorRawPp[eecPpBin]->GetXaxis()->FindBin(analysisDeltaR.first + epsilon);
        highAnalysisBin = energyEnergyCorrelatorRawPp[eecPpBin]->GetXaxis()->FindBin(analysisDeltaR.second - epsilon);
        lowPtIntegral = energyEnergyCorrelatorRawPp[std::make_tuple(iWeightExponent,iJetPt,trackPtBinsForDoubleRatio.first)]->Integral(lowAnalysisBin, highAnalysisBin, "width");
        highPtIntegral = energyEnergyCorrelatorRawPp[std::make_tuple(iWeightExponent,iJetPt,trackPtBinsForDoubleRatio.second)]->Integral(lowAnalysisBin, highAnalysisBin, "width");
        trackCorrelation = 1 - highPtIntegral / lowPtIntegral;

        // Now for the double ratio, we can combine the relevant uncertainties while scaling down the tracking related ones
        // by the expected overlap.
        systematicUncertaintyForDoubleRatioFromPp[eecPpBin] = (TH1D*) uncertaintyBackgroundSubtraction->Clone(Form("doubleRatioUncertaintyFromPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));
        for(int iBin = 1; iBin <= uncertaintyBackgroundSubtraction->GetNbinsX(); iBin++){
          systematicUncertaintyForDoubleRatioFromPp[eecPpBin]->SetBinError(iBin, TMath::Sqrt(TMath::Power(uncertaintyBackgroundSubtraction->GetBinError(iBin),2) + TMath::Power(uncertaintyTrackPairEfficiency->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintyTrackSelection->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintyMCnonclosure->GetBinError(iBin),2)));
        }

        // Read the histograms with shifted pT spectrum
        if(drawShiftedPtRatio){
          energyEnergyCorrelatorSignalShiftedPt[eecPpBin] = shiftedPtHistograms[iWeightExponent]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatchedShifted, iTrackPtMatchedShifted, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);
        }


        for(int iCentrality = firstDrawnCentralityBin[iWeightExponent]; iCentrality <= lastDrawnCentralityBin[iWeightExponent]; iCentrality++){
          iCentralityMatched = uncertaintyCard[kPbPb][iWeightExponent]->FindBinIndexCentrality(card[kPbPb][iWeightExponent]->GetBinBordersCentrality(iCentrality));
          eecPbPbBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt);
          eecPbPbUncorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kCorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);

          // Read the PbPb histograms
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin] = histograms[kPbPb][iWeightExponent]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin] = uncertainties[kPbPb][iWeightExponent]->GetUncorrelatedSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBin] = uncertainties[kPbPb][iWeightExponent]->GetCorrelatedSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty);

          // Read the relevant systematic uncertainties for double ratio from PbPb
          uncertaintyBackgroundSubtraction = uncertainties[kPbPb][iWeightExponent]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kBackgroundSubtraction);
          uncertaintyTrackPairEfficiency = uncertainties[kPbPb][iWeightExponent]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kTrackPairEfficiency);
          uncertaintyTrackSelection = uncertainties[kPbPb][iWeightExponent]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kTrackSelection);
          uncertaintyMCnonclosure = uncertainties[kPbPb][iWeightExponent]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kMonteCarloNonClosure);
          uncertaintySignalToBackgroundRatio = uncertainties[kPbPb][iWeightExponent]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kSignalToBackgroundRatio);

          // For double ratio, we can find the proper scaling factor for the tracking uncertainties by comparing integrals of raw
          // energy-energy correlators between different track pT bins. The higher track pT cuts are subsets of the lower cuts.
          // Thus integrals over the analysis region tell the degree of overlap of these regions.
          lowAnalysisBin = energyEnergyCorrelatorRawPbPb[eecPbPbBin]->GetXaxis()->FindBin(analysisDeltaR.first + epsilon);
          highAnalysisBin = energyEnergyCorrelatorRawPbPb[eecPbPbBin]->GetXaxis()->FindBin(analysisDeltaR.second - epsilon);
          lowPtIntegral = energyEnergyCorrelatorRawPbPb[std::make_tuple(iWeightExponent,iCentrality,iJetPt,trackPtBinsForDoubleRatio.first)]->Integral(lowAnalysisBin, highAnalysisBin, "width");
          highPtIntegral = energyEnergyCorrelatorRawPbPb[std::make_tuple(iWeightExponent,iCentrality,iJetPt,trackPtBinsForDoubleRatio.second)]->Integral(lowAnalysisBin, highAnalysisBin, "width");
          trackCorrelation = 1 - highPtIntegral / lowPtIntegral;

          // Now for the double ratio, we can combine the relevant uncertainties while scaling down the tracking related ones
          // by the expected overlap.
          systematicUncertaintyForDoubleRatioFromPbPb[eecPbPbBin] = (TH1D*) uncertaintyBackgroundSubtraction->Clone(Form("doubleRatioUncertaintyFromPbPb%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          for(int iBin = 1; iBin <= uncertaintyBackgroundSubtraction->GetNbinsX(); iBin++){
            systematicUncertaintyForDoubleRatioFromPbPb[eecPbPbBin]->SetBinError(iBin, TMath::Sqrt(TMath::Power(uncertaintyBackgroundSubtraction->GetBinError(iBin),2) + TMath::Power(uncertaintyTrackPairEfficiency->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintyTrackSelection->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintySignalToBackgroundRatio->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintyMCnonclosure->GetBinError(iBin),2)));
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
    eecPbPbBin = std::make_tuple(iWeightExponent, firstDrawnCentralityBin[iWeightExponent], firstDrawnJetPtBinEEC[iWeightExponent], firstDrawnTrackPtBinEEC[iWeightExponent]);
    lowAnalysisBin = energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->GetXaxis()->FindBin(analysisDeltaR.first + epsilon);
    highAnalysisBin = energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->GetXaxis()->FindBin(analysisDeltaR.second - epsilon);
    for(int iJetPt = firstDrawnJetPtBinEEC[iWeightExponent]; iJetPt <= lastDrawnJetPtBinEEC[iWeightExponent]; iJetPt++){
      for(int iTrackPt = lastDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt >= firstDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt--){
        if(!normalizeTo2GeV) referenceTrackPtBin = iTrackPt;
        eecPpBin = std::make_tuple(iWeightExponent, iJetPt, iTrackPt);
        eecPpReferenceBin = std::make_tuple(iWeightExponent, iJetPt, referenceTrackPtBin);
        eecPpUncorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kUncorrelatedUncertainty, iJetPt, iTrackPt);
        eecPpCorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kCorrelatedUncertainty, iJetPt, iTrackPt);

        energyEnergyCorrelatorSignalPp[eecPpBin]->Scale(1.0 / energyEnergyCorrelatorSignalPp[eecPpReferenceBin]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
        systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->Scale(energyEnergyCorrelatorSignalPp[eecPpBin]->GetBinContent(10) / systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetBinContent(10));
        systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->Scale(energyEnergyCorrelatorSignalPp[eecPpBin]->GetBinContent(10) / systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->GetBinContent(10));
        systematicUncertaintyForDoubleRatioFromPp[eecPpBin]->Scale(energyEnergyCorrelatorSignalPp[eecPpBin]->GetBinContent(10) / systematicUncertaintyForDoubleRatioFromPp[eecPpBin]->GetBinContent(10));

        // Shifted pT to regular pT ratios
        if(drawShiftedPtRatio){
          energyEnergyCorrelatorSignalShiftedPt[eecPpBin]->Scale(1.0 / energyEnergyCorrelatorSignalShiftedPt[eecPpReferenceBin]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

          energyEnergyCorrelatorShiftedPtRatio[eecPpBin] = (TH1D*) energyEnergyCorrelatorSignalShiftedPt[eecPpBin]->Clone(Form("shiftedPtRatio%d%d%d", iWeightExponent, iJetPt, iTrackPt));
          energyEnergyCorrelatorShiftedPtRatio[eecPpBin]->Divide(energyEnergyCorrelatorSignalPp[eecPpBin]);
        }

        for(int iCentrality = firstDrawnCentralityBin[iWeightExponent]; iCentrality <= lastDrawnCentralityBin[iWeightExponent]; iCentrality++){
          eecPbPbBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt);
          eecPbPbReferenceBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, referenceTrackPtBin);
          eecPbPbUncorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kCorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
          eecPpUncorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kUncorrelatedUncertainty, iJetPt, iTrackPt);
          eecPpCorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kCorrelatedUncertainty, iJetPt, iTrackPt);

          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->Scale(1.0 / energyEnergyCorrelatorSignalPbPb[eecPbPbReferenceBin]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->Scale(energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->GetBinContent(10) / systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetBinContent(10));
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBin]->Scale(energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->GetBinContent(10) / systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBin]->GetBinContent(10));
          systematicUncertaintyForDoubleRatioFromPbPb[eecPbPbBin]->Scale(energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->GetBinContent(10) / systematicUncertaintyForDoubleRatioFromPbPb[eecPbPbBin]->GetBinContent(10));

          energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin] = (TH1D*) energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->Clone(Form("energyEnergyCorrelatorRatio%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->Divide(energyEnergyCorrelatorSignalPp[eecPpBin]);

          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin] = (TH1D*) systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->Clone(Form("systematicUncertaintyRatioUncorrelated%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->Divide(systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]);

          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBin] = (TH1D*) systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBin]->Clone(Form("systematicUncertaintyRatioCorrelated%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBin]->Divide(systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]);
        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // Weight exponent loop

  // Reduce the statistical uncertainties by the overlapping statistics fraction for the double ratios
  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    for(int iJetPt = firstDrawnJetPtBinEEC[iWeightExponent]; iJetPt <= lastDrawnJetPtBinEEC[iWeightExponent]; iJetPt++){
      for(int iTrackPt = firstDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt <= lastDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt++){
        eecPpBin = std::make_tuple(iWeightExponent, iJetPt, iTrackPt);

        energyEnergyCorrelatorForDoubleRatioFromPp[eecPpBin] = (TH1D*) energyEnergyCorrelatorSignalPp[eecPpBin]->Clone(Form("energyEnergyCorrelatorForDoubleRatioFromPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));

        for(int iBin = 1; iBin <= energyEnergyCorrelatorForDoubleRatioFromPp[eecPpBin]->GetNbinsX(); iBin++){
          energyEnergyCorrelatorForDoubleRatioFromPp[eecPpBin]->SetBinError(iBin, energyEnergyCorrelatorForDoubleRatioFromPp[eecPpBin]->GetBinError(iBin) * trackCorrelation);
        }

        for(int iCentrality = firstDrawnCentralityBin[iWeightExponent]; iCentrality <= lastDrawnCentralityBin[iWeightExponent]; iCentrality++){
          eecPbPbBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt);
          energyEnergyCorrelatorForDoubleRatioFromPbPb[eecPbPbBin] = (TH1D*) energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->Clone(Form("energyEnergyCorrelatorForDoubleRatioFromPbPb%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));

          for(int iBin = 1; iBin <= energyEnergyCorrelatorForDoubleRatioFromPbPb[eecPbPbBin]->GetNbinsX(); iBin++){
            energyEnergyCorrelatorForDoubleRatioFromPbPb[eecPbPbBin]->SetBinError(iBin, energyEnergyCorrelatorForDoubleRatioFromPbPb[eecPbPbBin]->GetBinError(iBin) * trackCorrelation);
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
            eecPpBin = std::make_tuple(iWeightExponent, iJetPt, iTrackPt);
            eecPbPbBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt);

            // Calculate the single ratios with properly handled double ratio uncertainties
            energyEnergyCorrelatorForDoubleRatioFromPbPb[eecPbPbBin]->Divide(energyEnergyCorrelatorForDoubleRatioFromPp[eecPpBin]);
            systematicUncertaintyForDoubleRatioFromPbPb[eecPbPbBin]->Divide(systematicUncertaintyForDoubleRatioFromPp[eecPpBin]);
          } // Track pT loop
          
          eecDoubleRatioBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt);

          // Calculate the double ratios from the single ratios with properly handled uncertainties
          energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin] = (TH1D*) energyEnergyCorrelatorForDoubleRatioFromPbPb[std::make_tuple(iWeightExponent,iCentrality,iJetPt,trackPtBinsForDoubleRatio.second)]->Clone(Form("energyEnergyCorrelatorDoubleRatio%d%d%d", iWeightExponent, iCentrality, iJetPt));
          energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin]->Divide(energyEnergyCorrelatorForDoubleRatioFromPbPb[std::make_tuple(iWeightExponent,iCentrality,iJetPt,trackPtBinsForDoubleRatio.first)]);
          systematicUncertaintyDoubleRatio[eecDoubleRatioBin] = (TH1D*) systematicUncertaintyForDoubleRatioFromPbPb[std::make_tuple(iWeightExponent,iCentrality,iJetPt,trackPtBinsForDoubleRatio.second)]->Clone(Form("systematicUncertaintyDoubleRatio%d%d%d", iWeightExponent, iCentrality, iJetPt));
          systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->Divide(systematicUncertaintyForDoubleRatioFromPbPb[std::make_tuple(iWeightExponent,iCentrality,iJetPt,trackPtBinsForDoubleRatio.first)]);
        } // Jet pT loop
      } // Centrality llop
    } // Weight exponent loop
  }

  // For illustration purposes, create up and down shifted uncertainty bands
  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    for(int iJetPt = firstDrawnJetPtBinEEC[iWeightExponent]; iJetPt <= lastDrawnJetPtBinEEC[iWeightExponent]; iJetPt++){
      for(int iTrackPt = firstDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt <= lastDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt++){
        eecPpCorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kCorrelatedUncertainty, iJetPt, iTrackPt);
        eecPpCorrelatedSystematicsBinShapeUp = std::make_tuple(iWeightExponent, kCorrelatedUncertaintyShapeUp, iJetPt, iTrackPt);
        eecPpCorrelatedSystematicsBinShapeDown = std::make_tuple(iWeightExponent, kCorrelatedUncertaintyShapeDown, iJetPt, iTrackPt);

        systematicUncertaintyForPp[eecPpCorrelatedSystematicsBinShapeUp] = (TH1D*) systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->Clone(Form("correlatedUpBandPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));
        systematicUncertaintyForPp[eecPpCorrelatedSystematicsBinShapeDown] = (TH1D*) systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->Clone(Form("correlatedDownBandPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));
        calculateCorrelatedBands(systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin], systematicUncertaintyForPp[eecPpCorrelatedSystematicsBinShapeUp], systematicUncertaintyForPp[eecPpCorrelatedSystematicsBinShapeDown]);

        for(int iCentrality = firstDrawnCentralityBin[iWeightExponent]; iCentrality <= lastDrawnCentralityBin[iWeightExponent]; iCentrality++){
          eecPbPbCorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kCorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(iWeightExponent, kCorrelatedUncertaintyShapeUp, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(iWeightExponent, kCorrelatedUncertaintyShapeDown, iCentrality, iJetPt, iTrackPt);

          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp] = (TH1D*) systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBin]->Clone(Form("correlatedUpBand%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown] = (TH1D*) systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBin]->Clone(Form("correlatedDownBand%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          calculateCorrelatedBands(systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBin], systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp], systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]);

          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp] = (TH1D*) systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBin]->Clone(Form("correlatedUpBandForRatio%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown] = (TH1D*) systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBin]->Clone(Form("correlatedDownBandForRatio%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          calculateCorrelatedBands(systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBin], systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp], systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]);

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
  double bottomRowScale, bottomPadMargin, leftColumnScale, leftPadMargin;
  double leftMarginAdder, bottomMarginAdder;
  double thisPadScale;
  double firstRowScale;
  double cmsPosition;

  // Draw individual plots with all centralities mixed together in a single figure
  if(drawIndividualPlotsAllCentralities){

    for(int iJetPt = firstDrawnJetPtBinEEC[weightExponent-1]; iJetPt <= lastDrawnJetPtBinEEC[weightExponent-1]; iJetPt++){
      jetPtString = Form("%.0f < jet p_{T} < %.0f", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt));
      compactJetPtString = Form("_J=%.0f-%.0f", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt));

      for(int iTrackPt = firstDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt <= lastDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt++){
        trackPtString = Form("%.1f < track p_{T}", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt));
        compactTrackPtString = Form("_T>%.1f", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt));
        compactTrackPtString.ReplaceAll(".","v");

        eecPpBin = std::make_tuple(weightExponent-1, iJetPt, iTrackPt);
        eecPbPbReferenceBin = std::make_tuple(weightExponent-1, lastDrawnCentralityBin[weightExponent-1], iJetPt, iTrackPt);

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
        energyEnergyCorrelatorSignalPp[eecPpBin]->SetMarkerStyle(kFullDiamond);
        energyEnergyCorrelatorSignalPp[eecPpBin]->SetMarkerColor(kBlack);

        // Set the x-axis drawing range
        energyEnergyCorrelatorSignalPp[eecPpBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);

        // Draw the pp correlator to upper canves
        drawer->DrawHistogramToUpperPad(energyEnergyCorrelatorSignalPp[eecPpBin], "#Deltar", "EEC Signal", " ", "p");
        legend->AddEntry(energyEnergyCorrelatorSignalPp[eecPpBin], "pp", "p");

        // Draw the different centrality bins to the same plot
        for(int iCentrality = lastDrawnCentralityBin[weightExponent-1]; iCentrality >= firstDrawnCentralityBin[weightExponent-1]; iCentrality--) {
          eecPbPbBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt);
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->Draw("same,p");
          legend->AddEntry(energyEnergyCorrelatorSignalPbPb[eecPbPbBin], Form("PbPb %.0f-%.0f%%", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality)), "p");
        }

        // Draw the legend to the upper pad
        legend->Draw();

        // Linear scale for the ratio
        drawer->SetLogY(false);

        // Set the axis drawing ranges
        energyEnergyCorrelatorPbPbToPpRatio[eecPbPbReferenceBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
        energyEnergyCorrelatorPbPbToPpRatio[eecPbPbReferenceBin]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

        // Set the style for histograms
        for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++) {
          eecPbPbBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt);
          energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
        }

        drawer->SetGridY(true);
        drawer->DrawHistogramToLowerPad(energyEnergyCorrelatorPbPbToPpRatio[eecPbPbReferenceBin], "#Deltar", "#frac{PbPb}{pp}", " ");
        for(int iCentrality = lastDrawnCentralityBin[weightExponent-1] - 1; iCentrality >= firstDrawnCentralityBin[weightExponent-1]; iCentrality--) {
          eecPbPbBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt);
          energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->Draw("same,p");
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
    std::vector<SplitCanvas*> bigCanvas;
    bigCanvas.resize(nTrackPtBinsEEC);
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

          eecPpBin = std::make_tuple(weightExponent-1, iJetPt, iTrackPt);
          eecPpUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iJetPt, iTrackPt);
          eecPpCorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kCorrelatedUncertainty, iJetPt, iTrackPt);

          eecPbPbBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt);
          eecPbPbUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeUp, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeDown, iCentrality, iJetPt, iTrackPt);


          canvasIndex = (iCentrality - firstDrawnCentralityBin[weightExponent-1]) * (lastDrawnJetPtBinEEC[weightExponent-1] - firstDrawnJetPtBinEEC[weightExponent-1] + 1) + (iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]);
          bigCanvas[iTrackPt]->CD(canvasIndex);

          gPad->SetLogy();
          gPad->SetLogx();

          // The titles in the bottow row need to be scaled, since it has more margin than other rows
          thisPadScale = (iCentrality == lastDrawnCentralityBin[weightExponent-1]) ? bottomRowScale : 1;

          // Set the axis titles and labels
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitleOffset(0.85);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitleSize(0.2 * thisPadScale);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetLabelOffset(0.01);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetLabelSize(0.12 * thisPadScale);

          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitleOffset(0.5 / thisPadScale);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitleSize(0.2 * thisPadScale);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetLabelOffset(0.01 / thisPadScale);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetLabelSize(0.12 * thisPadScale);

          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitle("EEC");
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetStats(0);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetNdivisions(505);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->CenterTitle();
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitle("#Deltar");
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetNdivisions(505);

          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetTitle("");

          // Improve tick marks
          //systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTickLength(0.05);


          // Set the drawing ranges
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetRangeUser(yZoom1, yZoom2);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->GetXaxis()->SetRangeUser(0.006, 0.4);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->GetXaxis()->SetRangeUser(0.006, 0.4);
          systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(0.006, 0.4);

          // Set the drawing style for histograms
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetMarkerSize(1.2);
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerSize(1.2);

          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetMarkerStyle(markerStylePp);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetMarkerSize(1.2);
          energyEnergyCorrelatorSignalPp[eecPpBin]->SetMarkerStyle(markerStylePp);
          energyEnergyCorrelatorSignalPp[eecPpBin]->SetMarkerSize(1.2);

          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);

          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetFillColorAlpha(markerColorPp, 0.4);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetLineColor(markerColorPp);
          energyEnergyCorrelatorSignalPp[eecPpBin]->SetLineColor(markerColorPp);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetMarkerColor(markerColorPp);
          energyEnergyCorrelatorSignalPp[eecPpBin]->SetMarkerColor(markerColorPp);

          // Different drawing styles for correlated uncertainties
          systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->SetFillColorAlpha(markerColorPp, 0.4);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->SetFillColorAlpha(bandColorUpPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->SetFillColorAlpha(bandColorDownPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerSize(0);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerStyle(9);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerSize(0);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerStyle(9);

          // Draw the histograms
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->Draw("e2");
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->Draw("same,e2");
          systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->Draw("same,e3");
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->Draw("same,e3");
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->Draw("same,e3");
          energyEnergyCorrelatorSignalPp[eecPpBin]->Draw("same,p");
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->Draw("same,p");

          // Create a legend for jet pT to each pad
          leftMarginAdder = (iJetPt == firstDrawnJetPtBinEEC[weightExponent-1]) ? leftPadMargin : 0;
          bottomMarginAdder = (iCentrality == lastDrawnCentralityBin[weightExponent-1]) ? bottomPadMargin : 0;
          //legend = new TLegend(0.05 + leftMarginAdder, 0.7 + bottomMarginAdder*0.5, 0.5 / (1 - leftMarginAdder), 0.8 / (thisPadScale*0.81) - bottomMarginAdder*1.4);
          //legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
          
          //legend->Draw();

          // Create a legend for collision systems to each pad
          legend = new TLegend(0.03 + leftMarginAdder, 0.05 + bottomMarginAdder, 0.5 / (1 - leftMarginAdder), 0.45 / thisPadScale);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
          legend->AddEntry(systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin], Form("PbPb %.0f-%.0f%%", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality)), "p");
          legend->AddEntry(systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin], "pp", "p");
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
    std::vector<SplitCanvas*> bigRatioCanvas;
    bigRatioCanvas.resize(nTrackPtBinsEEC);
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

          eecPpBin = std::make_tuple(weightExponent-1, iJetPt, iTrackPt);
          eecPbPbBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt);
          eecPbPbUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeUp, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeDown, iCentrality, iJetPt, iTrackPt);

          canvasIndex = (iCentrality - firstDrawnCentralityBin[weightExponent-1]) * (lastDrawnJetPtBinEEC[weightExponent-1] - firstDrawnJetPtBinEEC[weightExponent-1] + 1) + (iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]);
          bigRatioCanvas[iTrackPt]->CD(canvasIndex);

          gPad->SetLogy(0); // Linear y-axis scale for ratios
          gPad->SetLogx();

          // The titles in the bottow row need to be scaled, since it has more margin than other rows
          thisPadScale = (iCentrality == lastDrawnCentralityBin[weightExponent-1]) ? bottomRowScale : 1;

          // Set the axis titles and labels
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitleOffset(0.85);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitleSize(0.2 * thisPadScale);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetLabelOffset(0.01);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetLabelSize(0.12 * thisPadScale);

          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitleOffset(0.5 / thisPadScale);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitleSize(0.2 * thisPadScale);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetLabelOffset(0.01 / thisPadScale);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetLabelSize(0.12 * thisPadScale);

          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitle("#frac{PbPb}{pp}");
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetStats(0);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetNdivisions(505);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitle("#Deltar");
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetNdivisions(505);

          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetTitle("");

          // Set the drawing ranges
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetRangeUser(0.2, 1.8);
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->GetXaxis()->SetRangeUser(0.006, 0.4);
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->GetXaxis()->SetRangeUser(0.006, 0.4);

          // Set the drawing style for histograms
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetMarkerSize(1.2);
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerSize(1.2);

          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);

          // Different drawing styles for correlated uncertainties
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->SetFillColorAlpha(bandColorUpPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->SetFillColorAlpha(bandColorDownPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerSize(0);
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerStyle(9);
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerSize(0);
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerStyle(9);

          // Draw the histograms
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->Draw("e2");
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->Draw("same,e3");
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->Draw("same,e3");
          energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->Draw("same,p");

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
          legend->AddEntry(systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin], Form("PbPb %.0f-%.0f%% / pp", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality)), "p,e2");
          legend->Draw();

          // Draw a line to one
          oneLine->Draw();

          // Draw the shifted pT to regular pT ratio
          if(drawShiftedPtRatio){
            energyEnergyCorrelatorShiftedPtRatio[eecPpBin]->SetLineStyle(2);
            energyEnergyCorrelatorShiftedPtRatio[eecPpBin]->SetLineColor(kBlack);
            energyEnergyCorrelatorShiftedPtRatio[eecPpBin]->Draw("same,C");
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

      eecDoubleRatioBin = std::make_tuple(weightExponent-1, iCentrality, doubleRatioJetPtBinIndex);

      doubleRatioCanvas->CD(canvasIndex++);

      gPad->SetLogy(false);
      gPad->SetLogx();

      // Set the axis titles and labels
      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetXaxis()->SetTitleOffset(1);
      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetXaxis()->SetTitleSize(0.09);
      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetXaxis()->SetLabelOffset(0.01);
      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetXaxis()->SetLabelSize(0.07);

      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetYaxis()->SetTitleOffset(1.4);
      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetYaxis()->SetTitleSize(0.08);
      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetYaxis()->SetLabelOffset(0.01);
      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetYaxis()->SetLabelSize(0.07);

      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetYaxis()->SetTitle(Form("#frac{PbPb/pp (p_{T}^{ch} > %.0f GeV)}{PbPb/pp (p_{T}^{ch} > %.0f GeV)}", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(trackPtBinsForDoubleRatio.second), card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(trackPtBinsForDoubleRatio.first)));
      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetXaxis()->CenterTitle();
      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->SetStats(0);
      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetYaxis()->SetNdivisions(505);
      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetYaxis()->CenterTitle();
      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetXaxis()->SetTitle("#Deltar");
      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetXaxis()->CenterTitle();
      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetXaxis()->SetNdivisions(505);

      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->SetTitle("");

      // Set the drawing ranges
      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetYaxis()->SetRangeUser(0.8, 1.2);

      // Set the drawing style for histograms
      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->SetMarkerSize(1.2);
      energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
      energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin]->SetMarkerSize(1.2);

      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
      energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
      energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);

      // Draw the histograms
      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->Draw("e2");
      energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin]->Draw("same,p");

      // Show the centrality bin in the legend
      leftMarginAdder = (iCentrality == doubleRatioCentralityBinIndex1) ? leftPadMargin-0.03 : 0;
      legend = new TLegend(0.05+leftMarginAdder, 0.3, 0.45 / (1-leftMarginAdder), 0.5);
      legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.07); legend->SetTextFont(62);
      legend->AddEntry((TObject*)0, Form("%.f < jet p_{T} < %.0f GeV",  card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(doubleRatioJetPtBinIndex), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(doubleRatioJetPtBinIndex)), "");
      legend->AddEntry(systematicUncertaintyDoubleRatio[eecDoubleRatioBin], Form("Centrality: %.0f-%.0f%%", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality)), "p");

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

    // If we are writing files for HepData, open a file for writing
    TFile *doubleRatioHepDataFile;
    if(saveHepDataFiles){
      doubleRatioHepDataFile = TFile::Open("hepdata_energyEnergyCorrelatorDoubleRatio_hin-23-004.root","RECREATE");
    }

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

        eecDoubleRatioBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt);

        // Set the axis titles and labels
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetXaxis()->SetTitleOffset(1);
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetXaxis()->SetTitleSize(0.09);
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetXaxis()->SetLabelOffset(0.01);
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetXaxis()->SetLabelSize(0.06);

        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetYaxis()->SetTitleOffset(1.2);
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetYaxis()->SetTitleSize(0.08);
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetYaxis()->SetLabelOffset(0.01);
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetYaxis()->SetLabelSize(0.06);

        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetYaxis()->SetTitle(Form("#frac{PbPb/pp (p_{T}^{ch} > %.0f GeV)}{PbPb/pp (p_{T}^{ch} > %.0f GeV)}", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(trackPtBinsForDoubleRatio.second), card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(trackPtBinsForDoubleRatio.first)));
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetXaxis()->CenterTitle();
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->SetStats(0);
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetYaxis()->SetNdivisions(505);
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetYaxis()->CenterTitle();
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetXaxis()->SetTitle("#Deltar");
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetXaxis()->CenterTitle();
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetXaxis()->SetNdivisions(505);

        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->SetTitle("");

        // Set the drawing ranges
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);

        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetYaxis()->SetRangeUser(1 - doubleRatioZoomMagnitude, 1 + doubleRatioZoomMagnitude);

        // Set the drawing styles 
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->SetMarkerSize(1.2);
        energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
        energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin]->SetMarkerSize(1.2);

        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
        energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
        energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
      }

      // Draw the histograms
      eecDoubleRatioBin = std::make_tuple(weightExponent-1, firstDrawnCentralityBin[weightExponent-1], iJetPt);
      systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->Draw("e2");

      for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]+1; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
        eecDoubleRatioBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt);
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->Draw("same,e2");
      }
      for(int iCentrality = lastDrawnCentralityBin[weightExponent-1]; iCentrality >= firstDrawnCentralityBin[weightExponent-1]; iCentrality--){
        eecDoubleRatioBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt);
        energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin]->Draw("same,p");
      }

      // Save the double ratio histograms for HepData
      if(saveHepDataFiles){

        for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){

          for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
            eecDoubleRatioBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt);

            // Choose the analysis region for the energy-energy correlator ratios and write them to file
            energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
            energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin]->Write(Form("energyEnergyCorrelatorDoubleRatio%s_C%.0f-%.0f_J%.0f-%.0f", energyWeightString[iWeightExponent].Data(), card[kPbPb][0]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][0]->GetHighBinBorderCentrality(iCentrality), card[kPbPb][0]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][0]->GetHighBinBorderJetPtEEC(iJetPt)));

            systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
            systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->Write(Form("energyEnergyCorrelatorDoubleRatio_systematicUncertainty%s_C%.0f-%.0f_J%.0f-%.0f", energyWeightString[iWeightExponent].Data(), card[kPbPb][0]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][0]->GetHighBinBorderCentrality(iCentrality), card[kPbPb][0]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][0]->GetHighBinBorderJetPtEEC(iJetPt)));

          } // Weight exponent loop
        } // Centrality loop
      } // Saving HepData files

      // Show the centrality bin in the legend
      legend = new TLegend(0.28, 0.22, 0.53, 0.5);
      legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.055); legend->SetTextFont(62);
      for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
        eecDoubleRatioBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt);
        legend->AddEntry(systematicUncertaintyDoubleRatio[eecDoubleRatioBin], Form("Centrality: %.0f-%.0f%%", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality)), "p");
      }

      legend->Draw();

      // Draw the line to one
      oneLine->Draw();

      mainTitle->SetTextFont(62);
      mainTitle->SetTextSize(0.08);

      if(tagSelection){
        mainTitle->DrawLatexNDC(0.29, 0.93, "CMS");

        mainTitle->SetTextFont(52);
        mainTitle->SetTextSize(0.05);
        mainTitle->DrawLatexNDC(0.27, 0.88, tagName[tagSelection]);
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

    // Close the HepData file if it was opened
    if(saveHepDataFiles){
      doubleRatioHepDataFile->Close();
    }

  } // If for drawing double ratios

  // ========================================================================================= //
  // ===== Drawing style with all energy weight exponent distributions in one big canvas ===== //
  // ========================================================================================= //
  
  if(drawBigCanvasAllDistributions){
    
    // Draw all the distributions to big canvases
    std::vector<SplitCanvas*> bigDualDistributionCanvas;
    bigDualDistributionCanvas.resize(nTrackPtBinsEEC);
    mainTitle = new TLatex();

    // If we are writing files for HepData, open a file for writing
    TFile *distributionHepDataFile;
    if(saveHepDataFiles){
      distributionHepDataFile = TFile::Open("hepdata_energyEnergyCorrelatorDistribution_hin-23-004.root","RECREATE");
    }

    for(int iTrackPt = firstDrawnTrackPtBinEEC[0]; iTrackPt <= lastDrawnTrackPtBinEEC[0]; iTrackPt++){
      compactTrackPtString = Form("_T>%.1f", card[kPbPb][0]->GetLowBinBorderTrackPtEEC(iTrackPt));
      compactTrackPtString.ReplaceAll(".","v");

      // Draw a big canvas and put all the plots in it
      bigDualDistributionCanvas[iTrackPt] = new SplitCanvas(Form("bigDualDistributionCanvas%d", iTrackPt), "", 1800, 1500);
      bigDualDistributionCanvas[iTrackPt]->SetMargin(0.08, 0.01, 0.08, 0.11);
      bigDualDistributionCanvas[iTrackPt]->DivideNeatly(4, 4);

      bottomRowScale = bigDualDistributionCanvas[iTrackPt]->GetBottomRowScale();
      bottomPadMargin = bigDualDistributionCanvas[iTrackPt]->GetBottomPadMargin();
      leftColumnScale = bigDualDistributionCanvas[iTrackPt]->GetLeftColumnScale();
      leftPadMargin = bigDualDistributionCanvas[iTrackPt]->GetLeftPadMargin();

      for(int iCentrality = firstDrawnCentralityBin[0]; iCentrality <= lastDrawnCentralityBin[0]; iCentrality++){
        for(int iJetPt = firstDrawnJetPtBinEEC[0]; iJetPt <= lastDrawnJetPtBinEEC[0]; iJetPt++) {
          jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", card[kPbPb][0]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][0]->GetHighBinBorderJetPtEEC(iJetPt));

          canvasIndex = (iCentrality - firstDrawnCentralityBin[0]) * (lastDrawnJetPtBinEEC[0] - firstDrawnJetPtBinEEC[0] + 1) + (iJetPt - firstDrawnJetPtBinEEC[0]);
          bigDualDistributionCanvas[iTrackPt]->CD(canvasIndex);

          gPad->SetLogy();
          gPad->SetLogx();

          // The titles in the bottow row need to be scaled, since it has more margin than other rows. Same for leftmost column
          thisPadScale = (iCentrality == lastDrawnCentralityBin[0]) ? bottomRowScale : 1;
          firstRowScale = (iJetPt == firstDrawnJetPtBinEEC[0]) ? leftColumnScale : 1;

          eecPbPbUncorrelatedSystematicsBin = std::make_tuple(0, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(0, kCorrelatedUncertaintyShapeUp, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(0, kCorrelatedUncertaintyShapeDown, iCentrality, iJetPt, iTrackPt);

          // Set the axis titles and labels
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitleOffset(0.85);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitleSize(0.2 * thisPadScale);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetLabelOffset(0.01);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetLabelSize(0.13 * thisPadScale);

          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitleOffset(0.5 / thisPadScale);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitleSize(0.2 * thisPadScale);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetLabelOffset(0.01 / thisPadScale);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetLabelSize(0.13 * thisPadScale);

          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitle("EEC");
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetStats(0);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetNdivisions(505);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->CenterTitle();
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitle("#Deltar");
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetNdivisions(505);

          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetTitle("");

          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTickLength(0.06 * thisPadScale / firstRowScale);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTickLength(0.05 * firstRowScale / thisPadScale);

          // Set the drawing ranges
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetRangeUser(0.04, 160);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->GetXaxis()->SetRangeUser(0.006, 0.4);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->GetXaxis()->SetRangeUser(0.006, 0.4);

          // Set the drawing style for histograms
          for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
            eecPpBin = std::make_tuple(iWeightExponent, iJetPt, iTrackPt);
            eecPpUncorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kUncorrelatedUncertainty, iJetPt, iTrackPt);
            eecPpCorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kCorrelatedUncertainty, iJetPt, iTrackPt);

            eecPbPbBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt);
            eecPbPbUncorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
            eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(iWeightExponent, kCorrelatedUncertaintyShapeUp, iCentrality, iJetPt, iTrackPt);
            eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(iWeightExponent, kCorrelatedUncertaintyShapeDown, iCentrality, iJetPt, iTrackPt);


            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetMarkerStyle(markerStylePbPb[iWeightExponent]);
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetMarkerSize(1.2);
            energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerStyle(markerStylePbPb[iWeightExponent]);
            energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerSize(1.2);

            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetFillColorAlpha(markerColorPbPb[iWeightExponent], 0.4);
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetLineColor(markerColorPbPb[iWeightExponent]);
            energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetLineColor(markerColorPbPb[iWeightExponent]);
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetMarkerColor(markerColorPbPb[iWeightExponent]);
            energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerColor(markerColorPbPb[iWeightExponent]);

            systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetMarkerStyle(markerStylePpEnergyWeight[iWeightExponent]);
            systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetMarkerSize(1.2);
            energyEnergyCorrelatorSignalPp[eecPpBin]->SetMarkerStyle(markerStylePpEnergyWeight[iWeightExponent]);
            energyEnergyCorrelatorSignalPp[eecPpBin]->SetMarkerSize(1.2);

            systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetFillColorAlpha(markerColorPpEnergyWeight[iWeightExponent], 0.4);
            systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetLineColor(markerColorPpEnergyWeight[iWeightExponent]);
            energyEnergyCorrelatorSignalPp[eecPpBin]->SetLineColor(markerColorPpEnergyWeight[iWeightExponent]);
            systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetMarkerColor(markerColorPpEnergyWeight[iWeightExponent]);
            energyEnergyCorrelatorSignalPp[eecPpBin]->SetMarkerColor(markerColorPpEnergyWeight[iWeightExponent]);          

            // Different drawing styles for correlated uncertainties
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->SetFillColorAlpha(bandColorUpPbPb[iWeightExponent], 0.4);
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->SetFillColorAlpha(bandColorDownPbPb[iWeightExponent], 0.4);
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerSize(0);
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerStyle(9);
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerSize(0);
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerStyle(9);

            systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->SetFillColorAlpha(markerColorPpEnergyWeight[iWeightExponent], 0.4);
          }

          // Draw the histograms
          systematicUncertaintyForPbPb[std::make_tuple(0,kUncorrelatedUncertainty,iCentrality,iJetPt,iTrackPt)]->Draw("e2");
          systematicUncertaintyForPbPb[std::make_tuple(0,kCorrelatedUncertaintyShapeUp,iCentrality,iJetPt,iTrackPt)]->Draw("same,e3");
          systematicUncertaintyForPbPb[std::make_tuple(0,kCorrelatedUncertaintyShapeDown,iCentrality,iJetPt,iTrackPt)]->Draw("same,e3");
          energyEnergyCorrelatorSignalPbPb[std::make_tuple(0,iCentrality,iJetPt,iTrackPt)]->Draw("same,p");
          systematicUncertaintyForPbPb[std::make_tuple(1,kUncorrelatedUncertainty,iCentrality,iJetPt,iTrackPt)]->Draw("same,e2");
          systematicUncertaintyForPbPb[std::make_tuple(1,kCorrelatedUncertaintyShapeUp,iCentrality,iJetPt,iTrackPt)]->Draw("same,e3");
          systematicUncertaintyForPbPb[std::make_tuple(1,kCorrelatedUncertaintyShapeDown,iCentrality,iJetPt,iTrackPt)]->Draw("same,e3");
          energyEnergyCorrelatorSignalPbPb[std::make_tuple(1,iCentrality,iJetPt,iTrackPt)]->Draw("same,p");

          systematicUncertaintyForPp[std::make_tuple(0,kUncorrelatedUncertainty,iJetPt,iTrackPt)]->Draw("same,e2");
          systematicUncertaintyForPp[std::make_tuple(0,kCorrelatedUncertainty,iJetPt,iTrackPt)]->Draw("same,e3");
          energyEnergyCorrelatorSignalPp[std::make_tuple(0,iJetPt,iTrackPt)]->Draw("same,p");
          systematicUncertaintyForPp[std::make_tuple(1,kUncorrelatedUncertainty,iJetPt,iTrackPt)]->Draw("same,e2");
          systematicUncertaintyForPp[std::make_tuple(1,kCorrelatedUncertainty,iJetPt,iTrackPt)]->Draw("same,e3");
          energyEnergyCorrelatorSignalPp[std::make_tuple(1,iJetPt,iTrackPt)]->Draw("same,p");

          // Save the relevant histograms for HepData submission
          if(saveHepDataFiles){

            hepDataTrackPtString = manipulator->StringifyNumber(card[kPbPb][0]->GetLowBinBorderTrackPtEEC(iTrackPt));
            
            // Only save the pp histograms once
            if(iCentrality == firstDrawnCentralityBin[0]){

              for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
                eecPpBin = std::make_tuple(iWeightExponent, iJetPt, iTrackPt);
                eecPpUncorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kUncorrelatedUncertainty, iJetPt, iTrackPt);
                eecPpCorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kCorrelatedUncertainty, iJetPt, iTrackPt);

                // Choose the analysis region for the energy-energy correlator histograms and write them to file
                energyEnergyCorrelatorSignalPp[eecPpBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
                energyEnergyCorrelatorSignalPp[eecPpBin]->Write(Form("energyEnergyCorrelator%s_pp_J%.0f-%.0f_T%s", energyWeightString[iWeightExponent].Data(), card[kPbPb][0]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][0]->GetHighBinBorderJetPtEEC(iJetPt), hepDataTrackPtString.Data()));

                systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
                systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->Write(Form("energyEnergyCorrelator_uncorrelatedUncertainty%s_pp_J%.0f-%.0f_T%s", energyWeightString[iWeightExponent].Data(), card[kPbPb][0]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][0]->GetHighBinBorderJetPtEEC(iJetPt), hepDataTrackPtString.Data()));

                systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
                systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->Write(Form("energyEnergyCorrelator_correlatedUncertainty%s_pp_J%.0f-%.0f_T%s", energyWeightString[iWeightExponent].Data(), card[kPbPb][0]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][0]->GetHighBinBorderJetPtEEC(iJetPt), hepDataTrackPtString.Data()));
              }
            }

            // Save also the PbPb histograms
            for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
              eecPbPbBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt);
              eecPbPbUncorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
              eecPbPbCorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kCorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);

              // Choose the analysis region for the energy-energy correlator histograms and write them to file
              energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
              energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->Write(Form("energyEnergyCorrelator%s_C%.0f-%.0f_J%.0f-%.0f_T%s", energyWeightString[iWeightExponent].Data(), card[kPbPb][0]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][0]->GetHighBinBorderCentrality(iCentrality), card[kPbPb][0]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][0]->GetHighBinBorderJetPtEEC(iJetPt), hepDataTrackPtString.Data()));

              systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
              systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->Write(Form("energyEnergyCorrelator_uncorrelatedUncertainty%s_C%.0f-%.0f_J%.0f-%.0f_T%s", energyWeightString[iWeightExponent].Data(), card[kPbPb][0]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][0]->GetHighBinBorderCentrality(iCentrality), card[kPbPb][0]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][0]->GetHighBinBorderJetPtEEC(iJetPt), hepDataTrackPtString.Data()));

              systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
              systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBin]->Write(Form("energyEnergyCorrelator_correlatedUncertainty%s_C%.0f-%.0f_J%.0f-%.0f_T%s", energyWeightString[iWeightExponent].Data(), card[kPbPb][0]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][0]->GetHighBinBorderCentrality(iCentrality), card[kPbPb][0]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][0]->GetHighBinBorderJetPtEEC(iJetPt), hepDataTrackPtString.Data()));
            }

          } // End for saving HepData files

          // Legend adders to put the legends to correct positions
          leftMarginAdder = (iJetPt == firstDrawnJetPtBinEEC[0]) ? leftPadMargin : 0;
          bottomMarginAdder = (iCentrality == lastDrawnCentralityBin[0]) ? bottomPadMargin : 0;

          // Add the jet pT binning to the top row legend
          if(iCentrality == firstDrawnCentralityBin[0]){
            legend = new TLegend(0.01 + leftMarginAdder/1.1, 0.85, 0.5 / (1 - leftMarginAdder), 0.95);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.12 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry((TObject*)0, Form("%s", jetPtString.Data()), "");
            legend->Draw();

          }

          // Add the centrality bin information to the leftmost column
          if(iJetPt == firstDrawnJetPtBinEEC[0]){
            legend = new TLegend(0.19, 0.085 + bottomMarginAdder*1.05, 0.7, 0.2 / (thisPadScale*0.78) + bottomMarginAdder*0.1 + 0.05);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.12 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry((TObject*)0, Form("PbPb %.0f-%.0f%%", card[kPbPb][0]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][0]->GetHighBinBorderCentrality(iCentrality)), "");
            legend->Draw();
          }

          // Add information about the nominal weight PbPb histogram to selected pad
          if(iJetPt == firstDrawnJetPtBinEEC[0]+1 && iCentrality == firstDrawnCentralityBin[0]+1){
            eecPbPbUncorrelatedSystematicsBin = std::make_tuple(0, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
            legend = new TLegend(0.4, 0.83, 0.8, 0.93);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.12 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry(systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin], "PbPb, n=1", "lpf");
            legend->Draw();
          }

          // Add information about the squared weight PbPb histogram to selected pad
          if(iJetPt == firstDrawnJetPtBinEEC[0]+2 && iCentrality == firstDrawnCentralityBin[0]+1){
            eecPbPbUncorrelatedSystematicsBin = std::make_tuple(1, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
            legend = new TLegend(0.4, 0.83, 0.8, 0.93);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.12 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry(systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin], "PbPb, n=2", "lpf");
            legend->Draw();
          }

          // Add information about the nominal weight PbPb histogram to selected pad
          if(iJetPt == firstDrawnJetPtBinEEC[0]+1 && iCentrality == firstDrawnCentralityBin[0]+2){
            eecPpUncorrelatedSystematicsBin = std::make_tuple(0, kUncorrelatedUncertainty, iJetPt, iTrackPt);
            legend = new TLegend(0.4, 0.83, 0.8, 0.93);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.12 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry(systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin], "pp, n=1", "lpf");
            legend->Draw();
          }

          // Add information about the squared weight PbPb histogram to selected pad
          if(iJetPt == firstDrawnJetPtBinEEC[0]+2 && iCentrality == firstDrawnCentralityBin[0]+2){
            eecPpUncorrelatedSystematicsBin = std::make_tuple(1, kUncorrelatedUncertainty, iJetPt, iTrackPt);
            legend = new TLegend(0.4, 0.83, 0.8, 0.93);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.12 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry(systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin], "pp, n=2", "lpf");
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

      if(tagSelection){
        mainTitle->SetTextSize(0.05);
        mainTitle->DrawLatexNDC(0.08, 0.945, "CMS");

        mainTitle->SetTextFont(52);
        mainTitle->SetTextSize(0.035);
        mainTitle->DrawLatexNDC(0.06, 0.91, tagName[tagSelection]);
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

    // Close the HepData file if it was opened
    if(saveHepDataFiles){
      distributionHepDataFile->Close();
    }

  } // If for drawing big canvases

  // ======================================================================================== //
  // ===== Drawing style with ratios from all energy weight exponents in one big canvas ===== //
  // ======================================================================================== //
  
  if(drawBigCanvasAllRatios){
    
    // Draw all the distributions to big canvases
    std::vector<SplitCanvas*> bigDualRatioCanvas;
    bigDualRatioCanvas.resize(nTrackPtBinsEEC);
    mainTitle = new TLatex();

    // If we are writing files for HepData, open a file for writing
    TFile *ratioHepDataFile;
    if(saveHepDataFiles){
      ratioHepDataFile = TFile::Open("hepdata_energyEnergyCorrelatorRatio_hin-23-004.root","RECREATE");
    }

    for(int iTrackPt = firstDrawnTrackPtBinEEC[0]; iTrackPt <= lastDrawnTrackPtBinEEC[0]; iTrackPt++){
      compactTrackPtString = Form("_T>%.1f", card[kPbPb][0]->GetLowBinBorderTrackPtEEC(iTrackPt));
      compactTrackPtString.ReplaceAll(".","v");

      // Draw a big canvas and put all the plots in it
      bigDualRatioCanvas[iTrackPt] = new SplitCanvas(Form("bigDualRatioCanvas%d", iTrackPt), "", 1800, 1500);
      bigDualRatioCanvas[iTrackPt]->SetMargin(0.08, 0.01, 0.08, 0.11);
      bigDualRatioCanvas[iTrackPt]->DivideNeatly(4, 4);

      bottomRowScale = bigDualRatioCanvas[iTrackPt]->GetBottomRowScale();
      bottomPadMargin = bigDualRatioCanvas[iTrackPt]->GetBottomPadMargin();
      leftColumnScale = bigDualRatioCanvas[iTrackPt]->GetLeftColumnScale();
      leftPadMargin = bigDualRatioCanvas[iTrackPt]->GetLeftPadMargin();

      for(int iCentrality = firstDrawnCentralityBin[0]; iCentrality <= lastDrawnCentralityBin[0]; iCentrality++){
        for(int iJetPt = firstDrawnJetPtBinEEC[0]; iJetPt <= lastDrawnJetPtBinEEC[0]; iJetPt++) {
          jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", card[kPbPb][0]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][0]->GetHighBinBorderJetPtEEC(iJetPt));

          eecPbPbUncorrelatedSystematicsBin = std::make_tuple(0, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(0, kCorrelatedUncertaintyShapeUp, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(0, kCorrelatedUncertaintyShapeDown, iCentrality, iJetPt, iTrackPt);

          canvasIndex = (iCentrality - firstDrawnCentralityBin[0]) * (lastDrawnJetPtBinEEC[0] - firstDrawnJetPtBinEEC[0] + 1) + (iJetPt - firstDrawnJetPtBinEEC[0]);
          bigDualRatioCanvas[iTrackPt]->CD(canvasIndex);

          gPad->SetLogy(0); // Linear y-axis scale for ratios
          gPad->SetLogx();

          // The titles in the bottow row need to be scaled, since it has more margin than other rows. Same for leftmost column
          thisPadScale = (iCentrality == lastDrawnCentralityBin[0]) ? bottomRowScale : 1;
          firstRowScale = (iJetPt == firstDrawnJetPtBinEEC[0]) ? leftColumnScale : 1;

          // Set the axis titles and labels
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitleOffset(0.85);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitleSize(0.2 * thisPadScale);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetLabelOffset(0.01);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetLabelSize(0.13 * thisPadScale);

          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitleOffset(0.5 / thisPadScale);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitleSize(0.2 * thisPadScale);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetLabelOffset(0.01 / thisPadScale);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetLabelSize(0.13 * thisPadScale);

          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitle("#frac{PbPb}{pp}");
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetStats(0);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetNdivisions(505);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitle("#Deltar");
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetNdivisions(505);

          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetTitle("");

          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTickLength(0.06 * thisPadScale / firstRowScale);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTickLength(0.05 * firstRowScale / thisPadScale);

          // Set the drawing ranges
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetRangeUser(0.2, 1.8);
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->GetXaxis()->SetRangeUser(0.006, 0.4);
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->GetXaxis()->SetRangeUser(0.006, 0.4);

          // Set the drawing style for histograms
          for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
            eecPbPbBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt);
            eecPbPbUncorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
            eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(iWeightExponent, kCorrelatedUncertaintyShapeUp, iCentrality, iJetPt, iTrackPt);
            eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(iWeightExponent, kCorrelatedUncertaintyShapeDown, iCentrality, iJetPt, iTrackPt);

            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetMarkerStyle(markerStylePbPb[iWeightExponent]);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetMarkerSize(1.2);
            energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerStyle(markerStylePbPb[iWeightExponent]);
            energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerSize(1.2);

            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetFillColorAlpha(markerColorPbPb[iWeightExponent], 0.4);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetLineColor(markerColorPbPb[iWeightExponent]);
            energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->SetLineColor(markerColorPbPb[iWeightExponent]);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetMarkerColor(markerColorPbPb[iWeightExponent]);
            energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->SetMarkerColor(markerColorPbPb[iWeightExponent]);

            // Different drawing styles for correlated uncertainties
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->SetFillColorAlpha(bandColorUpPbPb[iWeightExponent], 0.4);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->SetFillColorAlpha(bandColorDownPbPb[iWeightExponent], 0.4);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerSize(0);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerStyle(9);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerSize(0);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerStyle(9);
          }

          // Draw the histograms
          systematicUncertaintyPbPbToPpRatio[std::make_tuple(0,kUncorrelatedUncertainty,iCentrality,iJetPt,iTrackPt)]->Draw("e2");
          systematicUncertaintyPbPbToPpRatio[std::make_tuple(0,kCorrelatedUncertaintyShapeUp,iCentrality,iJetPt,iTrackPt)]->Draw("same,e3");
          systematicUncertaintyPbPbToPpRatio[std::make_tuple(0,kCorrelatedUncertaintyShapeDown,iCentrality,iJetPt,iTrackPt)]->Draw("same,e3");
          energyEnergyCorrelatorPbPbToPpRatio[std::make_tuple(0,iCentrality,iJetPt,iTrackPt)]->Draw("same,p");
          systematicUncertaintyPbPbToPpRatio[std::make_tuple(1,kUncorrelatedUncertainty,iCentrality,iJetPt,iTrackPt)]->Draw("same,e2");
          systematicUncertaintyPbPbToPpRatio[std::make_tuple(1,kCorrelatedUncertaintyShapeUp,iCentrality,iJetPt,iTrackPt)]->Draw("same,e3");
          systematicUncertaintyPbPbToPpRatio[std::make_tuple(1,kCorrelatedUncertaintyShapeDown,iCentrality,iJetPt,iTrackPt)]->Draw("same,e3");
          energyEnergyCorrelatorPbPbToPpRatio[std::make_tuple(1,iCentrality,iJetPt,iTrackPt)]->Draw("same,p");

          // Save the ratio histograms for HepData
          if(saveHepDataFiles){

            hepDataTrackPtString = manipulator->StringifyNumber(card[kPbPb][0]->GetLowBinBorderTrackPtEEC(iTrackPt));

            // Save also the PbPb histograms
            for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
              eecPbPbBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt);
              eecPbPbUncorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
              eecPbPbCorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kCorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);

              // Choose the analysis region for the energy-energy correlator ratios and write them to file
              energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
              energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->Write(Form("energyEnergyCorrelatorRatio%s_C%.0f-%.0f_J%.0f-%.0f_T%s", energyWeightString[iWeightExponent].Data(), card[kPbPb][0]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][0]->GetHighBinBorderCentrality(iCentrality), card[kPbPb][0]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][0]->GetHighBinBorderJetPtEEC(iJetPt), hepDataTrackPtString.Data()));

              systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
              systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->Write(Form("energyEnergyCorrelatorRatio_uncorrelatedUncertainty%s_C%.0f-%.0f_J%.0f-%.0f_T%s", energyWeightString[iWeightExponent].Data(), card[kPbPb][0]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][0]->GetHighBinBorderCentrality(iCentrality), card[kPbPb][0]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][0]->GetHighBinBorderJetPtEEC(iJetPt), hepDataTrackPtString.Data()));

              systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
              systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBin]->Write(Form("energyEnergyCorrelatorRatio_correlatedUncertainty%s_C%.0f-%.0f_J%.0f-%.0f_T%s", energyWeightString[iWeightExponent].Data(), card[kPbPb][0]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][0]->GetHighBinBorderCentrality(iCentrality), card[kPbPb][0]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][0]->GetHighBinBorderJetPtEEC(iJetPt), hepDataTrackPtString.Data()));
            }

          }

          // Legend adders to put the legends to correct positions
          leftMarginAdder = (iJetPt == firstDrawnJetPtBinEEC[0]) ? leftPadMargin : 0;
          bottomMarginAdder = (iCentrality == lastDrawnCentralityBin[0]) ? bottomPadMargin : 0;

          // Add the jet pT binning to the top row legend
          if(iCentrality == firstDrawnCentralityBin[0]){
            legend = new TLegend(0.01 + leftMarginAdder/1.1, 0.83, 0.5 / (1 - leftMarginAdder), 0.93);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.12 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry((TObject*)0, Form("%s", jetPtString.Data()), "");
            legend->Draw();

          }

          // Add the centrality bin information to the leftmost column
          if(iJetPt == firstDrawnJetPtBinEEC[0]){
            legend = new TLegend(0.19, 0.05 + bottomMarginAdder*1.1, 0.7, 0.2 / (thisPadScale*0.78) + bottomMarginAdder*0.1);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.12 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry((TObject*)0, Form("#frac{PbPb %.0f-%.0f%%}{pp}", card[kPbPb][0]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][0]->GetHighBinBorderCentrality(iCentrality)), "");
            legend->Draw();
          }

          // Add information about the nominal weight histogram to selected oad
          if(iJetPt == firstDrawnJetPtBinEEC[0]+1 && iCentrality == firstDrawnCentralityBin[0]+1){
            eecPbPbUncorrelatedSystematicsBin = std::make_tuple(0, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
            legend = new TLegend(0.1, 0.83, 0.5, 0.93);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.12 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry(systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin], "EEC ratio with n=1", "lpf");
            legend->Draw();
          }

          // Add information about the squared weight histogram to selected oad
          if(iJetPt == firstDrawnJetPtBinEEC[0]+2 && iCentrality == firstDrawnCentralityBin[0]+1){
            eecPbPbUncorrelatedSystematicsBin = std::make_tuple(1, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
            legend = new TLegend(0.1, 0.83, 0.5, 0.93);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.12 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry(systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin], "EEC ratio with n=2", "lpf");
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

      if(tagSelection){
        mainTitle->SetTextSize(0.05);
        mainTitle->DrawLatexNDC(0.08, 0.945, "CMS");

        mainTitle->SetTextFont(52);
        mainTitle->SetTextSize(0.035);
        mainTitle->DrawLatexNDC(0.06, 0.91, tagName[tagSelection]);
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

    // Close the HepData file if it was opened
    if(saveHepDataFiles){
      ratioHepDataFile->Close();
    }

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

          eecPbPbUncorrelatedSystematicsBin = std::make_tuple(0, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(0, kCorrelatedUncertaintyShapeUp, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(0, kCorrelatedUncertaintyShapeDown, iCentrality, iJetPt, iTrackPt);
          
          letterPaperDistributionCanvas->CD(canvasIndex++);

          gPad->SetLogy();
          gPad->SetLogx();

          // The titles in the bottow row need to be scaled, since it has more margin than other rows
          thisPadScale = (canvasIndex > 1) ? bottomRowScale : 1;

          // Set the axis titles and labels
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitleOffset(0.65);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitleSize(0.15 * thisPadScale);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetLabelOffset(0.0001);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetLabelSize(0.07 * thisPadScale);

          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitleOffset(0.6 / thisPadScale);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitleSize(0.15 * thisPadScale);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetLabelOffset(0.01 / thisPadScale);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetLabelSize(0.07 * thisPadScale);

          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitle("EEC");
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetStats(0);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetNdivisions(505);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->CenterTitle();
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitle("#Deltar");
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetNdivisions(505);

          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetTitle("");

          // Set the drawing ranges
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetRangeUser(0.04, 160);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->GetXaxis()->SetRangeUser(0.006, 0.4);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->GetXaxis()->SetRangeUser(0.006, 0.4);

          // Set the drawing style for histograms
          for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
            eecPpBin = std::make_tuple(iWeightExponent, iJetPt, iTrackPt);
            eecPpUncorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kUncorrelatedUncertainty, iJetPt, iTrackPt);
            eecPpCorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kCorrelatedUncertainty, iJetPt, iTrackPt);

            eecPbPbBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt);
            eecPbPbUncorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
            eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(iWeightExponent, kCorrelatedUncertaintyShapeUp, iCentrality, iJetPt, iTrackPt);
            eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(iWeightExponent, kCorrelatedUncertaintyShapeDown, iCentrality, iJetPt, iTrackPt);

            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetMarkerStyle(markerStylePbPb[iWeightExponent]);
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetMarkerSize(1.2);
            energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerStyle(markerStylePbPb[iWeightExponent]);
            energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerSize(1.2);

            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetFillColorAlpha(markerColorPbPb[iWeightExponent], 0.4);
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetLineColor(markerColorPbPb[iWeightExponent]);
            energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetLineColor(markerColorPbPb[iWeightExponent]);
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetMarkerColor(markerColorPbPb[iWeightExponent]);
            energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerColor(markerColorPbPb[iWeightExponent]);

            systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetMarkerStyle(markerStylePpEnergyWeight[iWeightExponent]);
            systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetMarkerSize(1.2);
            energyEnergyCorrelatorSignalPp[eecPpBin]->SetMarkerStyle(markerStylePpEnergyWeight[iWeightExponent]);
            energyEnergyCorrelatorSignalPp[eecPpBin]->SetMarkerSize(1.2);

            systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetFillColorAlpha(markerColorPpEnergyWeight[iWeightExponent], 0.4);
            systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetLineColor(markerColorPpEnergyWeight[iWeightExponent]);
            energyEnergyCorrelatorSignalPp[eecPpBin]->SetLineColor(markerColorPpEnergyWeight[iWeightExponent]);
            systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetMarkerColor(markerColorPpEnergyWeight[iWeightExponent]);
            energyEnergyCorrelatorSignalPp[eecPpBin]->SetMarkerColor(markerColorPpEnergyWeight[iWeightExponent]);          

            // Different drawing styles for correlated uncertainties
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->SetFillColorAlpha(bandColorUpPbPb[iWeightExponent], 0.4);
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->SetFillColorAlpha(bandColorDownPbPb[iWeightExponent], 0.4);
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerSize(0);
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerStyle(9);
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerSize(0);
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerStyle(9);

            systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->SetFillColorAlpha(markerColorPpEnergyWeight[iWeightExponent], 0.4);
          }

          // Draw the histograms
          systematicUncertaintyForPbPb[std::make_tuple(0,kUncorrelatedUncertainty,iCentrality,iJetPt,iTrackPt)]->Draw("e2");
          systematicUncertaintyForPbPb[std::make_tuple(0,kCorrelatedUncertaintyShapeUp,iCentrality,iJetPt,iTrackPt)]->Draw("same,e3");
          systematicUncertaintyForPbPb[std::make_tuple(0,kCorrelatedUncertaintyShapeDown,iCentrality,iJetPt,iTrackPt)]->Draw("same,e3");
          energyEnergyCorrelatorSignalPbPb[std::make_tuple(0,iCentrality,iJetPt,iTrackPt)]->Draw("same,p");
          systematicUncertaintyForPbPb[std::make_tuple(1,kUncorrelatedUncertainty,iCentrality,iJetPt,iTrackPt)]->Draw("same,e2");
          systematicUncertaintyForPbPb[std::make_tuple(1,kCorrelatedUncertaintyShapeUp,iCentrality,iJetPt,iTrackPt)]->Draw("same,e3");
          systematicUncertaintyForPbPb[std::make_tuple(1,kCorrelatedUncertaintyShapeDown,iCentrality,iJetPt,iTrackPt)]->Draw("same,e3");
          energyEnergyCorrelatorSignalPbPb[std::make_tuple(1,iCentrality,iJetPt,iTrackPt)]->Draw("same,p");

          systematicUncertaintyForPp[std::make_tuple(0,kUncorrelatedUncertainty,iJetPt,iTrackPt)]->Draw("same,e2");
          systematicUncertaintyForPp[std::make_tuple(0,kCorrelatedUncertainty,iJetPt,iTrackPt)]->Draw("same,e3");
          energyEnergyCorrelatorSignalPp[std::make_tuple(0,iJetPt,iTrackPt)]->Draw("same,p");
          systematicUncertaintyForPp[std::make_tuple(1,kUncorrelatedUncertainty,iJetPt,iTrackPt)]->Draw("same,e2");
          systematicUncertaintyForPp[std::make_tuple(1,kCorrelatedUncertainty,iJetPt,iTrackPt)]->Draw("same,e3");
          energyEnergyCorrelatorSignalPp[std::make_tuple(1,iJetPt,iTrackPt)]->Draw("same,p");

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
            anotherLegend->AddEntry(systematicUncertaintyForPbPb[std::make_tuple(0,kUncorrelatedUncertainty,iCentrality,iJetPt,iTrackPt)], "PbPb, n=1", "lpf");
            anotherLegend->AddEntry(systematicUncertaintyForPbPb[std::make_tuple(1,kUncorrelatedUncertainty,iCentrality,iJetPt,iTrackPt)], "PbPb, n=2", "lpf");
            anotherLegend->AddEntry(systematicUncertaintyForPp[std::make_tuple(0,kUncorrelatedUncertainty,iJetPt,iTrackPt)], "pp, n=1", "lpf");
            anotherLegend->AddEntry(systematicUncertaintyForPp[std::make_tuple(1,kUncorrelatedUncertainty,iJetPt,iTrackPt)], "pp, n=2", "lpf");
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

          eecPbPbUncorrelatedSystematicsBin = std::make_tuple(0, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(0, kCorrelatedUncertaintyShapeUp, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(0, kCorrelatedUncertaintyShapeDown, iCentrality, iJetPt, iTrackPt);
          
          letterPaperRatioCanvas->CD(canvasIndex++);

          gPad->SetLogy(0); // Linear y-axis scale for ratios
          gPad->SetLogx();

          // The titles in the bottow row need to be scaled, since it has more margin than other rows
          thisPadScale = (canvasIndex == 3) ? bottomRowScale : 1;


          // Set the axis titles and labels
          if(canvasIndex == 4){
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitleOffset(0.66);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetLabelOffset(-0.027);
          } else {
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitleOffset(0.88);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetLabelOffset(0.001);
          }
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitleSize(0.13 * thisPadScale);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetLabelSize(0.1 * thisPadScale);

          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitleOffset(0.63 / thisPadScale);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitleSize(0.13 * thisPadScale);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetLabelOffset(0.01 / thisPadScale);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetLabelSize(0.1 * thisPadScale);

          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitle("#frac{PbPb}{pp}");
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetStats(0);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetNdivisions(505);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitle("#Deltar");
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetNdivisions(505);

          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetTitle("");

          // Set the drawing ranges
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetRangeUser(0.2, 1.8);
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->GetXaxis()->SetRangeUser(0.006, 0.4);
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->GetXaxis()->SetRangeUser(0.006, 0.4);

          // Set the drawing style for histograms
          for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
            eecPbPbBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt);
            eecPbPbUncorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
            eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(iWeightExponent, kCorrelatedUncertaintyShapeUp, iCentrality, iJetPt, iTrackPt);
            eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(iWeightExponent, kCorrelatedUncertaintyShapeDown, iCentrality, iJetPt, iTrackPt);

            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetMarkerStyle(markerStylePbPb[iWeightExponent]);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetMarkerSize(1.2);
            energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerStyle(markerStylePbPb[iWeightExponent]);
            energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerSize(1.2);

            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetFillColorAlpha(markerColorPbPb[iWeightExponent], 0.4);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetLineColor(markerColorPbPb[iWeightExponent]);
            energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->SetLineColor(markerColorPbPb[iWeightExponent]);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetMarkerColor(markerColorPbPb[iWeightExponent]);
            energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->SetMarkerColor(markerColorPbPb[iWeightExponent]);

            // Different drawing styles for correlated uncertainties
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->SetFillColorAlpha(bandColorUpPbPb[iWeightExponent], 0.4);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->SetFillColorAlpha(bandColorDownPbPb[iWeightExponent], 0.4);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerSize(0);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerStyle(9);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerSize(0);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerStyle(9);
          }

          // Draw the histograms
          systematicUncertaintyPbPbToPpRatio[std::make_tuple(0,kUncorrelatedUncertainty,iCentrality,iJetPt,iTrackPt)]->Draw("e2");
          systematicUncertaintyPbPbToPpRatio[std::make_tuple(0,kCorrelatedUncertaintyShapeUp,iCentrality,iJetPt,iTrackPt)]->Draw("same,e3");
          systematicUncertaintyPbPbToPpRatio[std::make_tuple(0,kCorrelatedUncertaintyShapeDown,iCentrality,iJetPt,iTrackPt)]->Draw("same,e3");
          energyEnergyCorrelatorPbPbToPpRatio[std::make_tuple(0,iCentrality,iJetPt,iTrackPt)]->Draw("same,p");
          systematicUncertaintyPbPbToPpRatio[std::make_tuple(1,kUncorrelatedUncertainty,iCentrality,iJetPt,iTrackPt)]->Draw("same,e2");
          systematicUncertaintyPbPbToPpRatio[std::make_tuple(1,kCorrelatedUncertaintyShapeUp,iCentrality,iJetPt,iTrackPt)]->Draw("same,e3");
          systematicUncertaintyPbPbToPpRatio[std::make_tuple(1,kCorrelatedUncertaintyShapeDown,iCentrality,iJetPt,iTrackPt)]->Draw("same,e3");
          energyEnergyCorrelatorPbPbToPpRatio[std::make_tuple(1,iCentrality,iJetPt,iTrackPt)]->Draw("same,p");

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

            eecPbPbUncorrelatedSystematicsBin = std::make_tuple(0, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
            anotherLegend = new TLegend(0.28, 0.05, 0.88, 0.15);
            anotherLegend->SetFillStyle(0); anotherLegend->SetBorderSize(0); anotherLegend->SetTextSize(0.08 * thisPadScale); anotherLegend->SetTextFont(62);
            anotherLegend->AddEntry(systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin], "EEC ratio with n=1", "lpf");
            anotherLegend->Draw();
          }

          if(canvasIndex == 2){
            legend = new TLegend(0.4, 0.75, 0.7, 0.95);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.08 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry((TObject*)0, centralityString.Data(), "");
            legend->AddEntry((TObject*)0, trackPtString.Data(), "");
            legend->Draw();

            eecPbPbUncorrelatedSystematicsBin = std::make_tuple(1, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
            anotherLegend = new TLegend(0.08, 0.05, 0.83, 0.15);
            anotherLegend->SetFillStyle(0); anotherLegend->SetBorderSize(0); anotherLegend->SetTextSize(0.08 * thisPadScale); anotherLegend->SetTextFont(62);
            anotherLegend->AddEntry(systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin], "EEC ratio with n=2", "lpf");
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

      eecPpUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, firstJetPtIndex, iTrackPt);

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
      systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetTitleOffset(0.85);
      systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetTitleSize(0.08);
      systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetLabelOffset(0.001);
      systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetLabelSize(0.07);

      systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->SetTitleOffset(0.9);
      systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->SetTitleSize(0.08);
      systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->SetLabelOffset(0.01);
      systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->SetLabelSize(0.07);

      systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->SetTitle("EEC");
      systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
      systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetStats(0);
      systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->SetNdivisions(505);
      systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->CenterTitle();
      systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetTitle("#Deltar");
      systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
      systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetNdivisions(505);

      systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetTitle("");

      // Set the drawing ranges
      for(int iJetPt = firstDrawnJetPtBinEEC[weightExponent-1]; iJetPt <= lastDrawnJetPtBinEEC[weightExponent-1]; iJetPt++){
        eecPpBin = std::make_tuple(weightExponent-1, iJetPt, iTrackPt);
        eecPpUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iJetPt, iTrackPt);
        eecPpCorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kCorrelatedUncertainty, iJetPt, iTrackPt);

        systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
        systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->SetRangeUser(supplementaryDistributionZoomLow[weightExponent-1], supplementaryDistributionZoomHigh[weightExponent-1]);
        systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(0.006, 0.4);

        // Set the drawing style for histograms
        systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetMarkerStyle(markerStylePbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);
        systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetMarkerSize(1.2);
        energyEnergyCorrelatorSignalPp[eecPpBin]->SetMarkerStyle(markerStylePbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);
        energyEnergyCorrelatorSignalPp[eecPpBin]->SetMarkerSize(1.2);

        systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetFillColorAlpha(markerColorPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]], 0.4);
        systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetLineColor(markerColorPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);
        energyEnergyCorrelatorSignalPp[eecPpBin]->SetLineColor(markerColorPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);
        systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetMarkerColor(markerColorPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);
        energyEnergyCorrelatorSignalPp[eecPpBin]->SetMarkerColor(markerColorPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);

        // Different drawing styles for correlated uncertainties
        if(iJetPt % 2 == 0){
          systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->SetFillColorAlpha(bandColorDownPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]], 0.4);
        } else {
          systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->SetFillColorAlpha(bandColorUpPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]], 0.4);
        }
        systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->SetMarkerSize(0);
        systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->SetMarkerStyle(9);

        // Draw the systematic uncertainties
        if(iJetPt == firstDrawnJetPtBinEEC[weightExponent-1]){
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->Draw("e2");
        } else {
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->Draw("same,e2");
        }

        systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->Draw("same,e3");
      }

      // Make a legend for the jet pT bins
      legend = new TLegend(0.18, 0.18, 0.43, 0.46);
      legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(supplementaryLegendTextSize); legend->SetTextFont(62);

      // Draw the data points on top of the error bands
      for(int iJetPt = firstDrawnJetPtBinEEC[weightExponent-1]; iJetPt <= lastDrawnJetPtBinEEC[weightExponent-1]; iJetPt++){
        eecPpBin = std::make_tuple(weightExponent-1, iJetPt, iTrackPt);
        eecPpUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iJetPt, iTrackPt);
        energyEnergyCorrelatorSignalPp[eecPpBin]->Draw("same,p");
        legend->AddEntry(systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin], Form("%.0f < p_{T,jet} < %.0f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt)), "lpf");
      }
 
      legend->Draw();

      mainTitle->SetTextFont(62);
      mainTitle->SetTextSize(0.06);

      if(tagSelection){
        mainTitle->DrawLatexNDC(0.61+tagShift[tagSelection], 0.85, "CMS");

        mainTitle->SetTextFont(52);
        mainTitle->SetTextSize(0.045);
        mainTitle->DrawLatexNDC(0.745+tagShift[tagSelection], 0.85, tagName[tagSelection]);
      } else {
        mainTitle->DrawLatexNDC(0.78, 0.85, "CMS");
      }

      mainTitle->SetTextFont(42);
      mainTitle->SetTextSize(0.04);
      mainTitle->DrawLatexNDC(0.626, 0.95, "302 pb^{-1} pp (5.02 TeV)");

      mainTitle->SetTextFont(62);
      mainTitle->SetTextSize(supplementaryLegendTextSize);
      mainTitle->DrawLatexNDC(0.81, 0.78, Form("pp, n=%d", weightExponent));
      mainTitle->DrawLatexNDC(0.736, 0.71, trackPtString.Data());

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

        eecPbPbUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iCentrality, firstJetPtIndex, iTrackPt);
        
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
        systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitleOffset(0.85);
        systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitleSize(0.08);
        systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetLabelOffset(0.001);
        systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetLabelSize(0.07);

        systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitleOffset(0.9);
        systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitleSize(0.08);
        systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetLabelOffset(0.01);
        systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetLabelSize(0.07);

        systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitle("EEC");
        systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
        systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetStats(0);
        systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetNdivisions(505);
        systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->CenterTitle();
        systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitle("#Deltar");
        systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
        systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetNdivisions(505);

        systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetTitle("");

        // Set the drawing ranges
        for(int iJetPt = firstDrawnJetPtBinEEC[weightExponent-1]; iJetPt <= lastDrawnJetPtBinEEC[weightExponent-1]; iJetPt++){
          eecPbPbBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt);
          eecPbPbUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeUp, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeDown, iCentrality, iJetPt, iTrackPt);

          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetRangeUser(supplementaryDistributionZoomLow[weightExponent-1], supplementaryDistributionZoomHigh[weightExponent-1]);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->GetXaxis()->SetRangeUser(0.006, 0.4);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->GetXaxis()->SetRangeUser(0.006, 0.4);

          // Set the drawing style for histograms
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetMarkerStyle(markerStylePbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetMarkerSize(1.2);
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerStyle(markerStylePbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerSize(1.2);


          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetFillColorAlpha(markerColorPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]], 0.4);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetLineColor(markerColorPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetLineColor(markerColorPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetMarkerColor(markerColorPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerColor(markerColorPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]]);

          // Different drawing styles for correlated uncertainties
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->SetFillColorAlpha(bandColorUpPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]], 0.4);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->SetFillColorAlpha(bandColorDownPbPb[iJetPt - firstDrawnJetPtBinEEC[weightExponent-1]], 0.4);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerSize(0);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerStyle(9);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerSize(0);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerStyle(9);

          // Draw the systematic uncertainties
          if(iJetPt == firstDrawnJetPtBinEEC[weightExponent-1]){
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->Draw("e2");
          } else {
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->Draw("same,e2");
          }

          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->Draw("same,e3");
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->Draw("same,e3");
        }

        // Make a legend for the jet pT bins
        legend = new TLegend(0.18, 0.18, 0.43, 0.46);
        legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(supplementaryLegendTextSize); legend->SetTextFont(62);

        // Draw the data points on top of the error bands
        for(int iJetPt = firstDrawnJetPtBinEEC[weightExponent-1]; iJetPt <= lastDrawnJetPtBinEEC[weightExponent-1]; iJetPt++){
          eecPbPbBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt);
          eecPbPbUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->Draw("same,p");
          legend->AddEntry(systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin], Form("%.0f < p_{T,jet} < %.0f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt)), "lpf");
        }
          

        legend->Draw();

        mainTitle->SetTextFont(62);
        mainTitle->SetTextSize(0.06);

        if(tagSelection){
          mainTitle->DrawLatexNDC(0.61+tagShift[tagSelection], 0.85, "CMS");

          mainTitle->SetTextFont(52);
          mainTitle->SetTextSize(0.045);
          mainTitle->DrawLatexNDC(0.745+tagShift[tagSelection], 0.85, tagName[tagSelection]);
        } else {
          mainTitle->DrawLatexNDC(0.78, 0.85, "CMS");
        }

        mainTitle->SetTextFont(42);
        mainTitle->SetTextSize(0.04);
        mainTitle->DrawLatexNDC(0.568, 0.95, "1.70 nb^{-1} PbPb (5.02 TeV)");

        mainTitle->SetTextFont(62);
        mainTitle->SetTextSize(supplementaryLegendTextSize);
        if(iCentrality == 0){
          mainTitle->DrawLatexNDC(0.628, 0.78, Form("%s, n=%d", centralityString.Data(), weightExponent));
        } else {
          mainTitle->DrawLatexNDC(0.604, 0.78, Form("%s, n=%d", centralityString.Data(), weightExponent));
        }
        mainTitle->DrawLatexNDC(0.736, 0.71, trackPtString.Data());

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

          eecPpBin = std::make_tuple(weightExponent-1, iJetPt, iTrackPt);
          eecPpUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iJetPt, iTrackPt);
          eecPpCorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kCorrelatedUncertainty, iJetPt, iTrackPt);

          // Set the axis titles and labels
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetTitleOffset(0.85);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetTitleSize(0.08);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetLabelOffset(0.001);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetLabelSize(0.07);

          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->SetTitleOffset(0.9);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->SetTitleSize(0.08);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->SetLabelOffset(0.01);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->SetLabelSize(0.07);

          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->SetTitle("EEC");
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetStats(0);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->SetNdivisions(505);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->CenterTitle();
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetTitle("#Deltar");
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetNdivisions(505);

          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetTitle("");

          // Do the same for PbPb histograms
          for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
            eecPbPbUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);

            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitleOffset(0.85);
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitleSize(0.08);
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetLabelOffset(0.001);
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetLabelSize(0.07);

            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitleOffset(0.9);
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitleSize(0.08);
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetLabelOffset(0.01);
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetLabelSize(0.07);

            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitle("EEC");
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetStats(0);
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetNdivisions(505);
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->CenterTitle();
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitle("#Deltar");
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetNdivisions(505);

            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetTitle("");
          }

          // Set the drawing ranges
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->SetRangeUser(supplementaryDistributionZoomLow[weightExponent-1], supplementaryDistributionZoomHigh[weightExponent-1]);
          systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(0.006, 0.4);

          // Set the drawing style for pp
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetMarkerStyle(markerStylePp);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetMarkerSize(1.2);
          energyEnergyCorrelatorSignalPp[eecPpBin]->SetMarkerStyle(markerStylePp);
          energyEnergyCorrelatorSignalPp[eecPpBin]->SetMarkerSize(1.2);

          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetFillColorAlpha(markerColorPp, 0.4);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetLineColor(markerColorPp);
          energyEnergyCorrelatorSignalPp[eecPpBin]->SetLineColor(markerColorPp);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetMarkerColor(markerColorPp);
          energyEnergyCorrelatorSignalPp[eecPpBin]->SetMarkerColor(markerColorPp);

          systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->SetFillColorAlpha(markerColorPp, 0.4);
          systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->SetMarkerSize(0);
          systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->SetMarkerStyle(9);
      

          // Set the style also for PbPb plots
          for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
            eecPbPbBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt);
            eecPbPbUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
            eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeUp, iCentrality, iJetPt, iTrackPt);
            eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeDown, iCentrality, iJetPt, iTrackPt);

            // Set drawing ranges
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetRangeUser(supplementaryDistributionZoomLow[weightExponent-1], supplementaryDistributionZoomHigh[weightExponent-1]);
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->GetXaxis()->SetRangeUser(0.006, 0.4);
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->GetXaxis()->SetRangeUser(0.006, 0.4);

            // Set the drawing style for the PbPb histograms
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetMarkerSize(1.2);
            energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
            energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerSize(1.2);

            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
            energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
            energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);

            // Different drawing styles for correlated uncertainties
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->SetFillColorAlpha(bandColorUpPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->SetFillColorAlpha(bandColorDownPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerSize(0);
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerStyle(9);
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerSize(0);
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerStyle(9);
          }

          // First, draw all uncorrelated systematic uncertainties
          if(iDrawingStyle == 6){
            systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->Draw("e2");
            for(int iCentrality = lastDrawnCentralityBin[weightExponent-1]; iCentrality >= firstDrawnCentralityBin[weightExponent-1]; iCentrality--){
              eecPbPbUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
              systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->Draw("same,e2");
            }
          } else if (iDrawingStyle == 4){
            systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->Draw("e2");
          } else if (iDrawingStyle == 5){
            eecPbPbUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, 0, iJetPt, iTrackPt);
            systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->Draw("e2");
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->Draw("same,e2");
          } else {
            eecPbPbUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iDrawingStyle, iJetPt, iTrackPt);
            systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->Draw("e2");
          }

          // Then, draw all correlated uncertainties
          if(iDrawingStyle == 6){
            systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->Draw("same,e3");
            for(int iCentrality = lastDrawnCentralityBin[weightExponent-1]; iCentrality >= firstDrawnCentralityBin[weightExponent-1]; iCentrality--){
              eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeUp, iCentrality, iJetPt, iTrackPt);
              eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeDown, iCentrality, iJetPt, iTrackPt);
              systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->Draw("same,e3");
              systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->Draw("same,e3");
            }
          } else if(iDrawingStyle == 4){
            systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->Draw("same,e3");
          } else if(iDrawingStyle == 5){
            eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeUp, 0, iJetPt, iTrackPt);
            eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeDown, 0, iJetPt, iTrackPt);
            systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->Draw("same,e3");
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->Draw("same,e3");
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->Draw("same,e3");
          } else {
            eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeUp, iDrawingStyle, iJetPt, iTrackPt);
            eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeDown, iDrawingStyle, iJetPt, iTrackPt);
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->Draw("same,e3");
            systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->Draw("same,e3");
          }

          // Finally, draw the data points themselves
          if(iDrawingStyle == 6){
            for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
              eecPbPbBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt);
              energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->Draw("same,p");
            }
            energyEnergyCorrelatorSignalPp[eecPpBin]->Draw("same,p");
          } else if(iDrawingStyle == 4){
            energyEnergyCorrelatorSignalPp[eecPpBin]->Draw("same,p");
          } else if(iDrawingStyle == 5){
            eecPbPbBin = std::make_tuple(weightExponent-1, 0, iJetPt, iTrackPt);
            energyEnergyCorrelatorSignalPp[eecPpBin]->Draw("same,p");
            energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->Draw("same,p");
          } else {
            eecPbPbBin = std::make_tuple(weightExponent-1, iDrawingStyle, iJetPt, iTrackPt);
            energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->Draw("same,p");
          }

          // Make a legend for the jet pT bins
          if(physicsBriefingStyle){
            legend = new TLegend(0.18, 0.24, 0.47, 0.4);
          } else {
            legend = new TLegend(0.18, 0.18, 0.47, 0.52);
          }
          
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(supplementaryLegendTextSize); legend->SetTextFont(62);

          // Add legends from each distribution
          for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
            eecPbPbUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
            if(physicsBriefingStyle){
              if(iCentrality == 0){
                legend->AddEntry(systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin], "With quark-gluon plasma", "lpf");
              }
            } else {
              legend->AddEntry(systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin], Form("PbPb %.0f-%.0f%%", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality)), "lpf");
            }            
          }

          if(physicsBriefingStyle){
            legend->AddEntry(systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin], "Without quark-gluon plasma", "lpf");
          } else {
            legend->AddEntry(systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin], "pp", "lpf");
          }

          // Draw the legend
          legend->Draw();

          // Semi-transparent boxes for some drawing styles
          box = new TBox();
          box->SetFillColorAlpha(kWhite, 0.95);

          // No white boxes for physics briefing style
          if(!physicsBriefingStyle){

            // The pp text needs to be hidden for drawing styles 0-3
            if(iDrawingStyle < 4){
              if(weightExponent == 1){
                box->DrawBox(0.008, 0.49, 0.015, 0.7); // Hide the pp text (n=1)
              } else {
                box->DrawBox(0.008, 0.053, 0.015, 0.1); // Hide the pp text (n=2)
              }
            }

            // The 50-90% text needs to be drawn for styles 3 and 6
            if(iDrawingStyle != 3 && iDrawingStyle !=6){
              if(weightExponent == 1){
                box->DrawBox(0.008, 0.7, 0.041, 0.92); // Hide the 50-90% text (n=1)
              } else {
                box->DrawBox(0.008, 0.1, 0.041, 0.2); // Hide the 50-90% text (n=2)
              }
            }
          
            // The 30-50% text needs to be drawn for styles 2 and 6
            if(iDrawingStyle != 2 && iDrawingStyle != 6){
              if(weightExponent == 1){
                box->DrawBox(0.008, 0.92, 0.041, 1.33); // Hide the 30-50% text (n=1)
              } else {
                box->DrawBox(0.008, 0.2, 0.041, 0.35); // Hide the 30-50% text (n=2)
              }
            }
          
            // The 10-30% text needs to be drawn for styles 1 and 6
            if(iDrawingStyle != 1 && iDrawingStyle != 6){
              if(weightExponent == 1){
                box->DrawBox(0.008, 1.33, 0.041, 1.88); // Hide the 10-30% text (n=1)
              } else {
                box->DrawBox(0.008, 0.35, 0.041, 0.7); // Hide the 10-30% text (n=2)
              }
            }

            // The 0-10% text needs to be drawn for styles 0, 5 and 6
            if(iDrawingStyle > 0 && iDrawingStyle < 5){
              if(weightExponent == 1){
                box->DrawBox(0.008, 1.88, 0.04, 2.75); // Hide the 0-10% text (n=1)
              } else {
                box->DrawBox(0.008, 0.7, 0.04, 1.2); // Hide the 0-10% text (n=2)
              }
            }

          } // Physics briefing drawing style

          // Add some interesting information to the plot
          mainTitle->SetTextFont(62);
          mainTitle->SetTextSize(0.06);
          if(physicsBriefingStyle){
            mainTitle->DrawLatexNDC(0.185, 0.85, "CMS");
          } else {
            mainTitle->DrawLatexNDC(0.185, 0.55, "CMS");
          }

          if(tagSelection){
            mainTitle->SetTextFont(52);
            mainTitle->SetTextSize(0.045);
            if(physicsBriefingStyle){
              mainTitle->DrawLatexNDC(0.315, 0.85, tagName[tagSelection]);
            } else {
              mainTitle->DrawLatexNDC(0.315, 0.55, tagName[tagSelection]);
            }
          }

          mainTitle->SetTextFont(42);
          mainTitle->SetTextSize(0.04);
          mainTitle->DrawLatexNDC(0.165, 0.95, "1.70 nb^{-1} PbPb (5.02 TeV) + 302 pb^{-1} pp (5.02 TeV)");

          mainTitle->SetTextFont(62);
          mainTitle->SetTextSize(supplementaryLegendTextSize);

          if(!physicsBriefingStyle){
            mainTitle->DrawLatexNDC(0.522, 0.85, jetPtString.Data());
            mainTitle->DrawLatexNDC(0.675, 0.78, "anti-k_{T} R = 0.4");
            mainTitle->DrawLatexNDC(0.765, 0.71, "|#eta_{jet}| < 1.6");
            mainTitle->DrawLatexNDC(0.726, 0.64, trackPtString.Data());
            mainTitle->DrawLatexNDC(0.87, 0.57, Form("n=%d", weightExponent));
          }

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
            eecPbPbBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt);
            eecPbPbUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
            eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeUp, iCentrality, iJetPt, iTrackPt);
            eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeDown, iCentrality, iJetPt, iTrackPt);

            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitleOffset(0.85);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitleSize(0.08);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetLabelOffset(0.001);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetLabelSize(0.07);

            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitleOffset(0.9);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitleSize(0.08);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetLabelOffset(0.01);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetLabelSize(0.07);

            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitle("#frac{PbPb}{pp}");
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetStats(0);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetNdivisions(505);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->CenterTitle();
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitle("#Deltar");
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetNdivisions(505);

            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetTitle("");

            // Set drawing ranges
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetRangeUser(0.25, 1.75);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->GetXaxis()->SetRangeUser(0.006, 0.4);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->GetXaxis()->SetRangeUser(0.006, 0.4);

            // Set the drawing style
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetMarkerSize(1.2);
            energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
            energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->SetMarkerSize(1.2);

            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
            energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
            energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);

            // Different drawing styles for correlated uncertainties
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->SetFillColorAlpha(bandColorUpPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->SetFillColorAlpha(bandColorDownPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerSize(0);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerStyle(9);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerSize(0);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerStyle(9);
          }

          // First, draw all uncorrelated systematic uncertainties
          if(iDrawingStyle == 4){
            eecPbPbUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, lastDrawnCentralityBin[weightExponent-1], iJetPt, iTrackPt);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->Draw("e2");
            for(int iCentrality = lastDrawnCentralityBin[weightExponent-1]-1; iCentrality >= firstDrawnCentralityBin[weightExponent-1]; iCentrality--){
              eecPbPbUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
              systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->Draw("same,e2");
            }
          } else {
            eecPbPbUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iDrawingStyle, iJetPt, iTrackPt);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->Draw("e2");
          }

          // Then, draw all correlated uncertainties
          if(iDrawingStyle == 4){
            for(int iCentrality = lastDrawnCentralityBin[weightExponent-1]; iCentrality >= firstDrawnCentralityBin[weightExponent-1]; iCentrality--){
              eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeUp, iCentrality, iJetPt, iTrackPt);
              eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeDown, iCentrality, iJetPt, iTrackPt);
              systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->Draw("same,e3");
              systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->Draw("same,e3");
            }
          } else {
            eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeUp, iDrawingStyle, iJetPt, iTrackPt);
            eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeDown, iDrawingStyle, iJetPt, iTrackPt);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->Draw("same,e3");
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->Draw("same,e3");
          }

          // Finally, draw the data points themselves
          if(iDrawingStyle == 4){
            for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
              eecPbPbBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt);
              energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->Draw("same,p");
            }
          }  else {
            eecPbPbBin = std::make_tuple(weightExponent-1, iDrawingStyle, iJetPt, iTrackPt);
            energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->Draw("same,p");
          }

          // Make a legend for the jet pT bins
          legend = new TLegend(0.22, 0.18, 0.54, 0.45);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(supplementaryLegendTextSize); legend->SetTextFont(62);

          // Add legends from each distribution
          for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
            eecPbPbUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
            legend->AddEntry(systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin], Form("PbPb %.0f-%.0f%%", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality)), "lpf");
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
          cmsPosition = (tagSelection == kSupplementaryTag) ? 0.22 : 0.2;
          mainTitle->SetTextFont(62);
          mainTitle->SetTextSize(0.06);

          mainTitle->DrawLatexNDC(cmsPosition, 0.86, "CMS");

          if(tagSelection){
            mainTitle->SetTextFont(52);
            mainTitle->SetTextSize(0.045);
            mainTitle->DrawLatexNDC(0.4 - tagSelection*0.09, 0.92 - tagSelection*0.06, tagName[tagSelection]);
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

          eecPbPbUncorrelatedSystematicsBin = std::make_tuple(0, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);

          // Create a TCanvas for each centrality+jet pT+track pT bin
          theGreatCanvasOfRatioWeightExponentEdition = new SplitCanvas(Form("theGreatCanvasOfWeightExponentRatio%d%d%d", iTrackPt, iJetPt, iCentrality), "", 1300, 1100);
          theGreatCanvasOfRatioWeightExponentEdition->SetMargin(0.18, 0.01, 0.15, 0.07);
          theGreatCanvasOfRatioWeightExponentEdition->cd();

          gPad->SetLogy(false);
          gPad->SetLogx();

          mainTitle = new TLatex();

          // Set the axis titles and labels
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitleOffset(0.85);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitleSize(0.08);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetLabelOffset(0.001);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetLabelSize(0.07);

          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitleOffset(0.9);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitleSize(0.08);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetLabelOffset(0.01);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetLabelSize(0.07);

          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetTitle("#frac{PbPb}{pp}");
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetStats(0);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetNdivisions(505);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTitle("#Deltar");
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetNdivisions(505);

          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetTitle("");

          // Set drawing ranges
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetRangeUser(0.25, 1.75);
          for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
            eecPbPbBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt);
            eecPbPbUncorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
            eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(iWeightExponent, kCorrelatedUncertaintyShapeUp, iCentrality, iJetPt, iTrackPt);
            eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(iWeightExponent, kCorrelatedUncertaintyShapeDown, iCentrality, iJetPt, iTrackPt);

            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->GetXaxis()->SetRangeUser(0.006, 0.4);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->GetXaxis()->SetRangeUser(0.006, 0.4);


            // Set the drawing style
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetMarkerSize(1.2);
            energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->SetMarkerStyle(markerStylePbPb[iWeightExponent]);
            energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->SetMarkerSize(1.2);

            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetFillColorAlpha(markerColorPbPb[iWeightExponent], 0.4);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetLineColor(markerColorPbPb[iWeightExponent]);
            energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->SetLineColor(markerColorPbPb[iWeightExponent]);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetMarkerColor(markerColorPbPb[iWeightExponent]);
            energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->SetMarkerColor(markerColorPbPb[iWeightExponent]);

            // Different drawing styles for correlated uncertainties
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->SetFillColorAlpha(bandColorUpPbPb[iWeightExponent], 0.4);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->SetFillColorAlpha(bandColorDownPbPb[iWeightExponent], 0.4);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerSize(0);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerStyle(9);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerSize(0);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerStyle(9);
          }

          // First, draw all uncorrelated systematic uncertainties
          eecPbPbUncorrelatedSystematicsBin = std::make_tuple(0, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->Draw("e2");
          for(int iWeightExponent = 1; iWeightExponent < nWeightExponents; iWeightExponent++){
            eecPbPbUncorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
            systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->Draw("same,e2");
          }
          
          // Then, draw all correlated uncertainties
          for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
            eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(iWeightExponent, kCorrelatedUncertaintyShapeUp, iCentrality, iJetPt, iTrackPt);
            eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(iWeightExponent, kCorrelatedUncertaintyShapeDown, iCentrality, iJetPt, iTrackPt);
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->Draw("same,e3");
            systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->Draw("same,e3");
          }

          // Finally, draw the data points themselves
          for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
            eecPbPbBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt);
            energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->Draw("same,p");
          }

          // Make a legend for the jet pT bins
          legend = new TLegend(0.22, 0.2, 0.48, 0.34);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(supplementaryLegendTextSize); legend->SetTextFont(62);

          // Add legends from each distribution
          for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
            eecPbPbUncorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
            legend->AddEntry(systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin], Form("n = %d", iWeightExponent+1), "lpf");
          }

          // Draw the legend
          legend->Draw();

          // Draw a lone to one
          oneLine->Draw();

          // Add some interesting information to the plot
          cmsPosition = (tagSelection == kSupplementaryTag) ? 0.22 : 0.2;
          if(iTrackPt > 1) cmsPosition = 0.2;
          mainTitle->SetTextFont(62);
          mainTitle->SetTextSize(0.06);
          mainTitle->DrawLatexNDC(cmsPosition, 0.86, "CMS");

          if(tagSelection){
            mainTitle->SetTextFont(52);
            mainTitle->SetTextSize(0.045);
            if(iTrackPt > 1){
              mainTitle->DrawLatexNDC(0.31, 0.86, tagName[tagSelection]);
            } else {
              mainTitle->DrawLatexNDC(0.4 - tagSelection*0.09, 0.915 - tagSelection*0.055, tagName[tagSelection]);
            }
          }

          mainTitle->SetTextFont(42);
          mainTitle->SetTextSize(0.04);
          mainTitle->DrawLatexNDC(0.29, 0.95, "1.70 nb^{-1} PbPb (5.02 TeV) + 302 pb^{-1} pp (5.02 TeV)");

          mainTitle->SetTextFont(62);
          mainTitle->SetTextSize(supplementaryLegendTextSize);

          cmsPosition = 0;
          if(iTrackPt > 1) cmsPosition = tagShift[tagSelection];

          mainTitle->DrawLatexNDC(0.53 - cmsPosition, 0.86, jetPtString.Data());
          mainTitle->DrawLatexNDC(0.657 - cmsPosition, 0.79, "anti-k_{T} R = 0.4");
          mainTitle->DrawLatexNDC(0.735 - cmsPosition, 0.72, "|#eta_{jet}| < 1.6");
          mainTitle->DrawLatexNDC(0.701 - cmsPosition, 0.65, trackPtString.Data());

          mainTitle->DrawLatexNDC(0.22, 0.38, Form("#frac{%s}{pp}", centralityString.Data()));

          // Add some interesting information to the plot
          /*mainTitle->SetTextFont(62);
          mainTitle->SetTextSize(0.08);

          if(tagSelection){
            mainTitle->DrawLatexNDC(0.065, 0.93, "CMS");

            mainTitle->SetTextFont(52);
            mainTitle->SetTextSize(0.05);
            mainTitle->DrawLatexNDC(0.035, 0.88, tagName[tagSelection]);
          } else {
            mainTitle->DrawLatexNDC(0.02, 0.89, "CMS");
          }

          mainTitle->SetTextFont(42);
          mainTitle->SetTextSize(0.055);

          if(tagSelection){
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

  // Draw supplementary plots where all centrality bins are drawn into a single canvas
  if(drawSupplementaryEECShiftIllustration){

    TString shiftedJetPtString;
    TString extraString;

    TCanvas* theGreatCanvasOfShiftiness;

    for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){

      centralityString = Form("PbPb %.0f-%.0f%%", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality));
      compactCentralityString = Form("_C=%.0f-%.0f", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality));

      // Track pT loop
      for(int iTrackPt = firstDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt <= lastDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt++){

        trackPtString = Form("p_{T}^{ch} > %.0f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt));
        compactTrackPtString = Form("_T>%.1f", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt));
        compactTrackPtString.ReplaceAll(".","v");

        // Jet pT loop
        for(int iJetPt = firstDrawnJetPtBinEEC[weightExponent-1]; iJetPt < lastDrawnJetPtBinEEC[weightExponent-1]; iJetPt++){

          eecPpBin = std::make_tuple(weightExponent-1, iJetPt, iTrackPt);
          eecPpReferenceBin = std::make_tuple(weightExponent-1, iJetPt+1, iTrackPt);
          eecPpUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iJetPt, iTrackPt);
          eecPpCorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kCorrelatedUncertainty, iJetPt, iTrackPt);
          eecPpReferenceUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iJetPt+1, iTrackPt);
          eecPpReferenceCorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kCorrelatedUncertainty, iJetPt+1, iTrackPt);

          eecPbPbBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt);
          eecPbPbUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeUp, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeDown, iCentrality, iJetPt, iTrackPt);

          jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt));
          compactJetPtString = Form("_J=%.0f-%.0f", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt));
          shiftedJetPtString = Form("%.0f < jet p_{T} < %.0f GeV", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt+1), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt+1));

          // Create a TCanvas for each centrality+track/jet pT bin
          theGreatCanvasOfShiftiness = new SplitCanvas(Form("theGreatCanvasOfShiftiness%d%d%d", iTrackPt, iJetPt, iCentrality), "", 1100, 1100);
          theGreatCanvasOfShiftiness->SetMargin(0.15, 0.01, 0.15, 0.07);
          theGreatCanvasOfShiftiness->cd();

          gPad->SetLogy();
          gPad->SetLogx();

          mainTitle = new TLatex();

          // Set the axis titles and labels
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetTitleOffset(0.85);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetTitleSize(0.08);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetLabelOffset(0.001);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetLabelSize(0.07);

          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->SetTitleOffset(0.9);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->SetTitleSize(0.08);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->SetLabelOffset(0.01);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->SetLabelSize(0.07);

          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->SetTitle("EEC");
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetStats(0);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->SetNdivisions(505);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->CenterTitle();
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetTitle("#Deltar");
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->CenterTitle();
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetNdivisions(505);

          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetTitle("");

          // Set the drawing ranges
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->SetRangeUser(supplementaryDistributionZoomLow[weightExponent-1], supplementaryDistributionZoomHigh[weightExponent-1]);
          systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(0.006, 0.4);

          // Set the drawing style for pp in two jet pT bins
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetMarkerStyle(kFullCross);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetMarkerSize(1.2);
          energyEnergyCorrelatorSignalPp[eecPpBin]->SetMarkerStyle(kFullCross);
          energyEnergyCorrelatorSignalPp[eecPpBin]->SetMarkerSize(1.2);

          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetFillColorAlpha(kBlack, 0.4);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetLineColor(kBlack);
          energyEnergyCorrelatorSignalPp[eecPpBin]->SetLineColor(kBlack);
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetMarkerColor(kBlack);
          energyEnergyCorrelatorSignalPp[eecPpBin]->SetMarkerColor(kBlack);

          systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->SetFillColorAlpha(kBlack, 0.4);
          systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->SetMarkerSize(0);
          systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->SetMarkerStyle(9);

          systematicUncertaintyForPp[eecPpReferenceUncorrelatedSystematicsBin]->SetMarkerStyle(kFullCrossX);
          systematicUncertaintyForPp[eecPpReferenceUncorrelatedSystematicsBin]->SetMarkerSize(1.2);
          energyEnergyCorrelatorSignalPp[eecPpReferenceBin]->SetMarkerStyle(kFullCrossX);
          energyEnergyCorrelatorSignalPp[eecPpReferenceBin]->SetMarkerSize(1.2);

          systematicUncertaintyForPp[eecPpReferenceUncorrelatedSystematicsBin]->SetFillColorAlpha(kViolet+2, 0.4);
          systematicUncertaintyForPp[eecPpReferenceUncorrelatedSystematicsBin]->SetLineColor(kViolet+2);
          energyEnergyCorrelatorSignalPp[eecPpReferenceBin]->SetLineColor(kViolet+2);
          systematicUncertaintyForPp[eecPpReferenceUncorrelatedSystematicsBin]->SetMarkerColor(kViolet+2);
          energyEnergyCorrelatorSignalPp[eecPpReferenceBin]->SetMarkerColor(kViolet+2);

          systematicUncertaintyForPp[eecPpReferenceCorrelatedSystematicsBin]->SetFillColorAlpha(kViolet+2, 0.4);
          systematicUncertaintyForPp[eecPpReferenceCorrelatedSystematicsBin]->SetMarkerSize(0);
          systematicUncertaintyForPp[eecPpReferenceCorrelatedSystematicsBin]->SetMarkerStyle(9);
      

          // Set the style also for PbPb distributions
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetMarkerSize(1.2);
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerSize(1.2);

          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);

          // Different drawing styles for correlated uncertainties
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->SetFillColorAlpha(bandColorUpPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->SetFillColorAlpha(bandColorDownPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerSize(0);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerStyle(9);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerSize(0);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerStyle(9);

          // First, draw all uncorrelated systematic uncertainties
          systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->Draw("e2");
          if(!hideShiftInIllustration) systematicUncertaintyForPp[eecPpReferenceUncorrelatedSystematicsBin]->Draw("same,e2");
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->Draw("same,e2");

          // Then, draw all correlated uncertainties
          systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->Draw("same,e3");
          if(!hideShiftInIllustration) systematicUncertaintyForPp[eecPpReferenceCorrelatedSystematicsBin]->Draw("same,e3");
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->Draw("same,e3");
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->Draw("same,e3");

          // Finally, draw the data points themselves
          energyEnergyCorrelatorSignalPp[eecPpBin]->Draw("same,p");
          if(!hideShiftInIllustration) energyEnergyCorrelatorSignalPp[eecPpReferenceBin]->Draw("same,p");
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->Draw("same,p");

          // Make a legend for the jet pT bins
          legend = new TLegend(0.18, 0.18, 0.47, 0.41);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.04); legend->SetTextFont(62);

          // Add legends from each distribution
          legend->AddEntry(systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin], Form("pp, %s", jetPtString.Data()), "lpf");
          legend->AddEntry(systematicUncertaintyForPp[eecPpReferenceUncorrelatedSystematicsBin], Form("pp, %s", shiftedJetPtString.Data()), "lpf");
          legend->AddEntry(systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin], Form("%s, %s", centralityString.Data(), jetPtString.Data()), "lpf");

          // Draw the legend
          legend->Draw();

          // Semi-transparent boxes for some drawing styles
          box = new TBox();
          box->SetFillColorAlpha(kWhite, 0.95);

          extraString = "";
          if(hideShiftInIllustration){
            box->DrawBox(0.008, 0.7, 0.1, 1); // Hide shifted pT text
            extraString = "HiddenShift";
          }

          // Add some interesting information to the plot
          mainTitle->SetTextFont(62);
          mainTitle->SetTextSize(0.06);
          mainTitle->DrawLatexNDC(0.185, 0.85, "CMS");

          if(tagSelection){
            mainTitle->SetTextFont(52);
            mainTitle->SetTextSize(0.045);
            mainTitle->DrawLatexNDC(0.315, 0.85, tagName[tagSelection]);
          }

          mainTitle->SetTextFont(42);
          mainTitle->SetTextSize(0.04);
          mainTitle->DrawLatexNDC(0.165, 0.95, "1.70 nb^{-1} PbPb (5.02 TeV) + 302 pb^{-1} pp (5.02 TeV)");

          mainTitle->SetTextFont(62);
          mainTitle->SetTextSize(supplementaryLegendTextSize);

          //mainTitle->DrawLatexNDC(0.522, 0.85, jetPtString.Data());
          mainTitle->DrawLatexNDC(0.675, 0.85, "anti-k_{T} R = 0.4");
          mainTitle->DrawLatexNDC(0.765, 0.78, "|#eta_{jet}| < 1.6");
          mainTitle->DrawLatexNDC(0.726, 0.71, trackPtString.Data());
          mainTitle->DrawLatexNDC(0.87, 0.64, Form("n=%d", weightExponent));

          // Save the figures to file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/finalEECShiftIllustration%s%s%s%s%s.%s", extraString.Data(), saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
          }
        } // Drawing style loop
      } // Track pT loop
    } // Centrality loop
  } // If for drawing all centrality bins in one plot
}
