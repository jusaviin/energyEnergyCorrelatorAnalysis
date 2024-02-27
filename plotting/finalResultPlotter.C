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
  const int weightExponent = 1;

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
  inputFileName[kPbPb][0] = "data/eecAnalysis_akFlowJet_nominalEnergyWeight_optimizedUnfoldingBins_fixedCovarianceMatrix_updatedBackgroundSubtraction_processed_2024-02-23.root";
  inputFileName[kPbPb][1] = "data/eecAnalysis_akFlowJet_energyWeightSquared_optimizedUnfoldingBins_fixedCovarianceMatrix_updatedBackgroundSubtraction_processed_2024-02-23.root";
  inputFileName[kPp][0] = "data/ppData_pfJets_wtaAxis_nominalEnergyWeight_optimizedUnfoldingBins_fixedCovarianceMatrix_jet60or80triggers_unfoldingWithCovariance_processed_2024-01-23.root";
  inputFileName[kPp][1] = "data/ppData_pfJets_wtaAxis_energyWeightSquared_optimizedUnfoldingBins_fixedCovarianceMatrix_jet60or80triggers_unfoldingWithCovariance_processed_2024-01-23.root";
  TString uncertaintyFileName[kNDataTypes][nWeightExponents];
  uncertaintyFileName[kPbPb][0] = "systematicUncertainties/systematicUncertainties_PbPb_nominalEnergyWeight_includeMCnonClosure_2024-02-23.root";
  uncertaintyFileName[kPbPb][1] = "systematicUncertainties/systematicUncertainties_PbPb_energyWeightSquared_includeMCnonClosure_2024-02-23.root";
  uncertaintyFileName[kPp][0] = "systematicUncertainties/systematicUncertainties_pp_nominalEnergyWeight_includeMCnonClosure_2024-01-22.root";
  uncertaintyFileName[kPp][1] = "systematicUncertainties/systematicUncertainties_pp_energyWeightSquared_includeMCnonClosure_2024-01-22.root";

  TString shiftedPtFileName[nWeightExponents];
  shiftedPtFileName[0] = "data/ppData_pfJets_wtaAxis_shiftedJetPtBins_nominalEnergyWeight_jet60or80trigger_unfoldWithShiftedResponse_processed_2024-02-19.root";
  shiftedPtFileName[1] = "data/ppData_pfJets_wtaAxis_shiftedJetPtBins_energyWeightSquared_jet60or80trigger_unfoldWithShiftedResponse_processed_2024-02-19.root";

  // eecAnalysis_akFlowJet_nominalEnergyWeight_optimizedUnfoldingBins_fixedCovarianceMatrix_unfoldingWithCovariance_processed_2024-01-23.root
  // eecAnalysis_akFlowJet_nominalEnergyWeight_optimizedUnfoldingBins_unfoldingWithNominalSmear_processed_2024-01-17.root
  // eecAnalysis_akFlowJet_energyWeightSquared_optimizedUnfoldingBins_fixedCovarianceMatrix_unfoldingWithCovariance_processed_2024-01-23.root
  // eecAnalysis_akFlowJet_energyWeightSquared_optimizedUnfoldingBins_unfoldingWithNominalSmear_processed_2024-01-17.root
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
  bool drawDoubleRatioToSingleCanvas = true;
  bool drawBigCanvasAllRatios = false; // Draw ratios with all defined energy weight exponents to the same figure

  bool drawShiftedPtRatio = false;

  // Normalize all distributions to 2 GeV integral
  bool normalizeTo2GeV = false;
  int trackPtBinFor2GeV[nWeightExponents];
  std::pair<int, int> trackPtBinsForDoubleRatio[nWeightExponents];
  trackPtBinsForDoubleRatio[0] = std::make_pair(card[kPbPb][0]->GetBinIndexTrackPtEEC(2.0), card[kPbPb][0]->GetBinIndexTrackPtEEC(3.0));
  trackPtBinsForDoubleRatio[1] = std::make_pair(card[kPbPb][1]->GetBinIndexTrackPtEEC(1.0), card[kPbPb][1]->GetBinIndexTrackPtEEC(2.0));
  double trackPtForAllRatiosComparison = 2;

  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    trackPtBinFor2GeV[iWeightExponent] = card[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(2.0);

    if(drawBigCanvasAllRatios){
      firstDrawnTrackPtBinEEC[iWeightExponent] = card[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtForAllRatiosComparison);
      lastDrawnTrackPtBinEEC[iWeightExponent] = card[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtForAllRatiosComparison);
    }
  }

  // Select the bins to be drawn for double ratio plots
  std::pair<double, double> doubleRatioCentralityBin1 = std::make_pair(0.0,10.0);
  std::pair<double, double> doubleRatioCentralityBin2 = std::make_pair(10.0,30.0);
  std::pair<double, double> doubleRatioJetPtBin = std::make_pair(180,200);
  int doubleRatioCentralityBinIndex1;
  int doubleRatioCentralityBinIndex2;
  int doubleRatioJetPtBinIndex;
  TString doubleRatioJetPtString = Form("_%.0f<jetpt<%.0f", doubleRatioJetPtBin.first, doubleRatioJetPtBin.second);
  
  // Save the final plots
  const bool saveFigures = true;
  TString energyWeightString[nWeightExponents] = {"_nominalEnergyWeight", "_energyWeightSquared"};
  TString saveComment = energyWeightString[weightExponent-1] + "_optimizedUnfoldingBins_updatedBackground";
  TString figureFormat = "pdf";

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
    shiftedPtHistograms[iWeightExponent] = new EECHistogramManager(shiftedPtFile[iWeightExponent], shiftedPtCard[iWeightExponent]);
  
    // Load all unfolded energy-energy correlators
    shiftedPtHistograms[iWeightExponent]->SetLoadEnergyEnergyCorrelators(true);
    shiftedPtHistograms[iWeightExponent]->SetCentralityBinRange(0, shiftedPtCard[iWeightExponent]->GetNCentralityBins());
    shiftedPtHistograms[iWeightExponent]->SetTrackPtBinRangeEEC(0, shiftedPtCard[iWeightExponent]->GetNTrackPtBinsEEC());
    shiftedPtHistograms[iWeightExponent]->SetJetPtBinRangeEEC(0, shiftedPtCard[iWeightExponent]->GetNJetPtBinsEEC());

    // Load the histograms from the file
    shiftedPtHistograms[iWeightExponent]->LoadProcessedHistograms();

  } // Weight exponent loop
 
  // Energy-energy correlators and PbPb to pp ratios
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
        energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt] = NULL;
        energyEnergyCorrelatorSignalShiftedPt[iWeightExponent][iJetPt][iTrackPt] = NULL;
        energyEnergyCorrelatorShiftedPtRatio[iWeightExponent][iJetPt][iTrackPt] = NULL;

        energyEnergyCorrelatorForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt] = NULL;
        systematicUncertaintyForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt] = NULL;

        for(int iUncertainty = 0; iUncertainty < knSystematicUncertaintyTypes; iUncertainty++){
          systematicUncertaintyForPp[iWeightExponent][iUncertainty][iJetPt][iTrackPt] = NULL;
        } // Uncertainty type loop

        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
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

  // Define helper histograms to determine the uncertainties relevant for the double ratio
  TH1D* uncertaintyTrackSelection;
  TH1D* uncertaintyBackgroundSubtraction;
  TH1D* uncertaintyTrackPairEfficiency;
  TH1D* uncertaintyMCnonclosure;
  double trackCorrelation = 0.2;  // TODO: Check what is the actual number of overlapping tracks
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
      iJetPtMatchedShifted = shiftedPtCard[iWeightExponent]->FindBinIndexJetPtEEC(shiftedJetPtBin);
      for(int iTrackPt = firstDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt <= lastDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt++){
        iTrackPtMatchedPp = card[kPp][iWeightExponent]->FindBinIndexTrackPtEEC(card[kPbPb][iWeightExponent]->GetBinBordersTrackPtEEC(iTrackPt));
        iTrackPtMatchedPbPbUncertainty = uncertaintyCard[kPbPb][iWeightExponent]->FindBinIndexTrackPtEEC(card[kPbPb][iWeightExponent]->GetBinBordersTrackPtEEC(iTrackPt));
        iTrackPtMatchedPpUncertainty = uncertaintyCard[kPp][iWeightExponent]->FindBinIndexTrackPtEEC(card[kPbPb][iWeightExponent]->GetBinBordersTrackPtEEC(iTrackPt));
        iTrackPtMatchedShifted = shiftedPtCard[iWeightExponent]->FindBinIndexTrackPtEEC(card[kPbPb][iWeightExponent]->GetBinBordersTrackPtEEC(iTrackPt));

        // Read the pp histograms that do not have centrality binning
        energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt] = histograms[kPp][iWeightExponent]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatchedPp, iTrackPtMatchedPp, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);
        systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt] = uncertainties[kPp][iWeightExponent]->GetUncorrelatedSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty);
        systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt] = uncertainties[kPp][iWeightExponent]->GetCorrelatedSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty);

        // Read the relevant systematic uncertainties for double ratio from pp
        uncertaintyBackgroundSubtraction = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kBackgroundSubtraction);
        uncertaintyTrackPairEfficiency = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kTrackPairEfficiency);
        uncertaintyMCnonclosure = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kMonteCarloNonClosure);

        systematicUncertaintyForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt] = (TH1D*) uncertaintyBackgroundSubtraction->Clone(Form("doubleRatioUncertaintyFromPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));
        for(int iBin = 1; iBin <= uncertaintyBackgroundSubtraction->GetNbinsX(); iBin++){
          systematicUncertaintyForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt]->SetBinError(iBin, TMath::Sqrt(TMath::Power(uncertaintyBackgroundSubtraction->GetBinError(iBin),2) + TMath::Power(uncertaintyTrackPairEfficiency->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintyMCnonclosure->GetBinError(iBin),2)));
        }

        // Read the histograms with shifted pT spectrum
        energyEnergyCorrelatorSignalShiftedPt[iWeightExponent][iJetPt][iTrackPt] = shiftedPtHistograms[iWeightExponent]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatchedShifted, iTrackPtMatchedShifted, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);


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

          systematicUncertaintyForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt] = (TH1D*) uncertaintyBackgroundSubtraction->Clone(Form("doubleRatioUncertaintyFromPbPb%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          for(int iBin = 1; iBin <= uncertaintyBackgroundSubtraction->GetNbinsX(); iBin++){
            systematicUncertaintyForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetBinError(iBin, TMath::Sqrt(TMath::Power(uncertaintyBackgroundSubtraction->GetBinError(iBin),2) + TMath::Power(uncertaintyTrackPairEfficiency->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintyTrackSelection->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintyMCnonclosure->GetBinError(iBin),2)));
          }
        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // Weight exponent loop

  // ================================================ //
  //   Normalize the histograms and calculate ratios  //
  // ================================================ //

  double epsilon = 0.0001;
  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    int referenceTrackPtBin = trackPtBinFor2GeV[iWeightExponent];
    int lowAnalysisBin = energyEnergyCorrelatorSignalPbPb[iWeightExponent][firstDrawnCentralityBin[iWeightExponent]][firstDrawnJetPtBinEEC[iWeightExponent]][firstDrawnTrackPtBinEEC[iWeightExponent]]->GetXaxis()->FindBin(analysisDeltaR.first + epsilon);
    int highAnalysisBin = energyEnergyCorrelatorSignalPbPb[iWeightExponent][firstDrawnCentralityBin[iWeightExponent]][firstDrawnJetPtBinEEC[iWeightExponent]][firstDrawnTrackPtBinEEC[iWeightExponent]]->GetXaxis()->FindBin(analysisDeltaR.second - epsilon);
    for(int iJetPt = firstDrawnJetPtBinEEC[iWeightExponent]; iJetPt <= lastDrawnJetPtBinEEC[iWeightExponent]; iJetPt++){
      for(int iTrackPt = lastDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt >= firstDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt--){
        if(!normalizeTo2GeV) referenceTrackPtBin = iTrackPt;
        energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->Scale(1.0 / energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][referenceTrackPtBin]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
        systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt]->Scale(energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->GetBinContent(10) / systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetBinContent(10));
        systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt]->Scale(energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->GetBinContent(10) / systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt]->GetBinContent(10));
        systematicUncertaintyForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt]->Scale(energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->GetBinContent(10) / systematicUncertaintyForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt]->GetBinContent(10));

        // Shifted pT to regular pT ratios
        energyEnergyCorrelatorSignalShiftedPt[iWeightExponent][iJetPt][iTrackPt]->Scale(1.0 / energyEnergyCorrelatorSignalShiftedPt[iWeightExponent][iJetPt][referenceTrackPtBin]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

        energyEnergyCorrelatorShiftedPtRatio[iWeightExponent][iJetPt][iTrackPt] = (TH1D*) energyEnergyCorrelatorSignalShiftedPt[iWeightExponent][iJetPt][iTrackPt]->Clone(Form("shiftedPtRatio%d%d%d", iWeightExponent, iJetPt, iTrackPt));
        energyEnergyCorrelatorShiftedPtRatio[iWeightExponent][iJetPt][iTrackPt]->Divide(energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]);

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
          energyEnergyCorrelatorDoubleRatio[iWeightExponent][iCentrality][iJetPt] = (TH1D*) energyEnergyCorrelatorForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][trackPtBinsForDoubleRatio[iWeightExponent].second]->Clone(Form("energyEnergyCorrelatorDoubleRatio%d%d%d", iWeightExponent, iCentrality, iJetPt));
          energyEnergyCorrelatorDoubleRatio[iWeightExponent][iCentrality][iJetPt]->Divide(energyEnergyCorrelatorForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][trackPtBinsForDoubleRatio[iWeightExponent].first]);
          systematicUncertaintyDoubleRatio[iWeightExponent][iCentrality][iJetPt] = (TH1D*) systematicUncertaintyForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][trackPtBinsForDoubleRatio[iWeightExponent].second]->Clone(Form("systematicUncertaintyDoubleRatio%d%d%d", iWeightExponent, iCentrality, iJetPt));
          systematicUncertaintyDoubleRatio[iWeightExponent][iCentrality][iJetPt]->Divide(systematicUncertaintyForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][trackPtBinsForDoubleRatio[iWeightExponent].first]);
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
          lineDrawer->DrawLine(0.08, 0.8-(weightExponent-1)*0.2, 0.08, 1.2);
          lineDrawer->DrawLine(0.2, 0.8-(weightExponent-1)*0.2, 0.2, 1.2);

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

      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetTitle(Form("#frac{PbPb/pp (p_{T}^{ch} > %.0f GeV)}{PbPb/pp (p_{T}^{ch} > %.0f GeV)}", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(trackPtBinsForDoubleRatio[weightExponent-1].second), card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(trackPtBinsForDoubleRatio[weightExponent-1].first)));
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

    // For the double ratios, only draw one example bin
    doubleRatioJetPtBinIndex = card[kPbPb][weightExponent-1]->FindBinIndexJetPtEEC(doubleRatioJetPtBin);

    // Create a TCanvas for the double ratio
    TCanvas* theGreatCanvasOfDoubleRatio = new SplitCanvas("theGreatCanvasOfDoubleRatio", "", 1200, 800);
    theGreatCanvasOfDoubleRatio->SetMargin(0.24, 0.01, 0.2, 0.15);
    theGreatCanvasOfDoubleRatio->cd();

    gPad->SetLogy(false);
    gPad->SetLogx();

    double doubleRatioZoomMagnitude = 0.2;
    if(trackPtBinsForDoubleRatio[weightExponent-1].first < 2) doubleRatioZoomMagnitude = 0.3;
    if(weightExponent == 1) doubleRatioZoomMagnitude = 0.3;

    mainTitle = new TLatex();
    canvasIndex = 0;

    // Set the drawing and label styles
    for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){

      // Set the axis titles and labels
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetTitleOffset(1);
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetTitleSize(0.09);
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetLabelOffset(0.01);
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetXaxis()->SetLabelSize(0.07);

      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetTitleOffset(1.2);
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetTitleSize(0.08);
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetLabelOffset(0.01);
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetLabelSize(0.07);

      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetTitle(Form("#frac{PbPb/pp (p_{T}^{ch} > %.0f GeV)}{PbPb/pp (p_{T}^{ch} > %.0f GeV)}", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(trackPtBinsForDoubleRatio[weightExponent-1].second), card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(trackPtBinsForDoubleRatio[weightExponent-1].first)));
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

      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->GetYaxis()->SetRangeUser(1 - doubleRatioZoomMagnitude, 1 + doubleRatioZoomMagnitude);

      // Set the drawing styles 
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->SetMarkerSize(1.2);
      energyEnergyCorrelatorDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
      energyEnergyCorrelatorDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->SetMarkerSize(1.2);

      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]], 0.4);
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
      energyEnergyCorrelatorDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
      energyEnergyCorrelatorDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin[weightExponent-1]]);
    }

    // Draw the histograms
    systematicUncertaintyDoubleRatio[weightExponent-1][firstDrawnCentralityBin[weightExponent-1]][doubleRatioJetPtBinIndex]->Draw("e2");

    for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]+1; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
      systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->Draw("same,e2");
    }
    for(int iCentrality = lastDrawnCentralityBin[weightExponent-1]; iCentrality >= firstDrawnCentralityBin[weightExponent-1]; iCentrality--){
      energyEnergyCorrelatorDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex]->Draw("same,p");
    }

    // Show the centrality bin in the legend
    legend = new TLegend(0.28, 0.22, 0.53, 0.5);
    legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.055); legend->SetTextFont(62);
    for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
      legend->AddEntry(systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][doubleRatioJetPtBinIndex], Form("Centrality: %.0f-%.0f%%", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality)), "p");
    }

    legend->Draw();

    // Draw the line to one
    oneLine->Draw();

    mainTitle->SetTextFont(62);
    mainTitle->SetTextSize(0.08);
    mainTitle->DrawLatexNDC(0.3, 0.9, "CMS");

    mainTitle->SetTextFont(42);
    mainTitle->SetTextSize(0.055);
    mainTitle->DrawLatexNDC(0.47, 0.945, "PbPb #sqrt{s_{NN}} = 5.02 TeV, 1.70 nb^{-1}");
    mainTitle->DrawLatexNDC(0.47, 0.87, "pp #sqrt{s} = 5.02 TeV, 302 pb^{-1}");
    mainTitle->DrawLatexNDC(0.3, 0.68, "anti-k_{T} R = 0.4,");
    mainTitle->DrawLatexNDC(0.54, 0.68, Form("n=%d", weightExponent));
    mainTitle->DrawLatexNDC(0.3, 0.77, Form("%.f < p_{T,jet} < %.0f GeV, |#eta_{jet}| < 1.6",  card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(doubleRatioJetPtBinIndex), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(doubleRatioJetPtBinIndex)));

    // Save the figures to file
    if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/finalDoubleRatioSingleCanvas%s.%s", saveComment.Data(), figureFormat.Data()));
    }

  } // If for drawing double ratios

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
            legend = new TLegend(-0.03 + leftMarginAdder, 0.83, 0.5 / (1 - leftMarginAdder), 0.93);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry((TObject*)0, Form("#downarrow %s#downarrow", jetPtString.Data()), "");
            legend->Draw();

          }

          // Add the centrality bin information to the leftmost column
          if(iJetPt == firstDrawnJetPtBinEEC[0]){
            legend = new TLegend(0.19, 0.05 + bottomMarginAdder*1.1, 0.7, 0.2 / (thisPadScale*0.78) + bottomMarginAdder*0.1);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry((TObject*)0, Form("#frac{PbPb %.0f-%.0f%%}{pp} #rightarrow", card[kPbPb][0]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][0]->GetHighBinBorderCentrality(iCentrality)), "");
            legend->Draw();
          }

          // Add information about the nominal weight histogram to selected oad
          if(iJetPt == firstDrawnJetPtBinEEC[0]+1 && iCentrality == firstDrawnCentralityBin[0]+1){
            legend = new TLegend(0.1, 0.83, 0.5, 0.93);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry(systematicUncertaintyPbPbToPpRatio[0][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], "EEC ratio with n=1", "p,e2");
            legend->Draw();
          }

          // Add information about the squared weight histogram to selected oad
          if(iJetPt == firstDrawnJetPtBinEEC[0]+2 && iCentrality == firstDrawnCentralityBin[0]+1){
            legend = new TLegend(0.1, 0.83, 0.5, 0.93);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.1 * thisPadScale); legend->SetTextFont(62);
            legend->AddEntry(systematicUncertaintyPbPbToPpRatio[1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], "EEC ratio with n=2", "p,e2");
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
          lineDrawer->DrawLine(0.08, 0.6, 0.08, 1.2);
          lineDrawer->DrawLine(0.2, 0.6, 0.2, 1.2);

        }  // Jet pT loop
      }    // Centrality loop

      // Draw all necessary CMS text to the plot
      bigDualRatioCanvas[iTrackPt]->cd(0);

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
}
