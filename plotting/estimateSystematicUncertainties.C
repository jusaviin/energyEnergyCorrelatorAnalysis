#include "SystematicUncertaintyOrganizer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "EECHistogramManager.h"
#include "JDrawer.h"
#include "AlgorithmLibrary.h"

#include <algorithm>

// Function definitions. Implementations after the main macro.
TH1D* findTheDifference(TH1D* nominalResult, TH1D* variedResult[], const int nComparisonGraphs);
TH1D* findTheDifference(TH1D* nominalResult, TH1D* variedResult);
void drawIllustratingPlots(JDrawer* drawer, TH1D* nominalResult, TH1D* variedResult[], const int nVariations, const int iCentrality, const int iJetPt, const int iTrackPt, EECCard* card, TString comparisonLegend[], TString plotName, std::pair<double, double> drawingRange, std::pair<double, double> ratioZoom);
void drawIllustratingPlots(JDrawer* drawer, TH1D* nominalResult, TH1D* variedResult, const int iCentrality, const int iJetPt, const int iTrackPt, EECCard* card, TString comparisonLegend, TString plotName, std::pair<double, double> drawingRange, std::pair<double, double> ratioZoom);
void loadRelevantHistograms(EECHistogramManager* histograms);
void loadTrackingSystematicsHistograms(EECHistogramManager* histograms);

/*
 * Macro for estimating systematic uncertainties for energy-energy correlators
 *
 * The following sources of systematic uncertainty are estimated:
 *  - Jet energy resolition uncertainty (modifies unfolding response matrix, correlated)
 *  - Jet energy scale uncertainty (modifies unfolding response matrix, correlated)
 *  - Truth reference for jet pT unfolding (modifies unfolding response matrix, correlated)
 *  - Tracking efficiency (uncorrelated)
 *  - Track pair efficiency (uncorrelated)
 *  - Background subtraction (uncorrelated)
 *
 *  The results will be saved to a file and also slides with different contributions can be printed
 *
 * Arguments:
 *  const int weightExponent = Exponents for the energy weight in energy-energy correlators. Currently 1 and 2 are implemented.
 */
void estimateSystematicUncertainties(const int weightExponent = 1){

  // First, do a sanity check for the weight exponents. Currently only 1 and 2 are implemented.
  if(weightExponent < 1 || weightExponent > 2){
    cout << "ERROR! The weight exponent you gave has not been implemented!" << endl;
    cout << "Currently only values 1 and 2 are implemented." << endl;
    cout << "Please select one of these values for weight exponent." << endl;
    return;
  }
  
  // ==================================================================
  // ============================= Input ==============================
  // ==================================================================

  // Nominal results
  TString nominalResultFileName[2] = {"data/eecAnalysis_akFlowJet_wtaAxis_newTrackPairEfficiencySmoothed_unfoldingWithNominalSmear_processed_2023-07-13.root", "data/eecAnalysis_akFlowJet_wtaAxis_energyWeightSquared_firstFinalResultsWithFixedCard_processed_2023-10-23.root"};
  TFile* nominalResultFile = TFile::Open(nominalResultFileName[weightExponent-1]);
  EECCard* nominalResultCard = new EECCard(nominalResultFile);
  EECHistogramManager* nominalHistogramManager = new EECHistogramManager(nominalResultFile, nominalResultCard);
  loadRelevantHistograms(nominalHistogramManager);
  
  // Results unfolded with a response matrix smeared with jet energy resolution
  TString jetEnergyResolutionSmearDownFileName[2] = {"data/eecAnalysis_akFlowJet_wtaAxis_newTrackPairEfficiencySmoothed_unfoldingWithUncertaintySmearDown_processed_2023-07-13.root", "data/eecAnalysis_akFlowJet_wtaAxis_energyWeightSquared_unfoldingWithUncertaintySmearDown_processed_2023-10-23.root"};
  TString jetEnergyResolutionSmearUpFileName[2] = {"data/eecAnalysis_akFlowJet_wtaAxis_newTrackPairEfficiencySmoothed_unfoldingWithUncertaintySmearUp_processed_2023-07-13.root", "data/eecAnalysis_akFlowJet_wtaAxis_energyWeightSquared_unfoldingWithUncertaintySmearUp_processed_2023-10-23.root"};
  TFile* jetEnergyResolutionFile[2];
  jetEnergyResolutionFile[0] = TFile::Open(jetEnergyResolutionSmearDownFileName[weightExponent-1]);
  jetEnergyResolutionFile[1] = TFile::Open(jetEnergyResolutionSmearUpFileName[weightExponent-1]);
  EECCard* jetEnergyResolutionCard[2];
  EECHistogramManager* jetEnergyResolutionHistogramManager[2];
  for(int iJetEnergyResolutionFile = 0; iJetEnergyResolutionFile < 2; iJetEnergyResolutionFile++){
    jetEnergyResolutionCard[iJetEnergyResolutionFile] = new EECCard(jetEnergyResolutionFile[iJetEnergyResolutionFile]);
    jetEnergyResolutionHistogramManager[iJetEnergyResolutionFile] = new EECHistogramManager(jetEnergyResolutionFile[iJetEnergyResolutionFile], jetEnergyResolutionCard[iJetEnergyResolutionFile]);
    loadRelevantHistograms(jetEnergyResolutionHistogramManager[iJetEnergyResolutionFile]);
  }
  
  // Results unfolded with a response matrix smeared with jet energy scale
  TString jetEnergyScaleMinusFileName[2] = {"data/eecAnalysis_akFlowJet_wtaAxis_newTrackPairEfficiencySmoothed_unfoldingWithMinusJetEnergyScale_processed_2023-07-13.root", "data/eecAnalysis_akFlowJet_wtaAxis_energyWeightSquared_unfoldingWithMinusJetEnergyScale_processed_2023-10-23.root"};
  TString jetEnergyScalePlusFileName[2] = {"data/eecAnalysis_akFlowJet_wtaAxis_newTrackPairEfficiencySmoothed_unfoldingWithPlusJetEnergyScale_processed_2023-07-13.root", "data/eecAnalysis_akFlowJet_wtaAxis_energyWeightSquared_unfoldingWithPlusJetEnergyScale_processed_2023-10-23.root"};
  TFile* jetEnergyScaleFile[2];
  jetEnergyScaleFile[0] = TFile::Open(jetEnergyScaleMinusFileName[weightExponent-1]);
  jetEnergyScaleFile[1] = TFile::Open(jetEnergyScalePlusFileName[weightExponent-1]);
  EECCard* jetEnergyScaleCard[2];
  EECHistogramManager* jetEnergyScaleHistogramManager[2];
  for(int iJetEnergyScaleFile = 0; iJetEnergyScaleFile < 2; iJetEnergyScaleFile++){
    jetEnergyScaleCard[iJetEnergyScaleFile] = new EECCard(jetEnergyScaleFile[iJetEnergyScaleFile]);
    jetEnergyScaleHistogramManager[iJetEnergyScaleFile] = new EECHistogramManager(jetEnergyScaleFile[iJetEnergyScaleFile], jetEnergyScaleCard[iJetEnergyScaleFile]);
    loadRelevantHistograms(jetEnergyScaleHistogramManager[iJetEnergyScaleFile]);
  }

  // Results unfolded with a response matrix where jet pT spectrum is weighted to match the data
  TString jetPtPriorFileName[2] = {"data/eecAnalysis_akFlowJet_wtaAxis_newTrackPairEfficiencySmoothed_unfoldingWithModifiedPrior_processed_2023-07-13.root", "data/eecAnalysis_akFlowJet_wtaAxis_energyWeightSquared_unfoldingWithModifiedPrior_processed_2023-10-23.root"};
  TFile* jetPtPriorFile = TFile::Open(jetPtPriorFileName[weightExponent-1]);
  EECCard* jetPtPriorCard = new EECCard(jetPtPriorFile);
  EECHistogramManager* jetPtPriorHistogramManager = new EECHistogramManager(jetPtPriorFile, jetPtPriorCard);
  loadRelevantHistograms(jetPtPriorHistogramManager);

  // Results where background scaling factor is determined from 2% or 6% shifted simulation instead of nominal 4%
  TString backgroundSubtraction2pFileName[2] = {"data/eecAnalysis_akFlowJet_wtaAxis_newTrackPairEfficiencySmoothed_unfoldingWithNominalSmear_backgroundScaleUncertainty2pShift_processed_2023-07-13.root", "data/eecAnalysis_akFlowJet_wtaAxis_energyWeightSquared_unfoldingWithNominalSmear_backgroundScaleUncertainty2pShift_processed_2023-10-23.root"};
  TString backgroundSubtraction6pFileName[2] = {"data/eecAnalysis_akFlowJet_wtaAxis_newTrackPairEfficiencySmoothed_unfoldingWithNominalSmear_backgroundScaleUncertainty6pShift_processed_2023-07-13.root", "data/eecAnalysis_akFlowJet_wtaAxis_energyWeightSquared_unfoldingWithNominalSmear_backgroundScaleUncertainty6pShift_processed_2023-10-23.root"};
  TFile* backgroundSubtractionFile[2];
  backgroundSubtractionFile[0] = TFile::Open(backgroundSubtraction2pFileName[weightExponent-1]);
  backgroundSubtractionFile[1] = TFile::Open(backgroundSubtraction6pFileName[weightExponent-1]);
  EECCard* backgroundSubtractionCard[2];
  EECHistogramManager* backgroundSubtractionHistogramManager[2];
  for(int iBackgroundSubtractionFile = 0; iBackgroundSubtractionFile < 2; iBackgroundSubtractionFile++){
    backgroundSubtractionCard[iBackgroundSubtractionFile] = new EECCard(backgroundSubtractionFile[iBackgroundSubtractionFile]);
    backgroundSubtractionHistogramManager[iBackgroundSubtractionFile] = new EECHistogramManager(backgroundSubtractionFile[iBackgroundSubtractionFile], backgroundSubtractionCard[iBackgroundSubtractionFile]);
    loadRelevantHistograms(backgroundSubtractionHistogramManager[iBackgroundSubtractionFile]);
  }

  // Result with different track selections
  TString looseTrackSelectionFileName[2] = {"data/eecAnalysis_akFlowJet_wtaAxis_looseTrackSelection_newTrackPairEfficiencySmoothed_unfoldingWithNominalSmear_processed_2023-07-13.root", "data/eecAnalysis_akFlowJet_wtaAxis_energyWeightSquared_looseTrackCuts_processed_2023-10-23.root"};
  TString tightTrackSelectionFileName[2] = {"data/eecAnalysis_akFlowJet_wtaAxis_tightTrackSelection_newTrackPairEfficiencySmoothed_unfoldingWithNominalSmear_processed_2023-07-13.root", "data/eecAnalysis_akFlowJet_wtaAxis_energyWeightSquared_tightTrackCuts_processed_2023-10-23.root"};
  TFile* trackSelectionFile[2];
  trackSelectionFile[0] = TFile::Open(looseTrackSelectionFileName[weightExponent-1]);
  trackSelectionFile[1] = TFile::Open(tightTrackSelectionFileName[weightExponent-1]);
  EECCard* trackSelectionCard[2];
  EECHistogramManager* trackSelectionHistogramManager[2];
  for(int iTrackSelectionFile = 0; iTrackSelectionFile < 2; iTrackSelectionFile++){
    trackSelectionCard[iTrackSelectionFile] = new EECCard(trackSelectionFile[iTrackSelectionFile]);
    trackSelectionHistogramManager[iTrackSelectionFile] = new EECHistogramManager(trackSelectionFile[iTrackSelectionFile], trackSelectionCard[iTrackSelectionFile]);
    loadRelevantHistograms(trackSelectionHistogramManager[iTrackSelectionFile]);
  }

  // Results with varied single and pair track efficiency
  TString trackEfficiencyFileName[2] = {"data/eecAnalysis_akFlowJet_wtaAxis_newTrackPairEfficiencySmoothed_unfoldingWithNominalSmear_trackSystematics_processed_2023-07-13.root", "data/eecAnalysis_akFlowJet_wtaAxis_energyWeightSquared_trackingSystematics_processed_2023-10-23.root"};
  TFile* trackEfficiencyFile = TFile::Open(trackEfficiencyFileName[weightExponent-1]);
  EECCard* trackEfficiencyCard = new EECCard(trackEfficiencyFile);
  EECHistogramManager* trackEfficiencyHistogramManager = new EECHistogramManager(trackEfficiencyFile, trackEfficiencyCard);
  loadTrackingSystematicsHistograms(trackEfficiencyHistogramManager);

  // Results where the jet pT response matrix is determined from 2% or 6% shifted simulation instead of nominal 4%
  TString centralityShift2pFileName[2] = {"data/eecAnalysis_akFlowJet_wtaAxis_newTrackPairEfficiencySmoothed_unfoldingWith2pCentShift_processed_2023-07-13.root", "data/eecAnalysis_akFlowJet_wtaAxis_energyWeightSquared_unfoldingWith2pCentShift_processed_2023-10-23.root"};
  TString centralityShift6pFileName[2] = {"data/eecAnalysis_akFlowJet_wtaAxis_newTrackPairEfficiencySmoothed_unfoldingWith6pCentShift_processed_2023-07-13.root", "data/eecAnalysis_akFlowJet_wtaAxis_energyWeightSquared_unfoldingWith6pCentShift_processed_2023-10-23.root"};
  TFile* centralityShiftFile[2];
  centralityShiftFile[0] = TFile::Open(centralityShift2pFileName[weightExponent-1]);
  centralityShiftFile[1] = TFile::Open(centralityShift6pFileName[weightExponent-1]);
  EECCard* centralityShiftCard[2];
  EECHistogramManager* centralityShiftHistogramManager[2];
  for(int iCentralityShiftFile = 0; iCentralityShiftFile < 2; iCentralityShiftFile++){
    centralityShiftCard[iCentralityShiftFile] = new EECCard(centralityShiftFile[iCentralityShiftFile]);
    centralityShiftHistogramManager[iCentralityShiftFile] = new EECHistogramManager(centralityShiftFile[iCentralityShiftFile], centralityShiftCard[iCentralityShiftFile]);
    loadRelevantHistograms(centralityShiftHistogramManager[iCentralityShiftFile]);
  }

  // File containing relative uncertainties resulting from Monte Carlo non-closure
  TString monteCarloNonClosureFileName[2] = {"systematicUncertainties/monteCarloNonClosureRelative_PbPb_2023-07-16.root", "systematicUncertainties/monteCarloNonClosureRelative_PbPb_energyWeightSquared_2023-11-20.root"};
  TFile* monteCarloNonClosureFile = TFile::Open(monteCarloNonClosureFileName[weightExponent-1]);
  EECCard* monteCarloNonClosureCard = new EECCard(monteCarloNonClosureFile);
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = nominalResultCard->GetNCentralityBins();
  const int nJetPtBinsEEC = nominalResultCard->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = nominalResultCard->GetNTrackPtBinsEEC();
  
  // Estimate the uncertainties for all bins that have been unfolded
  int firstStudiedCentralityBin = nominalResultCard->GetFirstUnfoldedCentralityBin();
  int lastStudiedCentralityBin = nominalResultCard->GetLastUnfoldedCentralityBin();
  
  int firstStudiedJetPtBinEEC = nominalResultCard->GetFirstUnfoldedJetPtBin();
  int lastStudiedJetPtBinEEC = nominalResultCard->GetLastUnfoldedJetPtBin();
  
  int firstStudiedTrackPtBinEEC = nominalResultCard->GetFirstUnfoldedTrackPtBin();
  int lastStudiedTrackPtBinEEC = nominalResultCard->GetLastUnfoldedTrackPtBin();

  // Only draw example plots from selected subset of bins
  vector<int> drawnCentralityBins = {0,1,2,3};
  vector<int> drawnJetPtBins = {6,7,8,9};
  vector<int> drawnTrackPtBins = {3,4,5};
  
  const bool printUncertainties = false;
  
  std::pair<double, double> analysisDeltaR = std::make_pair(0.006, 0.39); // DeltaR span in which the analysis is done
  std::pair<double, double> ratioZoom = std::make_pair(0.8, 1.2);         // Y-axis zoom for ratios
  bool setAutomaticRatioZoom = true;                                      // If true, use predefined ratio zooms for systematic uncertainties
  
  TString nameAdder[] = {"","_energyWeightSquared"}; 
  TString outputFileName = Form("systematicUncertainties/systematicUncertainties_PbPb%s_includeMCnonClosure_2023-11-20.root", nameAdder[weightExponent-1].Data());
  
  // Option to skip evaluating some of the sources defined in SystematicUncertaintyOrganizer or not plotting examples of some
  bool skipUncertaintySource[SystematicUncertaintyOrganizer::knUncertaintySources];
  bool plotExample[SystematicUncertaintyOrganizer::knUncertaintySources];
  for(int iUncertainty = 0; iUncertainty < SystematicUncertaintyOrganizer::knUncertaintySources; iUncertainty++){
    skipUncertaintySource[iUncertainty] = false;
    plotExample[iUncertainty] = false;
  }
  //plotExample[SystematicUncertaintyOrganizer::kCentralityShift] = true;
  //skipUncertaintySource[SystematicUncertaintyOrganizer::kMonteCarloNonClosure] = true;
  
  // ==================================================================
  // ====================== Configuration done ========================
  // ==================================================================
  
  // Systematic uncertainty organizer can easily provide names for all your naming purposes
  SystematicUncertaintyOrganizer* nameGiver = new SystematicUncertaintyOrganizer();

  // Transformer used for Monte Carlo non-closure uncertainty
  AlgorithmLibrary* optimusPrimeTheTransformer = new AlgorithmLibrary();
  
  // Histograms for uncertainty estimation
  TH1D* nominalEnergyEnergyCorrelators[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* jetEnergyResolutionUncertaintyCorrelators[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][2];
  TH1D* jetEnergyScaleUncertaintyCorrelators[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][2];
  TH1D* jetPtPriorUncertaintyCorrelators[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* backgroundSubtractionUncertaintyCorrelators[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][2];
  TH1D* trackSelectionUncertaintyCorrelators[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][2];
  TH1D* singleTrackEfficiencyUncertaintyCorrelators[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][2];
  TH1D* trackPairEfficiencyUncertaintyCorrelators[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][2];
  TH1D* centralityShiftUncertaintyCorrelators[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][2];

  // Histograms to hold the systematic uncertainty results
  TH1D* energyEnergyCorrelatorSystematicUncertainties[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][SystematicUncertaintyOrganizer::knUncertaintySources];
  
  // Initialize all the histograms to NULL
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt] = NULL;
        jetPtPriorUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt] = NULL;
        for(int iFile = 0; iFile < 2; iFile++){
          jetEnergyResolutionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile] = NULL;
          jetEnergyScaleUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile] = NULL;
          backgroundSubtractionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile] = NULL;
          trackSelectionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile] = NULL;
          singleTrackEfficiencyUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile] = NULL;
          trackPairEfficiencyUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile] = NULL;
          centralityShiftUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile] = NULL;
        }
        for(int iUncertainty = 0; iUncertainty < SystematicUncertaintyOrganizer::knUncertaintySources; iUncertainty++){
          energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][iUncertainty] = NULL;
        }
      }
    }
  }

  // Read the histograms from the defined range from the file
  int iCentralityMatched;
  int iTrackPtMatched;
  int iJetPtMatched;
  double epsilon = 0.0001;
  int lowAnalysisBin, highAnalysisBin;
  std::pair<double,double> shiftedCentralityBin;
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

        // Read the nominal histograms
        nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt] = nominalHistogramManager->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

        // Normalize nominal histograms to one in the drawing range
        lowAnalysisBin = nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(analysisDeltaR.first + epsilon);
        highAnalysisBin = nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(analysisDeltaR.second - epsilon);
        nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt]->Scale(1.0 / nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

        // Read the jet pT prior histograms
        iCentralityMatched = jetPtPriorCard->FindBinIndexCentrality(nominalResultCard->GetBinBordersCentrality(iCentrality));
        iTrackPtMatched = jetPtPriorCard->FindBinIndexTrackPtEEC(nominalResultCard->GetBinBordersTrackPtEEC(iTrackPt));
        iJetPtMatched = jetPtPriorCard->FindBinIndexJetPtEEC(nominalResultCard->GetBinBordersJetPtEEC(iJetPt));
        jetPtPriorUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt] = jetPtPriorHistogramManager->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentralityMatched, iJetPtMatched, iTrackPtMatched, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

        // Normalize the jet pT prior histograms to one
        jetPtPriorUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt]->Scale(1.0 / jetPtPriorUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

        // Read the relative Monte Carlo non-closure histograms
        shiftedCentralityBin = nominalResultCard->GetBinBordersCentrality(iCentrality);
        iCentralityMatched = monteCarloNonClosureCard->FindBinIndexCentrality(shiftedCentralityBin.first+4, shiftedCentralityBin.second+4);
        iTrackPtMatched = monteCarloNonClosureCard->FindBinIndexTrackPtEEC(nominalResultCard->GetBinBordersTrackPtEEC(iTrackPt));
        iJetPtMatched = monteCarloNonClosureCard->FindBinIndexJetPtEEC(nominalResultCard->GetBinBordersJetPtEEC(iJetPt));

        energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kMonteCarloNonClosure] = (TH1D*) monteCarloNonClosureFile->Get(Form("relativeUncertaintyMonteCarloClosure_C%dJ%dT%d", iCentralityMatched, iJetPtMatched, iTrackPtMatched));

        // Setup the histogram manager for tracking efficiencies
        iCentralityMatched = trackEfficiencyCard->FindBinIndexCentrality(nominalResultCard->GetBinBordersCentrality(iCentrality));
        iTrackPtMatched = trackEfficiencyCard->FindBinIndexTrackPtEEC(nominalResultCard->GetBinBordersTrackPtEEC(iTrackPt));
        iJetPtMatched = trackEfficiencyCard->FindBinIndexJetPtEEC(nominalResultCard->GetBinBordersJetPtEEC(iJetPt));

        // Read the histograms for single track efficiency
        singleTrackEfficiencyUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][0] = trackEfficiencyHistogramManager->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationPlus, iCentralityMatched, iJetPtMatched, iTrackPtMatched, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);
        singleTrackEfficiencyUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][1] = trackEfficiencyHistogramManager->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationMinus, iCentralityMatched, iJetPtMatched, iTrackPtMatched, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

        // Read the histograms for track pair efficiency
        trackPairEfficiencyUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][0] = trackEfficiencyHistogramManager->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationPlus, iCentralityMatched, iJetPtMatched, iTrackPtMatched, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);
        trackPairEfficiencyUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][1] = trackEfficiencyHistogramManager->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationMinus, iCentralityMatched, iJetPtMatched, iTrackPtMatched, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

        for(int iFile = 0; iFile < 2; iFile++){

          // Normalize the single track efficiency histograms to one
          singleTrackEfficiencyUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile]->Scale(1.0 / singleTrackEfficiencyUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

          // Normalize the track pair efficiency histograms to one
          trackPairEfficiencyUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile]->Scale(1.0 / trackPairEfficiencyUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

          // Read the jet energy resolution histograms
          iCentralityMatched = jetEnergyResolutionCard[iFile]->FindBinIndexCentrality(nominalResultCard->GetBinBordersCentrality(iCentrality));
          iTrackPtMatched = jetEnergyResolutionCard[iFile]->FindBinIndexTrackPtEEC(nominalResultCard->GetBinBordersTrackPtEEC(iTrackPt));
          iJetPtMatched = jetEnergyResolutionCard[iFile]->FindBinIndexJetPtEEC(nominalResultCard->GetBinBordersJetPtEEC(iJetPt));
          jetEnergyResolutionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile] = jetEnergyResolutionHistogramManager[iFile]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentralityMatched, iJetPtMatched, iTrackPtMatched, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

          // Normalize the jet energy resolution histograms to one
          jetEnergyResolutionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile]->Scale(1.0 / jetEnergyResolutionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

          // Read the jet energy scale histograms
          iCentralityMatched = jetEnergyScaleCard[iFile]->FindBinIndexCentrality(nominalResultCard->GetBinBordersCentrality(iCentrality));
          iTrackPtMatched = jetEnergyScaleCard[iFile]->FindBinIndexTrackPtEEC(nominalResultCard->GetBinBordersTrackPtEEC(iTrackPt));
          iJetPtMatched = jetEnergyScaleCard[iFile]->FindBinIndexJetPtEEC(nominalResultCard->GetBinBordersJetPtEEC(iJetPt));
          jetEnergyScaleUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile] = jetEnergyScaleHistogramManager[iFile]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentralityMatched, iJetPtMatched, iTrackPtMatched, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

          // Normalize the jet energy scale histograms to one
          jetEnergyScaleUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile]->Scale(1.0 / jetEnergyScaleUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

          // Read the background subtraction uncertainty histograms
          iCentralityMatched = backgroundSubtractionCard[iFile]->FindBinIndexCentrality(nominalResultCard->GetBinBordersCentrality(iCentrality));
          iTrackPtMatched = backgroundSubtractionCard[iFile]->FindBinIndexTrackPtEEC(nominalResultCard->GetBinBordersTrackPtEEC(iTrackPt));
          iJetPtMatched = backgroundSubtractionCard[iFile]->FindBinIndexJetPtEEC(nominalResultCard->GetBinBordersJetPtEEC(iJetPt));
          backgroundSubtractionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile] = backgroundSubtractionHistogramManager[iFile]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentralityMatched, iJetPtMatched, iTrackPtMatched, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

          // Normalize the background subtraction histograms to one
          backgroundSubtractionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile]->Scale(1.0 / backgroundSubtractionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

          // Read the track selection uncertainty histograms
          iCentralityMatched = trackSelectionCard[iFile]->FindBinIndexCentrality(nominalResultCard->GetBinBordersCentrality(iCentrality));
          iTrackPtMatched = trackSelectionCard[iFile]->FindBinIndexTrackPtEEC(nominalResultCard->GetBinBordersTrackPtEEC(iTrackPt));
          iJetPtMatched = trackSelectionCard[iFile]->FindBinIndexJetPtEEC(nominalResultCard->GetBinBordersJetPtEEC(iJetPt));
          trackSelectionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile] = trackSelectionHistogramManager[iFile]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentralityMatched, iJetPtMatched, iTrackPtMatched, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

          // Normalize the background subtraction histograms to one
          trackSelectionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile]->Scale(1.0 / trackSelectionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

          // Read the centrality shift uncertainty histograms
          iCentralityMatched = centralityShiftCard[iFile]->FindBinIndexCentrality(nominalResultCard->GetBinBordersCentrality(iCentrality));
          iTrackPtMatched = centralityShiftCard[iFile]->FindBinIndexTrackPtEEC(nominalResultCard->GetBinBordersTrackPtEEC(iTrackPt));
          iJetPtMatched = centralityShiftCard[iFile]->FindBinIndexJetPtEEC(nominalResultCard->GetBinBordersJetPtEEC(iJetPt));
          centralityShiftUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile] = centralityShiftHistogramManager[iFile]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentralityMatched, iJetPtMatched, iTrackPtMatched, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

          // Normalize the background subtraction histograms to one
          centralityShiftUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile]->Scale(1.0 / centralityShiftUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

        }

      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
  
  // Drawer for illustrative plots
  JDrawer* drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  drawer->SetLogX(true);

  // Legends given to drawns graphs
  TString legendNames[2];
  
  // ================================================= //
  //   Uncertainty coming from jet energy resolution   //
  // ================================================= //
  
  if(!skipUncertaintySource[SystematicUncertaintyOrganizer::kJetEnergyResolution]){
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

          energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kJetEnergyResolution] = findTheDifference(nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt], jetEnergyResolutionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt], 2);

          // Draw example plots on how the uncertainty is obtained
          if(plotExample[SystematicUncertaintyOrganizer::kJetEnergyResolution] && std::binary_search(drawnCentralityBins.begin(), drawnCentralityBins.end(), iCentrality) && std::binary_search(drawnJetPtBins.begin(), drawnJetPtBins.end(), iJetPt) && std::binary_search(drawnTrackPtBins.begin(), drawnTrackPtBins.end(), iTrackPt)){
            legendNames[0] = "JER varied response down";
            legendNames[1] = "JER varied response up";

            // Set reasonable ratio zoom depending on centrality bin
            if(setAutomaticRatioZoom){
              ratioZoom = std::make_pair(0.9,1.1);
            }

            drawIllustratingPlots(drawer, nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt], jetEnergyResolutionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt], 2, iCentrality, iJetPt, iTrackPt, nominalResultCard, legendNames, nameGiver->GetSystematicUncertaintyName(SystematicUncertaintyOrganizer::kJetEnergyResolution), analysisDeltaR, ratioZoom);
          }
      
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Jet energy resolution uncertainty
  
  // ============================================ //
  //   Uncertainty coming from jet energy scale   //
  // ============================================ //

  if(!skipUncertaintySource[SystematicUncertaintyOrganizer::kJetEnergyScale]){
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

          energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kJetEnergyScale] = findTheDifference(nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt], jetEnergyScaleUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt], 2);

          // Draw example plots on how the uncertainty is obtained
          if(plotExample[SystematicUncertaintyOrganizer::kJetEnergyScale] && std::binary_search(drawnCentralityBins.begin(), drawnCentralityBins.end(), iCentrality) && std::binary_search(drawnJetPtBins.begin(), drawnJetPtBins.end(), iJetPt) && std::binary_search(drawnTrackPtBins.begin(), drawnTrackPtBins.end(), iTrackPt)){
            legendNames[0] = "JEC varied response down";
            legendNames[1] = "JEC varied response up";

            // Set reasonable ratio zoom
            if(setAutomaticRatioZoom){
              if(iCentrality == 0){
                ratioZoom = std::make_pair(0.7,1.3);
              } else if (iCentrality == 1){
                ratioZoom = std::make_pair(0.8,1.2);
              } else {
                ratioZoom = std::make_pair(0.85,1.15);
              }
              
            }

            drawIllustratingPlots(drawer, nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt], jetEnergyScaleUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt], 2, iCentrality, iJetPt, iTrackPt, nominalResultCard, legendNames, nameGiver->GetSystematicUncertaintyName(SystematicUncertaintyOrganizer::kJetEnergyScale), analysisDeltaR, ratioZoom);
          }
      
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Jet energy scale uncertainty

  // ==================================================== //
  //   Uncertainty coming from jet pT prior in unfolding  //
  // ==================================================== //

  if(!skipUncertaintySource[SystematicUncertaintyOrganizer::kUnfoldingTruth]){
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

          energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kUnfoldingTruth] = findTheDifference(nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt], jetPtPriorUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt]);

          // Draw example plots on how the uncertainty is obtained
          if(plotExample[SystematicUncertaintyOrganizer::kUnfoldingTruth] && std::binary_search(drawnCentralityBins.begin(), drawnCentralityBins.end(), iCentrality) && std::binary_search(drawnJetPtBins.begin(), drawnJetPtBins.end(), iJetPt) && std::binary_search(drawnTrackPtBins.begin(), drawnTrackPtBins.end(), iTrackPt)){
            legendNames[0] = "Weighted prior";

            // Set reasonable ratio zoom
            if(setAutomaticRatioZoom){
              ratioZoom = std::make_pair(0.9,1.1);
            }

            drawIllustratingPlots(drawer, nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt], jetPtPriorUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt], iCentrality, iJetPt, iTrackPt, nominalResultCard, legendNames[0], nameGiver->GetSystematicUncertaintyName(SystematicUncertaintyOrganizer::kUnfoldingTruth), analysisDeltaR, ratioZoom);
          }
      
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Jet energy scale uncertainty
  
  // ================================================== //
  //   Uncertainty coming from background subtraction   //
  // ================================================== //
  
  if(!skipUncertaintySource[SystematicUncertaintyOrganizer::kBackgroundSubtraction]){
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

          energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kBackgroundSubtraction] = findTheDifference(nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt], backgroundSubtractionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt], 2);

          // Draw example plots on how the uncertainty is obtained
          if(plotExample[SystematicUncertaintyOrganizer::kBackgroundSubtraction] && std::binary_search(drawnCentralityBins.begin(), drawnCentralityBins.end(), iCentrality) && std::binary_search(drawnJetPtBins.begin(), drawnJetPtBins.end(), iJetPt) && std::binary_search(drawnTrackPtBins.begin(), drawnTrackPtBins.end(), iTrackPt)){
            legendNames[0] = "Bg scale from 2% shift";
            legendNames[1] = "Bg scale from 6% shift";

            // Set reasonable ratio zoom
            if(setAutomaticRatioZoom){
              if(iCentrality == 0){
                ratioZoom = std::make_pair(0.85,1.15);
              } else {
                ratioZoom = std::make_pair(0.9,1.1);
              }
            }

            drawIllustratingPlots(drawer, nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt], backgroundSubtractionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt], 2, iCentrality, iJetPt, iTrackPt, nominalResultCard, legendNames, nameGiver->GetSystematicUncertaintyName(SystematicUncertaintyOrganizer::kBackgroundSubtraction), analysisDeltaR, ratioZoom);
          }
      
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Background subtraction uncertainty

  // =============================================== //
  //     Uncertainty coming from track selection     //
  // =============================================== //
  
  if(!skipUncertaintySource[SystematicUncertaintyOrganizer::kTrackSelection]){
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

          energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kTrackSelection] = findTheDifference(nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt], trackSelectionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt], 2);

          // Draw example plots on how the uncertainty is obtained
          if(plotExample[SystematicUncertaintyOrganizer::kTrackSelection] && std::binary_search(drawnCentralityBins.begin(), drawnCentralityBins.end(), iCentrality) && std::binary_search(drawnJetPtBins.begin(), drawnJetPtBins.end(), iJetPt) && std::binary_search(drawnTrackPtBins.begin(), drawnTrackPtBins.end(), iTrackPt)){
            legendNames[0] = "Loose track selection";
            legendNames[1] = "Tight track selection";

            // Set reasonable ratio zoom
            if(setAutomaticRatioZoom){
              ratioZoom = std::make_pair(0.8,1.2);
            }

            drawIllustratingPlots(drawer, nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt], trackSelectionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt], 2, iCentrality, iJetPt, iTrackPt, nominalResultCard, legendNames, nameGiver->GetSystematicUncertaintyName(SystematicUncertaintyOrganizer::kTrackSelection), analysisDeltaR, ratioZoom);
          }
      
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Track selection uncertainty

  // =================================================== //
  //   Uncertainty coming from single track efficiency   //
  // =================================================== //
  
  if(!skipUncertaintySource[SystematicUncertaintyOrganizer::kSingleTrackEfficiency]){
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

          energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kSingleTrackEfficiency] = findTheDifference(nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt], singleTrackEfficiencyUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt], 2);

          // Draw example plots on how the uncertainty is obtained
          if(plotExample[SystematicUncertaintyOrganizer::kSingleTrackEfficiency] && std::binary_search(drawnCentralityBins.begin(), drawnCentralityBins.end(), iCentrality) && std::binary_search(drawnJetPtBins.begin(), drawnJetPtBins.end(), iJetPt) && std::binary_search(drawnTrackPtBins.begin(), drawnTrackPtBins.end(), iTrackPt)){
            legendNames[0] = "Efficiency - 5%";
            legendNames[1] = "Efficiency + 5%";

            // Set reasonable ratio zoom
            if(setAutomaticRatioZoom){
              ratioZoom = std::make_pair(0.999999,1.000001);
            }

            drawIllustratingPlots(drawer, nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt], singleTrackEfficiencyUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt], 2, iCentrality, iJetPt, iTrackPt, nominalResultCard, legendNames, nameGiver->GetSystematicUncertaintyName(SystematicUncertaintyOrganizer::kSingleTrackEfficiency), analysisDeltaR, ratioZoom);
          }
      
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Track selection uncertainty

  // =================================================== //
  //    Uncertainty coming from track pair efficiency    //
  // =================================================== //
  
  if(!skipUncertaintySource[SystematicUncertaintyOrganizer::kTrackPairEfficiency]){
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

          energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kTrackPairEfficiency] = findTheDifference(nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt], trackPairEfficiencyUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt], 2);

          // Draw example plots on how the uncertainty is obtained
          if(plotExample[SystematicUncertaintyOrganizer::kTrackPairEfficiency] && std::binary_search(drawnCentralityBins.begin(), drawnCentralityBins.end(), iCentrality) && std::binary_search(drawnJetPtBins.begin(), drawnJetPtBins.end(), iJetPt) && std::binary_search(drawnTrackPtBins.begin(), drawnTrackPtBins.end(), iTrackPt)){
            legendNames[0] = "Pair efficiency - 10%";
            legendNames[1] = "Pair efficiency + 10%";

            // Set reasonable ratio zoom
            if(setAutomaticRatioZoom){
              ratioZoom = std::make_pair(0.8,1.2);
            }

            drawIllustratingPlots(drawer, nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt], trackPairEfficiencyUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt], 2, iCentrality, iJetPt, iTrackPt, nominalResultCard, legendNames, nameGiver->GetSystematicUncertaintyName(SystematicUncertaintyOrganizer::kTrackPairEfficiency), analysisDeltaR, ratioZoom);
          }
      
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Track selection uncertainty

  // ================================================ //
  //     Uncertainty coming from centrality shift     //
  // ================================================ //
  
  if(!skipUncertaintySource[SystematicUncertaintyOrganizer::kCentralityShift]){
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

          energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kCentralityShift] = findTheDifference(nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt], centralityShiftUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt], 2);

          // Draw example plots on how the uncertainty is obtained
          if(plotExample[SystematicUncertaintyOrganizer::kCentralityShift] && std::binary_search(drawnCentralityBins.begin(), drawnCentralityBins.end(), iCentrality) && std::binary_search(drawnJetPtBins.begin(), drawnJetPtBins.end(), iJetPt) && std::binary_search(drawnTrackPtBins.begin(), drawnTrackPtBins.end(), iTrackPt)){
            legendNames[0] = "Unfolding with 2% shift";
            legendNames[1] = "Unfolding with 6% shift";

            // Set reasonable ratio zoom
            if(setAutomaticRatioZoom){
              if(iCentrality == 0){
                ratioZoom = std::make_pair(0.85,1.15);
              } else {
                ratioZoom = std::make_pair(0.9,1.1);
              }
            }

            drawIllustratingPlots(drawer, nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt], centralityShiftUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt], 2, iCentrality, iJetPt, iTrackPt, nominalResultCard, legendNames, nameGiver->GetSystematicUncertaintyName(SystematicUncertaintyOrganizer::kCentralityShift), analysisDeltaR, ratioZoom);
          }
      
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Background subtraction uncertainty

  // ===================================================== //
  //   Uncertainty coming from non-closure in Monte Carlo  //
  // ===================================================== //

  if(!skipUncertaintySource[SystematicUncertaintyOrganizer::kMonteCarloNonClosure]){
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

          optimusPrimeTheTransformer->SuppressSingleBinFluctuations(energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kMonteCarloNonClosure], analysisDeltaR.first, analysisDeltaR.second, 3, 0.5);
          optimusPrimeTheTransformer->TransformToAbsoluteUncertainty(energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kMonteCarloNonClosure], nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt], false);

          // Note: Illustrating plots for this uncertainty can be plotter with fullAnalysisClosure.C macro
      
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Jet energy scale uncertainty
  
  // ============================================= //
  //   Add all uncertainty sources in quadrature   //
  // ============================================= //
  
  double sumOfSquares;
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

        // Clone the uncertainty histogram from the nominal one
        energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kAll] = (TH1D*) nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt]->Clone(Form("summedUncertainty%d%d%d", iCentrality, iJetPt, iTrackPt));
      
        // Calculate a sum of squared from all evaluated uncertainties
        for(int iBin = 1; iBin <= nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt]->GetNbinsX(); iBin++){
          sumOfSquares = 0;
          for(int iUncertainty = 0; iUncertainty < SystematicUncertaintyOrganizer::kAll; iUncertainty++){
            if(energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][iUncertainty] != NULL){
              sumOfSquares = sumOfSquares + TMath::Power(energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][iUncertainty]->GetBinError(iBin),2);
            }
          } // Uncertainty type loop
          energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kAll]->SetBinError(iBin, TMath::Sqrt(sumOfSquares));
        } // Bin loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
  
  // =========================================================== //
  //   Write the systematic uncertainty histograms into a file   //
  // =========================================================== //
  
  // If an output file name is given, put the total uncertainties to an output file
  if(outputFileName.EndsWith(".root")){
    TFile* outputFile = new TFile(outputFileName,"UPDATE");
    TString saveName;

    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
          for(int iUncertainty = 0; iUncertainty < SystematicUncertaintyOrganizer::knUncertaintySources; iUncertainty++){
            if(energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][iUncertainty] != NULL){
              saveName = Form("systematicUncertainty_%s_C%dJ%dT%d", nameGiver->GetSystematicUncertaintyName(iUncertainty).Data(), iCentrality, iJetPt, iTrackPt);
              energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][iUncertainty]->Write(saveName);
            }
          } // Uncertainty loop
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop

    // Write the card from the nominal file, since all indices are syncronized with that
    nominalResultCard->Write(outputFile);

    outputFile->Close();
  }
  
  // =========================================================================== //
  // Print slides summarizing all the different sources of systemtic uncertainty //
  // =========================================================================== //
  
  // if(printUncertainties){
    
  //   // Print a slide with uncertainties for each source and each centrality for each V
  //   char namer[100];
  //   for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
      
  //     sprintf(namer,"\\frametitle{Uncertainty for jet $v_{%d}$}", iFlow+1);
      
  //     cout << endl;
  //     cout << "\\begin{frame}" << endl;
  //     cout << namer << endl;
  //     cout << "\\begin{center}" << endl;
  //     cout << "  \\begin{tabular}{cccccc}" << endl;
  //     cout << "    \\toprule" << endl;
  //     cout << "    Source & C: 0-10 & C: 10-30 & C: 30-50 \\\\" << endl;
  //     cout << "    \\midrule" << endl;
      
  //     // Set the correct precision for printing floating point numbers
  //     cout << fixed << setprecision(4);
      
  //     // First, print the actual values of vn:s
  //     cout << "$v_{" << iFlow+1 << "}$ ";
      
  //     for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
  //       nominalResultGraph[iFlow]->GetPoint(iCentrality, nominalX, nominalY);
  //       cout << " & $" << nominalY << "$";
  //     }
  //     cout << " \\\\" << endl;
  //     cout << "    \\midrule" << endl;
      
  //     // Print the relative uncertainties
      
  //     for(int iUncertainty = 0; iUncertainty < SystematicUncertaintyOrganizer::knUncertaintySources; iUncertainty++){
  //       cout << uncertaintyNamer->GetLongRangeUncertaintyName(iUncertainty).Data();
  //       for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
  //         cout << " & $" << relativeUncertaintyTable[iUncertainty][iFlow][iCentrality] << "$";
  //       }
  //       cout << " \\\\" << endl;
  //     }
      
  //     cout << "    \\midrule" << endl;
      
  //     // Print the absolute uncertainties
      
  //     for(int iUncertainty = 0; iUncertainty < SystematicUncertaintyOrganizer::knUncertaintySources; iUncertainty++){
  //       cout << uncertaintyNamer->GetLongRangeUncertaintyName(iUncertainty).Data();
  //       for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
  //         cout << " & $" << absoluteUncertaintyTable[iUncertainty][iFlow][iCentrality] << "$";
  //       }
  //       cout << " \\\\" << endl;
  //     }
      
      
  //     cout << "    \\bottomrule" << endl;
  //     cout << "  \\end{tabular}" << endl;
  //     cout << "\\end{center}" << endl;
  //     cout << "\\begin{itemize}" << endl;
  //     cout << "  \\item Top: value. Middle: relative uncertainty. Bottom: absolute uncertainty." << endl;
  //     cout << "\\end{itemize}" << endl;
  //     cout << "\\end{frame}" << endl;
  //     cout << endl;
      
  //   } // vn loop
    
  // } // Printing uncertainties
}

/*
 * Function for finding the relative and absolute difference of points in two graphs
 *
 *  TH1D* nominalResult = Histogram containing nominal results
 *  TH1D* variedResult = Histograms containing varied results for systematic uncertainty estimation
 *  const int nVariations = Number of variations used to estimate the systematic uncertainties
 *
 *  return: Histogram where the error bars in each point correspond to estimated systematic uncertainty
 *
 */
TH1D* findTheDifference(TH1D *nominalResult, TH1D *variedResult[], const int nVariations){
  
  // Create an uncertainty histogram from the nominal results
  TH1D* uncertaintyHistogram = (TH1D*) nominalResult->Clone(Form("uncertaintyFor%s", variedResult[0]->GetName()));
  
  // Loop over the bins in the histogram, and in each bin assign the greates variation from the nominal result as uncertainty
  double biggestDifference, currentDifference;
  for(int iBin = 1; iBin <= nominalResult->GetNbinsX(); iBin++){
    biggestDifference = 0;
    for(int iVariation = 0; iVariation < nVariations; iVariation++){
      currentDifference = TMath::Abs(nominalResult->GetBinContent(iBin) - variedResult[iVariation]->GetBinContent(iBin));
      if(currentDifference > biggestDifference) biggestDifference = currentDifference;
    }
    uncertaintyHistogram->SetBinError(iBin, currentDifference);
  } // Histogram bin loop
  
  // Return the histogram where the errors represent the systematic uncertainties
  return uncertaintyHistogram;
  
}

/*
 * Function for finding the realtive and absolute difference of points in two graphs
 *
 *  TH1D* nominalResult = Histogram containing nominal results
 *  TH1D* variedResult = Histogram containing varied results for systematic uncertainty estimation
 *
 *  return: Histogram where the error bars in each point correspond to estimated systematic uncertainty
 *
 */
TH1D* findTheDifference(TH1D* nominalResult, TH1D* variedResult){
  
  // Enclose the single graph to an array and use the difference finder for graph array
  TH1D* comparisonArray[1] = {variedResult};
  return findTheDifference(nominalResult, comparisonArray, 1);
  
}

/*
 * Draw plots illustrating how systematic uncertainties are estimated form a certain source
 *
 *  JDrawer *drawer = JDrawer doing the dirty work in drawing
 *  TH1D* nominalResult = Histogram containing nominal results
 *  TH1D* variedResult[] = Array of histograms containing variation of results around the nominal one
 *  const int nVariations = Number of variations around the nominal result
 *  const int iCentrality = Centrality index of the considered bin
 *  const int iJetPt = Jet pT index of the considered bin
 *  const int iTrackPt = Track pT index of the considered bin
 *  EECCard* card = Card used to interpret the bin index
 *  TString comparisonLegend[] = An array of strings describing each result variation
 *  TString plotName = String added to saved plots. If left empty, the plots are not saved into files.
 *  std::pair<double, double> drawingRange = Drawing range for x-axis
 *  std::pair<double, double> ratioZoom = Y-axis zoom for the ratio
 */
void drawIllustratingPlots(JDrawer* drawer, TH1D* nominalResult, TH1D* variedResult[], const int nVariations, const int iCentrality, const int iJetPt, const int iTrackPt, EECCard* card, TString comparisonLegend[], TString plotName, std::pair<double, double> drawingRange, std::pair<double, double> ratioZoom){
  
  const int markers[] = {kFullDiamond, kFullDoubleDiamond, kFullCross, kFullFourTrianglesPlus, kFullSquare, kFullStar};
  const int colors[] = {kBlue, kRed, kGreen+3, kMagenta, kCyan, kViolet, kBlack};
  
  // Interpret the given binning information
  TString centralityString = Form("Cent: %.0f-%.0f%%", card->GetLowBinBorderCentrality(iCentrality), card->GetHighBinBorderCentrality(iCentrality));
  TString compactCentralityString = Form("_C=%.0f-%.0f", card->GetLowBinBorderCentrality(iCentrality), card->GetHighBinBorderCentrality(iCentrality));
  TString jetPtString = Form("%.0f < jet p_{T} < %.0f", card->GetLowBinBorderJetPtEEC(iJetPt), card->GetHighBinBorderJetPtEEC(iJetPt));
  TString compactJetPtString = Form("_J=%.0f-%.0f", card->GetLowBinBorderJetPtEEC(iJetPt), card->GetHighBinBorderJetPtEEC(iJetPt));
  TString trackPtString = Form("%.1f < track p_{T}",card->GetLowBinBorderTrackPtEEC(iTrackPt));
  TString compactTrackPtString = Form("_T>%.1f",card->GetLowBinBorderTrackPtEEC(iTrackPt));
  compactTrackPtString.ReplaceAll(".","v");

  // Calculate ratios between nominal and comparison graphs
  TH1D *ratioHistogram[nVariations];
  for(int iVariation = 0; iVariation < nVariations; iVariation++){
    ratioHistogram[iVariation] = (TH1D*) variedResult[iVariation]->Clone(Form("ratio_%s", comparisonLegend[iVariation].Data()));
    ratioHistogram[iVariation]->Divide(nominalResult);
  } // Varied results loop

  // Create a new canvas for the plot
  drawer->CreateSplitCanvas();

  // Use logarithmic axis for EEC
  drawer->SetLogY(true);

  // Setup the legend for plots
  TLegend *legend = new TLegend(0.23,0.3-nVariations*0.06,0.53,0.6);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  legend->AddEntry((TObject*) 0, Form("%s 5.02 TeV",card->GetAlternativeDataType().Data()), "");
  legend->AddEntry((TObject*) 0, centralityString.Data(),"");
  legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
  legend->AddEntry((TObject*) 0, trackPtString.Data(),"");

  // Set the x-axis drawing range
  nominalResult->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);

  // Draw the nominal result to the upper canvas
  nominalResult->SetLineColor(kBlack);
  drawer->DrawHistogramToUpperPad(nominalResult, "#Deltar", "EEC Signal", " ");
  legend->AddEntry(nominalResult, "Nominal", "l");

  // Draw the variations to the same plot
  for(int iVariation = 0; iVariation < nVariations; iVariation++) {
    variedResult[iVariation]->SetLineColor(colors[iVariation]);
    variedResult[iVariation]->Draw("same");
    legend->AddEntry(variedResult[iVariation], comparisonLegend[iVariation], "l");
  }

  // Draw the legends to the upper pad
  legend->Draw();

  // Linear scale for the ratio
  drawer->SetLogY(false);

  // Set the x-axis drawing range
  ratioHistogram[0]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);

  ratioHistogram[0]->SetLineColor(colors[0]);
  ratioHistogram[0]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
  drawer->SetGridY(true);
  drawer->DrawHistogramToLowerPad(ratioHistogram[0], "#Deltar", "#frac{Variation}{Nominal}", " ");
  for(int iVariation = 1; iVariation < nVariations; iVariation++) {
    ratioHistogram[iVariation]->SetLineColor(colors[iVariation]);
    ratioHistogram[iVariation]->Draw("same");
  }
  drawer->SetGridY(false);
  
  // If a plot name is given, save the plot in a file
  if(plotName != ""){
    gPad->GetCanvas()->SaveAs(Form("figures/systematicUncertainty_%s%s%s%s.pdf", plotName.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data()));
  }
}

/*
 * Draw plots illustrating how systematic uncertainties are estimated form a certain source
 *
 *  JDrawer* drawer = JDrawer doing the dirty work in drawing
 *  TH1D* nominalResult = Histogram containing nominal results
 *  TH1D* variedResult = Histogram containing variation of results around the nominal one
 *  const int iCentrality = Centrality index of the considered bin
 *  const int iJetPt = Jet pT index of the considered bin
 *  const int iTrackPt = Track pT index of the considered bin
 *  EECCard* card = Card used to interpret the binning information
 *  TString comparisonLegend = String describing the variation
 *  TString plotName = String added to saved plots. If left empty, the plots are not saved into files.
 *  std::pair<double, double> drawingRange = Drawing range for x-axis
 *  std::pair<double, double> ratioZoom = Y-axis zoom for the ratio
 */
void drawIllustratingPlots(JDrawer* drawer, TH1D* nominalResult, TH1D* variedResult, const int iCentrality, const int iJetPt, const int iTrackPt, EECCard* card, TString comparisonLegend, TString plotName, std::pair<double, double> drawingRange, std::pair<double, double> ratioZoom){
  TH1D* variedResultArray[1] = {variedResult};
  TString comparisonLegendArray[1] = {comparisonLegend};
  drawIllustratingPlots(drawer, nominalResult, variedResultArray, 1, iCentrality, iJetPt, iTrackPt, card, comparisonLegendArray, plotName, drawingRange, ratioZoom);
}

/*
 * Load histograms relevant for systematic uncertainty analysis from the histogram manager 
 */
void loadSelectedHistograms(EECHistogramManager* histograms, bool loadRegular, bool loadTrackSystematics){

  // Load energy-energy correlators
  histograms->SetLoadEnergyEnergyCorrelators(loadRegular);
  histograms->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationPlus(loadTrackSystematics);
  histograms->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationMinus(loadTrackSystematics);
  histograms->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationPlus(loadTrackSystematics);
  histograms->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationMinus(loadTrackSystematics);

  // Load all available bins
  histograms->SetCentralityBinRange(0, histograms->GetCard()->GetNCentralityBins() - 1);
  histograms->SetJetPtBinRangeEEC(0, histograms->GetCard()->GetNJetPtBinsEEC() - 1);
  histograms->SetTrackPtBinRangeEEC(0, histograms->GetCard()->GetNTrackPtBinsEEC() - 1);

  // Load the histograms from the file
  histograms->LoadProcessedHistograms();

}

/*
 * Load histograms relevant for systematic uncertainty analysis from the histogram manager 
 */
void loadRelevantHistograms(EECHistogramManager* histograms){

  // Load only regular energy-energy correlators
  loadSelectedHistograms(histograms, true, false);

}

/*
 * Load histograms related to tracking systematics from the histogram manager 
 */
void loadTrackingSystematicsHistograms(EECHistogramManager* histograms){

  // Load all energy-energy correlators
  loadSelectedHistograms(histograms, true, true);

}

