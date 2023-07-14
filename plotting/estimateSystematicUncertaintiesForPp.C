#include "SystematicUncertaintyOrganizer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "EECHistogramManager.h"
#include "JDrawer.h"

#include <algorithm>

// Function definitions. Implementations after the main macro.
TH1D* findTheDifference(TH1D* nominalResult, TH1D* variedResult[], const int nComparisonGraphs);
TH1D* findTheDifference(TH1D* nominalResult, TH1D* variedResult);
void drawIllustratingPlots(JDrawer* drawer, TH1D* nominalResult, TH1D* variedResult[], const int nVariations, const int iJetPt, const int iTrackPt, EECCard* card, TString comparisonLegend[], TString plotName, std::pair<double, double> drawingRange, std::pair<double, double> ratioZoom);
void drawIllustratingPlots(JDrawer* drawer, TH1D* nominalResult, TH1D* variedResult, const int iJetPt, const int iTrackPt, EECCard* card, TString comparisonLegend, TString plotName, std::pair<double, double> drawingRange, std::pair<double, double> ratioZoom);
void loadRelevantHistograms(EECHistogramManager* histograms);
void loadTrackingSystematicsHistograms(EECHistogramManager* histograms);

/*
 * Macro for estimating systematic uncertainties for energy-energy correlators in pp collisions
 * This is quite similar to the macro for PbPb, but some sources are treated slightly differently
 * and not all of the sources are the same, so I thought it might be simpler to have different 
 * macros for the different systems in this case.
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
 */
void estimateSystematicUncertaintiesForPp(){
  
  // ==================================================================
  // ============================= Input ==============================
  // ==================================================================

  // Nominal results
  TFile* nominalResultFile = TFile::Open("data/ppData_pfJets_wtaAxis_newTrackPairEfficiency_unfoldingWithNominalSmear_processed_2023-07-13.root");
  EECCard* nominalResultCard = new EECCard(nominalResultFile);
  EECHistogramManager* nominalHistogramManager = new EECHistogramManager(nominalResultFile, nominalResultCard);
  loadRelevantHistograms(nominalHistogramManager);
  
  // Results unfolded with a response matrix smeared with jet energy resolution
  TFile* jetEnergyResolutionFile[2];
  jetEnergyResolutionFile[0] = TFile::Open("data/ppData_pfJets_wtaAxis_newTrackPairEfficiency_unfoldingWithUncertaintySmearDown_processed_2023-07-13.root");
  jetEnergyResolutionFile[1] = TFile::Open("data/ppData_pfJets_wtaAxis_newTrackPairEfficiency_unfoldingWithUncertaintySmearUp_processed_2023-07-13.root");
  EECCard* jetEnergyResolutionCard[2];
  EECHistogramManager* jetEnergyResolutionHistogramManager[2];
  for(int iJetEnergyResolutionFile = 0; iJetEnergyResolutionFile < 2; iJetEnergyResolutionFile++){
    jetEnergyResolutionCard[iJetEnergyResolutionFile] = new EECCard(jetEnergyResolutionFile[iJetEnergyResolutionFile]);
    jetEnergyResolutionHistogramManager[iJetEnergyResolutionFile] = new EECHistogramManager(jetEnergyResolutionFile[iJetEnergyResolutionFile], jetEnergyResolutionCard[iJetEnergyResolutionFile]);
    loadRelevantHistograms(jetEnergyResolutionHistogramManager[iJetEnergyResolutionFile]);
  }
  
  // Results unfolded with a response matrix smeared with jet energy scale
  TFile* jetEnergyScaleFile[2];
  jetEnergyScaleFile[0] = TFile::Open("data/ppData_pfJets_wtaAxis_newTrackPairEfficiency_unfoldingWithMinusJetEnergyScale_processed_2023-07-13.root");
  jetEnergyScaleFile[1] = TFile::Open("data/ppData_pfJets_wtaAxis_newTrackPairEfficiency_unfoldingWithPlusJetEnergyScale_processed_2023-07-13.root");
  EECCard* jetEnergyScaleCard[2];
  EECHistogramManager* jetEnergyScaleHistogramManager[2];
  for(int iJetEnergyScaleFile = 0; iJetEnergyScaleFile < 2; iJetEnergyScaleFile++){
    jetEnergyScaleCard[iJetEnergyScaleFile] = new EECCard(jetEnergyScaleFile[iJetEnergyScaleFile]);
    jetEnergyScaleHistogramManager[iJetEnergyScaleFile] = new EECHistogramManager(jetEnergyScaleFile[iJetEnergyScaleFile], jetEnergyScaleCard[iJetEnergyScaleFile]);
    loadRelevantHistograms(jetEnergyScaleHistogramManager[iJetEnergyScaleFile]);
  }

  // Results unfolded with a response matrix where jet pT spectrum is weighted to match the data
  TFile* jetPtPriorFile = TFile::Open("data/ppData_pfJets_wtaAxis_newTrackPairEfficiency_unfoldingWithModifiedPrior_processed_2023-07-13.root");
  EECCard* jetPtPriorCard = new EECCard(jetPtPriorFile);
  EECHistogramManager* jetPtPriorHistogramManager = new EECHistogramManager(jetPtPriorFile, jetPtPriorCard);
  loadRelevantHistograms(jetPtPriorHistogramManager);

  // Results where background scaling factor is determined from peripheral Pythia+Hydjet instead of not subtracting background
  TFile* backgroundSubtractionFile;
  backgroundSubtractionFile = TFile::Open("data/ppData_pfJets_wtaAxis_newTrackPairEfficiency_unfoldingWithNominalSmear_backgroundSubtractionSystematics_processed_2023-07-13.root");
  EECCard* backgroundSubtractionCard = new EECCard(backgroundSubtractionFile);
  EECHistogramManager* backgroundSubtractionHistogramManager = new EECHistogramManager(backgroundSubtractionFile, backgroundSubtractionCard);
  loadRelevantHistograms(backgroundSubtractionHistogramManager);

  // Results with varied single and pair track efficiency
  TFile* trackEfficiencyFile = TFile::Open("data/ppData_pfJets_wtaAxis_trackSystematics_newTrackPairEfficiency_unfoldingWithNominalSmear_processed_2023-07-13.root");
  EECCard* trackEfficiencyCard = new EECCard(trackEfficiencyFile);
  EECHistogramManager* trackEfficiencyHistogramManager = new EECHistogramManager(trackEfficiencyFile, trackEfficiencyCard);
  loadTrackingSystematicsHistograms(trackEfficiencyHistogramManager);
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Find the number of bins from the card
  const int nJetPtBinsEEC = nominalResultCard->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = nominalResultCard->GetNTrackPtBinsEEC();
  
  // Estimate the uncertainties for all bins that have been unfolded
  int firstStudiedJetPtBinEEC = nominalResultCard->GetFirstUnfoldedJetPtBin();
  int lastStudiedJetPtBinEEC = nominalResultCard->GetLastUnfoldedJetPtBin();
  
  int firstStudiedTrackPtBinEEC = nominalResultCard->GetFirstUnfoldedTrackPtBin();
  int lastStudiedTrackPtBinEEC = nominalResultCard->GetLastUnfoldedTrackPtBin();

  // Only draw example plots from selected subset of bins
  vector<int> drawnJetPtBins = {2,3,4,5};
  vector<int> drawnTrackPtBins = {3,4,5};
  
  const bool printUncertainties = false;
  
  std::pair<double, double> analysisDeltaR = std::make_pair(0.006, 0.39); // DeltaR span in which the analysis is done
  std::pair<double, double> ratioZoom = std::make_pair(0.9, 1.1);         // Y-axis zoom for rations
  bool setAutomaticRatioZoom = true;                                      // If true, use predefined ratio zooms for systematic uncertainties
  
  TString outputFileName = "systematicUncertainties/systematicUncertaintiesForPp_jetMetUpdate_2023-07-14.root";
  
  // Option to skip evaluating some of the sources defined in SystematicUncertaintyOrganizer or not plotting examples of some
  bool skipUncertaintySource[SystematicUncertaintyOrganizer::knUncertaintySources];
  bool plotExample[SystematicUncertaintyOrganizer::knUncertaintySources];
  for(int iUncertainty = 0; iUncertainty < SystematicUncertaintyOrganizer::knUncertaintySources; iUncertainty++){
    skipUncertaintySource[iUncertainty] = false;
    plotExample[iUncertainty] = false;
  }
  //plotExample[SystematicUncertaintyOrganizer::kBackgroundSubtraction] = true;
  
  // ==================================================================
  // ====================== Configuration done ========================
  // ==================================================================
  
  // Systematic uncertainty organizer can easily provide names for all your naming purposes
  SystematicUncertaintyOrganizer* nameGiver = new SystematicUncertaintyOrganizer();
  
  // Histograms for uncertainty estimation
  TH1D* nominalEnergyEnergyCorrelators[nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* jetEnergyResolutionUncertaintyCorrelators[nJetPtBinsEEC][nTrackPtBinsEEC][2];
  TH1D* jetEnergyScaleUncertaintyCorrelators[nJetPtBinsEEC][nTrackPtBinsEEC][2];
  TH1D* jetPtPriorUncertaintyCorrelators[nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* backgroundSubtractionUncertaintyCorrelators[nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* singleTrackEfficiencyUncertaintyCorrelators[nJetPtBinsEEC][nTrackPtBinsEEC][2];
  TH1D* trackPairEfficiencyUncertaintyCorrelators[nJetPtBinsEEC][nTrackPtBinsEEC][2];

  // Histograms to hold the systematic uncertainty results
  TH1D* energyEnergyCorrelatorSystematicUncertainties[nJetPtBinsEEC][nTrackPtBinsEEC][SystematicUncertaintyOrganizer::knUncertaintySources];
  
  // Initialize all the histograms to NULL
  for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++) {
    for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++) {
      nominalEnergyEnergyCorrelators[iJetPt][iTrackPt] = NULL;
      jetPtPriorUncertaintyCorrelators[iJetPt][iTrackPt] = NULL;
      backgroundSubtractionUncertaintyCorrelators[iJetPt][iTrackPt] = NULL;
      for(int iFile = 0; iFile < 2; iFile++) {
        jetEnergyResolutionUncertaintyCorrelators[iJetPt][iTrackPt][iFile] = NULL;
        jetEnergyScaleUncertaintyCorrelators[iJetPt][iTrackPt][iFile] = NULL;
        singleTrackEfficiencyUncertaintyCorrelators[iJetPt][iTrackPt][iFile] = NULL;
        trackPairEfficiencyUncertaintyCorrelators[iJetPt][iTrackPt][iFile] = NULL;
      }
      for(int iUncertainty = 0; iUncertainty < SystematicUncertaintyOrganizer::knUncertaintySources; iUncertainty++) {
        energyEnergyCorrelatorSystematicUncertainties[iJetPt][iTrackPt][iUncertainty] = NULL;
      } // Uncertainty loop
    } // Track pT loop
  } // Jet pT loop

  // Read the histograms from the defined range from the file
  int iTrackPtMatched;
  int iJetPtMatched;
  double epsilon = 0.0001;
  int lowAnalysisBin, highAnalysisBin;

  for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++) {
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++) {

      // Read the nominal histograms
      nominalEnergyEnergyCorrelators[iJetPt][iTrackPt] = nominalHistogramManager->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPt, iTrackPt, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

      // Normalize nominal histograms to one in the drawing range
      lowAnalysisBin = nominalEnergyEnergyCorrelators[iJetPt][iTrackPt]->GetXaxis()->FindBin(analysisDeltaR.first + epsilon);
      highAnalysisBin = nominalEnergyEnergyCorrelators[iJetPt][iTrackPt]->GetXaxis()->FindBin(analysisDeltaR.second - epsilon);
      nominalEnergyEnergyCorrelators[iJetPt][iTrackPt]->Scale(1.0 / nominalEnergyEnergyCorrelators[iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

      for(int iFile = 0; iFile < 2; iFile++){

        // Read the jet energy resolution histograms
        iTrackPtMatched = jetEnergyResolutionCard[iFile]->FindBinIndexTrackPtEEC(nominalResultCard->GetBinBordersTrackPtEEC(iTrackPt));
        iJetPtMatched = jetEnergyResolutionCard[iFile]->FindBinIndexJetPtEEC(nominalResultCard->GetBinBordersJetPtEEC(iJetPt));
        jetEnergyResolutionUncertaintyCorrelators[iJetPt][iTrackPt][iFile] = jetEnergyResolutionHistogramManager[iFile]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatched, iTrackPtMatched, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

        // Normalize the jet energy resolution histograms to one
        jetEnergyResolutionUncertaintyCorrelators[iJetPt][iTrackPt][iFile]->Scale(1.0 / jetEnergyResolutionUncertaintyCorrelators[iJetPt][iTrackPt][iFile]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

        // Read the jet energy scale histograms
        iTrackPtMatched = jetEnergyScaleCard[iFile]->FindBinIndexTrackPtEEC(nominalResultCard->GetBinBordersTrackPtEEC(iTrackPt));
        iJetPtMatched = jetEnergyScaleCard[iFile]->FindBinIndexJetPtEEC(nominalResultCard->GetBinBordersJetPtEEC(iJetPt));
        jetEnergyScaleUncertaintyCorrelators[iJetPt][iTrackPt][iFile] = jetEnergyScaleHistogramManager[iFile]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatched, iTrackPtMatched, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

        // Normalize the jet energy scale histograms to one
        jetEnergyScaleUncertaintyCorrelators[iJetPt][iTrackPt][iFile]->Scale(1.0 / jetEnergyScaleUncertaintyCorrelators[iJetPt][iTrackPt][iFile]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
      }

      // Read the jet pT prior histograms
      iTrackPtMatched = jetPtPriorCard->FindBinIndexTrackPtEEC(nominalResultCard->GetBinBordersTrackPtEEC(iTrackPt));
      iJetPtMatched = jetPtPriorCard->FindBinIndexJetPtEEC(nominalResultCard->GetBinBordersJetPtEEC(iJetPt));
      jetPtPriorUncertaintyCorrelators[iJetPt][iTrackPt] = jetPtPriorHistogramManager->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatched, iTrackPtMatched, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

      // Normalize the jet energy resolution histograms to one
      jetPtPriorUncertaintyCorrelators[iJetPt][iTrackPt]->Scale(1.0 / jetPtPriorUncertaintyCorrelators[iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

      // Setup the histogram manager for tracking efficiencies
      iTrackPtMatched = trackEfficiencyCard->FindBinIndexTrackPtEEC(nominalResultCard->GetBinBordersTrackPtEEC(iTrackPt));
      iJetPtMatched = trackEfficiencyCard->FindBinIndexJetPtEEC(nominalResultCard->GetBinBordersJetPtEEC(iJetPt));

      // Read the histograms for single track efficiency
      singleTrackEfficiencyUncertaintyCorrelators[iJetPt][iTrackPt][0] = trackEfficiencyHistogramManager->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationPlus, 0, iJetPtMatched, iTrackPtMatched, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);
      singleTrackEfficiencyUncertaintyCorrelators[iJetPt][iTrackPt][1] = trackEfficiencyHistogramManager->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationMinus, 0, iJetPtMatched, iTrackPtMatched, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

      // Read the histograms for track pair efficiency
      trackPairEfficiencyUncertaintyCorrelators[iJetPt][iTrackPt][0] = trackEfficiencyHistogramManager->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationPlus, 0, iJetPtMatched, iTrackPtMatched, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);
      trackPairEfficiencyUncertaintyCorrelators[iJetPt][iTrackPt][1] = trackEfficiencyHistogramManager->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationMinus, 0, iJetPtMatched, iTrackPtMatched, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

      for(int iFile = 0; iFile < 2; iFile++) {
        // Normalize the single track efficiency histograms to one
        singleTrackEfficiencyUncertaintyCorrelators[iJetPt][iTrackPt][iFile]->Scale(1.0 / singleTrackEfficiencyUncertaintyCorrelators[iJetPt][iTrackPt][iFile]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

        // Normalize the track pair efficiency histograms to one
        trackPairEfficiencyUncertaintyCorrelators[iJetPt][iTrackPt][iFile]->Scale(1.0 / trackPairEfficiencyUncertaintyCorrelators[iJetPt][iTrackPt][iFile]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

      }

      // Read the background subtraction uncertainty histograms
      iTrackPtMatched = backgroundSubtractionCard->FindBinIndexTrackPtEEC(nominalResultCard->GetBinBordersTrackPtEEC(iTrackPt));
      iJetPtMatched = backgroundSubtractionCard->FindBinIndexJetPtEEC(nominalResultCard->GetBinBordersJetPtEEC(iJetPt));
      backgroundSubtractionUncertaintyCorrelators[iJetPt][iTrackPt] = backgroundSubtractionHistogramManager->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatched, iTrackPtMatched, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

      // Normalize the background subtraction histograms to one
      backgroundSubtractionUncertaintyCorrelators[iJetPt][iTrackPt]->Scale(1.0 / backgroundSubtractionUncertaintyCorrelators[iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

    }  // Track pT loop
  }    // Jet pT loop

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
    for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
        energyEnergyCorrelatorSystematicUncertainties[iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kJetEnergyResolution] = findTheDifference(nominalEnergyEnergyCorrelators[iJetPt][iTrackPt], jetEnergyResolutionUncertaintyCorrelators[iJetPt][iTrackPt], 2);

        // Draw example plots on how the uncertainty is obtained
        if(plotExample[SystematicUncertaintyOrganizer::kJetEnergyResolution] && std::binary_search(drawnJetPtBins.begin(), drawnJetPtBins.end(), iJetPt) && std::binary_search(drawnTrackPtBins.begin(), drawnTrackPtBins.end(), iTrackPt)) {
          legendNames[0] = "JER varied response down";
          legendNames[1] = "JER varied response up";

          // Set reasonable ratio zoom
          if(setAutomaticRatioZoom){
            ratioZoom = std::make_pair(0.95,1.05);
          }

          drawIllustratingPlots(drawer, nominalEnergyEnergyCorrelators[iJetPt][iTrackPt], jetEnergyResolutionUncertaintyCorrelators[iJetPt][iTrackPt], 2, iJetPt, iTrackPt, nominalResultCard, legendNames, nameGiver->GetSystematicUncertaintyName(SystematicUncertaintyOrganizer::kJetEnergyResolution), analysisDeltaR, ratioZoom);
        }

      } // Track pT loop
    } // Jet pT loop
  } // Jet energy resolution uncertainty

  // ============================================ //
  //   Uncertainty coming from jet energy scale   //
  // ============================================ //

  if(!skipUncertaintySource[SystematicUncertaintyOrganizer::kJetEnergyScale]){
    for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

        energyEnergyCorrelatorSystematicUncertainties[iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kJetEnergyScale] = findTheDifference(nominalEnergyEnergyCorrelators[iJetPt][iTrackPt], jetEnergyScaleUncertaintyCorrelators[iJetPt][iTrackPt], 2);

        // Draw example plots on how the uncertainty is obtained
        if(plotExample[SystematicUncertaintyOrganizer::kJetEnergyScale] && std::binary_search(drawnJetPtBins.begin(), drawnJetPtBins.end(), iJetPt) && std::binary_search(drawnTrackPtBins.begin(), drawnTrackPtBins.end(), iTrackPt)) {
          legendNames[0] = "JEC varied response down";
          legendNames[1] = "JEC varied response up";

          // Set reasonable ratio zoom
          if(setAutomaticRatioZoom){
            ratioZoom = std::make_pair(0.9,1.1);
          }

          drawIllustratingPlots(drawer, nominalEnergyEnergyCorrelators[iJetPt][iTrackPt], jetEnergyScaleUncertaintyCorrelators[iJetPt][iTrackPt], 2, iJetPt, iTrackPt, nominalResultCard, legendNames, nameGiver->GetSystematicUncertaintyName(SystematicUncertaintyOrganizer::kJetEnergyScale), analysisDeltaR, ratioZoom);
        }

      } // Track pT loop
    } // Jet pT loop
  } // Jet energy scale uncertainty

  // ==================================================== //
  //   Uncertainty coming from jet pT prior in unfolding  //
  // ==================================================== //

  if(!skipUncertaintySource[SystematicUncertaintyOrganizer::kUnfoldingTruth]){
    for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

        energyEnergyCorrelatorSystematicUncertainties[iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kUnfoldingTruth] = findTheDifference(nominalEnergyEnergyCorrelators[iJetPt][iTrackPt], jetPtPriorUncertaintyCorrelators[iJetPt][iTrackPt]);

        // Draw example plots on how the uncertainty is obtained
        if(plotExample[SystematicUncertaintyOrganizer::kUnfoldingTruth] && std::binary_search(drawnJetPtBins.begin(), drawnJetPtBins.end(), iJetPt) && std::binary_search(drawnTrackPtBins.begin(), drawnTrackPtBins.end(), iTrackPt)) {
          legendNames[0] = "Weighted prior";

          // Set reasonable ratio zoom
          if(setAutomaticRatioZoom){
            ratioZoom = std::make_pair(0.95,1.05);
          }

          drawIllustratingPlots(drawer, nominalEnergyEnergyCorrelators[iJetPt][iTrackPt], jetPtPriorUncertaintyCorrelators[iJetPt][iTrackPt], iJetPt, iTrackPt, nominalResultCard, legendNames[0], nameGiver->GetSystematicUncertaintyName(SystematicUncertaintyOrganizer::kUnfoldingTruth), analysisDeltaR, ratioZoom);
        }

      } // Track pT loop
    } // Jet pT loop
  } // Jet energy scale uncertainty

  // ================================================== //
  //   Uncertainty coming from background subtraction   //
  // ================================================== //
  
  if(!skipUncertaintySource[SystematicUncertaintyOrganizer::kBackgroundSubtraction]){
    for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

        energyEnergyCorrelatorSystematicUncertainties[iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kBackgroundSubtraction] = findTheDifference(nominalEnergyEnergyCorrelators[iJetPt][iTrackPt], backgroundSubtractionUncertaintyCorrelators[iJetPt][iTrackPt]);

        // Draw example plots on how the uncertainty is obtained
        if(plotExample[SystematicUncertaintyOrganizer::kBackgroundSubtraction] && std::binary_search(drawnJetPtBins.begin(), drawnJetPtBins.end(), iJetPt) && std::binary_search(drawnTrackPtBins.begin(), drawnTrackPtBins.end(), iTrackPt)) {
          legendNames[0] = "Bg from peripharal P+H";

          // Set reasonable ratio zoom
          if(setAutomaticRatioZoom){
            ratioZoom = std::make_pair(0.99,1.01);
          }

          drawIllustratingPlots(drawer, nominalEnergyEnergyCorrelators[iJetPt][iTrackPt], backgroundSubtractionUncertaintyCorrelators[iJetPt][iTrackPt], iJetPt, iTrackPt, nominalResultCard, legendNames[0], nameGiver->GetSystematicUncertaintyName(SystematicUncertaintyOrganizer::kBackgroundSubtraction), analysisDeltaR, ratioZoom);
        }

      } // Track pT loop
    } // Jet pT loop
  } // Background subtraction uncertainty

  // =================================================== //
  //   Uncertainty coming from single track efficiency   //
  // =================================================== //
  
  if(!skipUncertaintySource[SystematicUncertaintyOrganizer::kSingleTrackEfficiency]){
    for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

        energyEnergyCorrelatorSystematicUncertainties[iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kSingleTrackEfficiency] = findTheDifference(nominalEnergyEnergyCorrelators[iJetPt][iTrackPt], singleTrackEfficiencyUncertaintyCorrelators[iJetPt][iTrackPt], 2);        

        // Draw example plots on how the uncertainty is obtained
        if(plotExample[SystematicUncertaintyOrganizer::kSingleTrackEfficiency] && std::binary_search(drawnJetPtBins.begin(), drawnJetPtBins.end(), iJetPt) && std::binary_search(drawnTrackPtBins.begin(), drawnTrackPtBins.end(), iTrackPt)) {
          legendNames[0] = "Efficiency - 2.2%";
          legendNames[1] = "Efficiency + 2.2%";

          // Set reasonable ratio zoom
          if(setAutomaticRatioZoom){
            ratioZoom = std::make_pair(0.99999,1.00001);
          }

          drawIllustratingPlots(drawer, nominalEnergyEnergyCorrelators[iJetPt][iTrackPt], singleTrackEfficiencyUncertaintyCorrelators[iJetPt][iTrackPt], 2, iJetPt, iTrackPt, nominalResultCard, legendNames, nameGiver->GetSystematicUncertaintyName(SystematicUncertaintyOrganizer::kSingleTrackEfficiency), analysisDeltaR, ratioZoom);
        }

      } // Track pT loop
    } // Jet pT loop
  } // Track selection uncertainty

  // =================================================== //
  //    Uncertainty coming from track pair efficiency    //
  // =================================================== //

  if(!skipUncertaintySource[SystematicUncertaintyOrganizer::kTrackPairEfficiency]){
    for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

        energyEnergyCorrelatorSystematicUncertainties[iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kTrackPairEfficiency] = findTheDifference(nominalEnergyEnergyCorrelators[iJetPt][iTrackPt], trackPairEfficiencyUncertaintyCorrelators[iJetPt][iTrackPt], 2);

        // Draw example plots on how the uncertainty is obtained
        if(plotExample[SystematicUncertaintyOrganizer::kTrackPairEfficiency] && std::binary_search(drawnJetPtBins.begin(), drawnJetPtBins.end(), iJetPt) && std::binary_search(drawnTrackPtBins.begin(), drawnTrackPtBins.end(), iTrackPt)){
          legendNames[0] = "Pair efficiency - 4.4%";
          legendNames[1] = "Pair efficiency + 4.4%";

          // Set reasonable ratio zoom
          if(setAutomaticRatioZoom){
            ratioZoom = std::make_pair(0.95,1.05);
          }

          drawIllustratingPlots(drawer, nominalEnergyEnergyCorrelators[iJetPt][iTrackPt], trackPairEfficiencyUncertaintyCorrelators[iJetPt][iTrackPt], 2, iJetPt, iTrackPt, nominalResultCard, legendNames, nameGiver->GetSystematicUncertaintyName(SystematicUncertaintyOrganizer::kTrackPairEfficiency), analysisDeltaR, ratioZoom);
        }

      } // Track pT loop
    } // Jet pT loop
  } // Track selection uncertainty

  // ============================================= //
  //   Add all uncertainty sources in quadrature   //
  // ============================================= //
  
  double sumOfSquares;
  for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

      // Clone the uncertainty histogram from the nominal one
      energyEnergyCorrelatorSystematicUncertainties[iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kAll] = (TH1D*)nominalEnergyEnergyCorrelators[iJetPt][iTrackPt]->Clone(Form("summedUncertainty%d%d", iJetPt, iTrackPt));

      // Calculate a sum of squared from all evaluated uncertainties
      for(int iBin = 1; iBin <= nominalEnergyEnergyCorrelators[iJetPt][iTrackPt]->GetNbinsX(); iBin++){
        sumOfSquares = 0;
        for(int iUncertainty = 0; iUncertainty < SystematicUncertaintyOrganizer::kAll; iUncertainty++) {
          if(energyEnergyCorrelatorSystematicUncertainties[iJetPt][iTrackPt][iUncertainty] != NULL) {
            sumOfSquares = sumOfSquares + TMath::Power(energyEnergyCorrelatorSystematicUncertainties[iJetPt][iTrackPt][iUncertainty]->GetBinError(iBin), 2);
          }
        }  // Uncertainty type loop
        energyEnergyCorrelatorSystematicUncertainties[iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kAll]->SetBinError(iBin, TMath::Sqrt(sumOfSquares));
      } // Bin loop
    } // Track pT loop
  } // Jet pT loop

  // =========================================================== //
  //   Write the systematic uncertainty histograms into a file   //
  // =========================================================== //
  
  // If an output file name is given, put the total uncertainties to an output file
  if(outputFileName.EndsWith(".root")){
    TFile* outputFile = new TFile(outputFileName,"UPDATE");
    TString saveName;

    for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++) {
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++) {
        for(int iUncertainty = 0; iUncertainty < SystematicUncertaintyOrganizer::knUncertaintySources; iUncertainty++) {
          if(energyEnergyCorrelatorSystematicUncertainties[iJetPt][iTrackPt][iUncertainty] != NULL) {
            saveName = Form("systematicUncertainty_%s_C0J%dT%d", nameGiver->GetSystematicUncertaintyName(iUncertainty).Data(), iJetPt, iTrackPt);
            energyEnergyCorrelatorSystematicUncertainties[iJetPt][iTrackPt][iUncertainty]->Write(saveName);
          }
        } // Uncertainty loop
      } // Track pT loop
    } // Jet pT loop

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
 *  const int iJetPt = Jet pT index of the considered bin
 *  const int iTrackPt = Track pT index of the considered bin
 *  EECCard* card = Card used to interpret the bin index
 *  TString comparisonLegend[] = An array of strings describing each result variation
 *  TString plotName = String added to saved plots. If left empty, the plots are not saved into files.
 *  std::pair<double, double> drawingRange = Drawing range for x-axis
 *  std::pair<double, double> ratioZoom = Y-axis zoom for the ratio
 */
void drawIllustratingPlots(JDrawer* drawer, TH1D* nominalResult, TH1D* variedResult[], const int nVariations, const int iJetPt, const int iTrackPt, EECCard* card, TString comparisonLegend[], TString plotName, std::pair<double, double> drawingRange, std::pair<double, double> ratioZoom){
  
  const int markers[] = {kFullDiamond, kFullDoubleDiamond, kFullCross, kFullFourTrianglesPlus, kFullSquare, kFullStar};
  const int colors[] = {kBlue, kRed, kGreen+3, kMagenta, kCyan, kViolet, kBlack};
  
  // Interpret the given binning information
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
    gPad->GetCanvas()->SaveAs(Form("figures/systematicUncertainty_%s_pp%s%s.pdf", plotName.Data(), compactJetPtString.Data(), compactTrackPtString.Data()));
  }
}

/*
 * Draw plots illustrating how systematic uncertainties are estimated form a certain source
 *
 *  JDrawer* drawer = JDrawer doing the dirty work in drawing
 *  TH1D* nominalResult = Histogram containing nominal results
 *  TH1D* variedResult = Histogram containing variation of results around the nominal one
 *  const int iJetPt = Jet pT index of the considered bin
 *  const int iTrackPt = Track pT index of the considered bin
 *  EECCard* card = Card used to interpret the binning information
 *  TString comparisonLegend = String describing the variation
 *  TString plotName = String added to saved plots. If left empty, the plots are not saved into files.
 *  std::pair<double, double> drawingRange = Drawing range for x-axis
 *  std::pair<double, double> ratioZoom = Y-axis zoom for the ratio
 */
void drawIllustratingPlots(JDrawer* drawer, TH1D* nominalResult, TH1D* variedResult, const int iJetPt, const int iTrackPt, EECCard* card, TString comparisonLegend, TString plotName, std::pair<double, double> drawingRange, std::pair<double, double> ratioZoom){
  TH1D* variedResultArray[1] = {variedResult};
  TString comparisonLegendArray[1] = {comparisonLegend};
  drawIllustratingPlots(drawer, nominalResult, variedResultArray, 1, iJetPt, iTrackPt, card, comparisonLegendArray, plotName, drawingRange, ratioZoom);
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

