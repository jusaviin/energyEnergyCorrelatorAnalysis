#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Macro for studying the reflected cone background estimation using MC
 */
void studyReflectedConeBackground(){
  
  enum enumDataType{kPythiaHydjetSimulation, kMinBiasHydjetSimulation, knDataTypes};
  
  // Input files: index 0 = Pythia+Hydjet simulation, index 1 = minimum bias Hydjet simulation
  TString inputFileName[knDataTypes] = {"data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJet_wtaAxis_noTrigger_preprocessed_2022-10-14.root", "data/MinBiasHydjet_RecoGen_eecAnalysis_akFlowJet_firstMinBiasScan_noTrigger_preprocessed_2022-10-10.root"};
  // data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJet_wtaAxis_noTrigger_preprocessed_2022-10-14.root
  // data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJet_fakeFakeReflectedCone_noTrigger_preprocessed_2022-09-23.root
  // data/PbPbMC2018_GenGen_eecAnalysis_genJet_fakeFakeReflectedCone_noTrigger_preprocessed_2022-09-30.root
  // data/MinBiasHydjet_RecoGen_eecAnalysis_akFlowJet_firstMinBiasScan_noTrigger_preprocessed_2022-10-10.root
  
  // Open the input files and read the analysis cards
  TFile *inputFile[knDataTypes];
  EECCard *card[knDataTypes];
  for(int iDataType = 0; iDataType < knDataTypes; iDataType++){
    
    // Open the input file
    inputFile[iDataType] = TFile::Open(inputFileName[iDataType]);
    
    if(inputFile[iDataType] == NULL){
      cout << "Error! The file " << inputFileName[iDataType].Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }
    
    // Load the card from the file
    card[iDataType] = new EECCard(inputFile[iDataType]);
  }
  
  // Require matching centrality and track pT bins between Pythia+Hydjet and minimum bias Hydjet files
  const double epsilon = 0.00001;
  
  const int nCentralityBins = card[kPythiaHydjetSimulation]->GetNCentralityBins();
  if(nCentralityBins != card[kMinBiasHydjetSimulation]->GetNCentralityBins()){
    cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
    return;
  }
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    if(TMath::Abs(card[kPythiaHydjetSimulation]->GetLowBinBorderCentrality(iCentrality) - card[kMinBiasHydjetSimulation]->GetLowBinBorderCentrality(iCentrality)) > epsilon){
      cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
    if(TMath::Abs(card[kPythiaHydjetSimulation]->GetHighBinBorderCentrality(iCentrality) - card[kMinBiasHydjetSimulation]->GetHighBinBorderCentrality(iCentrality)) > epsilon){
      cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
  }
  
  const int nTrackPtBinsEEC = card[kPythiaHydjetSimulation]->GetNTrackPtBinsEEC();
  if(nTrackPtBinsEEC != card[kMinBiasHydjetSimulation]->GetNTrackPtBinsEEC()){
    cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
    return;
  }
  for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
    if(TMath::Abs(card[kPythiaHydjetSimulation]->GetLowBinBorderTrackPtEEC(iTrackPt) - card[kMinBiasHydjetSimulation]->GetLowBinBorderTrackPtEEC(iTrackPt)) > epsilon){
      cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
    if(TMath::Abs(card[kPythiaHydjetSimulation]->GetHighBinBorderTrackPtEEC(iTrackPt) - card[kMinBiasHydjetSimulation]->GetHighBinBorderTrackPtEEC(iTrackPt)) > epsilon){
      cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
  }
  
  // ====================================================
  //    Binning configuration for the background study
  // ====================================================
  
  // Find the number of bins from the card
  const int nJetPtBinsEEC[knDataTypes] = {card[kPythiaHydjetSimulation]->GetNJetPtBinsEEC(), card[kMinBiasHydjetSimulation]->GetNJetPtBinsEEC()};
  
  // Default binning ranges for reference
  // centrality = {0,10,30,50,90}
  // track pT =  {0.7,1,2,3,4,8,12,300}
  // jet pT in energy-energy correlators = {120,140,160,180,200,300,500,5030}
  // track pT in energy-energy correlators = {0.7,1,2,3,4,6,300}
  
  // Bin range to be integrated
  int firstStudiedCentralityBin = 0;
  int lastStudiedCentralityBin = 0;
  
  int firstStudiedJetPtBinEEC[knDataTypes] = {0,0};
  int lastStudiedJetPtBinEEC[knDataTypes] = {nJetPtBinsEEC[kPythiaHydjetSimulation], nJetPtBinsEEC[kMinBiasHydjetSimulation]}; // Note: Jets integrated over all pT ranges are in nJetPtBinsEEC bin
  
  int firstStudiedTrackPtBinEEC = 5;
  int lastStudiedTrackPtBinEEC = 5;
  
  // Select the types of energy-energy correlators are studied
  bool studyEnergyEnergyCorrelator[EECHistogramManager::knEnergyEnergyCorrelatorTypes];
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator] = true;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPt] = false;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorUncorrected] = false;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPtUncorrected] = false;
  
  // Logarithmic axes
  const bool logDeltaR = true;
  const bool logEEC = true;
  
  // Figure saving
  const bool saveFigures = false;  // Save figures
  const char* saveComment = "";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures
  
  // Instead of normalizing to the tail of the distribution, best match the background in the region where it is the most dominant
  bool optimalNormalization = true;
  
  // Create and setup a new histogram manager to project and handle the histograms
  EECHistogramManager *histograms[knDataTypes];
  for(int iDataType = 0; iDataType < knDataTypes; iDataType++){
    histograms[iDataType] = new EECHistogramManager(inputFile[iDataType],card[iDataType]);
    
    // Choose the energy-energy correlator types to load
    histograms[iDataType]->SetLoadEnergyEnergyCorrelators(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator]);
    histograms[iDataType]->SetLoadEnergyEnergyCorrelatorsJetPt(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPt]);
    histograms[iDataType]->SetLoadEnergyEnergyCorrelatorsUncorrected(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorUncorrected]);
    histograms[iDataType]->SetLoadEnergyEnergyCorrelatorsJetPtUncorrected(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPtUncorrected]);
    
    // Choose the bin ranges
    histograms[iDataType]->SetCentralityBinRange(firstStudiedCentralityBin,lastStudiedCentralityBin);
    histograms[iDataType]->SetJetPtBinRangeEEC(firstStudiedJetPtBinEEC[iDataType],lastStudiedJetPtBinEEC[iDataType]);
    histograms[iDataType]->SetTrackPtBinRangeEEC(firstStudiedTrackPtBinEEC,lastStudiedTrackPtBinEEC);
    
    // Load the histograms from the file
    histograms[iDataType]->LoadProcessedHistograms();
  }
  
  // Number of bins from the large deltaR used for the normalization of reflected cone background
  const int nBackgroundNormalizationBins = 9;
  
  // Indices for different types of energy-energy correlators
  enum enumCorrelatorPairingTypes{kTotalEEC, kSignalEEC, kSignalFakeEEC, kFakeFakeEEC, kBackgroundEEC, knPairingTypesEEC};
  enum enumReflectedConeType{kPairSignalReflectedCone, kPairOnlyReflectedCone, kCombineReflectedCone, kSubtractSignalReflectedCone, kSubtractCombined, knReflectedConeTypes};
  enum enumMinBiasHistogramStyle{kPureMinBias, kReflectedConeCombinedMinBias, kMinBiasBackgroundSubtracted, knMinBiasHistograms};
  enum enumMinBiasNormalizationStyle{kNormalizeToReflectedCone, kNormalizeToGeneratorLevel, knNormalizationStyles};
  
  TString matchString[knNormalizationStyles] = {"Ref cone match", "Gen level match"};
  
  // =========================================================================================================================== //
  // Select which ratios should be made from the energy-energy correlators and reflected cone histograms and their configuration //
  // =========================================================================================================================== //
  const int nRatioTypes = 6;
  std::pair<int,int> ratioIndex[nRatioTypes];
  bool drawComparisonType[nRatioTypes];
  const char* legendTextEnergyEnergyCorrelator[nRatioTypes];
  const char* legendTextReflectedCone[nRatioTypes];
  const char* ratioText[nRatioTypes];
  std::pair<double,double> yRange[nRatioTypes];
  std::pair<double,double> ratioZoom[nRatioTypes];
  const char* saveName[nRatioTypes];
  bool addSignalToTotalRatio[nRatioTypes];
  
  // Index 0: Compare signal+fake to corresponding reflected cone distribution
  drawComparisonType[0] = false;
  ratioIndex[0] = std::make_pair(kSignalFakeEEC, kPairSignalReflectedCone);
  legendTextEnergyEnergyCorrelator[0] = "Signal+fake pairs";
  legendTextReflectedCone[0] = "Jet+ref";
  ratioText[0] = "(Jet+ref) / BG";
  yRange[0] = std::make_pair(0.0005, 2);
  ratioZoom[0] = std::make_pair(0, 2);
  saveName[0] = "pythiaHydjetToSignalReflectedCone";
  addSignalToTotalRatio[0] = false;
  
  // Index 1: Compare fake+fake to corresponding reflected cone distribution
  drawComparisonType[1] = false;
  ratioIndex[1] = std::make_pair(kFakeFakeEEC, kPairOnlyReflectedCone);
  legendTextEnergyEnergyCorrelator[1] = "Fake+fake pairs";
  legendTextReflectedCone[1] = "Ref+ref";
  ratioText[1] = "(Ref+ref) / BG";
  yRange[1] = std::make_pair(0.0005, 2);
  ratioZoom[1] = std::make_pair(0, 50);
  saveName[1] = "hydjetHydjetToOnlyReflectedCone";
  addSignalToTotalRatio[1] = false;
  
  // Index 2: Compare total background to signal+fake reflected cone
  drawComparisonType[2] = false;
  ratioIndex[2] = std::make_pair(kBackgroundEEC, kPairSignalReflectedCone);
  legendTextEnergyEnergyCorrelator[2] = "All background pairs";
  legendTextReflectedCone[2] = "Jet+ref";
  ratioText[2] = "(Jet+ref) / BG";
  yRange[2] = std::make_pair(0.0005, 2);
  ratioZoom[2] = std::make_pair(0, 2);
  saveName[2] = "allBackgroundToSignalReflectedCone";
  addSignalToTotalRatio[2] = false;
  
  // Index 3: Compare total background to combined reflected cone
  drawComparisonType[3] = false;
  ratioIndex[3] = std::make_pair(kBackgroundEEC, kCombineReflectedCone);
  legendTextEnergyEnergyCorrelator[3] = "All background pairs";
  legendTextReflectedCone[3] = "Reflected";
  ratioText[3] = "Reflected / BG";
  yRange[3] = std::make_pair(0.0005, 2);
  ratioZoom[3] = std::make_pair(0, 2);
  saveName[3] = "allBackgroundToCombinedReflectedCone";
  addSignalToTotalRatio[3] = false;
  
  // Index 4: Compare signal to signal+fake reflected cone subtracted total distribution
  drawComparisonType[4] = false;
  ratioIndex[4] = std::make_pair(kSignalEEC, kSubtractSignalReflectedCone);
  legendTextEnergyEnergyCorrelator[4] = "Signal pairs";
  legendTextReflectedCone[4] = "Jet+ref";
  ratioText[4] = "BGsub / signal";
  yRange[4] = std::make_pair(0.001, 50);
  ratioZoom[4] = std::make_pair(0, 2);
  saveName[4] = "signalToSignalReflectedConeSubtracted";
  addSignalToTotalRatio[4] = false;
  
  // Index 5: Compare signal to combined reflected cone subtracted total distribution
  drawComparisonType[5] = true;
  ratioIndex[5] = std::make_pair(kSignalEEC, kSubtractCombined);
  legendTextEnergyEnergyCorrelator[5] = "Signal pairs";
  legendTextReflectedCone[5] = "BGsub";
  ratioText[5] = "BGsub / signal";
  yRange[5] = std::make_pair(0.001, 50);
  ratioZoom[5] = std::make_pair(0.5, 1.5);
  saveName[5] = "signalToReflectedConeSubtracted";
  addSignalToTotalRatio[5] = false;
  
  // ========================================================================================================================= //
  // Select which ratios should be made from the energy-energy correlators and minimum bias histograms and their configuration //
  // ========================================================================================================================= //
  const int nMinBiasRatioTypes = 3;
  std::pair<int,int> minBiasRatioIndex[nMinBiasRatioTypes];
  bool drawMinBiasComparisonType[nMinBiasRatioTypes];
  const char* minBiasLegendTextEnergyEnergyCorrelator[nMinBiasRatioTypes];
  const char* minBiasRatioText[nMinBiasRatioTypes];
  std::pair<double,double> yRangeMinBias[nMinBiasRatioTypes];
  std::pair<double,double> ratioZoomMinBias[nMinBiasRatioTypes];
  const char* saveNameMinBias[nMinBiasRatioTypes];
  bool addSignalToTotalRatioMinBias[nMinBiasRatioTypes];
  bool addMatchString[nMinBiasRatioTypes];
  double legendX1MinBias[nMinBiasRatioTypes];
  
  // Index 0: Compare fake+fake distribution to corresponding minimum bias distributions
  drawMinBiasComparisonType[0] = false;
  minBiasRatioIndex[0] = std::make_pair(kFakeFakeEEC, kPureMinBias);
  minBiasLegendTextEnergyEnergyCorrelator[0] = "Fake+fake pairs";
  minBiasRatioText[0] = "MinBias/FakeFake";
  yRangeMinBias[0] = std::make_pair(0.0005, 2);
  ratioZoomMinBias[0] = std::make_pair(0, 2);
  saveNameMinBias[0] = "pureMinBias";
  addSignalToTotalRatioMinBias[0] = false;
  addMatchString[0] = false;
  legendX1MinBias[0] = 0.27;
  
  // Index 1: Compare background distribution to one combined from reflected cone and minimum bias distributions
  drawMinBiasComparisonType[1] = false;
  minBiasRatioIndex[1] = std::make_pair(kBackgroundEEC, kReflectedConeCombinedMinBias);
  minBiasLegendTextEnergyEnergyCorrelator[1] = "Background pairs";
  minBiasRatioText[1] = "Combined/BG";
  yRangeMinBias[1] = std::make_pair(0.005, 5);
  ratioZoomMinBias[1] = std::make_pair(0, 2);
  saveNameMinBias[1] = "minBiasBackgroundCombined";
  addSignalToTotalRatioMinBias[1] = false;
  addMatchString[1] = true;
  legendX1MinBias[1] = 0.27;
  
  // Index 2: Compare fake+fake distribution to corresponding minimum bias distributions
  drawMinBiasComparisonType[2] = false;
  minBiasRatioIndex[2] = std::make_pair(kSignalEEC, kMinBiasBackgroundSubtracted);
  minBiasLegendTextEnergyEnergyCorrelator[2] = "Signal pairs";
  minBiasRatioText[2] = "BG Sub/Signal";
  yRangeMinBias[2] = std::make_pair(0.005, 20);
  ratioZoomMinBias[2] = std::make_pair(0, 2);
  saveNameMinBias[2] = "minBiasBackgrounsSubtracted";
  addSignalToTotalRatioMinBias[2] = false;
  addMatchString[2] = true;
  legendX1MinBias[2] = 0.15;
  
  // Energy-energy correlator histograms separated by subevents from the Pythia+Hydjet simulation
  TH1D* hEnergyEnergyCorrelator[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC[kPythiaHydjetSimulation]+1][nTrackPtBinsEEC][knPairingTypesEEC];
  
  // Reflected cone energy-energy correlators
  TH1D* hReflectedCone[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC[kPythiaHydjetSimulation]+1][nTrackPtBinsEEC][nBackgroundNormalizationBins+knPairingTypesEEC][knReflectedConeTypes];
  
  // Histograms for all different ratios
  TH1D* hRatio[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC[kPythiaHydjetSimulation]+1][nTrackPtBinsEEC][nBackgroundNormalizationBins+knPairingTypesEEC][nRatioTypes];
  TH1D* hSignalToTotalRatio[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC[kPythiaHydjetSimulation]+1][nTrackPtBinsEEC];
  
  // Energy-energy correlators from minimum bias Hydjet simulation
  TH1D* hMinimumBias[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC[kPythiaHydjetSimulation]+1][nJetPtBinsEEC[kMinBiasHydjetSimulation]+1][nTrackPtBinsEEC][knNormalizationStyles][knMinBiasHistograms];
  TH1D* hMinimumBiasRatio[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC[kPythiaHydjetSimulation]+1][nJetPtBinsEEC[kMinBiasHydjetSimulation]+1][nTrackPtBinsEEC][knNormalizationStyles][nMinBiasRatioTypes];
  
  // Initialize the energy-energy correlator histogram array to NULL
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC[kPythiaHydjetSimulation]+1; iJetPt++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
          hSignalToTotalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = NULL;
          for(int iPairingType = 0; iPairingType < knPairingTypesEEC; iPairingType++){
            hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iPairingType] = NULL;
          } // Pairing type loop
          for(int iNormalization = 0; iNormalization < nBackgroundNormalizationBins+knPairingTypesEEC; iNormalization++){
            for(int iReflectedConeType = 0; iReflectedConeType < knReflectedConeTypes; iReflectedConeType++){
              hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iReflectedConeType] = NULL;
            } // Reflected cone loop
            for(int iRatio = 0; iRatio < nRatioTypes; iRatio++){
              hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iRatio] = NULL;
            } // Ratio loop
          } // Normalization region loop
          for(int iJetPtMinBias = 0; iJetPtMinBias < nJetPtBinsEEC[kMinBiasHydjetSimulation]+1; iJetPtMinBias++){
            for(int iNormalization = 0; iNormalization < knNormalizationStyles; iNormalization++){
              for(int iMinBias = 0; iMinBias < knMinBiasHistograms; iMinBias++){
                hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalization][iMinBias] = NULL;
              }
              for(int iRatio = 0; iRatio < nMinBiasRatioTypes; iRatio++){
                hMinimumBiasRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalization][iRatio] = NULL;
              } // Ratio loop
            }
          } // Jet pT loop for minimum bias Hydjet
        } // Track pT loop
      } // Jet pT loop for Pythia+Hydjet
    } // Centrality loop
  } // Energy-energy correlator loop

  // Helper variables
  TH1D *helperHistogram, *helperHistogram2;
  double normalizationFactor;
  int lowIntegralBin, highIntegralBin;
  int studiedEnergyEnergyCorrelatorIndex = -1;
  int referenceIndexNormalization;
  double signalFaketoFakeFakeRatio[knNormalizationStyles];
  
  // Get the histograms from the histogram manager
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only read the selected energy-energy correlator types
    if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
    studiedEnergyEnergyCorrelatorIndex = iEnergyEnergyCorrelator;
    
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iJetPt = firstStudiedJetPtBinEEC[kPythiaHydjetSimulation]; iJetPt <= lastStudiedJetPtBinEEC[kPythiaHydjetSimulation]; iJetPt++){
        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
          
          // ================================================= //
          // Read the histograms from Pythia+Hydjet simulation //
          // ================================================= //
          
          // Histogram with all pair combinations. Normalize it to one
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC] = histograms[kPythiaHydjetSimulation]->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJetPair, EECHistograms::knSubeventTypes);
          normalizationFactor = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC]->Integral("width");
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC]->Scale(1/normalizationFactor);
          
          // Histogram with only signal. Normalize the sum of signal and background to total.
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kSignalEEC] = histograms[kPythiaHydjetSimulation]->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJetPair, EECHistograms::kPythiaPythia);
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kSignalEEC]->Scale(1/normalizationFactor);
          
          // Calculate the signal to total ratio
          hSignalToTotalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC]->Clone(Form("signalToTotalRatio%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt));
          hSignalToTotalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Divide(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kSignalEEC]);
          
          // Histograms with background contributions. Normalize to the total number of pairs.
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kSignalFakeEEC] = histograms[kPythiaHydjetSimulation]->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJetPair, EECHistograms::kPythiaHydjet);
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kBackgroundEEC] = (TH1D *) hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kSignalFakeEEC]->Clone(Form("energyEnergyCorrelatorBackground%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt));
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kFakeFakeEEC] = histograms[kPythiaHydjetSimulation]->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJetPair, EECHistograms::kHydjetHydjet);
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kBackgroundEEC]->Add(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kFakeFakeEEC]);
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kSignalFakeEEC]->Scale(1/normalizationFactor);
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kFakeFakeEEC]->Scale(1/normalizationFactor);
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kBackgroundEEC]->Scale(1/normalizationFactor);
          
          // Calculate the signal+fake to fake+fake ratio from generator level information
          signalFaketoFakeFakeRatio[kNormalizeToGeneratorLevel] = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kSignalFakeEEC]->Integral("width") / hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kFakeFakeEEC]->Integral("width");
          
          // The jet cone + reflected cone histogram.
          helperHistogram = histograms[kPythiaHydjetSimulation]->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSignalReflectedConePair, EECHistograms::knSubeventTypes);
          
          // The only reflected cone histogram.
          helperHistogram2 = histograms[kPythiaHydjetSimulation]->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kReflectedConePair, EECHistograms::knSubeventTypes);
          
          // Calculate the signal+fake to fake+fake ratio for the reflected cone histograms
          signalFaketoFakeFakeRatio[kNormalizeToReflectedCone] = helperHistogram->Integral("width") / helperHistogram2->Integral("width");
          
          // Normalize the reflected cone histogram to different DeltaR regions of the tail of the total distribution
          highIntegralBin = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC]->GetNbinsX();
          for(int iNormalization = 0; iNormalization < nBackgroundNormalizationBins+knPairingTypesEEC; iNormalization++){
            lowIntegralBin = highIntegralBin - (nBackgroundNormalizationBins - iNormalization - 1);
            if(iNormalization >= nBackgroundNormalizationBins) {
              lowIntegralBin = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC]->GetXaxis()->FindBin(0.3);
              highIntegralBin = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC]->GetXaxis()->FindBin(0.4);
            }
            
            // Normalization for jet cone + reflected cone histogram
            hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][kPairSignalReflectedCone] = (TH1D*) helperHistogram->Clone(Form("reflectedCone%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iNormalization));
            
            // Normalization for only reflected cone histogram
            hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][kPairOnlyReflectedCone] = (TH1D*) helperHistogram2->Clone(Form("reflectedConeOnly%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iNormalization));
            
            // Normalization for total reflected cone histogram, where all the pairs in went through once in each event
            hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][kCombineReflectedCone] = (TH1D*) helperHistogram->Clone(Form("reflectedConeCombined%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iNormalization));
            hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][kCombineReflectedCone]->Add(helperHistogram2);
            
            // Normalization for all reflected cone histograms
            referenceIndexNormalization = kTotalEEC;
            if(iNormalization >= nBackgroundNormalizationBins) referenceIndexNormalization = iNormalization - nBackgroundNormalizationBins;
            for(int iReflectedConeType = kPairSignalReflectedCone; iReflectedConeType <= kCombineReflectedCone; iReflectedConeType++){
              normalizationFactor = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][referenceIndexNormalization]->Integral(lowIntegralBin, highIntegralBin, "width") / hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iReflectedConeType]->Integral(lowIntegralBin, highIntegralBin, "width");
              hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iReflectedConeType]->Scale(normalizationFactor);
            }
            
            // The background subtracted histograms can be calculated from the normalized reflected cone histograms.
            
            // First clone the total distribution as the baseline
            for(int iReflectedConeType = kSubtractSignalReflectedCone; iReflectedConeType <= kSubtractCombined; iReflectedConeType++){
              hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iReflectedConeType] = (TH1D*) hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC]->Clone(Form("backgroundSubtracted%d%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iNormalization, iReflectedConeType));
            }
            
            // Then subtract the correct reflected cone histogram from the total histogram
            hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][kSubtractSignalReflectedCone]->Add(hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][kPairSignalReflectedCone],-1);
            hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][kSubtractCombined]->Add(hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][kCombineReflectedCone],-1);
            
            // After all the reflected cone histograms are ready, the defined ratios can be calculated
            for(int iRatio = 0; iRatio < nRatioTypes; iRatio++){
              hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iRatio] = (TH1D*) hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][ratioIndex[iRatio].second]->Clone(Form("ratio%d%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iNormalization, iRatio));
              hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iRatio]->Divide(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][ratioIndex[iRatio].first]);
            }
            
          } // Normalization region loop
          
          // ======================================================= //
          // Read the histograms from minimum bias Hydjet simulation //
          // ======================================================= //

          // Bins used for the integration in the total normalization
          lowIntegralBin = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC]->GetXaxis()->FindBin(0.3);
          highIntegralBin = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC]->GetXaxis()->FindBin(0.4);
          
          for(int iJetPtMinBias = firstStudiedJetPtBinEEC[kMinBiasHydjetSimulation]; iJetPtMinBias <= lastStudiedJetPtBinEEC[kMinBiasHydjetSimulation]; iJetPtMinBias++){
            for(int iNormalization = 0; iNormalization < knNormalizationStyles; iNormalization++){
              hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalization][kPureMinBias] = (TH1D*) histograms[kMinBiasHydjetSimulation]->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPtMinBias, iTrackPt, EECHistograms::kSameJetPair, EECHistograms::knSubeventTypes)->Clone(Form("pureMinBiasHistogram%d%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iJetPtMinBias, iTrackPt, iNormalization));
              
              // Normalize the minimum bias histograms to match the signal+fake to fake+fake ratio with different estimators
              normalizationFactor = hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0][kPairSignalReflectedCone]->Integral("width") / (hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalization][kPureMinBias]->Integral("width") * signalFaketoFakeFakeRatio[iNormalization]);
              hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalization][kPureMinBias]->Scale(normalizationFactor);
              hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalization][kReflectedConeCombinedMinBias] = (TH1D*) hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0][kPairSignalReflectedCone]->Clone(Form("combinedBackground%d%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iJetPtMinBias, iTrackPt, iNormalization));
              hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalization][kReflectedConeCombinedMinBias]->Add(hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalization][kPureMinBias]);
              
              // After the ratio is matched, normalize the distribution such that it matches the generator level background at 0.3 < DeltaR < 0.4
              normalizationFactor = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kBackgroundEEC]->Integral(lowIntegralBin, highIntegralBin, "width") / hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalization][kReflectedConeCombinedMinBias]->Integral(lowIntegralBin, highIntegralBin, "width");
              hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalization][kReflectedConeCombinedMinBias]->Scale(normalizationFactor);
              
              // For shape comparisons, also normalize the pure minium bias distribution such that it matches the fake+fake background at 0.3 < DeltaR < 0.4
              normalizationFactor = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kFakeFakeEEC]->Integral(lowIntegralBin, highIntegralBin, "width") / hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalization][kPureMinBias]->Integral(lowIntegralBin, highIntegralBin, "width");
              hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalization][kPureMinBias]->Scale(normalizationFactor);
              
              // Now that we have a normalized background estimate, we can do a background subtraction
              hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalization][kMinBiasBackgroundSubtracted] = (TH1D*) hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC]->Clone(Form("minBiasBackgroundSubtracted%d%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iJetPtMinBias, iTrackPt, iNormalization));
              hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalization][kMinBiasBackgroundSubtracted]->Add(hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalization][kReflectedConeCombinedMinBias], -1);
              
              // Calculate all defined ratios with minimum bias histograms
              for(int iRatio = 0; iRatio < nMinBiasRatioTypes; iRatio++){
                hMinimumBiasRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalization][iRatio] = (TH1D*) hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalization][minBiasRatioIndex[iRatio].second]->Clone(Form("minimumBiasRatio%d%d%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iJetPtMinBias, iTrackPt, iNormalization, iRatio));
                hMinimumBiasRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalization][iRatio]->Divide(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][minBiasRatioIndex[iRatio].first]);
              } // Ratio index loop
              
            } // Normalization style loop
          } // Jet pT loop for minimum bias hydjet simulation
          
        } // Track pT loop
      } // Jet pT loop for Pythia+Hydjet simulation
    } // Centrality loop
  } // Energy-energy correlator loop
  
  // Find the binning information for DeltaR bins from a histogram
  const int nDeltaRBins = hEnergyEnergyCorrelator[studiedEnergyEnergyCorrelatorIndex][firstStudiedCentralityBin][firstStudiedJetPtBinEEC[kPythiaHydjetSimulation]][firstStudiedTrackPtBinEEC][kTotalEEC]->GetNbinsX();
  double deltaRBinBorders[nDeltaRBins+1];
  deltaRBinBorders[0] = hEnergyEnergyCorrelator[studiedEnergyEnergyCorrelatorIndex][firstStudiedCentralityBin][firstStudiedJetPtBinEEC[kPythiaHydjetSimulation]][firstStudiedTrackPtBinEEC][kTotalEEC]->GetXaxis()->GetBinLowEdge(1);
  for(int iBin = 1; iBin <= nDeltaRBins; iBin++){
    deltaRBinBorders[iBin] = hEnergyEnergyCorrelator[studiedEnergyEnergyCorrelatorIndex][firstStudiedCentralityBin][firstStudiedJetPtBinEEC[kPythiaHydjetSimulation]][firstStudiedTrackPtBinEEC][kTotalEEC]->GetXaxis()->GetBinUpEdge(iBin);
  }
  
  
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  
  if(logDeltaR) drawer->SetLogX(true);
  
  TLegend *legend, *ptLegend;
  TString centralityString, trackPtString, jetPtString, jetPtStringMinBias;
  TString compactCentralityString, compactTrackPtString, compactJetPtString, compactJetPtStringMinBias;
  int color[9] = {kBlue, kRed, kGreen+2, kCyan, kMagenta, kOrange-1, kAzure+5, kOrange-2, kGray};
  int firstNormalizationIndex, lastNormalizationIndex, currentNormalizationIndex;
  
  // If optimal background matching is done, change the legend size
  double legendY1 = 0.09;
  double legendY2 = 0.94;
  if(optimalNormalization){
    legendY1 = 0.5;
    legendY2 = 0.86;
  }
  
  // ==========================================================
  // ===   Drawing Pythia-Hydjet to reflected cone ratios   ===
  // ==========================================================
  
  // Loop over all different ratio combinations
  for(int iRatio = 0; iRatio < nRatioTypes; iRatio++){
    
    // Only draw the selected ratio types
    if(!drawComparisonType[iRatio]) continue;
    
    // Loop over all selected histograms
    for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
      
      // Only read the selected energy-energy correlator types
      if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
      
      for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
        
        // Set the centrality information for legends and figure saving
        centralityString = Form("Cent: %.0f-%.0f%%",histograms[kPythiaHydjetSimulation]->GetCentralityBinBorder(iCentrality),histograms[kPythiaHydjetSimulation]->GetCentralityBinBorder(iCentrality+1));
        compactCentralityString = Form("_C=%.0f-%.0f",histograms[kPythiaHydjetSimulation]->GetCentralityBinBorder(iCentrality),histograms[kPythiaHydjetSimulation]->GetCentralityBinBorder(iCentrality+1));
        
        for(int iJetPt = firstStudiedJetPtBinEEC[kPythiaHydjetSimulation]; iJetPt <= lastStudiedJetPtBinEEC[kPythiaHydjetSimulation]; iJetPt++){
          
          // Set the jet pT information for legends and figure saving
          if(iJetPt == histograms[kPythiaHydjetSimulation]->GetNJetPtBinsEEC()){
            jetPtString = Form("Jet p_{T} > %.0f", histograms[kPythiaHydjetSimulation]->GetCard()->GetJetPtCut());
            compactJetPtString = "";
          } else {
            jetPtString = Form("%.0f < jet p_{T} < %.0f", histograms[kPythiaHydjetSimulation]->GetJetPtBinBorderEEC(iJetPt), histograms[kPythiaHydjetSimulation]->GetJetPtBinBorderEEC(iJetPt+1));
            compactJetPtString = Form("_J=%.0f-%.0f", histograms[kPythiaHydjetSimulation]->GetJetPtBinBorderEEC(iJetPt), histograms[kPythiaHydjetSimulation]->GetJetPtBinBorderEEC(iJetPt+1));
          }
          
          for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
            
            // Set the track pT information for legends and figure saving
            trackPtString = Form("%.1f < track p_{T}",histograms[kPythiaHydjetSimulation]->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString = Form("_T%.1f",histograms[kPythiaHydjetSimulation]->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString.ReplaceAll(".","v");
            
            // Create a legend for the figure
            legend = new TLegend(0.19,legendY1,0.39,legendY2);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->AddEntry((TObject*) 0, histograms[kPythiaHydjetSimulation]->GetCard()->GetAlternativeDataType().Data(), "");
            legend->AddEntry((TObject*) 0, centralityString.Data(),"");
            legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
            legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
            
            // Create a new canvas for the plot
            drawer->CreateSplitCanvas();
            
            // Logarithmic EEC axis
            if(logEEC) drawer->SetLogY(true);
            
            // For logarithmic drawing, cannot go down to zero
            if(logDeltaR){
              hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][ratioIndex[iRatio].first]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
            }
            
            // Set good y-ranges for plotting
            hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][ratioIndex[iRatio].first]->GetYaxis()->SetRangeUser(yRange[iRatio].first, yRange[iRatio].second);
            
            // Draw the background histogram to the upper canvas
            hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][ratioIndex[iRatio].first]->SetLineColor(kBlack);
            drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][ratioIndex[iRatio].first], "#Deltar", histograms[kPythiaHydjetSimulation]->GetEnergyEnergyCorrelatorAxisName(iEnergyEnergyCorrelator), " ");
            legend->AddEntry(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][ratioIndex[iRatio].first], legendTextEnergyEnergyCorrelator[iRatio], "l");
            
            // Option to draw the total distribution to the same figure
            firstNormalizationIndex = 0;
            if(addSignalToTotalRatio[iRatio]){
              firstNormalizationIndex = 1;
              
              hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC]->SetLineColor(color[0]);
              hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC]->Draw("same");
              legend->AddEntry(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC], "All pairs", "l");
              
            }
            
            // Option to normalize to the integral of the distribution we are comparing to
            lastNormalizationIndex = nBackgroundNormalizationBins-1;
            if(optimalNormalization) lastNormalizationIndex = firstNormalizationIndex;
            
            // Draw different normalization regions for the reflected cone histogram to the same plot
            for(int iNormalization = firstNormalizationIndex; iNormalization <= lastNormalizationIndex; iNormalization++){
              
              // Find the correct index for optimal normalization
              currentNormalizationIndex = iNormalization;
              if(optimalNormalization && (iNormalization == lastNormalizationIndex)) {
                currentNormalizationIndex = nBackgroundNormalizationBins + ratioIndex[iRatio].first;
                if(ratioIndex[iRatio].second == kSubtractSignalReflectedCone) currentNormalizationIndex = nBackgroundNormalizationBins + kSignalFakeEEC;
                if(ratioIndex[iRatio].second == kSubtractCombined) currentNormalizationIndex = nBackgroundNormalizationBins + kBackgroundEEC;
              }
                            
              // For logarithmic drawing, cannot go down to zero
              if(logDeltaR){
                hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][currentNormalizationIndex][ratioIndex[iRatio].second]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
              }
              
              hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][currentNormalizationIndex][ratioIndex[iRatio].second]->SetLineColor(color[iNormalization]);
              hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][currentNormalizationIndex][ratioIndex[iRatio].second]->Draw("same");
              
              if(optimalNormalization){
                legend->AddEntry(hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][currentNormalizationIndex][ratioIndex[iRatio].second], Form("%s, match yield", legendTextReflectedCone[iRatio]), "l");
              } else {
                legend->AddEntry(hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][currentNormalizationIndex][ratioIndex[iRatio].second], Form("%s, #Deltar > %.3f", legendTextReflectedCone[iRatio], deltaRBinBorders[nDeltaRBins-(nBackgroundNormalizationBins-iNormalization)]), "l");
              }
              
            } // Normalization region loop
            
            // Draw the legend to the upper pad
            legend->Draw();
            
            // Linear scale for the ratio
            drawer->SetLogY(false);
            
            // Draw the ratios to the lower pad
            for(int iNormalization = 0; iNormalization <= lastNormalizationIndex; iNormalization++){
              
              // Find the correct index for optimal normalization
              currentNormalizationIndex = iNormalization;
              if(optimalNormalization && (iNormalization == lastNormalizationIndex)) {
                currentNormalizationIndex = nBackgroundNormalizationBins + ratioIndex[iRatio].first;
                if(ratioIndex[iRatio].second == kSubtractSignalReflectedCone) currentNormalizationIndex = nBackgroundNormalizationBins + kSignalFakeEEC;
                if(ratioIndex[iRatio].second == kSubtractCombined) currentNormalizationIndex = nBackgroundNormalizationBins + kBackgroundEEC;
              }
              
              
              // For logarithmic drawing, cannot go down to zero
              if(logDeltaR){
                hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][currentNormalizationIndex][iRatio]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
              }
              
              hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][currentNormalizationIndex][iRatio]->SetLineColor(color[iNormalization]);
              if(iNormalization == 0){
                drawer->SetGridY(true);
                if(addSignalToTotalRatio[iRatio]){
                  
                  if(logDeltaR){
                    hSignalToTotalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
                  }
                  
                  hSignalToTotalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[0]);
                  hSignalToTotalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom[iRatio].first, ratioZoom[iRatio].second);
                  drawer->DrawHistogramToLowerPad(hSignalToTotalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt], "#Deltar", ratioText[iRatio], " ");
                } else {
                  hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][currentNormalizationIndex][iRatio]->GetYaxis()->SetRangeUser(ratioZoom[iRatio].first, ratioZoom[iRatio].second);
                  drawer->DrawHistogramToLowerPad(hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][currentNormalizationIndex][iRatio], "#Deltar", ratioText[iRatio], " ");
                }
                drawer->SetGridY(false);
              } else {
                hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][currentNormalizationIndex][iRatio]->Draw("same");
              }
              
            } // Normalization region loop
            
            // Save the figures to a file
            if(saveFigures){
              gPad->GetCanvas()->SaveAs(Form("figures/%sComparison%s_%s%s%s%s.%s", saveName[iRatio], saveComment, histograms[kPythiaHydjetSimulation]->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
            }
            
          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Energy-energy correlator loop
    
  } // Ratio type loop
 
  // ==========================================================================
  // ===   Drawing histograms including minimum bias Hydjet distributions   ===
  // ==========================================================================
  
  // Loop over all different minimum bias ratio combinations
  for(int iRatio = 0; iRatio < nMinBiasRatioTypes; iRatio++){
    
    // Only draw the selected ratio types
    if(!drawMinBiasComparisonType[iRatio]) continue;
    
    // Loop over all different normalization schemes
    for(int iNormalizationType = 0; iNormalizationType < knNormalizationStyles; iNormalizationType++){
      
      // There is no difference between the normalization schemes for pure hydjet
      if(minBiasRatioIndex[iRatio].second == kPureMinBias && iNormalizationType != 0) continue;
      
      // Loop over all selected histograms
      for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
        
        // Only read the selected energy-energy correlator types
        if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
        
        for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
          
          // Set the centrality information for legends and figure saving
          centralityString = Form("Cent: %.0f-%.0f%%",histograms[kPythiaHydjetSimulation]->GetCentralityBinBorder(iCentrality),histograms[kPythiaHydjetSimulation]->GetCentralityBinBorder(iCentrality+1));
          compactCentralityString = Form("_C=%.0f-%.0f",histograms[kPythiaHydjetSimulation]->GetCentralityBinBorder(iCentrality),histograms[kPythiaHydjetSimulation]->GetCentralityBinBorder(iCentrality+1));
          
          for(int iJetPt = firstStudiedJetPtBinEEC[kPythiaHydjetSimulation]; iJetPt <= lastStudiedJetPtBinEEC[kPythiaHydjetSimulation]; iJetPt++){
            
            // Set the jet pT information for legends and figure saving
            if(iJetPt == histograms[kPythiaHydjetSimulation]->GetNJetPtBinsEEC()){
              jetPtString = Form("Jet p_{T} > %.0f", histograms[kPythiaHydjetSimulation]->GetCard()->GetJetPtCut());
              compactJetPtString = "";
            } else {
              jetPtString = Form("%.0f < jet p_{T} < %.0f", histograms[kPythiaHydjetSimulation]->GetJetPtBinBorderEEC(iJetPt), histograms[kPythiaHydjetSimulation]->GetJetPtBinBorderEEC(iJetPt+1));
              compactJetPtString = Form("_J=%.0f-%.0f", histograms[kPythiaHydjetSimulation]->GetJetPtBinBorderEEC(iJetPt), histograms[kPythiaHydjetSimulation]->GetJetPtBinBorderEEC(iJetPt+1));
            }
            
            for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
              
              // Set the track pT information for legends and figure saving
              trackPtString = Form("%.1f < track p_{T}",histograms[kPythiaHydjetSimulation]->GetTrackPtBinBorderEEC(iTrackPt));
              compactTrackPtString = Form("_T%.1f",histograms[kPythiaHydjetSimulation]->GetTrackPtBinBorderEEC(iTrackPt));
              compactTrackPtString.ReplaceAll(".","v");
              
              // Create a legend for the figure
              legendY1 = 0.04;
              legendY2 = 0.28 + addMatchString[iRatio]*0.06;
              legend = new TLegend(legendX1MinBias[iRatio],legendY1,legendX1MinBias[iRatio]+0.2,legendY2);
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
              legend->AddEntry((TObject*) 0, histograms[kPythiaHydjetSimulation]->GetCard()->GetAlternativeDataType().Data(), "");
              legend->AddEntry((TObject*) 0, centralityString.Data(),"");
              legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
              legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
              if(addMatchString[iRatio]) legend->AddEntry((TObject*) 0, matchString[iNormalizationType].Data(),"");
              
              ptLegend = new TLegend(legendX1MinBias[iRatio]+0.33, legendY1, legendX1MinBias[iRatio]+0.53, legendY1+0.06*(lastStudiedJetPtBinEEC[kMinBiasHydjetSimulation] + 2 - firstStudiedJetPtBinEEC[kMinBiasHydjetSimulation]));
              ptLegend->SetFillStyle(0);ptLegend->SetBorderSize(0);ptLegend->SetTextSize(0.05);ptLegend->SetTextFont(62);
              
              // Create a new canvas for the plot
              drawer->CreateSplitCanvas();
              
              // Logarithmic EEC axis
              if(logEEC) drawer->SetLogY(true);
              
              // For logarithmic drawing, cannot go down to zero
              if(logDeltaR){
                hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][minBiasRatioIndex[iRatio].first]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
              }
              
              // Set good y-ranges for plotting
              hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][minBiasRatioIndex[iRatio].first]->GetYaxis()->SetRangeUser(yRangeMinBias[iRatio].first, yRangeMinBias[iRatio].second);
              
              // Draw the background histogram to the upper canvas
              hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][minBiasRatioIndex[iRatio].first]->SetLineColor(kBlack);
              drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][minBiasRatioIndex[iRatio].first], "#Deltar", histograms[kPythiaHydjetSimulation]->GetEnergyEnergyCorrelatorAxisName(iEnergyEnergyCorrelator), " ");
              ptLegend->AddEntry(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][minBiasRatioIndex[iRatio].first], minBiasLegendTextEnergyEnergyCorrelator[iRatio], "l");
              
              // Option to draw the total distribution to the same figure
              if(addSignalToTotalRatioMinBias[iRatio]){
                hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC]->SetLineColor(color[0]);
                hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC]->Draw("same");
                ptLegend->AddEntry(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC], "All pairs", "l");
                
              }
              
              // Draw background subtracted histograms with fake+fake backgrounds estimated from different jet pT ranges
              for(int iJetPtMinBias = firstStudiedJetPtBinEEC[kMinBiasHydjetSimulation]; iJetPtMinBias <= lastStudiedJetPtBinEEC[kMinBiasHydjetSimulation]; iJetPtMinBias++){
                
                // Set the jet pT information for legends and figure saving
                if(iJetPtMinBias == histograms[kMinBiasHydjetSimulation]->GetNJetPtBinsEEC()){
                  jetPtStringMinBias = Form("Jet p_{T} > %.0f", histograms[kMinBiasHydjetSimulation]->GetCard()->GetJetPtCut());
                  compactJetPtStringMinBias = "";
                } else {
                  jetPtStringMinBias = Form("%.0f < jet p_{T} < %.0f", histograms[kMinBiasHydjetSimulation]->GetJetPtBinBorderEEC(iJetPtMinBias), histograms[kMinBiasHydjetSimulation]->GetJetPtBinBorderEEC(iJetPtMinBias+1));
                  compactJetPtStringMinBias = Form("_J=%.0f-%.0f", histograms[kMinBiasHydjetSimulation]->GetJetPtBinBorderEEC(iJetPtMinBias), histograms[kMinBiasHydjetSimulation]->GetJetPtBinBorderEEC(iJetPtMinBias+1));
                }
                
                // For logarithmic drawing, cannot go down to zero
                if(logDeltaR){
                  hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalizationType][minBiasRatioIndex[iRatio].second]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
                }
                
                hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalizationType][minBiasRatioIndex[iRatio].second]->SetLineColor(color[iJetPtMinBias]);
                hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalizationType][minBiasRatioIndex[iRatio].second]->Draw("same");
                
                ptLegend->AddEntry(hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalizationType][minBiasRatioIndex[iRatio].second], Form("MB: %s", jetPtStringMinBias.Data()), "l");
                
              } // Min bias jet pT loop
              
              // Draw the legend to the upper pad
              legend->Draw();
              ptLegend->Draw();
              
              // Linear scale for the ratio
              drawer->SetLogY(false);
              
              // Draw the ratios to the lower pad
              for(int iJetPtMinBias = firstStudiedJetPtBinEEC[kMinBiasHydjetSimulation]; iJetPtMinBias <= lastStudiedJetPtBinEEC[kMinBiasHydjetSimulation]; iJetPtMinBias++){
                
                // For logarithmic drawing, cannot go down to zero
                if(logDeltaR){
                  hMinimumBiasRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalizationType][iRatio]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
                }
                
                hMinimumBiasRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalizationType][iRatio]->SetLineColor(color[iJetPtMinBias]);
                
                if(addSignalToTotalRatioMinBias[iRatio]){
                  
                  drawer->SetGridY(true);
                  if(logDeltaR){
                    hSignalToTotalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
                  }
                  
                  hSignalToTotalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[0]);
                  hSignalToTotalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoomMinBias[iRatio].first,ratioZoomMinBias[iRatio].second);
                  drawer->DrawHistogramToLowerPad(hSignalToTotalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt], "#Deltar", minBiasRatioText[iRatio], " ");
                  drawer->SetGridY(false);
                }
                
                if(iJetPtMinBias == firstStudiedJetPtBinEEC[kMinBiasHydjetSimulation] && !addSignalToTotalRatioMinBias[iRatio]){
                  drawer->SetGridY(true);
                  hMinimumBiasRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalizationType][iRatio]->GetYaxis()->SetRangeUser(ratioZoomMinBias[iRatio].first,ratioZoomMinBias[iRatio].second);
                  drawer->DrawHistogramToLowerPad(hMinimumBiasRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalizationType][iRatio], "#Deltar", minBiasRatioText[iRatio], " ");
                  drawer->SetGridY(false);
                } else {
                  hMinimumBiasRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iNormalizationType][iRatio]->Draw("same");
                }
                
              } // Min bias jet pT loop
              
              // Save the figures to a file
              if(saveFigures){
                gPad->GetCanvas()->SaveAs(Form("figures/%sComparison%s_%s%s%s%s.%s", saveNameMinBias[iRatio], saveComment, histograms[kPythiaHydjetSimulation]->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
              }
              
            } // Track pT loop
          } // Jet pT loop
        } // Centrality loop
      } // Energy-energy correlator loop
      
    } // Normalization type loop [Gen level/reflected cone]
  } // Ratio type loop
  
}
