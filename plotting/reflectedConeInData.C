#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Macro for studying the reflected cone background estimation using MC
 */
void reflectedConeInData(){
  
  enum enumDataType{kSignalEvent, kMinimumBias, knDataTypes};
  
  // Input files: index 0 = Pythia+Hydjet simulation, index 1 = minimum bias Hydjet simulation
  TString inputFileName[knDataTypes] = {"data/eecAnalysis_akFlowJets_fakeFakeReflectedCone_preprocessed_2022-09-23.root", "data/MinBiasHydjet_RecoGen_eecAnalysis_akFlowJet_firstMinBiasScan_noTrigger_preprocessed_2022-10-10.root"};
  // data/eecAnalysis_akFlowJets_fakeFakeReflectedCone_preprocessed_2022-09-23.root
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
  
  // Require matching centrality and track pT bins between signal event and minimum bias event files
  const double epsilon = 0.00001;
  
  const int nCentralityBins = card[kSignalEvent]->GetNCentralityBins();
  if(nCentralityBins != card[kMinimumBias]->GetNCentralityBins()){
    cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
    return;
  }
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    if(TMath::Abs(card[kSignalEvent]->GetLowBinBorderCentrality(iCentrality) - card[kMinimumBias]->GetLowBinBorderCentrality(iCentrality)) > epsilon){
      cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
    if(TMath::Abs(card[kSignalEvent]->GetHighBinBorderCentrality(iCentrality) - card[kMinimumBias]->GetHighBinBorderCentrality(iCentrality)) > epsilon){
      cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
  }
  
  const int nTrackPtBinsEEC = card[kSignalEvent]->GetNTrackPtBinsEEC();
  if(nTrackPtBinsEEC != card[kMinimumBias]->GetNTrackPtBinsEEC()){
    cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
    return;
  }
  for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
    if(TMath::Abs(card[kSignalEvent]->GetLowBinBorderTrackPtEEC(iTrackPt) - card[kMinimumBias]->GetLowBinBorderTrackPtEEC(iTrackPt)) > epsilon){
      cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
    if(TMath::Abs(card[kSignalEvent]->GetHighBinBorderTrackPtEEC(iTrackPt) - card[kMinimumBias]->GetHighBinBorderTrackPtEEC(iTrackPt)) > epsilon){
      cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
  }
  
  // ====================================================
  //    Binning configuration for the background study
  // ====================================================
  
  // Find the number of bins from the card
  const int nJetPtBinsEEC[knDataTypes] = {card[kSignalEvent]->GetNJetPtBinsEEC(), card[kMinimumBias]->GetNJetPtBinsEEC()};
  
  // Default binning ranges for reference
  // centrality = {0,10,30,50,90}
  // track pT =  {0.7,1,2,3,4,8,12,300}
  // jet pT in energy-energy correlators = {120,140,160,180,200,300,500,5030}
  // track pT in energy-energy correlators = {0.7,1,1.5,2,2.5,3,3.5,4,300}
  
  // Bin range to be integrated
  int firstStudiedCentralityBin = 0;
  int lastStudiedCentralityBin = 0;
  
  int firstStudiedJetPtBinEEC[knDataTypes] = {0,0};
  int lastStudiedJetPtBinEEC[knDataTypes] = {nJetPtBinsEEC[kSignalEvent], nJetPtBinsEEC[kMinimumBias]}; // Note: Jets integrated over all pT ranges are in nJetPtBinsEEC bin
  
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
  //enum enumCorrelatorPairingTypes{kTotalEEC, kSignalEEC, kSignalFakeEEC, kFakeFakeEEC, kBackgroundEEC, knPairingTypesEEC};
  enum enumReflectedConeType{kPairSignalReflectedCone, kPairOnlyReflectedCone, kCombineReflectedCone, kSubtractSignalReflectedCone, kSubtractCombined, knReflectedConeTypes};
  enum enumMinBiasHistogramStyle{kPureMinBias, kReflectedConeCombinedMinBias, kMinBiasBackgroundSubtracted, knMinBiasHistograms};
  enum enumNormalizationStyle{kNormalizeToReflectedCone, kNormalizeToTotal, knNormalizationStyles};
  
  TString matchString[knNormalizationStyles] = {"Ref cone match", "Gen level match"};
  
  // =========================================================================================================================== //
  // Select which ratios should be made from the energy-energy correlators and reflected cone histograms and their configuration //
  // =========================================================================================================================== //
  const int nRatioTypes = 1;
  bool addTotalDistribution[nRatioTypes];
  int ratioIndex[nRatioTypes];
  bool drawComparisonType[nRatioTypes];
  const char* legendTextReflectedCone[nRatioTypes];
  std::pair<double,double> yRange[nRatioTypes];
  std::pair<double,double> ratioZoom[nRatioTypes];
  const char* saveName[nRatioTypes];
  
  // Index 0: Compare signal+fake to corresponding reflected cone distribution
  drawComparisonType[0] = true;
  addTotalDistribution[0] = true;
  ratioIndex[0] = kPairSignalReflectedCone;
  legendTextReflectedCone[0] = "Jet+ref";
  yRange[0] = std::make_pair(0.0005, 2);
  ratioZoom[0] = std::make_pair(0, 2);
  saveName[0] = "reflectedConeNormalization";
  
  
  const int nMinBiasRatioTypes = 1;
  
  // Energy-energy correlator histograms separated by subevents from the Pythia+Hydjet simulation
  TH1D* hEnergyEnergyCorrelator[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC[kSignalEvent]+1][nTrackPtBinsEEC];
  
  // Reflected cone energy-energy correlators
  TH1D* hReflectedCone[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC[kSignalEvent]+1][nTrackPtBinsEEC][nBackgroundNormalizationBins][knReflectedConeTypes];
  
  // Histograms for all different ratios
  TH1D* hRatio[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC[kSignalEvent]+1][nTrackPtBinsEEC][nBackgroundNormalizationBins][nRatioTypes][2];
  
  // Energy-energy correlators from minimum bias Hydjet simulation
  TH1D* hMinimumBias[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC[kSignalEvent]+1][nJetPtBinsEEC[kMinimumBias]+1][nTrackPtBinsEEC][knMinBiasHistograms];
  TH1D* hMinimumBiasRatio[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC[kSignalEvent]+1][nJetPtBinsEEC[kMinimumBias]+1][nTrackPtBinsEEC][nMinBiasRatioTypes];
  
  // Initialize the energy-energy correlator histogram array to NULL
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC[kSignalEvent]+1; iJetPt++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = NULL;
          for(int iNormalization = 0; iNormalization < nBackgroundNormalizationBins; iNormalization++){
            for(int iReflectedConeType = 0; iReflectedConeType < knReflectedConeTypes; iReflectedConeType++){
              hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iReflectedConeType] = NULL;
            } // Reflected cone loop
            for(int iRatio = 0; iRatio < nRatioTypes; iRatio++){
              for(int iNormalizationType = 0; iNormalizationType < knNormalizationStyles; iNormalizationType++)
              hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iRatio][iNormalizationType] = NULL;
            } // Ratio loop
          } // Normalization region loop
          for(int iJetPtMinBias = 0; iJetPtMinBias < nJetPtBinsEEC[kMinimumBias]+1; iJetPtMinBias++){
            for(int iMinBias = 0; iMinBias < knMinBiasHistograms; iMinBias++){
              hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iMinBias] = NULL;
            }
            for(int iRatio = 0; iRatio < nMinBiasRatioTypes; iRatio++){
              hMinimumBiasRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iRatio] = NULL;
            } // Ratio loop
            
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
  double signalFaketoFakeFakeRatio;
  
  // Get the histograms from the histogram manager
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only read the selected energy-energy correlator types
    if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
    studiedEnergyEnergyCorrelatorIndex = iEnergyEnergyCorrelator;
    
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iJetPt = firstStudiedJetPtBinEEC[kSignalEvent]; iJetPt <= lastStudiedJetPtBinEEC[kSignalEvent]; iJetPt++){
        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
          
          // ================================================== //
          // ===   Read the histograms from signal events   === //
          // ================================================== //
          
          // Histogram with all pair combinations. Normalize it to one
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = histograms[kSignalEvent]->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJetPair);
          normalizationFactor = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Integral("width");
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Scale(1/normalizationFactor);
          
          // The jet cone + reflected cone histogram.
          helperHistogram = histograms[kSignalEvent]->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSignalReflectedConePair);
          
          // The only reflected cone histogram.
          helperHistogram2 = histograms[kSignalEvent]->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kReflectedConePair);
          
          // Calculate the signal+fake to fake+fake ratio for the reflected cone histograms
          signalFaketoFakeFakeRatio = helperHistogram->Integral("width") / helperHistogram2->Integral("width");
          
          // Normalize the reflected cone histogram to different DeltaR regions of the tail of the total distribution
          highIntegralBin = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetNbinsX();
          for(int iNormalization = 0; iNormalization < nBackgroundNormalizationBins; iNormalization++){
            lowIntegralBin = highIntegralBin - (nBackgroundNormalizationBins - iNormalization - 1);

            // Normalization for jet cone + reflected cone histogram
            hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][kPairSignalReflectedCone] = (TH1D*) helperHistogram->Clone(Form("reflectedCone%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iNormalization));
            
            // Normalization for only reflected cone histogram
            hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][kPairOnlyReflectedCone] = (TH1D*) helperHistogram2->Clone(Form("reflectedConeOnly%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iNormalization));
            
            // Normalization for total reflected cone histogram, where all the pairs in went through once in each event
            hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][kCombineReflectedCone] = (TH1D*) helperHistogram->Clone(Form("reflectedConeCombined%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iNormalization));
            hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][kCombineReflectedCone]->Add(helperHistogram2);
            
            // Normalization for all reflected cone histograms
            for(int iReflectedConeType = kPairSignalReflectedCone; iReflectedConeType <= kCombineReflectedCone; iReflectedConeType++){
              normalizationFactor = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Integral(lowIntegralBin, highIntegralBin, "width") / hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iReflectedConeType]->Integral(lowIntegralBin, highIntegralBin, "width");
              hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iReflectedConeType]->Scale(normalizationFactor);
            }
            
            // The background subtracted histograms can be calculated from the normalized reflected cone histograms.
            
            // First clone the total distribution as the baseline
            for(int iReflectedConeType = kSubtractSignalReflectedCone; iReflectedConeType <= kSubtractCombined; iReflectedConeType++){
              hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iReflectedConeType] = (TH1D*) hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Clone(Form("backgroundSubtracted%d%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iNormalization, iReflectedConeType));
            }
            
            // Then subtract the correct reflected cone histogram from the total histogram
            hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][kSubtractSignalReflectedCone]->Add(hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][kPairSignalReflectedCone],-1);
            hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][kSubtractCombined]->Add(hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][kCombineReflectedCone],-1);
            
            // After all the reflected cone histograms are ready, calculate the ratios to the first normalization bin
            for(int iRatio = 0; iRatio < nRatioTypes; iRatio++){
              hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iRatio][kNormalizeToReflectedCone] = (TH1D*) hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][ratioIndex[iRatio]]->Clone(Form("ratioToReflectedCone%d%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iNormalization, iRatio));
              hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iRatio][kNormalizeToReflectedCone]->Divide(hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0][ratioIndex[iRatio]]);
              
              hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iRatio][kNormalizeToTotal] = (TH1D*) hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][ratioIndex[iRatio]]->Clone(Form("ratioToTotal%d%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iNormalization, iRatio));
              hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iRatio][kNormalizeToTotal]->Divide(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]);
            }
            
          } // Normalization region loop
          
          // ====================================================== //
          // ===   Read the histograms from minimum bias data   === //
          // ====================================================== //
          
          /*
           TODO: Implementation for minimum bias
          for(int iJetPtMinBias = firstStudiedJetPtBinEEC[kMinimumBias]; iJetPtMinBias <= lastStudiedJetPtBinEEC[kMinimumBias]; iJetPtMinBias++){
            hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][kPureMinBias] = (TH1D*) histograms[kMinimumBias]->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPtMinBias, iTrackPt, EECHistograms::kSameJetPair, EECHistogramManager::knSubeventTypes)->Clone(Form("pureMinBiasHistogram%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iJetPtMinBias, iTrackPt));
            
            // Normalize the minimum bias histograms to match the signal+fake to fake+fake ratio with different estimators
            normalizationFactor = hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0][kPairSignalReflectedCone]->Integral("width") / (hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][kPureMinBias]->Integral("width") * signalFaketoFakeFakeRatio);
            hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][kPureMinBias]->Scale(normalizationFactor);
            hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][kReflectedConeCombinedMinBias] = (TH1D*) hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0][kPairSignalReflectedCone]->Clone(Form("combinedBackground%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iJetPtMinBias, iTrackPt));
            hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][kReflectedConeCombinedMinBias]->Add(hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][kPureMinBias]);
            
            // After the ratio is matched, normalize the distribution such that it matches the generator level background at 0.3 < DeltaR < 0.4
            // TODO: Normalization for the combined background.
            
            // Now that we have a normalized background estimate, we can do a background subtraction
            hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][kMinBiasBackgroundSubtracted] = (TH1D*) hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Clone(Form("minBiasBackgroundSubtracted%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iJetPtMinBias, iTrackPt));
            hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][kMinBiasBackgroundSubtracted]->Add(hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][kReflectedConeCombinedMinBias], -1);
            
            // Calculate all defined ratios with minimum bias histograms
            // TODO: Ratios also here need to be redefined
            for(int iRatio = 0; iRatio < nMinBiasRatioTypes; iRatio++){
              hMinimumBiasRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iRatio] = (TH1D*) hMinimumBias[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][minBiasRatioIndex[iRatio].second]->Clone(Form("minimumBiasRatio%d%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iJetPtMinBias, iTrackPt, iRatio));
              hMinimumBiasRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iRatio]->Divide(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]);
            } // Ratio index loop
            
          } // Jet pT loop for minimum bias hydjet simulation
          */
           
        } // Track pT loop
      } // Jet pT loop for Pythia+Hydjet simulation
    } // Centrality loop
  } // Energy-energy correlator loop
  
  // Find the binning information for DeltaR bins from a histogram
  const int nDeltaRBins = hEnergyEnergyCorrelator[studiedEnergyEnergyCorrelatorIndex][firstStudiedCentralityBin][firstStudiedJetPtBinEEC[kSignalEvent]][firstStudiedTrackPtBinEEC]->GetNbinsX();
  double deltaRBinBorders[nDeltaRBins+1];
  deltaRBinBorders[0] = hEnergyEnergyCorrelator[studiedEnergyEnergyCorrelatorIndex][firstStudiedCentralityBin][firstStudiedJetPtBinEEC[kSignalEvent]][firstStudiedTrackPtBinEEC]->GetXaxis()->GetBinLowEdge(1);
  for(int iBin = 1; iBin <= nDeltaRBins; iBin++){
    deltaRBinBorders[iBin] = hEnergyEnergyCorrelator[studiedEnergyEnergyCorrelatorIndex][firstStudiedCentralityBin][firstStudiedJetPtBinEEC[kSignalEvent]][firstStudiedTrackPtBinEEC]->GetXaxis()->GetBinUpEdge(iBin);
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
  int ratioType;
  TString ratioText;
  
  // Legend y-positions
  double legendY1 = 0.09;
  double legendY2 = 0.94;
  
  // ==========================================================
  // ===   Drawing Pythia-Hydjet to reflected cone ratios   ===
  // ==========================================================
  
  // Loop over all different ratio combinations
  for(int iRatio = 0; iRatio < nRatioTypes; iRatio++){
    
    // Only draw the selected ratio types
    if(!drawComparisonType[iRatio]) continue;
    
    // Make a very nice ratio text
    ratioText = Form("%s/#Deltar > %.3f", legendTextReflectedCone[iRatio], deltaRBinBorders[nDeltaRBins-nBackgroundNormalizationBins]);
    if(addTotalDistribution[iRatio]) ratioText = Form("%s/Total", legendTextReflectedCone[iRatio]);
    
    // Loop over all selected histograms
    for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
      
      // Only read the selected energy-energy correlator types
      if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
      
      for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
        
        // Set the centrality information for legends and figure saving
        centralityString = Form("Cent: %.0f-%.0f%%",histograms[kSignalEvent]->GetCentralityBinBorder(iCentrality),histograms[kSignalEvent]->GetCentralityBinBorder(iCentrality+1));
        compactCentralityString = Form("_C=%.0f-%.0f",histograms[kSignalEvent]->GetCentralityBinBorder(iCentrality),histograms[kSignalEvent]->GetCentralityBinBorder(iCentrality+1));
        
        for(int iJetPt = firstStudiedJetPtBinEEC[kSignalEvent]; iJetPt <= lastStudiedJetPtBinEEC[kSignalEvent]; iJetPt++){
          
          // Set the jet pT information for legends and figure saving
          if(iJetPt == histograms[kSignalEvent]->GetNJetPtBinsEEC()){
            jetPtString = Form("Jet p_{T} > %.0f", histograms[kSignalEvent]->GetCard()->GetJetPtCut());
            compactJetPtString = "";
          } else {
            jetPtString = Form("%.0f < jet p_{T} < %.0f", histograms[kSignalEvent]->GetJetPtBinBorderEEC(iJetPt), histograms[kSignalEvent]->GetJetPtBinBorderEEC(iJetPt+1));
            compactJetPtString = Form("_J=%.0f-%.0f", histograms[kSignalEvent]->GetJetPtBinBorderEEC(iJetPt), histograms[kSignalEvent]->GetJetPtBinBorderEEC(iJetPt+1));
          }
          
          for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
            
            // Set the track pT information for legends and figure saving
            trackPtString = Form("%.1f < track p_{T}",histograms[kSignalEvent]->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString = Form("_T%.1f",histograms[kSignalEvent]->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString.ReplaceAll(".","v");
            
            // Create a legend for the figure
            legend = new TLegend(0.19,legendY1,0.39,legendY2);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->AddEntry((TObject*) 0, histograms[kSignalEvent]->GetCard()->GetAlternativeDataType().Data(), "");
            legend->AddEntry((TObject*) 0, centralityString.Data(),"");
            legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
            legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
            
            // Create a new canvas for the plot
            drawer->CreateSplitCanvas();
            
            // Logarithmic EEC axis
            if(logEEC) drawer->SetLogY(true);
            
            // Draw the total distribution to canvas, if we select to do that:
            if(addTotalDistribution[iRatio]){
              
              // For logarithmic drawing, cannot go down to zero
              if(logDeltaR){
                hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
              }
              
              // Set good y-ranges for plotting
              hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(yRange[iRatio].first, yRange[iRatio].second);
              
              // Draw the background histogram to the upper canvas
              hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetLineColor(kBlack);
              drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt], "#Deltar", histograms[kSignalEvent]->GetEnergyEnergyCorrelatorAxisName(iEnergyEnergyCorrelator), " ");
              legend->AddEntry(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt], "All combinations", "l");
            } // Drawing total distribution to the canvas
            
            // Draw different normalization regions for the reflected cone histogram to the same plot
            for(int iNormalization = 0; iNormalization < nBackgroundNormalizationBins; iNormalization++){
              
              // For logarithmic drawing, cannot go down to zero
              if(logDeltaR){
                hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][ratioIndex[iRatio]]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
              }
              
              hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][ratioIndex[iRatio]]->SetLineColor(color[iNormalization]);
              
              // If the canvas has not been drawn before, do it in the first index
              if(!addTotalDistribution[iRatio] && iNormalization == 0){
                drawer->DrawHistogramToUpperPad(hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][ratioIndex[iRatio]], "#Deltar", histograms[kSignalEvent]->GetEnergyEnergyCorrelatorAxisName(iEnergyEnergyCorrelator), " ");
              } else {
                hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][ratioIndex[iRatio]]->Draw("same");
              }
              
              legend->AddEntry(hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][ratioIndex[iRatio]], Form("%s, #Deltar > %.3f", legendTextReflectedCone[iRatio], deltaRBinBorders[nDeltaRBins-(nBackgroundNormalizationBins-iNormalization)]), "l");
              
              
            } // Normalization region loop
            
            // Draw the legend to the upper pad
            legend->Draw();
            
            // Linear scale for the ratio
            drawer->SetLogY(false);
            
            // Choose the correct ratio to be drawn
            ratioType = kNormalizeToReflectedCone;
            if(addTotalDistribution[iRatio]) ratioType = kNormalizeToTotal;
            
            // Draw the ratios to the lower pad
            for(int iNormalization = 0; iNormalization < nBackgroundNormalizationBins; iNormalization++){
              
              // For logarithmic drawing, cannot go down to zero
              if(logDeltaR){
                hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iRatio][ratioType]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
              }
              
              hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iRatio][ratioType]->SetLineColor(color[iNormalization]);
              if(iNormalization == 0){
                drawer->SetGridY(true);
                hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iRatio][ratioType]->GetYaxis()->SetRangeUser(ratioZoom[iRatio].first, ratioZoom[iRatio].second);
                drawer->DrawHistogramToLowerPad(hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iRatio][ratioType], "#Deltar", ratioText.Data(), " ");
                drawer->SetGridY(false);
              } else {
                hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iRatio][ratioType]->Draw("same");
              }
              
            } // Normalization region loop
            
            // Save the figures to a file
            if(saveFigures){
              gPad->GetCanvas()->SaveAs(Form("figures/%sComparison%s_%s%s%s%s.%s", saveName[iRatio], saveComment, histograms[kSignalEvent]->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
            }
            
          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Energy-energy correlator loop
    
  } // Ratio type loop
 
  // ==========================================================================
  // ===   Drawing histograms including minimum bias Hydjet distributions   ===
  // ==========================================================================
  
  // TODO: Implementation
  /*
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
          centralityString = Form("Cent: %.0f-%.0f%%",histograms[kSignalEvent]->GetCentralityBinBorder(iCentrality),histograms[kSignalEvent]->GetCentralityBinBorder(iCentrality+1));
          compactCentralityString = Form("_C=%.0f-%.0f",histograms[kSignalEvent]->GetCentralityBinBorder(iCentrality),histograms[kSignalEvent]->GetCentralityBinBorder(iCentrality+1));
          
          for(int iJetPt = firstStudiedJetPtBinEEC[kSignalEvent]; iJetPt <= lastStudiedJetPtBinEEC[kSignalEvent]; iJetPt++){
            
            // Set the jet pT information for legends and figure saving
            if(iJetPt == histograms[kSignalEvent]->GetNJetPtBinsEEC()){
              jetPtString = Form("Jet p_{T} > %.0f", histograms[kSignalEvent]->GetCard()->GetJetPtCut());
              compactJetPtString = "";
            } else {
              jetPtString = Form("%.0f < jet p_{T} < %.0f", histograms[kSignalEvent]->GetJetPtBinBorderEEC(iJetPt), histograms[kSignalEvent]->GetJetPtBinBorderEEC(iJetPt+1));
              compactJetPtString = Form("_J=%.0f-%.0f", histograms[kSignalEvent]->GetJetPtBinBorderEEC(iJetPt), histograms[kSignalEvent]->GetJetPtBinBorderEEC(iJetPt+1));
            }
            
            for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
              
              // Set the track pT information for legends and figure saving
              trackPtString = Form("%.1f < track p_{T}",histograms[kSignalEvent]->GetTrackPtBinBorderEEC(iTrackPt));
              compactTrackPtString = Form("_T%.1f",histograms[kSignalEvent]->GetTrackPtBinBorderEEC(iTrackPt));
              compactTrackPtString.ReplaceAll(".","v");
              
              // Create a legend for the figure
              legendY1 = 0.04;
              legendY2 = 0.28 + addMatchString[iRatio]*0.06;
              legend = new TLegend(legendX1MinBias[iRatio],legendY1,legendX1MinBias[iRatio]+0.2,legendY2);
              legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
              legend->AddEntry((TObject*) 0, histograms[kSignalEvent]->GetCard()->GetAlternativeDataType().Data(), "");
              legend->AddEntry((TObject*) 0, centralityString.Data(),"");
              legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
              legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
              if(addMatchString[iRatio]) legend->AddEntry((TObject*) 0, matchString[iNormalizationType].Data(),"");
              
              ptLegend = new TLegend(legendX1MinBias[iRatio]+0.33, legendY1, legendX1MinBias[iRatio]+0.53, legendY1+0.06*(lastStudiedJetPtBinEEC[kMinimumBias] + 2 - firstStudiedJetPtBinEEC[kMinimumBias]));
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
              drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][minBiasRatioIndex[iRatio].first], "#Deltar", histograms[kSignalEvent]->GetEnergyEnergyCorrelatorAxisName(iEnergyEnergyCorrelator), " ");
              ptLegend->AddEntry(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][minBiasRatioIndex[iRatio].first], minBiasLegendTextEnergyEnergyCorrelator[iRatio], "l");
              
              // Option to draw the total distribution to the same figure
              if(addSignalToTotalRatioMinBias[iRatio]){
                hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC]->SetLineColor(color[0]);
                hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC]->Draw("same");
                ptLegend->AddEntry(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC], "All pairs", "l");
                
              }
              
              // Draw background subtracted histograms with fake+fake backgrounds estimated from different jet pT ranges
              for(int iJetPtMinBias = firstStudiedJetPtBinEEC[kMinimumBias]; iJetPtMinBias <= lastStudiedJetPtBinEEC[kMinimumBias]; iJetPtMinBias++){
                
                // Set the jet pT information for legends and figure saving
                if(iJetPtMinBias == histograms[kMinimumBias]->GetNJetPtBinsEEC()){
                  jetPtStringMinBias = Form("Jet p_{T} > %.0f", histograms[kMinimumBias]->GetCard()->GetJetPtCut());
                  compactJetPtStringMinBias = "";
                } else {
                  jetPtStringMinBias = Form("%.0f < jet p_{T} < %.0f", histograms[kMinimumBias]->GetJetPtBinBorderEEC(iJetPtMinBias), histograms[kMinimumBias]->GetJetPtBinBorderEEC(iJetPtMinBias+1));
                  compactJetPtStringMinBias = Form("_J=%.0f-%.0f", histograms[kMinimumBias]->GetJetPtBinBorderEEC(iJetPtMinBias), histograms[kMinimumBias]->GetJetPtBinBorderEEC(iJetPtMinBias+1));
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
              for(int iJetPtMinBias = firstStudiedJetPtBinEEC[kMinimumBias]; iJetPtMinBias <= lastStudiedJetPtBinEEC[kMinimumBias]; iJetPtMinBias++){
                
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
                
                if(iJetPtMinBias == firstStudiedJetPtBinEEC[kMinimumBias] && !addSignalToTotalRatioMinBias[iRatio]){
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
                gPad->GetCanvas()->SaveAs(Form("figures/%sComparison%s_%s%s%s%s.%s", saveNameMinBias[iRatio], saveComment, histograms[kSignalEvent]->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
              }
              
            } // Track pT loop
          } // Jet pT loop
        } // Centrality loop
      } // Energy-energy correlator loop
      
    } // Normalization type loop [Gen level/reflected cone]
  } // Ratio type loop
  */
}
