#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Draw comparison plots from
 */

/*
 * Macro for studying the reflected cone background estimation using MC
 */
void studyReflectedConeBackground(){

  // File from which the integrals are calculated
  TString inputFileName = "data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJet_fakeFakeReflectedCone_noTrigger_preprocessed_2022-09-23.root";
  // data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJet_fakeFakeReflectedCone_noTrigger_preprocessed_2022-09-23.root
  // data/PbPbMC2018_GenGen_eecAnalysis_genJet_fakeFakeReflectedCone_noTrigger_preprocessed_2022-09-30.root
  
  // Open the input file
  TFile *inputFile = TFile::Open(inputFileName);
  
  if(inputFile == NULL){
    cout << "Error! The file " << inputFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }
  
  // Load the card from the file and read the collision system
  EECCard *card = new EECCard(inputFile);
  TString collisionSystem = card->GetDataType();
  
  // ====================================================
  //    Binning configuration for the background study
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = card->GetNCentralityBins();
  const int nJetPtBinsEEC = card->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card->GetNTrackPtBinsEEC();
  
  // Default binning ranges for reference
  // centrality = {0,10,30,50,90}
  // track pT =  {0.7,1,2,3,4,8,12,300}
  // jet pT in energy-energy correlators = {120,140,160,180,200,300,500,5030}
  // track pT in energy-energy correlators = {0.7,1,2,3,4,6,300}
  
  // Bin range to be integrated
  int firstStudiedCentralityBin = 0;
  int lastStudiedCentralityBin = 0;
  
  int firstStudiedJetPtBinEEC = 0;
  int lastStudiedJetPtBinEEC = nJetPtBinsEEC; // Note: Jets integrated over all pT ranges are in nJetPtBinsEEC bin
  
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
  EECHistogramManager *histograms = new EECHistogramManager(inputFile,card);
  
  // Choose the energy-energy correlator types to load
  histograms->SetLoadEnergyEnergyCorrelators(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator]);
  histograms->SetLoadEnergyEnergyCorrelatorsJetPt(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPt]);
  histograms->SetLoadEnergyEnergyCorrelatorsUncorrected(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorUncorrected]);
  histograms->SetLoadEnergyEnergyCorrelatorsJetPtUncorrected(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPtUncorrected]);
  
  // Choose the bin ranges
  histograms->SetCentralityBinRange(firstStudiedCentralityBin,lastStudiedCentralityBin);
  histograms->SetJetPtBinRangeEEC(firstStudiedJetPtBinEEC,lastStudiedJetPtBinEEC);
  histograms->SetTrackPtBinRangeEEC(firstStudiedTrackPtBinEEC,lastStudiedTrackPtBinEEC);
  
  // Load the histograms from the file
  histograms->LoadProcessedHistograms();
  
  // Number of bins from the large deltaR used for the normalization of reflected cone background
  const int nBackgroundNormalizationBins = 9;
  
  // Indices for different types of energy-energy correlators
  enum enumCorrelatorPairingTypes{kTotalEEC, kSignalEEC, kSignalFakeEEC, kFakeFakeEEC, kBackgroundEEC, knPairingTypesEEC};
  enum enumReflectedConeType{kPairSignalReflectedCone, kPairOnlyReflectedCone, kCombineReflectedCone, kSubtractSignalReflectedCone, kSubtractCombined, knReflectedConeTypes};
  
  // Select which ratios should be made from the energy-energy correlators and reflected cone histograms and their configuration
  const int nRatioTypes = 5;
  std::pair<int,int> ratioIndex[nRatioTypes];
  bool drawComparisonType[nRatioTypes];
  const char* legendTextEnergyEnergyCorrelator[nRatioTypes];
  const char* legendTextReflectedCone[nRatioTypes];
  const char* ratioText[nRatioTypes];
  std::pair<double,double> yRange[nRatioTypes];
  std::pair<double,double> ratioZoom[nRatioTypes];
  const char* saveName[nRatioTypes];
  
  // Index 0: Compare signal+fake to corresponding reflected cone distribution
  drawComparisonType[0] = false;
  ratioIndex[0] = std::make_pair(kSignalFakeEEC, kPairSignalReflectedCone);
  legendTextEnergyEnergyCorrelator[0] = "Signal+fake pairs";
  legendTextReflectedCone[0] = "Jet+ref";
  ratioText[0] = "(Jet+ref) / BG";
  yRange[0] = std::make_pair(0.0005, 2);
  ratioZoom[0] = std::make_pair(0, 2);
  saveName[0] = "pythiaHydjetToSignalReflectedCone";
  
  // Index 1: Compare fake+fake to corresponding reflected cone distribution
  drawComparisonType[1] = false;
  ratioIndex[1] = std::make_pair(kFakeFakeEEC, kPairOnlyReflectedCone);
  legendTextEnergyEnergyCorrelator[1] = "Fake+fake pairs";
  legendTextReflectedCone[1] = "Ref+ref";
  ratioText[1] = "(Ref+ref) / BG";
  yRange[1] = std::make_pair(0.0005, 2);
  ratioZoom[1] = std::make_pair(0, 2);
  saveName[1] = "hydjetHydjetToOnlyReflectedCone";
  
  // Index 2: Compare total background to signal+fake reflected cone
  drawComparisonType[2] = false;
  ratioIndex[2] = std::make_pair(kBackgroundEEC, kPairSignalReflectedCone);
  legendTextEnergyEnergyCorrelator[2] = "All background pairs";
  legendTextReflectedCone[2] = "Jet+ref";
  ratioText[2] = "(Jet+ref) / BG";
  yRange[2] = std::make_pair(0.0005, 2);
  ratioZoom[2] = std::make_pair(0, 2);
  saveName[2] = "allBackgroundToSignalReflectedCone";
  
  // Index 3: Compare signal to signal+fake reflected cone subtracted total distribution
  drawComparisonType[3] = true;
  ratioIndex[3] = std::make_pair(kSignalEEC, kSubtractSignalReflectedCone);
  legendTextEnergyEnergyCorrelator[3] = "Signal pairs";
  legendTextReflectedCone[3] = "Jet+ref";
  ratioText[3] = "BGsub / signal";
  yRange[3] = std::make_pair(0.001, 50);
  ratioZoom[3] = std::make_pair(0, 2);
  saveName[3] = "signalToSignalReflectedConeSubtracted";
  
  // Index 4: Compare signal to combined reflected cone subtracted total distribution
  drawComparisonType[4] = false;
  ratioIndex[4] = std::make_pair(kSignalEEC, kSubtractCombined);
  legendTextEnergyEnergyCorrelator[4] = "Signal pairs";
  legendTextReflectedCone[4] = "Reflected";
  ratioText[4] = "BGsub / signal";
  yRange[4] = std::make_pair(0.001, 50);
  ratioZoom[4] = std::make_pair(0, 2);
  saveName[4] = "signalToReflectedConeSubtracted";
  
  // Energy-energy correlator histograms separated by subevents from the Pythia+Hydjet simulation
  TH1D* hEnergyEnergyCorrelator[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC][knPairingTypesEEC];
  
  // Reflected cone energy-energy correlators
  TH1D* hReflectedCone[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC][nBackgroundNormalizationBins][knReflectedConeTypes];
  
  // Histograms for all different ratios
  TH1D* hRatio[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC][nBackgroundNormalizationBins][nRatioTypes];
  
  // Initialize the energy-energy correlator histogram array to NULL
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC+1; iJetPt++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
          for(int iPairingType = 0; iPairingType < knPairingTypesEEC; iPairingType++){
            hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iPairingType] = NULL;
          } // Pairing type loop
          for(int iNormalization = 0; iNormalization < nBackgroundNormalizationBins; iNormalization++){
            for(int iReflectedConeType = 0; iReflectedConeType < knReflectedConeTypes; iReflectedConeType++){
              hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iReflectedConeType] = NULL;
            } // Reflected cone loop
            for(int iRatio = 0; iRatio < nRatioTypes; iRatio++){
              hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iRatio] = NULL;
            } // Ratio loop
          } // Normalization region loop
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Energy-energy correlator loop

  // Helper variables
  TH1D *helperHistogram, *helperHistogram2;
  double normalizationFactor;
  int lowIntegralBin, highIntegralBin;
  int studiedEnergyEnergyCorrelatorIndex = -1;
  
  // Get the histograms from the histogram manager
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only read the selected energy-energy correlator types
    if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
    studiedEnergyEnergyCorrelatorIndex = iEnergyEnergyCorrelator;
    
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
          
          // Histogram with all pair combinations. Normalize it to one
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC] = histograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJetPair, EECHistogramManager::knSubeventTypes);
          normalizationFactor = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC]->Integral("width");
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC]->Scale(1/normalizationFactor);
          
          // Histogram with only signal. Normalize the sum of signal and background to total.
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kSignalEEC] = histograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJetPair, EECHistogramManager::kPythiaPythia);
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kSignalEEC]->Scale(1/normalizationFactor);
          
          // Histograms with background contributions. Normalize to the total number of pairs.
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kSignalFakeEEC] = histograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJetPair, EECHistogramManager::kPythiaHydjet);
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kBackgroundEEC] = (TH1D *) hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kSignalFakeEEC]->Clone(Form("energyEnergyCorrelatorBackground%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt));
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kFakeFakeEEC] = histograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJetPair, EECHistogramManager::kHydjetHydjet);
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kBackgroundEEC]->Add(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kFakeFakeEEC]);
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kSignalFakeEEC]->Scale(1/normalizationFactor);
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kFakeFakeEEC]->Scale(1/normalizationFactor);
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kBackgroundEEC]->Scale(1/normalizationFactor);
          
          // The jet cone + reflected cone histogram.
          helperHistogram = histograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSignalReflectedConePair, EECHistogramManager::knSubeventTypes);
          
          // The only reflected cone histogram.
          helperHistogram2 = histograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kReflectedConePair, EECHistogramManager::knSubeventTypes);
          
          // Normalize the reflected cone histogram to different DeltaR regions of the tail of the total distribution
          highIntegralBin = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC]->GetNbinsX();
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
              normalizationFactor = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTotalEEC]->Integral(lowIntegralBin, highIntegralBin, "width") / hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iReflectedConeType]->Integral(lowIntegralBin, highIntegralBin, "width");
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
            
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Energy-energy correlator loop
  
  // Find the binning information for DeltaR bins from a histogram
  const int nDeltaRBins = hEnergyEnergyCorrelator[studiedEnergyEnergyCorrelatorIndex][firstStudiedCentralityBin][firstStudiedJetPtBinEEC][firstStudiedTrackPtBinEEC][kTotalEEC]->GetNbinsX();
  double deltaRBinBorders[nDeltaRBins+1];
  deltaRBinBorders[0] = hEnergyEnergyCorrelator[studiedEnergyEnergyCorrelatorIndex][firstStudiedCentralityBin][firstStudiedJetPtBinEEC][firstStudiedTrackPtBinEEC][kTotalEEC]->GetXaxis()->GetBinLowEdge(1);
  for(int iBin = 1; iBin <= nDeltaRBins; iBin++){
    deltaRBinBorders[iBin] = hEnergyEnergyCorrelator[studiedEnergyEnergyCorrelatorIndex][firstStudiedCentralityBin][firstStudiedJetPtBinEEC][firstStudiedTrackPtBinEEC][kTotalEEC]->GetXaxis()->GetBinUpEdge(iBin);
  }
  
  
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  
  if(logDeltaR) drawer->SetLogX(true);
  
  TLegend *legend;
  TString centralityString, trackPtString, jetPtString;
  TString compactCentralityString, compactTrackPtString, compactJetPtString;
  int color[9] = {kBlue, kRed, kGreen+2, kCyan, kMagenta, kOrange-1, kAzure+5, kOrange-2, kGray};
  double maximumY, currentY;
  
  // First set of plots: compare background histograms to reflected cone histograms with different normalization region
  for(int iRatio = 0; iRatio < nRatioTypes; iRatio++){
    
    // Only draw the selected ratio types
    if(!drawComparisonType[iRatio]) continue;
    
    // Loop over all selevted histograms
    for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
      
      // Only read the selected energy-energy correlator types
      if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
      
      for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
        
        // Set the centrality information for legends and figure saving
        centralityString = Form("Cent: %.0f-%.0f%%",histograms->GetCentralityBinBorder(iCentrality),histograms->GetCentralityBinBorder(iCentrality+1));
        compactCentralityString = Form("_C=%.0f-%.0f",histograms->GetCentralityBinBorder(iCentrality),histograms->GetCentralityBinBorder(iCentrality+1));
        
        for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
          
          // Set the jet pT information for legends and figure saving
          if(iJetPt == histograms->GetNJetPtBinsEEC()){
            jetPtString = Form("Jet p_{T} > %.0f", histograms->GetCard()->GetJetPtCut());
            compactJetPtString = "";
          } else {
            jetPtString = Form("%.0f < jet p_{T} < %.0f", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1));
            compactJetPtString = Form("_J=%.0f-%.0f", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1));
          }
          
          for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
            
            // Set the track pT information for legends and figure saving
            trackPtString = Form("%.1f < track p_{T}",histograms->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString = Form("_T%.1f",histograms->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString.ReplaceAll(".","v");
            
            // Create a legend for the figure
            legend = new TLegend(0.19,0.09,0.39,0.94);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->AddEntry((TObject*) 0, histograms->GetCard()->GetAlternativeDataType().Data(), "");
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
            drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][ratioIndex[iRatio].first], "#Deltar", histograms->GetEnergyEnergyCorrelatorAxisName(iEnergyEnergyCorrelator), " ");
            legend->AddEntry(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][ratioIndex[iRatio].first], legendTextEnergyEnergyCorrelator[iRatio], "l");
            
            // Draw different normalization regions for the reflected cone histogram to the same plot
            for(int iNormalization = 0; iNormalization < nBackgroundNormalizationBins; iNormalization++){
              
              // For logarithmic drawing, cannot go down to zero
              if(logDeltaR){
                hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][ratioIndex[iRatio].second]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
              }
              
              hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][ratioIndex[iRatio].second]->SetLineColor(color[iNormalization]);
              hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][ratioIndex[iRatio].second]->Draw("same");
              legend->AddEntry(hReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][ratioIndex[iRatio].second], Form("%s, #Deltar > %.3f", legendTextReflectedCone[iRatio], deltaRBinBorders[nDeltaRBins-(nBackgroundNormalizationBins-iNormalization)]), "l");
              
            } // Normalization region loop
            
            // Draw the legend to the upper pag
            legend->Draw();
            
            // Linear scale for the ratio
            drawer->SetLogY(false);
            
            // Draw the ratios to the lower pad
            for(int iNormalization = 0; iNormalization < nBackgroundNormalizationBins; iNormalization++){
              
              // For logarithmic drawing, cannot go down to zero
              if(logDeltaR){
                hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iRatio]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
              }
              
              hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iRatio]->SetLineColor(color[iNormalization]);
              if(iNormalization == 0){
                drawer->SetGridY(true);
                hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iRatio]->GetYaxis()->SetRangeUser(ratioZoom[iRatio].first, ratioZoom[iRatio].second);
                drawer->DrawHistogramToLowerPad(hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iRatio], "#Deltar", ratioText[iRatio], " ");
                drawer->SetGridY(false);
              } else {
                hRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization][iRatio]->Draw("same");
              }
              
            } // Normalization region loop
            
            // Save the figures to a file
            if(saveFigures){
              gPad->GetCanvas()->SaveAs(Form("figures/%sComparison%s_%s%s%s%s.%s", saveName[iRatio], saveComment, histograms->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
            }
            
          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Energy-energy correlator loop
    
  } // Ratio type loop
  
}
