#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"

/*
 * Macro for finding background normalization scale in MC. This is done by integrating the energy-energy correlator distributions for background within the signal cone, and in the reflected cone, and taking a ratio of these numbers.
 */
void findBackgroundNormalizationScale(){

  // File from which the integrals are calculated
  TString inputFileName = "data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_4pC_newSRC_cutBadPhiAndComb_wtaAxis_jetTrigger_preprocessed_2022-11-29.root";
  // PbPbMC2018_GenGen_eecAnalysis_akFlowJet_MnD_wtaAxis_noTrigger_preprocessed_2022-10-21.root
  // PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_removeFakeJets_wtaAxis_noTrigger_preprocessed_2022-11-22.root
  // PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_4pC_newSRC_cutBadPhiAndComb_wtaAxis_jetTrigger_preprocessed_2022-11-29.root
  // PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_4pC_nSRC_cutBadPhiAndCombM5_wtaAxis_jetTrigger_preprocessed_2022-11-29.root
  
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
  //  Binning configuration for the integral calculation
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
  
  // Select the types of energy-energy correlators to integrate
  bool integrateEnergyEnergyCorrelator[EECHistogramManager::knEnergyEnergyCorrelatorTypes];
  integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator] = true;
  integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPt] = false;
  integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorUncorrected] = false;
  integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPtUncorrected] = false;
  
  // Determine an inxed of a drawn energy-energy correlator
  int lowestEnergyEnergyCorrelatorIndex = -1;
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    if(integrateEnergyEnergyCorrelator[iEnergyEnergyCorrelator]){
      lowestEnergyEnergyCorrelatorIndex = iEnergyEnergyCorrelator;
      break;
    }
  }
    
  // Create and setup a new histogram manager to project and handle the histograms
  EECHistogramManager *histograms = new EECHistogramManager(inputFile,card);
  
  // Choose the energy-energy correlator types to load
  histograms->SetLoadEnergyEnergyCorrelators(integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator]);
  histograms->SetLoadEnergyEnergyCorrelatorsJetPt(integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPt]);
  histograms->SetLoadEnergyEnergyCorrelatorsUncorrected(integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorUncorrected]);
  histograms->SetLoadEnergyEnergyCorrelatorsJetPtUncorrected(integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPtUncorrected]);
  
  // Choose the bin ranges
  histograms->SetCentralityBinRange(0,nCentralityBins);
  histograms->SetJetPtBinRangeEEC(0,nJetPtBinsEEC);
  histograms->SetTrackPtBinRangeEEC(0,nTrackPtBinsEEC);
  
  // Load the histograms from the file
  histograms->LoadProcessedHistograms();
  
  // Initialize the energy-energy correlator histogram array to NULL
  TH1D* hEnergyEnergyCorrelator[EECHistogramManager::knEnergyEnergyCorrelatorTypes][EECHistograms::knPairingTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC][EECHistograms::knSubeventCombinations+1];
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    for(int iPairingType = 0; iPairingType < EECHistograms::knPairingTypes; iPairingType++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iJetPt = 0; iJetPt < nJetPtBinsEEC+1; iJetPt++){
          for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
            for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventCombinations+1; iSubevent++){
              hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iPairingType][iCentrality][iJetPt][iTrackPt][iSubevent] = NULL;
            } // Subevent loop
          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Jet cone pairing type loop
  } // Energy-energy correlator loop
  
  // Get the histograms from the histogram manager and calculate integrals
  double normalizationFactor;
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only read the selected energy-energy correlator types
    if(!integrateEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
    
    for(int iPairingType = 0; iPairingType < EECHistograms::knPairingTypes; iPairingType++){
      
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iJetPt = 0; iJetPt < nJetPtBinsEEC+1; iJetPt++){
          for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
            
            // For normalization, read the histogram with all pair combinations
            hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iPairingType][iCentrality][iJetPt][iTrackPt][EECHistograms::knSubeventCombinations] = histograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iPairingType, EECHistograms::knSubeventCombinations);
            
            for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventCombinations; iSubevent++){
              hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iPairingType][iCentrality][iJetPt][iTrackPt][iSubevent] = histograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iPairingType, iSubevent);
              
            } // Subevent loop
          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Jet cone pairing type loop
  } // Energy-energy correlator loop
  
  
  double signalFakeIntegral, fakeFakeIntegral, reflectedConeIntegral;
  TString centralityString, trackPtString, jetPtString;
  int lowIntegralBin = 1;
  int highIntegralBin = hEnergyEnergyCorrelator[lowestEnergyEnergyCorrelatorIndex][EECHistograms::kSameJetPair][0][0][0][EECHistograms::kPythiaHydjet]->FindBin(0.4);
  
  // Get the histograms from the histogram manager and calculate integrals
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only read the selected energy-energy correlator types
    if(!integrateEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
    
    // Print all the numbers to an array that can be copied to a post-processing code
    cout << " double " << histograms->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator) << "BackgroundScaleFactors[" << nCentralityBins << "][" << nJetPtBinsEEC+1 << "][" << nTrackPtBinsEEC << "] = {" << endl;
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      
      // Set the centrality information for legends and figure saving
      cout << " // " << Form("Centrality %.0f-%.0f%%", histograms->GetCentralityBinBorder(iCentrality), histograms->GetCentralityBinBorder(iCentrality+1)) << endl;
      cout << "{";
      
      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC+1; iJetPt++){
        
        cout << "{";
        
        for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
          
          // Calculate the integrals and the ratio
          signalFakeIntegral = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][EECHistograms::kSameJetPair][iCentrality][iJetPt][iTrackPt][EECHistograms::kPythiaHydjet]->Integral(lowIntegralBin, highIntegralBin, "width");
          fakeFakeIntegral = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][EECHistograms::kSameJetPair][iCentrality][iJetPt][iTrackPt][EECHistograms::kHydjetHydjet]->Integral(lowIntegralBin, highIntegralBin, "width");
          reflectedConeIntegral = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][EECHistograms::kSignalReflectedConePair][iCentrality][iJetPt][iTrackPt][EECHistograms::knSubeventCombinations]->Integral(lowIntegralBin, highIntegralBin, "width");
          
          // Print the number to the array:
          cout << reflectedConeIntegral / (signalFakeIntegral+fakeFakeIntegral);
          
          // Add the proper syntax
          if(iTrackPt < nTrackPtBinsEEC-1){
            cout << ",";
          } else {
            cout << "}";
          }
        } // Track pT loop
        
        if(iJetPt < nJetPtBinsEEC){
          cout << ", // " << Form("%.0f < jet pT < %.0f GeV", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1)) << endl;
        } else {
          cout << "}";
          if(iCentrality < nCentralityBins-1){
            cout << ", // " << "120 GeV < jet pT" << endl;
          } else {
            cout << " // " << "120 GeV < jet pT" << endl;
          }
        }
        
      } // Jet pT loop
    } // Centrality loop
    cout << "};" << endl;
  } // Energy-energy correlator loop
}
