#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"
#include "../src/EECHistograms.h"

/*
 * Macro for comparing histograms from two different sets of files
 */
void compareFileSets(){
  
  // Enumeration for different systems
  enum enumSystem{kPp, kPbPb, knSystems};

  // Enumeration for different predefined file comparison sets
  enum enumComparisonSet{kJetAxis, kRadiusWTA, kRadiusEScheme, kBigRadiusJetAxis, knComparisonSets};

  const bool useData = true; // Compare either pp and PbPb data or Pythia and Pythia+Hydjet simulations

  // Configuration on how many files are in sets, and how many file sets will be included in the configuration
  const int nFileSets = 2;
  const int nFilesInSet = 2;

  // For data, we can select a predefined comparison set. Set -1 for manual configuration
  const int comparisonSet = kRadiusEScheme;

  // Define the file names for the studied files
  TString fileName[knSystems][nFileSets][nFilesInSet];

  // Naming for axes in figures
  TString fileSetDescription[nFileSets];
  fileSetDescription[0] = "n = 1";
  fileSetDescription[1] = "n = 2";

  // The manual configuration here is overwritten for predefined sets
  TString fileInSetDescription[nFilesInSet];
  fileInSetDescription[0] = "E-scheme";
  fileInSetDescription[1] = "WTA";

  // The manual configuration here is overwritten for predefined sets
  

  if(useData){

    if(comparisonSet < 0){
      // Manual configuration

      // Files for pp
      fileName[kPp][0][0] = "data/eschemeAxis/ppData_pfJets_eschemeAxis_nominalEnergyWeight_processed_2025-02-21.root";
      fileName[kPp][1][0] = "data/eschemeAxis/ppData_pfJets_eschemeAxis_energyWeightSquared_processed_2025-02-21.root";
      fileName[kPp][0][1] = "data/eschemeAxis/ppData_pfJets_wtaAxis_nominalEnergyWeight_processed_2024-01-23.root";
      fileName[kPp][1][1] = "data/eschemeAxis/ppData_pfJets_wtaAxis_energyWeightSquared_processed_2024-01-23.root";

      // Files for PbPb
      fileName[kPbPb][0][0] = "data/eschemeAxis/eecAnalysis_akFlowJet_eschemeAxis_nominalEnergyWeight_mixedConeSubtracted_processed_2025-02-25.root";
      fileName[kPbPb][1][0] = "data/eschemeAxis/eecAnalysis_akFlowJet_eschemeAxis_energyWeightSquared_mixedConeSubtracted_processed_2025-02-25.root";
      fileName[kPbPb][0][1] = "data/eschemeAxis/eecAnalysis_akFlowJet_wtaAxis_nominalEnergyWeight_mixedConeSubtracted_processed_2024-03-31.root";
      fileName[kPbPb][1][1] = "data/eschemeAxis/eecAnalysis_akFlowJet_wtaAxis_energyWeightSquared_mixedConeSubtracted_processed_2024-03-31.root";

    } else if (comparisonSet == kJetAxis){
      // Comparison between different jet axes for R=0.4 cones

      // Files for pp
      fileName[kPp][0][0] = "data/eschemeAxis/ppData_pfJets_eschemeAxis_nominalEnergyWeight_processed_2025-02-21.root";
      fileName[kPp][1][0] = "data/eschemeAxis/ppData_pfJets_eschemeAxis_energyWeightSquared_processed_2025-02-21.root";
      fileName[kPp][0][1] = "data/eschemeAxis/ppData_pfJets_wtaAxis_nominalEnergyWeight_processed_2024-01-23.root";
      fileName[kPp][1][1] = "data/eschemeAxis/ppData_pfJets_wtaAxis_energyWeightSquared_processed_2024-01-23.root";

      // Files for PbPb
      fileName[kPbPb][0][0] = "data/eschemeAxis/eecAnalysis_akFlowJet_eschemeAxis_nominalEnergyWeight_mixedConeSubtracted_processed_2025-02-25.root";
      fileName[kPbPb][1][0] = "data/eschemeAxis/eecAnalysis_akFlowJet_eschemeAxis_energyWeightSquared_mixedConeSubtracted_processed_2025-02-25.root";
      fileName[kPbPb][0][1] = "data/eschemeAxis/eecAnalysis_akFlowJet_wtaAxis_nominalEnergyWeight_mixedConeSubtracted_processed_2024-03-31.root";
      fileName[kPbPb][1][1] = "data/eschemeAxis/eecAnalysis_akFlowJet_wtaAxis_energyWeightSquared_mixedConeSubtracted_processed_2024-03-31.root";

      fileInSetDescription[0] = "E-scheme";
      fileInSetDescription[1] = "WTA";

    } else if (comparisonSet == kRadiusWTA){
      // Comparison between R=0.4 and 0.8 cone radii for WTA axis

      // Files for pp
      fileName[kPp][0][0] = "data/eschemeAxis/ppData_pfJets_wtaAxis_radius8_nominalEnergyWeight_processed_2025-02-27.root";
      fileName[kPp][1][0] = "data/eschemeAxis/ppData_pfJets_wtaAxis_radius8_energyWeightSquared_processed_2025-02-27.root";
      fileName[kPp][0][1] = "data/eschemeAxis/ppData_pfJets_wtaAxis_nominalEnergyWeight_processed_2024-01-23.root";
      fileName[kPp][1][1] = "data/eschemeAxis/ppData_pfJets_wtaAxis_energyWeightSquared_processed_2024-01-23.root";

      // Files for PbPb
      fileName[kPbPb][0][0] = "data/eschemeAxis/eecAnalysis_akFlowJet_wtaAxis_nominalEnergyWeight_mixedConeSubtracted_radius8_processed_2025-02-27.root";
      fileName[kPbPb][1][0] = "data/eschemeAxis/eecAnalysis_akFlowJet_wtaAxis_energyWeightSquared_mixedConeSubtracted_radius8_processed_2025-02-27.root";
      fileName[kPbPb][0][1] = "data/eschemeAxis/eecAnalysis_akFlowJet_wtaAxis_nominalEnergyWeight_mixedConeSubtracted_processed_2024-03-31.root";
      fileName[kPbPb][1][1] = "data/eschemeAxis/eecAnalysis_akFlowJet_wtaAxis_energyWeightSquared_mixedConeSubtracted_processed_2024-03-31.root";

      fileInSetDescription[0] = "WTA, R=0.8";
      fileInSetDescription[1] = "WTA, R=0.4";

    } else if (comparisonSet == kRadiusEScheme){
      // Comparison between R=0.4 and 0.8 cone radii for E-scheme axis

      // Files for pp
      fileName[kPp][0][0] = "data/eschemeAxis/ppData_pfJets_eschemeAxis_radius8_nominalEnergyWeight_processed_2025-02-27.root";
      fileName[kPp][1][0] = "data/eschemeAxis/ppData_pfJets_eschemeAxis_radius8_energyWeightSquared_processed_2025-02-27.root";
      fileName[kPp][0][1] = "data/eschemeAxis/ppData_pfJets_eschemeAxis_nominalEnergyWeight_processed_2025-02-21.root";
      fileName[kPp][1][1] = "data/eschemeAxis/ppData_pfJets_eschemeAxis_energyWeightSquared_processed_2025-02-21.root";

      //fileName[kPp][0][0] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_nominalEnergyWeight_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_processed_2025-03-10.root";
      //fileName[kPp][1][0] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_energyWeightSquared_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_processed_2025-03-10.root";
      //fileName[kPp][0][1] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_nominalEnergyWeight_jetDeltaAxis_leadingParticleFlag_processed_2025-03-10.root";
      //fileName[kPp][1][1] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_energyWeightSquared_jetDeltaAxis_leadingParticleFlag_processed_2025-03-10.root";

      // Files for PbPb
      fileName[kPbPb][0][0] = "data/eschemeAxis/eecAnalysis_akFlowJet_eschemeAxis_nominalEnergyWeight_mixedConeSubtracted_radius8_processed_2025-02-27.root";
      fileName[kPbPb][1][0] = "data/eschemeAxis/eecAnalysis_akFlowJet_eschemeAxis_energyWeightSquared_mixedConeSubtracted_radius8_processed_2025-02-27.root";
      fileName[kPbPb][0][1] = "data/eschemeAxis/eecAnalysis_akFlowJet_eschemeAxis_nominalEnergyWeight_mixedConeSubtracted_processed_2025-02-25.root";
      fileName[kPbPb][1][1] = "data/eschemeAxis/eecAnalysis_akFlowJet_eschemeAxis_energyWeightSquared_mixedConeSubtracted_processed_2025-02-25.root";

      fileInSetDescription[0] = "E-scheme, R=0.8";
      fileInSetDescription[1] = "E-scheme, R=0.4";

    } else if (comparisonSet == kBigRadiusJetAxis){
      // Comparison between E-scheme and WTA axes for jet radius 0.8

      // Files for pp
      fileName[kPp][0][0] = "data/eschemeAxis/ppData_pfJets_eschemeAxis_radius8_nominalEnergyWeight_processed_2025-02-27.root";
      fileName[kPp][1][0] = "data/eschemeAxis/ppData_pfJets_eschemeAxis_radius8_energyWeightSquared_processed_2025-02-27.root";
      fileName[kPp][0][1] = "data/eschemeAxis/ppData_pfJets_wtaAxis_radius8_nominalEnergyWeight_processed_2025-02-27.root";
      fileName[kPp][1][1] = "data/eschemeAxis/ppData_pfJets_wtaAxis_radius8_energyWeightSquared_processed_2025-02-27.root";

      // Files for PbPb
      fileName[kPbPb][0][0] = "data/eschemeAxis/eecAnalysis_akFlowJet_eschemeAxis_nominalEnergyWeight_mixedConeSubtracted_radius8_processed_2025-02-27.root";
      fileName[kPbPb][1][0] = "data/eschemeAxis/eecAnalysis_akFlowJet_eschemeAxis_energyWeightSquared_mixedConeSubtracted_radius8_processed_2025-02-27.root";
      fileName[kPbPb][0][1] = "data/eschemeAxis/eecAnalysis_akFlowJet_wtaAxis_nominalEnergyWeight_mixedConeSubtracted_radius8_processed_2025-02-27.root";
      fileName[kPbPb][1][1] = "data/eschemeAxis/eecAnalysis_akFlowJet_wtaAxis_energyWeightSquared_mixedConeSubtracted_radius8_processed_2025-02-27.root";

      fileInSetDescription[0] = "E-scheme, R=0.8";
      fileInSetDescription[1] = "WTA, R=0.8";

    }


  } else {

    // Files for Pythia
    fileName[kPp][0][0] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_nominalEnergyWeight_processed_2025-02-19.root";
    fileName[kPp][1][0] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_energyWeightSquared_processed_2025-02-19.root";
    fileName[kPp][0][1] = "data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_optimizedUnfoldingBins_nominalSmear_truthReference_processed_2024-01-11.root";
    fileName[kPp][1][1] = "data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_optimizedUnfoldingBins_energyWeightSquared_nominalSmear_truthReference_processed_2024-01-10.root";

    // Files for Pythia+Hydjet
    fileName[kPbPb][0][0] = "data/eschemeAxis/PbPbMC2018_GenGen_akFlowJets_eschemeAxis_4pCentShift_cutBadPhi_nominalEnergyWeight_onlyReflectedCone_processed_2025-02-24.root";
    fileName[kPbPb][1][0] = "data/eschemeAxis/PbPbMC2018_GenGen_akFlowJets_eschemeAxis_4pCentShift_cutBadPhi_energyWeightSquared_onlyReflectedCone_processed_2025-02-24.root";
    fileName[kPbPb][0][1] = "data/PbPbMC2018_GenGen_akFlowJets_4pCentShift_cutBadPhi_optimizedUnfoldingBins_nominalSmear_truthReference_processed_2024-01-16.root";
    fileName[kPbPb][1][1] = "data/PbPbMC2018_GenGen_akFlowJets_4pCentShift_cutBadPhi_optimizedUnfoldingBins_energyWeightSquared_nominalSmear_truthReference_processed_2024-01-18.root";

  }

  // Give a descriptions for the ratios depending on what is defined in the sets
  TString ratioDescription = Form("#frac{%s}{%s}", fileInSetDescription[1].Data(), fileInSetDescription[0].Data());
  TString doubleRatioDescription = Form("frac{PbPb (%s/%s)}{pp (%s/%s)}", fileInSetDescription[1].Data(), fileInSetDescription[0].Data(), fileInSetDescription[1].Data(), fileInSetDescription[0].Data());
  
  // Open the files and check that they exist
  TFile* inputFile[knSystems][nFileSets][nFilesInSet];
  EECCard* card[knSystems][nFileSets][nFilesInSet];
  for(int iSystem = 0; iSystem < knSystems; iSystem++){
    for(int iFileSet = 0; iFileSet < nFileSets; iFileSet++){
      for(int iFileInSet = 0; iFileInSet < nFilesInSet; iFileInSet++){
    
        inputFile[iSystem][iFileSet][iFileInSet] = TFile::Open(fileName[iSystem][iFileSet][iFileInSet]);
    
        if(inputFile[iSystem][iFileSet][iFileInSet] == NULL){
          cout << "Error! The file " << fileName[iSystem][iFileSet][iFileInSet].Data() << " does not exist!" << endl;
          cout << "Maybe you forgot the data/ folder path?" << endl;
          cout << "Will not execute the code" << endl;
          return;
        }

        card[iSystem][iFileSet][iFileInSet]  = new EECCard(inputFile[iSystem][iFileSet][iFileInSet]);
      } // File in set loop
    } // File set loop
  } // System loop
  
  // ====================================================
  //               Binning configuration
  // ====================================================
  
  // Find the number of bins from the reference card
  const int nJetPtBinsEEC = card[kPbPb][0][0]->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card[kPbPb][0][0]->GetNTrackPtBinsEEC();
  const int nCentralityBins = card[kPbPb][0][0]->GetNCentralityBins();
  
  std::vector<std::pair<int,int>> comparedCentralityBin;
  comparedCentralityBin.push_back(std::make_pair(0,10));

  std::vector<std::pair<double,double>> comparedJetPtBin;
  comparedJetPtBin.push_back(std::make_pair(120,140));
  //comparedJetPtBin.push_back(std::make_pair(140,160));
  //comparedJetPtBin.push_back(std::make_pair(160,180));
  //comparedJetPtBin.push_back(std::make_pair(180,200));

  std::vector<double> comparedTrackPtBin;
  comparedTrackPtBin.push_back(1.0);
  //comparedTrackPtBin.push_back(2.0);
  //comparedTrackPtBin.push_back(3.0);

  // For Monte Carlo, a 4% centrality shift is applied
  if(card[kPbPb][0][0]->GetDataType().Contains("MC")){
    for(auto& centralityBin : comparedCentralityBin){
      centralityBin.first += 4;
      centralityBin.second += 4;
    }
  }

  // Choose the type of drawn energy-energy correlator
  // EECHistogramManager::kEnergyEnergyCorrelatorNormalized = Normalized energy-energy correlator
  // EECHistogramManager::kEnergyEnergyCorrelatorBackground = Estimated background
  // EECHistogramManager::kEnergyEnergyCorrelatorSignal = Background subtracted, but not unfolded energy-energy correlator
  // EECHistogramManager::kEnergyEnergyCorrelatorUnfolded = Unfolded energy-energy correlator
  // EECHistogramManager::kEnergyEnergyCorrelatorBackgroundAfterUnfolding = Estimated background after unfolding
  // EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal = Unfolded energy-energy correlator signal
  // EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels = Raw energy-energy correlator
  int drawnEnergyEnergyCorrelator = EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels;

  TString energyCorrelatorType[EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels+1];
  energyCorrelatorType[EECHistogramManager::kEnergyEnergyCorrelatorNormalized] = "Normalized";
  energyCorrelatorType[EECHistogramManager::kEnergyEnergyCorrelatorBackground] = "Background only";
  energyCorrelatorType[EECHistogramManager::kEnergyEnergyCorrelatorSignal] = "Signal";
  energyCorrelatorType[EECHistogramManager::kEnergyEnergyCorrelatorUnfolded] = "Unfolded";
  energyCorrelatorType[EECHistogramManager::kEnergyEnergyCorrelatorBackgroundAfterUnfolding] = "Unfolded background";
  energyCorrelatorType[EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal] = "Unfolded signal";
  energyCorrelatorType[EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels] = "Raw";

  // Select which pairing type and subevent to look at
  // These are only applied if drawnEnergyEnergyCorrelator is EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels

  // Available pairing types:
  //  kSameJetPair
  //  kSignalReflectedConePair
  //  kReflectedConePair
  //  kSignalMixedConePair
  //  kReflectedMixedConePair
  //  kMixedConePair
  //  kSignalSecondMixedConePair
  //  kReflectedSecondMixedConePair
  //  kMixedMixedConePair
  //  kSecondMixedConePair
  int pairingType = EECHistograms::kSameJetPair;

  // Available subevent combinations:
  //  kPythiaPythia
  //  kPythiaHydjet
  //  kHydjetPythia
  //  kHydjetHydjet
  //  knSubeventCombinations
  int subevent = EECHistograms::knSubeventCombinations;

  // The subevent decomposition is relevant only for MC
  if(useData) subevent = EECHistograms::knSubeventCombinations;

  // ====================================================
  //                Drawing configuration
  // ====================================================

  bool drawPp = true; // Include pp comparisons for single ratio histograms
  bool drawSingleRatios = true; // Draw the single ratio histograms
  bool drawDoubleRatios = false; // Draw the double ratio histograms
  
  // Figure saving
  const bool saveFigures = false;  // Save figures
  TString saveComment = "_jetAxisDifference_signal";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures

  // Drawing configuration
  std::pair<double, double> ratioZoom = std::make_pair(0.8, 1.2);
  std::pair<double, double> eecZoom = std::make_pair(0.2, 30);
  std::pair<double, double> doubleRatioZoom = std::make_pair(0.8, 1.2);
  const bool automaticZoom = true;

  // pp files do not have background, so do not draw pp if background is studied
  if(drawnEnergyEnergyCorrelator == EECHistogramManager::kEnergyEnergyCorrelatorBackground){
    drawPp = false;
    drawDoubleRatios = false;
  }

  // ====================================================
  //          Sanity check for unfolded signals
  // ====================================================

  // Sanity checks for unfolded distributions. Ensure that all the selected bins actually exist for unfolded distributions.
  if(drawnEnergyEnergyCorrelator != EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels){
    if(drawnEnergyEnergyCorrelator > EECHistogramManager::kEnergyEnergyCorrelatorSignal){
      for(int iSystem = 0; iSystem < knSystems; iSystem++){
        for(int iFileSet = 0; iFileSet < nFileSets; iFileSet++){
          for(int iFileInSet = 0; iFileInSet < nFilesInSet; iFileInSet++){

            // Sanity check for jet pT bins
            for(auto jetPtBin : comparedJetPtBin){
              if(card[iSystem][iFileSet][iFileInSet]->FindBinIndexJetPtEEC(jetPtBin) < card[iSystem][iFileSet][iFileInSet]->GetFirstUnfoldedJetPtBin() || card[iSystem][iFileSet][iFileInSet]->FindBinIndexJetPtEEC(jetPtBin) > card[iSystem][iFileSet][iFileInSet]->GetLastUnfoldedJetPtBin()){
                cout << "ERROR! Jet pT bin " << jetPtBin.first << "-" << jetPtBin.second << " does not exist in file " << fileName[iSystem][iFileSet][iFileInSet].Data() << endl;
                cout << "Please only choose jet pT bins that are included in the input files." << endl;
                return;
              }
            }

            // Sanity check for track pT bins
            for(auto trackPtBin : comparedTrackPtBin){
              if(card[iSystem][iFileSet][iFileInSet]->GetBinIndexTrackPtEEC(trackPtBin) < card[iSystem][iFileSet][iFileInSet]->GetFirstUnfoldedTrackPtBin() || card[iSystem][iFileSet][iFileInSet]->GetBinIndexTrackPtEEC(trackPtBin) > card[iSystem][iFileSet][iFileInSet]->GetLastUnfoldedTrackPtBin()){
                cout << "ERROR! Track pT cut > " << trackPtBin << " GeV does not exist in file " << fileName[iSystem][iFileSet][iFileInSet].Data() << endl;
                cout << "Please only choose track pT bins that are included in the input files." << endl;
                return;
              }
            } 
          } // File in set loop for sanity check
        } // File set loop for input sanity check
      } // System loop
    } // Only unfolded distributions
  } // Not raw distribution

  // ====================================================
  //                  Histogram loading
  // ====================================================
  
  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager* histograms[knSystems][nFileSets][nFilesInSet];
  for(int iSystem = 0; iSystem < knSystems; iSystem++){
    for(int iFileSet = 0; iFileSet < nFileSets; iFileSet++){
      for(int iFileInSet = 0; iFileInSet < nFilesInSet; iFileInSet++){

        histograms[iSystem][iFileSet][iFileInSet] = new EECHistogramManager(inputFile[iSystem][iFileSet][iFileInSet], card[iSystem][iFileSet][iFileInSet]);

        // Choose the energy-energy correlator types to load
        histograms[iSystem][iFileSet][iFileInSet]->SetLoadEnergyEnergyCorrelators(true);

        // Choose the bin ranges
        histograms[iSystem][iFileSet][iFileInSet]->SetCentralityBinRange(0, card[iSystem][iFileSet][iFileInSet]->GetNCentralityBins() - 1);
        histograms[iSystem][iFileSet][iFileInSet]->SetJetPtBinRangeEEC(0, card[iSystem][iFileSet][iFileInSet]->GetNJetPtBinsEEC() - 1);
        histograms[iSystem][iFileSet][iFileInSet]->SetTrackPtBinRangeEEC(0, card[iSystem][iFileSet][iFileInSet]->GetNTrackPtBinsEEC() - 1);

        // Load the histograms from the file
        histograms[iSystem][iFileSet][iFileInSet]->LoadProcessedHistograms();
      } // File in set loop for histogram loading
    } // File set loop for histogram loading
  } // System loop for histogram loading

  // Energy-energy correlator histograms
  TH1D* hEnergyEnergyCorrelator[nFileSets][nFilesInSet][nCentralityBins+1][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hEnergyEnergyCorrelatorRatio[nFileSets][nFilesInSet][nCentralityBins+1][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hEnergyEnergyCorrelatorDoubleRatio[nFileSets][nFilesInSet][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iFileSet = 0; iFileSet < nFileSets; iFileSet++){
    for(int iFileInSet = 0; iFileInSet < nFilesInSet; iFileInSet++){
      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
          for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
            hEnergyEnergyCorrelator[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt] = NULL;
            hEnergyEnergyCorrelatorRatio[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt] = NULL;
            hEnergyEnergyCorrelatorDoubleRatio[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt] = NULL;
          } // Centrality loop
          hEnergyEnergyCorrelator[iFileSet][iFileInSet][nCentralityBins][iJetPt][iTrackPt] = NULL;
          hEnergyEnergyCorrelatorRatio[iFileSet][iFileInSet][nCentralityBins][iJetPt][iTrackPt] = NULL;
        } // Track pT loop
      } // Jet pT loop
    } // File in set loop
  } // File set loop
  
  // Helper variables
  std::pair<double, double> drawingRange = std::make_pair(0.008, 0.39);
  double epsilon = 0.0001;
  int lowNormalizationBin, highNormalizationBin;
  int iTrackPt, iTrackPtReference;
  int iJetPt, iJetPtReference;
  int iCentrality, iCentralityReference;

  // Get the histograms from the histogram manager and normalize the signal histograms to one
  for(int iFileSet = 0; iFileSet < nFileSets; iFileSet++){
    for(int iFileInSet = 0; iFileInSet < nFilesInSet; iFileInSet++){
      for(auto jetPtBin : comparedJetPtBin){
        for(auto trackPtBin : comparedTrackPtBin){

          // Find the proper binning and express it in term of the bins in the first file
          iJetPt = card[kPp][iFileSet][iFileInSet]->FindBinIndexJetPtEEC(jetPtBin);
          iTrackPt = card[kPp][iFileSet][iFileInSet]->GetBinIndexTrackPtEEC(trackPtBin);
          iJetPtReference = card[kPbPb][0][0]->FindBinIndexJetPtEEC(jetPtBin);
          iTrackPtReference = card[kPbPb][0][0]->GetBinIndexTrackPtEEC(trackPtBin);

          // Load the selected energy-energy correlator histograms for pp
          if(drawnEnergyEnergyCorrelator == EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels){
            hEnergyEnergyCorrelator[iFileSet][iFileInSet][nCentralityBins][iJetPtReference][iTrackPtReference] = histograms[kPp][iFileSet][iFileInSet]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPt, iTrackPt, pairingType);
          } else {
            hEnergyEnergyCorrelator[iFileSet][iFileInSet][nCentralityBins][iJetPtReference][iTrackPtReference] = histograms[kPp][iFileSet][iFileInSet]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPt, iTrackPt, drawnEnergyEnergyCorrelator);
          }

          // Normalize the distributions to one in the drawingRange
          lowNormalizationBin = hEnergyEnergyCorrelator[iFileSet][iFileInSet][nCentralityBins][iJetPtReference][iTrackPtReference]->GetXaxis()->FindBin(drawingRange.first + epsilon);
          highNormalizationBin = hEnergyEnergyCorrelator[iFileSet][iFileInSet][nCentralityBins][iJetPtReference][iTrackPtReference]->GetXaxis()->FindBin(drawingRange.second - epsilon);

          hEnergyEnergyCorrelator[iFileSet][iFileInSet][nCentralityBins][iJetPtReference][iTrackPtReference]->Scale(1 / hEnergyEnergyCorrelator[iFileSet][iFileInSet][nCentralityBins][iJetPtReference][iTrackPtReference]->Integral(lowNormalizationBin, highNormalizationBin, "width"));

          // Do the same for all the selected centrality bins for PbPb
          iJetPt = card[kPbPb][iFileSet][iFileInSet]->FindBinIndexJetPtEEC(jetPtBin);
          iTrackPt = card[kPbPb][iFileSet][iFileInSet]->GetBinIndexTrackPtEEC(trackPtBin);
          for(auto centralityBin : comparedCentralityBin){
            iCentrality = card[kPbPb][iFileSet][iFileInSet]->FindBinIndexCentrality(centralityBin);
            iCentralityReference = card[kPbPb][0][0]->FindBinIndexCentrality(centralityBin);

            // Load the selected energy-energy correlator histograms for PbPb
            if(drawnEnergyEnergyCorrelator == EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels){
              hEnergyEnergyCorrelator[iFileSet][iFileInSet][iCentrality][iJetPtReference][iTrackPtReference] = histograms[kPbPb][iFileSet][iFileInSet]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, pairingType, subevent);
            } else {
              hEnergyEnergyCorrelator[iFileSet][iFileInSet][iCentrality][iJetPtReference][iTrackPtReference] = histograms[kPbPb][iFileSet][iFileInSet]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, drawnEnergyEnergyCorrelator);
            }

            if(hEnergyEnergyCorrelator[iFileSet][iFileInSet][iCentrality][iJetPtReference][iTrackPtReference] == NULL) cout << "We have a NULL" << endl;

            // Normalize the distributions to one in the drawingRange
            lowNormalizationBin = hEnergyEnergyCorrelator[iFileSet][iFileInSet][iCentrality][iJetPtReference][iTrackPtReference]->GetXaxis()->FindBin(drawingRange.first + epsilon);
            highNormalizationBin = hEnergyEnergyCorrelator[iFileSet][iFileInSet][iCentrality][iJetPtReference][iTrackPtReference]->GetXaxis()->FindBin(drawingRange.second - epsilon);

            hEnergyEnergyCorrelator[iFileSet][iFileInSet][iCentrality][iJetPtReference][iTrackPtReference]->Scale(1 / hEnergyEnergyCorrelator[iFileSet][iFileInSet][iCentrality][iJetPtReference][iTrackPtReference]->Integral(lowNormalizationBin, highNormalizationBin, "width"));
            
          } // Centrality loop
        } // Track pT loop
      } // Jet pT loop
    } // File in set loop
  } // File set loop

  // ====================================================
  //                 Calculating ratios
  // ====================================================

  // After all the histograms have been read, calculate the ratios within each set with respect to first file index in that set
  for(int iFileSet = 0; iFileSet < nFileSets; iFileSet++){
    for(int iFileInSet = 0; iFileInSet < nFilesInSet; iFileInSet++){
      for(auto jetPtBin : comparedJetPtBin){
        iJetPt = card[kPbPb][0][0]->FindBinIndexJetPtEEC(jetPtBin);
        for(auto trackPtBin : comparedTrackPtBin){
          iTrackPt = card[kPbPb][0][0]->GetBinIndexTrackPtEEC(trackPtBin);

          hEnergyEnergyCorrelatorRatio[iFileSet][iFileInSet][nCentralityBins][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelator[iFileSet][iFileInSet][nCentralityBins][iJetPt][iTrackPt]->Clone(Form("eecRatioPp%d%d%d%d", iFileSet, iFileInSet, iJetPt, iTrackPt));
          hEnergyEnergyCorrelatorRatio[iFileSet][iFileInSet][nCentralityBins][iJetPt][iTrackPt]->Divide(hEnergyEnergyCorrelator[iFileSet][0][nCentralityBins][iJetPt][iTrackPt]);

          for(auto centralityBin: comparedCentralityBin){
            iCentrality = card[kPbPb][0][0]->FindBinIndexCentrality(centralityBin);

            hEnergyEnergyCorrelatorRatio[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelator[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt]->Clone(Form("eecRatio%d%d%d%d%d", iFileSet, iFileInSet, iCentrality, iJetPt, iTrackPt));
            hEnergyEnergyCorrelatorRatio[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt]->Divide(hEnergyEnergyCorrelator[iFileSet][0][iCentrality][iJetPt][iTrackPt]);

          } // Centrality loop
        } // Track pT loop
      } // Jet pT loop
    } // File in set loop
  } // File set loop

  // We also want to calculate double ratios of PbPb to pp ratios of the single ratios
  for(int iFileSet = 0; iFileSet < nFileSets; iFileSet++){
    for(int iFileInSet = 0; iFileInSet < nFilesInSet; iFileInSet++){
      for(auto jetPtBin : comparedJetPtBin){
        iJetPt = card[kPbPb][0][0]->FindBinIndexJetPtEEC(jetPtBin);
        for(auto trackPtBin : comparedTrackPtBin){
          iTrackPt = card[kPbPb][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
          for(auto centralityBin: comparedCentralityBin){
            iCentrality = card[kPbPb][0][0]->FindBinIndexCentrality(centralityBin);

            hEnergyEnergyCorrelatorDoubleRatio[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelatorRatio[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt]->Clone(Form("eecDoubleRatio%d%d%d%d%d", iFileSet, iFileInSet, iCentrality, iJetPt, iTrackPt));
            hEnergyEnergyCorrelatorDoubleRatio[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt]->Divide(hEnergyEnergyCorrelatorRatio[iFileSet][iFileInSet][nCentralityBins][iJetPt][iTrackPt]);

          } // Centrality loop
        } // Track pT loop
      } // Jet pT loop
    } // File in set loop
  } // File set loop

  
  // ====================================================================
  //                Drawing the selected distributions
  // ====================================================================
  
  JDrawer* drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  
  drawer->SetLogX(true);

  TString compactJetPtString = "";
  TString compactTrackPtString = "";
  TString compactCentralityString = "";
  TString legendString;
  int markerStyle[] = {kFullSquare, kFullCircle, kOpenSquare, kOpenCircle, kFullCross};
  int color[] = {kBlack, kBlue, kGreen+3, kRed, kMagenta,kCyan,kOrange,kViolet+3,kPink-7,kSpring+3,kAzure-7};

  // Binning vectors
  std::vector<int> currentJetPtIndices;
  std::vector<int> currentTrackPtIndices;
  double minimumCandidate, maximumCandidate;
  TString systemForLegend = "";

  // Legend and line
  TLegend* legend;
  TLine* oneLine = new TLine(drawingRange.first, 1, drawingRange.second, 1);
  oneLine->SetLineStyle(2);

  // Draw the single ratio histograms
  if(drawSingleRatios){

    // For drawing file set comparisons, add flag for pp drawing to the centrality bin array
    std::vector<std::pair<int,int>> distributionCentralityBin;
    if(drawPp) distributionCentralityBin.push_back(std::make_pair(-1,-1));
    for(auto centralityBin : comparedCentralityBin){
      distributionCentralityBin.push_back(centralityBin);
    }

    for(auto jetPtBin : comparedJetPtBin){
      iJetPt = card[kPbPb][0][0]->FindBinIndexJetPtEEC(jetPtBin);
      for(auto trackPtBin : comparedTrackPtBin){
        iTrackPt = card[kPbPb][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
        for(auto centralityBin : distributionCentralityBin){
          if(centralityBin.first == -1){
            iCentrality = nCentralityBins;
            systemForLegend = card[kPp][0][0]->GetAlternativeDataType(false);
          } else {
            iCentrality = card[kPbPb][0][0]->FindBinIndexCentrality(centralityBin);
            systemForLegend = card[kPbPb][0][0]->GetAlternativeDataType(false);
          }
          
          // Create a new canvas for the plot
          drawer->CreateSplitCanvas();

          // Logarithmic EEC axis
          drawer->SetLogY(true);

          // Create the legend and add binning information to it
          legend = new TLegend(0.18,0.04,0.45,0.58);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);

          if(systemForLegend != "") legend->AddEntry((TObject*) 0, systemForLegend.Data(), "");
          if(centralityBin.first != -1) legend->AddEntry((TObject*) 0, Form("Centrality: %d-%d%%", centralityBin.first, centralityBin.second), "");
          legend->AddEntry((TObject*) 0, Form("%.0f < jet p_{T} < %.0f GeV", jetPtBin.first, jetPtBin.second), "");
          legend->AddEntry((TObject*) 0, Form("p_{T}^{ch} > %.1f GeV", trackPtBin), "");

          // Setup centrality, jet pT and track pT strings
          compactJetPtString = Form("_J=%.0f-%.0f", jetPtBin.first, jetPtBin.second);
          compactTrackPtString = Form("_T>%.1f",trackPtBin);
          compactTrackPtString.ReplaceAll(".","v");
          compactCentralityString = centralityBin.first == -1 ? "_pp" : Form("_C=%d-%d", centralityBin.first, centralityBin.second);


          // Set drawing style for all histograms
          for(int iFileSet = 0; iFileSet < nFileSets; iFileSet++){
            for(int iFileInSet = 0; iFileInSet < nFilesInSet; iFileInSet++){
              hEnergyEnergyCorrelator[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iFileInSet + nFilesInSet*iFileSet]);
              hEnergyEnergyCorrelator[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iFileInSet + nFilesInSet*iFileSet]);
              hEnergyEnergyCorrelator[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iFileInSet + nFilesInSet*iFileSet]); 
              hEnergyEnergyCorrelatorRatio[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iFileInSet + nFilesInSet*iFileSet]);
              hEnergyEnergyCorrelatorRatio[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iFileInSet + nFilesInSet*iFileSet]);
              hEnergyEnergyCorrelatorRatio[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iFileInSet + nFilesInSet*iFileSet]);
            } // File in set loop
          } // File set loop

          // Automatic zooming for the drawn histograms
          if(automaticZoom){
            eecZoom.first = 10000;
            eecZoom.second = 0;
            for(int iFileSet = 0; iFileSet < nFileSets; iFileSet++){
              for(int iFileInSet = 0; iFileInSet < nFilesInSet; iFileInSet++){
                hEnergyEnergyCorrelator[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
                minimumCandidate = hEnergyEnergyCorrelator[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt]->GetMinimum();
                maximumCandidate = hEnergyEnergyCorrelator[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt]->GetMaximum();
                if(minimumCandidate < eecZoom.first) eecZoom.first = minimumCandidate;
                if(maximumCandidate > eecZoom.second) eecZoom.second = maximumCandidate;
              } // File in set loop
            } // File set loop
            eecZoom.first = eecZoom.first / 2.0;
            eecZoom.second = eecZoom.second * 2.0;
          }

          // Set the x- and y-axis drawing ranges
          hEnergyEnergyCorrelator[0][0][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
          hEnergyEnergyCorrelator[0][0][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(eecZoom.first, eecZoom.second);

          // Draw the histograms to the upper canvas
          drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelator[0][0][iCentrality][iJetPt][iTrackPt], "#Deltar", Form("EEC (%s)", energyCorrelatorType[drawnEnergyEnergyCorrelator].Data()), " ", "p");

          for(int iFileSet = 0; iFileSet < nFileSets; iFileSet++){
            for(int iFileInSet = 0; iFileInSet < nFilesInSet; iFileInSet++){
              if(iFileSet == 0 && iFileInSet == 0) continue;
              hEnergyEnergyCorrelator[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt]->Draw("same,p");
            } // File in set loop
          } // File set loop

          // Add legends for drawn histograms
          for(int iFileSet = 0; iFileSet < nFileSets; iFileSet++){
            for(int iFileInSet = 0; iFileInSet < nFilesInSet; iFileInSet++){
              legend->AddEntry(hEnergyEnergyCorrelator[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt], Form("%s  %s", fileSetDescription[iFileSet].Data(), fileInSetDescription[iFileInSet].Data()), "p");
            } // File in set loop
          } // File set loop
  
          // Draw the legends to the upper pad
          legend->Draw();
          
          // Linear scale for the ratio
          drawer->SetLogY(false);
          
          // Set the drawing ranges
          hEnergyEnergyCorrelatorRatio[0][1][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
          hEnergyEnergyCorrelatorRatio[0][1][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

          // Draw the histograms
          drawer->SetGridY(true);

          drawer->DrawHistogramToLowerPad(hEnergyEnergyCorrelatorRatio[0][1][iCentrality][iJetPt][iTrackPt], "#Deltar", ratioDescription, " ", "p");
          for(int iFileSet = 0; iFileSet < nFileSets; iFileSet++){
            for(int iFileInSet = 1; iFileInSet < nFilesInSet; iFileInSet++){
              if(iFileSet == 0 && iFileInSet == 1) continue;
              hEnergyEnergyCorrelatorRatio[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt]->Draw("same");
            } // File in set loop
          } // File set loop
          drawer->SetGridY(false);
          
          // Save the figures to a file
          if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/eecComparison%s%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
          }
        } // Centrality looop
      } // Track pT loop
    } // Jet pT loop
  } // Drawing single ratios

  

  if(drawDoubleRatios){

    TString ratioDenominator;

    // To simplify ratio naming, cut the string from ',' and only use the second part in the ratio
    /*if(fileInSetDescription[0].Contains(",")){
      TObjArray* tokenizer = fileInSetDescription[0].Tokenize(",");

    } else {
      ratioDenominator = fileInSetDescription[0];
    }*/

    // Name for the double ratio axis
    TString doubleRatioAxis = Form("EEC (%s)  #frac{%s (%s / %s)}{%s (%s / %s)}", energyCorrelatorType[drawnEnergyEnergyCorrelator].Data(), card[kPbPb][0][0]->GetAlternativeDataType(false).Data(), fileInSetDescription[1].Data(), fileInSetDescription[0].Data(), card[kPp][0][0]->GetAlternativeDataType(false).Data(), fileInSetDescription[1].Data(), fileInSetDescription[0].Data());
    // "#frac{WTA (PbPb/pp)}{E-scheme (PbPb/pp)}";

    // Set a good drawing style for single canvas plots
    drawer->SetRelativeCanvasSize(1.1,1.1);
    drawer->SetLeftMargin(0.14);
    drawer->SetTopMargin(0.07);
    drawer->SetTitleOffsetY(1.7);
    drawer->SetTitleOffsetX(1.0);
    
    for(auto jetPtBin : comparedJetPtBin){
      iJetPt = card[kPbPb][0][0]->FindBinIndexJetPtEEC(jetPtBin);
      for(auto trackPtBin : comparedTrackPtBin){
        iTrackPt = card[kPbPb][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
        for(auto centralityBin : comparedCentralityBin){
          iCentrality = card[kPbPb][0][0]->FindBinIndexCentrality(centralityBin);
          
          // Linear scale for the double ratio
          drawer->SetLogY(false);

          // Create the legend and add binning information to it
          legend = new TLegend(0.18,0.62,0.45,0.9);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);

          legend->AddEntry((TObject*) 0, Form("Centrality: %d-%d%%", centralityBin.first, centralityBin.second), "");
          legend->AddEntry((TObject*) 0, Form("%.0f < jet p_{T} < %.0f GeV", jetPtBin.first, jetPtBin.second), "");
          legend->AddEntry((TObject*) 0, Form("p_{T}^{ch} > %.1f GeV", trackPtBin), "");

          // Setup centrality, jet pT and track pT strings
          compactJetPtString = Form("_J=%.0f-%.0f", jetPtBin.first, jetPtBin.second);
          compactTrackPtString = Form("_T>%.1f",trackPtBin);
          compactTrackPtString.ReplaceAll(".","v");
          compactCentralityString = Form("_C=%d-%d", centralityBin.first, centralityBin.second);

          // Set drawing style for all histograms
          for(int iFileSet = 0; iFileSet < nFileSets; iFileSet++){
            for(int iFileInSet = 0; iFileInSet < nFilesInSet; iFileInSet++){
              hEnergyEnergyCorrelatorDoubleRatio[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iFileInSet + nFilesInSet*iFileSet]);
              hEnergyEnergyCorrelatorDoubleRatio[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iFileInSet + nFilesInSet*iFileSet]);
              hEnergyEnergyCorrelatorDoubleRatio[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iFileInSet + nFilesInSet*iFileSet]);
            } // File in set loop
          } // File set loop

          // Set the x- and y-axis drawing ranges
          hEnergyEnergyCorrelatorDoubleRatio[1][1][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
          hEnergyEnergyCorrelatorDoubleRatio[1][1][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(doubleRatioZoom.first, doubleRatioZoom.second);

          // Draw the histograms
          drawer->DrawHistogram(hEnergyEnergyCorrelatorDoubleRatio[1][1][iCentrality][iJetPt][iTrackPt], "#Deltar", doubleRatioAxis.Data(), " ", "p");

          for(int iFileSet = 0; iFileSet < nFileSets; iFileSet++){
            for(int iFileInSet = 1; iFileInSet < nFilesInSet; iFileInSet++){
              if(iFileSet == 1 && iFileInSet == 1) continue;
              hEnergyEnergyCorrelatorDoubleRatio[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt]->Draw("same,p");
            } // File in set loop
          } // File set loop

          // Draw a line to one
          oneLine->Draw();

          // Add legends for drawn histograms
          for(int iFileSet = 0; iFileSet < nFileSets; iFileSet++){
            for(int iFileInSet = 1; iFileInSet < nFilesInSet; iFileInSet++){
              legend->AddEntry(hEnergyEnergyCorrelatorDoubleRatio[iFileSet][iFileInSet][iCentrality][iJetPt][iTrackPt], fileSetDescription[iFileSet].Data(), "p");
            } // File in set loop
          } // File set loop
  
          // Draw the legends to the upper pad
          legend->Draw();
          
          // Save the figures to a file
          if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/eecDoubleRatio%s%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
          }
        } // Centrality looop
      } // Track pT loop
    } // Jet pT loop
  } // Double ratio type loop

}
