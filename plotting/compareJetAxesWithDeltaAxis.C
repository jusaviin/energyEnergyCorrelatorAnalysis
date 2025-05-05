#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"
#include "../src/EECHistograms.h"

/*
 * Macro for comparing histograms from two different sets of files
 */
void compareJetAxesWithDeltaAxis(){
  
  // Enumeration for different jet axes
  enum enumJetAxis{kEscheme, kWTA, knJetAxes};

  // Enumeration for energy weight
  enum enumEnergyWeight{kNominalWeight, kSquaredWeight, knEnergyWeights};

  // Enumeration for jet radius
  enum enumJetRatios{kR0p4, kR0p8, knJetRadii};

  // Enumeration for vetoing jets in the region 0.4 < DeltaR < 0.8
  enum enumJetVeto{kNoVeto, kJetVeto, knVetoTypes};

  // Possibility to veto other jets within the extended radius around the jet axis
  const bool vetoCloseJets = false;

  // Define the file names for the studied files
  TString fileName[knVetoTypes][knJetRadii][knEnergyWeights][knJetAxes];

  // Naming for jet radii
  TString jetRadiusDescription[knJetRadii];
  jetRadiusDescription[kR0p4] = "R = 0.4";
  jetRadiusDescription[kR0p8] = "R = 0.8";

  TString jetRadiusSaveName[knJetRadii];
  jetRadiusSaveName[kR0p4] = "_R0p4";
  jetRadiusSaveName[kR0p8] = "_R0p8";

  // Naming for energy weights
  TString energyWeightDescription[knEnergyWeights];
  energyWeightDescription[kNominalWeight] = "n = 1";
  energyWeightDescription[kSquaredWeight] = "n = 2";

  TString compactEnergyWeightString[knEnergyWeights];
  compactEnergyWeightString[kNominalWeight] = "_nominalEnergyWeight";
  compactEnergyWeightString[kSquaredWeight] = "_energyWeightSquared";

  // Naming for jet axes
  TString jetAxisDescription[knJetAxes];
  jetAxisDescription[kEscheme] = "E-scheme";
  jetAxisDescription[kWTA] = "WTA";

  // Figure saving names for jet axes
  TString jetAxisSaveName[knJetAxes];
  jetAxisSaveName[kEscheme] = "_escheme";
  jetAxisSaveName[kWTA] = "_WTA";

  // Naming for veto types
  TString vetoDescription[knVetoTypes];
  vetoDescription[kNoVeto] = "";
  vetoDescription[kJetVeto] = "No other jets within 0.8";

  // Figure saving names for veto types
  TString vetoSaveName[knVetoTypes];
  vetoSaveName[kNoVeto] = "";
  vetoSaveName[kJetVeto] = "_vetoCloseJets";
  

  // Define files for each of these combinations. TODO: Update all files after merging in ACCRE is finished

  // Files for R = 0.4, no other jet veto
  fileName[kNoVeto][kR0p4][kNominalWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_nominalEnergyWeight_jetDeltaAxis_leadingParticleFlag_processed_2025-03-10.root";
  fileName[kNoVeto][kR0p4][kSquaredWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_energyWeightSquared_jetDeltaAxis_leadingParticleFlag_processed_2025-03-10.root";
  fileName[kNoVeto][kR0p4][kNominalWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_nominalEnergyWeight_jetDeltaAxis_leadingParticleFlag_processed_2025-03-10.root";
  fileName[kNoVeto][kR0p4][kSquaredWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_jetDeltaAxis_leadingParticleFlag_processed_2025-03-10.root";

  // Files for R = 0.8, no other jet veto
  fileName[kNoVeto][kR0p8][kNominalWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_nominalEnergyWeight_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_processed_2025-03-10.root";
  fileName[kNoVeto][kR0p8][kSquaredWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_energyWeightSquared_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_processed_2025-03-10.root";
  fileName[kNoVeto][kR0p8][kNominalWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_nominalEnergyWeight_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_processed_2025-03-10.root";
  fileName[kNoVeto][kR0p8][kSquaredWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_processed_2025-03-10.root";

  // Files for R = 0.4, veto if other jets within 0.8
  fileName[kJetVeto][kR0p4][kNominalWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_nominalEnergyWeight_jetDeltaAxis_jetRadius0p4_leadingParticleFlag_veto0p8Jets_processed_2025-03-16.root";
  fileName[kJetVeto][kR0p4][kSquaredWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_energyWeightSquared_jetDeltaAxis_jetRadius0p4_leadingParticleFlag_veto0p8Jets_processed_2025-03-16.root";
  fileName[kJetVeto][kR0p4][kNominalWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_nominalEnergyWeight_jetDeltaAxis_jetRadius0p4_leadingParticleFlag_veto0p8Jets_processed_2025-03-16.root";
  fileName[kJetVeto][kR0p4][kSquaredWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_jetDeltaAxis_jetRadius0p4_leadingParticleFlag_veto0p8Jets_processed_2025-03-16.root";

  // Files for R = 0.8, veto if other jets within 0.8
  fileName[kJetVeto][kR0p8][kNominalWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_nominalEnergyWeight_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_vetoCloseJets_processed_2025-03-14.root";
  fileName[kJetVeto][kR0p8][kSquaredWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_energyWeightSquared_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_vetoCloseJets_processed_2025-03-14.root";
  fileName[kJetVeto][kR0p8][kNominalWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_nominalEnergyWeight_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_vetoCloseJets_processed_2025-03-14.root";
  fileName[kJetVeto][kR0p8][kSquaredWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_vetoCloseJets_processed_2025-03-14.root";
  
  // Open the files and check that they exist
  TFile* inputFile[knVetoTypes][knJetRadii][knEnergyWeights][knJetAxes];
  EECCard* card[knVetoTypes][knJetRadii][knEnergyWeights][knJetAxes];
  for(int iVeto = 0; iVeto < knVetoTypes; iVeto++){
    for(int iJetRadius = 0; iJetRadius < knJetRadii; iJetRadius++){
      for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
        for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
    
          inputFile[iVeto][iJetRadius][iEnergyWeight][iJetAxis] = TFile::Open(fileName[iVeto][iJetRadius][iEnergyWeight][iJetAxis]);
    
          if(inputFile[iVeto][iJetRadius][iEnergyWeight][iJetAxis] == NULL){
            cout << "Error! The file " << fileName[iVeto][iJetRadius][iEnergyWeight][iJetAxis].Data() << " does not exist!" << endl;
            cout << "Maybe you forgot the data/ folder path?" << endl;
            cout << "Will not execute the code" << endl;
            return;
          }

          card[iVeto][iJetRadius][iEnergyWeight][iJetAxis]  = new EECCard(inputFile[iVeto][iJetRadius][iEnergyWeight][iJetAxis]);
        } // Jet axis loop
      } // Energy weight loop
    } // Jet radius loop
  } // Veto loop
  
  // Open the files and check that they exist
  TFile* inputFile[knSystems][knEnergyWeights][knJetAxes];
  EECCard* card[knSystems][knEnergyWeights][knJetAxes];
  for(int iSystem = 0; iSystem < knSystems; iSystem++){
    for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
      for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
    
        inputFile[iSystem][iEnergyWeight][iJetAxis] = TFile::Open(fileName[iSystem][iEnergyWeight][iJetAxis]);
    
        if(inputFile[iSystem][iEnergyWeight][iJetAxis] == NULL){
          cout << "Error! The file " << fileName[iSystem][iEnergyWeight][iJetAxis].Data() << " does not exist!" << endl;
          cout << "Maybe you forgot the data/ folder path?" << endl;
          cout << "Will not execute the code" << endl;
          return;
        }

        card[iSystem][iEnergyWeight][iJetAxis]  = new EECCard(inputFile[iSystem][iEnergyWeight][iJetAxis]);
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

  double maxJetDeltaAxis = 0.2;
  const int iMaxJetDeltaAxis = 100 * maxJetDeltaAxis;
  const double binTarget = 0.38; // For the trend, find this bin and use it to create the trend
  const bool normalizeDistributions = false; // Flag for normalizing distributions
  
  std::vector<std::pair<int,int>> comparedCentralityBin;
  //comparedCentralityBin.push_back(std::make_pair(0,10));

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

  // ====================================================
  //                Drawing configuration
  // ====================================================

  bool drawPp = true; // Include pp comparisons for single ratio histograms
  bool excludePbPb = comparedCentralityBin.size() == 0; // If no centrality bins are defined, exclude PbPb
  bool drawSingleRatios = false; // Draw the single ratio histograms
  bool drawDoubleRatios = false; // Draw the double ratio histograms
  bool drawTrends = false; // Draw edge effect trend as a function of DeltaR between E-scheme and WTA
  bool drawConsistencyCheck = false; // Compare summed DeltaR distributions to those projected without DeltaR cuts
  bool drawDeltaAxis = true; // Draw the DeltaR between E-scheme and WTA axis distributions themselves
  
  // Figure saving
  const bool saveFigures = false;  // Save figures
  TString saveComment = "_jetAxisDifference_signal";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures

  // Drawing configuration
  std::pair<double, double> ratioZoom = std::make_pair(0.8, 1.2);
  std::pair<double, double> eecZoom = std::make_pair(0.2, 30);
  std::pair<double, double> doubleRatioZoom = std::make_pair(0.8, 1.2);
  std::pair<double, double> trendZoom = std::make_pair(0.65,1.1);
  std::pair<double, double> consistencyCheckZoom = std::make_pair(0, 0.4);
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
      for(int iSystem = 0; iSystem < knSystems - excludePbPb; iSystem++){
        for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
          for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){

            // Sanity check for jet pT bins
            for(auto jetPtBin : comparedJetPtBin){
              if(card[iSystem][iEnergyWeight][iJetAxis]->FindBinIndexJetPtEEC(jetPtBin) < card[iSystem][iEnergyWeight][iJetAxis]->GetFirstUnfoldedJetPtBin() || card[iSystem][iEnergyWeight][iJetAxis]->FindBinIndexJetPtEEC(jetPtBin) > card[iSystem][iEnergyWeight][iJetAxis]->GetLastUnfoldedJetPtBin()){
                cout << "ERROR! Jet pT bin " << jetPtBin.first << "-" << jetPtBin.second << " does not exist in file " << fileName[iSystem][iEnergyWeight][iJetAxis].Data() << endl;
                cout << "Please only choose jet pT bins that are included in the input files." << endl;
                return;
              }
            }

            // Sanity check for track pT bins
            for(auto trackPtBin : comparedTrackPtBin){
              if(card[iSystem][iEnergyWeight][iJetAxis]->GetBinIndexTrackPtEEC(trackPtBin) < card[iSystem][iEnergyWeight][iJetAxis]->GetFirstUnfoldedTrackPtBin() || card[iSystem][iEnergyWeight][iJetAxis]->GetBinIndexTrackPtEEC(trackPtBin) > card[iSystem][iEnergyWeight][iJetAxis]->GetLastUnfoldedTrackPtBin()){
                cout << "ERROR! Track pT cut > " << trackPtBin << " GeV does not exist in file " << fileName[iSystem][iEnergyWeight][iJetAxis].Data() << endl;
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
  EECHistogramManager* histograms[knSystems][knEnergyWeights][knJetAxes];
  for(int iSystem = 0; iSystem < knSystems - excludePbPb; iSystem++){
    for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
      for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){

        histograms[iSystem][iEnergyWeight][iJetAxis] = new EECHistogramManager(inputFile[iSystem][iEnergyWeight][iJetAxis], card[iSystem][iEnergyWeight][iJetAxis]);

      } // File in set loop for histogram loading
    } // File set loop for histogram loading
  } // System loop for histogram loading

  // Energy-energy correlator histograms
  TH1D* hEnergyEnergyCorrelator[knEnergyWeights][knJetAxes][EECHistogramManager::kMaxJetDeltaAxisBins+2][nCentralityBins+1][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hEnergyEnergyCorrelatorRatio[knEnergyWeights][knJetAxes][EECHistogramManager::kMaxJetDeltaAxisBins+2][nCentralityBins+1][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hEnergyEnergyCorrelatorDoubleRatio[knEnergyWeights][knJetAxes][EECHistogramManager::kMaxJetDeltaAxisBins+2][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hEnergyEnergyConsistencyCheckRatio[knEnergyWeights][knJetAxes][EECHistogramManager::kMaxJetDeltaAxisBins+2][nCentralityBins+1][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hAxisRatioDeltaAxisTrend[knEnergyWeights][nCentralityBins+1][nJetPtBinsEEC][nTrackPtBinsEEC];

  TH1D* hJetAxisDeltaR[nCentralityBins+1][nJetPtBinsEEC];
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
    for(int iJetDeltaR = 0; iJetDeltaR < EECHistogramManager::kMaxJetDeltaAxisBins + 2; iJetDeltaR++){
      for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
        for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
          for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
            for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
              hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt] = NULL;
              hEnergyEnergyCorrelatorRatio[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt] = NULL;
              hEnergyEnergyCorrelatorDoubleRatio[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt] = NULL;
              hEnergyEnergyConsistencyCheckRatio[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt] = NULL;
            } // Centrality loop
            hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][nCentralityBins][iJetPt][iTrackPt] = NULL;
            hEnergyEnergyCorrelatorRatio[iEnergyWeight][iJetAxis][iJetDeltaR][nCentralityBins][iJetPt][iTrackPt] = NULL;
            hEnergyEnergyConsistencyCheckRatio[iEnergyWeight][iJetAxis][iJetDeltaR][nCentralityBins][iJetPt][iTrackPt] = NULL;
          } // Track pT loop
        } // Jet pT loop
      } // Jet axis loop
    } // DeltaR between E-scheme and WTA loop
    for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
          hAxisRatioDeltaAxisTrend[iEnergyWeight][iCentrality][iJetPt][iTrackPt] = NULL;
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Energy weight loop

  for(int iCentrality = 0; iCentrality <= nCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
      hJetAxisDeltaR[iCentrality][iJetPt] = NULL;
     } // Jet pT loop
  } // Centrality loop
  
  // Helper variables
  std::pair<double, double> drawingRange = std::make_pair(0.008, 0.39);
  double epsilon = 0.0001;
  int iTrackPt, iTrackPtReference;
  int iJetPt, iJetPtReference;
  int iCentrality, iCentralityReference;

  // Get the histograms from the histogram manager and normalize the signal histograms to one
  for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
    for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
      for(int iJetDeltaR = 0; iJetDeltaR < EECHistogramManager::kMaxJetDeltaAxisBins + 2; iJetDeltaR++){
        for(auto jetPtBin : comparedJetPtBin){
          for(auto trackPtBin : comparedTrackPtBin){

            // Find the proper binning and express it in term of the bins in the first file
            iJetPt = card[kPp][iEnergyWeight][iJetAxis]->FindBinIndexJetPtEEC(jetPtBin);
            iTrackPt = card[kPp][iEnergyWeight][iJetAxis]->GetBinIndexTrackPtEEC(trackPtBin);
            iJetPtReference = card[kPbPb][0][0]->FindBinIndexJetPtEEC(jetPtBin);
            iTrackPtReference = card[kPbPb][0][0]->GetBinIndexTrackPtEEC(trackPtBin);

            // For consistency check, sum all DeltaR between E-scheme and WTA axes bins together to make the total distribution
            if(iJetDeltaR == EECHistogramManager::kMaxJetDeltaAxisBins + 1){

              hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][nCentralityBins][iJetPtReference][iTrackPtReference] = (TH1D*) hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][0][nCentralityBins][iJetPtReference][iTrackPtReference]->Clone(Form("summedDistribution%d%d%d%d", iEnergyWeight, iJetAxis, iJetPtReference, iTrackPtReference));

              for(int jJetDeltaR = 1; jJetDeltaR < EECHistogramManager::kMaxJetDeltaAxisBins; jJetDeltaR++){

                hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][nCentralityBins][iJetPtReference][iTrackPtReference]->Add(hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][jJetDeltaR][nCentralityBins][iJetPtReference][iTrackPtReference]);

              } // Loop over all DeltaR between E-scheme and WTA axes bins

            } else {

              // Load the selected energy-energy correlator histograms for pp
              if(drawnEnergyEnergyCorrelator == EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels){
                hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][nCentralityBins][iJetPtReference][iTrackPtReference] = histograms[kPp][iEnergyWeight][iJetAxis]->GetHistogramEnergyEnergyCorrelatorJetDeltaAxis(EECHistogramManager::kEnergyEnergyCorrelator, iJetDeltaR, 0, iJetPt, iTrackPt, EECHistograms::kLeadingParticleTypes, pairingType);
              } else {
                // TODO: For processed histograms, I currently do not have the getter setup with JetDeltaR binning 
                hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][nCentralityBins][iJetPtReference][iTrackPtReference] = histograms[kPp][iEnergyWeight][iJetAxis]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPt, iTrackPt, drawnEnergyEnergyCorrelator);
              } 

            } // Regular histogram reading

            // Do the same for all the selected centrality bins for PbPb
            iJetPt = card[kPbPb][iEnergyWeight][iJetAxis]->FindBinIndexJetPtEEC(jetPtBin);
            iTrackPt = card[kPbPb][iEnergyWeight][iJetAxis]->GetBinIndexTrackPtEEC(trackPtBin);
            for(auto centralityBin : comparedCentralityBin){
              iCentrality = card[kPbPb][iEnergyWeight][iJetAxis]->FindBinIndexCentrality(centralityBin);
              iCentralityReference = card[kPbPb][0][0]->FindBinIndexCentrality(centralityBin);

              // For consistency check, sum all DeltaR between E-scheme and WTA axes bins together to make the total distribution
              if(iJetDeltaR == EECHistogramManager::kMaxJetDeltaAxisBins + 1){

                hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentralityReference][iJetPtReference][iTrackPtReference] = (TH1D*) hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][0][iCentralityReference][iJetPtReference][iTrackPtReference]->Clone(Form("summedDistribution%d%d%d%d%d", iEnergyWeight, iJetAxis, iCentralityReference, iJetPtReference, iTrackPtReference));

                for(int jJetDeltaR = 1; jJetDeltaR < EECHistogramManager::kMaxJetDeltaAxisBins; jJetDeltaR++){

                  hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentralityReference][iJetPtReference][iTrackPtReference]->Add(hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][jJetDeltaR][iCentralityReference][iJetPtReference][iTrackPtReference]);

                } // Loop over all DeltaR between E-scheme and WTA axes bins
              } else {

                // Load the selected energy-energy correlator histograms for PbPb
                if(drawnEnergyEnergyCorrelator == EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels){
                  hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPtReference][iTrackPtReference] = histograms[kPbPb][iEnergyWeight][iJetAxis]->GetHistogramEnergyEnergyCorrelatorJetDeltaAxis(EECHistogramManager::kEnergyEnergyCorrelator, iJetDeltaR, iCentrality, iJetPt, iTrackPt, EECHistograms::kLeadingParticleTypes, pairingType, subevent);
                } else {
                  // TODO: For processed histograms, I currently do not have the getter setup with JetDeltaR binning 
                  hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPtReference][iTrackPtReference] = histograms[kPbPb][iEnergyWeight][iJetAxis]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, drawnEnergyEnergyCorrelator);
                }

              } // Regular histogram reading
            
            } // Centrality loop
          } // Track pT loop
        } // Jet pT loop
      } // DeltaR between E-scheme and WTA axes loop
    } // Jet axis type loop
  } // Energy weight loop

  // Load the DeltaR between E-scheme and WTA aces histograms
  for(auto jetPtBin : comparedJetPtBin){
    iJetPt = card[kPp][0][0]->FindBinIndexJetPtEEC(jetPtBin);
    iJetPtReference = card[kPbPb][0][0]->FindBinIndexJetPtEEC(jetPtBin);

    hJetAxisDeltaR[nCentralityBins][iJetPtReference] = histograms[kPp][0][0]->GetHistogramJetDeltaAxis(0, iJetPt);
    
    for(auto centralityBin : comparedCentralityBin){
      iCentrality = card[kPbPb][0][0]->FindBinIndexCentrality(centralityBin);
      hJetAxisDeltaR[iCentrality][iJetPt] = histograms[kPbPb][0][0]->GetHistogramJetDeltaAxis(iCentrality, iJetPt);
    } // Centrality loop
  } // Jet pT loop

  // ====================================================
  //               Normalize distributions
  // ====================================================

  int lowNormalizationBin, highNormalizationBin;
  if(normalizeDistributions){
    for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
      for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
        for(int iJetDeltaR = 0; iJetDeltaR < EECHistogramManager::kMaxJetDeltaAxisBins + 2; iJetDeltaR++){
          for(auto jetPtBin : comparedJetPtBin){
            iJetPt = card[kPbPb][0][0]->FindBinIndexJetPtEEC(jetPtBin);
            for(auto trackPtBin : comparedTrackPtBin){
              iTrackPt = card[kPbPb][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
            
              // Normalize the distributions to one in the drawingRange
              lowNormalizationBin = hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][nCentralityBins][iJetPt][iTrackPt]->GetXaxis()->FindBin(drawingRange.first + epsilon);
              highNormalizationBin = hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][nCentralityBins][iJetPt][iTrackPt]->GetXaxis()->FindBin(drawingRange.second - epsilon);

             hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][nCentralityBins][iJetPt][iTrackPt]->Scale(1 / hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][nCentralityBins][iJetPt][iTrackPt]->Integral(lowNormalizationBin, highNormalizationBin, "width"));

              // Do the same for all the selected centrality bins for PbPb
              for(auto centralityBin : comparedCentralityBin){

                // Normalize the distributions to one in the drawingRange
                lowNormalizationBin = hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(drawingRange.first + epsilon);
                highNormalizationBin = hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(drawingRange.second - epsilon);

                hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->Scale(1 / hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->Integral(lowNormalizationBin, highNormalizationBin, "width"));
            
              } // Centrality loop
            } // Track pT loop
          } // Jet pT loop
        } // DeltaR between E-scheme and WTA axes loop
      } // Jet axis type loop
    } // Energy weight loop
  } // Normalizinf distributions

  // ====================================================
  //                 Calculating ratios
  // ====================================================

  // Make a centrality vector that contains also pp for easier looping over all results
  // For drawing file set comparisons, add flag for pp drawing to the centrality bin array
  std::vector<std::pair<int,int>> distributionCentralityBin;
  if(drawPp) distributionCentralityBin.push_back(std::make_pair(-1,-1));
  for(auto centralityBin : comparedCentralityBin){
    distributionCentralityBin.push_back(centralityBin);
  }

  // After all the histograms have been read, calculate the WTA to Escheme ratios
  for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
    for(int iJetDeltaR = 0; iJetDeltaR < EECHistogramManager::kMaxJetDeltaAxisBins + 2; iJetDeltaR++){
      for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
        for(auto jetPtBin : comparedJetPtBin){
          iJetPt = card[kPbPb][0][0]->FindBinIndexJetPtEEC(jetPtBin);
          for(auto trackPtBin : comparedTrackPtBin){
            iTrackPt = card[kPbPb][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
            for(auto centralityBin : distributionCentralityBin){

              // The pp centrality index is special, and cannot be read from the EECCard
              if(centralityBin.first < 0){
                iCentrality = nCentralityBins;
              } else {
                iCentrality = card[kPbPb][0][0]->FindBinIndexCentrality(centralityBin);
              }

                hEnergyEnergyCorrelatorRatio[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->Clone(Form("eecRatio%d%d%d%d%d%d", iEnergyWeight, iJetAxis, iJetDeltaR, iCentrality, iJetPt, iTrackPt));

              if(iJetDeltaR == EECHistogramManager::kMaxJetDeltaAxisBins + 1){
                // This is the ratio for consistency check. Instead of axis ratio, do a ratio between individual DeltaR between E-scheme and WTA bins axes summed together, and the one where these were never separated

                hEnergyEnergyCorrelatorRatio[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->Divide(hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][EECHistogramManager::kMaxJetDeltaAxisBins][iCentrality][iJetPt][iTrackPt]);
 
              } else {
                // In the default case, take a ratio with respect to E-scheme result
                
                hEnergyEnergyCorrelatorRatio[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->Divide(hEnergyEnergyCorrelator[iEnergyWeight][kEscheme][iJetDeltaR][iCentrality][iJetPt][iTrackPt]);
              }

            } // Centrality loop
          } // Track pT loop
        } // Jet pT loop
      } // Jet axis type loop
    } // DeltaR between E-scheme and WTA axes loop
  } // Energy weight loop

  // We also want to calculate double ratios of PbPb to pp ratios of the single ratios
  if(drawDoubleRatios){
    for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
      for(int iJetDeltaR = 0; iJetDeltaR <= EECHistogramManager::kMaxJetDeltaAxisBins; iJetDeltaR++){
        if(iJetDeltaR == iMaxJetDeltaAxis) iJetDeltaR = EECHistogramManager::kMaxJetDeltaAxisBins;
        for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
          for(auto jetPtBin : comparedJetPtBin){
            iJetPt = card[kPbPb][0][0]->FindBinIndexJetPtEEC(jetPtBin);
            for(auto trackPtBin : comparedTrackPtBin){
              iTrackPt = card[kPbPb][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
              for(auto centralityBin: comparedCentralityBin){
                iCentrality = card[kPbPb][0][0]->FindBinIndexCentrality(centralityBin);

                hEnergyEnergyCorrelatorDoubleRatio[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelatorRatio[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->Clone(Form("eecDoubleRatio%d%d%d%d%d%d", iEnergyWeight, iJetAxis, iJetDeltaR, iCentrality, iJetPt, iTrackPt));
                hEnergyEnergyCorrelatorDoubleRatio[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->Divide(hEnergyEnergyCorrelatorRatio[iEnergyWeight][iJetAxis][iJetDeltaR][nCentralityBins][iJetPt][iTrackPt]);

              } // Centrality loop
            } // Track pT loop
          } // Jet pT loop
        } // Jet axis loop
      } // DeltaR between E-scheme and WTA axes loop
    } // Energy weight loop
  }

  // For consistency check, calculate the ratio of each individual DeltaR between E-scheme and WTA axes bin to the total integrated one
  // After all the histograms have been read, calculate the WTA to Escheme ratios
  for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
    for(int iJetDeltaR = 0; iJetDeltaR < EECHistogramManager::kMaxJetDeltaAxisBins + 2; iJetDeltaR++){
      for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
        for(auto jetPtBin : comparedJetPtBin){
          iJetPt = card[kPbPb][0][0]->FindBinIndexJetPtEEC(jetPtBin);
          for(auto trackPtBin : comparedTrackPtBin){
            iTrackPt = card[kPbPb][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
            for(auto centralityBin : distributionCentralityBin){

              // The pp centrality index is special, and cannot be read from the EECCard
              if(centralityBin.first < 0){
                iCentrality = nCentralityBins;
              } else {
                iCentrality = card[kPbPb][0][0]->FindBinIndexCentrality(centralityBin);
              }

                hEnergyEnergyConsistencyCheckRatio[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->Clone(Form("eecConsistencyCheckRatio%d%d%d%d%d%d", iEnergyWeight, iJetAxis, iJetDeltaR, iCentrality, iJetPt, iTrackPt));
                hEnergyEnergyConsistencyCheckRatio[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->Divide(hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][EECHistogramManager::kMaxJetDeltaAxisBins][iCentrality][iJetPt][iTrackPt]);

            } // Centrality loop
          } // Track pT loop
        } // Jet pT loop
      } // Jet axis type loop
    } // DeltaR between E-scheme and WTA axes loop
  } // Energy weight loop

  // Once all the ratios have been calculated, calculate the WTA/E-scheme ratio trend as a function of DeltaR between the two axes in a specified DeltaR bin
  int targetBin;
  double targetContent, targetError;
  std::pair<double, double> trendBinBoundaries = std::make_pair(0,0);
  for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
    for(auto jetPtBin : comparedJetPtBin){
      iJetPt = card[kPbPb][0][0]->FindBinIndexJetPtEEC(jetPtBin);
        for(auto trackPtBin : comparedTrackPtBin){
          iTrackPt = card[kPbPb][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
          for(auto centralityBin : distributionCentralityBin){

            if(centralityBin.first == -1){
              iCentrality = nCentralityBins;
            } else {
              iCentrality = card[kPbPb][0][0]->FindBinIndexCentrality(centralityBin);
            }

            // Create a new histogram for the trend
            hAxisRatioDeltaAxisTrend[iEnergyWeight][iCentrality][iJetPt][iTrackPt] = new TH1D(Form("greatRatio%d%d%d%d", iEnergyWeight, iJetPt, iTrackPt, iCentrality), Form("greatRatio%d%d%d%d", iEnergyWeight, iJetPt, iTrackPt, iCentrality), iMaxJetDeltaAxis, 0, maxJetDeltaAxis);

            // Find a great bin to calculate the trend from
            targetBin = hEnergyEnergyCorrelatorRatio[iEnergyWeight][kWTA][0][iCentrality][iJetPt][iTrackPt]->FindBin(binTarget);
            trendBinBoundaries.first = hEnergyEnergyCorrelatorRatio[iEnergyWeight][kWTA][0][iCentrality][iJetPt][iTrackPt]->GetBinLowEdge(targetBin);
            trendBinBoundaries.second = trendBinBoundaries.first + hEnergyEnergyCorrelatorRatio[iEnergyWeight][kWTA][0][iCentrality][iJetPt][iTrackPt]->GetBinWidth(targetBin);

            for(int iJetDeltaR = 0; iJetDeltaR < iMaxJetDeltaAxis; iJetDeltaR++){
              targetContent = hEnergyEnergyCorrelatorRatio[iEnergyWeight][kWTA][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->GetBinContent(targetBin);
              targetError = hEnergyEnergyCorrelatorRatio[iEnergyWeight][kWTA][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->GetBinError(targetBin);

              hAxisRatioDeltaAxisTrend[iEnergyWeight][iCentrality][iJetPt][iTrackPt]->SetBinContent(iJetDeltaR+1, targetContent);
              hAxisRatioDeltaAxisTrend[iEnergyWeight][iCentrality][iJetPt][iTrackPt]->SetBinError(iJetDeltaR+1, targetError);
            }
        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // energy weight loop

  
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
  TString compactDeltaRString = "";
  TString deltaRString = "";
  TString legendString;
  int markerStyle[] = {kFullSquare, kFullCircle, kOpenSquare, kOpenCircle, kFullCross};
  int color[] = {kBlack, kBlue, kGreen+3, kRed, kMagenta,kCyan,kOrange,kViolet+3,kPink-7,kSpring+3,kAzure-7};

  // Binning vectors
  std::vector<int> currentJetPtIndices;
  std::vector<int> currentTrackPtIndices;
  double minimumCandidate, maximumCandidate;
  TString systemForLegend = "";
  TString consistencyString;
  TString consistencyRatio = "#frac{Summed}{Whole}";

  // Legend and line
  TLegend* legend;
  TLine* oneLine = new TLine(drawingRange.first, 1, drawingRange.second, 1);
  oneLine->SetLineStyle(2);

  // Draw consistency check to ensure adding together all DeltaR between E-scheme and WTA axis bins we recover the full distribution
  if(drawConsistencyCheck){

    for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
      for(int iJetDeltaR = 0; iJetDeltaR < EECHistogramManager::kMaxJetDeltaAxisBins + 2; iJetDeltaR++){
        if(iJetDeltaR == iMaxJetDeltaAxis) iJetDeltaR = EECHistogramManager::kMaxJetDeltaAxisBins+1;
        if(iJetDeltaR == EECHistogramManager::kMaxJetDeltaAxisBins+1){
          compactDeltaRString = "_unityCheck";
          deltaRString = "Summed";
        } else {
          compactDeltaRString = Form("_DA%d-%d", iJetDeltaR, iJetDeltaR+1);
          deltaRString = Form("%.2f < #Deltar_{j} < %.2f", iJetDeltaR/100.0, iJetDeltaR/100.0+0.01);
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
              legend->AddEntry((TObject*) 0, jetAxisDescription[iJetAxis].Data(), "");
              if(centralityBin.first != -1) legend->AddEntry((TObject*) 0, Form("Centrality: %d-%d%%", centralityBin.first, centralityBin.second), "");
              legend->AddEntry((TObject*) 0, Form("%.0f < jet p_{T} < %.0f GeV", jetPtBin.first, jetPtBin.second), "");
              legend->AddEntry((TObject*) 0, Form("p_{T}^{ch} > %.1f GeV", trackPtBin), "");

              // Setup centrality, jet pT and track pT strings
              compactJetPtString = Form("_J=%.0f-%.0f", jetPtBin.first, jetPtBin.second);
              compactTrackPtString = Form("_T>%.1f",trackPtBin);
              compactTrackPtString.ReplaceAll(".","v");
              compactCentralityString = centralityBin.first == -1 ? "_pp" : Form("_C=%d-%d", centralityBin.first, centralityBin.second);


              // Set drawing style for all histograms
              for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
                hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iEnergyWeight]);
                hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iEnergyWeight]);
                hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iEnergyWeight]);
                hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][EECHistogramManager::kMaxJetDeltaAxisBins][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iEnergyWeight + knEnergyWeights]);
                hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][EECHistogramManager::kMaxJetDeltaAxisBins][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iEnergyWeight + knEnergyWeights]);
                hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][EECHistogramManager::kMaxJetDeltaAxisBins][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iEnergyWeight + knEnergyWeights]); 
                hEnergyEnergyConsistencyCheckRatio[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iEnergyWeight]);
                hEnergyEnergyConsistencyCheckRatio[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iEnergyWeight]);
                hEnergyEnergyConsistencyCheckRatio[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iEnergyWeight]);
              } // Energy weight loop

              // Automatic zooming for the drawn histograms
              if(automaticZoom){
                hEnergyEnergyCorrelator[kNominalWeight][iJetAxis][EECHistogramManager::kMaxJetDeltaAxisBins][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
                hEnergyEnergyCorrelator[kSquaredWeight][iJetAxis][EECHistogramManager::kMaxJetDeltaAxisBins][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
                eecZoom.first = hEnergyEnergyCorrelator[kNominalWeight][iJetAxis][EECHistogramManager::kMaxJetDeltaAxisBins][iCentrality][iJetPt][iTrackPt]->GetMinimum();
                eecZoom.second = hEnergyEnergyCorrelator[kSquaredWeight][iJetAxis][EECHistogramManager::kMaxJetDeltaAxisBins][iCentrality][iJetPt][iTrackPt]->GetMaximum();
                for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
                  hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
                  minimumCandidate = hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->GetMinimum();
                  maximumCandidate = hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->GetMaximum();
                  if(minimumCandidate < eecZoom.first) eecZoom.first = minimumCandidate;
                  if(maximumCandidate > eecZoom.second) eecZoom.second = maximumCandidate;
                } // Energy weight loop
                eecZoom.first = eecZoom.first / 2.0;
                eecZoom.second = eecZoom.second * 2.0;
              }

              // Set the x- and y-axis drawing ranges
              hEnergyEnergyCorrelator[0][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
              hEnergyEnergyCorrelator[0][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(eecZoom.first, eecZoom.second);

              // Draw the histograms to the upper canvas
              drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelator[0][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt], "#Deltar", Form("EEC (%s)", energyCorrelatorType[drawnEnergyEnergyCorrelator].Data()), " ", "p");

              for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
                if(iEnergyWeight > 0) hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->Draw("same,p");
                hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][EECHistogramManager::kMaxJetDeltaAxisBins][iCentrality][iJetPt][iTrackPt]->Draw("same,p");
              } // Energy weight loop

              // Add legends for drawn histograms
              for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
                legend->AddEntry(hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt], Form("%s, %s", energyWeightDescription[iEnergyWeight].Data(), deltaRString.Data()), "p");
                legend->AddEntry(hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][EECHistogramManager::kMaxJetDeltaAxisBins][iCentrality][iJetPt][iTrackPt], Form("%s, #Deltar_{E-scheme}^{WTA} integrated", energyWeightDescription[iEnergyWeight].Data()), "p");
              } // Energy weight loop
  
              // Draw the legends to the upper pad
              legend->Draw();
          
              // Linear scale for the ratio
              drawer->SetLogY(false);
          
              // Set the drawing ranges
              hEnergyEnergyConsistencyCheckRatio[0][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
              if(iJetDeltaR < EECHistogramManager::kMaxJetDeltaAxisBins){
                hEnergyEnergyConsistencyCheckRatio[0][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(consistencyCheckZoom.first, consistencyCheckZoom.second);
              } else {
                hEnergyEnergyConsistencyCheckRatio[0][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.95, 1.05);
              }

              // Draw the histograms
              drawer->SetGridY(true);

              drawer->DrawHistogramToLowerPad(hEnergyEnergyConsistencyCheckRatio[0][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt], "#Deltar", Form("#frac{%s}{Integrated}", deltaRString.Data()), " ", "p");
              for(int iEnergyWeight = 1; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
                hEnergyEnergyConsistencyCheckRatio[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->Draw("same");
              } // Energy weight loop
              drawer->SetGridY(false);
          
              // Save the figures to a file
              if(saveFigures){
                gPad->GetCanvas()->SaveAs(Form("figures/eecConsistencyCheck%s%s%s%s%s%s.%s", saveComment.Data(), compactDeltaRString.Data(), jetAxisDescription[iJetAxis].Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
              }
            } // Centrality looop
          } // Track pT loop
        } // Jet pT loop
      } // DeltaR between E-scheme and WTA axes loop
    } // Jet axis loop
  }

  // Draw the single ratio histograms
  if(drawSingleRatios){

    for(int iJetDeltaR = 0; iJetDeltaR <= EECHistogramManager::kMaxJetDeltaAxisBins; iJetDeltaR++){
      if(iJetDeltaR == iMaxJetDeltaAxis) iJetDeltaR = EECHistogramManager::kMaxJetDeltaAxisBins;
      if(iJetDeltaR == EECHistogramManager::kMaxJetDeltaAxisBins) {
        compactDeltaRString = "";
      } else {
        compactDeltaRString = Form("_DA%d-%d", iJetDeltaR, iJetDeltaR+1);
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
            if(iJetDeltaR < EECHistogramManager::kMaxJetDeltaAxisBins) legend->AddEntry((TObject*) 0, Form("%.2f < #Deltar_{E-scheme}^{WTA} < %.2f", iJetDeltaR/100.0, iJetDeltaR/100.0+0.01), "");
            if(centralityBin.first != -1) legend->AddEntry((TObject*) 0, Form("Centrality: %d-%d%%", centralityBin.first, centralityBin.second), "");
            legend->AddEntry((TObject*) 0, Form("%.0f < jet p_{T} < %.0f GeV", jetPtBin.first, jetPtBin.second), "");
            legend->AddEntry((TObject*) 0, Form("p_{T}^{ch} > %.1f GeV", trackPtBin), "");

            // Setup centrality, jet pT and track pT strings
            compactJetPtString = Form("_J=%.0f-%.0f", jetPtBin.first, jetPtBin.second);
            compactTrackPtString = Form("_T>%.1f",trackPtBin);
            compactTrackPtString.ReplaceAll(".","v");
            compactCentralityString = centralityBin.first == -1 ? "_pp" : Form("_C=%d-%d", centralityBin.first, centralityBin.second);


            // Set drawing style for all histograms
            for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
              for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
                hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iJetAxis + knJetAxes*iEnergyWeight]);
                hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iJetAxis + knJetAxes*iEnergyWeight]);
                hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iJetAxis + knJetAxes*iEnergyWeight]); 
                hEnergyEnergyCorrelatorRatio[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iJetAxis + knJetAxes*iEnergyWeight]);
                hEnergyEnergyCorrelatorRatio[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iJetAxis + knJetAxes*iEnergyWeight]);
                hEnergyEnergyCorrelatorRatio[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iJetAxis + knJetAxes*iEnergyWeight]);
              } // File in set loop
            } // File set loop

            // Automatic zooming for the drawn histograms
            if(automaticZoom){
              eecZoom.first = 10000;
              eecZoom.second = 0;
              for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
                for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
                  hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
                  minimumCandidate = hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->GetMinimum();
                  maximumCandidate = hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->GetMaximum();
                  if(minimumCandidate < eecZoom.first) eecZoom.first = minimumCandidate;
                  if(maximumCandidate > eecZoom.second) eecZoom.second = maximumCandidate;
                } // Jet axis loop
              } // Energy weight loop
              eecZoom.first = eecZoom.first / 2.0;
              eecZoom.second = eecZoom.second * 2.0;
            }

            // Set the x- and y-axis drawing ranges
            hEnergyEnergyCorrelator[0][0][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
            hEnergyEnergyCorrelator[0][0][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(eecZoom.first, eecZoom.second);

            // Draw the histograms to the upper canvas
            drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelator[0][0][iJetDeltaR][iCentrality][iJetPt][iTrackPt], "#Deltar", Form("EEC (%s)", energyCorrelatorType[drawnEnergyEnergyCorrelator].Data()), " ", "p");

            for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
              for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
                if(iEnergyWeight == 0 && iJetAxis == 0) continue;
                hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->Draw("same,p");
              } // Jet axis loop
            } // Energy weight loop

            // Add legends for drawn histograms
            for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
              for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
                legend->AddEntry(hEnergyEnergyCorrelator[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt], Form("%s  %s", energyWeightDescription[iEnergyWeight].Data(), jetAxisDescription[iJetAxis].Data()), "p");
              } // Jet axis loop
            } // Energy weight loop
  
            // Draw the legends to the upper pad
            legend->Draw();
          
            // Linear scale for the ratio
            drawer->SetLogY(false);
          
            // Set the drawing ranges
            hEnergyEnergyCorrelatorRatio[0][1][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
            hEnergyEnergyCorrelatorRatio[0][1][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

            // Draw the histograms
            drawer->SetGridY(true);

            drawer->DrawHistogramToLowerPad(hEnergyEnergyCorrelatorRatio[0][1][iJetDeltaR][iCentrality][iJetPt][iTrackPt], "#Deltar", ratioDescription, " ", "p");
            for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
              for(int iJetAxis = 1; iJetAxis < knJetAxes; iJetAxis++){
                if(iEnergyWeight == 0 && iJetAxis == 1) continue;
                hEnergyEnergyCorrelatorRatio[iEnergyWeight][iJetAxis][iJetDeltaR][iCentrality][iJetPt][iTrackPt]->Draw("same");
              } // Jet axis loop
            } // Energy weight loop
            drawer->SetGridY(false);
          
            // Save the figures to a file
            if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/eecComparison%s%s%s%s%s.%s", saveComment.Data(), compactDeltaRString.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
            }
          } // Centrality looop
        } // Track pT loop
      } // Jet pT loop
    } // DeltaR between E-scheme and WTA axes loop
  } // Drawing single ratios

  

  if(drawDoubleRatios){

    TString ratioDenominator;

    // Name for the double ratio axis
    TString doubleRatioAxis = Form("EEC (%s)  #frac{%s (%s / %s)}{%s (%s / %s)}", energyCorrelatorType[drawnEnergyEnergyCorrelator].Data(), card[kPbPb][0][0]->GetAlternativeDataType(false).Data(), jetAxisDescription[kWTA].Data(), jetAxisDescription[kEscheme].Data(), card[kPp][0][0]->GetAlternativeDataType(false).Data(), jetAxisDescription[kWTA].Data(), jetAxisDescription[kEscheme].Data());
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
          for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
            for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
              hEnergyEnergyCorrelatorDoubleRatio[iEnergyWeight][iJetAxis][EECHistogramManager::kMaxJetDeltaAxisBins][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iJetAxis + knJetAxes*iEnergyWeight]);
              hEnergyEnergyCorrelatorDoubleRatio[iEnergyWeight][iJetAxis][EECHistogramManager::kMaxJetDeltaAxisBins][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iJetAxis + knJetAxes*iEnergyWeight]);
              hEnergyEnergyCorrelatorDoubleRatio[iEnergyWeight][iJetAxis][EECHistogramManager::kMaxJetDeltaAxisBins][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iJetAxis + knJetAxes*iEnergyWeight]);
            } // File in set loop
          } // File set loop

          // Set the x- and y-axis drawing ranges
          hEnergyEnergyCorrelatorDoubleRatio[1][1][EECHistogramManager::kMaxJetDeltaAxisBins][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
          hEnergyEnergyCorrelatorDoubleRatio[1][1][EECHistogramManager::kMaxJetDeltaAxisBins][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(doubleRatioZoom.first, doubleRatioZoom.second);

          // Draw the histograms
          drawer->DrawHistogram(hEnergyEnergyCorrelatorDoubleRatio[1][1][EECHistogramManager::kMaxJetDeltaAxisBins][iCentrality][iJetPt][iTrackPt], "#Deltar", doubleRatioAxis.Data(), " ", "p");

          for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
            for(int iJetAxis = 1; iJetAxis < knJetAxes; iJetAxis++){
              if(iEnergyWeight == 1 && iJetAxis == 1) continue;
              hEnergyEnergyCorrelatorDoubleRatio[iEnergyWeight][iJetAxis][EECHistogramManager::kMaxJetDeltaAxisBins][iCentrality][iJetPt][iTrackPt]->Draw("same,p");
            } // File in set loop
          } // File set loop

          // Draw a line to one
          oneLine->Draw();

          // Add legends for drawn histograms
          for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
            for(int iJetAxis = 1; iJetAxis < knJetAxes; iJetAxis++){
              legend->AddEntry(hEnergyEnergyCorrelatorDoubleRatio[iEnergyWeight][iJetAxis][EECHistogramManager::kMaxJetDeltaAxisBins][iCentrality][iJetPt][iTrackPt], energyWeightDescription[iEnergyWeight].Data(), "p");
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

  // Draw the single ratio histograms
  if(drawTrends){

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
        for(auto centralityBin : distributionCentralityBin){
          if(centralityBin.first == -1){
            iCentrality = nCentralityBins;
            systemForLegend = card[kPp][0][0]->GetAlternativeDataType(false);
          } else {
            iCentrality = card[kPbPb][0][0]->FindBinIndexCentrality(centralityBin);
            systemForLegend = card[kPbPb][0][0]->GetAlternativeDataType(false);
          }
          
          // Linear scale for ratio
          drawer->SetLogX(false);
          drawer->SetLogY(false);

          // Create the legend and add binning information to it
          legend = new TLegend(0.18,0.14,0.45,0.58);
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
          for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
            hAxisRatioDeltaAxisTrend[iEnergyWeight][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iEnergyWeight]);
            hAxisRatioDeltaAxisTrend[iEnergyWeight][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iEnergyWeight]);
            hAxisRatioDeltaAxisTrend[iEnergyWeight][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iEnergyWeight]); 
          } // Energy weight loop


          // Set the x- and y-axis drawing ranges
          hAxisRatioDeltaAxisTrend[0][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(trendZoom.first, trendZoom.second);

          // Draw the histograms to the upper canvas
          drawer->DrawHistogram(hAxisRatioDeltaAxisTrend[0][iCentrality][iJetPt][iTrackPt], "#Deltar_{E-scheme}^{WTA}", Form("EEC (%s) #frac{WTA}{E-scheme} (%.3f < #Deltar < %.3f)", energyCorrelatorType[drawnEnergyEnergyCorrelator].Data(), trendBinBoundaries.first, trendBinBoundaries.second), " ", "p");

          for(int iEnergyWeight = 1; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
            hAxisRatioDeltaAxisTrend[iEnergyWeight][iCentrality][iJetPt][iTrackPt]->Draw("same,p");
          } // Energy weight loop

          // Add legends for drawn histograms
          for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
            legend->AddEntry(hAxisRatioDeltaAxisTrend[iEnergyWeight][iCentrality][iJetPt][iTrackPt], Form("%s", energyWeightDescription[iEnergyWeight].Data()), "p");
          } // Energy weight loop
  
          // Draw the legends to the upper pad
          legend->Draw();
          
          // Save the figures to a file
          if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/eecRatioTrend%s%s%s%s%s.%s", saveComment.Data(), compactDeltaRString.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
          }
        } // Centrality looop
      } // Track pT loop
    } // Jet pT loop
  } // Drawing trends

  if(drawDeltaAxis){

    // Set a good drawing style for single canvas plots
    drawer->SetRelativeCanvasSize(1.1,1.1);
    drawer->SetLeftMargin(0.14);
    drawer->SetTopMargin(0.07);
    drawer->SetTitleOffsetY(1.7);
    drawer->SetTitleOffsetX(1.0);

    for(auto jetPtBin : comparedJetPtBin){
      iJetPt = card[kPbPb][0][0]->FindBinIndexJetPtEEC(jetPtBin);

      for(auto centralityBin : distributionCentralityBin){
        if(centralityBin.first == -1){
          iCentrality = nCentralityBins;
          systemForLegend = card[kPp][0][0]->GetAlternativeDataType(false);
        } else {
          iCentrality = card[kPbPb][0][0]->FindBinIndexCentrality(centralityBin);
          systemForLegend = card[kPbPb][0][0]->GetAlternativeDataType(false);
        }

        // Linear axis scales
        drawer->SetLogX(false);
        drawer->SetLogY(false);

        // Create the legend and add binning information to it
        legend = new TLegend(0.35,0.5,0.8,0.65);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);

        if(systemForLegend != "") legend->AddEntry((TObject*) 0, systemForLegend.Data(), "");
        if(centralityBin.first != -1) legend->AddEntry((TObject*) 0, Form("Centrality: %d-%d%%", centralityBin.first, centralityBin.second), "");
        legend->AddEntry((TObject*) 0, Form("%.0f < jet p_{T} < %.0f GeV", jetPtBin.first, jetPtBin.second), "");

        // Setup centrality, jet pT and track pT strings
        compactJetPtString = Form("_J=%.0f-%.0f", jetPtBin.first, jetPtBin.second);
        compactCentralityString = centralityBin.first == -1 ? "_pp" : Form("_C=%d-%d", centralityBin.first, centralityBin.second);

        // Set drawing style for the histograms
        hJetAxisDeltaR[iCentrality][iJetPt]->SetMarkerStyle(markerStyle[0]);
        hJetAxisDeltaR[iCentrality][iJetPt]->SetMarkerColor(color[0]);
        hJetAxisDeltaR[iCentrality][iJetPt]->SetLineColor(color[0]); 

        // Draw the histograms to the a canvas
        drawer->DrawHistogram(hJetAxisDeltaR[iCentrality][iJetPt], "#Deltar_{E-scheme}^{WTA}", "Counts", " ", "p");  

        // Draw the legend to the upper pad
        legend->Draw();   

        // Save the figures to a file
        if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/jetDeltaAxis%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), figureFormat));
        }   

      } // Centrality loop
    } // Jet pT loop

  }

}
