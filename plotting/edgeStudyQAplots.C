#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"
#include "../src/EECHistograms.h"

/*
 * Macro for comparing histograms from two different sets of files
 */
void edgeStudyQAplots(){
  
  // Enumeration for different jet axes
  enum enumJetAxis{kEscheme, kWTA, knJetAxes};

  // Enumeration for energy weight
  enum enumEnergyWeight{kNominalWeight, kSquaredWeight, knEnergyWeights};

  // Enumeration for jet radius
  enum enumJetRatios{kR0p4, kR0p8, knJetRadii};

  // Enumeration for vetoing jets in the region 0.4 < DeltaR < 0.8
  enum enumJetVeto{kNoVeto, kJetVeto, knVetoTypes};

  // Enumeration for type of simulation
  enum enumSimulation{kPythia, kHerwig, knSimulations};

  // Possibility to veto other jets within the extended radius around the jet axis
  const int vetoCloseJets = kNoVeto;

  // Select if you want to draw Pythia or Herwig QA plots
  const int iSimulation = kHerwig;

  // Define the file names for the studied files
  TString fileName[knSimulations][knVetoTypes][knJetRadii][knEnergyWeights][knJetAxes];

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

  // Naming for simulations
  TString simulationName[knSimulations];
  simulationName[kPythia] = "Pythia8";
  simulationName[kHerwig] = "Herwig7";

  // Figure saving names for simulations
  TString simulationSaveName[knSimulations];
  simulationSaveName[kPythia] = "_pythia";
  simulationSaveName[kHerwig] = "_herwig";

  // Naming for leading particle types
  TString leadingParticleName[EECHistograms::kLeadingParticleTypes+1];
  leadingParticleName[EECHistograms::kNotLeadingParticle] = "Pairs without leading particle";
  leadingParticleName[EECHistograms::kLeadingParticle] = "Pairs with leading particle";
  leadingParticleName[EECHistograms::kLeadingParticleTypes] = "All pairs";
  

  // Define files for each of these combinations.

  // ======================//
  //    Files for Pythia   //
  // ======================//

  // Files for R = 0.4, no other jet veto
  fileName[kPythia][kNoVeto][kR0p4][kNominalWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_nominalEnergyWeight_jetDeltaAxis_leadingParticleFlag_processed_2025-03-10.root";
  fileName[kPythia][kNoVeto][kR0p4][kSquaredWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_energyWeightSquared_jetDeltaAxis_leadingParticleFlag_processed_2025-03-10.root";
  fileName[kPythia][kNoVeto][kR0p4][kNominalWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_nominalEnergyWeight_jetDeltaAxis_leadingParticleFlag_processed_2025-03-10.root";
  fileName[kPythia][kNoVeto][kR0p4][kSquaredWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_jetDeltaAxis_leadingParticleFlag_processed_2025-03-10.root";

  // Files for R = 0.8, no other jet veto
  fileName[kPythia][kNoVeto][kR0p8][kNominalWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_nominalEnergyWeight_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_processed_2025-03-10.root";
  fileName[kPythia][kNoVeto][kR0p8][kSquaredWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_energyWeightSquared_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_processed_2025-03-10.root";
  fileName[kPythia][kNoVeto][kR0p8][kNominalWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_nominalEnergyWeight_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_processed_2025-03-10.root";
  fileName[kPythia][kNoVeto][kR0p8][kSquaredWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_processed_2025-03-10.root";

  // Files for R = 0.4, veto if other jets within 0.8
  fileName[kPythia][kJetVeto][kR0p4][kNominalWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_nominalEnergyWeight_jetDeltaAxis_jetRadius0p4_leadingParticleFlag_veto0p8Jets_processed_2025-03-16.root";
  fileName[kPythia][kJetVeto][kR0p4][kSquaredWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_energyWeightSquared_jetDeltaAxis_jetRadius0p4_leadingParticleFlag_veto0p8Jets_processed_2025-03-16.root";
  fileName[kPythia][kJetVeto][kR0p4][kNominalWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_nominalEnergyWeight_jetDeltaAxis_jetRadius0p4_leadingParticleFlag_veto0p8Jets_processed_2025-03-16.root";
  fileName[kPythia][kJetVeto][kR0p4][kSquaredWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_jetDeltaAxis_jetRadius0p4_leadingParticleFlag_veto0p8Jets_processed_2025-03-16.root";

  // Files for R = 0.8, veto if other jets within 0.8
  fileName[kPythia][kJetVeto][kR0p8][kNominalWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_nominalEnergyWeight_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_vetoCloseJets_processed_2025-03-14.root";
  fileName[kPythia][kJetVeto][kR0p8][kSquaredWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_energyWeightSquared_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_vetoCloseJets_processed_2025-03-14.root";
  fileName[kPythia][kJetVeto][kR0p8][kNominalWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_nominalEnergyWeight_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_vetoCloseJets_processed_2025-03-14.root";
  fileName[kPythia][kJetVeto][kR0p8][kSquaredWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_vetoCloseJets_processed_2025-03-14.root";

  // ======================//
  //    Files for Herwig   //
  // ======================//

  // Files for R = 0.4, no other jet veto
  fileName[kHerwig][kNoVeto][kR0p4][kNominalWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Herwig_pfJets_eschemeAxis_nominalEnergyWeight_jetDeltaAxis_jetRadius0p4_leadingParticleFlag_processed_2025-04-30.root";
  fileName[kHerwig][kNoVeto][kR0p4][kSquaredWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Herwig_pfJets_eschemeAxis_energyWeightSquared_jetDeltaAxis_jetRadius0p4_leadingParticleFlag_processed_2025-04-30.root";
  fileName[kHerwig][kNoVeto][kR0p4][kNominalWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_nominalEnergyWeight_jetDeltaAxis_jetRadius0p4_leadingParticleFlag_processed_2025-04-30.root";
  fileName[kHerwig][kNoVeto][kR0p4][kSquaredWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_energyWeightSquared_jetDeltaAxis_jetRadius0p4_leadingParticleFlag_processed_2025-04-30.root";

  // Files for R = 0.8, no other jet veto
  fileName[kHerwig][kNoVeto][kR0p8][kNominalWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Herwig_pfJets_eschemeAxis_nominalEnergyWeight_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_processed_2025-04-30.root";
  fileName[kHerwig][kNoVeto][kR0p8][kSquaredWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Herwig_pfJets_eschemeAxis_energyWeightSquared_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_processed_2025-04-30.root";
  fileName[kHerwig][kNoVeto][kR0p8][kNominalWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_nominalEnergyWeight_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_processed_2025-04-30.root";
  fileName[kHerwig][kNoVeto][kR0p8][kSquaredWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_energyWeightSquared_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_processed_2025-04-30.root";

  // Note: There are no files for Herwig with jet veto
  
  // Open the files and check that they exist
  TFile* inputFile[knJetRadii][knEnergyWeights][knJetAxes];
  EECCard* card[knJetRadii][knEnergyWeights][knJetAxes];
  for(int iJetRadius = 0; iJetRadius < knJetRadii; iJetRadius++){
    for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
      for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
    
        inputFile[iJetRadius][iEnergyWeight][iJetAxis] = TFile::Open(fileName[iSimulation][vetoCloseJets][iJetRadius][iEnergyWeight][iJetAxis]);
    
        if(inputFile[iJetRadius][iEnergyWeight][iJetAxis] == NULL){
          cout << "Error! The file " << fileName[iSimulation][vetoCloseJets][iJetRadius][iEnergyWeight][iJetAxis].Data() << " does not exist!" << endl;
          cout << "Maybe you forgot the data/ folder path?" << endl;
          cout << "Will not execute the code" << endl;
          return;
        }

        card[iJetRadius][iEnergyWeight][iJetAxis]  = new EECCard(inputFile[iJetRadius][iEnergyWeight][iJetAxis]);
      } // Jet axis loop
    } // Energy weight loop
  } // Jet radius loop
  
  // ====================================================
  //                    Configuration
  // ====================================================
  
  // Enumeration for possible normalizations for distributions
  enum normalizationType{kNoNormalization, kIntegralToOne, kPeakToOne, knNormlizationTypes};

  // Flag for normalizing the energy-energy correlator distributions to unity
  int normalizeDistributions = kPeakToOne;

  // Find the number of bins from the reference card
  const int nJetPtBinsEEC = card[0][0][0]->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card[0][0][0]->GetNTrackPtBinsEEC();
  const int nCentralityBins = card[0][0][0]->GetNCentralityBins();
  
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

  // If we are not looking at PbPb, clear the centrality vector and add a dummy bin there
  if(!card[0][0][0]->GetDataType().Contains("PbPb")){
    comparedCentralityBin.clear();
    comparedCentralityBin.push_back(std::make_pair(-10,-10));
  }

  // For Monte Carlo, a 4% centrality shift is applied
  if(card[0][0][0]->GetDataType().Contains("MC")){
    for(auto& centralityBin : comparedCentralityBin){
      centralityBin.first += 4;
      centralityBin.second += 4;
    }
  }

  // Select which pairing type and subevent to look at

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

  bool drawLeadingParticleCheck = true;
  
  // Figure saving
  const bool saveFigures = false;  // Save figures
  TString saveComment = "_jetAxisDifference_signal";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures

  // Drawing configuration
  std::pair<double, double> leadingParticleRatioZoom = std::make_pair(0, 1);
  std::pair<double, double> eecZoom = std::make_pair(0.2, 30);
  const bool automaticZoom = true;

  // ====================================================
  //                  Histogram loading
  // ====================================================
  
  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager* histograms[knJetRadii][knEnergyWeights][knJetAxes];
  for(int iJetRadius = 0; iJetRadius < knJetRadii; iJetRadius++){
    for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
      for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){

        histograms[iJetRadius][iEnergyWeight][iJetAxis] = new EECHistogramManager(inputFile[iJetRadius][iEnergyWeight][iJetAxis], card[iJetRadius][iEnergyWeight][iJetAxis]);
          
      } // Jet axis loop for histogram loading
    } // Energy weight loop for histogram loading
  } // Jet radius loop  for histogram loading

  // Energy-energy correlator histograms
  TH1D* hEnergyEnergyCorrelator[knEnergyWeights][knJetRadii][knJetAxes][EECHistograms::kLeadingParticleTypes+1][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hEnergyEnergyCorrelatorLeadingParticleRatio[knEnergyWeights][knJetRadii][knJetAxes][EECHistograms::kLeadingParticleTypes+1][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
    for(int iLeadingParticle = 0; iLeadingParticle < EECHistograms::kLeadingParticleTypes + 1; iLeadingParticle++){
      for(int iJetRadius = 0; iJetRadius < knJetRadii; iJetRadius++){
        for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
          for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
            for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
              for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
                hEnergyEnergyCorrelator[iEnergyWeight][iJetRadius][iJetAxis][iLeadingParticle][iCentrality][iJetPt][iTrackPt] = NULL;
                hEnergyEnergyCorrelatorLeadingParticleRatio[iEnergyWeight][iJetRadius][iJetAxis][iLeadingParticle][iCentrality][iJetPt][iTrackPt] = NULL;
              } // Centrality loop
            } // Track pT loop
          } // Jet pT loop
        } // Jet axis loop
      } // Jet radius loop
    } // Correlation including leading particle
  } // Energy weight loop

  
  // Helper variables
  std::pair<double, double> drawingRange = std::make_pair(0.008, 0.75);
  double epsilon = 0.0001;
  int iTrackPt, iTrackPtReference;
  int iJetPt, iJetPtReference;
  int iCentrality, iCentralityReference;

  // Get the histograms from the histogram manager and normalize the signal histograms to one
  for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
    for(int iJetRadius = 0; iJetRadius < knJetRadii; iJetRadius++){
      for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
        for(int iLeadingParticle = 0; iLeadingParticle < EECHistograms::kLeadingParticleTypes + 1; iLeadingParticle++){
          for(auto jetPtBin : comparedJetPtBin){
            for(auto trackPtBin : comparedTrackPtBin){
              for(auto centralityBin : comparedCentralityBin){

                // Find the proper binning and express it in term of the bins in the first file
                iJetPt = card[iJetRadius][iEnergyWeight][iJetAxis]->FindBinIndexJetPtEEC(jetPtBin);
                iTrackPt = card[iJetRadius][iEnergyWeight][iJetAxis]->GetBinIndexTrackPtEEC(trackPtBin);
                iJetPtReference = card[0][0][0]->FindBinIndexJetPtEEC(jetPtBin);
                iTrackPtReference = card[0][0][0]->GetBinIndexTrackPtEEC(trackPtBin);

                if(centralityBin.first < 0){
                  iCentrality = 0;
                  iCentralityReference = 0;
                } else {
                  iCentrality = card[iJetRadius][iEnergyWeight][iJetAxis]->FindBinIndexCentrality(centralityBin);
                  iCentralityReference = card[0][0][0]->FindBinIndexCentrality(centralityBin);
                }

                // Load the selected energy-energy correlator histograms
                hEnergyEnergyCorrelator[iEnergyWeight][iJetRadius][iJetAxis][iLeadingParticle][iCentrality][iJetPt][iTrackPt] = histograms[iJetRadius][iEnergyWeight][iJetAxis]->GetHistogramEnergyEnergyCorrelatorJetDeltaAxis(EECHistogramManager::kEnergyEnergyCorrelator, EECHistogramManager::kMaxJetDeltaAxisBins, iCentrality, iJetPt, iTrackPt, iLeadingParticle, pairingType, subevent);
            
              } // Centrality loop
            } // Track pT loop
          } // Jet pT loop
        } // DeltaR between E-scheme and WTA axes loop
      } // Jet axis type loop
    } // Jet radius loop
  } // Energy weight loop

  cout << "LOADED" << endl;

  // ====================================================
  //               Normalize distributions
  // ====================================================

  int lowNormalizationBin, highNormalizationBin;
  double normalizationFactor;
  if(normalizeDistributions){
    for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
      for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
        for(int iJetRadius = 0; iJetRadius < knJetRadii; iJetRadius++){
          for(auto jetPtBin : comparedJetPtBin){
            iJetPt = card[0][0][0]->FindBinIndexJetPtEEC(jetPtBin);
            for(auto trackPtBin : comparedTrackPtBin){
              iTrackPt = card[0][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
              for(auto centralityBin : comparedCentralityBin){
                iCentrality = centralityBin.first < 0 ? 0 : card[0][0][0]->FindBinIndexCentrality(centralityBin);

                if(hEnergyEnergyCorrelator[iEnergyWeight][iJetRadius][iJetAxis][EECHistograms::kLeadingParticleTypes][iCentrality][iJetPt][iTrackPt] == NULL){
                  cout << "NOOOOO" << endl;
                  cout << "Energy weight: " << energyWeightDescription[iEnergyWeight] << endl;
                  cout << "Jet radius: " << jetRadiusDescription[iJetRadius] << endl;
                  cout << "Jet axis: " << jetAxisDescription[iJetAxis] << endl;
                  cout << "Centrality: " << iCentrality << endl;
                  cout << "Jet pT: " << iJetPt << endl;
                  cout << "Track pT: " << iTrackPt << endl;
                }

                // Find the first and last bin in the region where we normalize the distribution
                lowNormalizationBin = hEnergyEnergyCorrelator[iEnergyWeight][iJetRadius][iJetAxis][EECHistograms::kLeadingParticleTypes][iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(drawingRange.first + epsilon);
                highNormalizationBin = hEnergyEnergyCorrelator[iEnergyWeight][iJetRadius][iJetAxis][EECHistograms::kLeadingParticleTypes][iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(drawingRange.second - epsilon);

                switch(normalizeDistributions){
                  case kIntegralToOne:
                    // Normalize the distributions to one in the drawingRange
                    normalizationFactor = hEnergyEnergyCorrelator[iEnergyWeight][iJetRadius][iJetAxis][EECHistograms::kLeadingParticleTypes][iCentrality][iJetPt][iTrackPt]->Integral(lowNormalizationBin, highNormalizationBin, "width");

                    //hEnergyEnergyCorrelator[iVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->Scale(1 / hEnergyEnergyCorrelator[iVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->Integral(lowNormalizationBin, highNormalizationBin, "width"));

                    break;

                  case kPeakToOne:
                    // Normalize the distributions such that the peak height is 1
                    normalizationFactor = 0;
                    for(int iBin = lowNormalizationBin; iBin <= highNormalizationBin; iBin++){
                      if(hEnergyEnergyCorrelator[iEnergyWeight][iJetRadius][iJetAxis][EECHistograms::kLeadingParticleTypes][iCentrality][iJetPt][iTrackPt]->GetBinContent(iBin) > normalizationFactor) normalizationFactor = hEnergyEnergyCorrelator[iEnergyWeight][iJetRadius][iJetAxis][EECHistograms::kLeadingParticleTypes][iCentrality][iJetPt][iTrackPt]->GetBinContent(iBin);
                    }

                    break;

                  default:
                    // Bad input, complain to user
                    cout << "You have selected an unsupperted normalization type!" << endl;
                    cout << "Please consider your actions and try again" << endl;
                    return; 
                }

                // Once the normalization factor is determined, normalize all leading particle histograms with the same factor in order to preserve the relative magnitudes of different histograms
                for(int iLeadingParticle = 0; iLeadingParticle < EECHistograms::kLeadingParticleTypes + 1; iLeadingParticle++){
                  hEnergyEnergyCorrelator[iEnergyWeight][iJetRadius][iJetAxis][iLeadingParticle][iCentrality][iJetPt][iTrackPt]->Scale(1.0 / normalizationFactor);
                }
            
              } // Centrality loop
            } // Track pT loop
          } // Jet pT loop
        } // Jet radius loop
      } // Jet axis type loop
    } // Energy weight loop
  } // Normalizing distributions

  cout << "NORMALIZED" << endl;

  // ====================================================
  //                 Calculating ratios
  // ====================================================

  // After all the histograms have been read, calculate ratios of different leading particle pairings with respect to the whole distribution
  for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
    for(int iLeadingParticle = 0; iLeadingParticle < EECHistograms::kLeadingParticleTypes + 1; iLeadingParticle++){
      for(int iJetRadius = 0; iJetRadius < knJetRadii; iJetRadius++){
        for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
          for(auto jetPtBin : comparedJetPtBin){
            iJetPt = card[0][0][0]->FindBinIndexJetPtEEC(jetPtBin);
            for(auto trackPtBin : comparedTrackPtBin){
              iTrackPt = card[0][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
              for(auto centralityBin : comparedCentralityBin){
                iCentrality = centralityBin.first < 0 ? 0 : card[0][0][0]->FindBinIndexCentrality(centralityBin);


                hEnergyEnergyCorrelatorLeadingParticleRatio[iEnergyWeight][iJetRadius][iJetAxis][iLeadingParticle][iCentrality][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelator[iEnergyWeight][iJetRadius][iJetAxis][iLeadingParticle][iCentrality][iJetPt][iTrackPt]->Clone(Form("leadingParticleRatio%d%d%d%d%d%d%d", iEnergyWeight, iLeadingParticle, iJetRadius, iJetAxis, iCentrality, iJetPt, iTrackPt));
                
                hEnergyEnergyCorrelatorLeadingParticleRatio[iEnergyWeight][iJetRadius][iJetAxis][iLeadingParticle][iCentrality][iJetPt][iTrackPt]->Divide(hEnergyEnergyCorrelator[iEnergyWeight][iJetRadius][iJetAxis][EECHistograms::kLeadingParticleTypes][iCentrality][iJetPt][iTrackPt]);

              } // Centrality loop
            } // Track pT loop
          } // Jet pT loop
        } // Jet axis type loop
      } // Jet radius loop
    } // Leading particle type loop
  } // Energy weight loop

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
  int markerStyle[] = {kFullCircle, kFullDiamond, kFullSquare, kOpenCircle, kFullCross};
  int color[] = {kBlue, kRed, kBlack, kGreen+3, kMagenta,kCyan,kOrange,kViolet+3,kPink-7,kSpring+3,kAzure-7};

  // Binning vectors
  std::vector<int> currentJetPtIndices;
  std::vector<int> currentTrackPtIndices;
  double minimumCandidate, maximumCandidate;

  // Legend and line
  TLegend* legend;
  TLine* oneLine = new TLine(drawingRange.first, 1, drawingRange.second, 1);
  oneLine->SetLineStyle(2);
  
  // Draw the single ratio histograms
  if(drawLeadingParticleCheck){

    for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
      for(int iJetRadius = 0; iJetRadius < knJetRadii; iJetRadius++){
        for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
          for(auto jetPtBin : comparedJetPtBin){
            iJetPt = card[0][0][0]->FindBinIndexJetPtEEC(jetPtBin);
            for(auto trackPtBin : comparedTrackPtBin){
              iTrackPt = card[0][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
              for(auto centralityBin : comparedCentralityBin){
                iCentrality = centralityBin.first < 0 ? 0 : card[0][0][0]->FindBinIndexCentrality(centralityBin);
          
                // Create a new canvas for the plot
                drawer->CreateSplitCanvas();

                // Logarithmic EEC axis
                drawer->SetLogY(true);

                // Create the legend and add binning information to it
                legend = new TLegend(0.18,0.04,0.45,0.58);
                legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);

                legend->AddEntry((TObject*) 0, simulationName[iSimulation].Data(), "");
                legend->AddEntry((TObject*) 0, Form("Jet axis: %s", jetAxisDescription[iJetAxis].Data()), "");
                legend->AddEntry((TObject*) 0, Form("%s, %s", jetRadiusDescription[iJetRadius].Data(), energyWeightDescription[iEnergyWeight].Data()), "");
                if(centralityBin.first >= 0) legend->AddEntry((TObject*) 0, Form("Centrality: %d-%d%%", centralityBin.first, centralityBin.second), "");
                legend->AddEntry((TObject*) 0, Form("%.0f < jet p_{T} < %.0f GeV", jetPtBin.first, jetPtBin.second), "");
                legend->AddEntry((TObject*) 0, Form("p_{T}^{ch} > %.1f GeV", trackPtBin), "");

                // Setup centrality, jet pT and track pT strings
                compactJetPtString = Form("_J=%.0f-%.0f", jetPtBin.first, jetPtBin.second);
                compactTrackPtString = Form("_T>%.1f",trackPtBin);
                compactTrackPtString.ReplaceAll(".","v");
                compactCentralityString = centralityBin.first == -1 ? "_pp" : Form("_C=%d-%d", centralityBin.first, centralityBin.second);


                // Set drawing style for all histograms
                for(int iLeadingParticle = 0; iLeadingParticle < EECHistograms::kLeadingParticleTypes + 1; iLeadingParticle++){
                  hEnergyEnergyCorrelator[iEnergyWeight][iJetRadius][iJetAxis][iLeadingParticle][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iLeadingParticle]);
                  hEnergyEnergyCorrelator[iEnergyWeight][iJetRadius][iJetAxis][iLeadingParticle][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iLeadingParticle]);
                  hEnergyEnergyCorrelator[iEnergyWeight][iJetRadius][iJetAxis][iLeadingParticle][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iLeadingParticle]); 
                  hEnergyEnergyCorrelatorLeadingParticleRatio[iEnergyWeight][iJetRadius][iJetAxis][iLeadingParticle][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iLeadingParticle]);
                  hEnergyEnergyCorrelatorLeadingParticleRatio[iEnergyWeight][iJetRadius][iJetAxis][iLeadingParticle][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iLeadingParticle]);
                  hEnergyEnergyCorrelatorLeadingParticleRatio[iEnergyWeight][iJetRadius][iJetAxis][iLeadingParticle][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iLeadingParticle]);
                } // Leading particle type loop

                // Automatic zooming for the drawn histograms
                if(automaticZoom){
                  eecZoom.first = 10000;
                  eecZoom.second = 0;
                  for(int iLeadingParticle = 0; iLeadingParticle < EECHistograms::kLeadingParticleTypes + 1; iLeadingParticle++){
                    hEnergyEnergyCorrelator[iEnergyWeight][iJetRadius][iJetAxis][iLeadingParticle][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
                    minimumCandidate = hEnergyEnergyCorrelator[iEnergyWeight][iJetRadius][iJetAxis][iLeadingParticle][iCentrality][iJetPt][iTrackPt]->GetMinimum();
                    maximumCandidate = hEnergyEnergyCorrelator[iEnergyWeight][iJetRadius][iJetAxis][iLeadingParticle][iCentrality][iJetPt][iTrackPt]->GetMaximum();
                    if(minimumCandidate < eecZoom.first) eecZoom.first = minimumCandidate;
                    if(maximumCandidate > eecZoom.second) eecZoom.second = maximumCandidate;
                  } // Leading particle type loop
                  eecZoom.first = eecZoom.first / 2.0;
                  eecZoom.second = eecZoom.second * 2.0;
                }

                // Set the x- and y-axis drawing ranges
                hEnergyEnergyCorrelator[iEnergyWeight][iJetRadius][iJetAxis][EECHistograms::kLeadingParticleTypes][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
                hEnergyEnergyCorrelator[iEnergyWeight][iJetRadius][iJetAxis][EECHistograms::kLeadingParticleTypes][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(eecZoom.first, eecZoom.second);

                // Draw the histograms to the upper canvas
                drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelator[iEnergyWeight][iJetRadius][iJetAxis][EECHistograms::kLeadingParticleTypes][iCentrality][iJetPt][iTrackPt], "#Deltar", "EEC", " ", "p");
                legend->AddEntry(hEnergyEnergyCorrelator[iEnergyWeight][iJetRadius][iJetAxis][EECHistograms::kLeadingParticleTypes][iCentrality][iJetPt][iTrackPt], leadingParticleName[EECHistograms::kLeadingParticleTypes].Data(), "p");

                for(int iLeadingParticle = 0; iLeadingParticle < EECHistograms::kLeadingParticleTypes; iLeadingParticle++){
                  hEnergyEnergyCorrelator[iEnergyWeight][iJetRadius][iJetAxis][iLeadingParticle][iCentrality][iJetPt][iTrackPt]->Draw("same,p");
                  legend->AddEntry(hEnergyEnergyCorrelator[iEnergyWeight][iJetRadius][iJetAxis][iLeadingParticle][iCentrality][iJetPt][iTrackPt], leadingParticleName[iLeadingParticle].Data(), "p");
                } // Leading particle type loop
  
                // Draw the legends to the upper pad
                legend->Draw();
          
                // Linear scale for the ratio
                drawer->SetLogY(false);
          
                // Set the drawing ranges
                hEnergyEnergyCorrelatorLeadingParticleRatio[iEnergyWeight][iJetRadius][iJetAxis][EECHistograms::kNotLeadingParticle][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
                hEnergyEnergyCorrelatorLeadingParticleRatio[iEnergyWeight][iJetRadius][iJetAxis][EECHistograms::kNotLeadingParticle][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(leadingParticleRatioZoom.first, leadingParticleRatioZoom.second);

                // Draw the histograms
                drawer->SetGridY(true);

                drawer->DrawHistogramToLowerPad(hEnergyEnergyCorrelatorLeadingParticleRatio[iEnergyWeight][iJetRadius][iJetAxis][EECHistograms::kNotLeadingParticle][iCentrality][iJetPt][iTrackPt], "#Deltar", "#frac{Selected pairs}{All pairs}", " ", "p");
                hEnergyEnergyCorrelatorLeadingParticleRatio[iEnergyWeight][iJetRadius][iJetAxis][EECHistograms::kLeadingParticle][iCentrality][iJetPt][iTrackPt]->Draw("p,same");

                drawer->SetGridY(false);
          
                // Save the figures to a file
                if(saveFigures){
                  gPad->GetCanvas()->SaveAs(Form("figures/eecLeadingParticleComponents%s%s%s%s%s%s%s%s.%s", saveComment.Data(), simulationSaveName[iSimulation].Data(), compactEnergyWeightString[iEnergyWeight].Data(), jetRadiusSaveName[iJetRadius].Data(), jetAxisSaveName[iJetAxis].Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
                } // Saving figures
              } // Centrality loop
            } // Track pT loop
          } // Jet pT loop
        } // Jet axis loop
      } // Jet radius loop
    } // Energy weight loop
  } // Drawing leading particle check

}
