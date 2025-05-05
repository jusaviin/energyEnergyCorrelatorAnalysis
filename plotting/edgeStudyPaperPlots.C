#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"
#include "../src/EECHistograms.h"

void drawVector(JDrawer* drawer, std::vector<std::pair<TH1D*,TString>> histogramVector, TString histogramLegendHeader, std::vector<TString> systemLegendVector, TString axisName, std::pair<double, double> xRange, std::pair<double, double> yRange, TString saveName);

/*
 * Macro producing the planned paper plots for the EEC edge study. These are
 *
 *  1) Average DeltaR between E-scheme and WTA axes as a function of jet pT
 *  2) Plot showing edge-loss from the R=0.4/R=0.8 ratio. The loss is difference of this integral from one in region ~0.2-0.4
 *  3) Plot of the edge-loss as a function of jet pT
 *  4) Plot of the edge loss as a function of average DeltaR between E-scheme and WTA axes corrending to different jet pT bins
 *
 * At least in the beginning, the study is only done in Pythia8. More systems might be added later.
 */
void edgeStudyPaperPlots(){

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
  const bool vetoCloseJets = false;

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
  TFile* inputFile[knSimulations][knVetoTypes][knJetRadii][knEnergyWeights][knJetAxes];
  EECCard* card[knSimulations][knVetoTypes][knJetRadii][knEnergyWeights][knJetAxes];
  for(int iSimulation = 0; iSimulation < knSimulations; iSimulation++){
    for(int iVeto = 0; iVeto < knVetoTypes; iVeto++){
      if(iSimulation == kHerwig && iVeto == kJetVeto) continue; // Note: We have no files with jet veto for Herwig
      for(int iJetRadius = 0; iJetRadius < knJetRadii; iJetRadius++){
        for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
          for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
    
            inputFile[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis] = TFile::Open(fileName[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis]);
    
            if(inputFile[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis] == NULL){
              cout << "Error! The file " << fileName[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis].Data() << " does not exist!" << endl;
              cout << "Maybe you forgot the data/ folder path?" << endl;
              cout << "Will not execute the code" << endl;
              return;
            }

            card[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis]  = new EECCard(inputFile[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis]);
          } // Jet axis loop
        } // Energy weight loop
      } // Jet radius loop
    } // Veto loop
  } // Simulation loop

  // ====================================================
  //         Configuration for derived histograms
  // ====================================================

  // Enumeration for possible normalizations for distributions
  enum normalizationType{kNoNormalization, kIntegralToOne, kPeakToOne, knNormlizationTypes};

  // Flag for normalizing the energy-energy correlator distributions to unity
  int normalizeDistributions = kPeakToOne;

  TString normalizationDescription[knNormlizationTypes];
  normalizationDescription[kNoNormalization] = "not normalized";
  normalizationDescription[kIntegralToOne] = "integral to one";
  normalizationDescription[kPeakToOne] = "peak height to one";

  // Range from which the edge loss is determined
  std::pair<double,double> edgeLossRegion = std::make_pair(0.05, 0.39);

  // Enumeration for different ways to quantify edge loss
  enum edgeLossTypeType{kRelativeEdgeLoss, kAbsoluteEdgeLoss, knEdgeLossMethods};

  // Select which edge loss method is used
  int edgeLossMethod = kRelativeEdgeLoss;

  TString edgeLossString[knEdgeLossMethods];
  edgeLossString[kRelativeEdgeLoss] = "Relative";
  edgeLossString[kAbsoluteEdgeLoss] = "Absolute";
  
  // ====================================================
  //               Binning configuration
  // ====================================================
  
  // Find the number of bins from the reference card
  const int nJetPtBinsEEC = card[0][0][0][0][0]->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card[0][0][0][0][0]->GetNTrackPtBinsEEC();
  const int nCentralityBins = card[0][0][0][0][0]->GetNCentralityBins();
  
  std::vector<std::pair<int,int>> comparedCentralityBin;
  comparedCentralityBin.push_back(std::make_pair(0,10));

  std::vector<std::pair<int,int>> comparedJetPtBin;
  comparedJetPtBin.push_back(std::make_pair(60,80));
  comparedJetPtBin.push_back(std::make_pair(80,100));
  comparedJetPtBin.push_back(std::make_pair(100,120));
  comparedJetPtBin.push_back(std::make_pair(120,140));
  comparedJetPtBin.push_back(std::make_pair(140,160));
  comparedJetPtBin.push_back(std::make_pair(160,180));
  comparedJetPtBin.push_back(std::make_pair(180,200));
  comparedJetPtBin.push_back(std::make_pair(200,220));
  comparedJetPtBin.push_back(std::make_pair(220,240));
  comparedJetPtBin.push_back(std::make_pair(240,260));
  comparedJetPtBin.push_back(std::make_pair(260,280));
  comparedJetPtBin.push_back(std::make_pair(280,300));

  const int nJetPtBins = comparedJetPtBin.size();

  std::vector<double> comparedTrackPtBin;
  comparedTrackPtBin.push_back(1.0);
  //comparedTrackPtBin.push_back(2.0);
  //comparedTrackPtBin.push_back(3.0);

  // We will use the comparedJetPtBin vector to construct a histogram axis, so the bins must be continuous. Ensure that this is the case
  for(int iJetPt = 0; iJetPt < nJetPtBins-1; iJetPt++){
    if(comparedJetPtBin.at(iJetPt).second != comparedJetPtBin.at(iJetPt+1).first){
      cout << "ERROR! The bins in comparedJetPtBin must be continuous!" << endl;
      cout << "However, you have configured bin " << iJetPt << " to be " << comparedJetPtBin.at(iJetPt).first << "-" << comparedJetPtBin.at(iJetPt).second << " and bin " << iJetPt+1 << " to be " << comparedJetPtBin.at(iJetPt+1).first << "-" << comparedJetPtBin.at(iJetPt+1).second << endl;
      cout << "This will cause a problem since these bins are used to construct an axis for a histogram!" << endl;
      cout << "Please fix this in the configuration" << endl;
      return;
    }
  }

  // If we are not looking at PbPb, clear the centrality vector and add a dummy bin there
  if(!card[0][0][0][0][0]->GetDataType().Contains("PbPb")){
    comparedCentralityBin.clear();
    comparedCentralityBin.push_back(std::make_pair(-10,-10));
  }

  // For Monte Carlo, a 4% centrality shift is applied
  if(card[0][0][0][0][0]->GetDataType().Contains("MC")){
    for(auto& centralityBin : comparedCentralityBin){
      centralityBin.first += 4;
      centralityBin.second += 4;
    }
  }

  // Choose which pairs are used for the energy-energy correlators

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

  // Paper plot selection
  bool drawDeltaJetAxisWithPt = false;        // Draw DeltaR between E-scheme and WTA jet axes as a function of pT
  bool fitOneOverPt = false;                  // Fit 1/pT function to average DeltaR between E-scheme and WTA jet axes as a function of pT histograms
  bool drawEdgeLossIllustration = false;      // Draw an illustration about edge loss
  bool drawEdgeLossWithPt = false;            // Draw the amount of edge-loss as a function of jet pT
  bool drawEdgeLossWithDeltaJetAxis = true;  // Draw the amount of edge-loss as a function of DeltaR between E-scheme and WTA jet axes 
  bool fitLinearToEdgeLoss = true;

  // QA plots
  bool drawJetVetoBiasIllustration = false;    // Draw plots illustrating how much vetoing close jets biases the distributions
  bool drawEdgeLossDefinition = false;         // Illustrate the definition of edge loss
  bool drawJetDeltaAxis = false;                // Draw the delta R distributions between jet axes
  
  // Figure saving
  const bool saveFigures = false;  // Save figures
  TString saveComment = "_jetAxisDifference_signal";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures

  // Drawing configuration
  std::pair<double, double> deltaJetAxisZoom = std::make_pair(0.02, 0.05);
  std::pair<double, double> edgeLossIllustrationZoom = std::make_pair(0.75, 1.25);
  std::pair<double, double> edgeLossZoom[knEdgeLossMethods];
  edgeLossZoom[kRelativeEdgeLoss] = std::make_pair(0.0, 0.03);
  edgeLossZoom[kAbsoluteEdgeLoss] = std::make_pair(0.0, 0.005);
  const bool automaticZoom = true;

  // ====================================================
  //                  Histogram loading
  // ====================================================
  
  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager* histograms[knSimulations][knVetoTypes][knJetRadii][knEnergyWeights][knJetAxes];
  for(int iSimulation = 0; iSimulation < knSimulations; iSimulation++){
    for(int iVeto = 0; iVeto < knVetoTypes; iVeto++){
      if(iSimulation == kHerwig && iVeto == kJetVeto) continue; // Note: We have no files with jet veto for Herwig
      for(int iJetRadius = 0; iJetRadius < knJetRadii; iJetRadius++){
        for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
          for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){

            histograms[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis] = new EECHistogramManager(inputFile[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis], card[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis]);
          
          } // Jet axis loop for histogram loading
        } // Energy weight loop for histogram loading
      } // Jet radius loop  for histogram loading
    } // Jet veto loop for histogram loading
  } // Simulation type loop for histogram loading

  // Energy-energy correlator histograms
  TH1D* hEnergyEnergyCorrelator[knSimulations][knVetoTypes][knJetRadii][knEnergyWeights][knJetAxes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hEnergyEnergyCorrelatorAxisRatio[knSimulations][knVetoTypes][knJetRadii][knEnergyWeights][knJetAxes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hEnergyEnergyCorrelatorRadiusRatio[knSimulations][knVetoTypes][knJetRadii][knEnergyWeights][knJetAxes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];

  // Histograms for DeltaR between E-scheme and WTA axes
  TH1D* hDeltaJetAxis[knSimulations][knJetAxes][nCentralityBins][nJetPtBinsEEC];
  TH1D* hAverageDeltaJetAxisVsPt[knSimulations][knJetAxes][nCentralityBins];

  // Histograms for edge loss amount as a function of jet pT and average DeltaR between E-scheme and WTA axes
  TH1D* hEdgeLossVsPt[knSimulations][knVetoTypes][knEnergyWeights][knJetAxes][nCentralityBins][nTrackPtBinsEEC];
  TGraphErrors* hEdgeLossVsDeltaJetAxis[knSimulations][knVetoTypes][knEnergyWeights][knJetAxes][nCentralityBins][nTrackPtBinsEEC];

  // Quantified edge loss for different kinematic selections
  double edgeLossAmount[knSimulations][knEdgeLossMethods][knVetoTypes][knEnergyWeights][knJetAxes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  double edgeLossError[knSimulations][knEdgeLossMethods][knVetoTypes][knEnergyWeights][knJetAxes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];

  // Fit functions
  TF1* fOneOverPtFit[knSimulations][knJetAxes][nCentralityBins];
  TF1* fLinearFit[knSimulations][knVetoTypes][knEnergyWeights][knJetAxes][nCentralityBins][nTrackPtBinsEEC];

  // QA histograms to study effect of certain selections
  TH1D* hEnergyEnergyCorrelatorVetoBiasCheck[knSimulations][knVetoTypes][knJetRadii][knEnergyWeights][knJetAxes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iSimulation = 0; iSimulation < knSimulations; iSimulation++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
        hAverageDeltaJetAxisVsPt[iSimulation][iJetAxis][iCentrality] = NULL;
        if(fitOneOverPt){
          fOneOverPtFit[iSimulation][iJetAxis][iCentrality] = new TF1(Form("oneOverPtFit%d%d", iCentrality, iJetAxis), "[0] + [1]/x", comparedJetPtBin.at(0).first, comparedJetPtBin.at(comparedJetPtBin.size()-1).second);
          fOneOverPtFit[iSimulation][iJetAxis][iCentrality]->SetParameters(1,1);
        } else {
          fOneOverPtFit[iSimulation][iJetAxis][iCentrality] = NULL;
        }
        for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
          hDeltaJetAxis[iSimulation][iJetAxis][iCentrality][iJetPt] = NULL;
          for(int iVeto = 0; iVeto < knVetoTypes; iVeto++){
            for(int iJetRadius = 0; iJetRadius < knJetRadii; iJetRadius++){
              for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
                for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){ 
                  hEnergyEnergyCorrelator[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt] = NULL;
                  hEnergyEnergyCorrelatorAxisRatio[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt] = NULL;
                  hEnergyEnergyCorrelatorRadiusRatio[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt] = NULL;
                  hEnergyEnergyCorrelatorVetoBiasCheck[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt] = NULL;
                } // Track pT loop
              } // Energy weight loop
            } // Jet radius loop
          } // Jet veto loop
        } // Jet pT loop
        for(int iVeto = 0; iVeto < knVetoTypes; iVeto++){
          for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
            for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
              hEdgeLossVsPt[iSimulation][iVeto][iEnergyWeight][iJetAxis][iCentrality][iTrackPt] = NULL;
              hEdgeLossVsDeltaJetAxis[iSimulation][iVeto][iEnergyWeight][iJetAxis][iCentrality][iTrackPt] = NULL;
              if(fitLinearToEdgeLoss){
                fLinearFit[iSimulation][iVeto][iEnergyWeight][iJetAxis][iCentrality][iTrackPt] = new TF1(Form("linearFit%d%d%d%d%d%d", iSimulation, iVeto, iEnergyWeight, iJetAxis, iCentrality, iTrackPt), "pol1", 0.025, 0.055);
                fLinearFit[iSimulation][iVeto][iEnergyWeight][iJetAxis][iCentrality][iTrackPt]->SetParameters(0,1);
              } else {
                fLinearFit[iSimulation][iVeto][iEnergyWeight][iJetAxis][iCentrality][iTrackPt] = NULL;
              }
            } // Track pT loop
          } // Energy weight loop
        } // Jet veto loop 
      } // Jet axis loop
    } // Centrality loop 

    // Initialize the edge loss amounts to 0
    for(int iEdgeLoss = 0; iEdgeLoss < knEdgeLossMethods; iEdgeLoss++){
      for(int iVeto = 0; iVeto < knVetoTypes; iVeto++){
        for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
          for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
            for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
              for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
                for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){ 
                  edgeLossAmount[iSimulation][iEdgeLoss][iVeto][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt] = 0;
                  edgeLossError[iSimulation][iEdgeLoss][iVeto][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt] = 0;
                } // Track pT loop
              } // Jet pT loop
            } // Centrality loop
          } // Jet axis type loop
        } // Energy weight loop
      } // Jet veto loop
    } // Used edge loss method loop
  } // Simulation loop
  
  // Helper variables
  std::pair<double, double> drawingRange = std::make_pair(0.008, 0.39);
  double epsilon = 0.0001;
  int iTrackPt, iTrackPtReference;
  int iJetPt, iJetPtReference;
  int iCentrality, iCentralityReference;

  // Load the energy-energy correlator histograms
  for(int iSimulation = 0; iSimulation < knSimulations; iSimulation++){
    for(int iVeto = 0; iVeto < knVetoTypes; iVeto++){
      if(iSimulation == kHerwig && iVeto == kJetVeto) continue; // Note: We have no files with jet veto for Herwig
      for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
        for(int iJetRadius = 0; iJetRadius < knJetRadii; iJetRadius++){
          for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
            for(auto centralityBin : comparedCentralityBin){
              for(auto jetPtBin : comparedJetPtBin){
                for(auto trackPtBin : comparedTrackPtBin){

                  // Find the proper binning and express it in term of the bins in the first file
                  iJetPt = card[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis]->FindBinIndexJetPtEEC(jetPtBin);
                  iTrackPt = card[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis]->GetBinIndexTrackPtEEC(trackPtBin);
                  iJetPtReference = card[0][0][0][0][0]->FindBinIndexJetPtEEC(jetPtBin);
                  iTrackPtReference = card[0][0][0][0][0]->GetBinIndexTrackPtEEC(trackPtBin);

                  if(centralityBin.first < 0){
                    iCentrality = 0;
                    iCentralityReference = 0;
                  } else {
                    iCentrality = card[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis]->FindBinIndexCentrality(centralityBin);
                    iCentralityReference = card[0][0][0][0][0]->FindBinIndexCentrality(centralityBin);
                  }

                  hEnergyEnergyCorrelator[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentralityReference][iJetPtReference][iTrackPtReference] = histograms[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, pairingType, subevent);
            
                } // Track pT loop
              } // Jet pT loop
            } // Centrality loop
          } // Jet axis type loop
        } // Jet radius loop
      } // Energy weight loop
    } // Jet veto loop 
  } // Simulation loop

  // Load the DeltaR between E-scheme and WTA axes histograms
  for(int iSimulation = 0; iSimulation < knSimulations; iSimulation++){
    for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
      for(auto jetPtBin : comparedJetPtBin){ 
        iJetPt = card[iSimulation][vetoCloseJets][kR0p4][kNominalWeight][iJetAxis]->FindBinIndexJetPtEEC(jetPtBin);
        iJetPtReference = card[0][0][0][0][0]->FindBinIndexJetPtEEC(jetPtBin);
        for(auto centralityBin : comparedCentralityBin){

          if(centralityBin.first < 0){
            iCentrality = 0;
            iCentralityReference = 0;
          } else {
            iCentrality = card[iSimulation][vetoCloseJets][kR0p4][kNominalWeight][iJetAxis]->FindBinIndexCentrality(centralityBin);
            iCentralityReference = card[0][0][0][0][0]->FindBinIndexCentrality(centralityBin);
          } 

          hDeltaJetAxis[iSimulation][iJetAxis][iCentralityReference][iJetPtReference] = histograms[iSimulation][vetoCloseJets][kR0p4][kNominalWeight][iJetAxis]->GetHistogramJetDeltaAxis(iCentrality, iJetPt);
        } // Centrality loop
      } // Jet pT loop
    } // Jet axis loop
  } // Simulation loop

  // ====================================================
  //               Normalize distributions
  // ====================================================

  int lowNormalizationBin, highNormalizationBin;
  double maxBinContent;
  if(normalizeDistributions){
    for(int iSimulation = 0; iSimulation < knSimulations; iSimulation++){
      for(int iVeto = 0; iVeto < knVetoTypes; iVeto++){
        if(iSimulation == kHerwig && iVeto == kJetVeto) continue; // Note: We have no files with jet veto for Herwig
        for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
          for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
            for(int iJetRadius = 0; iJetRadius < knJetRadii; iJetRadius++){
              for(auto jetPtBin : comparedJetPtBin){
                iJetPt = card[0][0][0][0][0]->FindBinIndexJetPtEEC(jetPtBin);
                for(auto trackPtBin : comparedTrackPtBin){
                  iTrackPt = card[0][0][0][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
                  // Do the same for all the selected centrality bins for PbPb
                  for(auto centralityBin : comparedCentralityBin){
                    iCentrality = centralityBin.first < 0 ? 0 : card[0][0][0][0][0]->FindBinIndexCentrality(centralityBin);

                    // Find the first and last bin in the region where we normalize the distribution
                    lowNormalizationBin = hEnergyEnergyCorrelator[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(drawingRange.first + epsilon);
                    highNormalizationBin = hEnergyEnergyCorrelator[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(drawingRange.second - epsilon);

                    switch(normalizeDistributions){
                    case kIntegralToOne:
                      // Normalize the distributions to one in the drawingRange
                      hEnergyEnergyCorrelator[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->Scale(1 / hEnergyEnergyCorrelator[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->Integral(lowNormalizationBin, highNormalizationBin, "width"));

                      break;

                    case kPeakToOne:
                      // Normalize the distributions such that the peak height is 1
                      maxBinContent = 0;
                      for(int iBin = lowNormalizationBin; iBin <= highNormalizationBin; iBin++){
                        if(hEnergyEnergyCorrelator[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->GetBinContent(iBin) > maxBinContent) maxBinContent = hEnergyEnergyCorrelator[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->GetBinContent(iBin);
                      }
                      hEnergyEnergyCorrelator[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->Scale(1 / maxBinContent);

                      break;

                    default:
                      // Bad input, complain to user
                      cout << "You have selected an unsupported normalization type!" << endl;
                      cout << "Please consider your actions and try again" << endl;
                      return; 
                    }
            
                  } // Centrality loop
                } // Track pT loop
              } // Jet pT loop
            } // Jet radius loop
          } // Jet axis type loop
        } // Energy weight loop
      } // Jet veto loop
    } // Simulation loop
  } // Normalizing distributions

  // ====================================================
  //                 Calculating ratios
  // ====================================================

  // Ratio calculation for WTA/E-scheme
  for(int iSimulation = 0; iSimulation < knSimulations; iSimulation++){
    for(int iVeto = 0; iVeto < knVetoTypes; iVeto++){
      if(iSimulation == kHerwig && iVeto == kJetVeto) continue; // Note: We have no files with jet veto for Herwig
      for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
        for(int iJetRadius = 0; iJetRadius < knJetRadii; iJetRadius++){
          for(auto jetPtBin : comparedJetPtBin){
            iJetPt = card[0][0][0][0][0]->FindBinIndexJetPtEEC(jetPtBin);
            for(auto trackPtBin : comparedTrackPtBin){
              iTrackPt = card[0][0][0][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
              for(auto centralityBin: comparedCentralityBin){
                iCentrality = centralityBin.first < 0 ? 0 : card[0][0][0][0][0]->FindBinIndexCentrality(centralityBin);
                for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){

                  hEnergyEnergyCorrelatorAxisRatio[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelator[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->Clone(Form("eecAxisRatio%d%d%d%d%d%d%d%d", iSimulation, iVeto, iJetRadius, iEnergyWeight, iJetAxis, iCentrality, iJetPt, iTrackPt));
                  hEnergyEnergyCorrelatorAxisRatio[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->Divide(hEnergyEnergyCorrelator[iSimulation][iVeto][iJetRadius][iEnergyWeight][kEscheme][iCentrality][iJetPt][iTrackPt]);
                } // Jet axis loop
              } // Centrality loop
            } // Track pT loop
          } // Jet pT loop
        } // Jet radius loop
      } // Energy weight loop
    } // Jet veto loop
  } // Simulation loop

  // Ratio calculation for R=0.4/R=0.8
  for(int iSimulation = 0; iSimulation < knSimulations; iSimulation++){
    for(int iVeto = 0; iVeto < knVetoTypes; iVeto++){
      if(iSimulation == kHerwig && iVeto == kJetVeto) continue; // Note: We have no files with jet veto for Herwig
      for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
        for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
          for(auto jetPtBin : comparedJetPtBin){
            iJetPt = card[0][0][0][0][0]->FindBinIndexJetPtEEC(jetPtBin);
            for(auto trackPtBin : comparedTrackPtBin){
              iTrackPt = card[0][0][0][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
              for(auto centralityBin: comparedCentralityBin){
                iCentrality = centralityBin.first < 0 ? 0 : card[0][0][0][0][0]->FindBinIndexCentrality(centralityBin);
                for(int iJetRadius = 0; iJetRadius < knJetRadii; iJetRadius++){

                  hEnergyEnergyCorrelatorRadiusRatio[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelator[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->Clone(Form("eecRadiusRatio%d%d%d%d%d%d%d%d", iSimulation, iVeto, iJetRadius, iEnergyWeight, iJetAxis, iCentrality, iJetPt, iTrackPt));
                  hEnergyEnergyCorrelatorRadiusRatio[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->Divide(hEnergyEnergyCorrelator[iSimulation][iVeto][kR0p8][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]);
                } // Jet axis loop
              } // Centrality loop
            } // Track pT loop
          } // Jet pT loop
        } // Jet radius loop
      } // Energy weight loop
    } // Jet veto loop
  } // Simulation loop

  // Ratio calculation for Jet veto/No veto
  for(int iSimulation = 0; iSimulation < knSimulations; iSimulation++){
    if(iSimulation == kHerwig) continue; // Note: We have no files with jet veto for Herwig
    for(int iVeto = 0; iVeto < knVetoTypes; iVeto++){
      for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
        for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
          for(auto jetPtBin : comparedJetPtBin){
            iJetPt = card[0][0][0][0][0]->FindBinIndexJetPtEEC(jetPtBin);
            for(auto trackPtBin : comparedTrackPtBin){
              iTrackPt = card[0][0][0][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
              for(auto centralityBin: comparedCentralityBin){
                iCentrality = centralityBin.first < 0 ? 0 : card[0][0][0][0][0]->FindBinIndexCentrality(centralityBin);
                for(int iJetRadius = 0; iJetRadius < knJetRadii; iJetRadius++){

                  hEnergyEnergyCorrelatorVetoBiasCheck[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelator[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->Clone(Form("eecVetoBiasCheck%d%d%d%d%d%d%d%d", iSimulation, iVeto, iJetRadius, iEnergyWeight, iJetAxis, iCentrality, iJetPt, iTrackPt));
                  hEnergyEnergyCorrelatorVetoBiasCheck[iSimulation][iVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->Divide(hEnergyEnergyCorrelator[iSimulation][kNoVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]);
                } // Jet axis loop
              } // Centrality loop
            } // Track pT loop
          } // Jet pT loop
        } // Jet radius loop
      } // Energy weight loop
    } // Jet veto loop
  } // Simulation loop

  // ====================================================
  //              Creating derived histograms
  // ====================================================

  // Create histograms with average DeltaR between E-scheme and WTA axes as a function of jet pT
  double jetPtBinBorderArray[nJetPtBins+1];
  double meanDeltaJetAxis, meanDeltaJetAxisError;
  int iBin = 1;

  // Collect the jet pT borders into an array that is used to define the x-axis for the average DeltaR between E-scheme and WTA axes as a function of jet pT histogram
  iJetPt = 0;
  for(auto jetPtBin : comparedJetPtBin){
    jetPtBinBorderArray[iJetPt++] = jetPtBin.first;
  }
  jetPtBinBorderArray[iJetPt] = comparedJetPtBin.at(iJetPt-1).second;

  // Create the histograms
  for(int iSimulation = 0; iSimulation < knSimulations; iSimulation++){
    for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
      for(auto centralityBin : comparedCentralityBin){
        iCentrality = centralityBin.first < 0 ? 0 : card[0][0][0][0][0]->FindBinIndexCentrality(centralityBin);
      
        // Make a new histogram with the x-axis defined by the determined jet pT bins
        hAverageDeltaJetAxisVsPt[iSimulation][iJetAxis][iCentrality] = new TH1D(Form("averageDeltaJetAxis%d%d%d", iSimulation, iJetAxis, iCentrality), Form("averageDeltaJetAxis%d%d%d", iSimulation, iJetAxis, iCentrality), nJetPtBins, jetPtBinBorderArray);
        hAverageDeltaJetAxisVsPt[iSimulation][iJetAxis][iCentrality]->Sumw2();

        // Once the histogram is created, fill each bin with the mean value from the DeltaJetAxis histograms
        iBin = 1;
        for(auto jetPtBin : comparedJetPtBin){
          iJetPt = card[0][0][0][0][0]->FindBinIndexJetPtEEC(jetPtBin);

          meanDeltaJetAxis = hDeltaJetAxis[iSimulation][iJetAxis][iCentrality][iJetPt]->GetMean();
          meanDeltaJetAxisError = hDeltaJetAxis[iSimulation][iJetAxis][iCentrality][iJetPt]->GetMeanError();

          hAverageDeltaJetAxisVsPt[iSimulation][iJetAxis][iCentrality]->SetBinContent(iBin, meanDeltaJetAxis);
          hAverageDeltaJetAxisVsPt[iSimulation][iJetAxis][iCentrality]->SetBinError(iBin++, meanDeltaJetAxisError);
        }

        if(fitOneOverPt){
          hAverageDeltaJetAxisVsPt[iSimulation][iJetAxis][iCentrality]->Fit(fOneOverPtFit[iSimulation][iJetAxis][iCentrality], "N");
        }
      } // Centrality loop
    } // Jet axis loop
  } // Simulation loop

  // Define the amount of edge loss from the integral of the ratio of different radii
  // Ratio calculation for R=0.4/R=0.8
  int firstEdgeLossBin, lastEdgeLossBin;
  double edgeLossLowerBoundary, edgeLossUpperBoundary;
  double totalArea, totalAreaError;
  double edgeLossIntegral, edgeLossIntegralError;
  for(int iSimulation = 0; iSimulation < knSimulations; iSimulation++){
    for(int iVeto = 0; iVeto < knVetoTypes; iVeto++){
      if(iSimulation == kHerwig && iVeto == kJetVeto) continue; // Note: We have no files with jet veto for Herwig
      for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
        for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
          for(auto jetPtBin : comparedJetPtBin){
            iJetPt = card[0][0][0][0][0]->FindBinIndexJetPtEEC(jetPtBin);
            for(auto trackPtBin : comparedTrackPtBin){
              iTrackPt = card[0][0][0][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
              for(auto centralityBin: comparedCentralityBin){
                iCentrality = centralityBin.first < 0 ? 0 : card[0][0][0][0][0]->FindBinIndexCentrality(centralityBin);
              
                // First calculate the relative edge loss

                // Find the bins corresponding to the defined region from which to calculate edge loss
                firstEdgeLossBin = hEnergyEnergyCorrelatorRadiusRatio[iSimulation][iVeto][kR0p4][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->FindBin(edgeLossRegion.first);
                lastEdgeLossBin = hEnergyEnergyCorrelatorRadiusRatio[iSimulation][iVeto][kR0p4][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->FindBin(edgeLossRegion.second);
                edgeLossLowerBoundary = hEnergyEnergyCorrelatorRadiusRatio[iSimulation][iVeto][kR0p4][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->GetXaxis()->GetBinLowEdge(firstEdgeLossBin);
                edgeLossUpperBoundary = hEnergyEnergyCorrelatorRadiusRatio[iSimulation][iVeto][kR0p4][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->GetXaxis()->GetBinUpEdge(lastEdgeLossBin);

                // Calculate the total area and the edge loss integral
                totalArea = edgeLossUpperBoundary - edgeLossLowerBoundary;
                edgeLossIntegral = hEnergyEnergyCorrelatorRadiusRatio[iSimulation][iVeto][kR0p4][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->IntegralAndError(firstEdgeLossBin, lastEdgeLossBin, edgeLossIntegralError, "width");

                // From these, we can determine the relative edge loss
                edgeLossAmount[iSimulation][kRelativeEdgeLoss][iVeto][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt] = totalArea - edgeLossIntegral;
                edgeLossError[iSimulation][kRelativeEdgeLoss][iVeto][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt] = edgeLossIntegralError;


                // Next, calculate the absolute edge loss
                totalArea = hEnergyEnergyCorrelator[iSimulation][iVeto][kR0p8][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->IntegralAndError(firstEdgeLossBin, lastEdgeLossBin, totalAreaError, "width");
                edgeLossIntegral = hEnergyEnergyCorrelator[iSimulation][iVeto][kR0p4][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->IntegralAndError(firstEdgeLossBin, lastEdgeLossBin, edgeLossIntegralError, "width");

                // Absolute energy loss is the difference of integrals between R=0.8 and R=0.4 in the defined region
                edgeLossAmount[iSimulation][kAbsoluteEdgeLoss][iVeto][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt] = totalArea - edgeLossIntegral;
                edgeLossError[iSimulation][kAbsoluteEdgeLoss][iVeto][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt] = TMath::Sqrt(TMath::Power(totalAreaError,2) + TMath::Power(edgeLossIntegralError,2));
              } // Centrality loop
            } // Track pT loop
          } // Jet pT loop
        } // Jet radius loop
      } // Energy weight loop
    } // Jet veto loop
  } // Simulation loop

  // After the edge losses have been calculated, create histograms with the edge loss as a function and jet pT, and average DeltaR between E-scheme and WTA jet axes
  double edgeLossForThisGraph[nJetPtBins];
  double edgeLossErrorForThisGraph[nJetPtBins];
  double deltaJetAxisForThisGraph[nJetPtBins];
  double deltaJetAxisErrorForThisGraph[nJetPtBins];
  for(int iSimulation = 0; iSimulation < knSimulations; iSimulation++){
    for(int iVeto = 0; iVeto < knVetoTypes; iVeto++){
      if(iSimulation == kHerwig && iVeto == kJetVeto) continue; // Note: We have no files with jet veto for Herwig
      for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
        for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
          for(auto trackPtBin : comparedTrackPtBin){
            iTrackPt = card[0][0][0][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
            for(auto centralityBin: comparedCentralityBin){
              iCentrality = centralityBin.first < 0 ? 0 : card[0][0][0][0][0]->FindBinIndexCentrality(centralityBin);

              // Make a new histogram with the x-axis defined by the determined jet pT bins
              hEdgeLossVsPt[iSimulation][iVeto][iEnergyWeight][iJetAxis][iCentrality][iTrackPt] = new TH1D(Form("edgeLoassVsPt%d%d%d%d%d%d", iSimulation, iVeto, iEnergyWeight, iJetAxis, iCentrality, iTrackPt), Form("edgeLossVsPt%d%d%d%d%d%d", iSimulation, iVeto, iEnergyWeight, iJetAxis, iCentrality, iTrackPt), nJetPtBins, jetPtBinBorderArray);
              hEdgeLossVsPt[iSimulation][iVeto][iEnergyWeight][iJetAxis][iCentrality][iTrackPt]->Sumw2();

              // Once the histogram is created, fill each bin with the edge loss amount
              iBin = 1;
              for(auto jetPtBin : comparedJetPtBin){
                iJetPt = card[0][0][0][0][0]->FindBinIndexJetPtEEC(jetPtBin);
                hEdgeLossVsPt[iSimulation][iVeto][iEnergyWeight][iJetAxis][iCentrality][iTrackPt]->SetBinContent(iBin, edgeLossAmount[iSimulation][edgeLossMethod][iVeto][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]);
                hEdgeLossVsPt[iSimulation][iVeto][iEnergyWeight][iJetAxis][iCentrality][iTrackPt]->SetBinError(iBin++, edgeLossError[iSimulation][edgeLossMethod][iVeto][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]);
              } 

              // Once we have the edge loss as a function of jet pT, we can map each jet pT bin onto an average DeltaR value between E-scheme and WTA axes and create a graph

              // Prepare the arrays for edge loss amount and DeltaR between different jet axes
              int iBin = 0;
              for(auto jetPtBin : comparedJetPtBin){
                iJetPt = card[0][0][0][0][0]->FindBinIndexJetPtEEC(jetPtBin);

                // First, collect the edge loss amount for this graph to an array
                edgeLossForThisGraph[iBin] = edgeLossAmount[iSimulation][edgeLossMethod][iVeto][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt];
                edgeLossErrorForThisGraph[iBin] = edgeLossError[iSimulation][edgeLossMethod][iVeto][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt];

                // Then, read the average DeltaR between E-scheme and WTA axes from the histograms
                deltaJetAxisForThisGraph[iBin] = hAverageDeltaJetAxisVsPt[iSimulation][iJetAxis][iCentrality]->GetBinContent(iBin+1);
                deltaJetAxisErrorForThisGraph[iBin] = hAverageDeltaJetAxisVsPt[iSimulation][iJetAxis][iCentrality]->GetBinError(iBin+1);

                // Increment the binning variable
                iBin++;

              }

              // After the arrays are prepared, we can create graphs using those
              hEdgeLossVsDeltaJetAxis[iSimulation][iVeto][iEnergyWeight][iJetAxis][iCentrality][iTrackPt] = new TGraphErrors(nJetPtBins, deltaJetAxisForThisGraph, edgeLossForThisGraph, deltaJetAxisErrorForThisGraph, edgeLossErrorForThisGraph);

              // Option to fit the edge loss vs. Delta jet axis distribution
              if(fitLinearToEdgeLoss){
                hEdgeLossVsDeltaJetAxis[iSimulation][iVeto][iEnergyWeight][iJetAxis][iCentrality][iTrackPt]->Fit(fLinearFit[iSimulation][iVeto][iEnergyWeight][iJetAxis][iCentrality][iTrackPt], "N");
              }

            } // Centrality loop
          } // Track pT loop
        } // Jet radius loop
      } // Energy weight loop
    } // Jet veto loop
  } // Simulation loop
  
  // ====================================================================
  //                Drawing the selected distributions
  // ====================================================================
  
  JDrawer* drawer = new JDrawer();
  drawer->SetRelativeCanvasSize(1.1,1.1);

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
  TLegend* edgeLossLegend;
  TLine* oneLine = new TLine(drawingRange.first, 1, drawingRange.second, 1);
  oneLine->SetLineStyle(2);

  // Variables needed to use the histogram drawing function
  std::vector<std::pair<TH1D*,TString>> histogramVector;
  TString histogramLegendHeader;
  std::vector<TString> systemLegendVector;
  TString saveName = "";

  // Draw the average DeltaR between E-scheme and WTA axes as a function of jet pT
  if(drawDeltaJetAxisWithPt){
    for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
      for(auto centralityBin : comparedCentralityBin){
        iCentrality = centralityBin.first < 0 ? 0 : card[0][0][0][0][0]->FindBinIndexCentrality(centralityBin);

        // Create the legend and add binning information to it
        legend = new TLegend(0.53,0.68,0.93,0.84);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);

        if(centralityBin.first >= 0) legend->AddEntry((TObject*) 0, Form("Centrality: %d-%d%%", centralityBin.first, centralityBin.second), "");
        //legend->AddEntry((TObject*) 0, Form("Jet selection: %s", jetAxisDescription[iJetAxis].Data()), "");

        // Setup centrality strings
        compactCentralityString = centralityBin.first == -1 ? "_pp" : Form("_C=%d-%d", centralityBin.first, centralityBin.second);

        // Set drawing style for the histograms
        for(int iSimulation = 0; iSimulation < knSimulations; iSimulation++){
          hAverageDeltaJetAxisVsPt[iSimulation][iJetAxis][iCentrality]->SetMarkerStyle(markerStyle[iSimulation]);
          hAverageDeltaJetAxisVsPt[iSimulation][iJetAxis][iCentrality]->SetMarkerColor(color[iSimulation]);
          hAverageDeltaJetAxisVsPt[iSimulation][iJetAxis][iCentrality]->SetLineColor(color[iSimulation]);
        }

        // Set drawing ranges for the axes
        hAverageDeltaJetAxisVsPt[kPythia][iJetAxis][iCentrality]->GetYaxis()->SetRangeUser(deltaJetAxisZoom.first, deltaJetAxisZoom.second);

        // Draw the histogram
        drawer->DrawHistogram(hAverageDeltaJetAxisVsPt[kPythia][iJetAxis][iCentrality], "Jet p_{T} (GeV)", "#LT#phi#GT", " ", "p");
        hAverageDeltaJetAxisVsPt[kHerwig][iJetAxis][iCentrality]->Draw("p,same");

        // Add legend for the histograms
        for(int iSimulation = 0; iSimulation < knSimulations; iSimulation++){
          legend->AddEntry(hAverageDeltaJetAxisVsPt[iSimulation][iJetAxis][iCentrality], simulationName[iSimulation].Data(), "p");
        }

        if(fitOneOverPt){
          for(int iSimulation = 0; iSimulation < knSimulations; iSimulation++){
            fOneOverPtFit[iSimulation][iJetAxis][iCentrality]->Draw("same");
            legend->AddEntry(fOneOverPtFit[iSimulation][iJetAxis][iCentrality], Form("%.2f + %.2f/p_{T}", fOneOverPtFit[iSimulation][iJetAxis][iCentrality]->GetParameter(0), fOneOverPtFit[iSimulation][iJetAxis][iCentrality]->GetParameter(1)), "l");
          }
        }

        // Draw the legend
        legend->Draw();

        // Save the figures to a file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/averageDeltaJetAxis%s%s%s.%s", saveComment.Data(), jetAxisSaveName[iJetAxis].Data(), compactCentralityString.Data(), figureFormat));
        }

      } // Centrality loop
    } // Jet axis type loop
  } // Drawing average DeltaR between E-scheme and WTA axes as a function of jet pT

  

  // Draw illustration plots for edge loss
  if(drawEdgeLossIllustration){

    for(int iSimulation = 0; iSimulation < knSimulations; iSimulation++){
      for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
        for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
          for(auto trackPtBin : comparedTrackPtBin){
            iTrackPt = card[0][0][0][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
            for(auto centralityBin : comparedCentralityBin){
              iCentrality = centralityBin.first < 0 ? 0 : card[0][0][0][0][0]->FindBinIndexCentrality(centralityBin);

              // Create vectors for histograms to be drawn
              histogramVector.clear();
              for(auto jetPtBin : comparedJetPtBin){
                iJetPt = card[0][0][0][0][0]->FindBinIndexJetPtEEC(jetPtBin);
                histogramVector.push_back(std::make_pair(hEnergyEnergyCorrelatorRadiusRatio[iSimulation][vetoCloseJets][kR0p4][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt], Form("%d-%d GeV", jetPtBin.first, jetPtBin.second)));
              }

              // Add the information we want to draw legends
              histogramLegendHeader = "Jet p_{T}";
              systemLegendVector.clear();
              systemLegendVector.push_back(simulationName[iSimulation]);
              systemLegendVector.push_back(Form("Jet axis: %s", jetAxisDescription[iJetAxis].Data()));
              if(vetoDescription[vetoCloseJets] != "") systemLegendVector.push_back(vetoDescription[vetoCloseJets]);
              if(centralityBin.first >= 0) systemLegendVector.push_back(Form("Centrality: %d-%d%%", centralityBin.first, centralityBin.second));
              systemLegendVector.push_back(Form("n = %d", iEnergyWeight+1));
              systemLegendVector.push_back(Form("p_{T}^{ch} > %.1f GeV", trackPtBin));

              // Create a good name if the histogram is saved
              if(saveFigures){
                // Setup centrality and track pT strings
                compactTrackPtString = Form("_T>%.1f",trackPtBin);
                compactTrackPtString.ReplaceAll(".","v");
                compactCentralityString = centralityBin.first == -1 ? "_pp" : Form("_C=%d-%d", centralityBin.first, centralityBin.second);

                saveName = Form("figures/edgeLossIllustration%s%s%s%s%s%s.%s", saveComment.Data(), simulationName[iSimulation].Data(), jetAxisSaveName[iJetAxis].Data(), compactEnergyWeightString[iEnergyWeight].Data(), compactCentralityString.Data(), compactTrackPtString.Data(), figureFormat);
              }

              // Now that we have histogram drawing information prepared, we can draw the histograms to canvas
              drawVector(drawer, histogramVector, histogramLegendHeader, systemLegendVector, "EEC #frac{R = 0.4}{R = 0.8}", drawingRange, edgeLossIllustrationZoom, saveName);

            } // Centrality looop
          } // Track pT loop
        } // Jet axis type loop
      } // Energy weight loop
    } // Simulation loop
  } // Drawing edge loss illustration

  // Draw illustration plots to show how vetoing jets with other close jets biases the distributions. Only done for Pythia
  if(drawJetVetoBiasIllustration){

    for(int iJetRadius = 0; iJetRadius < knJetRadii; iJetRadius++){
      for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
        for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
          for(auto trackPtBin : comparedTrackPtBin){
            iTrackPt = card[0][0][0][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
            for(auto centralityBin : comparedCentralityBin){
              iCentrality = centralityBin.first < 0 ? 0 : card[0][0][0][0][0]->FindBinIndexCentrality(centralityBin);

              // Create vectors for histograms to be drawn
              histogramVector.clear();
              for(auto jetPtBin : comparedJetPtBin){
                iJetPt = card[0][0][0][0][0]->FindBinIndexJetPtEEC(jetPtBin);
                histogramVector.push_back(std::make_pair(hEnergyEnergyCorrelatorVetoBiasCheck[kPythia][kJetVeto][iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt], Form("%d-%d GeV", jetPtBin.first, jetPtBin.second)));
              }

              // Add the information we want to draw legends
              histogramLegendHeader = "Jet p_{T}";
              systemLegendVector.clear();
              systemLegendVector.push_back(simulationName[kPythia].Data());
              systemLegendVector.push_back(Form("Jet axis: %s", jetAxisDescription[iJetAxis].Data()));
              if(centralityBin.first >= 0) systemLegendVector.push_back(Form("Centrality: %d-%d%%", centralityBin.first, centralityBin.second));
              systemLegendVector.push_back(Form("n = %d", iEnergyWeight+1));
              systemLegendVector.push_back(Form("p_{T}^{ch} > %.1f GeV", trackPtBin));

              // Create a good name if the histogram is saved
              if(saveFigures){
                // Setup centrality and track pT strings
                compactTrackPtString = Form("_T>%.1f",trackPtBin);
                compactTrackPtString.ReplaceAll(".","v");
                compactCentralityString = centralityBin.first == -1 ? "_pp" : Form("_C=%d-%d", centralityBin.first, centralityBin.second);

                saveName = Form("figures/vetoBiasCheck%s%s%s%s%s%s.%s", saveComment.Data(), jetAxisSaveName[iJetAxis].Data(), jetRadiusSaveName[iJetRadius].Data(), compactEnergyWeightString[iEnergyWeight].Data(), compactCentralityString.Data(), compactTrackPtString.Data(), figureFormat);
              }

              // Now that we have histogram drawing information prepared, we can draw the histograms to canvas
              drawVector(drawer, histogramVector, histogramLegendHeader, systemLegendVector, Form("EEC #frac{%s}{%s}", vetoDescription[kJetVeto].Data(), jetRadiusDescription[iJetRadius].Data()), drawingRange, edgeLossIllustrationZoom, saveName);

            } // Centrality looop
          } // Track pT loop
        } // Jet axis type loop
      } // Energy weight loop
    } // Jet radius loop
  } // Drawing bias check due to jet veto

  // Draw the edge loss as a function of jet pT
  if(drawEdgeLossWithPt){
    for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
      for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
        for(auto trackPtBin : comparedTrackPtBin){
          iTrackPt = card[0][0][0][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
          compactTrackPtString = Form("_T>%.1f",trackPtBin);
          compactTrackPtString.ReplaceAll(".","v");
          for(auto centralityBin: comparedCentralityBin){
            iCentrality = centralityBin.first < 0 ? 0 : card[0][0][0][0][0]->FindBinIndexCentrality(centralityBin);

            // Create the legend and add binning information to it
            legend = new TLegend(0.35,0.4,0.9,0.75);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.045);legend->SetTextFont(62);

            legend->AddEntry((TObject*) 0, Form("Jet selection: %s", jetAxisDescription[iJetAxis].Data()), "");
            if(vetoDescription[vetoCloseJets] != "") legend->AddEntry((TObject*) 0, vetoDescription[vetoCloseJets].Data(), "");
            if(centralityBin.first >= 0) legend->AddEntry((TObject*) 0, Form("Centrality: %d-%d%%", centralityBin.first, centralityBin.second), "");
            legend->AddEntry((TObject*) 0, Form("n = %d", iEnergyWeight+1), "");
            legend->AddEntry((TObject*) 0, Form("p_{T}^{ch} > %.1f GeV", trackPtBin), "");

            edgeLossLegend = new TLegend(0.15,0.8,0.55,0.92);
            edgeLossLegend->SetFillStyle(0);edgeLossLegend->SetBorderSize(0);
            edgeLossLegend->SetTextSize(0.045);edgeLossLegend->SetTextFont(62);
 
            edgeLossLegend->AddEntry((TObject*) 0, Form("Edge loss region: %.3f < #Deltar < %.3f", edgeLossLowerBoundary, edgeLossUpperBoundary), "");
            edgeLossLegend->AddEntry((TObject*) 0, Form("Normalization: %s", normalizationDescription[normalizeDistributions].Data()), "");

            // Setup centrality strings
            compactCentralityString = centralityBin.first == -1 ? "_pp" : Form("_C=%d-%d", centralityBin.first, centralityBin.second);

            // Set drawing style for the histograms
            for(int iSimulation = 0; iSimulation < knSimulations; iSimulation++){
              hEdgeLossVsPt[iSimulation][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt]->SetMarkerStyle(markerStyle[iSimulation]);
              hEdgeLossVsPt[iSimulation][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt]->SetMarkerColor(color[iSimulation]);
              hEdgeLossVsPt[iSimulation][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt]->SetLineColor(color[iSimulation]);
            }

            // Set drawing ranges for the axes
            hEdgeLossVsPt[kPythia][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt]->GetYaxis()->SetRangeUser(edgeLossZoom[edgeLossMethod].first, edgeLossZoom[edgeLossMethod].second);

            // Draw the histogram
            drawer->DrawHistogram(hEdgeLossVsPt[kPythia][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt], "Jet p_{T} (GeV)", Form("%s edge loss", edgeLossString[edgeLossMethod].Data()), " ", "p");
            hEdgeLossVsPt[kHerwig][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt]->Draw("same,p");

            // Add histograms to the legend
            for(int iSimulation = 0; iSimulation < knSimulations; iSimulation++){
              legend->AddEntry(hEdgeLossVsPt[iSimulation][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt], simulationName[iSimulation].Data(), "p");
            }

            // Draw the legend
            legend->Draw();
            edgeLossLegend->Draw();

            // Save the figures to a file
            if(saveFigures){
              gPad->GetCanvas()->SaveAs(Form("figures/edgeLossWithPt%s%s%s%s%s%s.%s", saveComment.Data(), vetoSaveName[vetoCloseJets].Data(), jetAxisSaveName[iJetAxis].Data(), energyWeightDescription[iEnergyWeight].Data(), compactCentralityString.Data(), compactTrackPtString.Data(), figureFormat));
            }
          } // Centrality loop
        } // Track pT loop
      } // Jet axis loop
    } // Energy weight loop
  } // Drawing the edge loss as a function of jet pT


  // Draw the edge loss as a function of DeltaR between E-scheme and WTA jet axes
  if(drawEdgeLossWithDeltaJetAxis){
    for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
      for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
        for(auto trackPtBin : comparedTrackPtBin){
          iTrackPt = card[0][0][0][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
          compactTrackPtString = Form("_T>%.1f",trackPtBin);
          compactTrackPtString.ReplaceAll(".","v");
          for(auto centralityBin: comparedCentralityBin){
            iCentrality = centralityBin.first < 0 ? 0 : card[0][0][0][0][0]->FindBinIndexCentrality(centralityBin);

            // Create the legend and add binning information to it
            legend = new TLegend(0.35,0.4,0.9,0.75);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.045);legend->SetTextFont(62);

            legend->AddEntry((TObject*) 0, Form("Jet selection: %s", jetAxisDescription[iJetAxis].Data()), "");
            if(vetoDescription[vetoCloseJets] != "") legend->AddEntry((TObject*) 0, vetoDescription[vetoCloseJets].Data(), "");
            if(centralityBin.first >= 0) legend->AddEntry((TObject*) 0, Form("Centrality: %d-%d%%", centralityBin.first, centralityBin.second), "");
            legend->AddEntry((TObject*) 0, Form("n = %d", iEnergyWeight+1), "");
            legend->AddEntry((TObject*) 0, Form("p_{T}^{ch} > %.1f GeV", trackPtBin), "");

            edgeLossLegend = new TLegend(0.15,0.8,0.55,0.92);
            edgeLossLegend->SetFillStyle(0);edgeLossLegend->SetBorderSize(0);
            edgeLossLegend->SetTextSize(0.045);edgeLossLegend->SetTextFont(62);
 
            edgeLossLegend->AddEntry((TObject*) 0, Form("Edge loss region: %.3f < #Deltar < %.3f", edgeLossLowerBoundary, edgeLossUpperBoundary), "");
            edgeLossLegend->AddEntry((TObject*) 0, Form("Normalization: %s", normalizationDescription[normalizeDistributions].Data()), "");

            // Setup centrality strings
            compactCentralityString = centralityBin.first == -1 ? "_pp" : Form("_C=%d-%d", centralityBin.first, centralityBin.second);

            // Set drawing style for the graphs
            for(int iSimulation = 0; iSimulation < knSimulations; iSimulation++){
              hEdgeLossVsDeltaJetAxis[iSimulation][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt]->SetMarkerStyle(markerStyle[iSimulation]);
              hEdgeLossVsDeltaJetAxis[iSimulation][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt]->SetMarkerColor(color[iSimulation]);
              hEdgeLossVsDeltaJetAxis[iSimulation][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt]->SetLineColor(color[iSimulation]);
              if(fitLinearToEdgeLoss){
                fLinearFit[iSimulation][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt]->SetLineColor(color[iSimulation]);
              }
            }

            // Draw the histogram
            drawer->DrawGraph(hEdgeLossVsDeltaJetAxis[kPythia][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt], 0.025, 0.055, 0, 0.03, "#LT#phi#GT", Form("%s edge loss", edgeLossString[edgeLossMethod].Data()), " ", "p");
            hEdgeLossVsDeltaJetAxis[kHerwig][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt]->Draw("same,p");

            // If reguested, draw the fits to the same plot
            if(fitLinearToEdgeLoss){
              for(int iSimulation = 0; iSimulation < knSimulations; iSimulation++){
                fLinearFit[iSimulation][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt]->Draw("same");
              }
            }

            // Add graphs to the legend
            for(int iSimulation = 0; iSimulation < knSimulations; iSimulation++){ 
              if(fitLinearToEdgeLoss){
                legend->AddEntry(hEdgeLossVsDeltaJetAxis[iSimulation][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt], Form("%s (#chi^{2}/ndf = %.3f)", simulationName[iSimulation].Data(), fLinearFit[iSimulation][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt]->GetChisquare() / fLinearFit[iSimulation][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt]->GetNDF()), "p");
                legend->AddEntry(fLinearFit[iSimulation][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt], Form("(%.4f #pm %.4f)", fLinearFit[iSimulation][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt]->GetParameter(0), fLinearFit[iSimulation][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt]->GetParError(0)), "l");
                legend->AddEntry(fLinearFit[iSimulation][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt], Form("+ (%.4f #pm %.4f) #LT#phi#GT", fLinearFit[iSimulation][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt]->GetParameter(1), fLinearFit[iSimulation][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt]->GetParError(1)), "l");
              } else {
                legend->AddEntry(hEdgeLossVsDeltaJetAxis[iSimulation][vetoCloseJets][iEnergyWeight][iJetAxis][iCentrality][iTrackPt], simulationName[iSimulation].Data(), "p");
              }
            }

            // Draw the legend
            legend->Draw();
            edgeLossLegend->Draw();

            // Save the figures to a file
            if(saveFigures){
              gPad->GetCanvas()->SaveAs(Form("figures/edgeLossWithDeltaJetAxis%s%s%s%s%s%s.%s", saveComment.Data(), vetoSaveName[vetoCloseJets].Data(), jetAxisSaveName[iJetAxis].Data(), energyWeightDescription[iEnergyWeight].Data(), compactCentralityString.Data(), compactTrackPtString.Data(), figureFormat));
            }
          } // Centrality loop
        } // Track pT loop
      } // Jet axis loop
    } // Energy weight loop
  } // Drawing the edge loss as a function of DeltaR between E-scheme and WTA jet axes


  // Draw an illustrative plot that shows how edge loss is defined
  if(drawEdgeLossDefinition){

    // Set a good drawing style for split canvas
    drawer->SetDefaultAppearanceSplitCanvas();
    drawer->SetRelativeCanvasSize(1.1,1.1);
    drawer->SetLeftMargin(0.14);
    drawer->SetTopMargin(0.07);
    drawer->SetTitleOffsetY(1.7);
    drawer->SetTitleOffsetX(1.0);
    drawer->SetLogY(true);

    // Use one example jet pT and track pT selection for the illustration
    std::pair<int,int> exampleJetPt = std::make_pair(120,140);
    double exampleTrackPt = 1;

    iJetPt = card[kPythia][vetoCloseJets][kR0p4][kNominalWeight][kEscheme]->FindBinIndexJetPtEEC(exampleJetPt);
    iTrackPt = card[kPythia][vetoCloseJets][kR0p4][kNominalWeight][kEscheme]->GetBinIndexTrackPtEEC(exampleTrackPt);

    // Create the legend and add binning information to it
    legend = new TLegend(0.18,0.04,0.45,0.58);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);

    legend->AddEntry((TObject*) 0, simulationName[kPythia].Data(), "");
    legend->AddEntry((TObject*) 0, "Jet axis: e-scheme", "");
    legend->AddEntry((TObject*) 0, Form("%d < jet p_{T} < %d GeV", exampleJetPt.first, exampleJetPt.second), "");
    legend->AddEntry((TObject*) 0, Form("p_{T}^{ch} > %.0f GeV", exampleTrackPt), "");
    legend->AddEntry((TObject*) 0, "n = 1", "");

    // Set drawing style for the histograms
    hEnergyEnergyCorrelator[kPythia][vetoCloseJets][kR0p4][kNominalWeight][kEscheme][0][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[0]);
    hEnergyEnergyCorrelator[kPythia][vetoCloseJets][kR0p4][kNominalWeight][kEscheme][0][iJetPt][iTrackPt]->SetMarkerColor(color[0]);
    hEnergyEnergyCorrelator[kPythia][vetoCloseJets][kR0p4][kNominalWeight][kEscheme][0][iJetPt][iTrackPt]->SetLineColor(color[0]);
    hEnergyEnergyCorrelator[kPythia][vetoCloseJets][kR0p8][kNominalWeight][kEscheme][0][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[1]);
    hEnergyEnergyCorrelator[kPythia][vetoCloseJets][kR0p8][kNominalWeight][kEscheme][0][iJetPt][iTrackPt]->SetMarkerColor(color[1]);
    hEnergyEnergyCorrelator[kPythia][vetoCloseJets][kR0p8][kNominalWeight][kEscheme][0][iJetPt][iTrackPt]->SetLineColor(color[1]);

    // Set the x- and y-axis drawing ranges
    hEnergyEnergyCorrelator[kPythia][vetoCloseJets][kR0p4][kNominalWeight][kEscheme][0][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.1, 0.39);

    // Draw the histograms to the upper pad
    drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelator[kPythia][vetoCloseJets][kR0p4][kNominalWeight][kEscheme][0][iJetPt][iTrackPt], "#Deltar", "EEC", " ", "p");
    hEnergyEnergyCorrelator[kPythia][vetoCloseJets][kR0p8][kNominalWeight][kEscheme][0][iJetPt][iTrackPt]->Draw("same,p");

    // Add histograms to the legend
    legend->AddEntry(hEnergyEnergyCorrelator[kPythia][vetoCloseJets][kR0p4][kNominalWeight][kEscheme][0][iJetPt][iTrackPt], "R = 0.4", "p");
    legend->AddEntry(hEnergyEnergyCorrelator[kPythia][vetoCloseJets][kR0p8][kNominalWeight][kEscheme][0][iJetPt][iTrackPt], "R = 0.8", "p");

    // Draw the legend
    legend->Draw();

    // Linear scale for the ratio
    drawer->SetLogY(false);
          
    // Set the drawing style ranges
    hEnergyEnergyCorrelatorRadiusRatio[kPythia][vetoCloseJets][kR0p4][kNominalWeight][kEscheme][0][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[0]);
    hEnergyEnergyCorrelatorRadiusRatio[kPythia][vetoCloseJets][kR0p4][kNominalWeight][kEscheme][0][iJetPt][iTrackPt]->SetMarkerColor(color[0]);
    hEnergyEnergyCorrelatorRadiusRatio[kPythia][vetoCloseJets][kR0p4][kNominalWeight][kEscheme][0][iJetPt][iTrackPt]->SetLineColor(color[0]);
    hEnergyEnergyCorrelatorRadiusRatio[kPythia][vetoCloseJets][kR0p4][kNominalWeight][kEscheme][0][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.1, 0.39);
    hEnergyEnergyCorrelatorRadiusRatio[kPythia][vetoCloseJets][kR0p4][kNominalWeight][kEscheme][0][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.8, 1.04);

    // Draw the histograms
    drawer->SetGridY(true);

    drawer->DrawHistogramToLowerPad(hEnergyEnergyCorrelatorRadiusRatio[kPythia][vetoCloseJets][kR0p4][kNominalWeight][kEscheme][0][iJetPt][iTrackPt], "#Deltar", "#frac{R = 0.4}{R = 0.8}", " ", "p");

    drawer->SetGridY(false);

  } // Drawing illustration of edge loss definition

  // Draw the edge loss as a function of DeltaR between E-scheme and WTA jet axes
  if(drawJetDeltaAxis){     
    for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
      for(auto centralityBin : comparedCentralityBin){
        iCentrality = centralityBin.first < 0 ? 0 : card[0][0][0][0][0]->FindBinIndexCentrality(centralityBin);
        for(auto jetPtBin : comparedJetPtBin){
          iJetPt = card[0][0][0][0][0]->FindBinIndexJetPtEEC(jetPtBin);

          // Create the legend and add binning information to it
          legend = new TLegend(0.35,0.4,0.9,0.75);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.045);legend->SetTextFont(62);

          //legend->AddEntry((TObject*) 0, Form("Jet selection: %s", jetAxisDescription[iJetAxis].Data()), "");
          if(vetoDescription[vetoCloseJets] != "") legend->AddEntry((TObject*) 0, vetoDescription[vetoCloseJets].Data(), "");
          if(centralityBin.first >= 0) legend->AddEntry((TObject*) 0, Form("Centrality: %d-%d%%", centralityBin.first, centralityBin.second), "");
          legend->AddEntry((TObject*) 0, Form("%d < jet p_{T} < %d", jetPtBin.first, jetPtBin.second), "");

          // Setup centrality strings
          compactCentralityString = centralityBin.first == -1 ? "_pp" : Form("_C=%d-%d", centralityBin.first, centralityBin.second);
          compactJetPtString = Form("_J=%d-%d", jetPtBin.first, jetPtBin.second);

          // Set drawing style for the histograms
          for(int iSimulation = 0; iSimulation < knSimulations; iSimulation++){
            hDeltaJetAxis[iSimulation][iJetAxis][iCentrality][iJetPt]->SetMarkerStyle(markerStyle[iSimulation]);
            hDeltaJetAxis[iSimulation][iJetAxis][iCentrality][iJetPt]->SetMarkerColor(color[iSimulation]);
            hDeltaJetAxis[iSimulation][iJetAxis][iCentrality][iJetPt]->SetLineColor(color[iSimulation]);

            // Normalize distributions to one
            hDeltaJetAxis[iSimulation][iJetAxis][iCentrality][iJetPt]->Scale(1.0 / hDeltaJetAxis[iSimulation][iJetAxis][iCentrality][iJetPt]->Integral("width"));

          }

          // Set drawing ranges for the axes
          //hDeltaJetAxis[kPythia][iJetAxis][iCentrality][iJetPt]->GetYaxis()->SetRangeUser(edgeLossZoom[edgeLossMethod].first, edgeLossZoom[edgeLossMethod].second);

          // Draw the histogram
          drawer->DrawHistogram(hDeltaJetAxis[kPythia][iJetAxis][iCentrality][iJetPt], "#phi", "#frac{dN}{d#phi}", " ", "p");
          hDeltaJetAxis[kHerwig][iJetAxis][iCentrality][iJetPt]->Draw("same,p");

          // Add histograms to the legend
          for(int iSimulation = 0; iSimulation < knSimulations; iSimulation++){
            legend->AddEntry(hDeltaJetAxis[iSimulation][iJetAxis][iCentrality][iJetPt], simulationName[iSimulation].Data(), "p");
          }

          // Draw the legend
          legend->Draw();

          // Save the figures to a file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/deltaJetAxis%s%s%s%s.%s", saveComment.Data(), vetoSaveName[vetoCloseJets].Data(), jetAxisSaveName[iJetAxis].Data(), compactCentralityString.Data(), figureFormat));
          }
        } // Jet pT loop
      } // Centrality loop
    } // Jet axis loop
  } // Drawing the edge loss as a function of DeltaR between E-scheme and WTA jet axes
  
}

// Common drawing macro for any kind of ratio plots to avoid copy-pasting code
//
//  JDrawer* drawer = Histogram drawing class;
//  std::vector<std::pair<TH1D*,TString>> histogramVector = Vector with all the histograms that we want to draw and a name for legend
//  TString histogramLegendHeader = Header given for the histogram legend
//  std::vector<TString> systemLegendVector = Everything to be added to the system legend of the plot
//  TString axisName = Name given to the y-axis of the plot
//  std::pair<double, double> xRange = Drawing range for the x-axis
//  std::pair<double, double> yRange = Drawing range for the y-axis
//  TString saveName = Name given to saved histograms
void drawVector(JDrawer* drawer, std::vector<std::pair<TH1D*,TString>> histogramVector, TString histogramLegendHeader, std::vector<TString> systemLegendVector, TString axisName, std::pair<double, double> xRange, std::pair<double, double> yRange, TString saveName){

  // Find the unmber of histograms to draw
  const int nHistograms = histogramVector.size();

  // Automatic color selection for all different jet pT bins
  gStyle->SetPalette(kBird);
  auto niceColors = TColor::GetPalette();
  double step = 254.0/nHistograms;
  int histogramColor[nHistograms];
  for(int iHistogram = 0; iHistogram < nHistograms; iHistogram++){
    histogramColor[iHistogram] = niceColors.At(TMath::Ceil(iHistogram*step));
  }

  // Logarithmic DeltaR axis
  drawer->SetLogX(true);

  // Create the legend and add binning information to it
  TLegend* legend = new TLegend(0.2,0.18,0.45,0.46);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.045);legend->SetTextFont(62);

  for(auto legendItem : systemLegendVector){
    legend->AddEntry((TObject*) 0, legendItem.Data(), "");
  }

  // Create a legend to hold jet pT bins
  TLegend* histogramLegend = new TLegend(0.2,0.6,0.92,0.9);
  histogramLegend->SetFillStyle(0);histogramLegend->SetBorderSize(0);
  histogramLegend->SetTextSize(0.03);histogramLegend->SetTextFont(62);
  histogramLegend->SetNColumns(3);
  histogramLegend->SetHeader(histogramLegendHeader.Data(),"C");

  // Set drawing style for all histograms
  int iBin = 0;
  for(auto histogramPair : histogramVector){
    histogramPair.first->SetLineColor(histogramColor[iBin++]);
  } // Jet pT loop

  // Set the x- and y-axis drawing ranges
  histogramVector.at(0).first->GetXaxis()->SetRangeUser(xRange.first, xRange.second);
  histogramVector.at(0).first->GetYaxis()->SetRangeUser(yRange.first, yRange.second);

  // Draw the histograms to canvas
  bool isFirst = true;
  for(auto histogramPair : histogramVector){
    if(isFirst){
      drawer->DrawHistogram(histogramPair.first, "#Deltar", axisName.Data(), " ", "l");
      isFirst = false;
    } else {
      histogramPair.first->Draw("same,l");
    }

    // Add a description of the histogram to the histogram legend
    histogramLegend->AddEntry(histogramPair.first, histogramPair.second.Data(), "l");
  }

  // Draw the legends
  legend->Draw();
  histogramLegend->Draw();
  
  // Save the figures
  if(saveName != ""){
    gPad->GetCanvas()->SaveAs(saveName.Data());
  }

}