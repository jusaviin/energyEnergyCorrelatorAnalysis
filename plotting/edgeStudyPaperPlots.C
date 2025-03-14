#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"
#include "../src/EECHistograms.h"

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

  // Define the file names for the studied files
  TString fileName[knJetRadii][knEnergyWeights][knJetAxes];

  // Naming for jet radii
  TString jetRadiusDescription[knJetRadii];
  jetRadiusDescription[kR0p4] = "R = 0.4";
  jetRadiusDescription[kR0p8] = "R = 0.8";

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
  

  // Define files for each of these combinations. TODO: Update all files after merging in ACCRE is finished

  // Files for R = 0.4
  fileName[kR0p4][kNominalWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_nominalEnergyWeight_jetDeltaAxis_leadingParticleFlag_processed_2025-03-10.root";
  fileName[kR0p4][kSquaredWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_energyWeightSquared_jetDeltaAxis_leadingParticleFlag_processed_2025-03-10.root";
  fileName[kR0p4][kNominalWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_nominalEnergyWeight_jetDeltaAxis_leadingParticleFlag_processed_2025-03-10.root";
  fileName[kR0p4][kSquaredWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_jetDeltaAxis_leadingParticleFlag_processed_2025-03-10.root";

  // Files for R = 0.8
  fileName[kR0p8][kNominalWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_nominalEnergyWeight_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_processed_2025-03-10.root";
  fileName[kR0p8][kSquaredWeight][kEscheme] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_eschemeAxis_energyWeightSquared_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_processed_2025-03-10.root";
  fileName[kR0p8][kNominalWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_nominalEnergyWeight_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_processed_2025-03-10.root";
  fileName[kR0p8][kSquaredWeight][kWTA] = "data/eschemeAxis/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_jetDeltaAxis_jetRadius0p8_leadingParticleFlag_processed_2025-03-10.root";
  
  // Open the files and check that they exist
  TFile* inputFile[knJetRadii][knEnergyWeights][knJetAxes];
  EECCard* card[knJetRadii][knEnergyWeights][knJetAxes];
  for(int iJetRadius = 0; iJetRadius < knJetRadii; iJetRadius++){
    for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
      for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
    
        inputFile[iJetRadius][iEnergyWeight][iJetAxis] = TFile::Open(fileName[iJetRadius][iEnergyWeight][iJetAxis]);
    
        if(inputFile[iJetRadius][iEnergyWeight][iJetAxis] == NULL){
          cout << "Error! The file " << fileName[iJetRadius][iEnergyWeight][iJetAxis].Data() << " does not exist!" << endl;
          cout << "Maybe you forgot the data/ folder path?" << endl;
          cout << "Will not execute the code" << endl;
          return;
        }

        card[iJetRadius][iEnergyWeight][iJetAxis]  = new EECCard(inputFile[iJetRadius][iEnergyWeight][iJetAxis]);
      } // Jet axis loop
    } // Energy weight loop
  } // System loop
  
  // ====================================================
  //               Binning configuration
  // ====================================================
  
  // Find the number of bins from the reference card
  const int nJetPtBinsEEC = card[0][0][0]->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card[0][0][0]->GetNTrackPtBinsEEC();
  const int nCentralityBins = card[0][0][0]->GetNCentralityBins();
  
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

  int nJetPtBins = comparedJetPtBin.size();

  std::vector<double> comparedTrackPtBin;
  comparedTrackPtBin.push_back(1.0);
  //comparedTrackPtBin.push_back(2.0);
  comparedTrackPtBin.push_back(3.0);

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

  // Flag for normalizing the energy-energy correlator distributions to unity
  bool normalizeDistributions = true;

  // ====================================================
  //                Drawing configuration
  // ====================================================

  // Paper plot selection
  bool drawDeltaJetAxisWithPt = false;        // Draw DeltaR between E-scheme and WTA jet axes as a function of pT
  bool fitOneOverPt = false;                  // Fit 1/pT function to average DeltaR between E-scheme and WTA jet axes as a function of pT histograms
  bool drawEdgeLossIllustration = true;      // Draw an illustration about edge loss
  bool drawEdgeLossWithPt = false;            // Draw the amount of edge-loss as a function of jet pT
  bool drawEdgeLossWithDelteJetAxis = false;  // Draw the amount of edge-loss as a function of DeltaR between E-scheme and WTA jet axes 
  
  // Figure saving
  const bool saveFigures = false;  // Save figures
  TString saveComment = "_jetAxisDifference_signal";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures

  // Drawing configuration
  std::pair<double, double> deltaJetAxisZoom = std::make_pair(0.02, 0.05);
  std::pair<double, double> edgeLossIllustrationZoom = std::make_pair(0.75, 1.25);
  std::pair<double, double> edgeLossZoom = std::make_pair(0.8, 1.2);
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

        // Choose the energy-energy correlator and jet histograms to load
        histograms[iJetRadius][iEnergyWeight][iJetAxis]->SetLoadJetHistograms(true);
        histograms[iJetRadius][iEnergyWeight][iJetAxis]->SetLoadEnergyEnergyCorrelators(true);

        // Choose the bin ranges
        histograms[iJetRadius][iEnergyWeight][iJetAxis]->SetCentralityBinRange(0, card[iJetRadius][iEnergyWeight][iJetAxis]->GetNCentralityBins() - 1);
        histograms[iJetRadius][iEnergyWeight][iJetAxis]->SetJetPtBinRangeEEC(0, card[iJetRadius][iEnergyWeight][iJetAxis]->GetNJetPtBinsEEC() - 1);
        histograms[iJetRadius][iEnergyWeight][iJetAxis]->SetTrackPtBinRangeEEC(0, card[iJetRadius][iEnergyWeight][iJetAxis]->GetNTrackPtBinsEEC() - 1);
        histograms[iJetRadius][iEnergyWeight][iJetAxis]->SetLoadedPairingType(pairingType, true);

        // Load the histograms from the file
        histograms[iJetRadius][iEnergyWeight][iJetAxis]->LoadProcessedHistograms();
      } // Jet axis loop for histogram loading
    } // Energy weight loop for histogram loading
  } // Jet radius loop  for histogram loading

  // Energy-energy correlator histograms
  TH1D* hEnergyEnergyCorrelator[knJetRadii][knEnergyWeights][knJetAxes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hEnergyEnergyCorrelatorAxisRatio[knJetRadii][knEnergyWeights][knJetAxes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hEnergyEnergyCorrelatorRadiusRatio[knJetRadii][knEnergyWeights][knJetAxes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];

  // Histograms for DeltaR between E-scheme and WTA axes
  TH1D* hDeltaJetAxis[knJetAxes][nCentralityBins][nJetPtBinsEEC];
  TH1D* hAverageDeltaJetAxisVsPt[knJetAxes][nCentralityBins];

  // Fit functions
  TF1* fOneOverPtFit[knJetAxes][nCentralityBins];
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
      hAverageDeltaJetAxisVsPt[iJetAxis][iCentrality] = NULL;
      fOneOverPtFit[iJetAxis][iCentrality] = new TF1(Form("oneOverPtFit%d%d", iCentrality, iJetAxis), "[0] + [1]/x", comparedJetPtBin.at(0).first, comparedJetPtBin.at(comparedJetPtBin.size()-1).second);
      fOneOverPtFit[iJetAxis][iCentrality]->SetParameters(1,1);
      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
        hDeltaJetAxis[iJetAxis][iCentrality][iJetPt] = NULL;
        for(int iJetRadius = 0; iJetRadius < knJetRadii; iJetRadius++){
          for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
            for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){ 
              hEnergyEnergyCorrelator[iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt] = NULL;
              hEnergyEnergyCorrelatorAxisRatio[iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt] = NULL;
              hEnergyEnergyCorrelatorRadiusRatio[iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt] = NULL;
            } // Track pT loop
          } // Energy weight loop
        } // Jet radius loop
      } // Jet pT loop
    } // Jet axis loop
  } // Centrality loop 
  
  // Helper variables
  std::pair<double, double> drawingRange = std::make_pair(0.008, 0.39);
  double epsilon = 0.0001;
  int iTrackPt, iTrackPtReference;
  int iJetPt, iJetPtReference;
  int iCentrality, iCentralityReference;

  // Load the energy-energy correlator histograms
  for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
    for(int iJetRadius = 0; iJetRadius < knJetRadii; iJetRadius++){
      for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
        for(auto centralityBin : comparedCentralityBin){
          for(auto jetPtBin : comparedJetPtBin){
            for(auto trackPtBin : comparedTrackPtBin){

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

              hEnergyEnergyCorrelator[iJetRadius][iEnergyWeight][iJetAxis][iCentralityReference][iJetPtReference][iTrackPtReference] = histograms[iJetRadius][iEnergyWeight][iJetAxis]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, pairingType, subevent);
            
            } // Track pT loop
          } // Jet pT loop
        } // Centrality loop
      } // Jet axis type loop
    } // Jet radius loop
  } // Energy weight loop

  // Load the DeltaR between E-scheme and WTA axes histograms
  for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
    for(auto jetPtBin : comparedJetPtBin){ 
      iJetPt = card[kR0p4][kNominalWeight][iJetAxis]->FindBinIndexJetPtEEC(jetPtBin);
      iJetPtReference = card[0][0][0]->FindBinIndexJetPtEEC(jetPtBin);
      for(auto centralityBin : comparedCentralityBin){

        if(centralityBin.first < 0){
          iCentrality = 0;
          iCentralityReference = 0;
        } else {
          iCentrality = card[kR0p4][kNominalWeight][iJetAxis]->FindBinIndexCentrality(centralityBin);
          iCentralityReference = card[0][0][0]->FindBinIndexCentrality(centralityBin);
        } 

        hDeltaJetAxis[iJetAxis][iCentralityReference][iJetPtReference] = histograms[kR0p4][kNominalWeight][iJetAxis]->GetHistogramJetDeltaAxis(iCentrality, iJetPt);
      } // Centrality loop
    } // Jet pT loop
  } // Jet axis loop

  // ====================================================
  //               Normalize distributions
  // ====================================================

  int lowNormalizationBin, highNormalizationBin;
  if(normalizeDistributions){
    for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
      for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
        for(int iJetRadius = 0; iJetRadius < knJetRadii; iJetRadius++){
          for(auto jetPtBin : comparedJetPtBin){
            iJetPt = card[0][0][0]->FindBinIndexJetPtEEC(jetPtBin);
            for(auto trackPtBin : comparedTrackPtBin){
              iTrackPt = card[0][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
              // Do the same for all the selected centrality bins for PbPb
              for(auto centralityBin : comparedCentralityBin){
                iCentrality = centralityBin.first < 0 ? 0 : card[0][0][0]->FindBinIndexJetPtEEC(centralityBin);

                // Normalize the distributions to one in the drawingRange
                lowNormalizationBin = hEnergyEnergyCorrelator[iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(drawingRange.first + epsilon);
                highNormalizationBin = hEnergyEnergyCorrelator[iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(drawingRange.second - epsilon);

                hEnergyEnergyCorrelator[iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->Scale(1 / hEnergyEnergyCorrelator[iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->Integral(lowNormalizationBin, highNormalizationBin, "width"));
            
              } // Centrality loop
            } // Track pT loop
          } // Jet pT loop
        } // Jet radius loop
      } // Jet axis type loop
    } // Energy weight loop
  } // Normalizing distributions

  // ====================================================
  //                 Calculating ratios
  // ====================================================

  // Ratio calculation for WTA/E-scheme
  for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
    for(int iJetRadius = 0; iJetRadius < knJetRadii; iJetRadius++){
      for(auto jetPtBin : comparedJetPtBin){
        iJetPt = card[0][0][0]->FindBinIndexJetPtEEC(jetPtBin);
        for(auto trackPtBin : comparedTrackPtBin){
          iTrackPt = card[0][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
          for(auto centralityBin: comparedCentralityBin){
            iCentrality = centralityBin.first < 0 ? 0 : card[0][0][0]->FindBinIndexJetPtEEC(centralityBin);
            for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){

              hEnergyEnergyCorrelatorAxisRatio[iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelator[iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->Clone(Form("eecAxisRatio%d%d%d%d%d%d", iJetRadius, iEnergyWeight, iJetAxis, iCentrality, iJetPt, iTrackPt));
              hEnergyEnergyCorrelatorAxisRatio[iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->Divide(hEnergyEnergyCorrelator[iJetRadius][iEnergyWeight][kEscheme][iCentrality][iJetPt][iTrackPt]);
            } // Jet axis loop
          } // Centrality loop
        } // Track pT loop
      } // Jet pT loop
    } // Jet radius loop
  } // Energy weight loop

  // Ratio calculation for R=0.4/R=0.8
  for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
    for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
      for(auto jetPtBin : comparedJetPtBin){
        iJetPt = card[0][0][0]->FindBinIndexJetPtEEC(jetPtBin);
        for(auto trackPtBin : comparedTrackPtBin){
          iTrackPt = card[0][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
          for(auto centralityBin: comparedCentralityBin){
            iCentrality = centralityBin.first < 0 ? 0 : card[0][0][0]->FindBinIndexJetPtEEC(centralityBin);
            for(int iJetRadius = 0; iJetRadius < knJetRadii; iJetRadius++){

              hEnergyEnergyCorrelatorRadiusRatio[iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelator[iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->Clone(Form("eecRadiusRatio%d%d%d%d%d%d", iJetRadius, iEnergyWeight, iJetAxis, iCentrality, iJetPt, iTrackPt));
              hEnergyEnergyCorrelatorRadiusRatio[iJetRadius][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->Divide(hEnergyEnergyCorrelator[kR0p8][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]);
            } // Jet axis loop
          } // Centrality loop
        } // Track pT loop
      } // Jet pT loop
    } // Jet radius loop
  } // Energy weight loop

  // ====================================================
  //              Creating derived histograms
  // ====================================================

  // Create histograms with average DeltaR between E-scheme and WTA axes as a function of jet pT
  const int maxJetPtBins = nJetPtBinsEEC;
  double jetPtBinBorderArray[maxJetPtBins+1];
  double meanDeltaJetAxis, meanDeltaJetAxisError;
  int iBin = 1;

  // Sanity check that we have allocated enough space for the double array
  if(nJetPtBins > maxJetPtBins){
    cout << "ERROR! You have more jet pT bins for average DeltaR between E-scheme and WTA axes than allowed by the code!" << endl;
    cout << "Please increase the value of maxJetPtBins variable in \"Creating derived histograms\"-section." << endl;
  }

  // Collect the jet pT borders into an array that is used to define the x-axis for the average DeltaR between E-scheme and WTA axes as a function of jet pT histogram
  iJetPt = 0;
  for(auto jetPtBin : comparedJetPtBin){
    jetPtBinBorderArray[iJetPt++] = jetPtBin.first;
  }
  jetPtBinBorderArray[iJetPt] = comparedJetPtBin.at(iJetPt-1).second;

  // Create the histograms
  for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
    for(auto centralityBin : comparedCentralityBin){
      iCentrality = centralityBin.first < 0 ? 0 : card[0][0][0]->FindBinIndexJetPtEEC(centralityBin);
      
      // Make a new histogram with the x-axis defined by the determined jet pT bins
      hAverageDeltaJetAxisVsPt[iJetAxis][iCentrality] = new TH1D(Form("averageDeltaJetAxis%d%d", iJetAxis, iCentrality), Form("averageDeltaJetAxis%d%d", iJetAxis, iCentrality), nJetPtBins, jetPtBinBorderArray);
      hAverageDeltaJetAxisVsPt[iJetAxis][iCentrality]->Sumw2();

      // Once the histogram is created, fill each bin with the mean value from the DeltaJetAxis histograms
      iBin = 1;
      for(auto jetPtBin : comparedJetPtBin){
        iJetPt = card[0][0][0]->FindBinIndexJetPtEEC(jetPtBin);

        meanDeltaJetAxis = hDeltaJetAxis[iJetAxis][iCentrality][iJetPt]->GetMean();
        meanDeltaJetAxisError = hDeltaJetAxis[iJetAxis][iCentrality][iJetPt]->GetMeanError();

        hAverageDeltaJetAxisVsPt[iJetAxis][iCentrality]->SetBinContent(iBin, meanDeltaJetAxis);
        hAverageDeltaJetAxisVsPt[iJetAxis][iCentrality]->SetBinError(iBin++, meanDeltaJetAxisError);
      }

      if(fitOneOverPt){
        hAverageDeltaJetAxisVsPt[iJetAxis][iCentrality]->Fit(fOneOverPtFit[iJetAxis][iCentrality], "N");
      }
    } // Centrality loop
  } // Jet axis loop
  
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

  // Automatic color selection for all different jet pT bins
  gStyle->SetPalette(kBird);
  auto niceColors = TColor::GetPalette();
  double step = 254.0/nJetPtBins;
  int jetPtColor[nJetPtBins];
  for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){
    jetPtColor[iJetPt] = niceColors.At(TMath::Ceil(iJetPt*step));
  }

  // Legend and line
  TLegend* legend;
  TLegend* jetPtLegend;
  TLine* oneLine = new TLine(drawingRange.first, 1, drawingRange.second, 1);
  oneLine->SetLineStyle(2);

  // Draw the average DeltaR between E-scheme and WTA axes as a function of jet pT
  if(drawDeltaJetAxisWithPt){
    for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
      for(auto centralityBin : comparedCentralityBin){
        iCentrality = centralityBin.first < 0 ? 0 : card[0][0][0]->FindBinIndexJetPtEEC(centralityBin);
        systemForLegend = card[kR0p4][kNominalWeight][iJetAxis]->GetAlternativeDataType(false);

        // Create the legend and add binning information to it
        legend = new TLegend(0.53,0.68,0.93,0.84);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);

        legend->AddEntry((TObject*) 0, systemForLegend.Data(), "");
        if(centralityBin.first >= 0) legend->AddEntry((TObject*) 0, Form("Centrality: %d-%d%%", centralityBin.first, centralityBin.second), "");
        //legend->AddEntry((TObject*) 0, Form("Jet selection: %s", jetAxisDescription[iJetAxis].Data()), "");

        // Setup centrality strings
        compactCentralityString = centralityBin.first == -1 ? "_pp" : Form("_C=%d-%d", centralityBin.first, centralityBin.second);

        // Set drawing style for the histograms
        hAverageDeltaJetAxisVsPt[iJetAxis][iCentrality]->SetMarkerStyle(markerStyle[0]);
        hAverageDeltaJetAxisVsPt[iJetAxis][iCentrality]->SetMarkerColor(color[0]);
        hAverageDeltaJetAxisVsPt[iJetAxis][iCentrality]->SetLineColor(color[0]); 

        // Set drawing ranges for the axes
        hAverageDeltaJetAxisVsPt[iJetAxis][iCentrality]->GetYaxis()->SetRangeUser(deltaJetAxisZoom.first, deltaJetAxisZoom.second);

        // Draw the histogram
        drawer->DrawHistogram(hAverageDeltaJetAxisVsPt[iJetAxis][iCentrality], "Jet p_{T} (GeV)", "#LT#phi#GT", " ", "p");

        if(fitOneOverPt){
          fOneOverPtFit[iJetAxis][iCentrality]->Draw("same");
          legend->AddEntry(fOneOverPtFit[iJetAxis][iCentrality], Form("%.2f + %.2f/p_{T}", fOneOverPtFit[iJetAxis][iCentrality]->GetParameter(0), fOneOverPtFit[iJetAxis][iCentrality]->GetParameter(1)), "l");
        }

        // Draw the legend
        legend->Draw();

        // Save the figures to a file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/averageDeltaJetAxis%s%s%s.%s", saveComment.Data(), jetAxisDescription[iJetAxis].Data(), compactCentralityString.Data(), figureFormat));
        }

      } // Centrality loop
    } // Jet axis type loop
  } // Drawing average DeltaR between E-scheme and WTA axes as a function of jet pT

  // Draw illustration plots for edge loss
  if(drawEdgeLossIllustration){

    for(int iEnergyWeight = 0; iEnergyWeight < knEnergyWeights; iEnergyWeight++){
      for(int iJetAxis = 0; iJetAxis < knJetAxes; iJetAxis++){
        for(auto trackPtBin : comparedTrackPtBin){
          iTrackPt = card[0][0][0]->GetBinIndexTrackPtEEC(trackPtBin);
          for(auto centralityBin : comparedCentralityBin){
            iCentrality = centralityBin.first < 0 ? 0 : card[0][0][0]->FindBinIndexJetPtEEC(centralityBin);
            systemForLegend = card[kR0p4][iEnergyWeight][iJetAxis]->GetAlternativeDataType(false);

            iJetPtReference = card[0][0][0]->FindBinIndexJetPtEEC(comparedJetPtBin.at(0));

            // Logarithmic DeltaR axis
            drawer->SetLogX(true);

            // Create the legend and add binning information to it
            legend = new TLegend(0.2,0.18,0.45,0.46);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);

            legend->AddEntry((TObject*) 0, systemForLegend.Data(), "");
            legend->AddEntry((TObject*) 0, Form("Jet axis: %s", jetAxisDescription[iJetAxis].Data()), "");
            if(centralityBin.first >= 0) legend->AddEntry((TObject*) 0, Form("Centrality: %d-%d%%", centralityBin.first, centralityBin.second), "");
            legend->AddEntry((TObject*) 0, Form("n = %d", iEnergyWeight+1), "");
            legend->AddEntry((TObject*) 0, Form("p_{T}^{ch} > %.1f GeV", trackPtBin), "");

            // Create a legend to hold jet pT bins
            jetPtLegend= new TLegend(0.2,0.6,0.92,0.9);
            jetPtLegend->SetFillStyle(0);jetPtLegend->SetBorderSize(0);
            jetPtLegend->SetTextSize(0.03);jetPtLegend->SetTextFont(62);
            jetPtLegend->SetNColumns(3);
            jetPtLegend->SetHeader("Jet p_{T}","C");

            // Setup centrality and track pT strings
            compactTrackPtString = Form("_T>%.1f",trackPtBin);
            compactTrackPtString.ReplaceAll(".","v");
            compactCentralityString = centralityBin.first == -1 ? "_pp" : Form("_C=%d-%d", centralityBin.first, centralityBin.second);


            // Set drawing style for all histograms
            iBin = 0;
            for(auto jetPtBin : comparedJetPtBin){
              iJetPt = card[0][0][0]->FindBinIndexJetPtEEC(jetPtBin);
              hEnergyEnergyCorrelatorRadiusRatio[kR0p4][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->SetLineColor(jetPtColor[iBin++]);
            } // Jet pT loop


            // Set the x- and y-axis drawing ranges
            hEnergyEnergyCorrelatorRadiusRatio[kR0p4][iEnergyWeight][iJetAxis][iCentrality][iJetPtReference][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
            hEnergyEnergyCorrelatorRadiusRatio[kR0p4][iEnergyWeight][iJetAxis][iCentrality][iJetPtReference][iTrackPt]->GetYaxis()->SetRangeUser(edgeLossIllustrationZoom.first, edgeLossIllustrationZoom.second);

            // Draw the histograms to the upper canvas
            drawer->DrawHistogram(hEnergyEnergyCorrelatorRadiusRatio[kR0p4][iEnergyWeight][iJetAxis][iCentrality][iJetPtReference][iTrackPt], "#Deltar", "EEC #frac{R = 0.4}{R = 0.8}", " ", "l");

            for(auto jetPtBin : comparedJetPtBin){
              iJetPt = card[0][0][0]->FindBinIndexJetPtEEC(jetPtBin);
              if(iJetPt == iJetPtReference) continue;
              hEnergyEnergyCorrelatorRadiusRatio[kR0p4][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt]->Draw("same,l");
            } // Jet pT loop

            // Add legends for drawn histograms TODO: These need to be in a separate small legend somewhere
            for(auto jetPtBin : comparedJetPtBin){
              iJetPt = card[0][0][0]->FindBinIndexJetPtEEC(jetPtBin);
              jetPtLegend->AddEntry(hEnergyEnergyCorrelatorRadiusRatio[kR0p4][iEnergyWeight][iJetAxis][iCentrality][iJetPt][iTrackPt], Form("%d-%d GeV", jetPtBin.first, jetPtBin.second), "l");
            } // Jet pT loop*/
  
            // Draw the legends
            legend->Draw();
            jetPtLegend->Draw();
          
            // Save the figures to a file
            if(saveFigures){
              gPad->GetCanvas()->SaveAs(Form("figures/edgeLossIllustration%s%s%s%s%s.%s", jetAxisDescription[iJetAxis].Data(), saveComment.Data(), compactEnergyWeightString[iEnergyWeight].Data(), compactCentralityString.Data(), compactTrackPtString.Data(), figureFormat));
            }
        } // Centrality looop
      } // Track pT loop
    } // Jet axis type loop
  } // Energy weight loop
  } // Drawing edge loss illustration
}
