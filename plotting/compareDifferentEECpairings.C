#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"
#include "AlgorithmLibrary.h"
#include "SystematicUncertaintyOrganizer.h"
#include "../src/EECHistograms.h"

/*
 * Macro for comparing final energy-energy correlators results
 */
void compareDifferentEECpairings(){
  
  // Predefined configuration for public figures
  bool truthToMixedConeBackgroundComparison = false;
  bool truthSignalToBackgroundSubtractedComparison = false;
  bool backgroundComponentsFromData = false;
  bool backgroundComponentsFromSimulation = false;

  // Sanity check for public figure drawing
  if(truthToMixedConeBackgroundComparison + truthSignalToBackgroundSubtractedComparison + backgroundComponentsFromData + backgroundComponentsFromSimulation > 1){
    cout << "ERROR! Only one predefined figure can be activated at a time." << endl;
    cout << "Please select only one figure to draw," << endl;
    return;
  } 

  // Studied file
  TString fileName = "veryNonShiftedData_processed.root";

  // eecAnalysis_akFlowJet_nominalEnergyWeight_combinedMixedConeBackground_processed_2024-04-25.root
  // eecAnalysis_akFlowJet_energyWeightSquared_combinedMixedConeBackground_processed_2024-05-02.root

  // eecAnalysis_akFlowJet_energyWeightSquared_allBackgrounds_includeSystematics_processed_2024-03-31.root
  // eecAnalysis_akFlowJet_nominalEnergyWeight_allBackgrounds_includeSystematics_processed_2024-03-31.root
  // eecAnalysis_akFlowJet_nominalEnergyWeight_allBackgrounds_noJetsInMidrapidity_lowStats_processed_2024-04-10.root
  // eecAnalysis_akFlowJet_nominalEnergyWeight_mixedCone_forwardRapidity_processed_2024-03-15.root
  // PbPbMC2018_GenGen_eecAnalysis_4pCentShift_cutBadPhi_nominalEnergyWeight_nominalSmear_mixedConeBackground_processed_2024-03-15.root
  // PbPbMC2018_GenGen_eecAnalysis_4pCentShift_cutBadPhi_nominalEnergyWeight_nominalSmear_onlyMixedConeBackground_processed_2024-03-15.root
  // PbPbMC2018_GenGen_eecAnalysis_4pCentShift_cutBadPhi_energyWeightSquared_allBackgrounds_noJetsInMidrapidity_mostStats_processed_2024-04-10.root
  // PbPbMC2018_GenGen_eecAnalysis_4pCentShift_cutBadPhi_nominalEnergyWeight_allBackgrounds_noJetsInMidrapidity_mostStats_processed_2024-04-10.root

  // data/PbPbMC2018_GenGen_eecAnalysis_4pCentShift_cutBadPhi_nominalEnergyWeight_allBackgrounds_someJobsMissing_processed_2024-04-01.root
  // data/PbPbMC2018_GenGen_eecAnalysis_4pCentShift_cutBadPhi_nominalEnergyWeight_allBackgrounds_matchMultiplicity_lowStats_processed_2024-04-12.root
  // data/PbPbMC2018_GenGen_eecAnalysis_4pCentShift_cutBadPhi_energyWeightSquared_allBackgrounds_matchMultiplicity_lowStats_processed_2024-04-12.root

  // Use predefined file for the comparison
  if(truthToMixedConeBackgroundComparison || truthSignalToBackgroundSubtractedComparison || backgroundComponentsFromSimulation){
    fileName = "data/PbPbMC2018_GenGen_eecAnalysis_4pCentShift_cutBadPhi_nominalEnergyWeight_allBackgrounds_matchMultiplicity_someMissing_processed_2024-04-24.root";
  }

  // Use predefined file when plotting 
  if(backgroundComponentsFromData){
    fileName = "data/eecAnalysis_akFlowJet_nominalEnergyWeight_combinedMixedConeBackground_unfoldingWithNominalSmear_processed_2024-05-28.root";
  }
  
  // Open the files and check that they exist
  TFile* inputFile = TFile::Open(fileName);

  if(inputFile == NULL){
    cout << "Error! The file " << fileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }

  EECCard* card = new EECCard(inputFile);

  bool isPbPb = card->GetDataType().Contains("PbPb");
  
  // ====================================================
  //               Binning configuration
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = card->GetNCentralityBins();
  const int nJetPtBinsEEC = card->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card->GetNTrackPtBinsEEC();

  // Define the bins that are compared
  std::vector<std::pair<double,double>> comparedCentralityBin;
  comparedCentralityBin.push_back(std::make_pair(0,10));
  //comparedCentralityBin.push_back(std::make_pair(10,30));
  //comparedCentralityBin.push_back(std::make_pair(30,50));
  //comparedCentralityBin.push_back(std::make_pair(50,90));
  
  std::vector<std::pair<double,double>> comparedJetPtBin;
  comparedJetPtBin.push_back(std::make_pair(30,40));
  //comparedJetPtBin.push_back(std::make_pair(140,160));
  //comparedJetPtBin.push_back(std::make_pair(160,180));
  //comparedJetPtBin.push_back(std::make_pair(180,200));

  std::vector<double> comparedTrackPtBin;
  comparedTrackPtBin.push_back(1.0);
  //comparedTrackPtBin.push_back(2.0);
  //comparedTrackPtBin.push_back(3.0);

  // For MC, just write 0-10 instead of 4-14 etc.
  int adjustCentralityPercentage = 1;

  TString systemString = card->GetAlternativeDataType(false);

  // Enumeration for different pairing types
  enum enumPairingType{kSameJetPair, kSignalReflectedConePair, kReflectedConePair, kSignalMixedConePair, kMixedConePair, kSignalSecondMixedConePair, kMixedMixedConePair, kSecondMixedConePair, kSignalPerpendicularConePair, kPerpendicularConePair, kSignalSecondPerpendicularConePair, kPerpendicularPerpendicularConePair, kSecondPerpendicularConePair, kCompiledBackgroundMixed, kCompiledBackgroundPerpendicularCone, kCompiledBackgroundTruthLevel, kCombinedSignalMixedConePair, kCombinedMixedConePair, kMixedBackgroundSubtracted, kPerpendicularBackgroundSubtracted, knPairingTypes};
  enum enumSubeventType{kPythiaPythia, kPythiaHydjet, kHydjetPythia, kHydjetHydjet, knSubeventCombinations, kAllBackground, knSubeventTypes};

  // Indices that go into pairing here:
  // First index = Jet cone pairing index
  //   EECHistograms::kSameJetPair
  //   EECHistograms::kSignalReflectedConePair
  //   EECHistograms::kReflectedConePair
  //   EECHistograms::kSignalMixedConePair
  //   EECHistograms::kReflectedMixedConePair
  //   EECHistograms::kMixedConePair
  //   EECHistograms::kSignalSecondMixedConePair
  //   EECHistograms::kReflectedSecondMixedConePair
  //   EECHistograms::kMixedMixedConePair
  //   EECHistograms::kSecondMixedConePair
  // Second index = Subevent pairing index
  //   EECHistograms::kPythiaPythia
  //   EECHistograms::kPythiaHydjet
  //   EECHistograms::kHydjetPythia
  //   EECHistograms::kHydjetHydjet
  //   EECHistograms::knSubeventCombinations (accept any subevent combination)
  std::vector<std::pair<int,int>> comparedEnergyEnergyCorrelatorPairings;
  comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kSameJetPair, kPythiaPythia));
  comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kMixedBackgroundSubtracted, knSubeventCombinations));
  comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kPerpendicularBackgroundSubtracted, knSubeventCombinations));
  
  //comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kSameJetPair, kPythiaHydjet));
  //comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kSameJetPair, kHydjetHydjet));
  //comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kMixedMixedConePair, knSubeventCombinations));
  //comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kCombinedMixedConePair, knSubeventCombinations));
  
  // Choose legen position for the plots
  enum enumLegendPosition{kBottomLeft, kBottomRight, kTopLeft, kPaperClosure};
  int legendPosition = kPaperClosure;

  // Option to manually provide legend text for the compared distributions
  bool useManualLegend = true;
  std::vector<TString> manualLegend;
  manualLegend.push_back("True signal");
  manualLegend.push_back("Mixed cone, plain HFSum");
  manualLegend.push_back("Perpendicular cone signal");
  //manualLegend.push_back("Pythia+EPOS");
  //manualLegend.push_back("EPOS+EPOS");

  // Flag for normalizing distributions to one over studied deltaR range
  bool normalizeDistributions = true;

  // Compare integrals of the distributions
  bool compareIntegrals = false;

  // It makes no sense to normalize distributions to 1 if you want to compare their integrals...
  if(compareIntegrals) normalizeDistributions = false;

  // If we are dealing with MC, shift the centrality by 4% as is done in order to match background energy density
  if(isPbPb){
    if(card->GetDataType().Contains("MC")){
      for(auto& centralityBin : comparedCentralityBin){
        centralityBin.first += 4;
        centralityBin.second += 4;
        systemString = "Pythia+Hydjet";
      }
    } else {
      adjustCentralityPercentage = 0; // Never adjust centrality values for data
    }
  } else {
    // If we are not doing PbPb collisions, ignore centrality
    comparedCentralityBin.clear();
    comparedCentralityBin.push_back(std::make_pair(-1,-1));
  }

  // Override the selected configuration if we are dealing with predefined figures
  if(truthToMixedConeBackgroundComparison){

    // Select the correct pairings for the background comparison
    comparedEnergyEnergyCorrelatorPairings.clear();
    comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kSameJetPair, kAllBackground));
    comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kCompiledBackgroundMixed, knSubeventCombinations));

    useManualLegend = true;
    manualLegend.clear();
    manualLegend.push_back("True background");
    manualLegend.push_back("Mixed cone background");

    // Do also other configuration
    legendPosition = kPaperClosure;
    normalizeDistributions = false;
    compareIntegrals = false;
 
  }

  if(truthSignalToBackgroundSubtractedComparison){

    // Select the correct pairings for the background comparison
    comparedEnergyEnergyCorrelatorPairings.clear();
    comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kSameJetPair, kPythiaPythia));
    comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kSameJetPair, knSubeventCombinations));
    comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kMixedBackgroundSubtracted, knSubeventCombinations));

    useManualLegend = true;
    manualLegend.clear();
    manualLegend.push_back("True signal");
    manualLegend.push_back("Signal + background");
    manualLegend.push_back("Extracted signal");

    // Do also other configuration
    legendPosition = kBottomLeft;
    normalizeDistributions = false;
    compareIntegrals = false;
 
  }

  if(backgroundComponentsFromData || backgroundComponentsFromSimulation){

    // Select the correct pairings for the background comparison
    comparedEnergyEnergyCorrelatorPairings.clear();
    comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kCompiledBackgroundMixed, knSubeventCombinations));
    comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kCombinedSignalMixedConePair, knSubeventCombinations));
    comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kCombinedMixedConePair, knSubeventCombinations));
    comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kMixedMixedConePair, knSubeventCombinations));

    useManualLegend = true;
    manualLegend.clear();
    manualLegend.push_back("Total background");
    manualLegend.push_back("Signal + mixed event");
    manualLegend.push_back("Mixed event 1 + Mixed event 1");
    manualLegend.push_back("Mixed event 1 + Mixed event 2");

    // Do also other configuration
    legendPosition = kBottomLeft;
    normalizeDistributions = false;
    compareIntegrals = false;
 
  }

  // ====================================================
  //                Drawing configuration
  // ====================================================
  
  // Figure saving
  const bool saveFigures = true;  // Save figures
  TString saveComment = "_plainHFMatch";   // Comment given for this specific file
  const char* figureFormat = "png"; // Format given for the figures

  int weightExponent = card->GetWeightExponent();
  if(weightExponent == 2){
    saveComment.Prepend("_energyWeightSquared");
  } else {
    saveComment.Prepend("_nominalEnergyWeight");
  }

  // Drawing configuration
  std::pair<double, double> ratioZoom = std::make_pair(0.85, 1.15);
  std::pair<double, double> eecZoom = std::make_pair(0.2, 30);
  bool logRatio = false;
  bool automaticZoom = true;

  // Configuration for predefined figures
  if(truthToMixedConeBackgroundComparison){
    ratioZoom = std::make_pair(0.91, 1.09);
    logRatio = false;
    automaticZoom = true;
  }

  if(truthSignalToBackgroundSubtractedComparison){
    //ratioZoom = std::make_pair(0, 14);
    ratioZoom = std::make_pair(0.8, 1.2);
    logRatio = false;
    automaticZoom = true;
  }

  // Configuration for predefined figures
  if(backgroundComponentsFromData || backgroundComponentsFromSimulation){
    ratioZoom = std::make_pair(0.01, 2);
    logRatio = true;
    automaticZoom = true;
  }

  // Sanity check for input. Ensure that all the selected bins actually exist in the input file.

  // Sanity check for centrality bins
  if(isPbPb){
    for(auto centralityBin : comparedCentralityBin){
      if(card->FindBinIndexCentrality(centralityBin) < 0){
        cout << "ERROR! Centrality bin " << centralityBin.first << "-" << centralityBin.second << " does not exist in file " << fileName.Data() << endl;
        cout << "Please only choose centrality bins that are included in the input files." << endl;
        return;
      }
    }
  }

  // Sanity check for jet pT bins
  for(auto jetPtBin : comparedJetPtBin){
    if(card->FindBinIndexJetPtEEC(jetPtBin) < 0){
      cout << "ERROR! Jet pT bin " << jetPtBin.first << "-" << jetPtBin.second << " does not exist in file " << fileName.Data() << endl;
      cout << "Please only choose jet pT bins that are included in the input files." << endl;
      return;
    }
  }

  // Sanity check for track pT bins
  for(auto trackPtBin : comparedTrackPtBin){
    if(card->GetBinIndexTrackPtEEC(trackPtBin) < 0){
      cout << "ERROR! Track pT cut > " << trackPtBin << " GeV does not exist in file " << fileName.Data() << endl;
      cout << "Please only choose track pT bins that are included in the input files." << endl;
      return;
    }
  } 
  
  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager* histograms = new EECHistogramManager(inputFile, card);

  // Load energy-energy correlators
  histograms->SetLoadEnergyEnergyCorrelators(true);

  // Energy-energy correlator histograms
  TH1D* hEnergyEnergyCorrelator[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][knPairingTypes][knSubeventTypes];
  TH1D* hEnergyEnergyCorrelatorRatio[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][knPairingTypes][knSubeventTypes];
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        for(int iPairing = 0; iPairing < knPairingTypes; iPairing++){
          for(int iSubevent = 0; iSubevent < knSubeventTypes; iSubevent++){
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][iPairing][iSubevent] = NULL;
            hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][iPairing][iSubevent] = NULL;
          } // Subevent type loop
        } // Pairing type loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
  
  // Helper variables
  std::pair<double, double> drawingRange = std::make_pair(0.008, 0.39);
  double epsilon = 0.0001;
  int lowNormalizationBin, highNormalizationBin;
  int iCentrality, iJetPt, iTrackPt;
  TH1D* helperHistogram;
  TH1D* backgroundHistogram;

  // Because we might define also derived histograms, we first need to load all possible pairing types and subevent combinations
  for(auto centralityBin : comparedCentralityBin){
    if(isPbPb){
      iCentrality = card->FindBinIndexCentrality(centralityBin);
    } else {
      iCentrality = 0;
    }
    for(auto jetPtBin : comparedJetPtBin){
      iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);
      for(auto trackPtBin : comparedTrackPtBin){
        iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);
        for(int iPairingType = 0; iPairingType < EECHistograms::knPairingTypes; iPairingType++){
          for(int iSubevent = 0; iSubevent <= EECHistograms::knSubeventCombinations; iSubevent++){

            // For data, there is no subevent division for the histograms
            if(!card->GetDataType().Contains("MC") && (iSubevent != EECHistograms::knSubeventCombinations)) continue;

            // Read the histogram corresponding to the defined bin
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent] = histograms->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iPairingType, iSubevent);

          } // Subevent loop
        } // Pairing type loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

  // If we have derived histograms, derive them before any normalizations take place
  for(auto pairingType : comparedEnergyEnergyCorrelatorPairings){
  
    // Compile total background contribution from different parts in truth level
    if(pairingType.first == kCompiledBackgroundTruthLevel){
      for(auto centralityBin : comparedCentralityBin){
        if(isPbPb){
          iCentrality = card->FindBinIndexCentrality(centralityBin);
        } else {
          iCentrality = 0;
        }
        for(auto jetPtBin : comparedJetPtBin){
          iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);
          for(auto trackPtBin : comparedTrackPtBin){
            iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);

            // Start by taking the signal-reflected cone pairs, which is the regular background estimate
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second] = (TH1D*) hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSignalReflectedConePair][knSubeventCombinations]->Clone(Form("compiledBackground%d%d%d%d%d", iCentrality, iJetPt, iTrackPt, pairingType.first, pairingType.second));

            // Subtract from the background fake-fake contribution between cones, since this is badly modeled
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSignalReflectedConePair][kHydjetHydjet], -1);

            // Add the real fake+fake contribution back to the background estimate
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSameJetPair][kHydjetHydjet]);

          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Compiling total background from different components

    // Compile total background contribution from different parts using jet cones in mixed events
    if(pairingType.first == kCompiledBackgroundMixed || pairingType.first == kCompiledBackgroundPerpendicularCone){

      int perpendicularShifter = 0;
      if(pairingType.first == kCompiledBackgroundPerpendicularCone){
        perpendicularShifter = 5;
      }

      for(auto centralityBin : comparedCentralityBin){
        if(isPbPb){
          iCentrality = card->FindBinIndexCentrality(centralityBin);
        } else {
          iCentrality = 0;
        }
        for(auto jetPtBin : comparedJetPtBin){
          iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);
          for(auto trackPtBin : comparedTrackPtBin){
            iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);

            // Start by taking the signal-mixed cone pairs, which is the basis of the background estimate
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second] = (TH1D*) hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSignalMixedConePair+perpendicularShifter][knSubeventCombinations]->Clone(Form("compiledBackground%d%d%d%d%d", iCentrality, iJetPt, iTrackPt, pairingType.first, pairingType.second));

            // Add the other signal-mixed cone component and normalize to keep the correct scaling
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSignalSecondMixedConePair+perpendicularShifter][knSubeventCombinations]);
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Scale(0.5);

            // Subtract from the background fake-fake contribution between cones, since this is badly modeled
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kMixedMixedConePair+perpendicularShifter][knSubeventCombinations], -1);

            // To improve statistics for the mixed cone only component, tae average of the two equivalent components
            helperHistogram = (TH1D*)hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kMixedConePair+perpendicularShifter][knSubeventCombinations]->Clone(Form("temporaryBackground%d%d%d%d%d", iCentrality, iJetPt, iTrackPt, pairingType.first, pairingType.second));
            helperHistogram->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSecondMixedConePair+perpendicularShifter][knSubeventCombinations]);
            helperHistogram->Scale(0.5);

            // Add the real fake+fake contribution back to the background estimate
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kMixedConePair+perpendicularShifter][knSubeventCombinations]);

          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Compiling total background from different components

    // Add together all background related subevents
    if(pairingType.second == kAllBackground){
      for(auto centralityBin : comparedCentralityBin){
        if(isPbPb){
          iCentrality = card->FindBinIndexCentrality(centralityBin);
        } else {
          iCentrality = 0;
        }
        for(auto jetPtBin : comparedJetPtBin){
          iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);
          for(auto trackPtBin : comparedTrackPtBin){
            iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);

            // Start from the signal+fake background
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second] = (TH1D*) hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][kPythiaHydjet]->Clone(Form("allBackground%d%d%d%d%d", iCentrality, iJetPt, iTrackPt, pairingType.first, pairingType.second));

            // If we are dealing with correlation between signal and refglected cones, add fake+signal component
            if(pairingType.first == kSignalReflectedConePair){
              hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][kHydjetPythia]);
            }

            // Next, add the fake+fake component to the total background
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][kHydjetHydjet]);

          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Adding together all background related subevents

    // Combine the two signal-mixed event components into one amalgamation
    if(pairingType.first == kCombinedSignalMixedConePair){
      for(auto centralityBin : comparedCentralityBin){
        if(isPbPb){
          iCentrality = card->FindBinIndexCentrality(centralityBin);
        } else {
          iCentrality = 0;
        }
        for(auto jetPtBin : comparedJetPtBin){
          iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);
          for(auto trackPtBin : comparedTrackPtBin){
            iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);

            // Start by taking the signal-mixed cone pairs
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second] = (TH1D*) hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSignalMixedConePair][knSubeventCombinations]->Clone(Form("combineSignalMixed%d%d%d%d%d", iCentrality, iJetPt, iTrackPt, pairingType.first, pairingType.second));

            // Add the signal-mixed cone pairs from the other mixed event
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSignalSecondMixedConePair][knSubeventCombinations]);

            // Take the average of these two to keep the scaling
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Scale(0.5);

          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Combining the two signal-mixed event components

    // Combine the two mixed event one components into one amalgamation
    if(pairingType.first == kCombinedMixedConePair){
      for(auto centralityBin : comparedCentralityBin){
        if(isPbPb){
          iCentrality = card->FindBinIndexCentrality(centralityBin);
        } else {
          iCentrality = 0;
        }
        for(auto jetPtBin : comparedJetPtBin){
          iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);
          for(auto trackPtBin : comparedTrackPtBin){
            iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);

            // Start by taking the signal-mixed cone pairs
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second] = (TH1D*) hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kMixedConePair][knSubeventCombinations]->Clone(Form("combineMixedOnly%d%d%d%d%d", iCentrality, iJetPt, iTrackPt, pairingType.first, pairingType.second));

            // Add the signal-mixed cone pairs from the other mixed event
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSecondMixedConePair][knSubeventCombinations]);

            // Take the average of these two to keep the scaling
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Scale(0.5);

          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Combining the two mixed event only components

    // Compile total background contribution from different parts using jet cones in mixed events. Then subtract it from the regular distribution
    if(pairingType.first == kMixedBackgroundSubtracted || pairingType.first == kPerpendicularBackgroundSubtracted){

      int perpendicularShifter = 0;
      if(pairingType.first == kPerpendicularBackgroundSubtracted){
        perpendicularShifter = 5;
      }

      for(auto centralityBin : comparedCentralityBin){
        if(isPbPb){
          iCentrality = card->FindBinIndexCentrality(centralityBin);
        } else {
          iCentrality = 0;
        }
        for(auto jetPtBin : comparedJetPtBin){
          iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);
          for(auto trackPtBin : comparedTrackPtBin){
            iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);

            // Start by taking the signal-mixed cone pairs
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second] = (TH1D*) hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSameJetPair][knSubeventCombinations]->Clone(Form("backgroundSubtracted%d%d%d%d%d", iCentrality, iJetPt, iTrackPt, pairingType.first, pairingType.second));

            // Start by taking the signal-mixed cone pairs, which is the basis of the background estimate
            backgroundHistogram = (TH1D*) hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSignalMixedConePair+perpendicularShifter][knSubeventCombinations]->Clone(Form("backgroundForSubtraction%d%d%d%d%d", iCentrality, iJetPt, iTrackPt, pairingType.first, pairingType.second));

            // Add the other signal-mixed cone component and normalize to keep the correct scaling
            backgroundHistogram->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSignalSecondMixedConePair+perpendicularShifter][knSubeventCombinations]);
            backgroundHistogram->Scale(0.5);

            // Subtract from the background fake-fake contribution between cones, since this is badly modeled
            backgroundHistogram->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kMixedMixedConePair+perpendicularShifter][knSubeventCombinations], -1);

            // To improve statistics for the mixed cone only component, tae average of the two equivalent components
            helperHistogram = (TH1D*)hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kMixedConePair+perpendicularShifter][knSubeventCombinations]->Clone(Form("temporaryBackgroundForSubtraction%d%d%d%d%d", iCentrality, iJetPt, iTrackPt, pairingType.first, pairingType.second));
            helperHistogram->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSecondMixedConePair+perpendicularShifter][knSubeventCombinations]);
            helperHistogram->Scale(0.5);

            // Add the real fake+fake contribution back to the background estimate
            backgroundHistogram->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kMixedConePair+perpendicularShifter][knSubeventCombinations]);

            // Finally subtract the background from the regular distribution
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(backgroundHistogram, -1);

          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Compiling total background from different components

  }

  // Go through all the histograms, normalize them, and calculate their ratios
  for(auto centralityBin : comparedCentralityBin){
    if(isPbPb){
      iCentrality = card->FindBinIndexCentrality(centralityBin);
    } else {
      iCentrality = 0;
    }
    for(auto jetPtBin : comparedJetPtBin){
      iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);
      for(auto trackPtBin : comparedTrackPtBin){
        iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);
        for(auto pairingType : comparedEnergyEnergyCorrelatorPairings){

          // Normalize the histogram to one in the studied range
          if(normalizeDistributions){
            lowNormalizationBin = hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetXaxis()->FindBin(drawingRange.first + epsilon);
            highNormalizationBin = hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetXaxis()->FindBin(drawingRange.second - epsilon);
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Scale(1 / hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Integral(lowNormalizationBin, highNormalizationBin, "width"));
          }

          // Calculate the ratio with respect to the first pairing type
          hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second] = (TH1D*) hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Clone(Form("pairingRatio%d%d%d%d%d", iCentrality, iJetPt, iTrackPt, pairingType.first, pairingType.second));
          hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Divide(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][comparedEnergyEnergyCorrelatorPairings.at(0).first][comparedEnergyEnergyCorrelatorPairings.at(0).second]);

        } // Pairing type loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
  
  // ==========================================================================
  //               Draw all the pairing types in the same figure
  // ==========================================================================
  
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
  TString energyWeightString = (weightExponent == 1) ? "n=1" : "n=2";
  TString legendString;
  int markerStyle[5] = {kFullCircle, kOpenSquare, kOpenCross, kFullStar, kFullCross};
  int color[] = {kBlack,kRed,kBlue,kGreen+3,kMagenta,kCyan,kOrange,kViolet+3,kPink-7,kSpring+3,kAzure-7};
  int styleIndex, legendIndex;
  bool firstLoop, secondLoop;

  double legendX1 = 0.18;
  double legendY1 = 0.04;
  double legendX2 = 0.45;
  double legendY2 = 0.58;

  if(legendPosition == kBottomRight){
    legendX1 = 0.58;
    legendY1 = 0.02;
    legendX2 = 0.85;
    legendY2 = 0.52;
  } else if (legendPosition == kTopLeft){
    legendX1 = 0.16;
    legendY1 = 0.43;
    legendX2 = 0.43;
    legendY2 = 0.95;
  } else if (legendPosition == kPaperClosure){
    legendX1 = 0.53;
    legendY1 = 0.03;
    legendX2 = 0.80;
    legendY2 = 0.46;
  }

  // Automatic zoom helper veriables
  double binContent;
  int firstBin, lastBin;
  TLegend* legend;
  TLegend* tagLegend;

  for(auto centralityBin : comparedCentralityBin){
    if(isPbPb){
      iCentrality = card->FindBinIndexCentrality(centralityBin);
      compactCentralityString = Form("_C=%.0f-%.0f", centralityBin.first, centralityBin.second);
    } else {
      iCentrality = 0;
      compactCentralityString = "";
    }
    for(auto jetPtBin : comparedJetPtBin){
      iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);
      compactJetPtString = Form("_J=%.0f-%.0f", jetPtBin.first, jetPtBin.second);
      for(auto trackPtBin : comparedTrackPtBin){
        iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);
        compactTrackPtString = Form("_T>%.1f",trackPtBin);
        compactTrackPtString.ReplaceAll(".","v");
          
        // Create a new canvas for the plot
        drawer->CreateSplitCanvas();

        // Logarithmic EEC axis
        drawer->SetLogY(true);

        // Create the legend and add jet and track pT information to it
        legend = new TLegend(legendX1, legendY1, legendX2, legendY2);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.06);legend->SetTextFont(62);

        if(isPbPb){
          legend->AddEntry((TObject*) 0, Form("%s %.0f-%.0f%%", systemString.Data(), centralityBin.first-4*adjustCentralityPercentage, centralityBin.second-4*adjustCentralityPercentage), "");
        }
        legend->AddEntry((TObject*) 0, Form("%.0f < jet p_{T} < %.0f GeV", jetPtBin.first, jetPtBin.second), "");
        legend->AddEntry((TObject*) 0, Form("p_{T}^{ch} > %.0f GeV, %s", trackPtBin, energyWeightString.Data()), "");

        // Set drawing style for all histograms
        styleIndex = 0;
        for(auto pairingType : comparedEnergyEnergyCorrelatorPairings){
          hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->SetMarkerStyle(markerStyle[styleIndex]);
          hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->SetMarkerColor(color[styleIndex]);
          hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->SetLineColor(color[styleIndex]); 
          hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->SetMarkerStyle(markerStyle[styleIndex]);
          hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->SetMarkerColor(color[styleIndex]);
          hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->SetLineColor(color[styleIndex++]);
        } // File loop

        // Automatic zooming for the drawn histograms
        if(automaticZoom){
          eecZoom.first = 1e12;
          eecZoom.second = 0;
          for(auto pairingType : comparedEnergyEnergyCorrelatorPairings){
            firstBin = hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetXaxis()->FindBin(drawingRange.first);
            lastBin = hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetXaxis()->FindBin(drawingRange.second);
            for(int iBin = firstBin; iBin <= lastBin; iBin++){
              binContent = hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetBinContent(iBin);
              if(binContent < eecZoom.first){
                if(binContent > 0){
                  eecZoom.first = binContent;
                }
              } 
              if(binContent > eecZoom.second) eecZoom.second = binContent;
            } // Bin loop
          } // File loop
          eecZoom.first = eecZoom.first / 2.0;
          eecZoom.second = eecZoom.second * 2.0;
        }

        firstLoop = true;
        legendIndex = 0;
        for(auto pairingType : comparedEnergyEnergyCorrelatorPairings){
          hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
          hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetYaxis()->SetRangeUser(eecZoom.first, eecZoom.second);

          if(firstLoop){
            drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second], "#Deltar", "EEC (A.U.)", " ");
            firstLoop = false;
          } else {
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Draw("same");
          }

          if(useManualLegend){
            legendString = manualLegend.at(legendIndex++);
          } else {
            legendString = histograms->GetPairingTypeSaveName(pairingType.first);
          }

          legend->AddEntry(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second], legendString.Data(), "p");
        }

        // Draw the legends to the upper pad
        legend->Draw();

        // Add CMS simulation tag to the canvas
        //if(legendPosition == kPaperClosure){
        //  tagLegend = new TLegend(0.16, 0.74, 0.43, 0.78);
        //  tagLegend->SetFillStyle(0);tagLegend->SetBorderSize(0);tagLegend->SetTextSize(0.07);tagLegend->SetTextFont(62);
        //  tagLegend->AddEntry((TObject*) 0, "CMS simulation", "");
        //  tagLegend->Draw();
        //}
          
        // Set linear scale for ratio, unless specifically asked to be logarithmic
        if(!logRatio) drawer->SetLogY(false);

        // Draw the histograms
        drawer->SetGridY(true);

        firstLoop = true;
        secondLoop = false;
        for(auto pairingType : comparedEnergyEnergyCorrelatorPairings){
          hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second); 
          hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
          if(firstLoop){
            firstLoop = false;
            secondLoop = true;
          } else if(secondLoop) {
            
            if(useManualLegend){
              if(truthToMixedConeBackgroundComparison){
                legendString = Form("#frac{Mixed cone bkg}{%s}", manualLegend.at(0).Data());
              } else {
                legendString = Form("#frac{Color}{%s}", manualLegend.at(0).Data());
              }
              
            } else {
              legendString = Form("#frac{Color}{%s}", histograms->GetPairingTypeSaveName(comparedEnergyEnergyCorrelatorPairings.at(0).first));
            }

            drawer->DrawHistogramToLowerPad(hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second], "#Deltar", legendString.Data(), " ");
            secondLoop = false;
          } else {
            hEnergyEnergyCorrelatorRatio[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Draw("same");
          }
        }

        drawer->SetGridY(false);
          
        // Save the figures to a file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/eecPairingTypeComparison%s%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
        }

        if(compareIntegrals){
          if(isPbPb){
            cout << "Bin: " << Form("Cent: %.0f-%.0f%%", centralityBin.first, centralityBin.second) << " " << Form("%.0f < jet p_{T} < %.0f", jetPtBin.first, jetPtBin.second) << " " << Form("%.1f < track p_{T}", trackPtBin) << endl;
          } else {
            cout << "Bin: " << Form("%.0f < jet p_{T} < %.0f", jetPtBin.first, jetPtBin.second) << " " << Form("%.1f < track p_{T}", trackPtBin) << endl;
          }
          legendIndex = 0;
          for(auto pairingType : comparedEnergyEnergyCorrelatorPairings){
            lowNormalizationBin = hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetXaxis()->FindBin(drawingRange.first + epsilon);
            highNormalizationBin = hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetXaxis()->FindBin(drawingRange.second - epsilon);
            cout << "Integral of " << manualLegend.at(legendIndex++) << " is " << hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Integral(lowNormalizationBin, highNormalizationBin, "width") << endl;
            cout << "Ratio to " << manualLegend.at(0) << " is " << hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Integral(lowNormalizationBin, highNormalizationBin, "width") / hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][comparedEnergyEnergyCorrelatorPairings.at(0).first][comparedEnergyEnergyCorrelatorPairings.at(0).second]->Integral(lowNormalizationBin, highNormalizationBin, "width") << endl;
          } // Loop over compared bins
        } // Compare integrals
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

}
