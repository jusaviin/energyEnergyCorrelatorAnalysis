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
  
  // Studied file
  TString fileName = "data/PbPbMC2018_GenGen_eecAnalysis_4pCentShift_cutBadPhi_energyWeightSquared_nominalSmear_onlyMixedConeBackground_processed_2024-03-20.root";

  // eecAnalysis_akFlowJet_nominalEnergyWeight_mixedConeBackground_processed_2024-03-13.root
  // eecAnalysis_akFlowJet_nominalEnergyWeight_mixedCone_midRapidity_processed_2024-03-15.root
  // eecAnalysis_akFlowJet_nominalEnergyWeight_mixedCone_forwardRapidity_processed_2024-03-15.root
  // PbPbMC2018_GenGen_eecAnalysis_4pCentShift_cutBadPhi_nominalEnergyWeight_nominalSmear_mixedConeBackground_processed_2024-03-15.root
  // PbPbMC2018_GenGen_eecAnalysis_4pCentShift_cutBadPhi_nominalEnergyWeight_nominalSmear_onlyMixedConeBackground_processed_2024-03-15.root
  // PbPbMC2018_GenGen_eecAnalysis_4pCentShift_cutBadPhi_energyWeightSquared_nominalSmear_onlyMixedConeBackground_processed_2024-03-20.root
  
  // Open the files and check that they exist
  TFile* inputFile = TFile::Open(fileName);

  if(inputFile == NULL){
    cout << "Error! The file " << fileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }

  EECCard* card = new EECCard(inputFile);
  
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
  comparedCentralityBin.push_back(std::make_pair(10,30));
  comparedCentralityBin.push_back(std::make_pair(30,50));
  comparedCentralityBin.push_back(std::make_pair(50,90));
  
  std::vector<std::pair<double,double>> comparedJetPtBin;
  comparedJetPtBin.push_back(std::make_pair(120,140));
  comparedJetPtBin.push_back(std::make_pair(140,160));
  comparedJetPtBin.push_back(std::make_pair(160,180));
  comparedJetPtBin.push_back(std::make_pair(180,200));

  std::vector<double> comparedTrackPtBin;
  comparedTrackPtBin.push_back(1.0);
  //comparedTrackPtBin.push_back(2.0);
  //comparedTrackPtBin.push_back(3.0);

  // Enumeration for different pairing types
  enum enumPairingType{kSameJetPair, kSignalReflectedConePair, kReflectedConePair, kSignalMixedConePair, kReflectedMixedConePair, kMixedConePair, kCompiledBackground, kCompiledBackgroundTruthLevel, knPairingTypes};
  enum enumSubeventType{kPythiaPythia, kPythiaHydjet, kHydjetPythia, kHydjetHydjet, knSubeventCombinations, kAllBackground, knSubeventTypes};

  // Indices that go into pairing here:
  // First index = Jet cone pairing index
  //   EECHistograms::kSameJetPair
  //   EECHistograms::kSignalReflectedConePair
  //   EECHistograms::kReflectedConePair
  //   EECHistograms::kSignalMixedConePair
  //   EECHistograms::kReflectedMixedConePair
  //   EECHistograms::kMixedConePair
  // Second index = Subevent pairing index
  //   EECHistograms::kPythiaPythia
  //   EECHistograms::kPythiaHydjet
  //   EECHistograms::kHydjetPythia
  //   EECHistograms::kHydjetHydjet
  //   EECHistograms::knSubeventCombinations (accept any subevent combination)
  std::vector<std::pair<int,int>> comparedEnergyEnergyCorrelatorPairings;
  comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kSameJetPair, kAllBackground));
  comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kCompiledBackground, knSubeventCombinations));
  //comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kReflectedConePair, kHydjetHydjet));
  //comparedEnergyEnergyCorrelatorPairings.push_back(std::make_pair(kMixedConePair, kHydjetHydjet));

  // Option to manually provide legend text for the compared distributions
  const bool useManualLegend = true;
  std::vector<TString> manualLegend;
  manualLegend.push_back("True background");
  manualLegend.push_back("Compiled background");
  //manualLegend.push_back("Mixed cone 2");

  // Flag for normalizing distributions to one over studied deltaR range
  bool normalizeDistributions = false;

  // Compare integrals of the distributions
  const bool compareIntegrals = true;

  // It makes no sense to normalize distributions to 1 if you want to compare their integrals...
  if(compareIntegrals) normalizeDistributions = false;

  // If we are dealing with MC, shift the centrality by 4% as is done in order to match background energy density
  if(card->GetDataType().Contains("MC")){
    for(auto& centralityBin : comparedCentralityBin){
      centralityBin.first += 4;
      centralityBin.second += 4;
    }
  }

  // ====================================================
  //                Drawing configuration
  // ====================================================
  
  // Figure saving
  const bool saveFigures = false;  // Save figures
  TString saveComment = "";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures

  int weightExponent = card->GetWeightExponent();
  if(weightExponent == 2){
    saveComment.Prepend("_energyWeightSquared");
  } else {
    saveComment.Prepend("_nominalEnergyWeight");
  }

  // Drawing configuration
  std::pair<double, double> ratioZoom = std::make_pair(0.8, 1.2);
  std::pair<double, double> eecZoom = std::make_pair(0.2, 30);
  const bool automaticZoom = true;

  // Sanity check for input. Ensure that all the selected bins actually exist in the input file.

  // Sanity check for centrality bins
  for(auto centralityBin : comparedCentralityBin){
    if(card->FindBinIndexCentrality(centralityBin) < 0){
      cout << "ERROR! Centrality bin " << centralityBin.first << "-" << centralityBin.second << " does not exist in file " << fileName.Data() << endl;
      cout << "Please only choose centrality bins that are included in the input files." << endl;
      return;
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

  // Choose the bin ranges
  histograms->SetCentralityBinRange(0, card->GetNCentralityBins() - 1);
  histograms->SetJetPtBinRangeEEC(0, card->GetNJetPtBinsEEC() - 1);
  histograms->SetTrackPtBinRangeEEC(0, card->GetNTrackPtBinsEEC() - 1);

  // Load the histograms from the file
  histograms->LoadProcessedHistograms();

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

  // Because we might define also derived histograms, we first need to load all possible pairing types and subevent combinations
  for(auto centralityBin : comparedCentralityBin){
    iCentrality = card->FindBinIndexCentrality(centralityBin);
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
        iCentrality = card->FindBinIndexCentrality(centralityBin);
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
    if(pairingType.first == kCompiledBackground){
      for(auto centralityBin : comparedCentralityBin){
        iCentrality = card->FindBinIndexCentrality(centralityBin);
        for(auto jetPtBin : comparedJetPtBin){
          iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);
          for(auto trackPtBin : comparedTrackPtBin){
            iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);

            // Start by taking the signal-reflected cone pairs, which is the regular background estimate
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second] = (TH1D*) hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kSignalReflectedConePair][knSubeventCombinations]->Clone(Form("compiledBackground%d%d%d%d%d", iCentrality, iJetPt, iTrackPt, pairingType.first, pairingType.second));

            // Subtract from the background fake-fake contribution between cones, since this is badly modeled
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kReflectedMixedConePair][knSubeventCombinations], -1);

            // Add the real fake+fake contribution back to the background estimate
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->Add(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][kReflectedConePair][knSubeventCombinations]);

          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Compiling total background from different components

    // Add together all background related subevents
    if(pairingType.second == kAllBackground){
      for(auto centralityBin : comparedCentralityBin){
        iCentrality = card->FindBinIndexCentrality(centralityBin);
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

  }

  // Go through all the histograms, normalize them, and calculate their ratios
  for(auto centralityBin : comparedCentralityBin){
    iCentrality = card->FindBinIndexCentrality(centralityBin);
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
  TString energyWeightString = (weightExponent == 1) ? "Nominal energy weight" : "Energy weight squared";
  TString legendString;
  int markerStyle[5] = {kFullCircle, kOpenSquare, kOpenCross, kFullStar, kFullCross};
  int color[] = {kBlack,kRed,kBlue,kGreen+3,kMagenta,kCyan,kOrange,kViolet+3,kPink-7,kSpring+3,kAzure-7};
  int styleIndex, legendIndex;
  bool firstLoop, secondLoop;

  // Automatic zoom helper veriables
  double minimumCandidate, maximumCandidate;

  for(auto centralityBin : comparedCentralityBin){
    iCentrality = card->FindBinIndexCentrality(centralityBin);
    compactCentralityString = Form("_C=%.0f-%.0f", centralityBin.first, centralityBin.second);
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
        TLegend* legend = new TLegend(0.18,0.04,0.45,0.58);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);

        legend->AddEntry((TObject*) 0, Form("Cent: %.0f-%.0f%%", centralityBin.first, centralityBin.second), "");
        legend->AddEntry((TObject*) 0, Form("%.0f < jet p_{T} < %.0f", jetPtBin.first, jetPtBin.second), "");
        legend->AddEntry((TObject*) 0, Form("%.1f < track p_{T}", trackPtBin), "");
        legend->AddEntry((TObject*) 0, energyWeightString, "");

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
          eecZoom.first = 10000;
          eecZoom.second = 0;
          for(auto pairingType : comparedEnergyEnergyCorrelatorPairings){
            hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
            minimumCandidate = hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetMinimum();
            maximumCandidate = hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second]->GetMaximum();
            if(minimumCandidate < eecZoom.first) eecZoom.first = minimumCandidate;
            if(maximumCandidate > eecZoom.second) eecZoom.second = maximumCandidate;
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
            drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt][pairingType.first][pairingType.second], "#Deltar", "EEC", " ");
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
          
        // Linear scale for the ratio
        drawer->SetLogY(false);

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
              legendString = Form("#frac{Color}{%s}", manualLegend.at(0).Data());
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
          cout << "Bin: " << Form("Cent: %.0f-%.0f%%", centralityBin.first, centralityBin.second) << " " << Form("%.0f < jet p_{T} < %.0f", jetPtBin.first, jetPtBin.second) << " " << Form("%.1f < track p_{T}", trackPtBin) << endl;
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
