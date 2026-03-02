#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Compare perpendicular cone background to mixed event background to understand the differences
 */
void perpendicularMixedBackground(){
  
  enum enumBackgroundType{kMixedCone, kPerpendicularCone, knBackgroundTypes};

  TString backgroundName[knBackgroundTypes];
  backgroundName[kMixedCone] = "Mixed cone";
  backgroundName[kPerpendicularCone] = "Perpendicular cone";
  
  // Input files: index 0 = Pythia+Hydjet simulation, index 1 = minimum bias Hydjet simulation
  TString inputFileName[knBackgroundTypes];
  //inputFileName[kMixedCone] = "data/PbPbMC2018_GenGen_eecAnalysis_4pCentShift_cutBadPhi_nominalEnergyWeight_mixedConeBackground_matchMultiplicity_semiOkStats_processed_2024-04-12.root";
  //inputFileName[kPerpendicularCone] = "data/pPb/PbPbMC2018_GenGen_akFlowJets_wtaAxis_4pCentShift_cutBadPhi_nominalEnergyWeight_perpendicularConeBackground_processed_2025-05-14.root";
  inputFileName[kMixedCone] = "data/eecAnalysis_akFlowJet_nominalEnergyWeight_combinedMixedConeBackground_unfoldingWithNominalSmear_processed_2024-05-28.root";
  inputFileName[kPerpendicularCone] = "data/perpendicularCone/eecAnalysis_akFlowJet_wtaAxis_nominalEnergyWeight_perpendicularCone_unfoldedResults_processed_2025-05-20.root";
  
  
  // Open the input files and read the analysis cards
  TFile* inputFile[knBackgroundTypes];
  EECCard* card[knBackgroundTypes];
  for(int iDataType = 0; iDataType < knBackgroundTypes; iDataType++){
    
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
  
  // ====================================================
  //    Binning configuration for the background study
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = card[0]->GetNCentralityBins();
  const int nJetPtBinsEEC = card[0]->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card[0]->GetNTrackPtBinsEEC();
  
  // Select explicitly which bins from the files are compared:
  std::vector<std::pair<double,double>> comparedCentralityBin;
  comparedCentralityBin.push_back(std::make_pair(0,10));
  //comparedCentralityBin.push_back(std::make_pair(10,30));
  //comparedCentralityBin.push_back(std::make_pair(30,50));
  //comparedCentralityBin.push_back(std::make_pair(50,90));

  std::vector<std::pair<double,double>> comparedJetPtBin;
  comparedJetPtBin.push_back(std::make_pair(120,140));
  //comparedJetPtBin.push_back(std::make_pair(40,50));
  //comparedJetPtBin.push_back(std::make_pair(50,60));
  //comparedJetPtBin.push_back(std::make_pair(60,80));

  std::vector<double> comparedTrackPtBin;
  comparedTrackPtBin.push_back(1.0);
  //comparedTrackPtBin.push_back(1.5);
  //comparedTrackPtBin.push_back(2.0);
  //comparedTrackPtBin.push_back(2.5);
  //comparedTrackPtBin.push_back(3.0);

  // For MC, shift the centrality bins by 4%
  bool isData = true;
  if(card[0]->GetDataType().Contains("MC")){
    for(auto& centralityBin : comparedCentralityBin){
      centralityBin.first += 4;
      centralityBin.second += 4;
    }
    isData = false;
  }
  
  // ====================================================
  //            Drawing style for the plots
  // ====================================================
  
  // Logarithmic axes
  const bool logDeltaR = true;
  const bool logEEC = true;
  
  // Figure saving
  const bool saveFigures = false;  // Save figures
  const char* saveComment = "_removeJetsFromReflectedCone";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures
  
  // Create and setup a new histogram manager to project and handle the histograms
  EECHistogramManager* histograms[knBackgroundTypes];
  for(int iDataType = 0; iDataType < knBackgroundTypes; iDataType++){
    histograms[iDataType] = new EECHistogramManager(inputFile[iDataType],card[iDataType]);
  }
  
  // Indices for different types of energy-energy correlators
  enum enumBackgroundComponents{kJetMixedCone, kSingleMixedCone, kMixedMixedCone, knBackgroundComponents};
  enum enumCorrelatorPairingTypes{kTotalEEC, kSignalEEC, kSignalFakeEEC, kFakeFakeEEC, kBackgroundEEC, knPairingTypesEEC};

  int backgroundSubeventCombination[knPairingTypesEEC];
  backgroundSubeventCombination[kTotalEEC] = EECHistograms::knSubeventCombinations;
  backgroundSubeventCombination[kSignalEEC] = EECHistograms::kPythiaPythia;
  backgroundSubeventCombination[kSignalFakeEEC] = EECHistograms::kPythiaHydjet;
  backgroundSubeventCombination[kFakeFakeEEC] = EECHistograms::kHydjetHydjet;
  backgroundSubeventCombination[kBackgroundEEC] = EECHistograms::knSubeventCombinations;

  int backgroundComponentMap[knBackgroundComponents][knBackgroundTypes];
  backgroundComponentMap[kJetMixedCone][kMixedCone] = EECHistograms::kSignalMixedConePair;
  backgroundComponentMap[kSingleMixedCone][kMixedCone] = EECHistograms::kMixedConePair;
  backgroundComponentMap[kMixedMixedCone][kMixedCone] = EECHistograms::kMixedMixedConePair;
  backgroundComponentMap[kJetMixedCone][kPerpendicularCone] = EECHistograms::kSignalPerpendicularConePair;
  backgroundComponentMap[kSingleMixedCone][kPerpendicularCone] = EECHistograms::kPerpendicularConePair;
  backgroundComponentMap[kMixedMixedCone][kPerpendicularCone] = EECHistograms::kPerpendicularPerpendicularConePair;

  TString backgroundComponentName[knBackgroundComponents];
  backgroundComponentName[kJetMixedCone] = "Jet-mixed cone pairs";
  backgroundComponentName[kSingleMixedCone] = "Mixed cone pairs";
  backgroundComponentName[kMixedMixedCone] = "Mixed-mixed cone pairs";

  TString backgroundComponentSaveName[knBackgroundComponents];
  backgroundComponentSaveName[kJetMixedCone] = "jetMixedConePair";
  backgroundComponentSaveName[kSingleMixedCone] = "mixedConePair";
  backgroundComponentSaveName[kMixedMixedCone] = "mixedMixedConePair";

  // ====================================================
  //                  Define histograms
  // ====================================================

  // Energy-energy correlator histograms separated by subevents from the Pythia+Hydjet simulation
  TH1D* hEnergyEnergyCorrelator[knBackgroundTypes][knPairingTypesEEC][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  
  // Background energy-energy correlators
  TH1D* hBackground[knBackgroundTypes][knPairingTypesEEC][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hBackgroundComponents[knBackgroundTypes][EECHistograms::knSubeventCombinations][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  
  // Ratios between different background methods
  TH1D* hBackgroundRatio[knBackgroundTypes][knPairingTypesEEC][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hBackgroundComponentRatio[knBackgroundTypes][EECHistograms::knSubeventCombinations][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  
  // Initialize the energy-energy correlator histogram array to NULL
  for(int iBackground = 0; iBackground < knBackgroundTypes; iBackground++){
    for(int iPairingType = 0; iPairingType < knPairingTypesEEC; iPairingType++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
          for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
            hEnergyEnergyCorrelator[iBackground][iPairingType][iCentrality][iJetPt][iTrackPt] = NULL;
            hBackground[iBackground][iPairingType][iCentrality][iJetPt][iTrackPt] = NULL;
            hBackgroundComponents[iBackground][iPairingType][iCentrality][iJetPt][iTrackPt] = NULL;
            hBackgroundRatio[iBackground][iPairingType][iCentrality][iJetPt][iTrackPt] = NULL;
            hBackgroundComponentRatio[iBackground][iPairingType][iCentrality][iJetPt][iTrackPt] = NULL;
          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Pairing type loop
  } // Background type loop

  // =========================================================
  //       Read histograms from files and normalize them
  // =========================================================

  // Helper variables
  std::pair<double, double> drawingRange = std::make_pair(0.008, 0.39);
  double epsilon = 0.0001;
  int lowNormalizationBin, highNormalizationBin;
  int iCentrality, iCentralityReference;
  int iTrackPt, iTrackPtReference;
  int iJetPt, iJetPtReference;
  double normalizationFactor;
  
  // Get the histograms from the histogram manager
  for(int iBackground = 0; iBackground < knBackgroundTypes; iBackground++){
    for(auto jetPtBin : comparedJetPtBin){
      for(auto trackPtBin : comparedTrackPtBin){
        for(auto centralityBin : comparedCentralityBin){

          // Find the proper binning and express it in term of the bins in the first file
          iJetPt = card[iBackground]->FindBinIndexJetPtEEC(jetPtBin);
          iTrackPt = card[iBackground]->GetBinIndexTrackPtEEC(trackPtBin);
          iCentrality = card[iBackground]->FindBinIndexCentrality(centralityBin);
          iJetPtReference = card[0]->FindBinIndexJetPtEEC(jetPtBin);
          iTrackPtReference = card[0]->GetBinIndexTrackPtEEC(trackPtBin);
          iCentralityReference = card[0]->FindBinIndexCentrality(centralityBin);
          
          // ================================================= //
          //             Read the jet cone histograms          //
          // ================================================= //
          
          // Histogram with all pair combinations. Normalize it to one
          hEnergyEnergyCorrelator[iBackground][kTotalEEC][iCentralityReference][iJetPtReference][iTrackPtReference] = histograms[iBackground]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJetPair, EECHistograms::knSubeventCombinations);
          lowNormalizationBin = hEnergyEnergyCorrelator[iBackground][kTotalEEC][iCentralityReference][iJetPtReference][iTrackPtReference]->GetXaxis()->FindBin(drawingRange.first);
          highNormalizationBin = hEnergyEnergyCorrelator[iBackground][kTotalEEC][iCentralityReference][iJetPtReference][iTrackPtReference]->GetXaxis()->FindBin(drawingRange.second);
          normalizationFactor = hEnergyEnergyCorrelator[iBackground][kTotalEEC][iCentralityReference][iJetPtReference][iTrackPtReference]->Integral(lowNormalizationBin, highNormalizationBin, "width");
          hEnergyEnergyCorrelator[iBackground][kTotalEEC][iCentralityReference][iJetPtReference][iTrackPtReference]->Scale(1/normalizationFactor);

          // Subevent decomposition can only be done for MC
          if(!isData){
          
            // Histogram with only signal. Normalize to the total integral
            hEnergyEnergyCorrelator[iBackground][kSignalEEC][iCentralityReference][iJetPtReference][iTrackPtReference] = histograms[iBackground]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJetPair, EECHistograms::kPythiaPythia);
            hEnergyEnergyCorrelator[iBackground][kSignalEEC][iCentralityReference][iJetPtReference][iTrackPtReference]->Scale(1/normalizationFactor);
          
            // Histograms with background contributions. Normalize to the total integral
            hEnergyEnergyCorrelator[iBackground][kSignalFakeEEC][iCentralityReference][iJetPtReference][iTrackPtReference] = histograms[iBackground]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJetPair, EECHistograms::kPythiaHydjet);
            hEnergyEnergyCorrelator[iBackground][kBackgroundEEC][iCentralityReference][iJetPtReference][iTrackPtReference] = (TH1D *) hEnergyEnergyCorrelator[iBackground][kSignalFakeEEC][iCentralityReference][iJetPtReference][iTrackPtReference]->Clone(Form("energyEnergyCorrelatorBackground%d%d%d%d", iBackground, iCentrality, iJetPt, iTrackPt));
            hEnergyEnergyCorrelator[iBackground][kFakeFakeEEC][iCentralityReference][iJetPtReference][iTrackPtReference] = histograms[iBackground]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJetPair, EECHistograms::kHydjetHydjet);
            hEnergyEnergyCorrelator[iBackground][kBackgroundEEC][iCentralityReference][iJetPtReference][iTrackPtReference]->Add(hEnergyEnergyCorrelator[iBackground][kFakeFakeEEC][iCentralityReference][iJetPtReference][iTrackPtReference]);

            hEnergyEnergyCorrelator[iBackground][kSignalFakeEEC][iCentralityReference][iJetPtReference][iTrackPtReference]->Scale(1/normalizationFactor);
            hEnergyEnergyCorrelator[iBackground][kFakeFakeEEC][iCentralityReference][iJetPtReference][iTrackPtReference]->Scale(1/normalizationFactor);
            hEnergyEnergyCorrelator[iBackground][kBackgroundEEC][iCentralityReference][iJetPtReference][iTrackPtReference]->Scale(1/normalizationFactor);
          }

          // =============================================================== //
          //       Properly normalize different background components        //
          // =============================================================== //    
          
          for(int iBackgroundComponent = 0; iBackgroundComponent < knBackgroundComponents; iBackgroundComponent++){
            
            // Start the background estimation from jet cone + mixed/perpendicular cone pairings
            hBackgroundComponents[iBackground][iBackgroundComponent][iCentralityReference][iJetPtReference][iTrackPtReference] = histograms[iBackground]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, backgroundComponentMap[iBackgroundComponent][iBackground], EECHistograms::knSubeventCombinations);

            hBackgroundComponents[iBackground][iBackgroundComponent][iCentralityReference][iJetPtReference][iTrackPtReference]->Scale(1/normalizationFactor);

            // Calculate the ratio with respect to mixed cone background components
            hBackgroundComponentRatio[iBackground][iBackgroundComponent][iCentralityReference][iJetPtReference][iTrackPtReference] = (TH1D*) hBackgroundComponents[iBackground][iBackgroundComponent][iCentralityReference][iJetPtReference][iTrackPtReference]->Clone(Form("backgroundComponentRatio%d%d%d%d%d", iBackground, iBackgroundComponent, iCentrality, iJetPt, iTrackPt));

            hBackgroundComponentRatio[iBackground][iBackgroundComponent][iCentralityReference][iJetPtReference][iTrackPtReference]->Divide(hBackgroundComponents[kMixedCone][iBackgroundComponent][iCentralityReference][iJetPtReference][iTrackPtReference]);
            
          } // Background component loop

          // ================================================================================ //
          //       Different subevent combinations for contructed background estimates        //
          // ================================================================================ //

          // TODO: Implementation  
          
        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // Background type loop
  
  // ==================== //
  // ===   Plotting   === //
  // ==================== //

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
  int color[9] = {kRed, kBlue, kGreen+2, kCyan, kMagenta, kOrange-1, kAzure+5, kOrange-2, kGray};
  int markerStyle[5] = {kFullCircle, kOpenSquare, kOpenCross, kFullStar, kFullCross};
  int firstNormalizationIndex, lastNormalizationIndex, currentNormalizationIndex;
  TH1D *referenceHistogram;
  
  // If optimal background matching is done, change the legend size
  double legendY1 = 0.09;
  double legendY2 = 0.94;
  
  // ================================================================================= //
  // ===   Draw comparison of different background components between estimators   === //
  // ================================================================================= //
  
  for(auto jetPtBin : comparedJetPtBin){
    jetPtString = Form("%.0f < jet p_{T} < %.0f", jetPtBin.first, jetPtBin.second);
    compactJetPtString = Form("_J=%.0f-%.0f", jetPtBin.first, jetPtBin.second);
    iJetPt = card[0]->FindBinIndexJetPtEEC(jetPtBin);       
    for(auto trackPtBin : comparedTrackPtBin){
      trackPtString = Form("%.1f < track p_{T}", trackPtBin);
      compactTrackPtString = Form("_T%.1f", trackPtBin);
      compactTrackPtString.ReplaceAll(".","v");
      iTrackPt = card[0]->GetBinIndexTrackPtEEC(trackPtBin);
      for(auto centralityBin : comparedCentralityBin){
        
        // Set the centrality information for legends and figure saving
        centralityString = Form("Cent: %.0f-%.0f%%", centralityBin.first, centralityBin.second);
        compactCentralityString = Form("_C=%.0f-%.0f", centralityBin.first, centralityBin.second);
        iCentrality = card[0]->FindBinIndexCentrality(centralityBin);

        // Loop over different background components
        for(int iBackgroundComponent = 0; iBackgroundComponent < knBackgroundComponents; iBackgroundComponent++){
            
          // Create a legend for the figure
          legend = new TLegend(0.19,legendY1,0.39,legendY2);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          legend->AddEntry((TObject*) 0, centralityString.Data(),"");
          legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
          legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
            
          // Create a new canvas for the plot
          drawer->CreateSplitCanvas();
            
          // Logarithmic EEC axis
          if(logEEC) drawer->SetLogY(true);

          // Set drawing style for background components
          for(int iBackground = 0; iBackground < knBackgroundTypes; iBackground++){
            hBackgroundComponents[iBackground][iBackgroundComponent][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iBackground]);
            hBackgroundComponents[iBackground][iBackgroundComponent][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iBackground]);
            hBackgroundComponents[iBackground][iBackgroundComponent][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iBackground]);
          }
            
          // Set the proper drawing range in DeltaR axis
          hBackgroundComponents[0][iBackgroundComponent][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
            
          // Draw the background components to the upper canvas
          for(int iBackground = 0; iBackground < knBackgroundTypes; iBackground++){
            if(iBackground == 0){
              drawer->DrawHistogramToUpperPad(hBackgroundComponents[iBackground][iBackgroundComponent][iCentrality][iJetPt][iTrackPt], "#Deltar", Form("EEC %s", backgroundComponentName[iBackgroundComponent].Data()), " ");
            } else {
              hBackgroundComponents[iBackground][iBackgroundComponent][iCentrality][iJetPt][iTrackPt]->Draw("same");
            }
            legend->AddEntry(hBackgroundComponents[iBackground][iBackgroundComponent][iCentrality][iJetPt][iTrackPt], backgroundName[iBackground].Data(), "p");
          }

          // Draw the legend to the upper pad
          legend->Draw();
            
          // Linear scale for the ratio
          drawer->SetLogY(false);

          // Set drawing style for the ratio histograms
          for(int iBackground = 0; iBackground < knBackgroundTypes; iBackground++){
            hBackgroundComponentRatio[iBackground][iBackgroundComponent][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iBackground]);
            hBackgroundComponentRatio[iBackground][iBackgroundComponent][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iBackground]);
            hBackgroundComponentRatio[iBackground][iBackgroundComponent][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iBackground]);
          }
  
          // Set the drawing range for the ratio
          hBackgroundComponentRatio[kPerpendicularCone][iBackgroundComponent][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);

          drawer->DrawHistogramToLowerPad(hBackgroundComponentRatio[kPerpendicularCone][iBackgroundComponent][iCentrality][iJetPt][iTrackPt], "#Deltar", "#frac{Perpendicular cone}{Mixed cone}", " ");
            
          // Save the figures to a file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/%sComparison%s_%s%s%s.%s", backgroundComponentSaveName[iBackgroundComponent].Data(), saveComment, compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
          }
            
        } // Background component loop
      } // Centrality loop
    } // Track pT loop
  } // Jet pT
  
}
