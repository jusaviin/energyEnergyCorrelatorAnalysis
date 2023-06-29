#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Macro for doing a bottom line test for unfolding.
 * The idea here is that unfolding should not enhance the discriminating power for models
 * To check this, we compare the unfolded data with Monte Carlo truth and not-unfolded data with reconstructed Monte Carlo.
 * In both cases, we should see similar values for the ratio.
 * If anything, the differences in unfolded-to-truth comparison should be larger.
 */
void bottomLineTestForUnfolding(){

  enum enumDataType{kPbPb, kPp, kNDataTypes};
  enum enumMCType{kReconstructed, kTruth, kNMCTypes};
  enum enumComparisonType{kData, kSimulation, kNComparisonTypes};

  // ============= //
  // Configuration //
  // ============= //
  
  // Input files
  TString dataFileName[kNDataTypes];
  dataFileName[kPbPb] = "data/eecAnalysis_akFlowJet_wtaAxis_unfoldingWithNominalSmear_processed_2023-06-25.root";
  dataFileName[kPp] = "data/ppData_pfJets_wtaAxis_unfoldingWithNominalSmear_processed_2023-06-25.root";

  TString simulationFileName[kNMCTypes][kNDataTypes];
  simulationFileName[kReconstructed][kPbPb] = "data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_nominalSmear_reconstructedReference_processed_2023-06-23.root";
  simulationFileName[kTruth][kPbPb] = "data/PbPbMC2018_GenGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_nominalSmear_truthReference_processed_2023-06-23.root";
  simulationFileName[kReconstructed][kPp] = "data/ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_32deltaRBins_nominalSmear_reconstructedReference_processed_2023-06-21.root";
  simulationFileName[kTruth][kPp] = "data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_32deltaRBins_nominalSmear_truthReference_processed_2023-06-21.root";

  
  TFile* dataFile[kNDataTypes];
  TFile* simulationFile[kNMCTypes][kNDataTypes];
  EECCard* dataCard[kNDataTypes];
  EECCard* simulationCard[kNMCTypes][kNDataTypes];
  for(int iFile = 0; iFile < kNDataTypes; iFile++){

    // Load the input file
    dataFile[iFile] = TFile::Open(dataFileName[iFile]);

    // Check that the input file exists
    if(dataFile[iFile] == NULL){
      cout << "Error! The file " << dataFileName[iFile].Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }
    
    // Load the card from the file
    dataCard[iFile] = new EECCard(dataFile[iFile]);

    for(int iMCType = 0; iMCType < kNMCTypes; iMCType++){

      // Load the input file
      simulationFile[iMCType][iFile] = TFile::Open(simulationFileName[iMCType][iFile]);

      // Check that the input file exists
      if(simulationFile[iMCType][iFile] == NULL){
        cout << "Error! The file " << simulationFileName[iMCType][iFile].Data() << " does not exist!" << endl;
        cout << "Maybe you forgot the data/ folder path?" << endl;
        cout << "Will not execute the code" << endl;
        return;
      }
    
      // Load the card from the file
      simulationCard[iMCType][iFile] = new EECCard(simulationFile[iMCType][iFile]);

    } // MC type loop
  } // File loop for opening files

  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Find the number of bins from the PbPb card
  const int nCentralityBins = dataCard[kPbPb]->GetNCentralityBins();
  const int nJetPtBinsEEC = dataCard[kPbPb]->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = dataCard[kPbPb]->GetNTrackPtBinsEEC();
  
  // Only compare the bins that are unfolded
  int firstDrawnCentralityBin = dataCard[kPbPb]->GetFirstUnfoldedCentralityBin();
  int lastDrawnCentralityBin = dataCard[kPbPb]->GetLastUnfoldedCentralityBin();
  
  int firstDrawnJetPtBinEEC = dataCard[kPbPb]->GetFirstUnfoldedJetPtBin();
  int lastDrawnJetPtBinEEC = dataCard[kPbPb]->GetLastUnfoldedJetPtBin();
  
  int firstDrawnTrackPtBinEEC = dataCard[kPbPb]->GetFirstUnfoldedTrackPtBin();
  int lastDrawnTrackPtBinEEC = dataCard[kPbPb]->GetLastUnfoldedTrackPtBin();
  
  // Save the final plots
  const bool saveFigures = false;
  TString saveComment = "_firstLook";

  // Ratio zoom settings
  std::pair<double, double> analysisDeltaR = std::make_pair(0.006, 0.39); // DeltaR span in which the analysis is done
  std::pair<double, double> ratioZoom = std::make_pair(0.4, 1.6);

  // =============================================== //
  // Read the histograms from the histogram managers //
  // =============================================== //

  // Create histogram managers for result files and systematic uncertainty organizers for systematic uncertainty files
  EECHistogramManager* dataHistograms[kNDataTypes];
  EECHistogramManager* simulationHistograms[kNMCTypes][kNDataTypes];
  int firstLoadedBin, lastLoadedBin;
  std::pair<double,double> shiftedCentralityBin;
  
  for(int iDataType = 0; iDataType < kNDataTypes; iDataType++){

    // Create a new histogram manager
    dataHistograms[iDataType] = new EECHistogramManager(dataFile[iDataType], dataCard[iDataType]);
  
    // Load all unfolded energy-energy correlators
    dataHistograms[iDataType]->SetLoadEnergyEnergyCorrelators(true);
    dataHistograms[iDataType]->SetCentralityBinRange(dataCard[iDataType]->GetFirstUnfoldedCentralityBin(), dataCard[iDataType]->GetLastUnfoldedCentralityBin());
    dataHistograms[iDataType]->SetTrackPtBinRangeEEC(dataCard[iDataType]->GetFirstUnfoldedTrackPtBin(), dataCard[iDataType]->GetLastUnfoldedTrackPtBin());
    dataHistograms[iDataType]->SetJetPtBinRangeEEC(dataCard[iDataType]->GetFirstUnfoldedJetPtBin(), dataCard[iDataType]->GetLastUnfoldedJetPtBin());

    // Load the histograms from the file
    dataHistograms[iDataType]->LoadProcessedHistograms();

    // Loop also over MC types
    for(int iMCType = 0; iMCType < kNMCTypes; iMCType++){

      // Create a new histogram manager
      simulationHistograms[iMCType][iDataType] = new EECHistogramManager(simulationFile[iMCType][iDataType], simulationCard[iMCType][iDataType]);

      // Load all unfolded energy-energy correlators for the same bins as done for data
      simulationHistograms[iMCType][iDataType]->SetLoadEnergyEnergyCorrelators(true);

      shiftedCentralityBin = dataCard[iDataType]->GetBinBordersCentrality(dataCard[iDataType]->GetFirstUnfoldedCentralityBin());
      shiftedCentralityBin.first += 4;
      shiftedCentralityBin.second += 4;
      firstLoadedBin = simulationCard[iMCType][iDataType]->FindBinIndexCentrality(shiftedCentralityBin);
      shiftedCentralityBin = dataCard[iDataType]->GetBinBordersCentrality(dataCard[iDataType]->GetLastUnfoldedCentralityBin());
      shiftedCentralityBin.first += 4;
      shiftedCentralityBin.second += 4;
      lastLoadedBin = simulationCard[iMCType][iDataType]->FindBinIndexCentrality(shiftedCentralityBin);
      simulationHistograms[iMCType][iDataType]->SetCentralityBinRange(firstLoadedBin, lastLoadedBin);

      firstLoadedBin = simulationCard[iMCType][iDataType]->FindBinIndexJetPtEEC(dataCard[iDataType]->GetBinBordersJetPtEEC(dataCard[iDataType]->GetFirstUnfoldedJetPtBin()));
      lastLoadedBin = simulationCard[iMCType][iDataType]->FindBinIndexJetPtEEC(dataCard[iDataType]->GetBinBordersJetPtEEC(dataCard[iDataType]->GetLastUnfoldedJetPtBin()));
      simulationHistograms[iMCType][iDataType]->SetJetPtBinRangeEEC(firstLoadedBin, lastLoadedBin);

      firstLoadedBin = simulationCard[iMCType][iDataType]->FindBinIndexTrackPtEEC(dataCard[iDataType]->GetBinBordersTrackPtEEC(dataCard[iDataType]->GetFirstUnfoldedTrackPtBin()));
      lastLoadedBin = simulationCard[iMCType][iDataType]->FindBinIndexTrackPtEEC(dataCard[iDataType]->GetBinBordersTrackPtEEC(dataCard[iDataType]->GetLastUnfoldedTrackPtBin()));
      simulationHistograms[iMCType][iDataType]->SetTrackPtBinRangeEEC(firstLoadedBin, lastLoadedBin);

      // Load the simulation histograms from the file
      simulationHistograms[iMCType][iDataType]->LoadProcessedHistograms();
      
    }

  }
 
  TH1D* energyEnergyCorrelatorUnfoldedPbPb[kNComparisonTypes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorRawPbPb[kNComparisonTypes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorUnfoldedPp[kNComparisonTypes][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorRawPp[kNComparisonTypes][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorUnfoldedPbPbRatio[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorRawPbPbRatio[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorUnfoldedPpRatio[nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorRawPpRatio[nJetPtBinsEEC][nTrackPtBinsEEC];

  // Initialize histograms to NULL
  for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
      energyEnergyCorrelatorUnfoldedPpRatio[iJetPt][iTrackPt] = NULL;
      energyEnergyCorrelatorRawPpRatio[iJetPt][iTrackPt] = NULL;
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        energyEnergyCorrelatorUnfoldedPbPbRatio[iCentrality][iJetPt][iTrackPt] = NULL;
        energyEnergyCorrelatorRawPbPbRatio[iCentrality][iJetPt][iTrackPt] = NULL;
        for(int iDataType = 0; iDataType < kNComparisonTypes; iDataType++){
          energyEnergyCorrelatorUnfoldedPbPb[iDataType][iCentrality][iJetPt][iTrackPt] = NULL;
          energyEnergyCorrelatorRawPbPb[iDataType][iCentrality][iJetPt][iTrackPt] = NULL;
        }
      } // Centrality loop
      for(int iDataType = 0; iDataType < kNComparisonTypes; iDataType++){
        energyEnergyCorrelatorUnfoldedPp[iDataType][iJetPt][iTrackPt] = NULL;
        energyEnergyCorrelatorRawPp[iDataType][iJetPt][iTrackPt] = NULL;
      }
    } // Track pT loop
  } // Jet pT loop

  // Read the histograms from managers
  int iTrackPtMatched;
  int iJetPtMatched;
  int iCentralityMatched;

  for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){

    for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){

      // Match the pp binning with PbPb binning
      iJetPtMatched = dataCard[kPp]->FindBinIndexJetPtEEC(dataCard[kPbPb]->GetBinBordersJetPtEEC(iJetPt));
      iTrackPtMatched = dataCard[kPp]->FindBinIndexTrackPtEEC(dataCard[kPbPb]->GetBinBordersTrackPtEEC(iTrackPt));

      // Read the pp histograms that do not have centrality binning
      energyEnergyCorrelatorUnfoldedPp[kData][iJetPt][iTrackPt] = dataHistograms[kPp]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatched, iTrackPtMatched, EECHistogramManager::kEnergyEnergyCorrelatorUnfolded);
      energyEnergyCorrelatorRawPp[kData][iJetPt][iTrackPt] = dataHistograms[kPp]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatched, iTrackPtMatched);

      // Match the truth level Pythia8 binning with PbPb binning
      iJetPtMatched = simulationCard[kTruth][kPp]->FindBinIndexJetPtEEC(dataCard[kPbPb]->GetBinBordersJetPtEEC(iJetPt));
      iTrackPtMatched = simulationCard[kTruth][kPp]->FindBinIndexTrackPtEEC(dataCard[kPbPb]->GetBinBordersTrackPtEEC(iTrackPt));

      // Read the histograms from truth level Pythia8 simulation
      energyEnergyCorrelatorUnfoldedPp[kSimulation][iJetPt][iTrackPt] = simulationHistograms[kTruth][kPp]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatched, iTrackPtMatched);

      // Match the reconstructed level Pythia8 binning with PbPb binning
      iJetPtMatched = simulationCard[kReconstructed][kPp]->FindBinIndexJetPtEEC(dataCard[kPbPb]->GetBinBordersJetPtEEC(iJetPt));
      iTrackPtMatched = simulationCard[kReconstructed][kPp]->FindBinIndexTrackPtEEC(dataCard[kPbPb]->GetBinBordersTrackPtEEC(iTrackPt));

      // Read the histograms from truth level Pythia8 simulation
      energyEnergyCorrelatorRawPp[kSimulation][iJetPt][iTrackPt] = simulationHistograms[kReconstructed][kPp]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatched, iTrackPtMatched);

      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){

        // Read the PbPb histograms
        energyEnergyCorrelatorUnfoldedPbPb[kData][iCentrality][iJetPt][iTrackPt] = dataHistograms[kPbPb]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistogramManager::kEnergyEnergyCorrelatorUnfolded);
        energyEnergyCorrelatorRawPbPb[kData][iCentrality][iJetPt][iTrackPt] = dataHistograms[kPbPb]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt);

        // Match the truth level Pythia+Hydjet binning with PbPb binning
        iJetPtMatched = simulationCard[kTruth][kPbPb]->FindBinIndexJetPtEEC(dataCard[kPbPb]->GetBinBordersJetPtEEC(iJetPt));
        iTrackPtMatched = simulationCard[kTruth][kPbPb]->FindBinIndexTrackPtEEC(dataCard[kPbPb]->GetBinBordersTrackPtEEC(iTrackPt));
        shiftedCentralityBin = dataCard[kPbPb]->GetBinBordersCentrality(iCentrality);
        shiftedCentralityBin.first += 4;
        shiftedCentralityBin.second += 4;
        iCentralityMatched = simulationCard[kTruth][kPbPb]->FindBinIndexCentrality(shiftedCentralityBin);

        // Read the histograms from truth level Pythia+Hydjet simulation
        energyEnergyCorrelatorUnfoldedPbPb[kSimulation][iCentrality][iJetPt][iTrackPt] = simulationHistograms[kTruth][kPbPb]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentralityMatched, iJetPtMatched, iTrackPtMatched);

        // Match the reconstructed level Pythia+Hydjet binning with PbPb binning
        iJetPtMatched = simulationCard[kReconstructed][kPbPb]->FindBinIndexJetPtEEC(dataCard[kPbPb]->GetBinBordersJetPtEEC(iJetPt));
        iTrackPtMatched = simulationCard[kReconstructed][kPbPb]->FindBinIndexTrackPtEEC(dataCard[kPbPb]->GetBinBordersTrackPtEEC(iTrackPt));
        shiftedCentralityBin = dataCard[kPbPb]->GetBinBordersCentrality(iCentrality);
        shiftedCentralityBin.first += 4;
        shiftedCentralityBin.second += 4;
        iCentralityMatched = simulationCard[kReconstructed][kPbPb]->FindBinIndexCentrality(shiftedCentralityBin);

        // Read the histograms from truth level Pythia+Hydjet simulation
        energyEnergyCorrelatorRawPbPb[kSimulation][iCentrality][iJetPt][iTrackPt] = simulationHistograms[kReconstructed][kPbPb]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentralityMatched, iJetPtMatched, iTrackPtMatched);

      } // Centrality loop
    } // Track pT loop
  } // Jet pT loop

  // ================================================ //
  //   Normalize the histograms and calculate ratios  //
  // ================================================ //

  double epsilon = 0.0001;
  int lowAnalysisBin = energyEnergyCorrelatorUnfoldedPbPb[kData][firstDrawnCentralityBin][firstDrawnJetPtBinEEC][firstDrawnTrackPtBinEEC]->GetXaxis()->FindBin(analysisDeltaR.first + epsilon);
  int highAnalysisBin = energyEnergyCorrelatorUnfoldedPbPb[kData][firstDrawnCentralityBin][firstDrawnJetPtBinEEC][firstDrawnTrackPtBinEEC]->GetXaxis()->FindBin(analysisDeltaR.second - epsilon);
  for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
    for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
      for(int iDataType = 0; iDataType < kNComparisonTypes; iDataType++){
        energyEnergyCorrelatorUnfoldedPp[iDataType][iJetPt][iTrackPt]->Scale(1.0 / energyEnergyCorrelatorUnfoldedPp[iDataType][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
        energyEnergyCorrelatorRawPp[iDataType][iJetPt][iTrackPt]->Scale(1.0 / energyEnergyCorrelatorRawPp[iDataType][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
      }

      // Ratio calculation for pp
      energyEnergyCorrelatorUnfoldedPpRatio[iJetPt][iTrackPt] = (TH1D*) energyEnergyCorrelatorUnfoldedPp[kData][iJetPt][iTrackPt]->Clone(Form("unfoldedPpRatio%d%d", iJetPt, iTrackPt));
      energyEnergyCorrelatorUnfoldedPpRatio[iJetPt][iTrackPt]->Divide(energyEnergyCorrelatorUnfoldedPp[kSimulation][iJetPt][iTrackPt]);

      energyEnergyCorrelatorRawPpRatio[iJetPt][iTrackPt] = (TH1D*) energyEnergyCorrelatorRawPp[kData][iJetPt][iTrackPt]->Clone(Form("rawPpRatio%d%d", iJetPt, iTrackPt));
      energyEnergyCorrelatorRawPpRatio[iJetPt][iTrackPt]->Divide(energyEnergyCorrelatorRawPp[kSimulation][iJetPt][iTrackPt]);

      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        for(int iDataType = 0; iDataType < kNComparisonTypes; iDataType++){
          energyEnergyCorrelatorUnfoldedPbPb[iDataType][iCentrality][iJetPt][iTrackPt]->Scale(1.0 / energyEnergyCorrelatorUnfoldedPbPb[iDataType][iCentrality][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
          energyEnergyCorrelatorRawPbPb[iDataType][iCentrality][iJetPt][iTrackPt]->Scale(1.0 / energyEnergyCorrelatorRawPbPb[iDataType][iCentrality][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
        }

        // Ratio calculation for PbPb
        energyEnergyCorrelatorUnfoldedPbPbRatio[iCentrality][iJetPt][iTrackPt] = (TH1D*) energyEnergyCorrelatorUnfoldedPbPb[kData][iCentrality][iJetPt][iTrackPt]->Clone(Form("unfoldedPbPbRatio%d%d%d", iCentrality, iJetPt, iTrackPt));
        energyEnergyCorrelatorUnfoldedPbPbRatio[iCentrality][iJetPt][iTrackPt]->Divide(energyEnergyCorrelatorUnfoldedPbPb[kSimulation][iCentrality][iJetPt][iTrackPt]);

        energyEnergyCorrelatorRawPbPbRatio[iCentrality][iJetPt][iTrackPt] = (TH1D*) energyEnergyCorrelatorRawPbPb[kData][iCentrality][iJetPt][iTrackPt]->Clone(Form("rawPbPbRatio%d%d%d", iCentrality, iJetPt, iTrackPt));
        energyEnergyCorrelatorRawPbPbRatio[iCentrality][iJetPt][iTrackPt]->Divide(energyEnergyCorrelatorRawPbPb[kSimulation][iCentrality][iJetPt][iTrackPt]);

      } // Centrality loop
    } // Track pT loop
  } // Jet pT loop

  // ========================================= //
  //   Draw the histograms in separate plots   //
  // ========================================= //

  JDrawer* drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1, 1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  drawer->SetLogX(true);  // Logarithmic deltaR axis

  TString centralityString;
  TString compactCentralityString;
  TString jetPtString;
  TString compactJetPtString;
  TString trackPtString;
  TString compactTrackPtString;

  // Common variables for different plots
  TLegend* legend;

  // Draw plots for all studied bins
  for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++) {
    jetPtString = Form("%.0f < jet p_{T} < %.0f", dataCard[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), dataCard[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));
    compactJetPtString = Form("_J=%.0f-%.0f", dataCard[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), dataCard[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));

    for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++) {
      trackPtString = Form("%.1f < track p_{T}", dataCard[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));
      compactTrackPtString = Form("_T>%.1f", dataCard[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));
      compactTrackPtString.ReplaceAll(".", "v");

      // Setup the legend for plots
      legend = new TLegend(0.53, 0.6, 0.83, 0.85);
      legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);
      legend->AddEntry((TObject*)0, "pp", "");
      legend->AddEntry((TObject*)0, jetPtString.Data(), "");
      legend->AddEntry((TObject*)0, trackPtString.Data(), "");

      // Linear scale for the ratio
      drawer->SetLogY(false);

      // Set the axis drawing ranges
      energyEnergyCorrelatorRawPpRatio[iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
      energyEnergyCorrelatorRawPpRatio[iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

      // Set the style for histograms
      energyEnergyCorrelatorRawPpRatio[iJetPt][iTrackPt]->SetLineColor(kBlack);
      energyEnergyCorrelatorUnfoldedPpRatio[iJetPt][iTrackPt]->SetLineColor(kRed);

      drawer->SetGridY(true);
      drawer->DrawHistogram(energyEnergyCorrelatorRawPpRatio[iJetPt][iTrackPt], "#Deltar", "#frac{Data}{simulation}", " ");
      energyEnergyCorrelatorUnfoldedPpRatio[iJetPt][iTrackPt]->Draw("same");
      drawer->SetGridY(false);

      // Add the histograms to the legend
      legend->AddEntry(energyEnergyCorrelatorRawPpRatio[iJetPt][iTrackPt], "Raw ratio", "l");
      legend->AddEntry(energyEnergyCorrelatorUnfoldedPpRatio[iJetPt][iTrackPt], "Unfolded ratio", "l");

      // Draw the legend
      legend->Draw();

      // If a plot name is given, save the plot in a file
      if(saveFigures) {
        gPad->GetCanvas()->SaveAs(Form("figures/bottomLineCheck%s_pp%s%s.pdf", saveComment.Data(), compactJetPtString.Data(), compactTrackPtString.Data()));
      }

      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++) {
        centralityString = Form("Cent: %.0f-%.0f%%", dataCard[kPbPb]->GetLowBinBorderCentrality(iCentrality), dataCard[kPbPb]->GetHighBinBorderCentrality(iCentrality));
        compactCentralityString = Form("_C=%.0f-%.0f", dataCard[kPbPb]->GetLowBinBorderCentrality(iCentrality), dataCard[kPbPb]->GetHighBinBorderCentrality(iCentrality));

        // Setup the legend for plots
        legend = new TLegend(0.53, 0.6, 0.83, 0.85);
        legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);
        legend->AddEntry((TObject*)0, centralityString.Data(), "");
        legend->AddEntry((TObject*)0, jetPtString.Data(), "");
        legend->AddEntry((TObject*)0, trackPtString.Data(), "");

        // Linear scale for the ratio
        drawer->SetLogY(false);

        // Set the axis drawing ranges
        energyEnergyCorrelatorRawPbPbRatio[iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
        energyEnergyCorrelatorRawPbPbRatio[iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

        // Set the style for histograms
        energyEnergyCorrelatorRawPbPbRatio[iCentrality][iJetPt][iTrackPt]->SetLineColor(kBlack);
        energyEnergyCorrelatorUnfoldedPbPbRatio[iCentrality][iJetPt][iTrackPt]->SetLineColor(kRed);

        drawer->SetGridY(true);
        drawer->DrawHistogram(energyEnergyCorrelatorRawPbPbRatio[iCentrality][iJetPt][iTrackPt], "#Deltar", "#frac{Data}{simulation}", " ");
        energyEnergyCorrelatorUnfoldedPbPbRatio[iCentrality][iJetPt][iTrackPt]->Draw("same");
        drawer->SetGridY(false);

        // Add the histograms to the legend
        legend->AddEntry(energyEnergyCorrelatorRawPbPbRatio[iCentrality][iJetPt][iTrackPt], "Raw ratio", "l");
        legend->AddEntry(energyEnergyCorrelatorUnfoldedPbPbRatio[iCentrality][iJetPt][iTrackPt], "Unfolded ratio", "l");

        // Draw the legend
        legend->Draw();

        // If a plot name is given, save the plot in a file
        if(saveFigures) {
          gPad->GetCanvas()->SaveAs(Form("figures/bottomLineCheck%s%s%s%s.pdf", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data()));
        }
      }  // Centrality loop
    }    // Track pT loop
  }      // Jet pT loop
}
