#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"
#include "AlgorithmLibrary.h"
#include "SystematicUncertaintyOrganizer.h"

/*
 * Macro for making closure plots for the analysis. It compares fully unfolded results to signal from MC truth.
 */
void fullAnalysisClosure(){

  // Enumeration for distribution type
  enum enumDistributionType{kMeasured, kTruth, kNDistributionTypes};
  bool isPbPbData = false;
  const int nSplits = isPbPbData ? 2 : 3;
  const int weightExponent = 1;

  // Ensure that a reasonable weight exponent is selected
  if(weightExponent < 1 || weightExponent > 2){
    cout << "ERROR! You have selected weightExponent = " << weightExponent << " while only 1 and 2 are implemented" << endl;
    cout << "Please select either 1 or 2 for the weightExponent!" << endl;
    return;
  }

  // Open the input files
  TString fileName[kNDistributionTypes][nSplits][2];

  if(isPbPbData){
    fileName[kMeasured][0][0] = "data/PbPbMC2018_RecoReco_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_nominalSmear_newTrackPairEfficiencySmoothed_part2_processed_2023-07-15.root";
    fileName[kMeasured][1][0] = "data/PbPbMC2018_RecoReco_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_nominalSmear_newTrackPairEfficiencySmoothed_part1_processed_2023-07-15.root";

    fileName[kMeasured][0][1] = "data/PbPbMC2018_RecoReco_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_energyWeightSquared_nominalSmear_part2_processed_2023-12-11.root";
    fileName[kMeasured][1][1] = "data/PbPbMC2018_RecoReco_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_energyWeightSquared_nominalSmear_part1_processed_2023-12-11.root";

    fileName[kTruth][0][0] = "data/PbPbMC2018_GenGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_nominalSmear_truthReference_part2_processed_2023-07-11.root";
    fileName[kTruth][1][0] = "data/PbPbMC2018_GenGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_nominalSmear_truthReference_part1_processed_2023-07-11.root";

    fileName[kTruth][0][1] = "data/PbPbMC2018_GenGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_energyWeightSquared_nominalSmear_truthReference_part2_processed_2023-10-24.root";
    fileName[kTruth][1][1] = "data/PbPbMC2018_GenGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_energyWeightSquared_nominalSmear_truthReference_part1_processed_2023-10-24.root";

  } else {
    fileName[kMeasured][0][0] = "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_32deltaRBins_newTrackPairEfficiency_nominalSmear_part2_processed_2023-07-15.root";
    fileName[kMeasured][1][0] = "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_32deltaRBins_newTrackPairEfficiency_nominalSmear_part1_processed_2023-07-15.root";
    fileName[kMeasured][2][0] = "data/ppMC2017_RecoReco_Herwig_pfJets_wtaAxis_32deltaRBins_newTrackPairEfficiency_nominalSmear_processed_2023-07-15.root";

    fileName[kMeasured][0][1] = "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_energyWeightSquared_nominalSmear_part2_processed_2023-11-10.root";
    fileName[kMeasured][1][1] = "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_energyWeightSquared_nominalSmear_part1_processed_2023-11-10.root";
    fileName[kMeasured][2][1] = "data/ppMC2017_RecoReco_Herwig_pfJets_wtaAxis_energyWeightSquared_nominalSmear_processed_2023-11-10.root";

    fileName[kTruth][0][0] = "data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_32deltaRBins_nominalSmear_truthReference_part2_processed_2023-06-21.root";
    fileName[kTruth][1][0] = "data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_32deltaRBins_nominalSmear_truthReference_part1_processed_2023-06-21.root";
    fileName[kTruth][2][0] = "data/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_32deltaRBins_firstTest_processed_2023-07-10.root";

    fileName[kTruth][0][1] = "data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_nominalSmear_truthReference_part2_processed_2023-10-30.root";
    fileName[kTruth][1][1] = "data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_nominalSmear_truthReference_part1_processed_2023-10-30.root";
    fileName[kTruth][2][1] = "data/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_energyWeightSquared_nominalSmear_processed_2023-11-10.root";
  }

  // Uncertainty file. First index is pp file, second PbPb file. The code will automatically select the correct one
  TString uncertaintyFileName[2][2] = {
    // Uncertainty files for pT1*pT2 weight
    {"systematicUncertainties/systematicUncertainties_pp_noMCnonClosure_2023-12-13.root", "systematicUncertainties/systematicUncertainties_jetMetUpdate_2023-07-14.root"},
    // Uncertainty files for pT1^{2}*pT2^{2} weight. TODO: Update the file names!
    {"systematicUncertainties/systematicUncertainties_pp_energyWeightSquared_noMCnonClosure_2023-12-13.root", "systematicUncertainties/systematicUncertainties_PbPb_energyWeightSquared_noMCnonClosure_2023-12-12.root"}
  };

  TFile* inputFile[kNDistributionTypes][nSplits];
  TFile* uncertaintyFile;
  EECCard* card[kNDistributionTypes][nSplits];
  EECCard* uncertaintyCard;

  for(int iFile = 0; iFile < kNDistributionTypes; iFile++){
    for(int iSplit = 0; iSplit < nSplits; iSplit++){
      inputFile[iFile][iSplit] = TFile::Open(fileName[iFile][iSplit][weightExponent-1]);
  
      // Check that the files exist
      if(inputFile[iFile][iSplit] == NULL){
        cout << "Error! The file " << fileName[iFile][iSplit][weightExponent-1].Data() << " does not exist!" << endl;
        cout << "Maybe you forgot the data/ folder path?" << endl;
        cout << "Will not execute the code" << endl;
        return;
      }

      card[iFile][iSplit] = new EECCard(inputFile[iFile][iSplit]);
    }
  }

  // File for systematic uncertainties. This is for drawing a band to ratio showing the relevant uncertainties for this comparison
  uncertaintyFile = TFile::Open(uncertaintyFileName[weightExponent-1][isPbPbData]);
  if(uncertaintyFile == NULL){
    cout << "Error! The file " << uncertaintyFileName[weightExponent-1][isPbPbData].Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the systematicUncertainties/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }
  uncertaintyCard = new EECCard(uncertaintyFile);
  
  // It is assumed that the different splits have the same binning. It might be worth implementing a check here to avoid bugs producing scary closures.

  // Declare algortihm library object for making stuff easier
  AlgorithmLibrary* optimusPrimeTheTransformer = new AlgorithmLibrary();
  TString today = optimusPrimeTheTransformer->GetToday();
  
  // ====================================================
  //               Binning configuration
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = card[kMeasured][0]->GetNCentralityBins();
  const int nJetPtBinsEEC = card[kMeasured][0]->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card[kMeasured][0]->GetNTrackPtBinsEEC();
  
  // Default binning ranges for reference
  // centrality = {0,10,30,50,90}
  // jet pT in energy-energy correlators = {120,140,160,180,200,300,500,5030}
  // track pT in energy-energy correlators = {0.7,1,1.5,2,2.5,3,3.5,4,300}
  
  // Draw the analysis closure for all bins that have been unfolded
  int firstDrawnCentralityBin = card[kMeasured][0]->GetFirstUnfoldedCentralityBin();
  int lastDrawnCentralityBin = card[kMeasured][0]->GetLastUnfoldedCentralityBin();
  
  int firstDrawnJetPtBinEEC = card[kMeasured][0]->GetFirstUnfoldedJetPtBin();
  int lastDrawnJetPtBinEEC = card[kMeasured][0]->GetLastUnfoldedJetPtBin();
  
  int firstDrawnTrackPtBinEEC = card[kMeasured][0]->GetFirstUnfoldedTrackPtBin();
  int lastDrawnTrackPtBinEEC = card[kMeasured][0]->GetLastUnfoldedTrackPtBin();
  
  // Select the types of energy-energy correlators are studied
  bool studyEnergyEnergyCorrelator[EECHistogramManager::knEnergyEnergyCorrelatorTypes];
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator] = true;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationPlus] = false;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationMinus] = false;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationPlus] = false;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationMinus] = false;
  
  // Logarithmic axes
  const bool logDeltaR = true;
  const bool logEEC = true;
  
  // Axis zooming
  std::pair<double,double> ratioZoom = std::make_pair(0.7, 1.3);

  // Legend text referring to splits
  TString splitLegendText[nSplits];
  if(isPbPbData){
    splitLegendText[0] = "split 1";
    splitLegendText[1] = "split 2";
  } else {
    splitLegendText[0] = "Pythia8, split 1";
    splitLegendText[1] = "Pythia8, split 2";
    splitLegendText[2] = "Herwig unfolded with Pythia";
  }
  
  // Figure saving
  const bool saveFigures = true;  // Save figures
  const char* saveComment = "_Pythia8_trackSelectionUpdate";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures

  // Save output file for Monte Carlo non-closure uncertainty
  const bool saveMonteCarloNonClosureFile = true;
  TString outputFileName[2][2] = {
    // Output file names for pT1*pT2 weight
    {Form("systematicUncertainties/monteCarloNonClosureRelative_pp_%s.root", today.Data()), Form("systematicUncertainties/monteCarloNonClosureRelative_PbPb_%s.root", today.Data())},
    // Output file names for pT1^{2}*pT2^{2} weight
    {Form("systematicUncertainties/monteCarloNonClosureRelative_pp_energyWeightSquared_%s.root", today.Data()), Form("systematicUncertainties/monteCarloNonClosureRelative_PbPb_energyWeightSquared_%s.root", today.Data())}
  };
  
  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager* histograms[kNDistributionTypes][nSplits];
  for(int iFile = 0; iFile < kNDistributionTypes; iFile++){
    for(int iSplit = 0; iSplit < nSplits; iSplit++){
      histograms[iFile][iSplit] = new EECHistogramManager(inputFile[iFile][iSplit], card[iFile][iSplit]);

      // Choose the energy-energy correlator types to load
      histograms[iFile][iSplit]->SetLoadEnergyEnergyCorrelators(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator]);
      histograms[iFile][iSplit]->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationPlus(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationPlus]);
      histograms[iFile][iSplit]->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationMinus(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationMinus]);
      histograms[iFile][iSplit]->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationPlus(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationPlus]);
      histograms[iFile][iSplit]->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationMinus(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationMinus]);

      // Choose the bin ranges
      histograms[iFile][iSplit]->SetCentralityBinRange(0, card[iFile][iSplit]->GetNCentralityBins() - 1);
      histograms[iFile][iSplit]->SetJetPtBinRangeEEC(0, card[iFile][iSplit]->GetNJetPtBinsEEC() - 1);
      histograms[iFile][iSplit]->SetTrackPtBinRangeEEC(0, card[iFile][iSplit]->GetNTrackPtBinsEEC() - 1);

      // Load the histograms from the file
      histograms[iFile][iSplit]->LoadProcessedHistograms();
    }
  }

  // Create systematic uncertainty organizer to illustrate the uncertainties on closures
  SystematicUncertaintyOrganizer* uncertaintyOrganizer = new SystematicUncertaintyOrganizer(uncertaintyFile);

  // Energy-energy correlator histograms
  TH1D* hEnergyEnergyCorrelatorSignal[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][kNDistributionTypes][nSplits]; // Last bin, true signal/extracted signal
  TH1D* hEnergyEnergyCorrelatorSignalRatio[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][nSplits]; // Ratio between true and extracted signal
  TH1D* hEnergyEnergyCorrelatorRelativeUncertainty[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC]; // Relative uncertainties for energy-energy correlators
  TH1D* hMonteCarloNonClosureUncertainty[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC]; // Relative Monte Carlo non-closure uncertaintiess estimated from difference between MC non-closure and relevant uncertainties
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
          hEnergyEnergyCorrelatorRelativeUncertainty[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = NULL;
          hMonteCarloNonClosureUncertainty[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = NULL;
          for(int iSplit = 0; iSplit < nSplits; iSplit++){
            for(int iSignalType = 0; iSignalType < kNDistributionTypes; iSignalType++){
              hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSignalType][iSplit] = NULL;
            } // Signal type loop
            hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSplit] = NULL;
          } // Split loop
        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
  } // Energy-energy correlator type loop
  
  // Helper histograms
  TH1D *helperHistogram;
  double normalizationFactor;
  std::pair<double,double> drawingRange = std::make_pair(0.008, 0.39);
  double epsilon = 0.0001;
  int lowSignalBin, highSignalBin;
  int iCentralityTruth, iJetPtTruth, iTrackPtTruth;
  int truthIndex = isPbPbData ? EECHistograms::kPythiaPythia : EECHistograms::knSubeventCombinations;
  
  // Get the histograms from the histogram manager and normalize the signal histograms to one
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only load the selected energy-energy correlator types
    if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;

    for(int iSplit = 0; iSplit < nSplits; iSplit++){

      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        iCentralityTruth = card[kTruth][iSplit]->FindBinIndexCentrality(card[kMeasured][iSplit]->GetBinBordersCentrality(iCentrality));
        for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
          iTrackPtTruth = card[kTruth][iSplit]->FindBinIndexTrackPtEEC(card[kMeasured][iSplit]->GetBinBordersTrackPtEEC(iTrackPt));
          for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
            iJetPtTruth = card[kTruth][iSplit]->FindBinIndexJetPtEEC(card[kMeasured][iSplit]->GetBinBordersJetPtEEC(iJetPt));

            // Read the signal from the unfolded measurement
            hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kMeasured][iSplit] = histograms[kMeasured][iSplit]->GetHistogramEnergyEnergyCorrelatorProcessed(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

            // Read the signal from MC truth
            hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTruth][iSplit] = histograms[kTruth][iSplit]->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentralityTruth, iJetPtTruth, iTrackPtTruth, EECHistograms::kSameJetPair, truthIndex);

            // Normalize the signal distributions to one in the drawingRange
            lowSignalBin = hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kMeasured][iSplit]->GetXaxis()->FindBin(drawingRange.first + epsilon);
            highSignalBin = hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kMeasured][iSplit]->GetXaxis()->FindBin(drawingRange.second - epsilon);
            for(int iSignalType = 0; iSignalType < kNDistributionTypes; iSignalType++){
              hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSignalType][iSplit]->Scale(1 / hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSignalType][iSplit]->Integral(lowSignalBin, highSignalBin, "width"));
            }

            // Calculate the extracted to true signal ratio
            hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSplit] = (TH1D*)hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kMeasured][iSplit]->Clone(Form("signalRatio%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iTrackPt, iJetPt, iSplit));

            hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSplit]->Divide(hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTruth][iSplit]);

          } // Jet pT loop
        } // Track pT loop
      } // Centrality loop
    } // Split loop
  } // Energy-energy correlator type loop

  TH1D* backgroundHelper = NULL;
  TH1D* trackPairHelper = NULL;
  TH1D* trackSelectionHelper = NULL;
  std::pair<double,double> shiftedCentralityBin;
  double sumOfSquares;

  // Read the uncertainties, add relevant sources in quadrature, and transform them into relative uncertainties
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only load the selected energy-energy correlator types
    if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;

    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      if(isPbPbData){
        // Take into account that MC has shifted centrality bins but data does not
        shiftedCentralityBin = card[kMeasured][0]->GetBinBordersCentrality(iCentrality);
        iCentralityTruth = uncertaintyCard->FindBinIndexCentrality(shiftedCentralityBin.first-4, shiftedCentralityBin.second-4);
      } else {
        iCentralityTruth = 0;
      }
      for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
        iTrackPtTruth = uncertaintyCard->FindBinIndexTrackPtEEC(card[kMeasured][0]->GetBinBordersTrackPtEEC(iTrackPt));
        for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
          iJetPtTruth = uncertaintyCard->FindBinIndexJetPtEEC(card[kMeasured][0]->GetBinBordersJetPtEEC(iJetPt));

          backgroundHelper = uncertaintyOrganizer->GetSystematicUncertainty(iCentralityTruth, iJetPtTruth, iTrackPtTruth, SystematicUncertaintyOrganizer::kBackgroundSubtraction);
          trackPairHelper = uncertaintyOrganizer->GetSystematicUncertainty(iCentralityTruth, iJetPtTruth, iTrackPtTruth, SystematicUncertaintyOrganizer::kTrackPairEfficiency);
          trackSelectionHelper = uncertaintyOrganizer->GetSystematicUncertainty(iCentralityTruth, iJetPtTruth, iTrackPtTruth, SystematicUncertaintyOrganizer::kTrackSelection);
          hEnergyEnergyCorrelatorRelativeUncertainty[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = uncertaintyOrganizer->GetSystematicUncertainty(iCentralityTruth, iJetPtTruth, iTrackPtTruth, SystematicUncertaintyOrganizer::kJetEnergyScale);

          // Add all the relevant uncertainties in quadrature
          for(int iBin = 1; iBin <= hEnergyEnergyCorrelatorRelativeUncertainty[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetNbinsX(); iBin++){
            sumOfSquares = 0;
            sumOfSquares += TMath::Power(hEnergyEnergyCorrelatorRelativeUncertainty[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetBinError(iBin),2);
            sumOfSquares += TMath::Power(backgroundHelper->GetBinError(iBin),2);
            sumOfSquares += TMath::Power(trackPairHelper->GetBinError(iBin),2);
            sumOfSquares += TMath::Power(trackSelectionHelper->GetBinError(iBin),2);
            hEnergyEnergyCorrelatorRelativeUncertainty[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetBinError(iBin, TMath::Sqrt(sumOfSquares));
          }

          // Transform the uncertainty histogram into relative uncertainties
          optimusPrimeTheTransformer->TransformToRelativeUncertainty(hEnergyEnergyCorrelatorRelativeUncertainty[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt], true);

        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop

  } // Energy-energy correlator type loop

  // Calculate the uncertainty for Monte Carlo non-closure
  TString histogramNamer;
  double averageNonClosure;
  double relevantUncertainty;
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only load the selected energy-energy correlator types
    if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;

    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
        for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){

          if(isPbPbData){
            histogramNamer = Form("relativeUncertaintyMonteCarloClosure_C%dJ%dT%d", iCentrality, iJetPt, iTrackPt);
          } else {
            histogramNamer = Form("relativeUncertaintyMonteCarloClosure_pp_J%dT%d", iJetPt, iTrackPt);
          }

          hMonteCarloNonClosureUncertainty[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelatorRelativeUncertainty[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Clone(histogramNamer);

          // Add all the relevant uncertainties in quadrature
          for(int iBin = 1; iBin <= hMonteCarloNonClosureUncertainty[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetNbinsX(); iBin++){
            averageNonClosure = hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0]->GetBinContent(iBin);
            averageNonClosure += hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][1]->GetBinContent(iBin);
            averageNonClosure /= 2;
            relevantUncertainty = hEnergyEnergyCorrelatorRelativeUncertainty[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetBinError(iBin);
            if(TMath::Abs(1-averageNonClosure) > relevantUncertainty){
              hMonteCarloNonClosureUncertainty[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetBinContent(iBin, TMath::Abs(1-averageNonClosure) - relevantUncertainty);
            } else {
              hMonteCarloNonClosureUncertainty[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetBinContent(iBin, 0);
            }
            // For selected bins, suppress obvious oscillation
            if(isPbPbData && iCentrality == 2 && iJetPt >= 8 && iBin == 23){
              hMonteCarloNonClosureUncertainty[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetBinContent(iBin, 0);
            }
            hMonteCarloNonClosureUncertainty[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetBinError(iBin, 0);
          }

        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
  } // Energy-energy correlator type loop

  // ==========================================================================
  //        Draw the selected energy-energy correlator signal ratios
  // ==========================================================================
  
  JDrawer* drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  
  if(logDeltaR) drawer->SetLogX(true);

  TLegend* legend;
  TString centralityString, trackPtString, jetPtString;
  TString compactCentralityString, compactTrackPtString, compactJetPtString;
  TString legendString;

  int recoSplitColor[] = {kRed, kBlue, kGreen+3};
  int truthSplitColor[] = {kMagenta, kCyan, kBlack};
  
  // Loop over all selected histograms
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only draw the selected energy-energy correlator types
    if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
    
    // Loop over centrality bins
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      
      // Set the centrality information for legends and figure saving
      if(isPbPbData){
        centralityString = Form("Cent: %.0f-%.0f%%", histograms[kMeasured][0]->GetCentralityBinBorder(iCentrality), histograms[kMeasured][0]->GetCentralityBinBorder(iCentrality+1));
        compactCentralityString = Form("_C=%.0f-%.0f", histograms[kMeasured][0]->GetCentralityBinBorder(iCentrality), histograms[kMeasured][0]->GetCentralityBinBorder(iCentrality+1));
      } else {
        centralityString = "";
        compactCentralityString = "";
      }
      
      // Loop over jet pT bins
      for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
        
        // Set the jet pT information for legends and figure saving
        if(iJetPt == histograms[kMeasured][0]->GetNJetPtBinsEEC()){
          jetPtString = Form("Jet p_{T} > %.0f", card[kMeasured][0]->GetJetPtCut());
          compactJetPtString = "";
        } else {
          jetPtString = Form("%.0f < jet p_{T} < %.0f", histograms[kMeasured][0]->GetJetPtBinBorderEEC(iJetPt), histograms[kMeasured][0]->GetJetPtBinBorderEEC(iJetPt+1));
          compactJetPtString = Form("_J=%.0f-%.0f", histograms[kMeasured][0]->GetJetPtBinBorderEEC(iJetPt), histograms[kMeasured][0]->GetJetPtBinBorderEEC(iJetPt+1));
        }
        
        // Loop over track pT bins
        for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
          
          trackPtString = Form("%.1f < track p_{T}",histograms[kMeasured][0]->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString = Form("_T>%.1f",histograms[kMeasured][0]->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString.ReplaceAll(".","v");
          
          // Create a legend for the figure
          legend = new TLegend(0.18,0.04,0.45,0.44+0.04*weightExponent);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          legend->AddEntry((TObject*) 0, Form("%s 5.02 TeV",card[kMeasured][0]->GetAlternativeDataType().Data()), "");
          if(weightExponent > 1) legend->AddEntry((TObject*) 0, "Energy weight squared","");
          if(isPbPbData) legend->AddEntry((TObject*) 0, centralityString.Data(),"");
          legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
          legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
          
          // Create a new canvas for the plot
          drawer->CreateSplitCanvas();
          
          // Logarithmic EEC axis
          if(logEEC) drawer->SetLogY(true);
          
          // Set the x-axis drawing range
          hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTruth][0]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
          
          // Draw the histograms to the upper canvas
          hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTruth][0]->SetLineColor(truthSplitColor[0]);
          drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTruth][0], "#Deltar", "EEC Signal", " ");
          
          hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kMeasured][0]->SetLineColor(recoSplitColor[0]);
          hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kMeasured][0]->Draw("same");

          for(int iSplit = 1; iSplit < nSplits; iSplit++){
            hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTruth][iSplit]->SetLineColor(truthSplitColor[iSplit]);
            hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTruth][iSplit]->Draw("same");
            hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kMeasured][iSplit]->SetLineColor(recoSplitColor[iSplit]);
            hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kMeasured][iSplit]->Draw("same");
          }

          for(int iSplit = 0; iSplit < nSplits; iSplit++){
            legend->AddEntry(hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTruth][iSplit], Form("True signal, %s", splitLegendText[iSplit].Data()), "l");
            legend->AddEntry(hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kMeasured][iSplit], Form("Unfolded signal, %s", splitLegendText[iSplit].Data()), "l");
          }
          
          // Draw the legends to the upper pad
          legend->Draw();
          
          // Linear scale for the ratio
          drawer->SetLogY(false);
          
          // For the ratio, first draw the systematic uncertainty band
          drawer->SetGridY(true);
          hEnergyEnergyCorrelatorRelativeUncertainty[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
          hEnergyEnergyCorrelatorRelativeUncertainty[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
          hEnergyEnergyCorrelatorRelativeUncertainty[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(kRed, 0.2);
          hEnergyEnergyCorrelatorRelativeUncertainty[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);
          hEnergyEnergyCorrelatorRelativeUncertainty[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
          drawer->DrawHistogramToLowerPad(hEnergyEnergyCorrelatorRelativeUncertainty[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt], "#Deltar", "#frac{Unfolded}{True}", " ","e2");

          // Then draw all the ratios to the same canvas
          for(int iSplit = 0; iSplit < nSplits; iSplit++){
            hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSplit]->SetLineColor(recoSplitColor[iSplit]);
            hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSplit]->Draw("same");
          }
          drawer->SetGridY(false);
          
          // Save the figures to a file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/%sFullAnalysisClosure%s%s%s%s.%s", histograms[kMeasured][0]->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator), saveComment, compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
          }
          
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Energy-energy correlator type loop

  // If we are writing a file evaluating MC non-closure uncertainties, do it
  if(saveMonteCarloNonClosureFile){
    cout << "Saving to file: " << outputFileName[weightExponent-1][isPbPbData].Data() << endl;
    TFile* outputFile = new TFile(outputFileName[weightExponent-1][isPbPbData],"UPDATE");

    for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){

      // Only save the selected energy-energy correlator types
      if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;

      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
          for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
            hMonteCarloNonClosureUncertainty[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Write();
          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Energy-energy correlator type loop

    // Write the card from the file with which all the bins have been synchronized
    card[kMeasured][0]->Write(outputFile);

    outputFile->Close();
  }
}
