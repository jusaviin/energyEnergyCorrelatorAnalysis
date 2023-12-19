#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "AlgorithmLibrary.h"

/*
 * Macro for finding background normalization scale in MC. This is done by integrating the energy-energy correlator distributions for background within the signal cone, and in the reflected cone, and taking a ratio of these numbers.
 */
void getEECaverageForCovariance(){

  // File from which the per jet average values in each bin are determined
  TString inputFileName = "data/ppData_pfJets_wtaAxis_jet60or80triggers_finalResults_processed_2023-08-07.root";
  // eecAnalysis_akFlowJet_wtaAxis_energyWeightSquared_nominalReflectedCone_processed_2023-12-01.root
  // eecAnalysis_akFlowJet_wtaAxis_newTrackPairEfficiencySmoothed_unfoldingWithNominalSmear_processed_2023-07-13.root
  // ppData_pfJets_wtaAxis_energyWeightSquared_jet60or80triggers_firstFinalResults_processed_2023-10-26.root
  // ppData_pfJets_wtaAxis_jet60or80triggers_finalResults_processed_2023-08-07.root
  
  // Open the input file
  TFile* inputFile = TFile::Open(inputFileName);
  
  if(inputFile == NULL){
    cout << "Error! The file " << inputFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }
  
  // Load the card from the file and read the collision system
  EECCard* card = new EECCard(inputFile);
  TString collisionSystem = card->GetDataType();
  const int iWeightExponent = card->GetWeightExponent();

  // Algorithm library to find today's date automagically
  AlgorithmLibrary* timeKeeper = new AlgorithmLibrary();
  TString today = timeKeeper->GetToday();
  TString energyWeightString[] = {"","","_energyWeightSquared"}; 

  // User configuration
  TString outputFileName = Form("averageEECfileForCovariance_%s%s_%s.root", collisionSystem.Data(), energyWeightString[iWeightExponent].Data(), today.Data());
  TString fileOption = "UPDATE";
  
  // ====================================================
  //               Binning configuration
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = card->GetNCentralityBins();
  const int nJetPtBinsEEC = card->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card->GetNTrackPtBinsEEC();
  
  // Default binning ranges for reference
  // centrality = {0,10,30,50,90}
  // track pT =  {0.7,1,2,3,4,8,12,300}
  // jet pT in energy-energy correlators = {120,140,160,180,200,300,500,5030}
  // track pT in energy-energy correlators = {0.7,1,2,3,4,6,300}
  
  // Select the types of energy-energy correlators to determine per jet average
  bool integrateEnergyEnergyCorrelator[EECHistogramManager::knEnergyEnergyCorrelatorTypes];
  integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator] = true;
  integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationPlus] = false;
  integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationMinus] = false;
  integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationPlus] = false;
  integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationMinus] = false;
  
  // Determine an inxed of a drawn energy-energy correlator
  int lowestEnergyEnergyCorrelatorIndex = -1;
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    if(integrateEnergyEnergyCorrelator[iEnergyEnergyCorrelator]){
      lowestEnergyEnergyCorrelatorIndex = iEnergyEnergyCorrelator;
      break;
    }
  }
    
  // Create and setup a new histogram manager to project and handle the histograms
  EECHistogramManager* histograms = new EECHistogramManager(inputFile,card);
  
  // Choose the energy-energy correlator types to load
  histograms->SetLoadEnergyEnergyCorrelators(integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator]);
  histograms->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationPlus(integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationPlus]);
  histograms->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationMinus(integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationMinus]);
  histograms->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationPlus(integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationPlus]);
  histograms->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationMinus(integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationMinus]);
  
  // Choose the bin ranges
  histograms->SetCentralityBinRange(0,nCentralityBins);
  histograms->SetJetPtBinRangeEEC(0,nJetPtBinsEEC);
  histograms->SetTrackPtBinRangeEEC(0,nTrackPtBinsEEC);
  
  // Load the histograms from the file
  histograms->LoadProcessedHistograms();
  
  // Initialize the energy-energy correlator histogram array to NULL
  TH1D* hEnergyEnergyCorrelator[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = NULL;
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Energy-energy correlator type loop
  
  // Get the histograms from the histogram manager and calculate per jet averages for enegy-energy correlators
  double normalizationFactor;
  std::pair<double,double> jetPtBinBorders;
  double originalBinContent, originalBinError, binWidth;
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only read the selected energy-energy correlator types
    if(!integrateEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){

        // For normalization, read the number of jets in the given centrality and pT bin
        jetPtBinBorders = card->GetBinBordersJetPtEEC(iJetPt);
        normalizationFactor = histograms->GetJetPtIntegral(iCentrality, jetPtBinBorders.first, jetPtBinBorders.second);

        for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
            
          // Read the energy-energy correlator from file
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = histograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt);

          // Undo the bin width normalization
          for(int iBin = 1; iBin <= hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetNbinsX(); iBin++){
            originalBinContent = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetBinContent(iBin);
            originalBinError = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetBinError(iBin);
            binWidth = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetBinWidth(iBin);
            hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetBinContent(iBin, originalBinContent*binWidth);
            hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetBinError(iBin, originalBinError*binWidth);
          }

          // After bin width normalization is removed, normalize the bin contents to the number of jets
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Scale(1.0 / normalizationFactor);
                          
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Energy-energy correlator loop

  // After the histograms have been normalized, save everything to file
  
  // Create the output file
  TString histogramNamer;
  TFile* outputFile = new TFile(outputFileName,fileOption);

  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only write the selected energy-energy correlator types
    if(!integrateEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
    
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
          histogramNamer = Form("eecAverageForCovariance_C%dJ%dT%d", iCentrality, iJetPt, iTrackPt);
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Write(histogramNamer.Data(), TObject::kOverwrite);
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Energy-energy correlator loop

  // Write the binning information using JCard
  card->Write(outputFile);

  // Close the file after everything is written
  outputFile->Close();
  
  // Delete the outputFile object
  delete outputFile;

}
