#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"

/*
 * Macro for calculating integrals of the EEC distributions.
 * Main idea is to check signal/background in large DeltaR in MC.
 */
void maximumPtWithinJetConeStudy(){

  // File from which the integrals are calculated
  TString inputFileName = "data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJet_addMaxParticlePtInJet_wtaAxis_noTrigger_preprocessed_2022-11-09.root";
  
  // Open the input file
  TFile *inputFile = TFile::Open(inputFileName);
  
  if(inputFile == NULL){
    cout << "Error! The file " << inputFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }
  
  // Load the card from the file and read the collision system
  EECCard *card = new EECCard(inputFile);
  TString collisionSystem = card->GetDataType();
  
  // ====================================================
  //  Binning configuration for the integral calculation
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = card->GetNCentralityBins();
  const int nJetPtBinsEEC = card->GetNJetPtBinsEEC();
  
  // Bin range to be integrated
  int firstIntegratedCentralityBin = 0;
  int lastIntegratedCentralityBin = 0;
  
  int firstIntegratedJetPtBinEEC = 4;
  int lastIntegratedJetPtBinEEC = 4; // Note: Jets integrated over all pT ranges are in nJetPtBinsEEC bin
    
  // Create and setup a new histogram manager to project and handle the histograms
  EECHistogramManager *histograms = new EECHistogramManager(inputFile,card);
  
  // Choose the energy-energy correlator types to load
  histograms->SetLoadMaxParticlePtWithinJetCone(true);
  
  // Choose the bin ranges
  histograms->SetCentralityBinRange(firstIntegratedCentralityBin,lastIntegratedCentralityBin);
  histograms->SetJetPtBinRangeEEC(firstIntegratedJetPtBinEEC,lastIntegratedJetPtBinEEC);
  
  // Load the histograms from the file
  histograms->LoadProcessedHistograms();
  
  // Initialize the energy-energy correlator histogram array to NULL
  const int nStudiedTrackPtBins = 3;
  TH1D* hMaxParticlePtWithinJetCone[nCentralityBins][nJetPtBinsEEC+1][nStudiedTrackPtBins+1];
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC+1; iJetPt++){
        for(int iTrackPt = 0; iTrackPt < nStudiedTrackPtBins+1; iTrackPt++){
          hMaxParticlePtWithinJetCone[iCentrality][iJetPt][iTrackPt] = NULL;
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  
  // Get the histograms from the histogram manager and calculate integrals
  const int nIntegralThresholds = 3;
  double integralThresholds[] = {10,15,20};
  double integralValues[nCentralityBins][nJetPtBinsEEC+1][nStudiedTrackPtBins+1][nIntegralThresholds+1];
  int lowBinIndex, highBinIndex;
    
  for(int iCentrality = firstIntegratedCentralityBin; iCentrality <= lastIntegratedCentralityBin; iCentrality++){
    for(int iJetPt = firstIntegratedJetPtBinEEC; iJetPt <= lastIntegratedJetPtBinEEC; iJetPt++){
      
      // Calculate the integrals without track pT selection
      hMaxParticlePtWithinJetCone[iCentrality][iJetPt][nStudiedTrackPtBins] = histograms->GetHistogramMaxSignalParticlePtInJetCone(iCentrality,iJetPt);
      integralValues[iCentrality][iJetPt][nStudiedTrackPtBins][nIntegralThresholds] = hMaxParticlePtWithinJetCone[iCentrality][iJetPt][nStudiedTrackPtBins]->Integral();
      for(int iThreshold = 0; iThreshold < nIntegralThresholds; iThreshold++){
        lowBinIndex = hMaxParticlePtWithinJetCone[iCentrality][iJetPt][nStudiedTrackPtBins]->GetXaxis()->FindBin(integralThresholds[iThreshold]+0.01);
        highBinIndex = hMaxParticlePtWithinJetCone[iCentrality][iJetPt][nStudiedTrackPtBins]->GetNbinsX();
        integralValues[iCentrality][iJetPt][nStudiedTrackPtBins][iThreshold] = hMaxParticlePtWithinJetCone[iCentrality][iJetPt][nStudiedTrackPtBins]->Integral(lowBinIndex, highBinIndex);
      }
      
      for(int iTrackPt = 0; iTrackPt < nStudiedTrackPtBins; iTrackPt++){
        
        // Calculate the integrals with the track pT cuts
        hMaxParticlePtWithinJetCone[iCentrality][iJetPt][iTrackPt] = histograms->GetHistogramMaxSignalParticlePtInJetConePtCut(iCentrality,iJetPt,iTrackPt);
        integralValues[iCentrality][iJetPt][iTrackPt][nIntegralThresholds] = hMaxParticlePtWithinJetCone[iCentrality][iJetPt][iTrackPt]->Integral();
        for(int iThreshold = 0; iThreshold < nIntegralThresholds; iThreshold++){
          lowBinIndex = hMaxParticlePtWithinJetCone[iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(integralThresholds[iThreshold]+0.01);
          highBinIndex = hMaxParticlePtWithinJetCone[iCentrality][iJetPt][iTrackPt]->GetNbinsX();
          integralValues[iCentrality][iJetPt][iTrackPt][iThreshold] = hMaxParticlePtWithinJetCone[iCentrality][iJetPt][iTrackPt]->Integral(lowBinIndex, highBinIndex);
        }
        
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
  
  // ================================================ //
  // Print to console information about the integrals //
  // ================================================ //
  
  // First, print the intagrals without requirement on background particle pT
  for(int iCentrality = firstIntegratedCentralityBin; iCentrality <= lastIntegratedCentralityBin; iCentrality++){
    for(int iJetPt = firstIntegratedJetPtBinEEC; iJetPt <= lastIntegratedJetPtBinEEC; iJetPt++){
      cout << "Centrality " << histograms->GetCentralityBinBorder(iCentrality) << "-" << histograms->GetCentralityBinBorder(iCentrality+1)  << "% ";
      cout << histograms->GetJetPtBinBorderEEC(iJetPt) << " < jet pT < " << histograms->GetJetPtBinBorderEEC(iJetPt+1) << " GeV" << endl;
      for(int iThreshold = 0; iThreshold < nIntegralThresholds; iThreshold++){
        cout << "Above " << integralThresholds[iThreshold] << " / Total: " << integralValues[iCentrality][iJetPt][nStudiedTrackPtBins][iThreshold] / integralValues[iCentrality][iJetPt][nStudiedTrackPtBins][nIntegralThresholds] << "  ";
      }
      cout << endl;
      
      // Loop over different maximum particle pT requirements for background particles
      for(int iTrackPt = 0; iTrackPt < nStudiedTrackPtBins; iTrackPt++){
        cout << "Background particle pT > " << histograms->GetMaxTrackPtWithinJetConeBinBorder(iTrackPt) << " GeV" << endl;
        for(int iThreshold = 0; iThreshold < nIntegralThresholds; iThreshold++){
          cout << "Above " << integralThresholds[iThreshold] << " / Total: " << integralValues[iCentrality][iJetPt][iTrackPt][iThreshold] / integralValues[iCentrality][iJetPt][iTrackPt][nIntegralThresholds] << "  ";
        }
        cout << endl;
      }
      
    } // Jet pT loop
  } // Centrality loop
  
}
