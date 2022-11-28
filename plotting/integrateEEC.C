#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"

/*
 * Macro for calculating integrals of the EEC distributions.
 * Main idea is to check signal/background in large DeltaR in MC.
 */
void integrateEEC(){

  // File from which the integrals are calculated
  TString inputFileName = "data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJet_updatedMultiplicityAndDensity_wtaAxis_noTrigger_preprocessed_2022-10-17.root";
  
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
  const int nTrackPtBinsEEC = card->GetNTrackPtBinsEEC();
  
  // Default binning ranges for reference
  // centrality = {0,10,30,50,90}
  // track pT =  {0.7,1,2,3,4,8,12,300}
  // jet pT in energy-energy correlators = {120,140,160,180,200,300,500,5030}
  // track pT in energy-energy correlators = {0.7,1,2,3,4,6,300}
  
  // Bin range to be integrated
  int firstIntegratedCentralityBin = 0;
  int lastIntegratedCentralityBin = 0;
  
  int firstIntegratedJetPtBinEEC = 0;
  int lastIntegratedJetPtBinEEC = nJetPtBinsEEC-1; // Note: Jets integrated over all pT ranges are in nJetPtBinsEEC bin
  
  int firstIntegratedTrackPtBinEEC = 7;
  int lastIntegratedTrackPtBinEEC = 7;
  
  // Select the types of energy-energy correlators to integrate
  bool integrateEnergyEnergyCorrelator[EECHistogramManager::knEnergyEnergyCorrelatorTypes];
  integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator] = true;
  integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPt] = false;
  integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorUncorrected] = false;
  integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPtUncorrected] = false;
  
  int integratedEnergyEnergyCorrelatorIndex = 0;
  
  // Create and setup a new histogram manager to project and handle the histograms
  EECHistogramManager *histograms = new EECHistogramManager(inputFile,card);
  
  // Choose the energy-energy correlator types to load
  histograms->SetLoadEnergyEnergyCorrelators(integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator]);
  histograms->SetLoadEnergyEnergyCorrelatorsJetPt(integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPt]);
  histograms->SetLoadEnergyEnergyCorrelatorsUncorrected(integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorUncorrected]);
  histograms->SetLoadEnergyEnergyCorrelatorsJetPtUncorrected(integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPtUncorrected]);
  
  // Choose the bin ranges
  histograms->SetCentralityBinRange(firstIntegratedCentralityBin,lastIntegratedCentralityBin);
  histograms->SetJetPtBinRangeEEC(firstIntegratedJetPtBinEEC,lastIntegratedJetPtBinEEC);
  histograms->SetTrackPtBinRangeEEC(firstIntegratedTrackPtBinEEC,lastIntegratedTrackPtBinEEC);
  
  // Load the histograms from the file
  histograms->LoadProcessedHistograms();
  
  // Initialize the energy-energy correlator histogram array to NULL
  TH1D* hEnergyEnergyCorrelator[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC][EECHistograms::knSubeventCombinations+1];
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC+1; iJetPt++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
          for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventCombinations+1; iSubevent++){
            hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSubevent] = NULL;
          } // Subevent loop
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Energy-energy correlator loop
  
  // Get the histograms from the histogram manager and calculate integrals
  double normalizationFactor;
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only read the selected energy-energy correlator types
    if(!integrateEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
    integratedEnergyEnergyCorrelatorIndex = iEnergyEnergyCorrelator;  // Remember at least one energy-energy correlator type that is loaded
    
    for(int iCentrality = firstIntegratedCentralityBin; iCentrality <= lastIntegratedCentralityBin; iCentrality++){
      for(int iJetPt = firstIntegratedJetPtBinEEC; iJetPt <= lastIntegratedJetPtBinEEC; iJetPt++){
        for(int iTrackPt = firstIntegratedTrackPtBinEEC; iTrackPt <= lastIntegratedTrackPtBinEEC; iTrackPt++){
          
          // For normalization, read the histogram with all pair combinations
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][EECHistograms::knSubeventCombinations] = histograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, 0, EECHistograms::knSubeventCombinations);
          normalizationFactor = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][EECHistograms::knSubeventCombinations]->Integral("width");
          hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][EECHistograms::knSubeventCombinations]->Scale(1/normalizationFactor);
          
          for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventCombinations; iSubevent++){
            hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSubevent] = histograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, 0, iSubevent);
            hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSubevent]->Scale(1/normalizationFactor);
            
          } // Subevent loop
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Energy-energy correlator loop
  
  // Find the binning information for DeltaR bins from a histogram
  const int nDeltaRBins = hEnergyEnergyCorrelator[integratedEnergyEnergyCorrelatorIndex][firstIntegratedCentralityBin][firstIntegratedJetPtBinEEC][firstIntegratedTrackPtBinEEC][0]->GetNbinsX();
  double deltaRBinBorders[nDeltaRBins+1];
  deltaRBinBorders[0] = hEnergyEnergyCorrelator[integratedEnergyEnergyCorrelatorIndex][firstIntegratedCentralityBin][firstIntegratedJetPtBinEEC][firstIntegratedTrackPtBinEEC][0]->GetXaxis()->GetBinLowEdge(1);
  for(int iBin = 1; iBin <= nDeltaRBins; iBin++){
    deltaRBinBorders[iBin] = hEnergyEnergyCorrelator[integratedEnergyEnergyCorrelatorIndex][firstIntegratedCentralityBin][firstIntegratedJetPtBinEEC][firstIntegratedTrackPtBinEEC][0]->GetXaxis()->GetBinUpEdge(iBin);
  }
  
  // Find the bin index for DeltaR = 0.41. We want to only integrate the high deltaR region and study the signal/background ratio there
  const int firstHighDeltaRBin = hEnergyEnergyCorrelator[integratedEnergyEnergyCorrelatorIndex][firstIntegratedCentralityBin][firstIntegratedJetPtBinEEC][firstIntegratedTrackPtBinEEC][0]->GetXaxis()->FindBin(0.41);
  
  // Do the integrals in the high pT region and print the results to console
  double signalIntegral, backgroundIntegral;
  const char* eecTitle[] = {"Energy-energy correlator", "Energy-energy correlator/jet $p_{\\mathrm{T}}$", "Energy-energy correlator UC", "Energy-energy correlator/jet $p_{\\mathrm{T}}$ UC"};
  
  
  // Get the histograms from the histogram manager and calculate integrals
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only read the selected energy-energy correlator types
    if(!integrateEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
    integratedEnergyEnergyCorrelatorIndex = iEnergyEnergyCorrelator;  // Remember at least one energy-energy correlator type that is loaded
    
    for(int iCentrality = firstIntegratedCentralityBin; iCentrality <= lastIntegratedCentralityBin; iCentrality++){
      for(int iJetPt = firstIntegratedJetPtBinEEC; iJetPt <= lastIntegratedJetPtBinEEC; iJetPt++){
        for(int iTrackPt = firstIntegratedTrackPtBinEEC; iTrackPt <= lastIntegratedTrackPtBinEEC; iTrackPt++){
          
          // Print the bin contents in a form of a beamer slide
          cout << endl;
          cout << "\\begin{frame}" << endl;
          cout << Form("\\frametitle{%s bin contents}", eecTitle[iEnergyEnergyCorrelator]) << endl;
          cout << "\\begin{itemize}" << endl;
          cout << "  \\item Signal: Pythia+Pythia, Background: Pythia+Hydjet, Hydjet+Hydjet" << endl;
          cout << "  \\item Centrality: " << histograms->GetCentralityBinBorder(iCentrality) << "-" << histograms->GetCentralityBinBorder(iCentrality+1) << "\\%" <<  endl;
          if(iJetPt == nJetPtBinsEEC){
            cout << "  \\item Track $p_{\\mathrm{T}} > \\qty{" << histograms->GetTrackPtBinBorderEEC(iTrackPt) << "}{GeV}$, Jet $p_{\\mathrm{T}}$ > \\qty{" << histograms->GetCard()->GetJetPtCut() << "}{GeV}$" << endl;
          } else {
            cout << "  \\item Track $p_{\\mathrm{T}} > \\qty{" << histograms->GetTrackPtBinBorderEEC(iTrackPt) << "}{GeV}$, $" << histograms->GetJetPtBinBorderEEC(iJetPt) << " < \\mathrm{Jet}\\;p_{\\mathrm{T}} < \\qty{" << histograms->GetJetPtBinBorderEEC(iJetPt+1) << "}{GeV}$" << endl;
          }
          cout << "\\end{itemize}" << endl;
          cout << "\\begin{center}" << endl;
          cout << "  \\begin{tabular}{cccc}" << endl;
          cout << "    \\toprule" << endl;
          cout << "    Bin & Signal & Background & Background/Signal \\\\" << endl;
          cout << "    \\midrule" << endl;
          
          // Set the correct precision for printing floating point numbers
          // cout << fixed << setprecision(4);
          
          for(int iBin = firstHighDeltaRBin; iBin <= nDeltaRBins; iBin++){
            
            // Calculate the integrals for signal and background
            signalIntegral = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0]->Integral(iBin, iBin, "width");
            backgroundIntegral = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][1]->Integral(iBin, iBin, "width");
            backgroundIntegral += hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][2]->Integral(iBin, iBin, "width");
            
            cout << "    $" << Form("%.3f",deltaRBinBorders[iBin-1]) << " < \\Delta r < " << Form("%.3f",deltaRBinBorders[iBin]) << "$ & $";
            cout << Form("%.7f",signalIntegral) << "$ & $";
            cout << Form("%.7f",backgroundIntegral) << "$ & $";
            cout << Form("%.2f",backgroundIntegral/signalIntegral) << "$ \\\\" << endl;
          }
          
          cout << "    \\bottomrule" << endl;
          cout << "  \\end{tabular}" << endl;
          cout << "\\end{center}" << endl;
          cout << "\\end{frame}" << endl;
          cout << endl;
          
          // Print the integral results in a form of a beamer slide
          cout << endl;
          cout << "\\begin{frame}" << endl;
          cout << Form("\\frametitle{%s integrals}", eecTitle[iEnergyEnergyCorrelator]) << endl;
          cout << "\\begin{itemize}" << endl;
          cout << "  \\item Signal: Pythia+Pythia, Background: Pythia+Hydjet, Hydjet+Hydjet" << endl;
          cout << "  \\item Centrality: " << histograms->GetCentralityBinBorder(iCentrality) << "-" << histograms->GetCentralityBinBorder(iCentrality+1) << "\\%" <<  endl;
          if(iJetPt == nJetPtBinsEEC){
            cout << "  \\item Track $p_{\\mathrm{T}} > \\qty{" << histograms->GetTrackPtBinBorderEEC(iTrackPt) << "}{GeV}$, Jet $p_{\\mathrm{T}}$ > \\qty{" << histograms->GetCard()->GetJetPtCut() << "}{GeV}$" << endl;
          } else {
            cout << "  \\item Track $p_{\\mathrm{T}} > \\qty{" << histograms->GetTrackPtBinBorderEEC(iTrackPt) << "}{GeV}$, $" << histograms->GetJetPtBinBorderEEC(iJetPt) << " < \\mathrm{Jet}\\;p_{\\mathrm{T}} < \\qty{" << histograms->GetJetPtBinBorderEEC(iJetPt+1) << "}{GeV}$" << endl;
          }
          cout << "\\end{itemize}" << endl;
          cout << "\\begin{center}" << endl;
          cout << "  \\begin{tabular}{cccc}" << endl;
          cout << "    \\toprule" << endl;
          cout << "    Region & Signal & Background & Background/Signal \\\\" << endl;
          cout << "    \\midrule" << endl;
          
          // Set the correct precision for printing floating point numbers
          // cout << fixed << setprecision(4);
          
          for(int iBin = firstHighDeltaRBin; iBin <= nDeltaRBins; iBin++){
            
            // Calculate the integrals for signal and background
            signalIntegral = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0]->Integral(iBin, nDeltaRBins, "width");
            backgroundIntegral = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][1]->Integral(iBin, nDeltaRBins, "width");
            backgroundIntegral += hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][2]->Integral(iBin, nDeltaRBins, "width");
            
            cout << "    $" << Form("%.3f",deltaRBinBorders[iBin-1]) << " < \\Delta r < " << deltaRBinBorders[nDeltaRBins] << "$ & $";
            cout << Form("%.7f",signalIntegral) << "$ & $";
            cout << Form("%.7f",backgroundIntegral) << "$ & $";
            cout << Form("%.2f",backgroundIntegral/signalIntegral) << "$ \\\\" << endl;
          }
          
          cout << "    \\bottomrule" << endl;
          cout << "  \\end{tabular}" << endl;
          cout << "\\end{center}" << endl;
          cout << "\\end{frame}" << endl;
          cout << endl;
          
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Energy-energy correlator loop
}
