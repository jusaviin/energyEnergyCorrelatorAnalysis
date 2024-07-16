#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"
#include "../src/EECHistograms.h"

/*
 * Macro for finding energy-energy correlator signal to background ratio from data
 */
void signalToBackgroundRatio(){
  
  // Input file
  TString inputFileName = "data/PbPbMC2018_GenGen_eecAnalysis_4pCentShift_cutBadPhi_nominalEnergyWeight_allBackgrounds_matchMultiplicity_someMissing_processed_2024-04-24.root";
  // eecAnalysis_akFlowJet_energyWeightSquared_optimizedUnfoldingBins_fixedCovarianceMatrix_unfoldingWithCovariance_processed_2024-01-23.root
  // eecAnalysis_akFlowJet_nominalEnergyWeight_optimizedUnfoldingBins_fixedCovarianceMatrix_unfoldingWithCovariance_processed_2024-01-23.root
  // PbPbMC2018_RecoReco_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_nominalEnergyWeight_optimizedUnfoldingBins_nominalSmear_processed_2024-01-19.root

  TFile* inputFile = TFile::Open(inputFileName);

  // Check that the input file exists
  if(inputFile == NULL){
    cout << "Error! The file " << inputFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }

  // Read EECCard from the file
  EECCard* card = new EECCard(inputFile);
  
  // ====================================================
  //               Binning configuration
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = card->GetNCentralityBins();
  const int nJetPtBinsEEC = card->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card->GetNTrackPtBinsEEC();
  
  // Select explicitly which bins from the files are compared:
  std::vector<std::pair<double,double>> comparedCentralityBin;
  comparedCentralityBin.push_back(std::make_pair(0,10));
  comparedCentralityBin.push_back(std::make_pair(10,30));
  comparedCentralityBin.push_back(std::make_pair(30,50));
  comparedCentralityBin.push_back(std::make_pair(50,90));

  std::vector<std::pair<double,double>> comparedJetPtBin;
  comparedJetPtBin.push_back(std::make_pair(100,120));
  comparedJetPtBin.push_back(std::make_pair(120,140));
  comparedJetPtBin.push_back(std::make_pair(140,160));
  comparedJetPtBin.push_back(std::make_pair(160,180));
  comparedJetPtBin.push_back(std::make_pair(180,200));
  comparedJetPtBin.push_back(std::make_pair(200,220));

  std::vector<double> comparedTrackPtBin;
  comparedTrackPtBin.push_back(1.0);
  comparedTrackPtBin.push_back(1.5);
  comparedTrackPtBin.push_back(2.0);
  comparedTrackPtBin.push_back(2.5);
  comparedTrackPtBin.push_back(3.0);

  // Fitting parameters
  const bool doFit = true;

  // Figure saving
  const bool saveFigures = false;  // Save figures
  const char* saveComment = "_includeFit";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures

  // If we are dealing with MC, add 4% centrality shift to centrality bins
  if(card->GetDataType().Contains("MC")){
    for(auto& centralityBin : comparedCentralityBin){
      centralityBin.first += 4;
      centralityBin.second += 4;
    }
  }

  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager* histograms;
  histograms = new EECHistogramManager(inputFile, card);

  // Choose the energy-energy correlator types to load
  histograms->SetLoadJetHistograms(true);
  histograms->SetLoadEnergyEnergyCorrelators(true);

  // Choose the bin ranges
  histograms->SetCentralityBinRange(0, card->GetNCentralityBins() - 1);
  histograms->SetJetPtBinRangeEEC(0, card->GetNJetPtBinsEEC() - 1);
  histograms->SetTrackPtBinRangeEEC(0, card->GetNTrackPtBinsEEC() - 1);

  // Load the histograms from the file
  histograms->LoadProcessedHistograms();

  // Energy-energy correlator histograms
  TH1D* hJetPt[nCentralityBins];
  TH1D* hEnergyEnergyCorrelator[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hReflectedCone[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    hJetPt[iCentrality] = NULL;
    for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt] = NULL;
        hReflectedCone[iCentrality][iJetPt][iTrackPt] = NULL;
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
  
  // Helper variables
  double integrationRange = 0.39;
  double epsilon = 0.001;
  int lowIntegralBin, highIntegralBin;
  double signalIntegral, backgroundIntegral;
  double signalIntegralError, backgroundIntegralError;
  int iCentrality;
  int iTrackPt;
  int iJetPt;
  int iJetMeanPt = 0;
  
  // Axes for mean jet pT vs. signal/background ratio plots 
  const int nAnalyzedJetPtBins = comparedJetPtBin.size();
  double jetMeanPtAxis[nCentralityBins][nTrackPtBinsEEC][nAnalyzedJetPtBins];
  double jetMeanPtAxisError[nCentralityBins][nTrackPtBinsEEC][nAnalyzedJetPtBins];
  double signalToBackgroundRatioAxis[nCentralityBins][nTrackPtBinsEEC][nAnalyzedJetPtBins];
  double signalToBackgroundRatioAxisError[nCentralityBins][nTrackPtBinsEEC][nAnalyzedJetPtBins];

  // Get the histograms from the histogram manager and calculate signal to background ratios and mean jet pT:s
  for(auto centralityBin : comparedCentralityBin){
    iCentrality = card->FindBinIndexCentrality(centralityBin);

    // Load the jet pT histogram for this centrality range
    hJetPt[iCentrality] = histograms->GetHistogramJetPt(iCentrality);

    for(auto trackPtBin : comparedTrackPtBin){
      iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);
      iJetMeanPt = 0;
      for(auto jetPtBin : comparedJetPtBin){
        iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);

        // Load the energy-energy correlator and reflected cone histograms
        hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt] = histograms->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt);
        hReflectedCone[iCentrality][iJetPt][iTrackPt] = histograms->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSignalReflectedConePair);

        // Calculate the signal and background integrals
        lowIntegralBin = 1;
        highIntegralBin = hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(integrationRange);
        signalIntegral = hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->IntegralAndError(lowIntegralBin, highIntegralBin, signalIntegralError, "width");
        backgroundIntegral = hReflectedCone[iCentrality][iJetPt][iTrackPt]->IntegralAndError(lowIntegralBin, highIntegralBin, backgroundIntegralError, "width");
       
        // Add the signal to background ratio and it's error to an array
        signalToBackgroundRatioAxis[iCentrality][iTrackPt][iJetMeanPt] = signalIntegral / backgroundIntegral;
        signalToBackgroundRatioAxisError[iCentrality][iTrackPt][iJetMeanPt] = TMath::Sqrt(TMath::Power(signalIntegralError / backgroundIntegral, 2) + TMath::Power(signalIntegral * backgroundIntegralError / backgroundIntegral / backgroundIntegral, 2));

        // Calculate the mean pT values in the studied bins and put them to an array
        lowIntegralBin = hJetPt[iCentrality]->FindBin(jetPtBin.first + epsilon);
        highIntegralBin = hJetPt[iCentrality]->FindBin(jetPtBin.second - epsilon);
        hJetPt[iCentrality]->GetXaxis()->SetRange(lowIntegralBin, highIntegralBin);
        jetMeanPtAxis[iCentrality][iTrackPt][iJetMeanPt] = hJetPt[iCentrality]->GetMean();
        jetMeanPtAxisError[iCentrality][iTrackPt][iJetMeanPt] = hJetPt[iCentrality]->GetMeanError();
        iJetMeanPt++;

      } // Jet pT loop
    } // Track pT loop
  } // Centrality loop

  // Fill the obtained information into TGraphErrors
  TGraphErrors* gSignalToBackgroundRatio[nCentralityBins][nTrackPtBinsEEC];

  // Initialize the signal to background graphs to NULL
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
      gSignalToBackgroundRatio[iCentrality][iTrackPt] = NULL;
    } // Track pT loop
  } // Centrality loop

  // Make new graphs for the selected bins
  for(auto centralityBin : comparedCentralityBin){
    iCentrality = card->FindBinIndexCentrality(centralityBin);
    for(auto trackPtBin : comparedTrackPtBin){
      iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);

      gSignalToBackgroundRatio[iCentrality][iTrackPt] = new TGraphErrors(nAnalyzedJetPtBins, jetMeanPtAxis[iCentrality][iTrackPt], signalToBackgroundRatioAxis[iCentrality][iTrackPt], jetMeanPtAxisError[iCentrality][iTrackPt], signalToBackgroundRatioAxisError[iCentrality][iTrackPt]);

    } // Track pT loop
  } // Centrality loop

  // Create fit functions for signal to background ratio
  TF1* fitToRatio[nCentralityBins][nTrackPtBinsEEC];

  // Initialize the fit functions
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){

      if(iCentrality == 3){
        fitToRatio[iCentrality][iTrackPt] = new TF1(Form("FitFunction%d%d", iCentrality, iTrackPt), "pol1", 120, 200);
        fitToRatio[iCentrality][iTrackPt]->SetParameter(0,0);
        fitToRatio[iCentrality][iTrackPt]->SetParameter(1,0.1);
      } else {
        fitToRatio[iCentrality][iTrackPt] = new TF1(Form("FitFunction%d%d", iCentrality, iTrackPt), "pol2", 120, 200);
        fitToRatio[iCentrality][iTrackPt]->SetParameter(0,0);
        fitToRatio[iCentrality][iTrackPt]->SetParameter(1,0.1);
        fitToRatio[iCentrality][iTrackPt]->SetParameter(2,0);
      }
    }
  }

  // Draw the analyzed bins to a canvas
  JDrawer* drawer = new JDrawer();
  drawer->SetDefaultAppearanceGraph();

  TLegend* legend;
  TString centralityString, compactCentralityString;
  TString trackPtString, compactTrackPtString;

  // Add a name describing the energy weight in the files
  TString energyWeightString, compactEnergyWeightString; 
  if(card->GetWeightExponent() == 2){
    energyWeightString = "Energy weight squared";
    compactEnergyWeightString = "_energyWeightSquared";
  } else {
    energyWeightString = "Nominal energy weight";
    compactEnergyWeightString = "_nominalEnergyWeight";
  }

  double minimumValue, maximumValue, drawMargin;

  // Draw each graph to separate canvas
  for(auto centralityBin : comparedCentralityBin){
    iCentrality = card->FindBinIndexCentrality(centralityBin);
    centralityString = Form("Pythia+Hydjet: %.0f-%.0f%%", centralityBin.first-4, centralityBin.second-4);
    compactCentralityString = Form("_C=%.0f-%.0f", centralityBin.first, centralityBin.second);
    for(auto trackPtBin : comparedTrackPtBin){
      iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);
      trackPtString = Form("p_{T}^{ch} > %.1f GeV", trackPtBin);
      compactTrackPtString = Form("_T>%.1f", trackPtBin);
      compactTrackPtString.ReplaceAll(".", "v");

      // Setup the legend for plots
      legend = new TLegend(0.18, 0.6, 0.58, 0.88);
      legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);
      legend->AddEntry((TObject*)0, energyWeightString.Data(), "");
      legend->AddEntry((TObject*)0, centralityString.Data(), "");
      legend->AddEntry((TObject*)0, trackPtString.Data(), "");

      // Find a good y-axis scale for the plots
      minimumValue = gSignalToBackgroundRatio[iCentrality][iTrackPt]->GetPointY(0);
      maximumValue = gSignalToBackgroundRatio[iCentrality][iTrackPt]->GetPointY(nAnalyzedJetPtBins-1);
      drawMargin = (maximumValue - minimumValue) * 0.2;

      // Set a fancy style for the graph
      gSignalToBackgroundRatio[iCentrality][iTrackPt]->SetMarkerStyle(kFullCircle);
      gSignalToBackgroundRatio[iCentrality][iTrackPt]->SetMarkerColor(kBlue);
      gSignalToBackgroundRatio[iCentrality][iTrackPt]->SetLineColor(kBlue);

      // Fit a line to the signal to background ratio
      if(doFit) gSignalToBackgroundRatio[iCentrality][iTrackPt]->Fit(fitToRatio[iCentrality][iTrackPt]);

      // Draw the graphs to canvas
      drawer->DrawGraphCustomAxes(gSignalToBackgroundRatio[iCentrality][iTrackPt], comparedJetPtBin.at(0).first, comparedJetPtBin.at(nAnalyzedJetPtBins-1).second, minimumValue - drawMargin, maximumValue + drawMargin, "Jet p_{T} (GeV)", "(Signal + BG) / BG", " ", "a,p");

      // Add the legend to the canvas
      legend->Draw();

      // Save the figures to a file
      if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/eecSignaltoBackgroundRatio%s%s%s%s.%s", saveComment, compactEnergyWeightString.Data(), compactCentralityString.Data(), compactTrackPtString.Data(), figureFormat));
      }
      

    } // Track pT loop
  } // Centrality loop

  // Print a table of mean pT:s to be copy-pasted into a corrector class
  cout << endl;
  cout << " double meanJetPt[" << comparedCentralityBin.size() << "][" << comparedJetPtBin.size() << "] = {" << endl;
  iJetPt = 0;
    
  for(auto centralityBin : comparedCentralityBin){
    iCentrality = card->FindBinIndexCentrality(centralityBin);
      
    cout << "{";

    for(auto jetPtBin : comparedJetPtBin){

      cout << jetMeanPtAxis[iCentrality][1][iJetPt];

      // Add the proper syntax
      if(jetPtBin.first != comparedJetPtBin.back().first){
        cout << ",";
      } else {
        cout << "}";
      }

      iJetPt++;
    } // Jet pT loop

    if(centralityBin.first != comparedCentralityBin.back().first){
      cout << ", // " << Form("Centrality %.0f-%.0f%%", centralityBin.first, centralityBin.second) << endl;
    } else {
      if(centralityBin.first != comparedCentralityBin.back().first){
        cout << ", // " << Form("Centrality %.0f-%.0f%%", centralityBin.first, centralityBin.second) << endl;
      } else {
        cout << " // " << Form("Centrality %.0f-%.0f%%", centralityBin.first, centralityBin.second) << endl;
      }
    }
  } // Centrality loop
  cout << "};" << endl;

  // Print a table of the fit parameters to be copy-pasted into a corrector class 
  cout << endl;
  cout << " double signalToBackgroundRatioFitParameters[" << comparedCentralityBin.size() << "][" << comparedTrackPtBin.size() << "][3] = {" << endl;

  for(auto centralityBin : comparedCentralityBin){
    iCentrality = card->FindBinIndexCentrality(centralityBin);
      
    // Set the centrality information for legends and figure saving
    cout << " // " << Form("Centrality %.0f-%.0f%%", centralityBin.first, centralityBin.second) << endl;
    cout << "{";
      
    for(auto trackPtBin : comparedTrackPtBin){
      iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);
        
      cout << "{";
        
      for(int iParameter = 0; iParameter < 3; iParameter++){
          
        // Print the parameter value to the array:
        if(iCentrality == 3 && iParameter == 2){
          cout << 0;
        } else {
          cout << fitToRatio[iCentrality][iTrackPt]->GetParameter(iParameter);
        }
          
        // Add the proper syntax
        if(iParameter < 2){
          cout << ",";
        } else {
          cout << "}";
        }
      } // Parameter pT loop
        
      if(trackPtBin != comparedTrackPtBin.back()){
        cout << ", // " << Form("track pT > %.1f GeV", trackPtBin) << endl;
      } else {
        cout << "}";
        if(trackPtBin != comparedTrackPtBin.back()){
          cout << ", // " << Form("track pT > %.1f GeV", trackPtBin) << endl;
        } else {
          cout << " // " << Form("track pT > %.1f GeV", trackPtBin) << endl;
        }
      }
        
    } // Track pT loop
  } // Centrality loop
  cout << "};" << endl;
  

}
