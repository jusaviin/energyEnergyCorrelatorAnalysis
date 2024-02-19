#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Macro for finding energy-energy correlator signal to background ratio from data
 */
void findMeanPtForUnfoldedEnergyEnergyCorrelators(){
  
  // Input file
  TString inputFileName = "data/eecAnalysis_akFlowJet_nominalEnergyWeight_optimizedUnfoldingBins_manyUnfoldingBinsForScaling_processed_2024-01-17.root";
  // eecAnalysis_akFlowJet_energyWeightSquared_optimizedUnfoldingBins_fixedCovarianceMatrix_unfoldingWithCovariance_processed_2024-01-23.root
  // eecAnalysis_akFlowJet_nominalEnergyWeight_optimizedUnfoldingBins_fixedCovarianceMatrix_unfoldingWithCovariance_processed_2024-01-23.root
  // eecAnalysis_akFlowJet_energyWeightSquared_optimizedUnfoldingBins_manyUnfoldingBinsForScaling_processed_2024-01-17.root
  // eecAnalysis_akFlowJet_nominalEnergyWeight_optimizedUnfoldingBins_manyUnfoldingBinsForScaling_processed_2024-01-17.root

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
  const int weightExponent = card->GetWeightExponent();
  
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
  bool individualCentrality = true; // True = make different figure for each bin. False = plot all centrality bin to the same figure.

  std::vector<std::pair<double,double>> comparedJetPtBin;
  comparedJetPtBin.push_back(std::make_pair(100,120));
  comparedJetPtBin.push_back(std::make_pair(120,140));
  comparedJetPtBin.push_back(std::make_pair(140,160));
  comparedJetPtBin.push_back(std::make_pair(160,180));
  comparedJetPtBin.push_back(std::make_pair(180,200));
  comparedJetPtBin.push_back(std::make_pair(200,220));
  comparedJetPtBin.push_back(std::make_pair(220,240));
  bool individualJetPt = true; // True = make different figure for each bin. False = plot all jet pT bin to the same figure.

  std::vector<double> comparedTrackPtBin;
  if(weightExponent > 1){
    comparedTrackPtBin.push_back(1.0);
    comparedTrackPtBin.push_back(1.5);
  }
  comparedTrackPtBin.push_back(2.0);
  comparedTrackPtBin.push_back(2.5);
  comparedTrackPtBin.push_back(3.0);
  bool individualTrackPt = true; // True = make different figure for each bin. False = plot all track pT bin to the same figure.

  // Drawing options
  const bool drawRawEnergyEnergyCorrelatorFits = false;
  const bool drawUnfoldedEnergyEnergyCorrelatorFits = false;

  // If we are dealing with MC, add 4% centrality shift to centrality bins
  if(card->GetDataType().Contains("MC")){
    for(auto& centralityBin : comparedCentralityBin){
      centralityBin.first += 4;
      centralityBin.second += 4;
    }
  }

  // Only allow one variable for which all bins are plotted to the same figure
  if(individualCentrality + individualJetPt + individualTrackPt < 2){
    cout << "You are tring to plot too many bins to the same figure!" << endl;
    cout << "This macro can only plot all the bins from one variable." << endl;
    cout << "Please check your configuration!" << endl; 
    return;
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
  TH1D* hUnfoldedEnergyEnergyCorrelator[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    hJetPt[iCentrality] = NULL;
    for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt] = NULL;
        hUnfoldedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt] = NULL;
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
  int iJetUnfoldedPt = 0;
  std::vector<bool> skipUnfolded;
  bool firstLoop = true;
  
  // Axes for mean jet pT vs. signal/background ratio plots 
  const int nAnalyzedJetPtBins = comparedJetPtBin.size();
  double jetMeanPtAxis[nCentralityBins][nTrackPtBinsEEC][nAnalyzedJetPtBins];
  double jetMeanPtAxisError[nCentralityBins][nTrackPtBinsEEC][nAnalyzedJetPtBins];
  double deltaRpeakRaw[nCentralityBins][nTrackPtBinsEEC][nAnalyzedJetPtBins];
  double deltaRpeakRawError[nCentralityBins][nTrackPtBinsEEC][nAnalyzedJetPtBins];

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

        // Load raw and unfolded energy-energy correlator histograms and normalize them to make scales more reasonable
        hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt] = histograms->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt);
        hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Scale(1.0/hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Integral("width"));
        hUnfoldedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt] = histograms->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistogramManager::kEnergyEnergyCorrelatorUnfolded);

        // Skip unfolded correlators in jet pT bins in which the unfolding is not performed
        if(firstLoop){
          skipUnfolded.push_back(hUnfoldedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt] == NULL);
        }

        if(!skipUnfolded.back()){
          hUnfoldedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Scale(1.0/hUnfoldedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Integral("width"));
        }

        // Calculate the mean pT values in the studied bins and put them to an array
        lowIntegralBin = hJetPt[iCentrality]->FindBin(jetPtBin.first + epsilon);
        highIntegralBin = hJetPt[iCentrality]->FindBin(jetPtBin.second - epsilon);
        hJetPt[iCentrality]->GetXaxis()->SetRange(lowIntegralBin, highIntegralBin);
        jetMeanPtAxis[iCentrality][iTrackPt][iJetMeanPt] = hJetPt[iCentrality]->GetMean();
        jetMeanPtAxisError[iCentrality][iTrackPt][iJetMeanPt] = hJetPt[iCentrality]->GetMeanError();
        iJetMeanPt++;

      } // Jet pT loop
      firstLoop = false;
    } // Track pT loop
  } // Centrality loop

  // Make the unfolded peak position arrays taking into account that some bins are skipped
  int totalUnfoldedBins = 0;
  for(auto isSkipped : skipUnfolded){
    totalUnfoldedBins = totalUnfoldedBins + 1 - isSkipped;
  }
  const int nUnfoldedJetPtBins = totalUnfoldedBins;
  double deltaRpeakUnfolded[nCentralityBins][nTrackPtBinsEEC][nUnfoldedJetPtBins];
  double deltaRpeakUnfoldedError[nCentralityBins][nTrackPtBinsEEC][nUnfoldedJetPtBins];
  double jetMeanPtUnfoldedAxis[nCentralityBins][nTrackPtBinsEEC][nUnfoldedJetPtBins];
  double jetMeanPtUnfoldedAxisError[nCentralityBins][nTrackPtBinsEEC][nUnfoldedJetPtBins];

  // Copy the jet mean pT values from regular arrays to unfolded arrays
  for(auto centralityBin : comparedCentralityBin){
    iCentrality = card->FindBinIndexCentrality(centralityBin);
    for(auto trackPtBin : comparedTrackPtBin){
      iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);
      iJetMeanPt = 0;
      iJetUnfoldedPt = 0;
      for(auto jetPtBin : comparedJetPtBin){

        if(skipUnfolded.at(iJetMeanPt)){
          iJetMeanPt++;
        } else {
          jetMeanPtUnfoldedAxis[iCentrality][iTrackPt][iJetUnfoldedPt] = jetMeanPtAxis[iCentrality][iTrackPt][iJetMeanPt];
          jetMeanPtUnfoldedAxisError[iCentrality][iTrackPt][iJetUnfoldedPt] = jetMeanPtAxisError[iCentrality][iTrackPt][iJetMeanPt];
          iJetMeanPt++;
          iJetUnfoldedPt++;
        }

      } // Jet pT loop
    } // Track pT loop
  } // Centrality loop


  // =================================================================== //
  // Fit the energy-energy correlator distributions and see what happens //
  // =================================================================== //

  // Create fit functions for signal to background ratio
  TF1* fitForPeak[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TF1* fitForUnfoldedPeak[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];

  // Initialize the fit functions
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){

        fitForPeak[iCentrality][iJetPt][iTrackPt] =  new TF1(Form("peakFitFunction%d%d%d", iCentrality, iJetPt, iTrackPt), "pol3", 0, 0.8);
        fitForUnfoldedPeak[iCentrality][iJetPt][iTrackPt] =  new TF1(Form("unfoldedPeakFitFunction%d%d%d", iCentrality, iJetPt, iTrackPt), "pol3", 0, 0.8);

      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

  // Define jet pT dependent fitting range for the distributions
  std::vector<std::pair<std::pair<double,double>,std::pair<double,double>>> fitRangeForJetPt;
  fitRangeForJetPt.push_back(std::make_pair(std::make_pair(100,120),std::make_pair(0.02,0.05)));
  fitRangeForJetPt.push_back(std::make_pair(std::make_pair(120,140),std::make_pair(0.01,0.04)));
  fitRangeForJetPt.push_back(std::make_pair(std::make_pair(140,160),std::make_pair(0.01,0.04)));
  fitRangeForJetPt.push_back(std::make_pair(std::make_pair(160,180),std::make_pair(0.008,0.03)));
  fitRangeForJetPt.push_back(std::make_pair(std::make_pair(180,200),std::make_pair(0.008,0.028)));
  fitRangeForJetPt.push_back(std::make_pair(std::make_pair(200,220),std::make_pair(0.008,0.028)));
  fitRangeForJetPt.push_back(std::make_pair(std::make_pair(220,240),std::make_pair(0.006,0.025)));

  if(weightExponent == 2){
    fitRangeForJetPt.clear();
    fitRangeForJetPt.push_back(std::make_pair(std::make_pair(100,120),std::make_pair(0.008,0.03)));
    fitRangeForJetPt.push_back(std::make_pair(std::make_pair(120,140),std::make_pair(0.008,0.028)));
    fitRangeForJetPt.push_back(std::make_pair(std::make_pair(140,160),std::make_pair(0.006,0.025)));
    fitRangeForJetPt.push_back(std::make_pair(std::make_pair(160,180),std::make_pair(0.005,0.02)));
    fitRangeForJetPt.push_back(std::make_pair(std::make_pair(180,200),std::make_pair(0.005,0.02)));
    fitRangeForJetPt.push_back(std::make_pair(std::make_pair(200,220),std::make_pair(0.004,0.015)));
    fitRangeForJetPt.push_back(std::make_pair(std::make_pair(220,240),std::make_pair(0.004,0.015)));
  }

  // Make a drawer to check that the fits are reasonable
  JDrawer* drawer = new JDrawer();
  TLegend* legend;
  TString centralityString, compactCentralityString;
  TString trackPtString, compactTrackPtString;
  TString jetPtString, compactJetPtString;

  // Loop over all defined bins and make the fits
  std::pair<double,double> thisFitRange;
  bool fitRangeFound;
  for(auto centralityBin : comparedCentralityBin){
    iCentrality = card->FindBinIndexCentrality(centralityBin);
    centralityString = Form("Cent: %.0f-%.0f%%", centralityBin.first, centralityBin.second);
    compactCentralityString = Form("_C=%.0f-%.0f", centralityBin.first, centralityBin.second);
    for(auto trackPtBin : comparedTrackPtBin){
      iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);
      trackPtString = Form("%.1f < track p_{T}", trackPtBin);
      compactTrackPtString = Form("_T>%.1f", trackPtBin);
      compactTrackPtString.ReplaceAll(".", "v");
      iJetMeanPt = 0;
      iJetUnfoldedPt = 0;
      for(auto jetPtBin : comparedJetPtBin){
        iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);
        jetPtString = Form("%.0f < jet p_{T} < %.0f", jetPtBin.first, jetPtBin.second);
        compactJetPtString = Form("_J=%.0f-%.0f", jetPtBin.first, jetPtBin.second);

        // Find the fit range for the defined jet pT bin
        fitRangeFound = false;
        for(auto fitRangeForThisJetPt : fitRangeForJetPt){
          if(fitRangeForThisJetPt.first == jetPtBin){
            thisFitRange = fitRangeForThisJetPt.second;
            fitRangeFound = true;
            break;
          }
        }

        // If we did not find fit range, crash the code
        if(!fitRangeFound){
          cout << "ERROR! No fit range defined for " << jetPtBin.first << " < jet pT < " << jetPtBin.second << endl;
          cout << "Please check your code! Cannot do fitting without a proper range!" << endl;
          return;
        }

        // Once the fitting range is found, do the fitting with the fit function
        hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Fit(fitForPeak[iCentrality][iJetPt][iTrackPt], "0", "", thisFitRange.first, thisFitRange.second);

        // Draw the fits if specified
        if(drawRawEnergyEnergyCorrelatorFits){

          // After fitting is done, draw the distribution an fit to canvas to see how the fit looks like
          hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.005,0.08);
          drawer->DrawHistogram(hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt], "#Deltar", "EEC", " ");
          fitForPeak[iCentrality][iJetPt][iTrackPt]->SetRange(thisFitRange.first, thisFitRange.second);
          fitForPeak[iCentrality][iJetPt][iTrackPt]->Draw("same");

          legend = new TLegend(0.5, 0.25, 0.8, 0.55);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);
          legend->AddEntry((TObject*)0, "Raw correlator", "");
          legend->AddEntry((TObject*)0, centralityString.Data(), "");
          legend->AddEntry((TObject*)0, jetPtString.Data(), "");
          legend->AddEntry((TObject*)0, trackPtString.Data(), "");
          legend->Draw();

        } // Drawing fits to energy-energy correlators

        // Find the maximum from the fit function and save the value to the fit function array
        deltaRpeakRaw[iCentrality][iTrackPt][iJetMeanPt] = fitForPeak[iCentrality][iJetPt][iTrackPt]->GetMaximumX(thisFitRange.first, thisFitRange.second);
        deltaRpeakRawError[iCentrality][iTrackPt][iJetMeanPt] = 0.001; // Maybe do some more fancy error determination here?

        // There are some bins which are not unfolded. We cannot do fitting in those bins
        if(!skipUnfolded.at(iJetMeanPt)){

          // Do the fitting also for the unfolded energy-energy correlators
          hUnfoldedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Fit(fitForUnfoldedPeak[iCentrality][iJetPt][iTrackPt], "0", "", thisFitRange.first, thisFitRange.second);

          // Draw the unfolded fits if specified
          if(drawUnfoldedEnergyEnergyCorrelatorFits){

            // After fitting is done, draw the distribution an fit to canvas to see how the fit looks like
            hUnfoldedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.005,0.08);
            drawer->DrawHistogram(hUnfoldedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt], "#Deltar", "EEC", " ");
            fitForUnfoldedPeak[iCentrality][iJetPt][iTrackPt]->SetRange(thisFitRange.first, thisFitRange.second);
            fitForUnfoldedPeak[iCentrality][iJetPt][iTrackPt]->Draw("same");

            legend = new TLegend(0.5, 0.25, 0.8, 0.55);
            legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);
            legend->AddEntry((TObject*)0, "Unfolded correlator", "");
            legend->AddEntry((TObject*)0, centralityString.Data(), "");
            legend->AddEntry((TObject*)0, jetPtString.Data(), "");
            legend->AddEntry((TObject*)0, trackPtString.Data(), "");
            legend->Draw();

          } // Drawing fits to unfolded energy-energy correlators

          // Find the maximum from the fit function and save the value to the fit function array
          deltaRpeakUnfolded[iCentrality][iTrackPt][iJetUnfoldedPt] = fitForUnfoldedPeak[iCentrality][iJetPt][iTrackPt]->GetMaximumX(thisFitRange.first, thisFitRange.second);
          deltaRpeakUnfoldedError[iCentrality][iTrackPt][iJetUnfoldedPt] = 0.001; // Maybe do some more fancy error determination here?

          iJetUnfoldedPt++;

        }

        // Increment the graph jet pT index
        iJetMeanPt++;

      } // Jet pT loop
    } // Track pT loop
  } // Centrality loop
  
  // ========================================================================= //
  // Put the determined peak positions to graphs and find the trend in mean pT //
  // ========================================================================= //

  // Fill the obtained information into TGraphErrors
  TGraphErrors* gDeltaRPeakRaw[nCentralityBins][nTrackPtBinsEEC];
  TGraphErrors* gDeltaRPeakUnfolded[nCentralityBins][nTrackPtBinsEEC];

  // Initialize the signal to background graphs to NULL
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
      gDeltaRPeakRaw[iCentrality][iTrackPt] = NULL;
      gDeltaRPeakUnfolded[iCentrality][iTrackPt] = NULL;
    } // Track pT loop
  } // Centrality loop

  // Make new graphs for the selected bins
  for(auto centralityBin : comparedCentralityBin){
    iCentrality = card->FindBinIndexCentrality(centralityBin);
    for(auto trackPtBin : comparedTrackPtBin){
      iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);

      gDeltaRPeakRaw[iCentrality][iTrackPt] = new TGraphErrors(nAnalyzedJetPtBins, jetMeanPtAxis[iCentrality][iTrackPt], deltaRpeakRaw[iCentrality][iTrackPt], jetMeanPtAxisError[iCentrality][iTrackPt], deltaRpeakRawError[iCentrality][iTrackPt]);
      gDeltaRPeakUnfolded[iCentrality][iTrackPt] = new TGraphErrors(nUnfoldedJetPtBins, jetMeanPtUnfoldedAxis[iCentrality][iTrackPt], deltaRpeakUnfolded[iCentrality][iTrackPt], jetMeanPtUnfoldedAxisError[iCentrality][iTrackPt], deltaRpeakUnfoldedError[iCentrality][iTrackPt]);

    } // Track pT loop
  } // Centrality loop

  // Create functions to fit the peak position graphs
  TF1* fitToRawDeltaRPeak[nCentralityBins][nTrackPtBinsEEC];
  TF1* fitToUnfoldedDeltaRPeak[nCentralityBins][nTrackPtBinsEEC];

  // Initialize the fit functions
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        fitToRawDeltaRPeak[iCentrality][iTrackPt] = new TF1(Form("deltaRpeakFit%d%d", iCentrality, iTrackPt), "pol2", comparedJetPtBin.at(0).first, comparedJetPtBin.back().second);
        fitToRawDeltaRPeak[iCentrality][iTrackPt]->SetParameter(0,0);
        fitToRawDeltaRPeak[iCentrality][iTrackPt]->SetParameter(1,0.1);
        fitToRawDeltaRPeak[iCentrality][iTrackPt]->SetParameter(2,0);

        fitToUnfoldedDeltaRPeak[iCentrality][iTrackPt] = new TF1(Form("unfoldedDeltaRpeakFit%d%d", iCentrality, iTrackPt), "pol2", comparedJetPtBin.at(0).first, comparedJetPtBin.back().second);
        fitToUnfoldedDeltaRPeak[iCentrality][iTrackPt]->SetParameter(0,0);
        fitToUnfoldedDeltaRPeak[iCentrality][iTrackPt]->SetParameter(1,0.1);
        fitToUnfoldedDeltaRPeak[iCentrality][iTrackPt]->SetParameter(2,0);
    }
  }

  // Set drawing style for graphs
  drawer->SetDefaultAppearanceGraph();

  // Helper variables
  double minimumValue, maximumValue, drawMargin;

  // Draw each graph to separate canvas
  for(auto centralityBin : comparedCentralityBin){
    iCentrality = card->FindBinIndexCentrality(centralityBin);
    centralityString = Form("Cent: %.0f-%.0f%%", centralityBin.first, centralityBin.second);
    compactCentralityString = Form("_C=%.0f-%.0f", centralityBin.first, centralityBin.second);
    for(auto trackPtBin : comparedTrackPtBin){
      iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);
      trackPtString = Form("%.1f < track p_{T}", trackPtBin);
      compactTrackPtString = Form("_T>%.1f", trackPtBin);
      compactTrackPtString.ReplaceAll(".", "v");

      // Setup the legend for plots
      legend = new TLegend(0.5, 0.7, 0.9, 0.85);
      legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);
      legend->AddEntry((TObject*)0, centralityString.Data(), "");
      legend->AddEntry((TObject*)0, trackPtString.Data(), "");

      // Find a good y-axis scale for the plots
      maximumValue = gDeltaRPeakRaw[iCentrality][iTrackPt]->GetPointY(0);
      minimumValue = gDeltaRPeakRaw[iCentrality][iTrackPt]->GetPointY(nAnalyzedJetPtBins-1);
      drawMargin = (maximumValue - minimumValue) * 0.2;

      // Set a fancy style for the graph
      gDeltaRPeakRaw[iCentrality][iTrackPt]->SetMarkerStyle(kFullCircle);
      gDeltaRPeakRaw[iCentrality][iTrackPt]->SetMarkerColor(kBlue);
      gDeltaRPeakRaw[iCentrality][iTrackPt]->SetLineColor(kBlue);

      // Fit a line to the raw peak position
      gDeltaRPeakRaw[iCentrality][iTrackPt]->Fit(fitToRawDeltaRPeak[iCentrality][iTrackPt], "0");

      // Draw the graphs to canvas
      drawer->DrawGraphCustomAxes(gDeltaRPeakRaw[iCentrality][iTrackPt], comparedJetPtBin.at(0).first, comparedJetPtBin.at(nAnalyzedJetPtBins-1).second, minimumValue - drawMargin, maximumValue + drawMargin, "Jet p_{T}", "Peak #Deltar", " ", "a,p");

      // Add the legend to the canvas
      legend->Draw();

      // Draw the fit function
      fitToRawDeltaRPeak[iCentrality][iTrackPt]->SetLineColor(kBlue);
      fitToRawDeltaRPeak[iCentrality][iTrackPt]->Draw("same");

      // Also fit the unfolded peak position
      gDeltaRPeakUnfolded[iCentrality][iTrackPt]->Fit(fitToUnfoldedDeltaRPeak[iCentrality][iTrackPt], "0");

      // Draw the unfolded values also to the graph
      gDeltaRPeakUnfolded[iCentrality][iTrackPt]->SetMarkerStyle(kFullCross);
      gDeltaRPeakUnfolded[iCentrality][iTrackPt]->SetMarkerColor(kRed);
      gDeltaRPeakUnfolded[iCentrality][iTrackPt]->SetLineColor(kRed);

      gDeltaRPeakUnfolded[iCentrality][iTrackPt]->Draw("same,p");

      // Draw the fit function
      fitToUnfoldedDeltaRPeak[iCentrality][iTrackPt]->SetLineColor(kRed);
      fitToUnfoldedDeltaRPeak[iCentrality][iTrackPt]->Draw("same");
      

    } // Track pT loop
  } // Centrality loop

  // Determine the upshift be the difference of the two fits
  double meanJetPt;
  double targetValue;
  double shiftedPt;
  double averageRelativeShift[nCentralityBins];
  int nBinsForAverage[nCentralityBins];
  for(auto centralityBin : comparedCentralityBin){
    iCentrality = card->FindBinIndexCentrality(centralityBin);
    averageRelativeShift[iCentrality] = 0;
    nBinsForAverage[iCentrality] = 0;
    for(auto trackPtBin : comparedTrackPtBin){
      iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);
      iJetPt = 0;
      for(auto jetPtBin : comparedJetPtBin){
        meanJetPt = jetMeanPtAxis[iCentrality][iTrackPt][iJetPt];

        // Only print the results for analysis region
        if(meanJetPt < 120 || meanJetPt > 200) {
          iJetPt++; 
          continue;
        }

        targetValue = fitToUnfoldedDeltaRPeak[iCentrality][iTrackPt]->Eval(meanJetPt);
        shiftedPt = fitToRawDeltaRPeak[iCentrality][iTrackPt]->GetX(targetValue, comparedJetPtBin.at(0).first, comparedJetPtBin.back().second);

        // Print some information to estimate what is going on
        cout << "Centrality: " << centralityBin.first << "-" << centralityBin.second << "  Track pT > " << trackPtBin << " GeV" << endl;
        cout << "Jet pT " << meanJetPt << " is shifted to " << shiftedPt << " with an upshift of " << shiftedPt-meanJetPt << " which is " <<  ((shiftedPt-meanJetPt)/meanJetPt)*100 << "%% of the original." << endl;

        // Calculate the average relative shift for each centrality bin
        averageRelativeShift[iCentrality] += (shiftedPt-meanJetPt)/meanJetPt;
        nBinsForAverage[iCentrality]++;
   
        iJetPt++;
      } // Jet pT loop
    } // Track pT loop
  } // Centrality loop

  cout << endl;
  for(auto centralityBin : comparedCentralityBin){
    iCentrality = card->FindBinIndexCentrality(centralityBin);

    cout << "Average relative shift for centrality " << centralityBin.first << "-" << centralityBin.second << " is " << averageRelativeShift[iCentrality]/nBinsForAverage[iCentrality] << endl;
  } // Centrality loop  

  /*
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

  */
  

}
