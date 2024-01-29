#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Macro for comparing energy-energy correlators in custom bins
 */
void compareRelativeUncertainties(){
  
  // Files for comparison
  TString fileName = "data/eecAnalysis_akFlowJet_energyWeightSquared_optimizedUnfoldingBins_unfoldingWithNominalSmear_processed_2024-01-17.root";
  // eecAnalysis_akFlowJet_wtaAxis_energyWeightSquared_covarianceMatrix_unfoldingWithCovarianceMatrix_processed_2023-12-14.root
  // eecAnalysis_akFlowJet_wtaAxis_energyWeightSquared_covarianceMatrix_unfoldingOriginal_processed_2023-12-14.root
  // ppData_pfJets_wtaAxis_energyWeightSquared_covarianceMatrix_jet60or80triggers_unfoldingWithCovarianceMatrix_processed_2023-12-14.root
  
  // Open the file and check that is exists
  TFile* inputFile = TFile::Open(fileName);

  if(inputFile == NULL){
    cout << "Error! The file " << fileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }

  EECCard* card = new EECCard(inputFile);
  
  // Check if we are using PbPb or pp data
  TString collisionSystem = card->GetDataType();
  bool isPbPbData = collisionSystem.Contains("PbPb");
  
  // ====================================================
  //               Binning configuration
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = card->GetNCentralityBins();
  const int nJetPtBinsEEC = card->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card->GetNTrackPtBinsEEC();
  
  // Select explicitly which bins from the files are compared:
  std::vector<std::pair<double,double>> comparedCentralityBin;
  //comparedCentralityBin.push_back(std::make_pair(0,10));
  //comparedCentralityBin.push_back(std::make_pair(10,30));
  //comparedCentralityBin.push_back(std::make_pair(30,50));
  comparedCentralityBin.push_back(std::make_pair(50,90));
  bool individualCentrality = true; // True = make different figure for each bin. False = plot all centrality bin to the same figure.

  if(!isPbPbData){
    comparedCentralityBin.clear();
    comparedCentralityBin.push_back(std::make_pair(0,100));
    individualCentrality = true;
  }

  std::vector<std::pair<double,double>> comparedJetPtBin;
  //comparedJetPtBin.push_back(std::make_pair(120,140));
  //comparedJetPtBin.push_back(std::make_pair(140,160));
  //comparedJetPtBin.push_back(std::make_pair(160,180));
  comparedJetPtBin.push_back(std::make_pair(180,200));
  bool individualJetPt = true; // True = make different figure for each bin. False = plot all jet pT bin to the same figure.

  std::vector<double> comparedTrackPtBin;
  comparedTrackPtBin.push_back(1.0);
  //comparedTrackPtBin.push_back(1.5);
  comparedTrackPtBin.push_back(2.0);
  //comparedTrackPtBin.push_back(2.5);
  comparedTrackPtBin.push_back(3.0);
  bool individualTrackPt = false; // True = make different figure for each bin. False = plot all track pT bin to the same figure.

  // ====================================================
  //                    Configuration
  // ====================================================
  
  // Zooming
  std::pair<double, double> drawingRange = std::make_pair(0.006, 0.39);
  std::pair<double, double> eecZoom = std::make_pair(0,0.5);

  // Figure saving
  const bool saveFigures = true;  // Save figures
  const char* saveComment = "_pythia8";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures

  // Drawing configuration
  std::pair<double, double> ratioZoom = std::make_pair(0.7, 1.3);

  // Only allow one variable for which all bins are plotted to the same figure
  if(individualCentrality + individualJetPt + individualTrackPt < 2){
    cout << "You are tring to plot too many bins to the same figure!" << endl;
    cout << "This macro can only plot all the bins from one variable." << endl;
    cout << "Please check your configuration!" << endl; 
    return;
  }

  // Based on the input from user, make a construct that can be handled later in the code to plot the selected bins
  enum enumTupleDecoder{kCentrality, kJetPt, kTrackPt}; // Components of the n-tuple in binning vector
  std::vector<std::tuple<std::vector<std::pair<double,double>>, std::vector<std::pair<double,double>>, std::vector<double>>> binningInformation;

  // Separate plots for each centrality bin
  std::vector<std::pair<double,double>> temporaryCentralityVector;
  std::vector<std::pair<double,double>> temporaryJetPtVector;
  std::vector<double> temporaryTrackPtVector;

  if(individualCentrality){
    for(auto centralityBin : comparedCentralityBin){
      temporaryCentralityVector.clear();
      temporaryCentralityVector.push_back(centralityBin);

      // Separate plots for each jet pT bin
      if(individualJetPt){

        for(auto jetPtBin : comparedJetPtBin){
          temporaryJetPtVector.clear();
          temporaryJetPtVector.push_back(jetPtBin);

          // Separate track pT bins for each bin
          if(individualTrackPt){

            for(auto trackPtBin : comparedTrackPtBin){
              temporaryTrackPtVector.clear();
              temporaryTrackPtVector.push_back(trackPtBin);

              binningInformation.push_back(std::make_tuple(temporaryCentralityVector, temporaryJetPtVector, temporaryTrackPtVector));

            }

          // All defined track pT bins in a single plot
          } else {

            binningInformation.push_back(std::make_tuple(temporaryCentralityVector, temporaryJetPtVector, comparedTrackPtBin));

          } // Track pT binning if-else

        } // Loop over ket pT bins

      // All defined jet pT bins in a single plot
      } else {

        // Separate track pT bins for each bin
        if(individualTrackPt){

          for(auto trackPtBin : comparedTrackPtBin){
            temporaryTrackPtVector.clear();
            temporaryTrackPtVector.push_back(trackPtBin);

            binningInformation.push_back(std::make_tuple(temporaryCentralityVector, comparedJetPtBin, temporaryTrackPtVector));

          }

          // All defined track pT bins in a single plot
          } else {

            binningInformation.push_back(std::make_tuple(temporaryCentralityVector, comparedJetPtBin, comparedTrackPtBin));

          } // Track pT binning if-else

      } // Jet pT binning if-else

    } // Loop over centrality bins

  // All defined centrality bins in a single plot 
  } else {

    // Separate plots for each jet pT bin
    if(individualJetPt){

      for(auto jetPtBin : comparedJetPtBin){
        temporaryJetPtVector.clear();
        temporaryJetPtVector.push_back(jetPtBin);

        // Separate track pT bins for each bin
        if(individualTrackPt){

          for(auto trackPtBin : comparedTrackPtBin){
            temporaryTrackPtVector.clear();
            temporaryTrackPtVector.push_back(trackPtBin);

            binningInformation.push_back(std::make_tuple(comparedCentralityBin, temporaryJetPtVector, temporaryTrackPtVector));

          }

        // All defined track pT bins in a single plot
        } else {

          binningInformation.push_back(std::make_tuple(comparedCentralityBin, temporaryJetPtVector, comparedTrackPtBin));

        } // Track pT binning if-else

      } // Loop over ket pT bins

    // All defined jet pT bins in a single plot
    } else {

      // Separate track pT bins for each bin
      if(individualTrackPt){

        for(auto trackPtBin : comparedTrackPtBin){
          temporaryTrackPtVector.clear();
          temporaryTrackPtVector.push_back(trackPtBin);

          binningInformation.push_back(std::make_tuple(comparedCentralityBin, comparedJetPtBin, temporaryTrackPtVector));

        }

        // All defined track pT bins in a single plot
        } else {

          binningInformation.push_back(std::make_tuple(comparedCentralityBin, comparedJetPtBin, comparedTrackPtBin));

        } // Track pT binning if-else

    } // Jet pT binning if-else
    
  } // Centrality binning if-else
  

  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager* histograms;
  histograms = new EECHistogramManager(inputFile, card);

  // Choose the energy-energy correlator types to load
  histograms->SetLoadEnergyEnergyCorrelators(true);

  // Choose the bin ranges
  histograms->SetCentralityBinRange(0, card->GetNCentralityBins() - 1);
  histograms->SetJetPtBinRangeEEC(0, card->GetNJetPtBinsEEC() - 1);
  histograms->SetTrackPtBinRangeEEC(0, card->GetNTrackPtBinsEEC() - 1);

  // Load the histograms from the file
  histograms->LoadProcessedHistograms();

  // Energy-energy correlator histograms
  TH1D* hEnergyEnergyCorrelator[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];  // Raw energy-energy correlators
  TH1D* hEnergyEnergyCorrelatorUnfolded[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];  // Unfolded energy-energy correlators
  TH1D* hEnergyEnergyCorrelatorRelativeUncertainty[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];  // Relative uncertainty in each bin
  TH1D* hEnergyEnergyCorrelatorUnfoldedRelativeUncertainty[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];  // Relative uncertainty in each bin
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt] = NULL;
        hEnergyEnergyCorrelatorUnfolded[iCentrality][iJetPt][iTrackPt] = NULL;
        hEnergyEnergyCorrelatorRelativeUncertainty[iCentrality][iJetPt][iTrackPt] = NULL;
        hEnergyEnergyCorrelatorUnfoldedRelativeUncertainty[iCentrality][iJetPt][iTrackPt] = NULL;
      } // Centrality loop
    } // Track pT loop
  } // Jet pT loop

  // Helper histograms
  double epsilon = 0.0001;
  int lowNormalizationBin, highNormalizationBin;
  int iCentrality;
  int iTrackPt;
  int iJetPt;
  double binContent, binError;

  // Get the histograms from the histogram manager and normalize the signal histograms to one
  for(auto jetPtBin : comparedJetPtBin){
    for(auto trackPtBin : comparedTrackPtBin){
      for(auto centralityBin : comparedCentralityBin){

        // Binning indices
        iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);
        iTrackPt = card->GetBinIndexTrackPtEEC(trackPtBin);
        iCentrality = isPbPbData ? card->FindBinIndexCentrality(centralityBin) : 0;

        // Load the raw energy-energy correlator histograms

        hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt] = histograms->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt);

        // Normalize the distributions to one in the drawingRange
        lowNormalizationBin = hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(drawingRange.first + epsilon);
        highNormalizationBin = hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(drawingRange.second - epsilon);

        hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Scale(1 / hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Integral(lowNormalizationBin, highNormalizationBin, "width"));

        // Calculate the relative uncertainty in each bin
        hEnergyEnergyCorrelatorRelativeUncertainty[iCentrality][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Clone(Form("relativeUncertainty%d%d%d", iCentrality, iJetPt, iTrackPt));
        for(int iBin = 1; iBin <= hEnergyEnergyCorrelatorRelativeUncertainty[iCentrality][iJetPt][iTrackPt]->GetNbinsX(); iBin++){
          binContent = hEnergyEnergyCorrelatorRelativeUncertainty[iCentrality][iJetPt][iTrackPt]->GetBinContent(iBin);
          binError = hEnergyEnergyCorrelatorRelativeUncertainty[iCentrality][iJetPt][iTrackPt]->GetBinError(iBin);
          hEnergyEnergyCorrelatorRelativeUncertainty[iCentrality][iJetPt][iTrackPt]->SetBinContent(iBin, binError / binContent);
          hEnergyEnergyCorrelatorRelativeUncertainty[iCentrality][iJetPt][iTrackPt]->SetBinError(iBin, 0);
        }

        // Load the unfolded energy-energy correlator histograms

        hEnergyEnergyCorrelatorUnfolded[iCentrality][iJetPt][iTrackPt] = histograms->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

        // Normalize the distributions to one in the drawingRange
        lowNormalizationBin = hEnergyEnergyCorrelatorUnfolded[iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(drawingRange.first + epsilon);
        highNormalizationBin = hEnergyEnergyCorrelatorUnfolded[iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(drawingRange.second - epsilon);

        hEnergyEnergyCorrelatorUnfolded[iCentrality][iJetPt][iTrackPt]->Scale(1 / hEnergyEnergyCorrelatorUnfolded[iCentrality][iJetPt][iTrackPt]->Integral(lowNormalizationBin, highNormalizationBin, "width"));

        // Calculate the relative uncertainty in each bin
        hEnergyEnergyCorrelatorUnfoldedRelativeUncertainty[iCentrality][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelatorUnfolded[iCentrality][iJetPt][iTrackPt]->Clone(Form("relativeUncertaintyUnfolded%d%d%d", iCentrality, iJetPt, iTrackPt));
        for(int iBin = 1; iBin <= hEnergyEnergyCorrelatorUnfoldedRelativeUncertainty[iCentrality][iJetPt][iTrackPt]->GetNbinsX(); iBin++){
          binContent = hEnergyEnergyCorrelatorUnfoldedRelativeUncertainty[iCentrality][iJetPt][iTrackPt]->GetBinContent(iBin);
          binError = hEnergyEnergyCorrelatorUnfoldedRelativeUncertainty[iCentrality][iJetPt][iTrackPt]->GetBinError(iBin);
          hEnergyEnergyCorrelatorUnfoldedRelativeUncertainty[iCentrality][iJetPt][iTrackPt]->SetBinContent(iBin, binError / binContent);
          hEnergyEnergyCorrelatorUnfoldedRelativeUncertainty[iCentrality][iJetPt][iTrackPt]->SetBinError(iBin, 0);
        }

      } // Centrality loop
    } // Track pT loop
  } // Jet pT loop

  
  // ==========================================================================
  //                    All the ratios in the same figure
  // ==========================================================================
  
  JDrawer* drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  
  drawer->SetLogX(true);

  TString compactCentralityString = "";
  TString compactJetPtString = "";
  TString compactTrackPtString = "";
  TString comparedVariableString = "";
  TString legendString;
  int markerStyle[2] = {kFullCircle, kOpenSquare};
  int color[] = {kBlack,kRed,kBlue,kGreen+3,kMagenta,kCyan,kOrange,kViolet+3,kPink-7,kSpring+3,kAzure-7};

  // Binning vectors
  std::vector<int> currentCentralityIndices;
  std::vector<int> currentJetPtIndices;
  std::vector<int> currentTrackPtIndices;
  bool colorWithCentrality = false;
  bool colorWithJetPt = false;
  bool colorWithTrackPt = false;
  int colorFinder = 0;
  int firstCentralityBin = 0;
  int firstTrackPtBin = 0;
  int firstJetPtBin = 0;
  int legendCentralityIndex = 0;
  int legendJetPtIndex = 0;
  int legendTrackPtIndex = 0;
  TString individualLegend;

  for(auto plottedBin : binningInformation){

    // Find the indices for each drawn bin
    currentCentralityIndices.clear();
    if(isPbPbData){
      for(auto centralityBin : std::get<kCentrality>(plottedBin)){
        currentCentralityIndices.push_back(card->FindBinIndexCentrality(centralityBin));
      }
    } else {
      currentCentralityIndices.push_back(0);
    }

    currentJetPtIndices.clear();
    for(auto jetPtBin : std::get<kJetPt>(plottedBin)){
      currentJetPtIndices.push_back(card->FindBinIndexJetPtEEC(jetPtBin));
    }

    currentTrackPtIndices.clear();
    for(auto trackPtBin : std::get<kTrackPt>(plottedBin)){
      currentTrackPtIndices.push_back(card->GetBinIndexTrackPtEEC(trackPtBin));
    }

    // If the lenght of the indices array is larger than 1, we are plotting all these bins to the same plot
    colorWithCentrality = (currentCentralityIndices.size() > 1);
    colorWithTrackPt = (currentTrackPtIndices.size() > 1);
    colorWithJetPt = (currentJetPtIndices.size() > 1);

    // Remember the first drawn bin index
    firstCentralityBin = currentCentralityIndices.at(0);
    firstTrackPtBin = currentTrackPtIndices.at(0);
    firstJetPtBin = currentJetPtIndices.at(0);

    TLegend* legend = new TLegend(0.4,0.56,0.7,0.9);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);

    // Add common legend variables and define figure naming in case figures are saved
    if(isPbPbData){
      if(!colorWithCentrality){ 
        legend->AddEntry((TObject*) 0, Form("Cent: %.0f-%.0f%%", std::get<kCentrality>(plottedBin).at(0).first, std::get<kCentrality>(plottedBin).at(0).second), "");
        compactCentralityString = Form("_C=%.0f-%.0f", std::get<kCentrality>(plottedBin).at(0).first, std::get<kCentrality>(plottedBin).at(0).second);
      } else {
        compactCentralityString = "";
        comparedVariableString = "_centralityComparison";
      }
    }

    if(!colorWithJetPt) {
      legend->AddEntry((TObject*) 0, Form("%.0f < jet p_{T} < %.0f", std::get<kJetPt>(plottedBin).at(0).first, std::get<kJetPt>(plottedBin).at(0).second), "");
      compactJetPtString = Form("_J=%.0f-%.0f", std::get<kJetPt>(plottedBin).at(0).first, std::get<kJetPt>(plottedBin).at(0).second);
    } else {
      compactJetPtString = "";
      comparedVariableString = "_jetPtComparison";
    }

    if(!colorWithTrackPt){ 
      legend->AddEntry((TObject*) 0, Form("%.1f < track p_{T}", std::get<kTrackPt>(plottedBin).at(0)), "");
      compactTrackPtString = Form("_T>%.1f",std::get<kTrackPt>(plottedBin).at(0));
      compactTrackPtString.ReplaceAll(".","v");
    } else {
      compactTrackPtString = "";
      comparedVariableString = "_trackPtComparison";
    }

    // Set drawing style for all histograms
    colorFinder = 0;
    for(int iJetPt : currentJetPtIndices){
      for(int iTrackPt : currentTrackPtIndices){
        for(int iCentrality : currentCentralityIndices){
          hEnergyEnergyCorrelatorRelativeUncertainty[iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[0]);
          hEnergyEnergyCorrelatorRelativeUncertainty[iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[colorFinder]);
          hEnergyEnergyCorrelatorRelativeUncertainty[iCentrality][iJetPt][iTrackPt]->SetLineColor(color[colorFinder]);
          hEnergyEnergyCorrelatorUnfoldedRelativeUncertainty[iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[0]);
          hEnergyEnergyCorrelatorUnfoldedRelativeUncertainty[iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[colorFinder]);
          hEnergyEnergyCorrelatorUnfoldedRelativeUncertainty[iCentrality][iJetPt][iTrackPt]->SetLineColor(color[colorFinder]); 
        } // Centrality binning
        if(colorWithTrackPt) colorFinder++;
      } // Track pT binning 
      if(colorWithJetPt) colorFinder++;
    } // Jet pT binning

    // Set the x-axis drawing range
    hEnergyEnergyCorrelatorRelativeUncertainty[firstCentralityBin][firstJetPtBin][firstTrackPtBin]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
    hEnergyEnergyCorrelatorRelativeUncertainty[firstCentralityBin][firstJetPtBin][firstTrackPtBin]->GetYaxis()->SetRangeUser(eecZoom.first, eecZoom.second);
          
    // Draw the histograms to the upper canvas
    drawer->DrawHistogram(hEnergyEnergyCorrelatorRelativeUncertainty[firstCentralityBin][firstJetPtBin][firstTrackPtBin], "#Deltar", "EEC relative error", " ");

    for(int iCentrality : currentCentralityIndices){
      for(int iJetPt : currentJetPtIndices){
        for(int iTrackPt : currentTrackPtIndices){
          if(iCentrality == firstCentralityBin && iJetPt == firstJetPtBin && iTrackPt == firstTrackPtBin) continue;
          hEnergyEnergyCorrelatorRelativeUncertainty[iCentrality][iJetPt][iTrackPt]->Draw("same");
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop

    // Add legends for drawn histograms
    legendCentralityIndex = 0;
    legendJetPtIndex = 0;
    legendTrackPtIndex = 0;
    individualLegend = "";
    for(int iJetPt : currentJetPtIndices){
      if(colorWithJetPt) individualLegend = Form(" %.0f < jet p_{T} < %.0f", std::get<kJetPt>(plottedBin).at(legendJetPtIndex).first, std::get<kJetPt>(plottedBin).at(legendJetPtIndex).second);
      legendJetPtIndex++;
      for(int iTrackPt : currentTrackPtIndices){
        if(colorWithTrackPt) individualLegend = Form(" %.1f < track p_{T}", std::get<kTrackPt>(plottedBin).at(legendTrackPtIndex++));
        for(int iCentrality : currentCentralityIndices){
          if(colorWithCentrality && isPbPbData) individualLegend = Form(" Cent: %.0f-%.0f%%", std::get<kCentrality>(plottedBin).at(legendCentralityIndex).first, std::get<kCentrality>(plottedBin).at(legendCentralityIndex).second);
          legendCentralityIndex++;
          legend->AddEntry(hEnergyEnergyCorrelatorRelativeUncertainty[iCentrality][iJetPt][iTrackPt], Form("%s %s", collisionSystem.Data(), individualLegend.Data()), "l");
        } // Centrality loop 
      } // Track pT loop
    } // Jet pT loop
  
    // Draw the legends to the upper pad
    legend->Draw();
          
    // Save the figures to a file
    if(saveFigures){
    gPad->GetCanvas()->SaveAs(Form("figures/eecRelativeUncertainty%s%s%s%s%s.%s", saveComment, comparedVariableString.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
    }

    // Draw also the relative uncertainties after unfolding

    legend = new TLegend(0.4,0.56,0.7,0.9);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);

    // Add common legend variables and define figure naming in case figures are saved
    if(isPbPbData){
      if(!colorWithCentrality){ 
        legend->AddEntry((TObject*) 0, Form("Cent: %.0f-%.0f%%", std::get<kCentrality>(plottedBin).at(0).first, std::get<kCentrality>(plottedBin).at(0).second), "");
        compactCentralityString = Form("_C=%.0f-%.0f", std::get<kCentrality>(plottedBin).at(0).first, std::get<kCentrality>(plottedBin).at(0).second);
      } else {
        compactCentralityString = "";
        comparedVariableString = "_centralityComparison";
      }
    }

    if(!colorWithJetPt) {
      legend->AddEntry((TObject*) 0, Form("%.0f < jet p_{T} < %.0f", std::get<kJetPt>(plottedBin).at(0).first, std::get<kJetPt>(plottedBin).at(0).second), "");
      compactJetPtString = Form("_J=%.0f-%.0f", std::get<kJetPt>(plottedBin).at(0).first, std::get<kJetPt>(plottedBin).at(0).second);
    } else {
      compactJetPtString = "";
      comparedVariableString = "_jetPtComparison";
    }

    if(!colorWithTrackPt){ 
      legend->AddEntry((TObject*) 0, Form("%.1f < track p_{T}", std::get<kTrackPt>(plottedBin).at(0)), "");
      compactTrackPtString = Form("_T>%.1f",std::get<kTrackPt>(plottedBin).at(0));
      compactTrackPtString.ReplaceAll(".","v");
    } else {
      compactTrackPtString = "";
      comparedVariableString = "_trackPtComparison";
    }

    // Set the x-axis drawing range
    hEnergyEnergyCorrelatorUnfoldedRelativeUncertainty[firstCentralityBin][firstJetPtBin][firstTrackPtBin]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
    hEnergyEnergyCorrelatorUnfoldedRelativeUncertainty[firstCentralityBin][firstJetPtBin][firstTrackPtBin]->GetYaxis()->SetRangeUser(eecZoom.first, eecZoom.second);
          
    // Draw the histograms to the upper canvas
    drawer->DrawHistogram(hEnergyEnergyCorrelatorUnfoldedRelativeUncertainty[firstCentralityBin][firstJetPtBin][firstTrackPtBin], "#Deltar", "Unfolded signal EEC relative error", " ");

    for(int iCentrality : currentCentralityIndices){
      for(int iJetPt : currentJetPtIndices){
        for(int iTrackPt : currentTrackPtIndices){
          if(iCentrality == firstCentralityBin && iJetPt == firstJetPtBin && iTrackPt == firstTrackPtBin) continue;
          hEnergyEnergyCorrelatorUnfoldedRelativeUncertainty[iCentrality][iJetPt][iTrackPt]->Draw("same");
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop

    // Add legends for drawn histograms
    legendCentralityIndex = 0;
    legendJetPtIndex = 0;
    legendTrackPtIndex = 0;
    individualLegend = "";
    for(int iJetPt : currentJetPtIndices){
      if(colorWithJetPt) individualLegend = Form(" %.0f < jet p_{T} < %.0f", std::get<kJetPt>(plottedBin).at(legendJetPtIndex).first, std::get<kJetPt>(plottedBin).at(legendJetPtIndex).second);
      legendJetPtIndex++;
      for(int iTrackPt : currentTrackPtIndices){
        if(colorWithTrackPt) individualLegend = Form(" %.1f < track p_{T}", std::get<kTrackPt>(plottedBin).at(legendTrackPtIndex++));
        for(int iCentrality : currentCentralityIndices){
          if(colorWithCentrality && isPbPbData) individualLegend = Form(" Cent: %.0f-%.0f%%", std::get<kCentrality>(plottedBin).at(legendCentralityIndex).first, std::get<kCentrality>(plottedBin).at(legendCentralityIndex).second);
          legendCentralityIndex++;
          legend->AddEntry(hEnergyEnergyCorrelatorUnfoldedRelativeUncertainty[iCentrality][iJetPt][iTrackPt], Form("%s %s", collisionSystem.Data(), individualLegend.Data()), "l");
        } // Centrality loop 
      } // Track pT loop
    } // Jet pT loop
  
    // Draw the legends to the upper pad
    legend->Draw();
          
    // Save the figures to a file
    if(saveFigures){
    gPad->GetCanvas()->SaveAs(Form("figures/unfoldedEecRelativeUncertainty%s%s%s%s%s.%s", saveComment, comparedVariableString.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
    }
  }

}
