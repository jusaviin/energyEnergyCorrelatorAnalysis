#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Macro for comparing final energy-energy correlators results
 */
void compareFinalEECResults(){
  
  // Files for comparison
  const int nComparisonFiles = 1;
  TString pbPbFileName[nComparisonFiles];
  //pbPbFileName[0] = "data/eecAnalysis_akFlowJet_wtaAxis_newTrackPairEfficiencySmoothed_unfoldingWithNominalSmear_processed_2023-07-13.root";
  pbPbFileName[0] = "data/eecAnalysis_akFlowJet_wtaAxis_energyWeightSquared_firstFinalResultsWithFixedCard_processed_2023-10-23.root";

  TString ppFileName[nComparisonFiles];
  //ppFileName[0] = "data/ppData_pfJets_wtaAxis_newTrackPairEfficiency_unfoldingWithNominalSmear_processed_2023-07-13.root";
  ppFileName[0] = "data/ppData_pfJets_wtaAxis_energyWeightSquared_jet60or80triggers_firstFinalResults_processed_2023-10-26.root";
  
  
  // Open the files and check that they exist
  TFile* pbPbInputFile[nComparisonFiles];
  EECCard* pbPbCard[nComparisonFiles];
  TFile* ppInputFile[nComparisonFiles];
  EECCard* ppCard[nComparisonFiles];
  TString energyWeightString[nComparisonFiles];
  int weightExponent;
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    
    pbPbInputFile[iFile] = TFile::Open(pbPbFileName[iFile]);
    
    if(pbPbInputFile[iFile] == NULL){
      cout << "Error! The file " << pbPbFileName[iFile].Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }

    pbPbCard[iFile] = new EECCard(pbPbInputFile[iFile]);

    ppInputFile[iFile] = TFile::Open(ppFileName[iFile]);
    
    if(ppInputFile[iFile] == NULL){
      cout << "Error! The file " << ppFileName[iFile].Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }

    ppCard[iFile] = new EECCard(ppInputFile[iFile]);

    // Check that the energy weights from PbPb and pp match. Otherwise taking ratio is meaningless
    weightExponent = pbPbCard[iFile]->GetWeightExponent();
    if(weightExponent != ppCard[iFile]->GetWeightExponent()){
      cout << "ERROR! You are trying to take a PbPb to pp ratio with different energy weight exponents!" << endl;
      cout << "I cannot allow you to do this. The ratio would not be what you want." << endl;
      cout << "Please check that you are comparing files with the same energy weight exponents." << endl;
      return;
    }

    // Put the correct energy weight to a string
    if(weightExponent == 1){
      energyWeightString[iFile] = "p_{T,1}p_{T,2}";
    } else {
      energyWeightString[iFile] = Form("p_{T,1}^{%d}p_{T,2}^{%d}", weightExponent, weightExponent);
    }
  }
  
  // ====================================================
  //               Binning configuration
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = pbPbCard[0]->GetNCentralityBins();
  const int nJetPtBinsEEC = pbPbCard[0]->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = pbPbCard[0]->GetNTrackPtBinsEEC();
  
  // Select explicitly which bins from the files are compared:
  std::vector<std::pair<double,double>> comparedCentralityBin;
  comparedCentralityBin.push_back(std::make_pair(0,10));
  //comparedCentralityBin.push_back(std::make_pair(10,30));
  //comparedCentralityBin.push_back(std::make_pair(30,50));
  //comparedCentralityBin.push_back(std::make_pair(50,90));
  bool individualCentrality = true; // True = make different figure for each bin. False = plot all centrality bin to the same figure.

  std::vector<std::pair<double,double>> comparedJetPtBin;
  comparedJetPtBin.push_back(std::make_pair(120,140));
  comparedJetPtBin.push_back(std::make_pair(140,160));
  comparedJetPtBin.push_back(std::make_pair(160,180));
  comparedJetPtBin.push_back(std::make_pair(180,200));
  bool individualJetPt = true; // True = make different figure for each bin. False = plot all jet pT bin to the same figure.

  std::vector<double> comparedTrackPtBin;
  comparedTrackPtBin.push_back(1.0);
  comparedTrackPtBin.push_back(1.5);
  comparedTrackPtBin.push_back(2.0);
  comparedTrackPtBin.push_back(2.5);
  comparedTrackPtBin.push_back(3.0);
  bool individualTrackPt = false; // True = make different figure for each bin. False = plot all track pT bin to the same figure.

  // ====================================================
  //                Drawing configuration
  // ====================================================
  
  // Figure saving
  const bool saveFigures = false;  // Save figures
  const char* saveComment = "_energyWeightComparison";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures

  // Drawing configuration
  std::pair<double, double> ratioZoom = std::make_pair(0.4, 1.6);
  std::pair<double, double> eecZoom = std::make_pair(0.05, 50);

  // Sanity checks for input. Ensure that all the selected bins actually exist in the input files
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){

    // Sanity check for centrality bins
    for(auto centralityBin : comparedCentralityBin){
      if(pbPbCard[iFile]->FindBinIndexCentrality(centralityBin) < pbPbCard[iFile]->GetFirstUnfoldedCentralityBin() || pbPbCard[iFile]->FindBinIndexCentrality(centralityBin) > pbPbCard[iFile]->GetLastUnfoldedCentralityBin()){
        cout << "ERROR! Centrality bin " << centralityBin.first << "-" << centralityBin.second << " does not exist in file " << pbPbFileName[iFile].Data() << endl;
        cout << "Please only choose centrality bins that are included in the input files." << endl;
        return;
      } 
    }

    // Sanity check for jet pT bins
    for(auto jetPtBin : comparedJetPtBin){
      if(pbPbCard[iFile]->FindBinIndexJetPtEEC(jetPtBin) < pbPbCard[iFile]->GetFirstUnfoldedJetPtBin() || pbPbCard[iFile]->FindBinIndexJetPtEEC(jetPtBin) > pbPbCard[iFile]->GetLastUnfoldedJetPtBin()){
        cout << "ERROR! Jet pT bin " << jetPtBin.first << "-" << jetPtBin.second << " does not exist in file " << pbPbFileName[iFile].Data() << endl;
        cout << "Please only choose jet pT bins that are included in the input files." << endl;
        return;
      }
      if(ppCard[iFile]->FindBinIndexJetPtEEC(jetPtBin) < ppCard[iFile]->GetFirstUnfoldedJetPtBin() || ppCard[iFile]->FindBinIndexJetPtEEC(jetPtBin) > ppCard[iFile]->GetLastUnfoldedJetPtBin()){
        cout << "ERROR! Jet pT bin " << jetPtBin.first << "-" << jetPtBin.second << " does not exist in file " << ppFileName[iFile].Data() << endl;
        cout << "Please only choose jet pT bins that are included in the input files." << endl;
        return;
      } 
    }

    // Sanity check for track pT bins
    for(auto trackPtBin : comparedTrackPtBin){
      if(pbPbCard[iFile]->GetBinIndexTrackPtEEC(trackPtBin) < pbPbCard[iFile]->GetFirstUnfoldedTrackPtBin() || pbPbCard[iFile]->GetBinIndexTrackPtEEC(trackPtBin) > pbPbCard[iFile]->GetLastUnfoldedTrackPtBin()){
        cout << "ERROR! Track pT cut > " << trackPtBin << " GeV does not exist in file " << pbPbFileName[iFile].Data() << endl;
        cout << "Please only choose track pT bins that are included in the input files." << endl;
        return;
      }
      if(ppCard[iFile]->GetBinIndexTrackPtEEC(trackPtBin) < ppCard[iFile]->GetFirstUnfoldedTrackPtBin() || ppCard[iFile]->GetBinIndexTrackPtEEC(trackPtBin) > ppCard[iFile]->GetLastUnfoldedTrackPtBin()){
        cout << "ERROR! Track pT cut > " << trackPtBin << " GeV does not exist in file " << ppFileName[iFile].Data() << endl;
        cout << "Please only choose track pT bins that are included in the input files." << endl;
        return;
      }
    } 
  } // File loop for input sanity check

  // Only allow one variable for which all bins are plotted to the same figure
  if(individualCentrality + individualJetPt + individualTrackPt < 2){
    cout << "You are tring to plot too many bins to the same figure!" << endl;
    cout << "This macro can only plot all the bins from one variable." << endl;
    cout << "Please check your configuration!" << endl; 
    return;
  }

  // Do not allow several bins plotted to same figure if several files are compared
  if(nComparisonFiles > 1 && (!individualTrackPt || !individualJetPt || !individualCentrality)){
    cout << "If you compare files, you cannot draw several bins from one variable to single figure." << endl;
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
  EECHistogramManager* pbPbHistograms[nComparisonFiles];
  EECHistogramManager* ppHistograms[nComparisonFiles];
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    pbPbHistograms[iFile] = new EECHistogramManager(pbPbInputFile[iFile], pbPbCard[iFile]);

    // Choose the energy-energy correlator types to load
    pbPbHistograms[iFile]->SetLoadEnergyEnergyCorrelators(true);

    // Choose the bin ranges
    pbPbHistograms[iFile]->SetCentralityBinRange(0, pbPbCard[iFile]->GetNCentralityBins() - 1);
    pbPbHistograms[iFile]->SetJetPtBinRangeEEC(0, pbPbCard[iFile]->GetNJetPtBinsEEC() - 1);
    pbPbHistograms[iFile]->SetTrackPtBinRangeEEC(0, pbPbCard[iFile]->GetNTrackPtBinsEEC() - 1);

    // Load the histograms from the file
    pbPbHistograms[iFile]->LoadProcessedHistograms();

    ppHistograms[iFile] = new EECHistogramManager(ppInputFile[iFile], ppCard[iFile]);

    // Choose the energy-energy correlator types to load
    ppHistograms[iFile]->SetLoadEnergyEnergyCorrelators(true);

    // Choose the bin ranges
    ppHistograms[iFile]->SetCentralityBinRange(0, 1);
    ppHistograms[iFile]->SetJetPtBinRangeEEC(0, ppCard[iFile]->GetNJetPtBinsEEC() - 1);
    ppHistograms[iFile]->SetTrackPtBinRangeEEC(0, ppCard[iFile]->GetNTrackPtBinsEEC() - 1);

    // Load the histograms from the file
    ppHistograms[iFile]->LoadProcessedHistograms();
  }

  // Energy-energy correlator histograms
  TH1D* hEnergyEnergyCorrelatorPbPb[nComparisonFiles][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];  // Energy-energy correlators for PbPb
  TH1D* hEnergyEnergyCorrelatorPp[nComparisonFiles][nJetPtBinsEEC][nTrackPtBinsEEC];    // Energy-energy correlators for pp
  TH1D* hEnergyEnergyCorrelatorRatio[nComparisonFiles][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC]; // Ratio between PbPb and pp distributions
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        hEnergyEnergyCorrelatorPp[iFile][iJetPt][iTrackPt] = NULL;
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          hEnergyEnergyCorrelatorPbPb[iFile][iCentrality][iJetPt][iTrackPt] = NULL;
          hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt] = NULL;
        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // File loop
  
  // Helper histograms
  std::pair<double, double> drawingRange = std::make_pair(0.006, 0.39);
  double epsilon = 0.0001;
  int lowNormalizationBin, highNormalizationBin;
  int iCentrality, iCentralityReference;
  int iTrackPt, iTrackPtReference;
  int iJetPt, iJetPtReference;

  // Get the histograms from the histogram manager and normalize the signal histograms to one
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    for(auto jetPtBin : comparedJetPtBin){
      for(auto trackPtBin : comparedTrackPtBin){

        // Load the final unfolded signal histograms for pp
        iCentrality = 0;
        iJetPt = ppCard[iFile]->FindBinIndexJetPtEEC(jetPtBin);
        iTrackPt = ppCard[iFile]->GetBinIndexTrackPtEEC(trackPtBin);
        iJetPtReference = pbPbCard[0]->FindBinIndexJetPtEEC(jetPtBin);
        iTrackPtReference = pbPbCard[0]->GetBinIndexTrackPtEEC(trackPtBin);

        // Read the histogram using pp indices, fill it with PbPb indices from file 0 to facilitate drawing later
        hEnergyEnergyCorrelatorPp[iFile][iJetPtReference][iTrackPtReference] = ppHistograms[iFile]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

        // Normalize the distributions to one in the drawingRange
        lowNormalizationBin = hEnergyEnergyCorrelatorPp[iFile][iJetPtReference][iTrackPtReference]->GetXaxis()->FindBin(drawingRange.first + epsilon);
        highNormalizationBin = hEnergyEnergyCorrelatorPp[iFile][iJetPtReference][iTrackPtReference]->GetXaxis()->FindBin(drawingRange.second - epsilon);

        hEnergyEnergyCorrelatorPp[iFile][iJetPtReference][iTrackPtReference]->Scale(1 / hEnergyEnergyCorrelatorPp[iFile][iJetPtReference][iTrackPtReference]->Integral(lowNormalizationBin, highNormalizationBin, "width"));

        for(auto centralityBin : comparedCentralityBin){

          // Load the final unfolded signal histograms for PbPb
          iJetPt = pbPbCard[iFile]->FindBinIndexJetPtEEC(jetPtBin);
          iTrackPt = pbPbCard[iFile]->GetBinIndexTrackPtEEC(trackPtBin);
          iCentrality = pbPbCard[iFile]->FindBinIndexCentrality(centralityBin);
          iCentralityReference = pbPbCard[0]->FindBinIndexCentrality(centralityBin);

          hEnergyEnergyCorrelatorPbPb[iFile][iCentralityReference][iJetPtReference][iTrackPtReference] = pbPbHistograms[iFile]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

          // Normalize the distributions to one in the drawingRange
          lowNormalizationBin = hEnergyEnergyCorrelatorPbPb[iFile][iCentralityReference][iJetPtReference][iTrackPtReference]->GetXaxis()->FindBin(drawingRange.first + epsilon);
          highNormalizationBin = hEnergyEnergyCorrelatorPbPb[iFile][iCentralityReference][iJetPtReference][iTrackPtReference]->GetXaxis()->FindBin(drawingRange.second - epsilon);

          hEnergyEnergyCorrelatorPbPb[iFile][iCentralityReference][iJetPtReference][iTrackPtReference]->Scale(1 / hEnergyEnergyCorrelatorPbPb[iFile][iCentralityReference][iJetPtReference][iTrackPtReference]->Integral(lowNormalizationBin, highNormalizationBin, "width"));

          // Calculate the PbPb to pp ratio
          hEnergyEnergyCorrelatorRatio[iFile][iCentralityReference][iJetPtReference][iTrackPtReference] = (TH1D*)hEnergyEnergyCorrelatorPbPb[iFile][iCentralityReference][iJetPtReference][iTrackPtReference]->Clone(Form("eecRatio%d", iFile));
          hEnergyEnergyCorrelatorRatio[iFile][iCentralityReference][iJetPtReference][iTrackPtReference]->Divide(hEnergyEnergyCorrelatorPp[iFile][iJetPtReference][iTrackPtReference]);

        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // File loop
  
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
    for(auto centralityBin : std::get<kCentrality>(plottedBin)){
      currentCentralityIndices.push_back(pbPbCard[0]->FindBinIndexCentrality(centralityBin));
    }

    currentJetPtIndices.clear();
    for(auto jetPtBin : std::get<kJetPt>(plottedBin)){
      currentJetPtIndices.push_back(pbPbCard[0]->FindBinIndexJetPtEEC(jetPtBin));
    }

    currentTrackPtIndices.clear();
    for(auto trackPtBin : std::get<kTrackPt>(plottedBin)){
      currentTrackPtIndices.push_back(pbPbCard[0]->GetBinIndexTrackPtEEC(trackPtBin));
    }

    // If the lenght of the indices array is larger than 1, we are plotting all these bins to the same plot
    colorWithCentrality = (currentCentralityIndices.size() > 1);
    colorWithTrackPt = (currentTrackPtIndices.size() > 1);
    colorWithJetPt = (currentJetPtIndices.size() > 1);

    // Remember the first drawn bin index
    firstCentralityBin = currentCentralityIndices.at(0);
    firstTrackPtBin = currentTrackPtIndices.at(0);
    firstJetPtBin = currentJetPtIndices.at(0);

    // Create a new canvas for the plot
    drawer->CreateSplitCanvas();
          
    // Logarithmic EEC axis
    drawer->SetLogY(true);

    TLegend* legend = new TLegend(0.18,0.04,0.45,0.58);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);

    // Add common legend variables and define figure naming in case figures are saved
    if(!colorWithCentrality){ 
      legend->AddEntry((TObject*) 0, Form("Cent: %.0f-%.0f%%", std::get<kCentrality>(plottedBin).at(0).first, std::get<kCentrality>(plottedBin).at(0).second), "");
      compactCentralityString = Form("_C=%.0f-%.0f", std::get<kCentrality>(plottedBin).at(0).first, std::get<kCentrality>(plottedBin).at(0).second);
    } else {
      compactCentralityString = "";
      comparedVariableString = "_centralityComparison";
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
    for(int iFile = 0; iFile < nComparisonFiles; iFile++){
      for(int iJetPt : currentJetPtIndices){
        for(int iTrackPt : currentTrackPtIndices){
          hEnergyEnergyCorrelatorPp[iFile][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[1]);
          hEnergyEnergyCorrelatorPp[iFile][iJetPt][iTrackPt]->SetMarkerColor(color[iFile+colorFinder]);
          hEnergyEnergyCorrelatorPp[iFile][iJetPt][iTrackPt]->SetLineColor(color[iFile+colorFinder]);
          for(int iCentrality : currentCentralityIndices){
            if(colorWithCentrality) colorFinder++;
            hEnergyEnergyCorrelatorPbPb[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[0]);
            hEnergyEnergyCorrelatorPbPb[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iFile+colorFinder]);
            hEnergyEnergyCorrelatorPbPb[iFile][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iFile+colorFinder]); 
            hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[0]);
            hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iFile+colorFinder]);
            hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iFile+colorFinder]);
          } // Centrality binning
          if(colorWithTrackPt) colorFinder++;
        } // Track pT binning 
        if(colorWithJetPt) colorFinder++;
      } // Jet pT binning
    } // File loop

    // Set the x-axis drawing range
    hEnergyEnergyCorrelatorPbPb[0][firstCentralityBin][firstJetPtBin][firstTrackPtBin]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
    hEnergyEnergyCorrelatorPbPb[0][firstCentralityBin][firstJetPtBin][firstTrackPtBin]->GetYaxis()->SetRangeUser(eecZoom.first, eecZoom.second);
          
    // Draw the histograms to the upper canvas
    drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelatorPbPb[0][firstCentralityBin][firstJetPtBin][firstTrackPtBin], "#Deltar", "EEC", " ");

    for(int iCentrality : currentCentralityIndices){
      for(int iJetPt : currentJetPtIndices){
        for(int iTrackPt : currentTrackPtIndices){
          for(int iFile = 0; iFile < nComparisonFiles; iFile++){
            if(iFile == 0 && iCentrality == firstCentralityBin && iJetPt == firstJetPtBin && iTrackPt == firstTrackPtBin) continue;
            hEnergyEnergyCorrelatorPbPb[iFile][iCentrality][iJetPt][iTrackPt]->Draw("same");
          } // File loop
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
    for(int iJetPt : currentJetPtIndices){
      for(int iTrackPt : currentTrackPtIndices){
        for(int iFile = 0; iFile < nComparisonFiles; iFile++){
          hEnergyEnergyCorrelatorPp[iFile][iJetPt][iTrackPt]->Draw("same");
        } // File loop
      } // Track pT loop
    } // Jet pT loop

    // Add legends for drawn histograms
    legendCentralityIndex = 0;
    legendJetPtIndex = 0;
    legendTrackPtIndex = 0;
    for(int iFile = 0; iFile < nComparisonFiles; iFile++){
      individualLegend = "";
      for(int iJetPt : currentJetPtIndices){
        if(colorWithJetPt) individualLegend = Form(" %.0f < jet p_{T} < %.0f", std::get<kJetPt>(plottedBin).at(legendJetPtIndex).first, std::get<kJetPt>(plottedBin).at(legendJetPtIndex).second);
        legendJetPtIndex++;
        for(int iTrackPt : currentTrackPtIndices){
          if(colorWithTrackPt) individualLegend = Form(" %.1f < track p_{T}", std::get<kTrackPt>(plottedBin).at(legendTrackPtIndex++));
          for(int iCentrality : currentCentralityIndices){
            if(colorWithCentrality) individualLegend = Form(" Cent: %.0f-%.0f%%", std::get<kCentrality>(plottedBin).at(legendCentralityIndex).first, std::get<kCentrality>(plottedBin).at(legendCentralityIndex).second);
            legendCentralityIndex++;
            legend->AddEntry(hEnergyEnergyCorrelatorPbPb[iFile][iCentrality][iJetPt][iTrackPt], Form("PbPb %s%s", energyWeightString[iFile].Data(), individualLegend.Data()), "p");
          } // Centrality loop 
          if(colorWithCentrality) individualLegend = "";
          legend->AddEntry(hEnergyEnergyCorrelatorPp[iFile][iJetPt][iTrackPt], Form("pp %s%s", energyWeightString[iFile].Data(), individualLegend.Data()), "p");
        } // Track pT loop
      } // Jet pT loop
    } // File loop
  
    // Draw the legends to the upper pad
    legend->Draw();
          
    // Linear scale for the ratio
    drawer->SetLogY(false);
          
    // Set the drawing ranges
    hEnergyEnergyCorrelatorRatio[0][firstCentralityBin][firstJetPtBin][firstTrackPtBin]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);          
    hEnergyEnergyCorrelatorRatio[0][firstCentralityBin][firstJetPtBin][firstTrackPtBin]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

    // Draw the histograms
    drawer->SetGridY(true);
    drawer->DrawHistogramToLowerPad(hEnergyEnergyCorrelatorRatio[0][firstCentralityBin][firstJetPtBin][firstTrackPtBin], "#Deltar", "#frac{PbPb}{pp}", " ");
    for(int iCentrality : currentCentralityIndices){
      for(int iJetPt : currentJetPtIndices){
        for(int iTrackPt : currentTrackPtIndices){
          for(int iFile = 0; iFile < nComparisonFiles; iFile++){
            if(iFile == 0 && iCentrality == firstCentralityBin && iJetPt == firstJetPtBin && iTrackPt == firstTrackPtBin) continue;
            hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->Draw("same");
          } // File loop
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
    drawer->SetGridY(false);
          
    // Save the figures to a file
    if(saveFigures){
    gPad->GetCanvas()->SaveAs(Form("figures/eecFinalResults%s%s%s%s%s.%s", saveComment, comparedVariableString.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
    }
  }

}
