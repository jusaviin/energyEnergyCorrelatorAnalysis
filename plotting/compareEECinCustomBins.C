#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"
#include "AlgorithmLibrary.h"

/*
 * Macro for comparing energy-energy correlators in custom bins
 *
 *  const int presetComparison = Preset comparison style for plots. Can be set to 0 to manually define everything
 *  const double lowDrawRange = Lowest drawn DeltaR value in the figure. The default 0.008 is used in data analysis.
 */
void compareEECinCustomBins(const int presetComparison = 0, const double lowDrawRange = 0.008){

  // Enumeration for cumulant scaling variants
  enum enumCumulantVariations{kJetShape, kCorrelator, kAdvancedCorrelator, knCumulants};

  TString cumulantName[knCumulants] = {"jetShape", "correlator", "advancedCorrelator"};
  TString legendText[knCumulants] = {"E1C cumulant scale", "E2C cumulant scale (old)", "E2C cumulant scale (new)"};
  
  // Files for comparison
  std::vector<TString> fileName;
  TString presetSaveComment = "";

  // Use this region to manually define the comparison files. Other files for preset comparisons are not meant to be changed.
  if(presetComparison == 0){
    fileName.push_back("data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_nominalEnergyWeight_includeJetShape_processed_2024-07-17.root");
    fileName.push_back("data/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_nominalEnergyWeight_jetShape_processed_2024-07-18.root");
  }

  if(presetComparison == 1){
    // Preset value 1: 5 GeV jet pT shift in Pythia8, nominal energy weight
    fileName.push_back("data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_nominalEnergyWeight_includeJetShape_processed_2024-07-17.root");
    fileName.push_back("data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_nominalEnergyWeight_includeJetShape_5GeVshiftedPt_processed_2024-07-18.root");
    presetSaveComment = "_pythia5GeVShift";
  } else if (presetComparison == 2){
    // Preset value 2: 15 GeV jet pT shift in Pythia8, nominal energy weight
    fileName.push_back("data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_nominalEnergyWeight_includeJetShape_processed_2024-07-17.root");
    fileName.push_back("data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_nominalEnergyWeight_includeJetShape_shiftedJetPt_processed_2024-07-17.root");
    presetSaveComment = "_pythia15GeVShift";
  } else if (presetComparison == 3){
    // Preset value 3: 5 GeV jet pT shift in Herwig7, nominal energy weight
    fileName.push_back("data/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_nominalEnergyWeight_jetShape_processed_2024-07-18.root");
    fileName.push_back("data/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_nominalEnergyWeight_jetShape_5GeVshiftedPt_processed_2024-07-18.root");
    presetSaveComment = "_herwig5GeVShift";
  } else if (presetComparison == 4){
    // Preset value 4: 15 GeV jet pT shift in Herwig7, nominal energy weight
    fileName.push_back("data/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_nominalEnergyWeight_jetShape_processed_2024-07-18.root");
    fileName.push_back("data/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_nominalEnergyWeight_jetShape_15GeVshiftedPt_processed_2024-07-18.root");
    presetSaveComment = "_herwig15GeVShift";
  } else if (presetComparison == 5){
    // Preset value 5: Nominal Pythia8 compared with 5 GeV shifted Herwig7, nominal energy weight
    fileName.push_back("data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_nominalEnergyWeight_includeJetShape_processed_2024-07-17.root");
    fileName.push_back("data/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_nominalEnergyWeight_jetShape_5GeVshiftedPt_processed_2024-07-18.root");
    presetSaveComment = "_pythiaHerwigCombined5GeVShift";
  } else if (presetComparison == 6){
    // Preset value 6: Nominal Pythia8 compared with 15 GeV shifted Herwig7, nominal energy weight
    fileName.push_back("data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_nominalEnergyWeight_includeJetShape_processed_2024-07-17.root");
    fileName.push_back("data/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_nominalEnergyWeight_jetShape_15GeVshiftedPt_processed_2024-07-18.root");
    presetSaveComment = "_pythiaHerwigCombined15GeVShift";
  } else if(presetComparison == 7){
    // Preset value 7: 5 GeV jet pT shift in Pythia8, squared energy weight
    fileName.push_back("data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_includeJetShape_processed_2024-07-17.root");
    fileName.push_back("data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_includeJetShape_5GeVshiftedPt_processed_2024-07-18.root");
    presetSaveComment = "_pythia5GeVShift";
  } else if(presetComparison == 8){
    // Preset value 8: 15 GeV jet pT shift in Pythia8, squared energy weight
    fileName.push_back("data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_includeJetShape_processed_2024-07-17.root");
    fileName.push_back("data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_includeJetShape_shiftedJetPt_processed_2024-07-17.root");
    presetSaveComment = "_pythia15GeVShift";
  } else if(presetComparison == 9){
    // Preset value 9: 5 GeV jet pT shift in Herwig7, squared energy weight
    fileName.push_back("data/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_energyWeightSquared_jetShape_processed_2024-07-18.root");
    fileName.push_back("data/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_energyWeightSquared_jetShape_5GeVshiftedPt_processed_2024-07-18.root");
    presetSaveComment = "_herwig5GeVShift";
  } else if(presetComparison == 10){
    // Preset value 10: 15 GeV jet pT shift in Herwig7, squared energy weight
    fileName.push_back("data/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_energyWeightSquared_jetShape_processed_2024-07-18.root");
    fileName.push_back("data/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_energyWeightSquared_jetShape_15GeVshiftedPt_processed_2024-07-18.root");
    presetSaveComment = "_herwig15GeVShift";
  } else if(presetComparison == 11){
    // Preset value 11: Nominal Pythia8 compared with 5 GeV shifted Herwig7, squared energy weight
    fileName.push_back("data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_includeJetShape_processed_2024-07-17.root");
    fileName.push_back("data/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_energyWeightSquared_jetShape_5GeVshiftedPt_processed_2024-07-18.root");
    presetSaveComment = "_pythiaHerwigCombined5GeVShift";
  } else if(presetComparison == 12){
    // Preset value 12: Nominal Pythia8 compared with 15 GeV shifted Herwig7, squared energy weight
    fileName.push_back("data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_includeJetShape_processed_2024-07-17.root");
    fileName.push_back("data/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_energyWeightSquared_jetShape_15GeVshiftedPt_processed_2024-07-18.root");
    presetSaveComment = "_pythiaHerwigCombined15GeVShift";
  } else if(presetComparison == 13){
    // Preset value 13: 5 GeV shift in Pythia8 and Herwig7 drawn in the same plot
    fileName.push_back("data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_nominalEnergyWeight_includeJetShape_processed_2024-07-17.root");
    fileName.push_back("data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_nominalEnergyWeight_includeJetShape_5GeVshiftedPt_processed_2024-07-18.root");
    fileName.push_back("data/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_nominalEnergyWeight_jetShape_processed_2024-07-18.root");
    fileName.push_back("data/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_nominalEnergyWeight_jetShape_5GeVshiftedPt_processed_2024-07-18.root");
    presetSaveComment = "_pythiaHerwig5GeVShift";
  } else if(presetComparison == 14){
    // Preset value 14: 15 GeV shift in Pythia8 and Herwig7 drawn in the same plot
    fileName.push_back("data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_nominalEnergyWeight_includeJetShape_processed_2024-07-17.root");
    fileName.push_back("data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_nominalEnergyWeight_includeJetShape_shiftedJetPt_processed_2024-07-17.root");
    fileName.push_back("data/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_nominalEnergyWeight_jetShape_processed_2024-07-18.root");
    fileName.push_back("data/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_nominalEnergyWeight_jetShape_15GeVshiftedPt_processed_2024-07-18.root");
    presetSaveComment = "_pythiaHerwig15GeVShift";
  }

  // Find the number of comparison files from the vector size
  const int nComparisonFiles = fileName.size();
  
  // Open the files and check that they exist
  TFile* inputFile[nComparisonFiles];
  EECCard* card[nComparisonFiles];
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    
    inputFile[iFile] = TFile::Open(fileName.at(iFile));
    
    if(inputFile[iFile] == NULL){
      cout << "Error! The file " << fileName.at(iFile).Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }

    card[iFile] = new EECCard(inputFile[iFile]);
  }
  
  // Check if we are using PbPb or pp data
  TString collisionSystem = card[0]->GetDataType();
  bool isPbPbData = collisionSystem.Contains("PbPb");

  // Option to define system string for the legend
  std::vector<TString> systemForLegend;

  // Automatically determine from file name which files are Pythia8 ans which are Herwig7
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    if(fileName.at(iFile).Contains("Pythia8")){
      systemForLegend.push_back("Pythia8");
    } else {
      systemForLegend.push_back("Herwig7");
    }
  }

  // If every entry in the vector is the same, only keep the first one
  if(std::count(systemForLegend.begin(), systemForLegend.end(), systemForLegend.at(0)) == systemForLegend.size()){
    systemForLegend.resize(1);
  }
  
  
  // ====================================================
  //               Binning configuration
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = card[0]->GetNCentralityBins();
  const int nJetPtBinsEEC = card[0]->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card[0]->GetNTrackPtBinsEEC();

  // For legends, collect common and not common selections
  std::vector<TString> commonLegend;
  std::vector<TString> individualLegend;
  std::vector<TString> legendComment;

  // Define how much the jet pT is shifted in each file
  std::vector<double> jetPtShift;

  if(presetComparison > 0){
    // Determine preset comparison values form the input number
    jetPtShift.push_back(0); // The first file in the list defines the unshifted bin
    jetPtShift.push_back(15 - (10 * (presetComparison % 2) ) ); // Odd numbers have 5 GeV shift and even 15 GeV shift

    // Both Pythia and Herwig are drawn in the same figure for high preset comparison numbers
    if(presetComparison > 12){
      jetPtShift.push_back(0); // The first file in the list defines the unshifted bin
      jetPtShift.push_back(15 - (10 * (presetComparison % 2) ) ); // Odd numbers have 5 GeV shift and even 15 GeV shift
    }

  } else {
    // Manually defined shifts in the jet pT
    jetPtShift.push_back(0); // The first file in the list defines the unshifted bin
    jetPtShift.push_back(0); // The other files can contain results shifted with respect to the first file
  }  

  legendComment.push_back("(pp, no E-loss)");
  legendComment.push_back("(PbPb, before E-loss)");

  // Check that each file has a shift definition
  if(jetPtShift.size() < nComparisonFiles){
    cout << "ERROR! You have not defined a jet pT shift for every file!" << endl;
    cout << "Cannot run the code. Please define jet pT shifts for each compared file." << endl;
    return;
  }
  
  // Select explicitly which bins from the files are compared:
  std::vector<std::pair<double,double>> comparedCentralityBin;
  comparedCentralityBin.push_back(std::make_pair(0,100));

  std::vector<std::pair<double,double>> comparedJetPtBin;
  comparedJetPtBin.push_back(std::make_pair(120,140));
  //comparedJetPtBin.push_back(std::make_pair(140,160));
  //comparedJetPtBin.push_back(std::make_pair(160,180));
  //comparedJetPtBin.push_back(std::make_pair(180,200));

  std::vector<double> comparedTrackPtBin;
  comparedTrackPtBin.push_back(1);

  // Determine the used energy weight from the card
  double energyWeight = card[0]->GetWeightExponent();

  // Flag for drawing plots only with ratios
  bool drawDistributionsAndRatios = true;
  bool drawOnlyRatios = false;

  // Flag for drawing the jet shape scale
  bool drawJetShapeCumulantScale = false;
  std::vector<int> drawnCumulantIndex;  
  drawnCumulantIndex.push_back(kCorrelator);
  drawnCumulantIndex.push_back(kAdvancedCorrelator);

  if(drawnCumulantIndex.size() > 1 &&  fileName.size() > 2){
    cout << "ERROR! Please only draw one cumulant style if you add more than two files!" << endl;
    cout << "Now, the execution of the code will stop. I cannot let you make crowded figures." << endl;
    return;
  }

  // Create an object to easilty manipulate histograms
  AlgorithmLibrary* optimusPrimeTheTransformer = new AlgorithmLibrary();

  // ====================================================
  //                    Configuration
  // ====================================================
  
  // Figure saving
  const bool saveFigures = true;  // Save figures
  const char* figureFormat = "pdf"; // Format given for the figures

  TString saveComment;
  if(presetComparison > 0){
    saveComment = presetSaveComment; // Read the predefined value for the save comment
  } else {
    saveComment = "_pythiaHerwig5GeVShift"; // Manually defined save comment
  }

  // Add the energy weight string to the save comment
  if(energyWeight == 1){
    saveComment.Append("_nominalEnergyWeight");
  } else {
    saveComment.Append("_energyWeightSquared");
  }

  // Drawing configuration
  std::pair<double, double> ratioZoom = std::make_pair(0.75, 1.25); // For only ratio: 0.86, 1.72
  
  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager* histograms[nComparisonFiles];
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    histograms[iFile] = new EECHistogramManager(inputFile[iFile], card[iFile]);

    // Choose the energy-energy correlator types to load
    histograms[iFile]->SetLoadEnergyEnergyCorrelators(true);

    // Choose the bin ranges
    histograms[iFile]->SetCentralityBinRange(0, card[iFile]->GetNCentralityBins() - 1);
    histograms[iFile]->SetJetPtBinRangeEEC(0, card[iFile]->GetNJetPtBinsEEC() - 1);
    histograms[iFile]->SetTrackPtBinRangeEEC(0, card[iFile]->GetNTrackPtBinsEEC() - 1);

    // Load the histograms from the file
    histograms[iFile]->LoadProcessedHistograms();
  }

  // Energy-energy correlator, jet shape, and cumulant histograms
  TH1D* hEnergyEnergyCorrelator[nComparisonFiles][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hEnergyEnergyCorrelatorRatio[nComparisonFiles][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hJetShape[nComparisonFiles][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hCumulant[knCumulants][nComparisonFiles][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hCumulantScaledDistribution[knCumulants][nComparisonFiles][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hCumulantScaledRatio[knCumulants][nComparisonFiles][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];

  // Jet pT spectra for determining mean jet pT in each analysis bin
  TH1D* hJetPt[nComparisonFiles][nCentralityBins];
  double meanJetPt[nComparisonFiles][nCentralityBins][nJetPtBinsEEC];
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for (int iFile = 0; iFile < nComparisonFiles; iFile++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      hJetPt[iFile][iCentrality] = NULL;
      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
        meanJetPt[iFile][iCentrality][iJetPt] = 0;
        for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
          hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt] = NULL;
          hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt] = NULL;
          hJetShape[iFile][iCentrality][iJetPt][iTrackPt] = NULL;
          for(int iCumulant = 0; iCumulant < knCumulants; iCumulant++){
            hCumulant[iCumulant][iFile][iCentrality][iJetPt][iTrackPt] = NULL;
            hCumulantScaledDistribution[iCumulant][iFile][iCentrality][iJetPt][iTrackPt] = NULL;
            hCumulantScaledRatio[iCumulant][iFile][iCentrality][iJetPt][iTrackPt] = NULL;
          }
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // File loop
  
  // Helper histograms
  std::pair<double, double> drawingRange = std::make_pair(lowDrawRange, 0.39);
  double epsilon = 0.0001;
  int lowNormalizationBin, highNormalizationBin;
  int iCentrality, iCentralityMatched;
  int iTrackPt, iTrackPtMatched;
  int iJetPt, iJetPtMatched;
  TH1D* helperHistogram;

  // If drawing range is increased from original, add this to the saved name
  if(drawingRange.first < 0.008){
    saveComment.Append("_includeLowAngles");
  }

  // Get the histograms from the histogram manager and normalize the signal histograms to one
  int referenceFile = 0;
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){

    // Use even iFile indices as reference files with respect to which the ratios are taken
    if(iFile%2 == 0) referenceFile = iFile;

    for(auto centralityBin : comparedCentralityBin){
      iCentrality = isPbPbData ? card[0]->FindBinIndexCentrality(centralityBin) : 0;
      iCentralityMatched = isPbPbData ? card[iFile]->FindBinIndexCentrality(centralityBin) : 0;

      // Read the jet pT histograms
      hJetPt[iFile][iCentrality] = histograms[iFile]->GetHistogramJetPt(iCentralityMatched);

      for(auto jetPtBin : comparedJetPtBin){
        iJetPt = card[0]->FindBinIndexJetPtEEC(jetPtBin);
        iJetPtMatched = card[iFile]->FindBinIndexJetPtEEC(jetPtBin.first + jetPtShift.at(iFile), jetPtBin.second + jetPtShift.at(iFile));

        // Set the range of x-axis to correspond to the current bin range
        hJetPt[iFile][iCentrality]->GetXaxis()->SetRangeUser(jetPtBin.first + jetPtShift.at(iFile) + epsilon, jetPtBin.second + jetPtShift.at(iFile) - epsilon);

        // Get the mean jet pT in the defined range from the histogram
        meanJetPt[iFile][iCentrality][iJetPt] = hJetPt[iFile][iCentrality]->GetMean();

        for(auto trackPtBin : comparedTrackPtBin){
          iTrackPt = card[0]->GetBinIndexTrackPtEEC(trackPtBin);
          iTrackPtMatched = card[iFile]->GetBinIndexTrackPtEEC(trackPtBin);
      
          // Read the energy-energy correlator histograms
          hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt] = histograms[iFile]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentralityMatched, iJetPtMatched, iTrackPtMatched);

          // Normalize the distributions to one in the drawingRange
          lowNormalizationBin = hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(drawingRange.first + epsilon);
          highNormalizationBin = hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(drawingRange.second - epsilon);

          hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->Scale(1 / hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->Integral(lowNormalizationBin, highNormalizationBin, "width"));

          // Calculate the ratio to reference index
          hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt] = (TH1D*)hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->Clone(Form("eecRatio%d", iFile));
          hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->Divide(hEnergyEnergyCorrelator[referenceFile][iCentrality][iJetPt][iTrackPt]);

          // Read the jet shape histograms
          hJetShape[iFile][iCentrality][iJetPt][iTrackPt] = histograms[iFile]->GetHistogramJetShape(iCentralityMatched, iJetPtMatched, iTrackPtMatched);

          // Normalize the distributions such that the peak amplitude is one
          hJetShape[iFile][iCentrality][iJetPt][iTrackPt]->Scale(1 / hJetShape[iFile][iCentrality][iJetPt][iTrackPt]->GetMaximum());

          // Square the jet shape histograms
          optimusPrimeTheTransformer->SquareHistogram(hJetShape[iFile][iCentrality][iJetPt][iTrackPt]);

          // For the squared histogram, calculate the cumulant in the full DeltaR region and in the analysis region
          hCumulant[kJetShape][iFile][iCentrality][iJetPt][iTrackPt] = optimusPrimeTheTransformer->GetCumulant(hJetShape[iFile][iCentrality][iJetPt][iTrackPt], lowNormalizationBin);
          hCumulant[kJetShape][iFile][iCentrality][iJetPt][iTrackPt]->SetName(Form("jetShapeCumulant%d%d%d%d", iFile, iCentrality, iJetPt, iTrackPt));

          // Take the squared cumulant also for energy-energy correlators
          helperHistogram = (TH1D*) hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->Clone(Form("eecHelper%d%d%d%d", iFile, iCentrality, iJetPt, iTrackPt));
          helperHistogram->Scale(1 / helperHistogram->GetMaximum());
          optimusPrimeTheTransformer->SquareHistogram(helperHistogram);

          // Calculate the cumulant from the squared energy-energy correlator distribution
          hCumulant[kCorrelator][iFile][iCentrality][iJetPt][iTrackPt] = optimusPrimeTheTransformer->GetCumulant(helperHistogram, lowNormalizationBin);
          hCumulant[kCorrelator][iFile][iCentrality][iJetPt][iTrackPt]->SetName(Form("correlatorCumulant%d%d%d%d", iFile, iCentrality, iJetPt, iTrackPt));

          // For advanced cumulant scaling, use the regular cumulant histograms as a base
          hCumulant[kAdvancedCorrelator][iFile][iCentrality][iJetPt][iTrackPt] = (TH1D*) hCumulant[kCorrelator][iFile][iCentrality][iJetPt][iTrackPt]->Clone(Form("advancedCorrelatorCumulant%d%d%d%d", iFile, iCentrality, iJetPt, iTrackPt));

        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop

  } // File loop

  // Helper histogram
  TH1D* scaleHistogram;
  double distributionScale;
  double currentBinContent;
  double binCenter;
  double energyLoss;
  double advancedScaleFactor;

  // Scale the ratio histogram with the ratio of cumulants of squared jet shapes
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){

    // Use even iFile indices as reference files with respect to which the ratios are taken
    if(iFile%2 == 0) referenceFile = iFile;

    for(auto centralityBin : comparedCentralityBin){
      iCentrality = isPbPbData ? card[0]->FindBinIndexCentrality(centralityBin) : 0;
      for(auto jetPtBin : comparedJetPtBin){
        iJetPt = card[0]->FindBinIndexJetPtEEC(jetPtBin);
        energyLoss = meanJetPt[iFile][iCentrality][iJetPt] - meanJetPt[0][iCentrality][iJetPt];
        for(auto trackPtBin : comparedTrackPtBin){
          iTrackPt = card[0]->GetBinIndexTrackPtEEC(trackPtBin);

          for(int iCumulant : drawnCumulantIndex){

            advancedScaleFactor = 1;

            hCumulantScaledDistribution[iCumulant][iFile][iCentrality][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->Clone(Form("%sCumulantScaledDistribution%d%d%d%d", cumulantName[iCumulant].Data(), iFile, iCentrality, iJetPt, iTrackPt));
            scaleHistogram = (TH1D*) hCumulant[iCumulant][iFile][iCentrality][iJetPt][iTrackPt]->Clone(Form("scaleFor%sCumulant%d%d%d%d", cumulantName[iCumulant].Data(), iFile, iCentrality, iJetPt, iTrackPt));
            scaleHistogram->Divide(hCumulant[iCumulant][referenceFile][iCentrality][iJetPt][iTrackPt]);

            if(iCumulant == kAdvancedCorrelator) {
              optimusPrimeTheTransformer->ExponentiateHistogram(scaleHistogram, 2.0/3.0);
              advancedScaleFactor = 2.0/3.0;
            }

            // For the scale histogram, add the energy loss term to the cumulant ratio
            for(int iBin = 1; iBin <= scaleHistogram->GetNbinsX(); iBin++){
              currentBinContent = scaleHistogram->GetBinContent(iBin);
              binCenter = scaleHistogram->GetBinCenter(iBin);
              scaleHistogram->SetBinContent(iBin, currentBinContent - (energyLoss / meanJetPt[iFile][iCentrality][iJetPt] * (1 - advancedScaleFactor * binCenter)));
            }

            hCumulantScaledDistribution[iCumulant][iFile][iCentrality][iJetPt][iTrackPt]->Divide(scaleHistogram);
            distributionScale = hCumulantScaledDistribution[iCumulant][iFile][iCentrality][iJetPt][iTrackPt]->Integral(lowNormalizationBin, highNormalizationBin, "width");
            hCumulantScaledDistribution[iCumulant][iFile][iCentrality][iJetPt][iTrackPt]->Scale(1.0 / distributionScale);

            hCumulantScaledRatio[iCumulant][iFile][iCentrality][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->Clone(Form("%sCumulantScaledRatio%d%d%d%d", cumulantName[iCumulant].Data(), iFile, iCentrality, iJetPt, iTrackPt));
            hCumulantScaledRatio[iCumulant][iFile][iCentrality][iJetPt][iTrackPt]->Divide(scaleHistogram);
            hCumulantScaledRatio[iCumulant][iFile][iCentrality][iJetPt][iTrackPt]->Scale(1.0 / distributionScale);

          } // Variant loop for cumulants

        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // File loop
  
  
  // ==========================================================================
  //                    All the ratios in the same figure
  // ==========================================================================
  
  JDrawer* drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.57);
  drawer->SetTitleOffsetX(1.0);
  
  drawer->SetLogX(true);

  double legendXadder = 0;
  if(energyWeight == 1 && drawingRange.first < 0.008){
    legendXadder = 0.1;
  }

  TString compactCentralityString, compactTrackPtString, compactJetPtString;
  TString legendString;
  TString scaleString;
  int markerStyle[] = {kOpenSquare, kOpenCircle, kOpenCross, kOpenDoubleDiamond, kOpenDiamond, kOpenStar};
  //int markerStyle[] = {53, 54, 57, 65, 56, 58};
  int color[] = {kBlack, kRed, kBlue, kGreen+3, kMagenta, kCyan};
  int cumulantMarkerStyle[] = {kOpenDiamond, kOpenDiamond, kOpenCrossX, kOpenCrossX, kOpenDoubleDiamond, kOpenDoubleDiamond};
  //int cumulantMarkerStyle[] = {56, 56, 67, 67, 65, 65};
  int cumulantColor[] = {kMagenta, kMagenta, kCyan, kCyan, kGreen+3, kGreen+3};
  int styleIndex;
  TLegend* legend;
  TLegend* ratioLegend[4];
  TLatex* missingTextFiller = new TLatex();
  TCanvas* thisCanvas;
  TPad* overlay;
  TLine* lineDrawer = new TLine();
  lineDrawer->SetLineStyle(2);
  std::pair<double, double> histogramMinMax;

  if(drawDistributionsAndRatios){

    for(auto centralityBin : comparedCentralityBin){
      iCentrality = isPbPbData ? card[0]->FindBinIndexCentrality(centralityBin) : 0;

      if(isPbPbData){
        compactCentralityString = Form("_C=%.0f-%.0f", centralityBin.first, centralityBin.second);
      } else {
      compactCentralityString = "_pp";  
      }

      for(auto jetPtBin : comparedJetPtBin){
        iJetPt = card[0]->FindBinIndexJetPtEEC(jetPtBin);
        compactJetPtString = Form("_J=%.0f-%.0f", jetPtBin.first, jetPtBin.second);
        for(auto trackPtBin : comparedTrackPtBin){
          iTrackPt = card[0]->GetBinIndexTrackPtEEC(trackPtBin);
          compactTrackPtString = Form("_T>%.1f",trackPtBin);
          compactTrackPtString.ReplaceAll(".","v");

          // Create a new canvas for the plot
          drawer->CreateSplitCanvas();
          
          // Logarithmic EEC axis
          drawer->SetLogY(true);

          legend = new TLegend(0.16+legendXadder, 0.04, 0.43+legendXadder, 0.52);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);

          // If we have manually determined system string, use it. Otherwise read the system from card
          if(systemForLegend.size() == 1){
            legend->AddEntry((TObject*) 0, Form("%s 5.02 TeV", systemForLegend.at(0).Data()), "");
          } else if (systemForLegend.size() > 1){
            legend->AddEntry((TObject*) 0, "pp simulation 5.02 TeV", "");
          } else {
            legend->AddEntry((TObject*) 0, Form("%s 5.02 TeV",card[0]->GetAlternativeDataType(false).Data()), "");
          } 

          // Add centrality and charged particle pT information
          if(isPbPbData){
            legend->AddEntry((TObject*) 0, Form("Cent: %.0f-%.0f%%", centralityBin.first, centralityBin.second), "");
          }
          legend->AddEntry((TObject*) 0, Form("p_{T}^{ch} > %.1f GeV, n=%.0f", trackPtBin, energyWeight), "");

          // Set drawing style for all histograms
          for(int iFile = 0; iFile < nComparisonFiles; iFile++){
            hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iFile]);
            hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iFile]);
            hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iFile]);
            //hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->SetLineWidth(2);
            hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iFile]);
            hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iFile]);
            hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iFile]);
            //hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetLineWidth(2);

            // Different styles for cumulant depending if jet shape cumulant are drawn or not
            if(drawnCumulantIndex.size() > 1){
              // Several cumulants from one shifted system are drawn to the same sigure

              styleIndex = 0;
              for(int iCumulant : drawnCumulantIndex){

                hCumulantScaledDistribution[iCumulant][iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(cumulantMarkerStyle[styleIndex]);
                hCumulantScaledDistribution[iCumulant][iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(cumulantColor[styleIndex]);
                hCumulantScaledDistribution[iCumulant][iFile][iCentrality][iJetPt][iTrackPt]->SetLineColor(cumulantColor[styleIndex]);
                //hCumulantScaledDistribution[iCumulant][iFile][iCentrality][iJetPt][iTrackPt]->SetLineWidth(2);

                hCumulantScaledRatio[iCumulant][iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(cumulantMarkerStyle[styleIndex]);
                hCumulantScaledRatio[iCumulant][iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(cumulantColor[styleIndex]);
                hCumulantScaledRatio[iCumulant][iFile][iCentrality][iJetPt][iTrackPt]->SetLineColor(cumulantColor[styleIndex]);
                //hCumulantScaledRatio[iCumulant][iFile][iCentrality][iJetPt][iTrackPt]->SetLineWidth(2);

                styleIndex = styleIndex + 2;

              }

            } else {
              // Single cumulant for several files are drawn to the same figure

              hCumulantScaledDistribution[drawnCumulantIndex.at(0)][iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(cumulantMarkerStyle[iFile]);
              hCumulantScaledDistribution[drawnCumulantIndex.at(0)][iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(cumulantColor[iFile]);
              hCumulantScaledDistribution[drawnCumulantIndex.at(0)][iFile][iCentrality][iJetPt][iTrackPt]->SetLineColor(cumulantColor[iFile]);

              hCumulantScaledRatio[drawnCumulantIndex.at(0)][iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(cumulantMarkerStyle[iFile]);
              hCumulantScaledRatio[drawnCumulantIndex.at(0)][iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(cumulantColor[iFile]);
              hCumulantScaledRatio[drawnCumulantIndex.at(0)][iFile][iCentrality][iJetPt][iTrackPt]->SetLineColor(cumulantColor[iFile]);

            }
          } // Setting drawing style for histograms

          // Scale some histograms for better visual clarity, and find minimum and maximum values from all the histograms
          histogramMinMax = std::make_pair(10e10, -10e10);
          distributionScale = 1;
          for(int iFile = 0; iFile < nComparisonFiles; iFile++){

            hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->Scale(distributionScale);
            for(int iCumulant : drawnCumulantIndex){
              hCumulantScaledDistribution[iCumulant][iFile][iCentrality][iJetPt][iTrackPt]->Scale(distributionScale);
            }
  
            if(iFile%2 == 1) distributionScale *=2;

            histogramMinMax = optimusPrimeTheTransformer->FindHistogramMinMax(hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt], histogramMinMax, drawingRange);

          }

          // Add some margin to minimum and maximum values
          histogramMinMax.first = histogramMinMax.first / 2.0;
          histogramMinMax.second = histogramMinMax.second * 1.5;

          // Set the axis drawing ranges
          hEnergyEnergyCorrelator[0][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
          hEnergyEnergyCorrelator[0][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(histogramMinMax.first, histogramMinMax.second);
          
          // Draw the histograms to the upper canvas
          drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelator[0][iCentrality][iJetPt][iTrackPt], "#Deltar", "EEC", " ");

          for(int iFile = 1; iFile < nComparisonFiles; iFile++){

            // Draw the other energy-energy correlator distirbutions
            hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->Draw("same");

            // Draw also the cumulant scaled distributions to the same plot
            if(iFile%2 == 1){
              for(int iCumulant : drawnCumulantIndex){
                hCumulantScaledDistribution[iCumulant][iFile][iCentrality][iJetPt][iTrackPt]->Draw("same");
              }
            }
          }

          distributionScale = 1;
          for(int iFile = 0; iFile < nComparisonFiles; iFile++){
            scaleString = "";
            if(iFile > 1){
              scaleString = Form(" #times %.0f", distributionScale);
            }
            if(systemForLegend.size() > 1){
              legend->AddEntry(hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt], Form("%s, %.0f < p_{T,jet} < %.0f GeV%s", systemForLegend.at(iFile).Data(), jetPtBin.first + jetPtShift.at(iFile), jetPtBin.second + jetPtShift.at(iFile), scaleString.Data()), "p");

              // If jet shape cumulant scales are not drawn, add cumulants to main legend
              if(iFile%2 == 1){
                for(int iCumulant : drawnCumulantIndex){
                  legend->AddEntry(hCumulantScaledDistribution[iCumulant][iFile][iCentrality][iJetPt][iTrackPt], Form("%s, %.0f < p_{T,jet} < %.0f GeV, %s%s", systemForLegend.at(iFile).Data(), jetPtBin.first + jetPtShift.at(iFile), jetPtBin.second + jetPtShift.at(iFile), legendText[iCumulant].Data(), scaleString.Data()), "p");
                }
                distributionScale *= 2;
              }
            } else {
              legend->AddEntry(hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt], Form("%.0f < jet p_{T} < %.0f GeV %s", jetPtBin.first + jetPtShift.at(iFile), jetPtBin.second + jetPtShift.at(iFile), legendComment.at(iFile).Data()), "p");

              if(iFile%2 == 1){
                for(int iCumulant : drawnCumulantIndex){
                  legend->AddEntry(hCumulantScaledDistribution[iCumulant][iFile][iCentrality][iJetPt][iTrackPt], Form("%.0f < p_{T,jet} < %.0f GeV, %s%s", jetPtBin.first + jetPtShift.at(iFile), jetPtBin.second + jetPtShift.at(iFile), legendText[iCumulant].Data(), scaleString.Data()), "p");
                } // Cumulant style loop
              } // If for odd file indices
            } // If for length of systemForLegend array
          } // File loop
  
          // Draw the legends to the upper pad
          legend->Draw();
          
          // Linear scale for the ratio
          drawer->SetLogY(false);
          
          // Set the drawing ranges
          hEnergyEnergyCorrelatorRatio[1][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);          
          hEnergyEnergyCorrelatorRatio[1][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

          // Draw the histograms
          drawer->SetGridY(true);
          drawer->DrawHistogramToLowerPad(hEnergyEnergyCorrelatorRatio[1][iCentrality][iJetPt][iTrackPt], "#Deltar", Form("#frac{%.0f < p_{T,jet} < %.0f GeV}{%.0f < p_{T,jet} < %.0f GeV}", jetPtBin.first + jetPtShift.at(1), jetPtBin.second + jetPtShift.at(1), jetPtBin.first, jetPtBin.second), " ");

          for(int iFile = 3; iFile < nComparisonFiles; iFile = iFile+2){
            hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->Draw("same");
          }
          for(int iFile = 1; iFile < nComparisonFiles; iFile = iFile+2){
            for(int iCumulant : drawnCumulantIndex){
              hCumulantScaledRatio[iCumulant][iFile][iCentrality][iJetPt][iTrackPt]->Draw("same");
            }
          }

          drawer->SetGridY(false);

          // The eV in GeV in the lower pad title is not drawn correctly because it gets behing the upper pad Add it manually:
          thisCanvas = drawer->GetCanvas();
          thisCanvas->cd();

          overlay = new TPad("overlay", "Overlay", 0, 0, 1, 1);
          overlay->SetFillStyle(4000); // Transparent fill
          overlay->Draw();
          overlay->cd();

          missingTextFiller->SetTextFont(42);
          missingTextFiller->SetTextSize(0.038);
          missingTextFiller->SetTextAngle(90);
          missingTextFiller->DrawLatexNDC(0.083, 0.3827, "eV");
          missingTextFiller->DrawLatexNDC(0.026, 0.3827, "eV");
          
          // Save the figures to a file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/eecCustomBinComparison%s%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
          } // Save figures
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Drawing distributions together with ratios

  // Draw only the ratio plots with all jet pT bins drawn in the same plot
  if(drawOnlyRatios){

    drawer->SetRelativeCanvasSize(1,1.5);
    drawer->SetLeftMargin(0.1);

    int ratioMarkerStyle[] = {kFullSquare, kFullCircle, kFullCross, kFullCrossX, kOpenSquare, kOpenCircle, kOpenCross, kOpenCrossX};
    int ratioColor[] = {kBlack, kRed, kBlue, kGreen+3, kBlack, kRed, kBlue, kGreen+3};
    int ratioIndex;

    // Here probably need to setup new settings for JDrawer

    for(auto centralityBin : comparedCentralityBin){
      iCentrality = isPbPbData ? card[0]->FindBinIndexCentrality(centralityBin) : 0;

      if(isPbPbData){
        compactCentralityString = Form("_C=%.0f-%.0f", centralityBin.first, centralityBin.second);
      } else {
        compactCentralityString = "_pp";  
      }

      for(auto trackPtBin : comparedTrackPtBin){
        iTrackPt = card[0]->GetBinIndexTrackPtEEC(trackPtBin);
        compactTrackPtString = Form("_T>%.1f",trackPtBin);
        compactTrackPtString.ReplaceAll(".","v");
          
        // Use linear axis scale since we are only drawing ratios
        drawer->SetLogY(false);

        // Create two distinct legends
        ratioLegend[0] = new TLegend(0.122, 0.4, 0.43, 0.45);
        ratioLegend[1] = new TLegend(0.1, 0.3, 0.9, 0.4); 
        ratioLegend[2] = new TLegend(0.122, 0.8, 0.43, 0.9);
        ratioLegend[3] = new TLegend(0.1, 0.7, 0.9, 0.8);
        

        for(int iLegend = 0; iLegend < 4; iLegend++){
          ratioLegend[iLegend]->SetFillStyle(0);ratioLegend[iLegend]->SetBorderSize(0);
          ratioLegend[iLegend]->SetTextSize(0.037);ratioLegend[iLegend]->SetTextFont(62);
        }

        ratioLegend[2]->AddEntry((TObject*) 0, Form("pp simulation 5.02 TeV, p_{T}^{ch} > %.1f GeV, n=%.0f", trackPtBin, energyWeight), "");
        ratioLegend[2]->AddEntry((TObject*) 0, "E2C Cumulant Scale #times 1.4", "");
        ratioLegend[0]->AddEntry((TObject*) 0, "E2C Cumulant Scale", "");

        ratioLegend[1]->SetNColumns(2);
        ratioLegend[3]->SetNColumns(2);

        // Set drawing style for correlator cumulant scaled ratios
        ratioIndex = 0;
        for(int iFile = 1; iFile < nComparisonFiles; iFile = iFile + 2){
          for(auto jetPtBin : comparedJetPtBin){
            iJetPt = card[0]->FindBinIndexJetPtEEC(jetPtBin);

            hCumulantScaledRatio[drawnCumulantIndex.at(0)][iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(ratioMarkerStyle[ratioIndex]);
            hCumulantScaledRatio[drawnCumulantIndex.at(0)][iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(ratioColor[ratioIndex]);
            hCumulantScaledRatio[drawnCumulantIndex.at(0)][iFile][iCentrality][iJetPt][iTrackPt]->SetLineColor(ratioColor[ratioIndex++]);
          } // Jet pT loop
        } // Setting drawing style for histograms

        // Draw all the hisrograms to the canvas
        ratioIndex = 0;
        distributionScale = 1;
        scaleString = "";
        for(int iFile = 1; iFile < nComparisonFiles; iFile = iFile + 2){
          distributionScale = 1 + (iFile-1)*0.2;
          for(auto jetPtBin : comparedJetPtBin){
            iJetPt = card[0]->FindBinIndexJetPtEEC(jetPtBin);

            // Scale distributions from different files to keep the drawing cleaner
            hCumulantScaledRatio[drawnCumulantIndex.at(0)][iFile][iCentrality][iJetPt][iTrackPt]->Scale(distributionScale);

            // Draw all the ratios to the same canvas
            if(ratioIndex++ == 0){
              // For the first histogram, set correct drawing range and draw it to the canvas
              hCumulantScaledRatio[drawnCumulantIndex.at(0)][iFile][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
              hCumulantScaledRatio[drawnCumulantIndex.at(0)][iFile][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
              drawer->DrawHistogram(hCumulantScaledRatio[drawnCumulantIndex.at(0)][iFile][iCentrality][iJetPt][iTrackPt], "#Deltar", Form("Ratio to %.0f GeV lower p_{T,jet} bin", jetPtShift.at(1)), " ");
            } else {
              hCumulantScaledRatio[drawnCumulantIndex.at(0)][iFile][iCentrality][iJetPt][iTrackPt]->Draw("same");
            }

            // Add things to legend
            if(iFile > 1){
              scaleString = Form(" #times %.1f", distributionScale);
            }

            if(systemForLegend.size() > 1){
              ratioLegend[iFile]->AddEntry(hCumulantScaledRatio[drawnCumulantIndex.at(0)][iFile][iCentrality][iJetPt][iTrackPt], Form("%s, %.0f < p_{T,jet} < %.0f GeV", systemForLegend.at(iFile).Data(), jetPtBin.first + jetPtShift.at(iFile), jetPtBin.second + jetPtShift.at(iFile)), "p");
            } else {
              ratioLegend[iFile]->AddEntry(hCumulantScaledRatio[drawnCumulantIndex.at(0)][iFile][iCentrality][iJetPt][iTrackPt], Form("%.0f < jet p_{T} < %.0f GeV %s", jetPtBin.first + jetPtShift.at(iFile), jetPtBin.second + jetPtShift.at(iFile), legendComment.at(iFile).Data()), "p");
            }

            // Draw lines to where different distributions are scaled
            if(jetPtBin == comparedJetPtBin.at(0)){
              lineDrawer->DrawLine(drawingRange.first, distributionScale, drawingRange.second, distributionScale);
            }
          } // Jet pT loop
        } // File loop
  
        // Draw the legends
        for(int iLegend = 0; iLegend < 4; iLegend++){
          ratioLegend[iLegend]->Draw();
        }
          
        // Save the figures to a file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/eecCustomBinRatioComparison%s%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
        } // Saving the file
      } // Track pT loop
    } // Centrality loop
  } // Drawing only ratios
}
