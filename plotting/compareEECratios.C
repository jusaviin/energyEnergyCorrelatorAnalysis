#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"
#include "../src/EECHistograms.h"
#include "AlgorithmLibrary.h"
#include "SystematicUncertaintyOrganizer.h"

/*
 * Macro for comparing final energy-energy correlators results
 */
void compareEECratios(){
  
  // Files for comparison
  std::vector<std::pair<TString,TString>> fileName;
  //fileName.push_back(std::make_pair("data/pPb/ppData_pfJets_eschemeAxis_nominalEnergyWeight_jetEtaCMcut_noBackgroundSubtraction_processed_2025-06-30.root", "data/pPb/pPb_5TeV_pToMinusEta_pfJets_eschemeAxis_nominalEnergyWeight_minimumBias_jetEtaMCcut_noSubtraction_processed_2025-06-05.root"));
  fileName.push_back(std::make_pair("data/pPb/ppData_pfJets_eschemeAxis_nominalEnergyWeight_jetEtaCMcut_mixedConeSubtractedHFShift25_processed_2025-06-30.root", "data/pPb/pPb_5TeV_pToMinusEta_pfJets_eschemeAxis_nominalEnergyWeight_minimumBias_jetEtaMCcut_mixedConeHFshift28_processed_2025-06-30.root"));
  fileName.push_back(std::make_pair("data/pPb/ppData_pfJets_eschemeAxis_nominalEnergyWeight_jetEtaCMcut_perpendicularConeSubtracted_processed_2025-06-30.root", "data/pPb/pPb_5TeV_pToMinusEta_pfJets_eschemeAxis_nominalEnergyWeight_minimumBias_jetEtaMCcut_perpendicularConeSubtracted_processed_2025-06-05.root"));

  const int nComparisonFiles = fileName.size();

  std::vector<TString> fileDescription;
  //fileDescription.push_back("pPb/pp 5 TeV, no bg subtraction");
  fileDescription.push_back("pPb/pp 5 TeV, mixed cone sub");
  fileDescription.push_back("pPb/pp 5 TeV, perp cone sub");

  // Check that a description exists for each file
  if(fileDescription.size() < fileName.size()){
    cout << "ERROR! Not enough file descriptions given. Please give a description for all your files!" << endl;
    return;
  }
  
  // Open the files and check that they exist
  std::vector<std::pair<TFile*,TFile*>> inputFile;
  std::vector<std::pair<EECCard*, EECCard*>> card;
  TFile* openedFile1;
  TFile* openedFile2;
  for(auto theseFiles : fileName){
    
    openedFile1 = TFile::Open(theseFiles.first);
    
    if(openedFile1 == NULL){
      cout << "Error! The file " << theseFiles.first.Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }

    openedFile2 = TFile::Open(theseFiles.second);
    
    if(openedFile2 == NULL){
      cout << "Error! The file " << theseFiles.second.Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }

    inputFile.push_back(std::make_pair(openedFile1, openedFile2));
    card.push_back(std::make_pair(new EECCard(openedFile1), new EECCard(openedFile2)));

  }

  // Use the first non-pp card as the referencea
  EECCard* referenceCard = card.at(0).second;

  // ====================================================
  //               Binning configuration
  // ====================================================

  // Find the number of bins from the reference card
  const int nCentralityBins = referenceCard->GetNCentralityBins();
  const int nJetPtBinsEEC = referenceCard->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = referenceCard->GetNTrackPtBinsEEC();
  
  
  // Select explicitly which bins from the files are compared:
  std::vector<std::pair<double,double>> comparedCentralityBin;
  comparedCentralityBin.push_back(std::make_pair(0,10));
  //comparedCentralityBin.push_back(std::make_pair(10,30));
  //comparedCentralityBin.push_back(std::make_pair(30,50));
  //comparedCentralityBin.push_back(std::make_pair(50,90));
  bool individualCentrality = true; // True = make different figure for each bin. False = plot all centrality bin to the same figure.

  std::vector<std::pair<double,double>> comparedJetPtBin;
  comparedJetPtBin.push_back(std::make_pair(30,40));
  comparedJetPtBin.push_back(std::make_pair(40,50));
  comparedJetPtBin.push_back(std::make_pair(50,60));
  comparedJetPtBin.push_back(std::make_pair(60,80));
  bool individualJetPt = true; // True = make different figure for each bin. False = plot all jet pT bin to the same figure.

  std::vector<double> comparedTrackPtBin;
  comparedTrackPtBin.push_back(1.0);
  //comparedTrackPtBin.push_back(1.5);
  //comparedTrackPtBin.push_back(2.0);
  //comparedTrackPtBin.push_back(2.5);
  //comparedTrackPtBin.push_back(3.0);
  bool individualTrackPt = true; // True = make different figure for each bin. False = plot all track pT bin to the same figure.

  // Different normalization options
  enum enumNormalization{kNoNormalization, kNormalizeToPairs, kNormalizeToJets};
  const int normalizeDistributions = kNormalizeToJets;

  // Choose the type of draw energy-energy correlator
  // EECHistogramManager::kEnergyEnergyCorrelatorNormalized = Normalized energy-energy correlator
  // EECHistogramManager::kEnergyEnergyCorrelatorBackground = Estimated background
  // EECHistogramManager::kEnergyEnergyCorrelatorSignal = Background subtracted, but not unfolded energy-energy correlator
  // EECHistogramManager::kEnergyEnergyCorrelatorUnfolded = Unfolded energy-energy correlator
  // EECHistogramManager::kEnergyEnergyCorrelatorBackgroundAfterUnfolding = Estimated background after unfolding
  // EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal = Unfolded energy-energy correlator signal
  // EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels = Raw energy-energy correlator
  int drawnEnergyEnergyCorrelator = EECHistogramManager::kEnergyEnergyCorrelatorSignal;

  // Choose the pairing type if raw energy-energy correlator is drawn
  // EECHistograms::kSameJetPair;
  // EECHistograms::kSignalReflectedConePair;
  // EECHistograms::kReflectedConePair;
  // EECHistograms::kSignalMixedConePair;
  // EECHistograms::kReflectedMixedConePair;
  // EECHistograms::kMixedConePair;
  // EECHistograms::kSignalSecondMixedConePair; 
  // EECHistograms::kReflectedSecondMixedConePair; 
  // EECHistograms::kMixedMixedConePair; 
  // EECHistograms::kSecondMixedConePair;
  int iPairingType = EECHistograms::kSameJetPair;

  // No centrality selection if there is no PbPb data
  std::vector<bool> skipCentrality;
  for(auto thisCard : card){
    skipCentrality.push_back(!thisCard.second->GetDataType().Contains("PbPb"));
  }

  // If we are dealing with MC, shift the centrality by 4% as is done in order to match background energy density
  //if(card[0]->GetDataType().Contains("MC")){
  //  for(auto& centralityBin : comparedCentralityBin){
  //    centralityBin.first += 4;
  //    centralityBin.second += 4;
  //  }
  //}

  // ====================================================
  //                Drawing configuration
  // ====================================================
  
  // Figure saving
  const bool saveFigures = true;  // Save figures
  const char* saveComment = "_aliceAspectRatioJetNormalized";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures

  // Drawing configuration
  std::pair<double, double> ratioZoom = std::make_pair(0.7, 1.35);

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
  std::vector<std::pair<EECHistogramManager*, EECHistogramManager*>> histograms;
  for(int iFile = 0; iFile < inputFile.size(); iFile++){
    histograms.push_back(std::make_pair(new EECHistogramManager(inputFile.at(iFile).first, card.at(iFile).first), new EECHistogramManager(inputFile.at(iFile).second, card.at(iFile).second)));
  }

  // Energy-energy correlator histograms
  TH1D* hEnergyEnergyCorrelatorReference[nComparisonFiles][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hEnergyEnergyCorrelatorMedium[nComparisonFiles][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hEnergyEnergyCorrelatorRatio[nComparisonFiles][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        hEnergyEnergyCorrelatorReference[nComparisonFiles][nJetPtBinsEEC][nTrackPtBinsEEC] = NULL;
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          hEnergyEnergyCorrelatorMedium[iFile][iCentrality][iJetPt][iTrackPt] = NULL;
          hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt] = NULL;
        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // File loop
  
  
  // Helper histograms
  std::pair<double, double> drawingRange = std::make_pair(0.008, 0.39);
  double epsilon = 0.0001;
  int lowNormalizationBin, highNormalizationBin;
  int iCentrality, iCentralityReference;
  int iTrackPt, iTrackPtReference;
  int iJetPt, iJetPtReference;
  std::pair<double,double> referenceCentralityBin;

  // Get the histograms from the histogram manager and normalize the signal histograms to one
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    for(auto jetPtBin : comparedJetPtBin){
      for(auto trackPtBin : comparedTrackPtBin){

        // Read the reference energy-energy correlators (pp)
        iJetPt = card.at(iFile).first->FindBinIndexJetPtEEC(jetPtBin);
        iTrackPt = card.at(iFile).first->GetBinIndexTrackPtEEC(trackPtBin);
        iJetPtReference = referenceCard->FindBinIndexJetPtEEC(jetPtBin);
        iTrackPtReference = referenceCard->GetBinIndexTrackPtEEC(trackPtBin);

        // Load the selected energy-energy correlator histogram
        if(drawnEnergyEnergyCorrelator == EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels){
          hEnergyEnergyCorrelatorReference[iFile][iJetPtReference][iTrackPtReference] = histograms.at(iFile).first->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPt, iTrackPt, iPairingType);
        } else {
          hEnergyEnergyCorrelatorReference[iFile][iJetPtReference][iTrackPtReference] = histograms.at(iFile).first->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPt, iTrackPt, drawnEnergyEnergyCorrelator);
        }

        // Normalize the distributions to one in the drawingRange
        if(normalizeDistributions == kNormalizeToPairs){
          lowNormalizationBin = hEnergyEnergyCorrelatorReference[iFile][iJetPtReference][iTrackPtReference]->GetXaxis()->FindBin(drawingRange.first + epsilon);
          highNormalizationBin = hEnergyEnergyCorrelatorReference[iFile][iJetPtReference][iTrackPtReference]->GetXaxis()->FindBin(drawingRange.second - epsilon);

          hEnergyEnergyCorrelatorReference[iFile][iJetPtReference][iTrackPtReference]->Scale(1.0 / hEnergyEnergyCorrelatorReference[iFile][iJetPtReference][iTrackPtReference]->Integral(lowNormalizationBin, highNormalizationBin, "width"));
        } else if (normalizeDistributions == kNormalizeToJets){
          hEnergyEnergyCorrelatorReference[iFile][iJetPtReference][iTrackPtReference]->Scale(1.0 / histograms.at(iFile).first->GetJetPtIntegral(0, jetPtBin.first, jetPtBin.second));
        }

        for(auto centralityBin : comparedCentralityBin){

          // Read the medium modified energy-energy correlators (pPb or PbPb)

          // For MC, shift the centrality bin borders
          referenceCentralityBin = centralityBin;
          if(referenceCard->GetDataType().Contains("MC")){
            referenceCentralityBin.first += 4;
            referenceCentralityBin.second += 4;
          }
          if(card.at(iFile).second->GetDataType().Contains("MC")){
            centralityBin.first += 4;
            centralityBin.second += 4;
          }

          // Find the proper binning and express it in term of the bins in the first file
          iJetPt = card.at(iFile).second->FindBinIndexJetPtEEC(jetPtBin);
          iTrackPt = card.at(iFile).second->GetBinIndexTrackPtEEC(trackPtBin);

          if(skipCentrality.at(iFile)){
            iCentrality = 0;
            iCentralityReference = 0;
          } else {
            iCentrality = card.at(iFile).second->FindBinIndexCentrality(centralityBin);
            iCentralityReference = referenceCard->FindBinIndexCentrality(referenceCentralityBin);
          }

          // Load the selected energy-energy correlator histogram
          if(drawnEnergyEnergyCorrelator == EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels){
            hEnergyEnergyCorrelatorMedium[iFile][iCentralityReference][iJetPtReference][iTrackPtReference] = histograms.at(iFile).second->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iPairingType);
          } else {
            hEnergyEnergyCorrelatorMedium[iFile][iCentralityReference][iJetPtReference][iTrackPtReference] = histograms.at(iFile).second->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, drawnEnergyEnergyCorrelator);
          }

          // Normalize the distributions to one in the drawingRange
          if(normalizeDistributions == kNormalizeToPairs){
            lowNormalizationBin = hEnergyEnergyCorrelatorMedium[iFile][iCentralityReference][iJetPtReference][iTrackPtReference]->GetXaxis()->FindBin(drawingRange.first + epsilon);
            highNormalizationBin = hEnergyEnergyCorrelatorMedium[iFile][iCentralityReference][iJetPtReference][iTrackPtReference]->GetXaxis()->FindBin(drawingRange.second - epsilon);

            hEnergyEnergyCorrelatorMedium[iFile][iCentralityReference][iJetPtReference][iTrackPtReference]->Scale(1 / hEnergyEnergyCorrelatorMedium[iFile][iCentralityReference][iJetPtReference][iTrackPtReference]->Integral(lowNormalizationBin, highNormalizationBin, "width"));
          } else if (normalizeDistributions == kNormalizeToJets){
          hEnergyEnergyCorrelatorMedium[iFile][iCentralityReference][iJetPtReference][iTrackPtReference]->Scale(1.0 / histograms.at(iFile).second->GetJetPtIntegral(iCentrality, jetPtBin.first, jetPtBin.second));
          }

        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // File loop

  // After all the histograms have been read, calculate the ratios
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    for(auto jetPtBin : comparedJetPtBin){
      iJetPt = referenceCard->FindBinIndexJetPtEEC(jetPtBin);
      for(auto trackPtBin : comparedTrackPtBin){
        iTrackPt = referenceCard->GetBinIndexTrackPtEEC(trackPtBin);
        for(auto centralityBin : comparedCentralityBin){

          if(referenceCard->GetDataType().Contains("MC")){
            centralityBin.first += 4;
            centralityBin.second += 4;
          }

          if(skipCentrality.at(iFile)){
            iCentrality = 0;
          } else {
            iCentrality = referenceCard->FindBinIndexCentrality(centralityBin);
          }

          hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelatorMedium[iFile][iCentrality][iJetPt][iTrackPt]->Clone(Form("eecRatio%d%d%d%d", iFile, iCentrality, iJetPt, iTrackPt));
          hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->Divide(hEnergyEnergyCorrelatorReference[iFile][iJetPt][iTrackPt]);

        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // File loop

  // ==========================================================================
  //                  Draw all the ratios in the same figure
  // ==========================================================================
  
  JDrawer* drawer = new JDrawer();
  //drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,0.94);
  //drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.055);
  drawer->SetBottomMargin(0.155);
  //drawer->SetTitleOffsetY(1.7);
  //drawer->SetTitleOffsetX(1.0);
  
  drawer->SetLogX(true);

  TLine* oneLine = new TLine(drawingRange.first, 1, drawingRange.second, 1);
  oneLine->SetLineStyle(2);

  TString compactCentralityString = "";
  TString compactJetPtString = "";
  TString compactTrackPtString = "";
  TString comparedVariableString = "";
  TString ratioName = "";
  TString energyWeightString = (referenceCard->GetWeightExponent() == 1) ? "Nominal energy weight" : "Energy weight squared";
  energyWeightString = "Per jet normalized";
  TString legendString;
  int markerStyle[5] = {kFullCircle, kOpenSquare, kOpenCross, kFullStar, kFullCross};
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
  double minimumCandidate, maximumCandidate;
  TString individualLegend;
  TString yAxisName = "EEC (raw)";

  for(auto plottedBin : binningInformation){

    // Find the indices for each drawn bin
    currentCentralityIndices.clear();
    for(auto centralityBin : std::get<kCentrality>(plottedBin)){

      if(referenceCard->GetDataType().Contains("MC")){
        centralityBin.first += 4;
        centralityBin.second += 4;
      }
      if(referenceCard->GetDataType().Contains("pp") || referenceCard->GetDataType().Contains("pPb")){
        currentCentralityIndices.push_back(0);
      } else {
        currentCentralityIndices.push_back(referenceCard->FindBinIndexCentrality(centralityBin));
      }
    }

    currentJetPtIndices.clear();
    for(auto jetPtBin : std::get<kJetPt>(plottedBin)){
      currentJetPtIndices.push_back(referenceCard->FindBinIndexJetPtEEC(jetPtBin));
    }

    currentTrackPtIndices.clear();
    for(auto trackPtBin : std::get<kTrackPt>(plottedBin)){
      currentTrackPtIndices.push_back(referenceCard->GetBinIndexTrackPtEEC(trackPtBin));
    }

    // If the lenght of the indices array is larger than 1, we are plotting all these bins to the same plot
    colorWithCentrality = (currentCentralityIndices.size() > 1);
    colorWithTrackPt = (currentTrackPtIndices.size() > 1);
    colorWithJetPt = (currentJetPtIndices.size() > 1);

    // Remember the first drawn bin index
    firstCentralityBin = currentCentralityIndices.at(0);
    firstTrackPtBin = currentTrackPtIndices.at(0);
    firstJetPtBin = currentJetPtIndices.at(0);

    // Linear scale for the ratio
    drawer->SetLogY(false);

    TLegend* legend = new TLegend(0.18,0.74,0.45,0.98);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);

    TLegend* anotherLegend = new TLegend(0.18,0.14,0.45,0.44);
    anotherLegend->SetFillStyle(0);anotherLegend->SetBorderSize(0);anotherLegend->SetTextSize(0.04);anotherLegend->SetTextFont(62);

    // Add the file description to the legend if coloring with any variable
    if(colorWithCentrality || colorWithJetPt || colorWithTrackPt){
      legend->AddEntry((TObject*) 0, fileDescription.at(0), "");
    }

    // Add common legend variables and define figure naming in case figures are saved
    if(!colorWithCentrality){ 
      if(referenceCard->GetDataType().Contains("pp") || referenceCard->GetDataType().Contains("pPb")){
        compactCentralityString = "";
      } else {
        legend->AddEntry((TObject*) 0, Form("Cent: %.0f-%.0f%%", std::get<kCentrality>(plottedBin).at(0).first, std::get<kCentrality>(plottedBin).at(0).second), "");
        compactCentralityString = Form("_C=%.0f-%.0f", std::get<kCentrality>(plottedBin).at(0).first, std::get<kCentrality>(plottedBin).at(0).second);
      }
    } else {
      compactCentralityString = "";
      comparedVariableString = "_centralityComparison";
    }

    if(!colorWithJetPt) {
      legend->AddEntry((TObject*) 0, Form("%.0f < jet p_{T} < %.0f GeV", std::get<kJetPt>(plottedBin).at(0).first, std::get<kJetPt>(plottedBin).at(0).second), "");
      compactJetPtString = Form("_J=%.0f-%.0f", std::get<kJetPt>(plottedBin).at(0).first, std::get<kJetPt>(plottedBin).at(0).second);
    } else {
      compactJetPtString = "";
      comparedVariableString = "_jetPtComparison";
    }

    if(!colorWithTrackPt){ 
      legend->AddEntry((TObject*) 0, Form("p_{T}^{ch} > %.1f GeV", std::get<kTrackPt>(plottedBin).at(0)), "");
      compactTrackPtString = Form("_T>%.1f",std::get<kTrackPt>(plottedBin).at(0));
      compactTrackPtString.ReplaceAll(".","v");
    } else {
      compactTrackPtString = "";
      comparedVariableString = "_trackPtComparison";
    }

    legend->AddEntry((TObject*) 0, energyWeightString, "");

    // Set drawing style for all histograms
    colorFinder = 0;
    for(int iFile = 0; iFile < nComparisonFiles; iFile++){
      for(int iJetPt : currentJetPtIndices){
        for(int iTrackPt : currentTrackPtIndices){
          for(int iCentrality : currentCentralityIndices){
            if(colorWithCentrality) colorFinder++;
            hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iFile+colorFinder]);
            hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iFile+colorFinder]);
            hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iFile+colorFinder]);
          } // Centrality binning
          if(colorWithTrackPt) colorFinder++;
        } // Track pT binning 
        if(colorWithJetPt) colorFinder++;
      } // Jet pT binning
    } // File loop

    // Set the x- and y-axis drawing ranges
    hEnergyEnergyCorrelatorRatio[0][firstCentralityBin][firstJetPtBin][firstTrackPtBin]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
    hEnergyEnergyCorrelatorRatio[0][firstCentralityBin][firstJetPtBin][firstTrackPtBin]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
          
    // Draw the histograms to the upper canvas
    if(drawnEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels){
      yAxisName = Form("EEC %s", histograms.at(0).first->GetEnergyEnergyCorrelatorProcessSaveName(drawnEnergyEnergyCorrelator));
    }
    drawer->DrawHistogram(hEnergyEnergyCorrelatorRatio[0][firstCentralityBin][firstJetPtBin][firstTrackPtBin], "#Deltar", yAxisName, " ");

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

    // Add legends for drawn histograms
    legendCentralityIndex = 0;
    legendJetPtIndex = 0;
    legendTrackPtIndex = 0;
    for(int iFile = 0; iFile < nComparisonFiles; iFile++){
      if(nComparisonFiles > 1){
        individualLegend = fileDescription.at(iFile);
      } else {
        individualLegend = "";
      }
      for(int iJetPt : currentJetPtIndices){
        if(colorWithJetPt) individualLegend = Form(" %.0f < jet p_{T} < %.0f GeV", std::get<kJetPt>(plottedBin).at(legendJetPtIndex).first, std::get<kJetPt>(plottedBin).at(legendJetPtIndex).second);
        legendJetPtIndex++;
        for(int iTrackPt : currentTrackPtIndices){
          if(colorWithTrackPt) individualLegend = Form(" track p_{T} > %.1f GeV", std::get<kTrackPt>(plottedBin).at(legendTrackPtIndex++));
          for(int iCentrality : currentCentralityIndices){
            if(colorWithCentrality) individualLegend = Form(" Cent: %.0f-%.0f%%", std::get<kCentrality>(plottedBin).at(legendCentralityIndex).first, std::get<kCentrality>(plottedBin).at(legendCentralityIndex).second);
            legendCentralityIndex++;
            anotherLegend->AddEntry(hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt], individualLegend.Data(), "p");
          } // Centrality loop 
        } // Track pT loop
      } // Jet pT loop
    } // File loop
  
    // Draw the legend
    legend->Draw();
    anotherLegend->Draw();

    // Add a line to one
    oneLine->Draw();
          
    // Save the figures to a file
    if(saveFigures){
    gPad->GetCanvas()->SaveAs(Form("figures/eecRatioComparison%s%s%s%s%s.%s", saveComment, comparedVariableString.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
    }
  }

}
