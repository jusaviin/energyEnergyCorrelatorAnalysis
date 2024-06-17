#include "SystematicUncertaintyOrganizer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECHistogramManager.h"
#include "EECCard.h"
#include "JDrawer.h"
#include "AlgorithmLibrary.h"
#include "HybridModelHistogramManager.h"

/*
 * Determine up and down uncertainty bands for shape correlated uncertainties
 */
void calculateCorrelatedBands(TH1D* correlatedHistogram, TH1D* correlatedBandUpHistogram, TH1D* correlatedBandDownHistogram){
  double binContent, binError;
  int criticalBin = correlatedHistogram->GetXaxis()->FindBin(0.055);
  for(int iBin = 1; iBin <= criticalBin; iBin++){
    binContent = correlatedHistogram->GetBinContent(iBin);
    binError = correlatedHistogram->GetBinError(iBin);
    correlatedBandUpHistogram->SetBinContent(iBin, binContent+binError/2.0);
    correlatedBandUpHistogram->SetBinError(iBin, binError/2.0);
    correlatedBandDownHistogram->SetBinContent(iBin, binContent-binError/2.0);
    correlatedBandDownHistogram->SetBinError(iBin, binError/2.0);
  }
  for(int iBin = criticalBin+1; iBin <= correlatedHistogram->GetNbinsX(); iBin++){
    binContent = correlatedHistogram->GetBinContent(iBin);
    binError = correlatedHistogram->GetBinError(iBin);
    correlatedBandUpHistogram->SetBinContent(iBin, binContent-binError/2.0);
    correlatedBandUpHistogram->SetBinError(iBin, binError/2.0);
    correlatedBandDownHistogram->SetBinContent(iBin, binContent+binError/2.0);
    correlatedBandDownHistogram->SetBinError(iBin, binError/2.0);
  }
}

/*
 * Macro for making final result plots comparing energy-energy correlators between pp and PbPb
 */
void modelComparison(){

  enum enumDataType{kPbPb, kPp, kNDataTypes};
  enum enumSystematicUncertaintyType{kUncorrelatedUncertainty, kCorrelatedUncertainty, kCorrelatedUncertaintyShapeUp, kCorrelatedUncertaintyShapeDown, knSystematicUncertaintyTypes};

  // ============= //
  // Configuration //
  // ============= //

  // Select the weight exponent that is used for the figures
  int weightExponent = 1;

  // Check that the selected weight exponent is reasonable
  if(weightExponent < 1 || weightExponent > 2){
    cout << "ERROR! You selected " << weightExponent << " for the value for weight exponent." << endl;
    cout << "However, only 1 and 2 are implemented." << endl;
    cout << "Please select one of the implemented values." << endl;
    return;
  }
  
  const int nWeightExponents = 2;

  // Input files
  TString inputFileName[kNDataTypes][nWeightExponents];
  inputFileName[kPbPb][0] = "data/eecAnalysis_akFlowJet_nominalEnergyWeight_combinedMixedConeBackground_fixedCovarianceMatrix_processed_2024-05-28.root";
  inputFileName[kPbPb][1] = "data/eecAnalysis_akFlowJet_energyWeightSquared_combinedMixedConeBackground_fixedCovarianceMatrix_processed_2024-05-28.root";
  inputFileName[kPp][0] = "data/ppData_pfJets_wtaAxis_nominalEnergyWeight_optimizedUnfoldingBins_fixedCovarianceMatrix_jet60or80triggers_addLowPtBins_processed_2024-04-18.root";
  inputFileName[kPp][1] = "data/ppData_pfJets_wtaAxis_energyWeightSquared_optimizedUnfoldingBins_fixedCovarianceMatrix_jet60or80triggers_unfoldingWithCovariance_processed_2024-01-23.root";
  TString uncertaintyFileName[kNDataTypes][nWeightExponents];
  uncertaintyFileName[kPbPb][0] = "systematicUncertainties/systematicUncertainties_PbPb_nominalEnergyWeight_combinedMixedConeBackground_noMCnonClosure_2024-05-28.root";
  uncertaintyFileName[kPbPb][1] = "systematicUncertainties/systematicUncertainties_PbPb_energyWeightSquared_combinedMixedConeBackground_noMCnonClosure_2024-05-28.root";
  uncertaintyFileName[kPp][0] = "systematicUncertainties/systematicUncertainties_pp_nominalEnergyWeight_noMCnonClosure_2024-05-02.root";
  uncertaintyFileName[kPp][1] = "systematicUncertainties/systematicUncertainties_pp_energyWeightSquared_noMCnonClosure_2024-05-02.root";

  // String pointing to the folder where the hybrid model predictions are located
  TString hybridModelFolder = "theoryComparison/hybridModel/data";
  
  TFile* inputFile[kNDataTypes][nWeightExponents];
  TFile* uncertaintyFile[kNDataTypes][nWeightExponents];
  EECCard* card[kNDataTypes][nWeightExponents];
  EECCard* uncertaintyCard[kNDataTypes][nWeightExponents];
  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    for(int iFile = 0; iFile < kNDataTypes; iFile++){

      // Load the input file
      inputFile[iFile][iWeightExponent] = TFile::Open(inputFileName[iFile][iWeightExponent]);

      // Check that the input file exists
      if(inputFile[iFile][iWeightExponent] == NULL){
        cout << "Error! The file " << inputFileName[iFile][iWeightExponent].Data() << " does not exist!" << endl;
        cout << "Maybe you forgot the data/ folder path?" << endl;
        cout << "Will not execute the code" << endl;
        return;
      }
    
      // Load the card from the file and read the collision system
      card[iFile][iWeightExponent] = new EECCard(inputFile[iFile][iWeightExponent]);

      // Load the uncertainty file
      uncertaintyFile[iFile][iWeightExponent] = TFile::Open(uncertaintyFileName[iFile][iWeightExponent]);

      // Check that the uncertianty file exists
      if(uncertaintyFile[iFile][iWeightExponent] == NULL){
        cout << "Error! The file " << uncertaintyFileName[iFile][iWeightExponent].Data() << " does not exist!" << endl;
        cout << "Maybe you forgot the systematicUncertainties/ folder path?" << endl;
        cout << "Will not execute the code" << endl;
        return;
      }
    
      // Load the card from the file and read the collision system
      uncertaintyCard[iFile][iWeightExponent] = new EECCard(uncertaintyFile[iFile][iWeightExponent]);
    }  // File loop

  } // Weight exponent loop

  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Find the number of bins from the PbPb card
  const int nCentralityBins = card[kPbPb][0]->GetNCentralityBins();
  const int nJetPtBinsEEC = card[kPbPb][0]->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card[kPbPb][0]->GetNTrackPtBinsEEC();
  
  // Select explicitly which bins from the files are compared:
  std::vector<std::pair<double,double>> drawnCentralityBin;
  drawnCentralityBin.push_back(std::make_pair(0,10));
  drawnCentralityBin.push_back(std::make_pair(10,30));
  drawnCentralityBin.push_back(std::make_pair(30,50));
  drawnCentralityBin.push_back(std::make_pair(50,90));
  
  std::vector<std::pair<double,double>> drawnJetPtBin;
  drawnJetPtBin.push_back(std::make_pair(120,140));
  drawnJetPtBin.push_back(std::make_pair(140,160));
  drawnJetPtBin.push_back(std::make_pair(160,180));
  drawnJetPtBin.push_back(std::make_pair(180,200));

  std::vector<double> drawnTrackPtBin;
  drawnTrackPtBin.push_back(1.0);
  //drawnTrackPtBin.push_back(1.5);
  drawnTrackPtBin.push_back(2.0);
  //drawnTrackPtBin.push_back(2.5);
  //drawnTrackPtBin.push_back(3.0);

  // Choose which plots to draw
  bool drawDistributionDataToTheoryComparison = true;
  bool drawRatioDataToTheoryComparison = true;

  bool drawVerticalLines = false; // Draw illustrative vertical lines

  // Save the final plots
  const bool saveFigures = true;
  TString energyWeightString[nWeightExponents] = {"_nominalEnergyWeight", "_energyWeightSquared"};
  TString energyWeightLegend[nWeightExponents] = {"n=1", "n=2"};
  TString saveComment =  "_noSystemticsForData";
  TString figureFormat = "pdf";
  saveComment.Prepend(energyWeightString[weightExponent-1]);

  // Ratio zoom settings
  std::pair<double, double> ratioZoom = std::make_pair(0.4, 1.6);
  std::pair<double, double> analysisDeltaR = std::make_pair(0.008, 0.39); // DeltaR span in which the analysis is done
  
  // Marker colors and styles
  int markerStylePbPb[] = {kFullSquare, kFullCircle, kFullCross, kFullFourTrianglesPlus};
  int markerStylePp = kFullDiamond;
  int markerColorPbPb[] = {kRed, kBlue, kMagenta, kGreen+3};
  int markerColorPp = kBlack;
  int bandColorUpPbPb[] = {kOrange+7, kViolet-3, kPink-3, kOrange-3};
  int bandColorDownPbPb[] = {kPink+9, kAzure+8, kViolet+6, kSpring};
  int markerStylePpEnergyWeight[] = {kFullDiamond, kFullDoubleDiamond};
  int markerColorPpEnergyWeight[] = {kViolet-2, kCyan+1};

  TLine* lineDrawer = new TLine();
  lineDrawer->SetLineStyle(2);
  lineDrawer->SetLineColor(kBlack);

  // =============================================== //
  // Read the histograms from the histogram managers //
  // =============================================== //

  // Create histogram managers for result files and systematic uncertainty organizers for systematic uncertainty files
  EECHistogramManager* histograms[kNDataTypes][nWeightExponents];
  SystematicUncertaintyOrganizer* uncertainties[kNDataTypes][nWeightExponents];
  
  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    for(int iDataType = 0; iDataType < kNDataTypes; iDataType++){

      // Create a new histogram manager
      histograms[iDataType][iWeightExponent] = new EECHistogramManager(inputFile[iDataType][iWeightExponent], card[iDataType][iWeightExponent]);
  
      // Load all unfolded energy-energy correlators
      histograms[iDataType][iWeightExponent]->SetLoadEnergyEnergyCorrelators(true);
      histograms[iDataType][iWeightExponent]->SetCentralityBinRange(0, card[iDataType][iWeightExponent]->GetNCentralityBins());
      histograms[iDataType][iWeightExponent]->SetTrackPtBinRangeEEC(0, card[iDataType][iWeightExponent]->GetNTrackPtBinsEEC());
      histograms[iDataType][iWeightExponent]->SetJetPtBinRangeEEC(0, card[iDataType][iWeightExponent]->GetNJetPtBinsEEC());

      // Load the histograms from the file
      histograms[iDataType][iWeightExponent]->LoadProcessedHistograms();

      // Create a new systematic uncertainty organizer
      uncertainties[iDataType][iWeightExponent] = new SystematicUncertaintyOrganizer(uncertaintyFile[iDataType][iWeightExponent]);

    } // Data type loop

  } // Weight exponent loop
 
  // Energy-energy correlators and PbPb to pp ratios
  TH1D* energyEnergyCorrelatorRawPbPb[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorRawPp[nWeightExponents][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorSignalPbPb[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorSignalPp[nWeightExponents][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorPbPbToPpRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyForPbPb[nWeightExponents][knSystematicUncertaintyTypes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyForPp[nWeightExponents][knSystematicUncertaintyTypes][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyPbPbToPpRatio[nWeightExponents][knSystematicUncertaintyTypes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];

  // Energy-energy correlators for Hybrid model
  TGraphErrors* energyEnergyCorrelatorHybridModelPp[nWeightExponents][nJetPtBinsEEC][nTrackPtBinsEEC][2];
  TGraphErrors* energyEnergyCorrelatorHybridModelPbPb[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][2];
  TGraphErrors* energyEnergyCorrelatorHybridModelPbPbToPpRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][2];
  TH1D* histogrammifiedHybridModelPp[nWeightExponents][nJetPtBinsEEC][nTrackPtBinsEEC][2];
  TH1D* histogrammifiedHybridModelPbPb[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][2];
  TH1D* histogrammifiedHybridModelPbPbToPpRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][2];
  TH1D* hybridModelToDataRatioPp[nWeightExponents][nJetPtBinsEEC][nTrackPtBinsEEC][2];
  TH1D* hybridModelToDataRatioPbPb[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][2];
  TH1D* hybridModelToDataRatioPbPbToPpRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][2];

  // Initialize histograms to NULL
  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        energyEnergyCorrelatorRawPp[iWeightExponent][iJetPt][iTrackPt] = NULL;
        energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt] = NULL;

        for(int iUncertainty = 0; iUncertainty < knSystematicUncertaintyTypes; iUncertainty++){
          systematicUncertaintyForPp[iWeightExponent][iUncertainty][iJetPt][iTrackPt] = NULL;
        } // Uncertainty type loop

        for(int iWake = 0; iWake < 2; iWake++){
          energyEnergyCorrelatorHybridModelPp[iWeightExponent][iJetPt][iTrackPt][iWake] = NULL;
          histogrammifiedHybridModelPp[iWeightExponent][iJetPt][iTrackPt][iWake] = NULL;
          hybridModelToDataRatioPp[iWeightExponent][iJetPt][iTrackPt][iWake] = NULL;
        } // Wake type loop

        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          energyEnergyCorrelatorRawPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt] = NULL;
          energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt] = NULL;
          energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt] = NULL;

          for(int iUncertainty = 0; iUncertainty < knSystematicUncertaintyTypes; iUncertainty++){
            systematicUncertaintyForPbPb[iWeightExponent][iUncertainty][iCentrality][iJetPt][iTrackPt] = NULL;
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][iUncertainty][iCentrality][iJetPt][iTrackPt] = NULL;
          }

          for(int iWake = 0; iWake < 2; iWake++){
            energyEnergyCorrelatorHybridModelPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = NULL;
            energyEnergyCorrelatorHybridModelPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = NULL;
            histogrammifiedHybridModelPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = NULL;
            histogrammifiedHybridModelPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = NULL;
            hybridModelToDataRatioPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = NULL;
            hybridModelToDataRatioPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = NULL;
          } // Wake type loop

        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // Weight exponent loop

  // Read the histograms from managers
  int iTrackPt, iTrackPtMatchedPp, iTrackPtMatchedPbPbUncertainty, iTrackPtMatchedPpUncertainty;
  int iJetPt, iJetPtMatchedPp, iJetPtMatchedPbPbUncertainty, iJetPtMatchedPpUncertainty;
  int iCentrality, iCentralityMatched;
  int referenceTrackPtBin, lowAnalysisBin, highAnalysisBin;
  double lowPtIntegral, highPtIntegral;
  double epsilon = 0.0001;

  // Define helper histograms to determine the uncertainties relevant for the double ratio
  TH1D* uncertaintyTrackSelection;
  TH1D* uncertaintyBackgroundSubtraction;
  TH1D* uncertaintyTrackPairEfficiency;
  TH1D* uncertaintyMCnonclosure;
  TH1D* uncertaintySignalToBackgroundRatio;

  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    for(auto jetPtBin : drawnJetPtBin){

      // Match jet pT bin indices
      iJetPt = card[kPbPb][iWeightExponent]->FindBinIndexJetPtEEC(jetPtBin);
      iJetPtMatchedPp = card[kPp][iWeightExponent]->FindBinIndexJetPtEEC(jetPtBin);
      iJetPtMatchedPbPbUncertainty = uncertaintyCard[kPbPb][iWeightExponent]->FindBinIndexJetPtEEC(jetPtBin);
      iJetPtMatchedPpUncertainty = uncertaintyCard[kPp][iWeightExponent]->FindBinIndexJetPtEEC(jetPtBin);

      // Before the main track pT loop, loop over all raw energy-energy correlators to have the information available to scale the tracking related uncertainties in double ratios.
      for(auto trackPtBin : drawnTrackPtBin){
        iTrackPt = card[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);
        iTrackPtMatchedPp = card[kPp][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);

        // Read the raw pp energy-energy correlator distributions
        energyEnergyCorrelatorRawPp[iWeightExponent][iJetPt][iTrackPt] = histograms[kPp][iWeightExponent]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatchedPp, iTrackPtMatchedPp);

        for(auto centralityBin : drawnCentralityBin){

          iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);

          // Read the raw PbPb energy-energy correlator distributions
          energyEnergyCorrelatorRawPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt] = histograms[kPbPb][iWeightExponent]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt);
        }
      }

      // Then go to the main track pT loop
      for(auto trackPtBin : drawnTrackPtBin){
        iTrackPt = card[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);
        iTrackPtMatchedPp = card[kPp][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);
        iTrackPtMatchedPbPbUncertainty = uncertaintyCard[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);
        iTrackPtMatchedPpUncertainty = uncertaintyCard[kPp][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);

        // Read the pp histograms that do not have centrality binning
        energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt] = histograms[kPp][iWeightExponent]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatchedPp, iTrackPtMatchedPp, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);
        systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt] = uncertainties[kPp][iWeightExponent]->GetUncorrelatedSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty);
        systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt] = uncertainties[kPp][iWeightExponent]->GetCorrelatedSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty);

        // Read the relevant systematic uncertainties for double ratio from pp
        uncertaintyBackgroundSubtraction = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kBackgroundSubtraction);
        uncertaintyTrackPairEfficiency = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kTrackPairEfficiency);
        uncertaintyMCnonclosure = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kMonteCarloNonClosure);
        uncertaintyTrackSelection = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kTrackSelection);

        for(auto centralityBin : drawnCentralityBin){

          iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);
          iCentralityMatched = uncertaintyCard[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);

          // Read the PbPb histograms
          energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt] = histograms[kPbPb][iWeightExponent]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);
          systematicUncertaintyForPbPb[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt] = uncertainties[kPbPb][iWeightExponent]->GetUncorrelatedSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty);
          systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt] = uncertainties[kPbPb][iWeightExponent]->GetCorrelatedSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty);

          // Read the relevant systematic uncertainties for double ratio from PbPb
          uncertaintyBackgroundSubtraction = uncertainties[kPbPb][iWeightExponent]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kBackgroundSubtraction);
          uncertaintyTrackPairEfficiency = uncertainties[kPbPb][iWeightExponent]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kTrackPairEfficiency);
          uncertaintyTrackSelection = uncertainties[kPbPb][iWeightExponent]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kTrackSelection);
          uncertaintyMCnonclosure = uncertainties[kPbPb][iWeightExponent]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kMonteCarloNonClosure);
          uncertaintySignalToBackgroundRatio = uncertainties[kPbPb][iWeightExponent]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kSignalToBackgroundRatio);
        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // Weight exponent loop

  // ================================================ //
  //   Normalize the histograms and calculate ratios  //
  // ================================================ //

  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(drawnCentralityBin.at(0));
    iJetPt = card[kPbPb][iWeightExponent]->FindBinIndexJetPtEEC(drawnJetPtBin.at(0));
    iTrackPt = card[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(drawnTrackPtBin.at(0));
    lowAnalysisBin = energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(analysisDeltaR.first + epsilon);
    highAnalysisBin = energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(analysisDeltaR.second - epsilon);
    for(auto jetPtBin : drawnJetPtBin){
      iJetPt = card[kPbPb][iWeightExponent]->FindBinIndexJetPtEEC(jetPtBin);
      for(auto trackPtBin : drawnTrackPtBin){
        iTrackPt = card[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);
        energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->Scale(1.0 / energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
        systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt]->Scale(energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->GetBinContent(10) / systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetBinContent(10));
        systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt]->Scale(energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->GetBinContent(10) / systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt]->GetBinContent(10));

        for(auto centralityBin : drawnCentralityBin){

          iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);

          energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Scale(1.0 / energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
          systematicUncertaintyForPbPb[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Scale(energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->GetBinContent(10) / systematicUncertaintyForPbPb[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetBinContent(10));
          systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Scale(energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->GetBinContent(10) / systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetBinContent(10));

          energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt] = (TH1D*) energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Clone(Form("energyEnergyCorrelatorRatio%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Divide(energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]);

          systematicUncertaintyPbPbToPpRatio[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPbPb[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("systematicUncertaintyRatioUncorrelated%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          systematicUncertaintyPbPbToPpRatio[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Divide(systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt]);

          systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("systematicUncertaintyRatioCorrelated%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Divide(systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt]);
        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // Weight exponent loop

  // For illustration purposes, create up and down shifted uncertainty bands
  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    for(auto jetPtBin : drawnJetPtBin){
      iJetPt = card[kPbPb][iWeightExponent]->FindBinIndexJetPtEEC(jetPtBin);
      for(auto trackPtBin : drawnTrackPtBin){
        iTrackPt = card[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);

        systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertaintyShapeUp][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt]->Clone(Form("correlatedUpBandPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));
        systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertaintyShapeDown][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt]->Clone(Form("correlatedDownBandPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));
        calculateCorrelatedBands(systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt], systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertaintyShapeUp][iJetPt][iTrackPt], systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertaintyShapeDown][iJetPt][iTrackPt]);

        for(auto centralityBin : drawnCentralityBin){

          iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);

          systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("correlatedUpBand%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("correlatedDownBand%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          calculateCorrelatedBands(systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt], systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]);

          systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("correlatedUpBandForRatio%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("correlatedDownBandForRatio%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          calculateCorrelatedBands(systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt], systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]);

        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // Weight exponent loop

  // Create an object for shaking histograms out from graphs
  AlgorithmLibrary* graphShaker = new AlgorithmLibrary();

  // Then read the histograms for the hybrid model
  HybridModelHistogramManager* hybridHistograms = new HybridModelHistogramManager(hybridModelFolder);
  double weightExponentHybrid;
  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    weightExponentHybrid = iWeightExponent+1;
    for(auto jetPtBin : drawnJetPtBin){
      iJetPt = card[kPbPb][iWeightExponent]->FindBinIndexJetPtEEC(jetPtBin);
      for(auto trackPtBin : drawnTrackPtBin){
        iTrackPt = card[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);
        for(int iWake = 0; iWake < 2; iWake++){

          // Load the pp histograms
          energyEnergyCorrelatorHybridModelPp[iWeightExponent][iJetPt][iTrackPt][iWake] = hybridHistograms->GetEnergyEnergyCorrelatorPp(jetPtBin, trackPtBin, weightExponentHybrid, iWake);

          // Create a histogrammified version of the TGraphErrors
          histogrammifiedHybridModelPp[iWeightExponent][iJetPt][iTrackPt][iWake] = graphShaker->HistogrammifyWithErrors(energyEnergyCorrelatorHybridModelPp[iWeightExponent][iJetPt][iTrackPt][iWake], energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]);

          // Now we can take a ratio between hybrid model prediction and data
          hybridModelToDataRatioPp[iWeightExponent][iJetPt][iTrackPt][iWake] = (TH1D*) histogrammifiedHybridModelPp[iWeightExponent][iJetPt][iTrackPt][iWake]->Clone(Form("hybridModelRatioToDataPp%d%d%d%d", iWeightExponent, iJetPt, iTrackPt, iWake));
          hybridModelToDataRatioPp[iWeightExponent][iJetPt][iTrackPt][iWake]->Divide(energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]);

          for(auto centralityBin : drawnCentralityBin){
            iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);

            // Load the PbPb and ratio histograms
            energyEnergyCorrelatorHybridModelPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = hybridHistograms->GetEnergyEnergyCorrelatorPbPb(centralityBin, jetPtBin, trackPtBin, weightExponentHybrid, iWake);
            energyEnergyCorrelatorHybridModelPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = hybridHistograms->GetEnergyEnergyCorrelatorPbPbToPpRatio(centralityBin, jetPtBin, trackPtBin, weightExponentHybrid, iWake);

            // Create a histogrammified versions of the TGraphErrors
            histogrammifiedHybridModelPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = graphShaker->HistogrammifyWithErrors(energyEnergyCorrelatorHybridModelPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake], energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]);
            histogrammifiedHybridModelPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = graphShaker->HistogrammifyWithErrors(energyEnergyCorrelatorHybridModelPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake], energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt]);

            // Now we can take a ratio between hybrid model prediction and data
            hybridModelToDataRatioPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = (TH1D*) histogrammifiedHybridModelPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake]->Clone(Form("hybridModelRatioToDataPbPb%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iWake));
            hybridModelToDataRatioPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake]->Divide(energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]);

            hybridModelToDataRatioPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = (TH1D*) histogrammifiedHybridModelPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake]->Clone(Form("hybridModelRatioToDataPbPbToPpRatio%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iWake));
            hybridModelToDataRatioPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake]->Divide(energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt]);
        
          } // Centrality loop
        } // Wake loop
      } // Track pT loop
    } // Jet pT loop
  } // Weight exponent loop

  // ========================================= //
  //   Draw the histograms in separate plots   //
  // ========================================= //

  JDrawer* drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  drawer->SetLogX(true); // Logarithmic deltaR axis

  TString centralityString;
  TString compactCentralityString;
  TString jetPtString;
  TString compactJetPtString;
  TString trackPtString;
  TString compactTrackPtString;

  // Common variables for different plots
  TLegend* legend;
  TLegend* anotherLegend;
  TLatex* mainTitle;
  TLine* oneLine = new TLine(analysisDeltaR.first, 1, analysisDeltaR.second, 1);
  TBox* box;
  oneLine->SetLineColor(kBlack);
  oneLine->SetLineStyle(2);
  int canvasIndex;
  double bottomRowScale, bottomPadMargin, leftPadMargin;
  double leftMarginAdder, bottomMarginAdder;
  double thisPadScale;

  // Normalization and style for theory predictions
  int color[9] = {kRed, kBlue, kMagenta, kCyan, kGreen+3, kOrange+7, kViolet-3, kPink-3, kOrange-3};
  int iPrediction;

  TString wakeString[2] = {"Hybrid, no wake", "Hybrid, with wake"};

  // Draw individual plots with all centralities mixed together in a single figure
  if(drawDistributionDataToTheoryComparison){

    for(auto jetPtBin : drawnJetPtBin){
      iJetPt = card[kPbPb][weightExponent-1]->FindBinIndexJetPtEEC(jetPtBin);

      jetPtString = Form("%.0f < jet p_{T} < %.0f", jetPtBin.first, jetPtBin.second);
      compactJetPtString = Form("_J=%.0f-%.0f", jetPtBin.first, jetPtBin.second);

      for(auto trackPtBin : drawnTrackPtBin){
        iTrackPt = card[kPbPb][weightExponent-1]->GetBinIndexTrackPtEEC(trackPtBin);

        trackPtString = Form("%.1f < track p_{T}", trackPtBin);
        compactTrackPtString = Form("_T>%.1f", trackPtBin);
        compactTrackPtString.ReplaceAll(".","v");

        // Before going to centrality loop, draw the pp distribution

        // Create a new canvas for the plot
        drawer->CreateSplitCanvas();

        // Use logarithmic axis for EEC
        drawer->SetLogY(true);

        // Setup the legend for plots
        legend = new TLegend(0.23, 0.05, 0.53, 0.6);
        legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);
        legend->AddEntry((TObject*)0, jetPtString.Data(), "");
        legend->AddEntry((TObject*)0, trackPtString.Data(), "");
        legend->AddEntry((TObject*)0, energyWeightLegend[weightExponent-1].Data(), "");

        // Set the drawing style for pp histogram
        energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->SetMarkerStyle(kFullSquare);
        energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->SetMarkerColor(kBlack);

        // Set the x-axis drawing range
        energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);

        // Draw the data correlator to upper canves
        drawer->DrawHistogramToUpperPad(energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt], "#Deltar", "EEC", " ", "p");
        legend->AddEntry(energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt], "pp data", "p");

        // Compare the prediction with and without wake
        for(int iWake = 0; iWake < 2; iWake++){

          // There are some bins for which the prediction does not exist
          if(energyEnergyCorrelatorHybridModelPp[weightExponent-1][iJetPt][iTrackPt][iWake] == NULL) continue;

          // Give some nice styles for the predictions
          energyEnergyCorrelatorHybridModelPp[weightExponent-1][iJetPt][iTrackPt][iWake]->SetLineColor(color[iWake]);
          energyEnergyCorrelatorHybridModelPp[weightExponent-1][iJetPt][iTrackPt][iWake]->SetMarkerColor(color[iWake]);
          energyEnergyCorrelatorHybridModelPp[weightExponent-1][iJetPt][iTrackPt][iWake]->SetMarkerStyle(kFullCircle);
          energyEnergyCorrelatorHybridModelPp[weightExponent-1][iJetPt][iTrackPt][iWake]->SetFillColorAlpha(color[iWake], 0.4);

          // Draw the prediction to the same canvas as the data
          energyEnergyCorrelatorHybridModelPp[weightExponent-1][iJetPt][iTrackPt][iWake]->Draw("3,same");

          // Add a legend for the theory prediction
          legend->AddEntry(energyEnergyCorrelatorHybridModelPp[weightExponent-1][iJetPt][iTrackPt][iWake], wakeString[iWake], "f");
        }

        // Draw the legend to the upper pad
        legend->Draw();

        // Linear scale for the ratio
        drawer->SetLogY(false);

        // Set the axis drawing ranges for ratio
        hybridModelToDataRatioPp[weightExponent-1][iJetPt][iTrackPt][0]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
        hybridModelToDataRatioPp[weightExponent-1][iJetPt][iTrackPt][0]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

        // Set the style for histograms
        for(int iWake = 0; iWake < 2; iWake++){
          hybridModelToDataRatioPp[weightExponent-1][iJetPt][iTrackPt][iWake]->SetMarkerColor(color[iWake]);
          hybridModelToDataRatioPp[weightExponent-1][iJetPt][iTrackPt][iWake]->SetMarkerStyle(kFullCircle);
          hybridModelToDataRatioPp[weightExponent-1][iJetPt][iTrackPt][iWake]->SetMarkerSize(0);
          hybridModelToDataRatioPp[weightExponent-1][iJetPt][iTrackPt][iWake]->SetFillColorAlpha(color[iWake], 0.4);
        }

        drawer->SetGridY(true);
        drawer->DrawHistogramToLowerPad(hybridModelToDataRatioPp[weightExponent-1][iJetPt][iTrackPt][0], "#Deltar", "#frac{Theory}{Data}", " ", "e3");
        for(int iWake = 1; iWake < 2; iWake++){
          hybridModelToDataRatioPp[weightExponent-1][iJetPt][iTrackPt][iWake]->Draw("same,e3");
        }
        drawer->SetGridY(false);

        // If a plot name is given, save the plot in a file
        if(saveFigures) {
          gPad->GetCanvas()->SaveAs(Form("figures/energyEnergyCorrelator_hybridModelDistribution%s_pp%s%s.%s", saveComment.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
        }


        for(auto centralityBin : drawnCentralityBin){
          iCentrality = card[kPbPb][weightExponent-1]->FindBinIndexCentrality(centralityBin);

          centralityString = Form("PbPb %.0f-%.0f%%", centralityBin.first, centralityBin.second);
          compactCentralityString = Form("_C=%.0f-%.0f", centralityBin.first, centralityBin.second);

          // Create a new canvas for the plot
          drawer->CreateSplitCanvas();

          // Use logarithmic axis for EEC
          drawer->SetLogY(true);

          // Setup the legend for plots
          legend = new TLegend(0.23, 0.05, 0.53, 0.6);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);
          legend->AddEntry((TObject*)0, centralityString.Data(), "");
          legend->AddEntry((TObject*)0, jetPtString.Data(), "");
          legend->AddEntry((TObject*)0, trackPtString.Data(), "");
          legend->AddEntry((TObject*)0, energyWeightLegend[weightExponent-1].Data(), "");

          // Set the drawing style for PbPb histogram
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(kFullSquare);
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(kBlack);

          // Set the x-axis drawing range
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);

          // Draw the data correlator to upper canves
          drawer->DrawHistogramToUpperPad(energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt], "#Deltar", "EEC", " ", "p");
          legend->AddEntry(energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt], "PbPb data", "p");

          // Compare the prediction with and without wake
          for(int iWake = 0; iWake < 2; iWake++){

            // There are some bins for which the prediction does not exist
            if(energyEnergyCorrelatorHybridModelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake] == NULL) continue;

            // Give some nice styles for the predictions
            energyEnergyCorrelatorHybridModelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetLineColor(color[iWake]);
            energyEnergyCorrelatorHybridModelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetMarkerColor(color[iWake]);
            energyEnergyCorrelatorHybridModelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetMarkerStyle(kFullCircle);
            energyEnergyCorrelatorHybridModelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetFillColorAlpha(color[iWake], 0.4);

            // Draw the prediction to the same canvas as the data
            energyEnergyCorrelatorHybridModelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->Draw("3,same");

            // Add a legend for the theory prediction
            legend->AddEntry(energyEnergyCorrelatorHybridModelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake], wakeString[iWake], "f");
          }

          // Draw the legend to the upper pad
          legend->Draw();

          // Linear scale for the ratio
          drawer->SetLogY(false);

          // Set the axis drawing ranges for ratio
          hybridModelToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][0]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          hybridModelToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][0]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

          // Set the style for histograms
          for(int iWake = 0; iWake < 2; iWake++){
            hybridModelToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetMarkerColor(color[iWake]);
            hybridModelToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetMarkerStyle(kFullCircle);
            hybridModelToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetMarkerSize(0);
            hybridModelToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetFillColorAlpha(color[iWake], 0.4);
          }

          drawer->SetGridY(true);
          drawer->DrawHistogramToLowerPad(hybridModelToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][0], "#Deltar", "#frac{Theory}{Data}", " ", "e3");
          for(int iWake = 1; iWake < 2; iWake++){
            hybridModelToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->Draw("same,e3");
          }
          drawer->SetGridY(false);

          // If a plot name is given, save the plot in a file
          if(saveFigures) {
            gPad->GetCanvas()->SaveAs(Form("figures/energyEnergyCorrelator_hybridModelDistribution%s%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
          }

        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Comparing theory with the PbPb distribution

  // Draw individual plots with all centralities mixed together in a single figure
  if(drawRatioDataToTheoryComparison){

    for(auto centralityBin : drawnCentralityBin){
      iCentrality = card[kPbPb][weightExponent-1]->FindBinIndexCentrality(centralityBin);

      centralityString = Form("PbPb %.0f-%.0f%%", centralityBin.first, centralityBin.second);
      compactCentralityString = Form("_C=%.0f-%.0f", centralityBin.first, centralityBin.second);

      for(auto jetPtBin : drawnJetPtBin){
        iJetPt = card[kPbPb][weightExponent-1]->FindBinIndexJetPtEEC(jetPtBin);

        jetPtString = Form("%.0f < jet p_{T} < %.0f", jetPtBin.first, jetPtBin.second);
        compactJetPtString = Form("_J=%.0f-%.0f", jetPtBin.first, jetPtBin.second);

        for(auto trackPtBin : drawnTrackPtBin){
          iTrackPt = card[kPbPb][weightExponent-1]->GetBinIndexTrackPtEEC(trackPtBin);

          trackPtString = Form("%.1f < track p_{T}", trackPtBin);
          compactTrackPtString = Form("_T>%.1f", trackPtBin);
          compactTrackPtString.ReplaceAll(".","v");

          // Create a new canvas for the plot
          drawer->CreateSplitCanvas();

          // Use logarithmic axis for EEC
          drawer->SetLogY(true);

          // Setup the legend for plots
          legend = new TLegend(0.23, 0.05, 0.53, 0.6);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);
          legend->AddEntry((TObject*)0, centralityString.Data(), "");
          legend->AddEntry((TObject*)0, jetPtString.Data(), "");
          legend->AddEntry((TObject*)0, trackPtString.Data(), "");
          legend->AddEntry((TObject*)0, energyWeightLegend[weightExponent-1].Data(), "");

          // Set the drawing style for pp histogram
          energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(kFullSquare);
          energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(kBlack);

          // Set the x-axis drawing range
          energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);

          // Draw the pp correlator to upper canves
          drawer->DrawHistogramToUpperPad(energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt], "#Deltar", "#frac{PbPb}{pp}", " ", "p");
          legend->AddEntry(energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt], "Data", "p");

          // Draw the hybrid predictions with and without wake to the same plot
          for(int iWake = 0; iWake < 2; iWake++){

            // There are some bins for which the prediction does not exist
            if(energyEnergyCorrelatorHybridModelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake] == NULL) continue;

            // Give some nice styles for the predictions
            energyEnergyCorrelatorHybridModelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetLineColor(color[iWake]);
            energyEnergyCorrelatorHybridModelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetMarkerColor(color[iWake]);
            energyEnergyCorrelatorHybridModelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetMarkerStyle(kFullCircle);
            energyEnergyCorrelatorHybridModelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetFillColorAlpha(color[iWake], 0.4);

            // Draw the prediction to the same canvas as the data
            energyEnergyCorrelatorHybridModelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->Draw("3,same");

            // Add a legend for the theory prediction
            legend->AddEntry(energyEnergyCorrelatorHybridModelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake], wakeString[iWake], "f");
          }

          // Draw the legend to the upper pad
          legend->Draw();

          // Draw a line to one
          oneLine->Draw();

          // Linear scale for the ratio
          drawer->SetLogY(false);

          // Set the axis drawing ranges
          hybridModelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][0]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          hybridModelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][0]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

          // Set the style for histograms
          for(int iWake = 0; iWake < 2; iWake++){
            hybridModelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetMarkerColor(color[iPrediction]);
            hybridModelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetMarkerStyle(kFullCircle);
            hybridModelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetMarkerSize(0);
            hybridModelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetFillColorAlpha(color[iWake], 0.4);
          }

          drawer->SetGridY(true);
          drawer->DrawHistogramToLowerPad(hybridModelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][0], "#Deltar", "#frac{Theory}{Data}", " ", "e3");
          for(int iWake = 1; iWake < 2; iWake++){
            hybridModelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->Draw("same,e3");
          }
          drawer->SetGridY(false);

          // If a plot name is given, save the plot in a file
          if(saveFigures) {
            gPad->GetCanvas()->SaveAs(Form("figures/energyEnergyCorrelator_hybridModelRatio%s%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
          }

        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Comparing theory with the PbPb distribution
}
