#include "SystematicUncertaintyOrganizer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECHistogramManager.h"
#include "EECCard.h"
#include "JDrawer.h"
#include "AlgorithmLibrary.h"

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
void theoryFitter(){

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

  // eecAnalysis_akFlowJet_nominalEnergyWeight_optimizedUnfoldingBins_fixedCovarianceMatrix_updatedBackgroundSubtraction_processed_2024-02-23.root
  // eecAnalysis_akFlowJet_energyWeightSquared_optimizedUnfoldingBins_fixedCovarianceMatrix_updatedBackgroundSubtraction_processed_2024-02-23.root
  // ppData_pfJets_wtaAxis_nominalEnergyWeight_optimizedUnfoldingBins_fixedCovarianceMatrix_jet60or80triggers_unfoldingWithCovariance_processed_2024-01-23.root
  // ppData_pfJets_wtaAxis_nominalEnergyWeight_optimizedUnfoldingBins_jet60or80triggers_unfoldingWithNominalSmear_processed_2024-01-17.root
  // ppData_pfJets_wtaAxis_energyWeightSquared_optimizedUnfoldingBins_fixedCovarianceMatrix_jet60or80triggers_unfoldingWithCovariance_processed_2024-01-23.root
  // ppData_pfJets_wtaAxis_energyWeightSquared_optimizedUnfoldingBins_jet60or80triggers_unfoldingWithNominalSmear_processed_2024-01-17.root
  
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
  
  // The final results are available for all the bins that are unfolded
  int firstDrawnCentralityBin[nWeightExponents];
  int lastDrawnCentralityBin[nWeightExponents];
  
  int firstDrawnJetPtBinEEC[nWeightExponents];
  int lastDrawnJetPtBinEEC[nWeightExponents];
  
  int firstDrawnTrackPtBinEEC[nWeightExponents];
  int lastDrawnTrackPtBinEEC[nWeightExponents];

  // Select the bins for which the theory predictions are available
  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    firstDrawnCentralityBin[iWeightExponent] = 0;
    lastDrawnCentralityBin[iWeightExponent] = 0;
    firstDrawnJetPtBinEEC[iWeightExponent] = 6;
    lastDrawnJetPtBinEEC[iWeightExponent] = 6;
    firstDrawnTrackPtBinEEC[iWeightExponent] = 1;
    lastDrawnTrackPtBinEEC[iWeightExponent] = 1;
  }

  // Choose which plots to draw
  bool drawDistributionDataToTheoryComparison = true;
  bool drawRatioDataToTheoryComparison = true;

  bool drawVerticalLines = false; // Draw illustrative vertical lines

  // Normalize all distributions to 2 GeV integral
  bool normalizeTo2GeV = false;
  int trackPtBinFor2GeV[nWeightExponents];

  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    trackPtBinFor2GeV[iWeightExponent] = card[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(2.0);
  }
  
  // Save the final plots
  const bool saveFigures = true;
  TString energyWeightString[nWeightExponents] = {"_nominalEnergyWeight", "_energyWeightSquared"};
  TString saveComment =  "_newColorScheme";
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

  // Initialize histograms to NULL
  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        energyEnergyCorrelatorRawPp[iWeightExponent][iJetPt][iTrackPt] = NULL;
        energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt] = NULL;

        for(int iUncertainty = 0; iUncertainty < knSystematicUncertaintyTypes; iUncertainty++){
          systematicUncertaintyForPp[iWeightExponent][iUncertainty][iJetPt][iTrackPt] = NULL;
        } // Uncertainty type loop

        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          energyEnergyCorrelatorRawPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt] = NULL;
          energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt] = NULL;
          energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt] = NULL;

          for(int iUncertainty = 0; iUncertainty < knSystematicUncertaintyTypes; iUncertainty++){
            systematicUncertaintyForPbPb[iWeightExponent][iUncertainty][iCentrality][iJetPt][iTrackPt] = NULL;
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][iUncertainty][iCentrality][iJetPt][iTrackPt] = NULL;
          }
        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // Weight exponent loop

  // Read the histograms from managers
  int iTrackPtMatchedPp, iTrackPtMatchedPbPbUncertainty, iTrackPtMatchedPpUncertainty;
  int iJetPtMatchedPp, iJetPtMatchedPbPbUncertainty, iJetPtMatchedPpUncertainty;
  int iCentralityMatched;
  int referenceTrackPtBin, lowAnalysisBin, highAnalysisBin;
  double lowPtIntegral, highPtIntegral;
  double epsilon = 0.0001;

  // Define helper histograms to determine the uncertainties relevant for the double ratio
  TH1D* uncertaintyTrackSelection;
  TH1D* uncertaintyBackgroundSubtraction;
  TH1D* uncertaintyTrackPairEfficiency;
  TH1D* uncertaintyMCnonclosure;
  TH1D* uncertaintySignalToBackgroundRatio;
  double trackCorrelation;
  std::pair<double,double> currentJetPtBin;

  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    for(int iJetPt = firstDrawnJetPtBinEEC[iWeightExponent]; iJetPt <= lastDrawnJetPtBinEEC[iWeightExponent]; iJetPt++){

      // Find the current jet pT bin
      currentJetPtBin = card[kPbPb][iWeightExponent]->GetBinBordersJetPtEEC(iJetPt);

      iJetPtMatchedPp = card[kPp][iWeightExponent]->FindBinIndexJetPtEEC(currentJetPtBin);
      iJetPtMatchedPbPbUncertainty = uncertaintyCard[kPbPb][iWeightExponent]->FindBinIndexJetPtEEC(currentJetPtBin);
      iJetPtMatchedPpUncertainty = uncertaintyCard[kPp][iWeightExponent]->FindBinIndexJetPtEEC(currentJetPtBin);

      // Before the main track pT loop, loop over all raw energy-energy correlators to have the information available to scale the tracking related uncertainties in double ratios.
      for(int iTrackPt = firstDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt <= lastDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt++){
        iTrackPtMatchedPp = card[kPp][iWeightExponent]->FindBinIndexTrackPtEEC(card[kPbPb][iWeightExponent]->GetBinBordersTrackPtEEC(iTrackPt));

        // Read the raw pp energy-energy correlator distributions
        energyEnergyCorrelatorRawPp[iWeightExponent][iJetPt][iTrackPt] = histograms[kPp][iWeightExponent]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatchedPp, iTrackPtMatchedPp);

        for(int iCentrality = firstDrawnCentralityBin[iWeightExponent]; iCentrality <= lastDrawnCentralityBin[iWeightExponent]; iCentrality++){

          // Read the raw PbPb energy-energy correlator distributions
          energyEnergyCorrelatorRawPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt] = histograms[kPbPb][iWeightExponent]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt);
        }
      }

      // Then go to the main track pT loop
      for(int iTrackPt = firstDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt <= lastDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt++){
        iTrackPtMatchedPp = card[kPp][iWeightExponent]->FindBinIndexTrackPtEEC(card[kPbPb][iWeightExponent]->GetBinBordersTrackPtEEC(iTrackPt));
        iTrackPtMatchedPbPbUncertainty = uncertaintyCard[kPbPb][iWeightExponent]->FindBinIndexTrackPtEEC(card[kPbPb][iWeightExponent]->GetBinBordersTrackPtEEC(iTrackPt));
        iTrackPtMatchedPpUncertainty = uncertaintyCard[kPp][iWeightExponent]->FindBinIndexTrackPtEEC(card[kPbPb][iWeightExponent]->GetBinBordersTrackPtEEC(iTrackPt));

        // Read the pp histograms that do not have centrality binning
        energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt] = histograms[kPp][iWeightExponent]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatchedPp, iTrackPtMatchedPp, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);
        systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt] = uncertainties[kPp][iWeightExponent]->GetUncorrelatedSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty);
        systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt] = uncertainties[kPp][iWeightExponent]->GetCorrelatedSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty);

        // Read the relevant systematic uncertainties for double ratio from pp
        uncertaintyBackgroundSubtraction = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kBackgroundSubtraction);
        uncertaintyTrackPairEfficiency = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kTrackPairEfficiency);
        uncertaintyMCnonclosure = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kMonteCarloNonClosure);
        uncertaintyTrackSelection = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kTrackSelection);

        for(int iCentrality = firstDrawnCentralityBin[iWeightExponent]; iCentrality <= lastDrawnCentralityBin[iWeightExponent]; iCentrality++){
          iCentralityMatched = uncertaintyCard[kPbPb][iWeightExponent]->FindBinIndexCentrality(card[kPbPb][iWeightExponent]->GetBinBordersCentrality(iCentrality));

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
    referenceTrackPtBin = trackPtBinFor2GeV[iWeightExponent];
    lowAnalysisBin = energyEnergyCorrelatorSignalPbPb[iWeightExponent][firstDrawnCentralityBin[iWeightExponent]][firstDrawnJetPtBinEEC[iWeightExponent]][firstDrawnTrackPtBinEEC[iWeightExponent]]->GetXaxis()->FindBin(analysisDeltaR.first + epsilon);
    highAnalysisBin = energyEnergyCorrelatorSignalPbPb[iWeightExponent][firstDrawnCentralityBin[iWeightExponent]][firstDrawnJetPtBinEEC[iWeightExponent]][firstDrawnTrackPtBinEEC[iWeightExponent]]->GetXaxis()->FindBin(analysisDeltaR.second - epsilon);
    for(int iJetPt = firstDrawnJetPtBinEEC[iWeightExponent]; iJetPt <= lastDrawnJetPtBinEEC[iWeightExponent]; iJetPt++){
      for(int iTrackPt = lastDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt >= firstDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt--){
        if(!normalizeTo2GeV) referenceTrackPtBin = iTrackPt;
        energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->Scale(1.0 / energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][referenceTrackPtBin]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
        systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt]->Scale(energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->GetBinContent(10) / systematicUncertaintyForPp[iWeightExponent][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetBinContent(10));
        systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt]->Scale(energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->GetBinContent(10) / systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt]->GetBinContent(10));

        for(int iCentrality = firstDrawnCentralityBin[iWeightExponent]; iCentrality <= lastDrawnCentralityBin[iWeightExponent]; iCentrality++){
          energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Scale(1.0 / energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][referenceTrackPtBin]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
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
    for(int iJetPt = firstDrawnJetPtBinEEC[iWeightExponent]; iJetPt <= lastDrawnJetPtBinEEC[iWeightExponent]; iJetPt++){
      for(int iTrackPt = firstDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt <= lastDrawnTrackPtBinEEC[iWeightExponent]; iTrackPt++){

        systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertaintyShapeUp][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt]->Clone(Form("correlatedUpBandPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));
        systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertaintyShapeDown][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt]->Clone(Form("correlatedDownBandPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));
        calculateCorrelatedBands(systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertainty][iJetPt][iTrackPt], systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertaintyShapeUp][iJetPt][iTrackPt], systematicUncertaintyForPp[iWeightExponent][kCorrelatedUncertaintyShapeDown][iJetPt][iTrackPt]);

        for(int iCentrality = firstDrawnCentralityBin[iWeightExponent]; iCentrality <= lastDrawnCentralityBin[iWeightExponent]; iCentrality++){

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

  // After the data distributions are ready, get the theory predictions
  AlgorithmLibrary *theoryProvider = new AlgorithmLibrary();
  auto theoryCurvesForDistribution = theoryProvider->GetGraphsFromDatFile("theoryComparison/ForFittingpt140-160b3p2GLV.dat");
  auto theoryCurvesForRatio = theoryProvider->GetGraphsFromDatFile("theoryComparison/ForFittingpt140-160b3p2GLVratioAA_pp.dat");

  // Data to theory ratio histograms
  const int nPredictions = theoryCurvesForDistribution.size();
  TH1D* histogrammifiedTheoryDistribution[nPredictions];
  TH1D* histogrammifiedTheoryRatio[nPredictions];
  TH1D* distributionToTheoryRatio[nPredictions];
  TH1D* ratioToTheoryRatio[nPredictions];
  for(int iPrediction = 0; iPrediction < nPredictions; iPrediction++){
    histogrammifiedTheoryDistribution[iPrediction] = NULL;
    histogrammifiedTheoryRatio[iPrediction] = NULL;
    distributionToTheoryRatio[iPrediction] = NULL;
    ratioToTheoryRatio[iPrediction] = NULL;
  }

  // Find normalization by matching the integrals of data and histogrammified theory in the free quark/gluon region
  double dataIntegralDistribution, theoryIntegralDistribution;
  double dataIntegralRatio, theoryIntegralRatio;
  double scaleFactorDistribution[nPredictions];
  double scaleFactorRatio[nPredictions];
  double normalizationRegionDeltaR;
  for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
    for(int iJetPt = firstDrawnJetPtBinEEC[weightExponent-1]; iJetPt <= lastDrawnJetPtBinEEC[weightExponent-1]; iJetPt++){
      for(int iTrackPt = firstDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt <= lastDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt++){

        // Calculate the integrals in the data
        dataIntegralDistribution = energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->Integral(13, 27, "width");
        dataIntegralRatio = energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]->Integral(13, 27, "width");
        normalizationRegionDeltaR = energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->GetXaxis()->GetBinLowEdge(13);

        // Loop over all predictions and find scaling factors
        for(int iPrediction = 0; iPrediction < nPredictions; iPrediction++){

          // Find scaling factor for distributions
          histogrammifiedTheoryDistribution[iPrediction] = theoryProvider->Histogrammify(theoryCurvesForDistribution.at(iPrediction).second, energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]);
          theoryIntegralDistribution = histogrammifiedTheoryDistribution[iPrediction]->Integral(13, 27, "width");
          scaleFactorDistribution[iPrediction] = dataIntegralDistribution/theoryIntegralDistribution;

          // Find scaling factors for ratios
          histogrammifiedTheoryRatio[iPrediction] = theoryProvider->Histogrammify(theoryCurvesForRatio.at(iPrediction).second, energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]);
          theoryIntegralRatio = histogrammifiedTheoryRatio[iPrediction]->Integral(13, 27, "width");
          scaleFactorRatio[iPrediction] = dataIntegralRatio/theoryIntegralRatio;

        }

      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop


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

  // Draw individual plots with all centralities mixed together in a single figure
  if(drawDistributionDataToTheoryComparison){

    for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
      centralityString = Form("PbPb %.0f-%.0f%%", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality));
      compactCentralityString = Form("_C=%.0f-%.0f", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality));

      for(int iJetPt = firstDrawnJetPtBinEEC[weightExponent-1]; iJetPt <= lastDrawnJetPtBinEEC[weightExponent-1]; iJetPt++){
        jetPtString = Form("%.0f < jet p_{T} < %.0f", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt));
        compactJetPtString = Form("_J=%.0f-%.0f", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt));

        for(int iTrackPt = firstDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt <= lastDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt++){
          trackPtString = Form("%.1f < track p_{T}", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt));
          compactTrackPtString = Form("_T>%.1f", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt));
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

          // Set the drawing style for pp histogram
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(kFullSquare);
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(kBlack);

          // Set the x-axis drawing range
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);

          // Draw the pp correlator to upper canves
          drawer->DrawHistogramToUpperPad(energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt], "#Deltar", "EEC", " ", "p");
          legend->AddEntry(energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt], "PbPb data", "p");

          // Draw the different centrality bins to the same plot
          iPrediction = 0;
          for(auto thisPrediction : theoryCurvesForDistribution){
            thisPrediction.second->Scale(scaleFactorDistribution[iPrediction]);

            // Give some nice styles for the predictions
            thisPrediction.second->SetLineColor(color[iPrediction]);
            thisPrediction.second->SetMarkerColor(color[iPrediction]);
            thisPrediction.second->SetMarkerStyle(kFullCircle);

            // Draw the prediction to the same canvas as the data
            thisPrediction.second->Draw("pl,same");

            // Add a legend for the theory prediction
            legend->AddEntry(thisPrediction.second, thisPrediction.first, "p");

            // Increment the counter for styles
            iPrediction++;
          }

          // Draw the legend to the upper pad
          legend->Draw();

          // Draw a line showing the normalization region for predictions
          lineDrawer->DrawLine(normalizationRegionDeltaR, 1, normalizationRegionDeltaR, 20);

          // Linear scale for the ratio
          drawer->SetLogY(false);

          // Take a ratio between histogrammified graph and data
          iPrediction = 0;
          for(auto thisPrediction : theoryCurvesForDistribution){
            distributionToTheoryRatio[iPrediction] = theoryProvider->Histogrammify(thisPrediction.second, energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]);
            distributionToTheoryRatio[iPrediction]->Divide(energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]);
            iPrediction++;
          }

          // Set the axis drawing ranges
          distributionToTheoryRatio[0]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          distributionToTheoryRatio[0]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

          // Set the style for histograms
          for(int iPrediction = 0; iPrediction < nPredictions; iPrediction++){
            distributionToTheoryRatio[iPrediction]->SetMarkerColor(color[iPrediction]);
            distributionToTheoryRatio[iPrediction]->SetMarkerStyle(kFullCircle);
          }

          drawer->SetGridY(true);
          drawer->DrawHistogramToLowerPad(distributionToTheoryRatio[0], "#Deltar", "#frac{Theory}{Data}", " ", "p");
          for(int iPrediction = 1; iPrediction < nPredictions; iPrediction++){
            distributionToTheoryRatio[iPrediction]->Draw("same,p");
          }
          drawer->SetGridY(false);

          // Draw a line showing the normalization region for predictions
          lineDrawer->DrawLine(normalizationRegionDeltaR, 0.5, normalizationRegionDeltaR, 1.5);

          // If a plot name is given, save the plot in a file
          if(saveFigures) {
            gPad->GetCanvas()->SaveAs(Form("figures/energyEnergyCorrelator_distributionTheoryFit%s%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
          }

        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Comparing theory with the PbPb distribution

  // Draw individual plots with all centralities mixed together in a single figure
  if(drawRatioDataToTheoryComparison){

    for(int iCentrality = firstDrawnCentralityBin[weightExponent-1]; iCentrality <= lastDrawnCentralityBin[weightExponent-1]; iCentrality++){
      centralityString = Form("PbPb %.0f-%.0f%%", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality));
      compactCentralityString = Form("_C=%.0f-%.0f", card[kPbPb][weightExponent-1]->GetLowBinBorderCentrality(iCentrality), card[kPbPb][weightExponent-1]->GetHighBinBorderCentrality(iCentrality));

      for(int iJetPt = firstDrawnJetPtBinEEC[weightExponent-1]; iJetPt <= lastDrawnJetPtBinEEC[weightExponent-1]; iJetPt++){
        jetPtString = Form("%.0f < jet p_{T} < %.0f", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt));
        compactJetPtString = Form("_J=%.0f-%.0f", card[kPbPb][weightExponent-1]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb][weightExponent-1]->GetHighBinBorderJetPtEEC(iJetPt));

        for(int iTrackPt = firstDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt <= lastDrawnTrackPtBinEEC[weightExponent-1]; iTrackPt++){
          trackPtString = Form("%.1f < track p_{T}", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt));
          compactTrackPtString = Form("_T>%.1f", card[kPbPb][weightExponent-1]->GetLowBinBorderTrackPtEEC(iTrackPt));
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

          // Set the drawing style for pp histogram
          energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(kFullSquare);
          energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(kBlack);

          // Set the x-axis drawing range
          energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);

          // Draw the pp correlator to upper canves
          drawer->DrawHistogramToUpperPad(energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt], "#Deltar", "#frac{PbPb}{pp}", " ", "p");
          legend->AddEntry(energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt], "Data", "p");

          // Draw the different centrality bins to the same plot
          iPrediction = 0;
          for(auto thisPrediction : theoryCurvesForRatio){
            thisPrediction.second->Scale(scaleFactorRatio[iPrediction]);

            // Give some nice styles for the predictions
            thisPrediction.second->SetLineColor(color[iPrediction]);
            thisPrediction.second->SetMarkerColor(color[iPrediction]);
            thisPrediction.second->SetMarkerStyle(kFullCircle);

            // Draw the prediction to the same canvas as the data
            thisPrediction.second->Draw("pl,same");

            // Add a legend for the theory prediction
            legend->AddEntry(thisPrediction.second, thisPrediction.first, "p");

            // Increment the counter for styles
            iPrediction++;
          }

          // Draw the legend to the upper pad
          legend->Draw();

          // Draw a line showing the normalization region for predictions
          lineDrawer->DrawLine(normalizationRegionDeltaR, 0.5, normalizationRegionDeltaR, 2);

          // Linear scale for the ratio
          drawer->SetLogY(false);

          // Take a ratio between histogrammified graph and data
          iPrediction = 0;
          for(auto thisPrediction : theoryCurvesForRatio){
            ratioToTheoryRatio[iPrediction] = theoryProvider->Histogrammify(thisPrediction.second, energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]);
            ratioToTheoryRatio[iPrediction]->Divide(energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]);
            iPrediction++;
          }

          // Set the axis drawing ranges
          ratioToTheoryRatio[0]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          ratioToTheoryRatio[0]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

          // Set the style for histograms
          for(int iPrediction = 0; iPrediction < nPredictions; iPrediction++){
            ratioToTheoryRatio[iPrediction]->SetMarkerColor(color[iPrediction]);
            ratioToTheoryRatio[iPrediction]->SetMarkerStyle(kFullCircle);
          }

          drawer->SetGridY(true);
          drawer->DrawHistogramToLowerPad(ratioToTheoryRatio[0], "#Deltar", "#frac{Theory}{Data}", " ", "p");
          for(int iPrediction = 1; iPrediction < nPredictions; iPrediction++){
            ratioToTheoryRatio[iPrediction]->Draw("same,p");
          }
          drawer->SetGridY(false);

          // Draw a line showing the normalization region for predictions
          lineDrawer->DrawLine(normalizationRegionDeltaR, 0.5, normalizationRegionDeltaR, 1.5);

          // If a plot name is given, save the plot in a file
          if(saveFigures) {
            gPad->GetCanvas()->SaveAs(Form("figures/energyEnergyCorrelator_ratioTheoryFit%s%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
          }

        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Comparing theory with the PbPb distribution
}
