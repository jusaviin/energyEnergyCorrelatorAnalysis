#include "SystematicUncertaintyOrganizer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECHistogramManager.h"
#include "EECCard.h"
#include "JDrawer.h"
#include "AlgorithmLibrary.h"
#include "HybridModelHistogramManager.h"
#include "HolguinHistogramManager.h"
#include "CoLBTHistogramManager.h"
#include "JewelHistogramManager.h"

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
 *  Arguments:
 *   int weightExponent = Drawn energy weight exponent
 *   int theoryComparisonIndex  = Index for the selected theory predictions
 *       0 = Hybrid model with different wake configuration
 *       1 = Holguin perturbative calculations with different k-values
 *       2 = CoLBT predictions with different q-values
 *       3 = Selection of Hybrid, Holguin and CoLBT predictions
 *       4 = Best k from Holguin and CoLBT
 *       5 = Only JEWEL
 *       6 = Best k from Holguin and JEWEL
 *   int drawnPlots = Index to select which plots are drawn to the canvas
 *       0 = Compare energy-energy correlator distributions to theory predictions
 *       1 = Compare PbPb/pp ratios to theory predictions
 *       2 = Compare double ratios to theory predictions
 *   bool paperPlorMode = Use default configuration for drawing all the plots included in the paper
 *       true = Force default paper plot configuration
 *       false = Allow user to make their own configuration
 */
void modelComparison(int weightExponent = 1, int theoryComparisonIndex = 0, int drawnPlots = 0, bool paperPlotMode = false){

  enum enumDataType{kPbPb, kPp, kNDataTypes};
  enum enumSystematicUncertaintyType{kUncorrelatedUncertainty, kCorrelatedUncertainty, kCorrelatedUncertaintyShapeUp, kCorrelatedUncertaintyShapeDown, knSystematicUncertaintyTypes};
  enum enumrelativeUncertaintyType{kRelativeUncertaintyStatisticalUp, kRelativeUncertaintyStatisticalDown, kRelativeUncertaintySystematic, knRelativeUncertaintyTypes};
  enum enumPpMCTypes{kPythia, kHerwig, knPpMCTypes};
  enum enumPpModelComparison{kPpCompareHybrid, kPpComparePythia, kPpCompareHerwig, knPpModelComparisons};

  // ============= //
  // Configuration //
  // ============= //

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
  uncertaintyFileName[kPbPb][0] = "systematicUncertainties/systematicUncertainties_PbPb_nominalEnergyWeight_combinedMxedConeBackground_includeMCstats_2025-04-18.root";
  uncertaintyFileName[kPbPb][1] = "systematicUncertainties/systematicUncertainties_PbPb_energyWeightSquared_combinedMxedConeBackground_includeMCstats_2025-04-18.root";
  uncertaintyFileName[kPp][0] = "systematicUncertainties/systematicUncertainties_pp_nominalEnergyWeight_includeMCstats_2025-04-21.root";
  uncertaintyFileName[kPp][1] = "systematicUncertainties/systematicUncertainties_pp_energyWeightSquared_includeMCstats_2025-04-21.root";

  // Input files for pp MC results
  TString ppMCFileName[knPpMCTypes][nWeightExponents];
  ppMCFileName[kPythia][0] = "data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_optimizedUnfoldingBins_nominalSmear_truthReference_processed_2024-01-11.root";
  ppMCFileName[kPythia][1] = "data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_optimizedUnfoldingBins_energyWeightSquared_nominalSmear_truthReference_processed_2024-01-10.root";
  ppMCFileName[kHerwig][0] = "data/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_nominalEnergyWeight_officialSimulation_nominalSmear_processed_2024-06-24.root";
  ppMCFileName[kHerwig][1] = "data/ppMC2017_GenGen_Herwig_pfJets_wtaAxis_energyWeightSquared_officialSimulation_nominalSmear_processed_2024-06-24.root";


  // String pointing to the folder where the hybrid model predictions are located
  TString hybridModelFolder = "theoryComparison/hybridModel/data";
  TString holguinDataFolder = "theoryComparison/holguin";
  TString coLBTFolder = "theoryComparison/colbt";
  TString jewelDataFolder = "theoryComparison/jewel";
  
  TFile* inputFile[kNDataTypes][nWeightExponents];
  TFile* uncertaintyFile[kNDataTypes][nWeightExponents];
  TFile* ppMCFile[knPpMCTypes][nWeightExponents];
  EECCard* card[kNDataTypes][nWeightExponents];
  EECCard* uncertaintyCard[kNDataTypes][nWeightExponents];
  EECCard* ppMCCard[knPpMCTypes][nWeightExponents];
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

    for(int iFile = 0; iFile < knPpMCTypes; iFile++){

      // Load the pp MC file
      ppMCFile[iFile][iWeightExponent] = TFile::Open(ppMCFileName[iFile][iWeightExponent]);

      // Check that the input file exists
      if(ppMCFile[iFile][iWeightExponent] == NULL){
        cout << "Error! The file " << ppMCFileName[iFile][iWeightExponent].Data() << " does not exist!" << endl;
        cout << "Maybe you forgot the data/ folder path?" << endl;
        cout << "Will not execute the code" << endl;
        return;
      }
    
      // Load the card from the file and read the collision system
      ppMCCard[iFile][iWeightExponent] = new EECCard(ppMCFile[iFile][iWeightExponent]);

    }

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

  // Always include all bins in the paper plot mode
  if(paperPlotMode){
    drawnCentralityBin.clear();
    drawnCentralityBin.push_back(std::make_pair(0,10));
    drawnCentralityBin.push_back(std::make_pair(10,30));
    drawnCentralityBin.push_back(std::make_pair(30,50));
    drawnCentralityBin.push_back(std::make_pair(50,90));
 
    drawnJetPtBin.clear();
    drawnJetPtBin.push_back(std::make_pair(120,140));
    drawnJetPtBin.push_back(std::make_pair(140,160));
    drawnJetPtBin.push_back(std::make_pair(160,180));
    drawnJetPtBin.push_back(std::make_pair(180,200));

    drawnTrackPtBin.clear();
    drawnTrackPtBin.push_back(1.0);
    drawnTrackPtBin.push_back(2.0);
  }

  // Define which k-values to study for Holguin predictions
  std::vector<double> holguinKValue;
  holguinKValue.push_back(0.1);
  holguinKValue.push_back(0.3);
  holguinKValue.push_back(0.5);

  // Define the q-values to study for CoLBT predictions
  std::vector<double> coLBTqValue;
  coLBTqValue.push_back(0.5);
  coLBTqValue.push_back(1);

  // Index for the selected theory predictions
  // 0 = Hybrid model with different wake configuration
  // 1 = Holguin perturbative calculations with different k-values
  // 2 = CoLBT predictions with different q-values
  // 3 = Selection of Hybrid, Holguin and CoLBT predictions
  // 4 = Best k from Holguin and CoLBT
  // 5 = Only JEWEL
  // 6 = Best k from Holguin and JEWEL
  TString theorySaveName[7] = {"hybridModel", "holguin", "colbt", "threeModels", "holguinAndColbt", "jewel", "holguinAndJewel"};


  // Binning for double ratio plots
  std::pair<double, double> trackPtCutsForDoubleDatio = std::make_pair(1.0, 2.0);
  std::pair<int, int> trackPtBinsForDoubleRatio = std::make_pair(card[kPbPb][0]->GetBinIndexTrackPtEEC(trackPtCutsForDoubleDatio.first), card[kPbPb][0]->GetBinIndexTrackPtEEC(trackPtCutsForDoubleDatio.second));

  // Choose which plots to draw
  bool drawDistributionDataToTheoryComparison = (drawnPlots == 0);
  bool drawRatioDataToTheoryComparison = (drawnPlots == 1);
  bool drawDoubleRatioDataToTheoryComparison = (drawnPlots == 2);

  // Flag for adding preliminary tag to the figures
  bool addPreliminaryTag = false;

  // Legend mode where the differences in plots shown side-by-side in the main text are highlighted
  bool mainTextHighlightMode = false;

  // Save the final plots
  const bool saveFigures = true;
  TString energyWeightString[nWeightExponents] = {"_nominalEnergyWeight", "_energyWeightSquared"};
  TString energyWeightLegend[nWeightExponents] = {"n=1", "n=2"};
  TString saveComment =  "";
  TString figureFormat = "pdf";

  if(paperPlotMode){
    saveComment =  "";
    figureFormat = "pdf";
    addPreliminaryTag = false;
    mainTextHighlightMode = false;
  }

  saveComment.Prepend(energyWeightString[weightExponent-1]);

  // Ratio zoom settings
  std::pair<double, double> distributionZoom[nWeightExponents] = {std::make_pair(0.12, 25), std::make_pair(0.04, 60)};
  std::pair<double, double> ratioZoomPpDistribution = std::make_pair(0.71, 1.29);
  std::pair<double, double> ratioZoom = std::make_pair(0.4, 1.6);
  std::pair<double, double> analysisDeltaR = std::make_pair(0.008, 0.39); // DeltaR span in which the analysis is done
  std::pair<double, double> pbpbToPpRatioZoom = std::make_pair(0.16, 1.84);
  std::pair<double, double> doubleRatioZoom = std::make_pair(0.5, 1.5);
  
  // Marker colors and styles
  int markerStylePbPb[] = {kFullSquare, kFullCircle, kFullCross, kFullFourTrianglesPlus};
  int markerStylePp = kFullDiamond;
  int markerStyleTheory[] = {kFullCircle, kFullCross, kFullCrossX};
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
  EECHistogramManager* ppMCHistograms[knPpMCTypes][nWeightExponents];
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

    for(int iMCType = 0; iMCType < knPpMCTypes; iMCType++){

      // Create a new histogram manager
      ppMCHistograms[iMCType][iWeightExponent] = new EECHistogramManager(ppMCFile[iMCType][iWeightExponent], ppMCCard[iMCType][iWeightExponent]);
  
      // Load all unfolded energy-energy correlators
      ppMCHistograms[iMCType][iWeightExponent]->SetLoadEnergyEnergyCorrelators(true);
      ppMCHistograms[iMCType][iWeightExponent]->SetCentralityBinRange(0, ppMCCard[iMCType][iWeightExponent]->GetNCentralityBins());
      ppMCHistograms[iMCType][iWeightExponent]->SetTrackPtBinRangeEEC(0, ppMCCard[iMCType][iWeightExponent]->GetNTrackPtBinsEEC());
      ppMCHistograms[iMCType][iWeightExponent]->SetJetPtBinRangeEEC(0, ppMCCard[iMCType][iWeightExponent]->GetNJetPtBinsEEC());

      // Load the histograms from the file
      ppMCHistograms[iMCType][iWeightExponent]->LoadProcessedHistograms();

    }

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

  // Histograms for relative uncertainties
  TH1D* hRelativeUncertaintyPp[nWeightExponents][nJetPtBinsEEC][nTrackPtBinsEEC][knRelativeUncertaintyTypes];
  TH1D* hRelativeUncertaintyPbPb[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][knRelativeUncertaintyTypes];
  TH1D* hRelativeUncertaintyPbPbToPpRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][knRelativeUncertaintyTypes];

  // Double ratios from energy-energy correlator histograms
  TH1D* energyEnergyCorrelatorForDoubleRatioFromPbPb[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorForDoubleRatioFromPp[nWeightExponents][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorDoubleRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC];
  TH1D* systematicUncertaintyForDoubleRatioFromPbPb[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyForDoubleRatioFromPp[nWeightExponents][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyDoubleRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC];

  // Relative uncertainties for double ratios
  TH1D* hRelativeUncertaintyDoubleRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC][knRelativeUncertaintyTypes];

  // Energy-energy correlators from pp MC models
  TH1D* hEnergyEnergyCorrelatorPpMC[nWeightExponents][nJetPtBinsEEC][nTrackPtBinsEEC][knPpMCTypes];
  TH1D* hEnergyEnergyCorrelatorPpMCToDataRatio[nWeightExponents][nJetPtBinsEEC][nTrackPtBinsEEC][knPpMCTypes];

  // Energy-energy correlators for Hybrid model
  TGraphErrors* energyEnergyCorrelatorHybridModelPp[nWeightExponents][nJetPtBinsEEC][nTrackPtBinsEEC][HybridModelHistogramManager::kWakeConfigurations];
  TGraphErrors* energyEnergyCorrelatorHybridModelPbPb[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][HybridModelHistogramManager::kWakeConfigurations];
  TGraphErrors* energyEnergyCorrelatorHybridModelPbPbToPpRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][HybridModelHistogramManager::kWakeConfigurations];
  TH1D* histogrammifiedHybridModelPp[nWeightExponents][nJetPtBinsEEC][nTrackPtBinsEEC][HybridModelHistogramManager::kWakeConfigurations];
  TH1D* histogrammifiedHybridModelPbPb[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][HybridModelHistogramManager::kWakeConfigurations];
  TH1D* histogrammifiedHybridModelPbPbToPpRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][HybridModelHistogramManager::kWakeConfigurations];
  TH1D* hybridModelToDataRatioPp[nWeightExponents][nJetPtBinsEEC][nTrackPtBinsEEC][HybridModelHistogramManager::kWakeConfigurations];
  TH1D* hybridModelToDataRatioPbPb[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][HybridModelHistogramManager::kWakeConfigurations];
  TH1D* hybridModelToDataRatioPbPbToPpRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][HybridModelHistogramManager::kWakeConfigurations];

   // Double ratio histograms for the Hybrid model
  TH1D* histogrammifiedHybridModelDoubleRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC][HybridModelHistogramManager::kWakeConfigurations];
  TH1D* hybridModelToDataRatioDoubleRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC][HybridModelHistogramManager::kWakeConfigurations];

  // Energy-energy correlators for perturbative calculations by Holguin et. al.
  TGraph* energyEnergyCorrelatorHolguinPbPb[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][HolguinHistogramManager::kMaxKValues];
  TGraph* energyEnergyCorrelatorHolguinPbPbToPpRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][HolguinHistogramManager::kMaxKValues];
  TH1D* histogrammifiedHolguinPbPb[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][HolguinHistogramManager::kMaxKValues];
  TH1D* histogrammifiedHolguinPbPbToPpRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][HolguinHistogramManager::kMaxKValues];
  TH1D* holguinToDataRatioPbPb[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][HolguinHistogramManager::kMaxKValues];
  TH1D* holguinToDataRatioPbPbToPpRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][HolguinHistogramManager::kMaxKValues];

  // Energy-energy correlators for CoLBT predictions
  TGraphErrors* energyEnergyCorrelatorCoLBTPbPbToPpRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][CoLBTHistogramManager::kQValues];
  TH1D* histogrammifiedCoLBTPbPbToPpRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][CoLBTHistogramManager::kQValues];
  TH1D* coLBTToDataRatioPbPbToPpRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][CoLBTHistogramManager::kQValues];

  // Energy-energy correlators for JEWEL predictions
  TH1D* energyEnergyCorrelatorJewelPp[nWeightExponents][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorJewelPbPb[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][JewelHistogramManager::kRecoilSettings];
  TH1D* energyEnergyCorrelatorJewelPbPbToPpRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][JewelHistogramManager::kRecoilSettings];

  TH1D* jewelToDataRatioPp[nWeightExponents][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* jewelToDataRatioPbPb[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][JewelHistogramManager::kRecoilSettings];
  TH1D* jewelToDataRatioPbPbToPpRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][JewelHistogramManager::kRecoilSettings];
  TH1D* jewelDoubleRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC][JewelHistogramManager::kRecoilSettings];
  TH1D* jewelToDataRatioDoubleRatio[nWeightExponents][nCentralityBins][nJetPtBinsEEC][JewelHistogramManager::kRecoilSettings];

  // Initialize histograms to NULL
  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        energyEnergyCorrelatorDoubleRatio[iWeightExponent][iCentrality][iJetPt] = NULL;
        systematicUncertaintyDoubleRatio[iWeightExponent][iCentrality][iJetPt] = NULL;

        for(int iWake = 0; iWake < HybridModelHistogramManager::kWakeConfigurations; iWake++){
          histogrammifiedHybridModelDoubleRatio[iWeightExponent][iCentrality][iJetPt][iWake] = NULL;
          hybridModelToDataRatioDoubleRatio[iWeightExponent][iCentrality][iJetPt][iWake] = NULL;
        }

        for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){
          jewelDoubleRatio[iWeightExponent][iCentrality][iJetPt][iRecoil] = NULL;
          jewelToDataRatioDoubleRatio[iWeightExponent][iCentrality][iJetPt][iRecoil] = NULL;
        }

        for(int iUncertainty = 0; iUncertainty < knRelativeUncertaintyTypes; iUncertainty++){
          hRelativeUncertaintyDoubleRatio[iWeightExponent][iCentrality][iJetPt][iUncertainty] = NULL;
        }
      }

      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        energyEnergyCorrelatorRawPp[iWeightExponent][iJetPt][iTrackPt] = NULL;
        energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt] = NULL;

        energyEnergyCorrelatorForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt] = NULL;
        systematicUncertaintyForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt] = NULL;

        energyEnergyCorrelatorJewelPp[iWeightExponent][iJetPt][iTrackPt] = NULL;
        jewelToDataRatioPp[iWeightExponent][iJetPt][iTrackPt] = NULL;

        for(int iUncertainty = 0; iUncertainty < knRelativeUncertaintyTypes; iUncertainty++){
          hRelativeUncertaintyPp[iWeightExponent][iJetPt][iTrackPt][iUncertainty] = NULL;
        }

        for(int iUncertainty = 0; iUncertainty < knSystematicUncertaintyTypes; iUncertainty++){
          systematicUncertaintyForPp[iWeightExponent][iUncertainty][iJetPt][iTrackPt] = NULL;
        } // Uncertainty type loop

        for(int iMCType = 0; iMCType < knPpMCTypes; iMCType++){
          hEnergyEnergyCorrelatorPpMC[iWeightExponent][iJetPt][iTrackPt][iMCType] = NULL;
          hEnergyEnergyCorrelatorPpMCToDataRatio[iWeightExponent][iJetPt][iTrackPt][iMCType] = NULL;
        }

        for(int iWake = 0; iWake < HybridModelHistogramManager::kWakeConfigurations; iWake++){
          energyEnergyCorrelatorHybridModelPp[iWeightExponent][iJetPt][iTrackPt][iWake] = NULL;
          histogrammifiedHybridModelPp[iWeightExponent][iJetPt][iTrackPt][iWake] = NULL;
          hybridModelToDataRatioPp[iWeightExponent][iJetPt][iTrackPt][iWake] = NULL;
        } // Wake type loop

        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          energyEnergyCorrelatorRawPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt] = NULL;
          energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt] = NULL;
          energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt] = NULL;

          energyEnergyCorrelatorForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt] = NULL;
          systematicUncertaintyForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt] = NULL;

          for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){
            energyEnergyCorrelatorJewelPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iRecoil] = NULL;
            energyEnergyCorrelatorJewelPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iRecoil] = NULL;
            jewelToDataRatioPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iRecoil] = NULL;
            jewelToDataRatioPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iRecoil] = NULL;
          }

          for(int iUncertainty = 0; iUncertainty < knRelativeUncertaintyTypes; iUncertainty++){
            hRelativeUncertaintyPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iUncertainty] = NULL;
            hRelativeUncertaintyPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iUncertainty] = NULL;
          }

          for(int iUncertainty = 0; iUncertainty < knSystematicUncertaintyTypes; iUncertainty++){
            systematicUncertaintyForPbPb[iWeightExponent][iUncertainty][iCentrality][iJetPt][iTrackPt] = NULL;
            systematicUncertaintyPbPbToPpRatio[iWeightExponent][iUncertainty][iCentrality][iJetPt][iTrackPt] = NULL;
          }

          for(int iWake = 0; iWake < HybridModelHistogramManager::kWakeConfigurations; iWake++){
            energyEnergyCorrelatorHybridModelPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = NULL;
            energyEnergyCorrelatorHybridModelPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = NULL;
            histogrammifiedHybridModelPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = NULL;
            histogrammifiedHybridModelPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = NULL;
            hybridModelToDataRatioPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = NULL;
            hybridModelToDataRatioPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = NULL;
          } // Wake type loop

          for(int iKValue = 0; iKValue < HolguinHistogramManager::kMaxKValues; iKValue++){
            energyEnergyCorrelatorHolguinPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue] = NULL;
            energyEnergyCorrelatorHolguinPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue] = NULL;
            histogrammifiedHolguinPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue] = NULL;
            histogrammifiedHolguinPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue] = NULL;
            holguinToDataRatioPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue] = NULL;
            holguinToDataRatioPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue] = NULL;
          }

          for(int iQValue = 0; iQValue < CoLBTHistogramManager::kQValues; iQValue++){
            energyEnergyCorrelatorCoLBTPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iQValue] = NULL;
            histogrammifiedCoLBTPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iQValue] = NULL;
            coLBTToDataRatioPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iQValue] = NULL;
          }

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
  double trackCorrelation;
  double epsilon = 0.0001;

  // Define helper histograms to determine the uncertainties relevant for the double ratio
  TH1D* uncertaintyTrackSelection;
  TH1D* uncertaintyBackgroundSubtraction;
  TH1D* uncertaintyTrackPairEfficiency;
  TH1D* uncertaintyMCnonclosure;
  TH1D* uncertaintySignalToBackgroundRatio;

  // Create a prime transformer for all your histogram manipulation needs
  AlgorithmLibrary* optimusPrimeTheTransformer = new AlgorithmLibrary();
  double binError;
  double otherBinError;

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

        // For relative uncertainties, read the total systematic uncertainties for pp collisions and clone the signal distribution
        hRelativeUncertaintyPp[iWeightExponent][iJetPt][iTrackPt][kRelativeUncertaintySystematic] = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty);
        hRelativeUncertaintyPp[iWeightExponent][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalUp] = (TH1D*) energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->Clone(Form("relativeStatisticalUncertaintyUpPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));
        hRelativeUncertaintyPp[iWeightExponent][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalDown] = (TH1D*) energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->Clone(Form("relativeStatisticalUncertaintyDownPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));

        // Transform the uncertainties into relative uncertainties
        for(int iUncertainty = 0; iUncertainty < knRelativeUncertaintyTypes; iUncertainty++){
          optimusPrimeTheTransformer->TransformToRelativeUncertainty(hRelativeUncertaintyPp[iWeightExponent][iJetPt][iTrackPt][iUncertainty], true);
        }

        // Move the statistical uncertainty bin contents up and down by the uncertainties
        for(int iBin = 1; iBin <= hRelativeUncertaintyPp[iWeightExponent][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalUp]->GetNbinsX(); iBin++){
          binError = hRelativeUncertaintyPp[iWeightExponent][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalUp]->GetBinError(iBin);
          hRelativeUncertaintyPp[iWeightExponent][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalUp]->SetBinContent(iBin, 1 + binError);
          hRelativeUncertaintyPp[iWeightExponent][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalDown]->SetBinContent(iBin, 1 - binError);
        }

        if(drawDoubleRatioDataToTheoryComparison){

          // Read the relevant systematic uncertainties for double ratio from pp
          uncertaintyBackgroundSubtraction = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kBackgroundSubtraction);
          uncertaintyTrackPairEfficiency = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kTrackPairEfficiency);
          uncertaintyMCnonclosure = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kMonteCarloNonClosure);
          uncertaintyTrackSelection = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty, SystematicUncertaintyOrganizer::kTrackSelection);

          // For double ratio, we can find the proper scaling factor for the tracking uncertainties by comparing integrals of raw
          // energy-energy correlators between different track pT bins. The higher track pT cuts are subsets of the lower cuts.
          // Thus integrals over the analysis region tell the degree of overlap of these regions.
          lowAnalysisBin = energyEnergyCorrelatorRawPp[iWeightExponent][iJetPt][iTrackPt]->GetXaxis()->FindBin(analysisDeltaR.first + epsilon);
          highAnalysisBin = energyEnergyCorrelatorRawPp[iWeightExponent][iJetPt][iTrackPt]->GetXaxis()->FindBin(analysisDeltaR.second - epsilon);
          lowPtIntegral = energyEnergyCorrelatorRawPp[iWeightExponent][iJetPt][trackPtBinsForDoubleRatio.first]->Integral(lowAnalysisBin, highAnalysisBin, "width");
          highPtIntegral = energyEnergyCorrelatorRawPp[iWeightExponent][iJetPt][trackPtBinsForDoubleRatio.second]->Integral(lowAnalysisBin, highAnalysisBin, "width");
          trackCorrelation = 1 - highPtIntegral / lowPtIntegral;

          // Now for the double ratio, we can combine the relevant uncertainties while scaling down the tracking related ones
          // by the expected overlap.
          systematicUncertaintyForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt] = (TH1D*) uncertaintyBackgroundSubtraction->Clone(Form("doubleRatioUncertaintyFromPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));
          for(int iBin = 1; iBin <= uncertaintyBackgroundSubtraction->GetNbinsX(); iBin++){
            systematicUncertaintyForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt]->SetBinError(iBin, TMath::Sqrt(TMath::Power(uncertaintyBackgroundSubtraction->GetBinError(iBin),2) + TMath::Power(uncertaintyTrackPairEfficiency->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintyTrackSelection->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintyMCnonclosure->GetBinError(iBin),2)));
          }

        } // Histograms relevant to double ratio


        for(auto centralityBin : drawnCentralityBin){

          iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);
          iCentralityMatched = uncertaintyCard[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);

          // Read the PbPb histograms
          energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt] = histograms[kPbPb][iWeightExponent]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);
          systematicUncertaintyForPbPb[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt] = uncertainties[kPbPb][iWeightExponent]->GetUncorrelatedSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty);
          systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt] = uncertainties[kPbPb][iWeightExponent]->GetCorrelatedSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty);

          // For relative uncertainties, read the total systematic uncertainties for PbPb collisions and clone the signal distribution
          hRelativeUncertaintyPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintySystematic] = uncertainties[kPbPb][iWeightExponent]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty);
          hRelativeUncertaintyPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalUp] = (TH1D*) energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Clone(Form("relativeStatisticalUncertaintyUpPbPb%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          hRelativeUncertaintyPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalDown] = (TH1D*) energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Clone(Form("relativeStatisticalUncertaintyDownPbPb%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));

          // Transform the uncertainties into relative uncertainties
          for(int iUncertainty = 0; iUncertainty < knRelativeUncertaintyTypes; iUncertainty++){
            optimusPrimeTheTransformer->TransformToRelativeUncertainty(hRelativeUncertaintyPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iUncertainty], true);
          }

          // Move the statistical uncertainty bin contents up and down by the uncertainties
          for(int iBin = 1; iBin <= hRelativeUncertaintyPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalUp]->GetNbinsX(); iBin++){
            binError = hRelativeUncertaintyPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalUp]->GetBinError(iBin);
            hRelativeUncertaintyPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalUp]->SetBinContent(iBin, 1 + binError);
            hRelativeUncertaintyPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalDown]->SetBinContent(iBin, 1 - binError);
          }

          if(drawDoubleRatioDataToTheoryComparison){

            // Read the relevant systematic uncertainties for double ratio from PbPb
            uncertaintyBackgroundSubtraction = uncertainties[kPbPb][iWeightExponent]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kBackgroundSubtraction);
            uncertaintyTrackPairEfficiency = uncertainties[kPbPb][iWeightExponent]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kTrackPairEfficiency);
            uncertaintyTrackSelection = uncertainties[kPbPb][iWeightExponent]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kTrackSelection);
            uncertaintyMCnonclosure = uncertainties[kPbPb][iWeightExponent]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kMonteCarloNonClosure);
            uncertaintySignalToBackgroundRatio = uncertainties[kPbPb][iWeightExponent]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty, SystematicUncertaintyOrganizer::kSignalToBackgroundRatio);

            // For double ratio, we can find the proper scaling factor for the tracking uncertainties by comparing integrals of raw
            // energy-energy correlators between different track pT bins. The higher track pT cuts are subsets of the lower cuts.
            // Thus integrals over the analysis region tell the degree of overlap of these regions.
            lowAnalysisBin = energyEnergyCorrelatorRawPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(analysisDeltaR.first + epsilon);
            highAnalysisBin = energyEnergyCorrelatorRawPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(analysisDeltaR.second - epsilon);
            lowPtIntegral = energyEnergyCorrelatorRawPbPb[iWeightExponent][iCentrality][iJetPt][trackPtBinsForDoubleRatio.first]->Integral(lowAnalysisBin, highAnalysisBin, "width");
            highPtIntegral = energyEnergyCorrelatorRawPbPb[iWeightExponent][iCentrality][iJetPt][trackPtBinsForDoubleRatio.second]->Integral(lowAnalysisBin, highAnalysisBin, "width");
            trackCorrelation = 1 - highPtIntegral / lowPtIntegral;

            // Now for the double ratio, we can combine the relevant uncertainties while scaling down the tracking related ones
            // by the expected overlap.
            systematicUncertaintyForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt] = (TH1D*) uncertaintyBackgroundSubtraction->Clone(Form("doubleRatioUncertaintyFromPbPb%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
            for(int iBin = 1; iBin <= uncertaintyBackgroundSubtraction->GetNbinsX(); iBin++){
              systematicUncertaintyForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetBinError(iBin, TMath::Sqrt(TMath::Power(uncertaintyBackgroundSubtraction->GetBinError(iBin),2) + TMath::Power(uncertaintyTrackPairEfficiency->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintyTrackSelection->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintySignalToBackgroundRatio->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintyMCnonclosure->GetBinError(iBin),2)));
            }

          } // Doing stuff relevant to double ratio
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

        if(drawDoubleRatioDataToTheoryComparison){
          systematicUncertaintyForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt]->Scale(energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->GetBinContent(10) / systematicUncertaintyForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt]->GetBinContent(10));
        } // Stuff relevant for double ratios

        for(auto centralityBin : drawnCentralityBin){

          iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);

          energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Scale(1.0 / energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
          systematicUncertaintyForPbPb[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Scale(energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->GetBinContent(10) / systematicUncertaintyForPbPb[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetBinContent(10));
          systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Scale(energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->GetBinContent(10) / systematicUncertaintyForPbPb[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetBinContent(10));
          
          if(drawDoubleRatioDataToTheoryComparison){
            systematicUncertaintyForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Scale(energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->GetBinContent(10) / systematicUncertaintyForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->GetBinContent(10));
          } // Studd relevant for double ratios

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

  // Reduce the statistical uncertainties by the overlapping statistics fraction for the double ratios
  if(drawDoubleRatioDataToTheoryComparison){
    for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
      for(auto jetPtBin : drawnJetPtBin){
        iJetPt = card[kPbPb][iWeightExponent]->FindBinIndexJetPtEEC(jetPtBin);
        for(auto trackPtBin : drawnTrackPtBin){
          iTrackPt = card[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);

          energyEnergyCorrelatorForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt] = (TH1D*) energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->Clone(Form("energyEnergyCorrelatorForDoubleRatioFromPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));

          for(int iBin = 1; iBin <= energyEnergyCorrelatorForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt]->GetNbinsX(); iBin++){
            energyEnergyCorrelatorForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt]->SetBinError(iBin, energyEnergyCorrelatorForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt]->GetBinError(iBin) * trackCorrelation);
          }

          for(auto centralityBin : drawnCentralityBin){
            iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);

            energyEnergyCorrelatorForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt] = (TH1D*) energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Clone(Form("energyEnergyCorrelatorForDoubleRatioFromPbPb%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));

            for(int iBin = 1; iBin <= energyEnergyCorrelatorForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->GetNbinsX(); iBin++){
              energyEnergyCorrelatorForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->SetBinError(iBin, energyEnergyCorrelatorForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->GetBinError(iBin) * trackCorrelation);
            }

          } // Centrality loop
        } // Track pT loop
      } // Jet pT loop
    } // Weight exponent loop
  

    // After regular ratios have been calculated, proceed to calculating double ratio. Here we need to use different histograms as above to properly take into account systematic uncertainty cancellation due to correlated datasets.
    for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
      for(auto centralityBin : drawnCentralityBin){
        iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);
        for(auto jetPtBin : drawnJetPtBin){
        iJetPt = card[kPbPb][iWeightExponent]->FindBinIndexJetPtEEC(jetPtBin);
          for(auto trackPtBin : drawnTrackPtBin){
            iTrackPt = card[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);

            // Calculate the single ratios with properly handled double ratio uncertainties
            energyEnergyCorrelatorForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Divide(energyEnergyCorrelatorForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt]);
            systematicUncertaintyForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Divide(systematicUncertaintyForDoubleRatioFromPp[iWeightExponent][iJetPt][iTrackPt]);
          } // Track pT loop
          // Calculate the double ratios from the single ratios with properly handled uncertainties
          energyEnergyCorrelatorDoubleRatio[iWeightExponent][iCentrality][iJetPt] = (TH1D*) energyEnergyCorrelatorForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][trackPtBinsForDoubleRatio.second]->Clone(Form("energyEnergyCorrelatorDoubleRatio%d%d%d", iWeightExponent, iCentrality, iJetPt));
          energyEnergyCorrelatorDoubleRatio[iWeightExponent][iCentrality][iJetPt]->Divide(energyEnergyCorrelatorForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][trackPtBinsForDoubleRatio.first]);
          systematicUncertaintyDoubleRatio[iWeightExponent][iCentrality][iJetPt] = (TH1D*) systematicUncertaintyForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][trackPtBinsForDoubleRatio.second]->Clone(Form("systematicUncertaintyDoubleRatio%d%d%d", iWeightExponent, iCentrality, iJetPt));
          systematicUncertaintyDoubleRatio[iWeightExponent][iCentrality][iJetPt]->Divide(systematicUncertaintyForDoubleRatioFromPbPb[iWeightExponent][iCentrality][iJetPt][trackPtBinsForDoubleRatio.first]);
        } // Jet pT loop
      } // Centrality llop
    } // Weight exponent loop

  } // Reducing statistical uncertainty by overlapping statistics for double ratios and then calculating them

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

  // Create histograms for relative data uncertainties for PbPb to pp ratio
  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    for(auto jetPtBin : drawnJetPtBin){
      iJetPt = card[kPbPb][iWeightExponent]->FindBinIndexJetPtEEC(jetPtBin);
      for(auto trackPtBin : drawnTrackPtBin){
        iTrackPt = card[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);
        for(auto centralityBin : drawnCentralityBin){
          iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);

          // For relative systematic uncertainty, add the correlated and uncorrelated systematic uncertainties in quadrature
          hRelativeUncertaintyPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintySystematic] = (TH1D*) systematicUncertaintyPbPbToPpRatio[iWeightExponent][kCorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("relativeSystematicUncertaintyForPbPbToPpRatio%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));

          for(int iBin = 1; iBin <= hRelativeUncertaintyPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintySystematic]->GetNbinsX(); iBin++){
            binError = hRelativeUncertaintyPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintySystematic]->GetBinError(iBin);
            otherBinError = systematicUncertaintyPbPbToPpRatio[iWeightExponent][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetBinError(iBin);
            hRelativeUncertaintyPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintySystematic]->SetBinError(iBin, TMath::Sqrt(binError*binError + otherBinError*otherBinError));
          }

          // Relative statistical uncertainties we can directly get from the ratio
          hRelativeUncertaintyPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalUp] = (TH1D*) energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Clone(Form("relativeStatisticalUncertaintyUpForPbPbToPpRatio%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          hRelativeUncertaintyPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalDown] = (TH1D*) energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Clone(Form("relativeStatisticalUncertaintyDownForPbPbToPpRatio%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));

          // Transform the uncertainties into relative uncertainties
          for(int iUncertainty = 0; iUncertainty < knRelativeUncertaintyTypes; iUncertainty++){
            optimusPrimeTheTransformer->TransformToRelativeUncertainty(hRelativeUncertaintyPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iUncertainty], true);
          }

          // Move the statistical uncertainty bin contents up and down by the uncertainties
          for(int iBin = 1; iBin <= hRelativeUncertaintyPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalUp]->GetNbinsX(); iBin++){
            binError = hRelativeUncertaintyPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalUp]->GetBinError(iBin);
            hRelativeUncertaintyPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalUp]->SetBinContent(iBin, 1 + binError);
            hRelativeUncertaintyPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalDown]->SetBinContent(iBin, 1 - binError);
          }

        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // Weight exponent loop

  // Calculate the relative uncertainties also for double ratios
  if(drawDoubleRatioDataToTheoryComparison){
    for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
      for(auto centralityBin : drawnCentralityBin){
        iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);
        for(auto jetPtBin : drawnJetPtBin){
          iJetPt = card[kPbPb][iWeightExponent]->FindBinIndexJetPtEEC(jetPtBin);

          // For relative uncertainties, read the total systematic uncertainties
          hRelativeUncertaintyDoubleRatio[iWeightExponent][iCentrality][iJetPt][kRelativeUncertaintySystematic] = (TH1D*) systematicUncertaintyDoubleRatio[iWeightExponent][iCentrality][iJetPt]->Clone(Form("relativeSystematicUncertaintyDoubleRatio%d%d%d", iWeightExponent, iCentrality, iJetPt));
          hRelativeUncertaintyDoubleRatio[iWeightExponent][iCentrality][iJetPt][kRelativeUncertaintyStatisticalUp] = (TH1D*) energyEnergyCorrelatorDoubleRatio[iWeightExponent][iCentrality][iJetPt]->Clone(Form("relativeStatisticalUncertaintyUpDoubleRatio%d%d%d", iWeightExponent, iCentrality, iJetPt));
          hRelativeUncertaintyDoubleRatio[iWeightExponent][iCentrality][iJetPt][kRelativeUncertaintyStatisticalDown] = (TH1D*) energyEnergyCorrelatorDoubleRatio[iWeightExponent][iCentrality][iJetPt]->Clone(Form("relativeStatisticalUncertaintyDownDoubleRatio%d%d%d", iWeightExponent, iCentrality, iJetPt));

          // Transform the uncertainties into relative uncertainties
          for(int iUncertainty = 0; iUncertainty < knRelativeUncertaintyTypes; iUncertainty++){
            optimusPrimeTheTransformer->TransformToRelativeUncertainty(hRelativeUncertaintyDoubleRatio[iWeightExponent][iCentrality][iJetPt][iUncertainty], true);
          }

          // Move the statistical uncertainty bin contents up and down by the uncertainties
          for(int iBin = 1; iBin <= hRelativeUncertaintyDoubleRatio[iWeightExponent][iCentrality][iJetPt][kRelativeUncertaintyStatisticalUp]->GetNbinsX(); iBin++){
            binError = hRelativeUncertaintyDoubleRatio[iWeightExponent][iCentrality][iJetPt][kRelativeUncertaintyStatisticalUp]->GetBinError(iBin);
            hRelativeUncertaintyDoubleRatio[iWeightExponent][iCentrality][iJetPt][kRelativeUncertaintyStatisticalUp]->SetBinContent(iBin, 1 + binError);
            hRelativeUncertaintyDoubleRatio[iWeightExponent][iCentrality][iJetPt][kRelativeUncertaintyStatisticalDown]->SetBinContent(iBin, 1 - binError);
          }
        } // Jet pT loop
      } // Centrality loop
    } // Weight exponent loop
  } // Relative uncertainties for double ratios

  // Read the histograms for pp MC
  TH1D* errorlessData;
  for(int iMCType = 0; iMCType < knPpMCTypes; iMCType++){
    for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
      for(auto jetPtBin : drawnJetPtBin){
      iJetPt = card[kPbPb][iWeightExponent]->FindBinIndexJetPtEEC(jetPtBin);
      iJetPtMatchedPp = ppMCCard[iMCType][iWeightExponent]->FindBinIndexJetPtEEC(jetPtBin);
        for(auto trackPtBin : drawnTrackPtBin){
          iTrackPt = card[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);
          iTrackPtMatchedPp = ppMCCard[iMCType][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);

          // Ask the histogram managers to deliver the proper histograms
          hEnergyEnergyCorrelatorPpMC[iWeightExponent][iJetPt][iTrackPt][iMCType] = ppMCHistograms[iMCType][iWeightExponent]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatchedPp, iTrackPtMatchedPp, EECHistograms::kSameJetPair);

          // Normalize the histogram to the analysis region
          lowAnalysisBin = hEnergyEnergyCorrelatorPpMC[iWeightExponent][iJetPt][iTrackPt][iMCType]->GetXaxis()->FindBin(analysisDeltaR.first + epsilon);
          highAnalysisBin = hEnergyEnergyCorrelatorPpMC[iWeightExponent][iJetPt][iTrackPt][iMCType]->GetXaxis()->FindBin(analysisDeltaR.second - epsilon);
          hEnergyEnergyCorrelatorPpMC[iWeightExponent][iJetPt][iTrackPt][iMCType]->Scale(1.0 / hEnergyEnergyCorrelatorPpMC[iWeightExponent][iJetPt][iTrackPt][iMCType]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

          // Take a ratio between data distribution and MC distribution disregarding the uncertainty in the data
          hEnergyEnergyCorrelatorPpMCToDataRatio[iWeightExponent][iJetPt][iTrackPt][iMCType] = (TH1D*) hEnergyEnergyCorrelatorPpMC[iWeightExponent][iJetPt][iTrackPt][iMCType]->Clone(Form("monteCarloRatioToDataPp%d%d%d%d", iWeightExponent, iJetPt, iTrackPt, iMCType));
          errorlessData = (TH1D*) energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->Clone(Form("errorlessComparisonToMC%d%d%d%d", iWeightExponent, iJetPt, iTrackPt, iMCType));
          optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
          hEnergyEnergyCorrelatorPpMCToDataRatio[iWeightExponent][iJetPt][iTrackPt][iMCType]->Divide(errorlessData);

        } // Track pT loop
      } // Jet pT loop
    } // Weight exponent loop
  } // pp MC type loop

  // Then read the histograms for the hybrid model
  HybridModelHistogramManager* hybridHistograms = new HybridModelHistogramManager(hybridModelFolder);
  double weightExponentHybrid;
  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    weightExponentHybrid = iWeightExponent+1;
    for(auto jetPtBin : drawnJetPtBin){
      iJetPt = card[kPbPb][iWeightExponent]->FindBinIndexJetPtEEC(jetPtBin);
      for(auto trackPtBin : drawnTrackPtBin){
        iTrackPt = card[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);
        for(int iWake = 0; iWake < HybridModelHistogramManager::kWakeConfigurations; iWake++){

          // Load the pp histograms
          energyEnergyCorrelatorHybridModelPp[iWeightExponent][iJetPt][iTrackPt][iWake] = hybridHistograms->GetEnergyEnergyCorrelatorPp(jetPtBin, trackPtBin, weightExponentHybrid, iWake);

          // Create a histogrammified version of the TGraphErrors
          histogrammifiedHybridModelPp[iWeightExponent][iJetPt][iTrackPt][iWake] = optimusPrimeTheTransformer->HistogrammifyWithErrors(energyEnergyCorrelatorHybridModelPp[iWeightExponent][iJetPt][iTrackPt][iWake], energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]);

          // Now we can take a ratio between hybrid model prediction and data
          hybridModelToDataRatioPp[iWeightExponent][iJetPt][iTrackPt][iWake] = (TH1D*) histogrammifiedHybridModelPp[iWeightExponent][iJetPt][iTrackPt][iWake]->Clone(Form("hybridModelRatioToDataPp%d%d%d%d", iWeightExponent, iJetPt, iTrackPt, iWake));
          errorlessData = (TH1D*) energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->Clone(Form("errorlessHybrid%d%d%d%d", iWeightExponent, iJetPt, iTrackPt, iWake));
          optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
          hybridModelToDataRatioPp[iWeightExponent][iJetPt][iTrackPt][iWake]->Divide(errorlessData);

          for(auto centralityBin : drawnCentralityBin){
            iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);

            // Load the PbPb and ratio histograms
            energyEnergyCorrelatorHybridModelPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = hybridHistograms->GetEnergyEnergyCorrelatorPbPb(centralityBin, jetPtBin, trackPtBin, weightExponentHybrid, iWake);
            energyEnergyCorrelatorHybridModelPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = hybridHistograms->GetEnergyEnergyCorrelatorPbPbToPpRatio(centralityBin, jetPtBin, trackPtBin, weightExponentHybrid, iWake);

            // Create a histogrammified versions of the TGraphErrors
            histogrammifiedHybridModelPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = optimusPrimeTheTransformer->HistogrammifyWithErrors(energyEnergyCorrelatorHybridModelPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake], energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]);
            histogrammifiedHybridModelPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = optimusPrimeTheTransformer->HistogrammifyWithErrors(energyEnergyCorrelatorHybridModelPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake], energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt]);

            // Now we can take a ratio between hybrid model prediction and data
            hybridModelToDataRatioPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = (TH1D*) histogrammifiedHybridModelPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake]->Clone(Form("hybridModelRatioToDataPbPb%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iWake));
            errorlessData = (TH1D*) energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Clone(Form("errorlessHybridPbPb%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iWake));
            optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
            hybridModelToDataRatioPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake]->Divide(errorlessData);

            hybridModelToDataRatioPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake] = (TH1D*) histogrammifiedHybridModelPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake]->Clone(Form("hybridModelRatioToDataPbPbToPpRatio%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iWake));
            errorlessData = (TH1D*) energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Clone(Form("errorlessHybridRatio%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iWake));
            optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
            hybridModelToDataRatioPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iWake]->Divide(errorlessData);
        
          } // Centrality loop
        } // Wake loop
      } // Track pT loop
    } // Jet pT loop
  } // Weight exponent loop

  // After all the Hybrid model histograms have been read, calculate double ratios from them
  if(drawDoubleRatioDataToTheoryComparison){
    for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
      weightExponentHybrid = iWeightExponent+1;
      for(auto centralityBin : drawnCentralityBin){
        iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);
        for(auto jetPtBin : drawnJetPtBin){
          iJetPt = card[kPbPb][iWeightExponent]->FindBinIndexJetPtEEC(jetPtBin);
          for(int iWake = 0; iWake < HybridModelHistogramManager::kWakeConfigurations; iWake++){
 
            // First, calculate the double ratio in Hybrid model
            histogrammifiedHybridModelDoubleRatio[iWeightExponent][iCentrality][iJetPt][iWake] = (TH1D*) histogrammifiedHybridModelPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][trackPtBinsForDoubleRatio.second][iWake]->Clone(Form("energyEnergyCorrelatorDoubleRatioHybrid%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iWake));
            histogrammifiedHybridModelDoubleRatio[iWeightExponent][iCentrality][iJetPt][iWake]->Divide(histogrammifiedHybridModelPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][trackPtBinsForDoubleRatio.first][iWake]);

            // Then, take a ratio between the double ratio in Hybrid model and double ratio in data
            hybridModelToDataRatioDoubleRatio[iWeightExponent][iCentrality][iJetPt][iWake] = (TH1D*) histogrammifiedHybridModelDoubleRatio[iWeightExponent][iCentrality][iJetPt][iWake]->Clone(Form("energyEnergyCorrelatorDoubleRatioDataHybridComparison%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iWake));
            errorlessData = (TH1D*) energyEnergyCorrelatorDoubleRatio[iWeightExponent][iCentrality][iJetPt]->Clone(Form("errorlessComparisonToHybridDoubleRatio%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iWake));
            optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
            hybridModelToDataRatioDoubleRatio[iWeightExponent][iCentrality][iJetPt][iWake]->Divide(errorlessData);

          } // Wake loop
        } // Jet pT loop
      } // Centrality loop
    } // Weight exponent loop
  } // Only calculate double ratios from Hybrid model if that is drawn

  // Helper veriables for normalizing the perturbative calculation
  double normalizationRegionDeltaRLow, normalizationRegionDeltaRHigh;
  double dataIntegralDistribution, theoryIntegralDistribution;
  double dataIntegralRatio, theoryIntegralRatio;
  int lowNormalizationBin = 13;
  int highNormalizationBin = 19;

  // Next, read the predictions from the perturbative calculation by Holguin et. al.
  HolguinHistogramManager* holguinHistograms = new HolguinHistogramManager(holguinDataFolder);
  int iKValue;
  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    for(auto jetPtBin : drawnJetPtBin){
      iJetPt = card[kPbPb][iWeightExponent]->FindBinIndexJetPtEEC(jetPtBin);
      for(auto trackPtBin : drawnTrackPtBin){
        iTrackPt = card[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);
        for(auto myKValue : holguinKValue){
          iKValue = holguinHistograms->FindKValueIndex(myKValue);
          if(iKValue < 0) continue;
          for(auto centralityBin : drawnCentralityBin){
            iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);

            // Load the PbPb and ratio histograms
            energyEnergyCorrelatorHolguinPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue] = holguinHistograms->GetEnergyEnergyCorrelatorPbPb(centralityBin, jetPtBin, trackPtBin, iWeightExponent+1, myKValue);
            energyEnergyCorrelatorHolguinPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue] = holguinHistograms->GetEnergyEnergyCorrelatorPbPbToPpRatio(centralityBin, jetPtBin, trackPtBin, iWeightExponent+1, myKValue);

            // There are only a few bins for which predictions exist. Before continuing, check if the histograms are NULL
            if(energyEnergyCorrelatorHolguinPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue] == NULL) continue;

            // Create a histogrammified versions of the TGraphs
            histogrammifiedHolguinPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue] = optimusPrimeTheTransformer->Histogrammify(energyEnergyCorrelatorHolguinPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue], energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]);
            histogrammifiedHolguinPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue] = optimusPrimeTheTransformer->Histogrammify(energyEnergyCorrelatorHolguinPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue], energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt]);

            // The absolute normalization for these theory predictions needs to be obtained from data
            // This is done by matching the integral in region 0.0419615 < DeltaR < 0.125882

            // Calculate the integrals in the data
            dataIntegralDistribution = energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Integral(lowNormalizationBin, highNormalizationBin, "width");
            dataIntegralRatio = energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Integral(lowNormalizationBin, highNormalizationBin, "width");

            // Calculate the integrals for the histogrammified distributions
            theoryIntegralDistribution = histogrammifiedHolguinPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue]->Integral(lowNormalizationBin, highNormalizationBin, "width");
            theoryIntegralRatio = histogrammifiedHolguinPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue]->Integral(lowNormalizationBin, highNormalizationBin, "width");

            // With these numbers, we can normalize the perturbative calculations
            energyEnergyCorrelatorHolguinPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue]->Scale(dataIntegralDistribution/theoryIntegralDistribution);
            histogrammifiedHolguinPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue]->Scale(dataIntegralDistribution/theoryIntegralDistribution);
            energyEnergyCorrelatorHolguinPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue]->Scale(dataIntegralRatio/theoryIntegralRatio);
            histogrammifiedHolguinPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue]->Scale(dataIntegralRatio/theoryIntegralRatio);

            // Now we can take a ratio between Holguin prediction and data
            holguinToDataRatioPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue] = (TH1D*) histogrammifiedHolguinPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue]->Clone(Form("holguinRatioToDataPbPb%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iKValue));
            errorlessData = (TH1D*) energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Clone(Form("errorlessHolguinPbPb%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iKValue));
            optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
            holguinToDataRatioPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue]->Divide(errorlessData);

            holguinToDataRatioPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue] = (TH1D*) histogrammifiedHolguinPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue]->Clone(Form("holguinRatioToDataPbPbToPpRatio%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iKValue));
            errorlessData = (TH1D*) energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Clone(Form("errorlessHolguinRatio%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iKValue));
            optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
            holguinToDataRatioPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iKValue]->Divide(errorlessData);
        
          } // Centrality loop
        } // Wake loop
      } // Track pT loop
    } // Jet pT loop
  } // Weight exponent loop

   // Then read the histograms for CoLBT model
  CoLBTHistogramManager* coLBTHistograms = new CoLBTHistogramManager(coLBTFolder);
  int iQValue;
  for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
    for(auto jetPtBin : drawnJetPtBin){
      iJetPt = card[kPbPb][iWeightExponent]->FindBinIndexJetPtEEC(jetPtBin);
      for(auto trackPtBin : drawnTrackPtBin){
        iTrackPt = card[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);
        for(auto centralityBin : drawnCentralityBin){
          iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);
          for(auto qValue : coLBTqValue){
            iQValue = coLBTHistograms->FindQValueIndex(qValue);
            if(iQValue < 0) continue;

            // Load the PbPb to pp ratio histograms
            energyEnergyCorrelatorCoLBTPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iQValue] = coLBTHistograms->GetEnergyEnergyCorrelatorPbPbToPpRatio(centralityBin, jetPtBin, trackPtBin, iWeightExponent+1, qValue);

            // There are only a few bins for which predictions exist. Before continuing, check if the histograms are NULL
            if(energyEnergyCorrelatorCoLBTPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iQValue] == NULL) continue;

            // Create a histogrammified versions of the TGraphErrors
            histogrammifiedCoLBTPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iQValue] = optimusPrimeTheTransformer->HistogrammifyWithErrors(energyEnergyCorrelatorCoLBTPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iQValue], energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt]);

            // Now we can take a ratio between CoLBT prediction and data
            coLBTToDataRatioPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iQValue] = (TH1D*) histogrammifiedCoLBTPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iQValue]->Clone(Form("coLBTRatioToDataPbPbToPpRatio%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iQValue));
            errorlessData = (TH1D*) energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Clone(Form("errorlessCoLBTRatio%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iQValue));
            optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
            coLBTToDataRatioPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iQValue]->Divide(errorlessData);
        
          } // Q-vaule loop
        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // Weight exponent loop

  // Read the histograms for JEWEL
  JewelHistogramManager* jewelHistograms = new JewelHistogramManager(jewelDataFolder);
  if(theoryComparisonIndex == 5 || theoryComparisonIndex == 6 || drawnPlots == 0){
    for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
      for(auto jetPtBin : drawnJetPtBin){
        iJetPt = card[kPbPb][iWeightExponent]->FindBinIndexJetPtEEC(jetPtBin);
        for(auto trackPtBin : drawnTrackPtBin){
          iTrackPt = card[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);

          // Load the pp histograms
          energyEnergyCorrelatorJewelPp[iWeightExponent][iJetPt][iTrackPt] = jewelHistograms->GetEnergyEnergyCorrelatorPp(jetPtBin, trackPtBin, iWeightExponent+1);

          // There are only a few bins for which predictions exist. Before continuing, check if the histograms are NULL
          if(energyEnergyCorrelatorJewelPp[iWeightExponent][iJetPt][iTrackPt] == NULL) continue;

          // Now we can take a ratio between JEWEL prediction and data
          jewelToDataRatioPp[iWeightExponent][iJetPt][iTrackPt] = (TH1D*) energyEnergyCorrelatorJewelPp[iWeightExponent][iJetPt][iTrackPt]->Clone(Form("jewelRatioToDataPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));
          errorlessData = (TH1D*) energyEnergyCorrelatorSignalPp[iWeightExponent][iJetPt][iTrackPt]->Clone(Form("errorlessJewel%d%d%d", iWeightExponent, iJetPt, iTrackPt));
          optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
          jewelToDataRatioPp[iWeightExponent][iJetPt][iTrackPt]->Divide(errorlessData);

          for(auto centralityBin : drawnCentralityBin){
            iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);

            for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){

              // Load the PbPb and ratio histograms
              energyEnergyCorrelatorJewelPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iRecoil] = jewelHistograms->GetEnergyEnergyCorrelatorPbPb(centralityBin, jetPtBin, trackPtBin, iWeightExponent+1, iRecoil);
              energyEnergyCorrelatorJewelPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iRecoil] = jewelHistograms->GetEnergyEnergyCorrelatorPbPbToPpRatio(centralityBin, jetPtBin, trackPtBin, iWeightExponent+1, iRecoil);

              // The prediction might not exist in all bins
              if(energyEnergyCorrelatorJewelPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iRecoil] == NULL) continue;

              // Take a ratio between JEWEL prediction and data
              jewelToDataRatioPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iRecoil] = (TH1D*) energyEnergyCorrelatorJewelPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iRecoil]->Clone(Form("jewelRatioToDataPbPb%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iRecoil));
              errorlessData = (TH1D*) energyEnergyCorrelatorSignalPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Clone(Form("errorlessJewelPbPb%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iRecoil));
              optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
              jewelToDataRatioPbPb[iWeightExponent][iCentrality][iJetPt][iTrackPt][iRecoil]->Divide(errorlessData);

              jewelToDataRatioPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iRecoil] = (TH1D*) energyEnergyCorrelatorJewelPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iRecoil]->Clone(Form("jewelRatioToDataPbPbToPpRatio%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iRecoil));
              errorlessData = (TH1D*) energyEnergyCorrelatorPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt]->Clone(Form("errorlessJewelRatio%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iRecoil));
              optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
              jewelToDataRatioPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][iTrackPt][iRecoil]->Divide(errorlessData);
            } // Recoil loop
          } // Centrality loop
        } // Track pT loop
      } // Jet pT loop
    } // Weight exponent loop

    // After all the JEWEL histograms have been read, calculate double ratios from them
    if(drawDoubleRatioDataToTheoryComparison){
      for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
        for(auto centralityBin : drawnCentralityBin){
          if(jewelHistograms->FindCentralityBinIndex(centralityBin) < 0) continue; // Skip bins tht contain no data
          iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);
          for(auto jetPtBin : drawnJetPtBin){
            iJetPt = card[kPbPb][iWeightExponent]->FindBinIndexJetPtEEC(jetPtBin);
            for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){
 
              // First, calculate the double ratio in JEWEL
              jewelDoubleRatio[iWeightExponent][iCentrality][iJetPt][iRecoil] = (TH1D*) energyEnergyCorrelatorJewelPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][trackPtBinsForDoubleRatio.second][iRecoil]->Clone(Form("energyEnergyCorrelatorDoubleRatioJewel%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iRecoil));
              jewelDoubleRatio[iWeightExponent][iCentrality][iJetPt][iRecoil]->Divide(energyEnergyCorrelatorJewelPbPbToPpRatio[iWeightExponent][iCentrality][iJetPt][trackPtBinsForDoubleRatio.first][iRecoil]);

              // Then, take a ratio between the double ratio in JEWEL and double ratio in data
              jewelToDataRatioDoubleRatio[iWeightExponent][iCentrality][iJetPt][iRecoil] = (TH1D*) jewelDoubleRatio[iWeightExponent][iCentrality][iJetPt][iRecoil]->Clone(Form("energyEnergyCorrelatorDoubleRatioDataJewelComparison%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iRecoil));
              errorlessData = (TH1D*) energyEnergyCorrelatorDoubleRatio[iWeightExponent][iCentrality][iJetPt]->Clone(Form("errorlessComparisonToJewelDoubleRatio%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iRecoil));
              optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
              jewelToDataRatioDoubleRatio[iWeightExponent][iCentrality][iJetPt][iRecoil]->Divide(errorlessData);

            } // Wake loop
          } // Jet pT loop
        } // Centrality loop
      } // Weight exponent loop
    } // Only calculate double ratios from Hybrid model if that is drawn
  }


  // ========================================= //
  //   Draw the histograms in separate plots   //
  // ========================================= //

  JDrawer* drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.05);
  drawer->SetLeftMargin(0.16);
  drawer->SetRightMargin(0.01);
  drawer->SetTopMargin(0.08);
  drawer->SetBottomMargin(0.09);

  // For published plots, apparently all the text needs to be bigger than default in JDrawer:
  drawer->SetTitleSizeX(40);
  drawer->SetTitleSizeY(35);
  drawer->SetLabelSizeX(25);
  drawer->SetLabelSizeY(25);

  // The offsets for titles and labels
  drawer->SetTitleOffsetY(1.25);
  drawer->SetTitleOffsetX(0.6);
  drawer->SetLabelOffsetY(0.006);

  // Logarithmic deltaR axis
  drawer->SetLogX(true); 

  drawer->SetTick(0,0);

  TString centralityString;
  TString compactCentralityString;
  TString jetPtString;
  TString compactJetPtString;
  TString trackPtString;
  TString compactTrackPtString;
  TString distributionString;

  // Common variables for different plots
  TLegend* legend;
  TLegend* anotherLegend;
  TLatex* mrLatexer = new TLatex();
  TLine* oneLine = new TLine(analysisDeltaR.first, 1, analysisDeltaR.second, 1);
  TBox* box;
  oneLine->SetLineColor(kBlack);
  oneLine->SetLineStyle(3);
  int styleIndex;
  double bottomRowScale, bottomPadMargin, leftPadMargin;
  double leftMarginAdder, bottomMarginAdder;
  double thisPadScale;

  // Legend position variables
  double legendX1, legendX2, legendY1, legendY2;

  // Legend text size for easy tuning
  double legendTextUpperCanvasSize = 0.07;
  double legendTextLowerCanvasSize = 0.1;

  // Adjust tick sizes for x-axis
  double tickSizeX = 0.05;

  // Normalization and style for theory predictions
  int color[9] = {kRed, kBlue, kGreen+3, kMagenta, kCyan, kOrange+7, kViolet-3, kPink-3, kOrange-3};

  // Change opacity for filled areas in order to make things visible in grauscale
  double jetWakeOpacity[HybridModelHistogramManager::kWakeConfigurations] = {0.4, 0.6, 0.5};
  double colbtOpacity[CoLBTHistogramManager::kQValues] = {0.4, 0.6};

  int jewelMarkerStyle[JewelHistogramManager::kRecoilSettings];
  jewelMarkerStyle[0] = kFullCross;
  jewelMarkerStyle[1] = kFullCrossX;

  // Colors for model comparisons in pp data
  int ppModelColor[knPpModelComparisons];
  ppModelColor[kPpCompareHybrid] = kRed;
  ppModelColor[kPpComparePythia] = kBlue;
  ppModelColor[kPpCompareHerwig] = kGreen+3;

  int ppMCIndexToModelColorMap[knPpMCTypes] = {1, 2};

  int ppMCMarkerStyle[knPpMCTypes];
  ppMCMarkerStyle[kPythia] = kFullDiamond;
  ppMCMarkerStyle[kHerwig] = kFullCircle;

  // Draw comparison of energy-energy correlator distributions between data and theory
  if(drawDistributionDataToTheoryComparison){

    for(auto jetPtBin : drawnJetPtBin){
      iJetPt = card[kPbPb][weightExponent-1]->FindBinIndexJetPtEEC(jetPtBin);

      jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", jetPtBin.first, jetPtBin.second);
      compactJetPtString = Form("_J=%.0f-%.0f", jetPtBin.first, jetPtBin.second);

      for(auto trackPtBin : drawnTrackPtBin){
        iTrackPt = card[kPbPb][weightExponent-1]->GetBinIndexTrackPtEEC(trackPtBin);

        trackPtString = Form("p_{T}^{ch} > %.0f GeV", trackPtBin);
        compactTrackPtString = Form("_T>%.1f", trackPtBin);
        compactTrackPtString.ReplaceAll(".","v");

        // Before going to centrality loop, draw the pp distribution

        // Create a new canvas for the plot
        drawer->CreateSplitCanvas();

        // Use logarithmic axis for EEC
        drawer->SetLogY(true);

        // Setup the legend for the plots
        legendX1 = 0.2; legendX2 = 0.5;
        legendY1 = 0.05; legendY2 = 0.55;

        // Make the legend slightly smaller in the main text highlight mode
        if(mainTextHighlightMode) legendY2 = 0.5;

        legend = new TLegend(legendX1, legendY1, legendX2, legendY2);
        legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(legendTextUpperCanvasSize); legend->SetTextFont(62);

        // Setup the another legend for the plots
        legendX1 = 0.5; legendX2 = 0.8;
        legendY1 = 0.05; legendY2 = 0.25;

        // Fow main text highlight mode, the legend needs to be more to the right and up
        if(mainTextHighlightMode){
          legendX1 = 0.58; legendY1 = 0.08;
        }

        anotherLegend = new TLegend(legendX1, legendY1, legendX2, legendY2);
        anotherLegend->SetFillStyle(0); anotherLegend->SetBorderSize(0); anotherLegend->SetTextSize(legendTextUpperCanvasSize); anotherLegend->SetTextFont(62);
        if(mainTextHighlightMode) anotherLegend->SetTextSize(legendTextUpperCanvasSize*2);

        if(theoryComparisonIndex > 5){
          anotherLegend->AddEntry((TObject*)0, energyWeightLegend[weightExponent-1].Data(), "");
          anotherLegend->AddEntry((TObject*)0, trackPtString.Data(), "");
        } else {
          if(mainTextHighlightMode){
            anotherLegend->AddEntry((TObject*)0, energyWeightLegend[weightExponent-1].Data(), "");
          } else {
            legend->AddEntry((TObject*)0, energyWeightLegend[weightExponent-1].Data(), "");
          }
          legend->AddEntry((TObject*)0, trackPtString.Data(), "");
        }

        // Set the drawing style for pp data and uncertainty histograms
        energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->SetMarkerStyle(kFullSquare);
        energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->SetMarkerColor(kBlack);

        systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetFillColorAlpha(kBlack, 0.4);
        systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetLineColor(kBlack);
        systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerColor(kBlack);
        systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->SetMarkerStyle(kFullSquare);

        systematicUncertaintyForPp[weightExponent-1][kCorrelatedUncertainty][iJetPt][iTrackPt]->SetFillColorAlpha(kBlack, 0.4);

        // Set the x-axis drawing range
        systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
        systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(distributionZoom[weightExponent-1].first, distributionZoom[weightExponent-1].second);

        // Adjust tick size for x-axis
        systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt]->GetXaxis()->SetTickLength(tickSizeX);

        // Draw the data correlator to upper canves
        drawer->DrawHistogramToUpperPad(systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt], "#Deltar", "EEC", " ", "e2");
        systematicUncertaintyForPp[weightExponent-1][kCorrelatedUncertainty][iJetPt][iTrackPt]->Draw("same,e3");
        energyEnergyCorrelatorSignalPp[weightExponent-1][iJetPt][iTrackPt]->Draw("same,p");

        legend->AddEntry(systematicUncertaintyForPp[weightExponent-1][kUncorrelatedUncertainty][iJetPt][iTrackPt], "pp data", "lpf");

        // Pythia+Herwig+Hybrid

        // Add legend for Pythia8 and Herwig7, unless we are drawing only JEWEL
        if(theoryComparisonIndex != 5){
          legend->AddEntry(hEnergyEnergyCorrelatorPpMC[weightExponent-1][iJetPt][iTrackPt][kPythia], "Pythia8 CP5", "pl");
          legend->AddEntry(hEnergyEnergyCorrelatorPpMC[weightExponent-1][iJetPt][iTrackPt][kHerwig], "Herwig7 CH3", "pl");
    
          // There is no wake in pp, so do only comparison with Hybrid without wake
          if(energyEnergyCorrelatorHybridModelPp[weightExponent-1][iJetPt][iTrackPt][0] != NULL){

            // Give some nice styles for the predictions
            energyEnergyCorrelatorHybridModelPp[weightExponent-1][iJetPt][iTrackPt][0]->SetLineColor(ppModelColor[kPpCompareHybrid]);
            energyEnergyCorrelatorHybridModelPp[weightExponent-1][iJetPt][iTrackPt][0]->SetLineWidth(0);
            energyEnergyCorrelatorHybridModelPp[weightExponent-1][iJetPt][iTrackPt][0]->SetMarkerColor(ppModelColor[kPpCompareHybrid]);
            energyEnergyCorrelatorHybridModelPp[weightExponent-1][iJetPt][iTrackPt][0]->SetMarkerStyle(kFullCircle);
            energyEnergyCorrelatorHybridModelPp[weightExponent-1][iJetPt][iTrackPt][0]->SetFillColorAlpha(ppModelColor[kPpCompareHybrid], 0.4);

            // Draw the prediction to the same canvas as the data
            energyEnergyCorrelatorHybridModelPp[weightExponent-1][iJetPt][iTrackPt][0]->Draw("3,same");

            // Add a legend for the theory prediction
            legend->AddEntry(energyEnergyCorrelatorHybridModelPp[weightExponent-1][iJetPt][iTrackPt][0], "Hybrid model", "f");
          }

          // After hybrid model, add also Pythia8 and Herwig7 predictions
          for(int iMCType = 0; iMCType < knPpMCTypes; iMCType++){

            // Give some nice styles for the predictions
            hEnergyEnergyCorrelatorPpMC[weightExponent-1][iJetPt][iTrackPt][iMCType]->SetLineColor(ppModelColor[ppMCIndexToModelColorMap[iMCType]]);
            hEnergyEnergyCorrelatorPpMC[weightExponent-1][iJetPt][iTrackPt][iMCType]->SetMarkerColor(ppModelColor[ppMCIndexToModelColorMap[iMCType]]);
            hEnergyEnergyCorrelatorPpMC[weightExponent-1][iJetPt][iTrackPt][iMCType]->SetMarkerStyle(ppMCMarkerStyle[iMCType]);
            hEnergyEnergyCorrelatorPpMC[weightExponent-1][iJetPt][iTrackPt][iMCType]->SetFillColorAlpha(ppModelColor[ppMCIndexToModelColorMap[iMCType]], 0.4);

            // Draw the prediction to the same canvas as the data
            hEnergyEnergyCorrelatorPpMC[weightExponent-1][iJetPt][iTrackPt][iMCType]->Draw("same,p");
          }
        }

        if(theoryComparisonIndex >= 5){
          // Option to also draw JEWEL in addition to other vacuum models

          // Set a nice style for the JEWEL histograms
          if(theoryComparisonIndex == 5){
            energyEnergyCorrelatorJewelPp[weightExponent-1][iJetPt][iTrackPt]->SetLineColor(kRed);
            energyEnergyCorrelatorJewelPp[weightExponent-1][iJetPt][iTrackPt]->SetMarkerColor(kRed);
            energyEnergyCorrelatorJewelPp[weightExponent-1][iJetPt][iTrackPt]->SetMarkerStyle(kFullCircle);
          } else {
            energyEnergyCorrelatorJewelPp[weightExponent-1][iJetPt][iTrackPt]->SetLineColor(kMagenta);
            energyEnergyCorrelatorJewelPp[weightExponent-1][iJetPt][iTrackPt]->SetMarkerColor(kMagenta);
            energyEnergyCorrelatorJewelPp[weightExponent-1][iJetPt][iTrackPt]->SetMarkerStyle(kFullCross);
          }

          // Draw the jewel histogram and add it to the legend
          energyEnergyCorrelatorJewelPp[weightExponent-1][iJetPt][iTrackPt]->Draw("same,pl");
          legend->AddEntry(energyEnergyCorrelatorJewelPp[weightExponent-1][iJetPt][iTrackPt], "JEWEL", "pl");
        }

        // Draw the legend to the upper pad
        legend->Draw();
        if(theoryComparisonIndex >= 5 || mainTextHighlightMode) anotherLegend->Draw();

        // Draw latex messages to the plots
        mrLatexer->SetTextFont(62);
        mrLatexer->SetTextSize(0.1);
        mrLatexer->DrawLatexNDC(0.16, 0.9, "CMS");

        if(addPreliminaryTag){
          mrLatexer->SetTextFont(52);
          mrLatexer->SetTextSize(0.075);
          mrLatexer->DrawLatexNDC(0.28, 0.9, "Preliminary");
        }

        // Luminosity
        mrLatexer->SetTextFont(42);
        mrLatexer->SetTextSize(0.065);
        mrLatexer->DrawLatexNDC(0.655, 0.9, "302 pb^{-1} pp (5.02 TeV)");

        // Binning
        mrLatexer->SetTextFont(62);
        mrLatexer->SetTextSize(0.07);
        mrLatexer->DrawLatexNDC(0.58, 0.775, jetPtString.Data());
        mrLatexer->DrawLatexNDC(0.712, 0.68, "anti-k_{T} R = 0.4");
        mrLatexer->DrawLatexNDC(0.795, 0.59, "|#eta_{jet}| < 1.6");

        // Linear scale for the ratio
        drawer->SetLogY(false);

        // Set the axis drawing ranges for ratio
        hRelativeUncertaintyPp[weightExponent-1][iJetPt][iTrackPt][kRelativeUncertaintySystematic]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
        hRelativeUncertaintyPp[weightExponent-1][iJetPt][iTrackPt][kRelativeUncertaintySystematic]->GetYaxis()->SetRangeUser(ratioZoomPpDistribution.first, ratioZoomPpDistribution.second);

        // Adjust tick sizes for the x-axis
        hRelativeUncertaintyPp[weightExponent-1][iJetPt][iTrackPt][kRelativeUncertaintySystematic]->GetXaxis()->SetTickLength(tickSizeX * 1.5);

        // Set the style for uncertainty bands for systematic and statistical uncertainties from data
        for(int iUncertainty = 0; iUncertainty < knRelativeUncertaintyTypes; iUncertainty++){
          hRelativeUncertaintyPp[weightExponent-1][iJetPt][iTrackPt][iUncertainty]->SetLineColor(kBlack);
          hRelativeUncertaintyPp[weightExponent-1][iJetPt][iTrackPt][iUncertainty]->SetLineStyle(9);
          hRelativeUncertaintyPp[weightExponent-1][iJetPt][iTrackPt][iUncertainty]->SetMarkerStyle(kFullCircle);
          hRelativeUncertaintyPp[weightExponent-1][iJetPt][iTrackPt][iUncertainty]->SetMarkerSize(0);
        }
        hRelativeUncertaintyPp[weightExponent-1][iJetPt][iTrackPt][kRelativeUncertaintySystematic]->SetFillColorAlpha(kBlack, 0.4);
        hRelativeUncertaintyPp[weightExponent-1][iJetPt][iTrackPt][kRelativeUncertaintySystematic]->SetLineWidth(0);

        // Create a new legend to show the different data uncertainty bands
        legendX1 = 0.25; legendX2 = 0.95;
        legendY1 = 0.86; legendY2 = 0.96;
        legend = new TLegend(legendX1, legendY1, legendX2, legendY2);
        legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(legendTextLowerCanvasSize); legend->SetTextFont(62);
        legend->SetNColumns(2);

        // Draw the error bars from data
        drawer->DrawHistogramToLowerPad(hRelativeUncertaintyPp[weightExponent-1][iJetPt][iTrackPt][kRelativeUncertaintySystematic], "#Deltar", "#frac{Theory}{Data}", " ", "e3");
        hRelativeUncertaintyPp[weightExponent-1][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalUp]->Draw("same,HIST,C");
        hRelativeUncertaintyPp[weightExponent-1][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalDown]->Draw("same,HIST,C");

        // Add the different uncertainties to the legend
        legend->AddEntry(hRelativeUncertaintyPp[weightExponent-1][iJetPt][iTrackPt][kRelativeUncertaintySystematic], "Data syst. unc.", "f");
        legend->AddEntry(hRelativeUncertaintyPp[weightExponent-1][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalUp], "Data stat. unc.", "l");

        // Draw Pythia, Herwig, and Hybrid unless we require only JEWEL to be drawn
        if(theoryComparisonIndex != 5){

          // Set the style for MC to data comparison histograms
          hybridModelToDataRatioPp[weightExponent-1][iJetPt][iTrackPt][0]->SetMarkerStyle(kFullCircle);
          hybridModelToDataRatioPp[weightExponent-1][iJetPt][iTrackPt][0]->SetMarkerSize(0);
          hybridModelToDataRatioPp[weightExponent-1][iJetPt][iTrackPt][0]->SetFillColorAlpha(ppModelColor[kPpCompareHybrid], 0.4);

          // Set the style for pp MC simulation to data ratios
          for(int iMCType = 0; iMCType < knPpMCTypes; iMCType++){
            hEnergyEnergyCorrelatorPpMCToDataRatio[weightExponent-1][iJetPt][iTrackPt][iMCType]->SetMarkerStyle(ppMCMarkerStyle[iMCType]);
            hEnergyEnergyCorrelatorPpMCToDataRatio[weightExponent-1][iJetPt][iTrackPt][iMCType]->SetMarkerColor(ppModelColor[ppMCIndexToModelColorMap[iMCType]]);
            hEnergyEnergyCorrelatorPpMCToDataRatio[weightExponent-1][iJetPt][iTrackPt][iMCType]->SetLineColor(ppModelColor[ppMCIndexToModelColorMap[iMCType]]);
            hEnergyEnergyCorrelatorPpMCToDataRatio[weightExponent-1][iJetPt][iTrackPt][iMCType]->SetFillColorAlpha(ppModelColor[ppMCIndexToModelColorMap[iMCType]], 0.4);
          }

          // Draw the ratio with respect to hybrid model
          hybridModelToDataRatioPp[weightExponent-1][iJetPt][iTrackPt][0]->Draw("same,e3");

          // Draw the ratio with respect to Pythia and Herwig simulations
          hEnergyEnergyCorrelatorPpMCToDataRatio[weightExponent-1][iJetPt][iTrackPt][kPythia]->Draw("same,p");
          hEnergyEnergyCorrelatorPpMCToDataRatio[weightExponent-1][iJetPt][iTrackPt][kHerwig]->Draw("same,p");

        }

        if(theoryComparisonIndex >= 5){
          // Theory comparison index above 5: add JEWEL vacuum

          // Set the style for the histograms
          if(theoryComparisonIndex == 5){
            jewelToDataRatioPp[weightExponent-1][iJetPt][iTrackPt]->SetMarkerStyle(kFullCircle);
            jewelToDataRatioPp[weightExponent-1][iJetPt][iTrackPt]->SetMarkerColor(kRed);
            jewelToDataRatioPp[weightExponent-1][iJetPt][iTrackPt]->SetLineColor(kRed);
          } else {
            jewelToDataRatioPp[weightExponent-1][iJetPt][iTrackPt]->SetMarkerStyle(kFullCross);
            jewelToDataRatioPp[weightExponent-1][iJetPt][iTrackPt]->SetMarkerColor(kMagenta);
            jewelToDataRatioPp[weightExponent-1][iJetPt][iTrackPt]->SetLineColor(kMagenta);
          }

          // Draw the histogram to the same canvas as data
          jewelToDataRatioPp[weightExponent-1][iJetPt][iTrackPt]->Draw("same,pl");
        }

        // Draw the legend to the lower pad
        legend->Draw();

        // Draw a line at one
        oneLine->Draw();

        distributionString = theoryComparisonIndex < 5 ? "Model" : "ExtendedModel";

        // If a plot name is given, save the plot in a file
        if(saveFigures) {
          gPad->GetCanvas()->SaveAs(Form("figures/energyEnergyCorrelator_distribution%sComparison%s_pp%s%s.%s", distributionString.Data(), saveComment.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
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
          legendX1 = 0.21; legendX2 = 0.51;
          legendY1 = 0.05; legendY2 = 0.6;
          legend = new TLegend(legendX1, legendY1, legendX2, legendY2);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(legendTextUpperCanvasSize); legend->SetTextFont(62);
          legend->AddEntry((TObject*)0, energyWeightLegend[weightExponent-1].Data(), "");
          legend->AddEntry((TObject*)0, trackPtString.Data(), "");

          // Set the drawing style for PbPb histogram
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(kFullSquare);
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(kBlack);
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetLineColor(kBlack);

          // Set drawing style for systematic uncertainties
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(kBlack, 0.4);
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetLineColor(kBlack);
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(kBlack);
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(kFullSquare);

          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(kGray+3, 0.4);
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(kBlue+4, 0.4);
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);

          // Set the axis drawing ranges
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(distributionZoom[weightExponent-1].first, distributionZoom[weightExponent-1].second);

          // Adjust tick size for x-axis
          systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTickLength(tickSizeX);

          // Draw the systematic uncertainties to the upper canvas
          drawer->DrawHistogramToUpperPad(systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], "#Deltar", "EEC", " ", "e2");
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          systematicUncertaintyForPbPb[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");

          // Draw the data points to the same canvas and add the histogram to the legend
          energyEnergyCorrelatorSignalPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt]->Draw("same,p");

          legend->AddEntry(systematicUncertaintyForPbPb[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], centralityString.Data(), "lpf");

          // Choose a set of theory predictions to draw to the figure
          if(theoryComparisonIndex == 0){
            // Theory comparison index 0, Hybrid model with different wake settings

            // Compare the prediction with and without wake
            for(int iWake = 0; iWake < HybridModelHistogramManager::kWakeConfigurations; iWake++){

              // There are some bins for which the prediction does not exist
              if(energyEnergyCorrelatorHybridModelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake] == NULL) continue;

              // Give some nice styles for the predictions
              energyEnergyCorrelatorHybridModelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetLineColor(color[iWake]);
              energyEnergyCorrelatorHybridModelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetMarkerColor(color[iWake]);
              energyEnergyCorrelatorHybridModelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetMarkerStyle(kFullCircle);
              energyEnergyCorrelatorHybridModelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetFillColorAlpha(color[iWake], jetWakeOpacity[iWake]);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorHybridModelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->Draw("3,same");

              // Add a legend for the theory prediction
              legend->AddEntry(energyEnergyCorrelatorHybridModelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake], hybridHistograms->GetWakeName(iWake), "f");
            }
          } else if(theoryComparisonIndex == 1){
            // Theory comparison index 1, perturbative calculations from Holguin et. al. with different k-values

            // Compare the prediction with different k-values
            styleIndex = 0;
            for(auto myKValue : holguinKValue){
              iKValue = holguinHistograms->FindKValueIndex(myKValue);

              // There are some bins for which the prediction does not exist
              if(energyEnergyCorrelatorHolguinPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue] == NULL) continue;

              // Give some nice styles for the predictions
              energyEnergyCorrelatorHolguinPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetLineColor(color[iKValue]);
              energyEnergyCorrelatorHolguinPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerColor(color[iKValue]);
              energyEnergyCorrelatorHolguinPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerStyle(markerStyleTheory[styleIndex++]);
              energyEnergyCorrelatorHolguinPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetFillColorAlpha(color[iKValue], 0.4);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorHolguinPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->Draw("pl,same");

              // Add a legend for the theory prediction
              legend->AddEntry(energyEnergyCorrelatorHolguinPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue], Form("Holguin, k=%.1f", myKValue), "pl");
            }
          } else if(theoryComparisonIndex == 5){
            // Theory comparison index 5, comparison to JEWEL

            // Compare the prediction with and without recoils
            for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){

              // There might be some bins for which the prediction does not exist
              if(energyEnergyCorrelatorJewelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil] == NULL) continue;

              // Give some nice styles for the predictions
              energyEnergyCorrelatorJewelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->SetLineColor(color[iRecoil]);
              energyEnergyCorrelatorJewelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->SetMarkerColor(color[iRecoil]);
              energyEnergyCorrelatorJewelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->SetMarkerStyle(jewelMarkerStyle[iRecoil]);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorJewelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->Draw("pl,same");

              // Add a legend for the theory prediction
              legend->AddEntry(energyEnergyCorrelatorJewelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil], jewelHistograms->GetRecoilName(iRecoil), "pl");
            }
          } else if(theoryComparisonIndex == 6){
            // Theory comparison index 6, comparison to best k-value from Holguin and JEWEL

            // Add the perturbative calculation with k=0.3
            iKValue = holguinHistograms->FindKValueIndex(0.3);

            // There are some bins for which the prediction does not exist
            if(energyEnergyCorrelatorHolguinPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue] != NULL){

              // Give some nice styles for the predictions
              energyEnergyCorrelatorHolguinPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetLineColor(color[2]);
              energyEnergyCorrelatorHolguinPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerColor(color[2]);
              energyEnergyCorrelatorHolguinPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerStyle(kFullCircle);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorHolguinPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->Draw("pl,same");

              // Add a legend for the theory prediction
              legend->AddEntry(energyEnergyCorrelatorHolguinPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue], "Holguin, k=0.3", "pl");
            }

            // Add the JEWEL predictions with and without recoil to the same plot

            for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){

              // There might be some bins for which the prediction does not exist
              if(energyEnergyCorrelatorJewelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil] == NULL) continue;

              // Give some nice styles for the predictions
              energyEnergyCorrelatorJewelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->SetLineColor(color[iRecoil]);
              energyEnergyCorrelatorJewelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->SetMarkerColor(color[iRecoil]);
              energyEnergyCorrelatorJewelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->SetMarkerStyle(jewelMarkerStyle[iRecoil]);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorJewelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->Draw("pl,same");

              // Add a legend for the theory prediction
              legend->AddEntry(energyEnergyCorrelatorJewelPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil], jewelHistograms->GetRecoilName(iRecoil), "pl");
            }
          }

          // Draw the legend to the upper pad
          legend->Draw();

          // Draw latex messages to the plots
          mrLatexer->SetTextFont(62);
          mrLatexer->SetTextSize(0.1);
          mrLatexer->DrawLatexNDC(0.16, 0.9, "CMS");

          if(addPreliminaryTag){
            mrLatexer->SetTextFont(52);
            mrLatexer->SetTextSize(0.075);
            mrLatexer->DrawLatexNDC(0.28, 0.9, "Preliminary");
          }
          
          // Luminosity
          mrLatexer->SetTextFont(42);
          mrLatexer->SetTextSize(0.065);
          mrLatexer->DrawLatexNDC(0.6, 0.9, "1.70 nb^{-1} PbPb (5.02 TeV)");

          // Binning
          mrLatexer->SetTextFont(62);
          mrLatexer->SetTextSize(0.07);
          mrLatexer->DrawLatexNDC(0.58, 0.775, jetPtString.Data());
          mrLatexer->DrawLatexNDC(0.712, 0.68, "anti-k_{T} R = 0.4");
          mrLatexer->DrawLatexNDC(0.795, 0.59, "|#eta_{jet}| < 1.6");

          // Linear scale for the ratio
          drawer->SetLogY(false);

          // Set the axis drawing ranges for ratio
          hRelativeUncertaintyPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintySystematic]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          hRelativeUncertaintyPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintySystematic]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

          // Adjust tick size for x-axis
          hRelativeUncertaintyPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintySystematic]->GetXaxis()->SetTickLength(tickSizeX*1.5);

          // Set the style for uncertainty bands for systematic and statistical uncertainties from data
          for(int iUncertainty = 0; iUncertainty < knRelativeUncertaintyTypes; iUncertainty++){
            hRelativeUncertaintyPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iUncertainty]->SetLineColor(kBlack);
            hRelativeUncertaintyPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iUncertainty]->SetLineStyle(9);
            hRelativeUncertaintyPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iUncertainty]->SetMarkerStyle(kFullCircle);
            hRelativeUncertaintyPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iUncertainty]->SetMarkerSize(0);
          }
          hRelativeUncertaintyPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintySystematic]->SetFillColorAlpha(kBlack, 0.4);
          hRelativeUncertaintyPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintySystematic]->SetLineWidth(0);


          // Create a new legend to show the different data uncertainty bands
          legendX1 = 0.25; legendX2 = 0.95;
          legendY1 = 0.86; legendY2 = 0.96;
          legend = new TLegend(legendX1, legendY1, legendX2, legendY2);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(legendTextLowerCanvasSize); legend->SetTextFont(62);
          legend->SetNColumns(2);

          // Draw the error bars from the data
          drawer->DrawHistogramToLowerPad(hRelativeUncertaintyPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintySystematic], "#Deltar", "#frac{Theory}{Data}", " ", "e3");
          hRelativeUncertaintyPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalUp]->Draw("same,HIST,C");
          hRelativeUncertaintyPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalDown]->Draw("same,HIST,C");

          // Add the different uncertainties to the legend
          legend->AddEntry(hRelativeUncertaintyPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintySystematic], "Data syst. unc.", "f");
          legend->AddEntry(hRelativeUncertaintyPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalUp], "Data stat. unc.", "l");

          // Choose a set of theory predictions to draw to the figure
          if(theoryComparisonIndex == 0){
            // Theory comparison index 0, Hybrid model with different wake settings
            for(int iWake = 0; iWake < HybridModelHistogramManager::kWakeConfigurations; iWake++){
              hybridModelToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetMarkerColor(color[iWake]);
              hybridModelToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetMarkerStyle(kFullCircle);
              hybridModelToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetMarkerSize(0);
              hybridModelToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetFillColorAlpha(color[iWake], jetWakeOpacity[iWake]);
              hybridModelToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->Draw("same,e3");
            }
          } else if(theoryComparisonIndex == 1){
            // Theory comparison index 1, perturbative calculations from Holguin et. al. with different k-values
            styleIndex = 0;
            for(auto myKValue : holguinKValue){
              iKValue = holguinHistograms->FindKValueIndex(myKValue);
              if(holguinToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue] == NULL) continue;
              holguinToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerColor(color[iKValue]);
              holguinToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetLineWidth(3);
              holguinToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetLineColor(color[iKValue]);
              holguinToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerStyle(markerStyleTheory[styleIndex++]);
              holguinToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->Draw("same,HIST,C");
            }
          } else if(theoryComparisonIndex == 5){
            // Theory comparison index 5, JEWEL prediction

            for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){
              if(jewelToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil] == NULL) continue;
              jewelToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->SetMarkerColor(color[iRecoil]);
              jewelToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->SetMarkerStyle(jewelMarkerStyle[iRecoil]);
              jewelToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->Draw("same,pl");
            }
          } else if(theoryComparisonIndex == 6){
            // Theory comparison index 6, best Holguin k-value and JEWEL prediction

            // Draw the perturbative calculation
            iKValue = holguinHistograms->FindKValueIndex(0.3);
            if(holguinToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue] != NULL){
              holguinToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerColor(color[2]);
              holguinToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetLineWidth(3);
              holguinToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetLineColor(color[2]);
              holguinToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerStyle(kFullCircle);
              holguinToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->Draw("same,HIST,C");
            }

            // Add the JEWEL predictions to the same figure
            for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){
              if(jewelToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil] == NULL) continue;
              jewelToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->SetMarkerColor(color[iRecoil]);
              jewelToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->SetMarkerStyle(jewelMarkerStyle[iRecoil]);
              jewelToDataRatioPbPb[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->Draw("same,pl");
            }
          }

          // Draw the legend to the lower pad
          legend->Draw();

          // Draw a line at one
          oneLine->Draw();

          // If a plot name is given, save the plot in a file
          if(saveFigures) {
            gPad->GetCanvas()->SaveAs(Form("figures/energyEnergyCorrelator_%sDistribution%s%s%s%s.%s", theorySaveName[theoryComparisonIndex].Data(), saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
          }

        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Comparing theory with the PbPb distribution

  // Draw data/theory comparison for PbPb to pp ratios
  if(drawRatioDataToTheoryComparison){

    // Change the zoom setting if we are drawing JEWEL, in order for points not to overlap with legends
    if(theoryComparisonIndex == 5){
      pbpbToPpRatioZoom.first = 0.1;
      pbpbToPpRatioZoom.second = 2.7;
      ratioZoom.first = 0.0;
      ratioZoom.second = 2.1;
    }

    if(theoryComparisonIndex == 6){
      pbpbToPpRatioZoom.first = 0.1;
      pbpbToPpRatioZoom.second = 2.3;
    }

    for(auto centralityBin : drawnCentralityBin){
      iCentrality = card[kPbPb][weightExponent-1]->FindBinIndexCentrality(centralityBin);

      centralityString = Form("PbPb %.0f-%.0f%%", centralityBin.first, centralityBin.second);
      compactCentralityString = Form("_C=%.0f-%.0f", centralityBin.first, centralityBin.second);

      for(auto jetPtBin : drawnJetPtBin){
        iJetPt = card[kPbPb][weightExponent-1]->FindBinIndexJetPtEEC(jetPtBin);

        jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", jetPtBin.first, jetPtBin.second);
        compactJetPtString = Form("_J=%.0f-%.0f", jetPtBin.first, jetPtBin.second);

        for(auto trackPtBin : drawnTrackPtBin){
          iTrackPt = card[kPbPb][weightExponent-1]->GetBinIndexTrackPtEEC(trackPtBin);

          trackPtString = Form("p_{T}^{ch} > %.0f GeV", trackPtBin);
          compactTrackPtString = Form("_T>%.1f", trackPtBin);
          compactTrackPtString.ReplaceAll(".","v");

          // Create a new canvas for the plot
          drawer->CreateSplitCanvas();

          // No logarithmic drawing for ratio
          drawer->SetLogY(false);

          // Setup the legend for plots
          legendX1 = 0.54; legendX2 = 0.84;
          legendY1 = 0.04; legendY2 = 0.22;
          if(theoryComparisonIndex == 2){
            legendX1 = 0.19; legendX2 = 0.49;
          } else if(theoryComparisonIndex == 5 && weightExponent == 2){
            legendX1 = 0.46; legendX2 = 0.8;
            legendY1 = 0.55; legendY2 = 0.46;
          } else if(theoryComparisonIndex >= 5){
            legendX1 = 0.58; legendX2 = 0.88;
            legendY1 = 0.54; legendY2 = 0.36;
          }
          legend = new TLegend(legendX1, legendY1, legendX2, legendY2);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(legendTextUpperCanvasSize); legend->SetTextFont(62);
          if(theoryComparisonIndex == 5 && weightExponent == 2){
            legend->AddEntry((TObject*)0, Form("%s, %s", trackPtString.Data(), energyWeightLegend[weightExponent-1].Data()), "");
          } else {
            legend->AddEntry((TObject*)0, trackPtString.Data(), "");
            legend->AddEntry((TObject*)0, energyWeightLegend[weightExponent-1].Data(), "");
          }

          // Make another legend to which all the different histograms are collected.
          legendX1 = 0.19; legendX2 = 0.49;
          legendY1 = 0.04; legendY2 = 0.31;
          if(theoryComparisonIndex == 2){
            legendX1 = 0.54; legendX2 = 0.84;
          } else if(theoryComparisonIndex >= 5){
            legendX1 = 0.5; legendX2 = 0.9;
            legendY1 = 0.56; legendY2 = 0.83;
          }
          anotherLegend = new TLegend(legendX1, legendY1, legendX2, legendY2);
          anotherLegend->SetFillStyle(0); anotherLegend->SetBorderSize(0); anotherLegend->SetTextSize(legendTextUpperCanvasSize); anotherLegend->SetTextFont(62);

          // Set the drawing style for PbPb to pp ratio histograms
          energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(kFullSquare);
          energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(kBlack);
          energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]->SetLineColor(kBlack);

          // Set drawing style for systematic uncertainties
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(kBlack, 0.4);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetLineColor(kBlack);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(kBlack);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(kFullSquare);

          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(kGray+3, 0.4);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(kBlue+4, 0.4);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);

          // Set the x-axis drawing range
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(pbpbToPpRatioZoom.first, pbpbToPpRatioZoom.second);

          // Adjust tick size for x-axis
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTickLength(tickSizeX);

          // Draw first the systematic uncertainties to the upper canves
          drawer->DrawHistogramToUpperPad(systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], "#Deltar", "#frac{PbPb}{pp}", " ", "e2");
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeUp][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");
          systematicUncertaintyPbPbToPpRatio[weightExponent-1][kCorrelatedUncertaintyShapeDown][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");

          // Then draw the PbPb to pp ratio and add a legend for it
          energyEnergyCorrelatorPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt]->Draw("same,p");
          anotherLegend->AddEntry(systematicUncertaintyPbPbToPpRatio[weightExponent-1][kUncorrelatedUncertainty][iCentrality][iJetPt][iTrackPt], Form("%s / pp", centralityString.Data()), "lpf");

          // Choose a set of theory predictions to draw to the figure
          if(theoryComparisonIndex == 0){
            // Theory comparison index 0, Hybrid model with different wake settings

            for(int iWake = 0; iWake < HybridModelHistogramManager::kWakeConfigurations; iWake++){

              // There are some bins for which the prediction does not exist
              if(energyEnergyCorrelatorHybridModelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake] == NULL) continue;

              // Give some nice styles for the predictions
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetLineColor(color[iWake]);
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetMarkerColor(color[iWake]);
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetMarkerStyle(kFullCircle);
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetFillColorAlpha(color[iWake], jetWakeOpacity[iWake]);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->Draw("3,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorHybridModelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake], hybridHistograms->GetWakeName(iWake), "f");
            }
          } else if (theoryComparisonIndex == 1){
            // Theory comparison index 1, perturbative calculations from Holguin et. al. with different k-values

            styleIndex = 0;
            for(auto myKValue : holguinKValue){
              iKValue = holguinHistograms->FindKValueIndex(myKValue);

              // There are some bins for which the prediction does not exist
              if(energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue] == NULL) continue;

              // Give some nice styles for the predictions
              energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetLineColor(color[iKValue]);
              energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerColor(color[iKValue]);
              energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerStyle(markerStyleTheory[styleIndex++]);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->Draw("pl,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue], Form("Holguin, k=%.1f", myKValue), "pl");
            }

          } else if(theoryComparisonIndex == 2){
            // Theory comparison index 2, CoLBT prediction

            for(auto myQValue : coLBTqValue){
              iQValue = coLBTHistograms->FindQValueIndex(myQValue);

              // There are some bins for which the prediction does not exist
              if(energyEnergyCorrelatorCoLBTPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue] == NULL) continue;

              // Give some nice styles for the predictions
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetLineColor(color[iQValue]);
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetMarkerColor(color[iQValue]);
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetMarkerStyle(kFullCircle);
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetFillColorAlpha(color[iQValue], colbtOpacity[iQValue]);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->Draw("3,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorCoLBTPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue], Form("CoLBT, q=%.1f", myQValue), "f");
            }
          } else if(theoryComparisonIndex == 3){
            // Theory comparison index 3, comparison between Hybrid model, Holguin calculation and CoLBT

            // First, set the style and draw the Hybrid model prediction with full wake

            // There are some bins for which the prediction does not exist
            if(energyEnergyCorrelatorHybridModelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][HybridModelHistogramManager::kFullWake] != NULL){

              // Give some nice styles for the predictions
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][HybridModelHistogramManager::kFullWake]->SetLineColor(color[2]);
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][HybridModelHistogramManager::kFullWake]->SetMarkerColor(color[2]);
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][HybridModelHistogramManager::kFullWake]->SetMarkerStyle(kFullCircle);
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][HybridModelHistogramManager::kFullWake]->SetFillColorAlpha(color[2], 0.4);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][HybridModelHistogramManager::kFullWake]->Draw("3,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorHybridModelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][HybridModelHistogramManager::kFullWake], hybridHistograms->GetWakeName(HybridModelHistogramManager::kFullWake), "f");
              }

            // Next, add the perturbative calculation with k=0.3
            iKValue = holguinHistograms->FindKValueIndex(0.3);

            // There are some bins for which the prediction does not exist
            if(energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue] != NULL){

              // Give some nice styles for the predictions
              energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetLineColor(color[1]);
              energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerColor(color[1]);
              energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerStyle(kFullCircle);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->Draw("pl,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue], "Holguin, k=0.3", "pl");
            }

            // Finally, add CoLBT prediction with q=1
            iQValue = coLBTHistograms->FindQValueIndex(1);

            // There are some bins for which the prediction does not exist
            if(energyEnergyCorrelatorCoLBTPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue] != NULL){

              // Give some nice styles for the predictions
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetLineColor(color[0]);
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetMarkerColor(color[0]);
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetMarkerStyle(kFullCircle);
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetFillColorAlpha(color[0], 0.4);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->Draw("3,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorCoLBTPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue], "CoLBT, q=1", "f");
            }

          } else if (theoryComparisonIndex == 4){
            // Theory comparison index 4, comparison between best k in Holguin calculation and CoLBT

            // Add the perturbative calculation with k=0.3
            iKValue = holguinHistograms->FindKValueIndex(0.3);

            // There are some bins for which the prediction does not exist
            if(energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue] != NULL){

              // Give some nice styles for the predictions
              energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetLineColor(kMagenta);
              energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerColor(kMagenta);
              energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerStyle(kFullCircle);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->Draw("pl,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue], "Holguin, k=0.3", "pl");
            } 

            // Then add all different q-values for CoLBT
            for(auto myQValue : coLBTqValue){
              iQValue = coLBTHistograms->FindQValueIndex(myQValue);

              // There are some bins for which the prediction does not exist
              if(energyEnergyCorrelatorCoLBTPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue] == NULL) continue;

              // Give some nice styles for the predictions
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetLineColor(color[iQValue]);
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetMarkerColor(color[iQValue]);
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetMarkerStyle(kFullCircle);
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetFillColorAlpha(color[iQValue], 0.4);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->Draw("3,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorCoLBTPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue], Form("CoLBT, q=%.1f", myQValue), "f");
            }

          } else if (theoryComparisonIndex == 5){
            // Theory comparison index 5, comparison to JEWEL

            // Compare distributions with and without recoil
            for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){

              // There might be some bins for which the prediction does not exist
              if(energyEnergyCorrelatorJewelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil] == NULL) continue;

              // Give some nice styles for the predictions
              energyEnergyCorrelatorJewelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->SetLineColor(color[iRecoil]);
              energyEnergyCorrelatorJewelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->SetMarkerColor(color[iRecoil]);
              energyEnergyCorrelatorJewelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->SetMarkerStyle(jewelMarkerStyle[iRecoil]);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorJewelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->Draw("pl,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorJewelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil], jewelHistograms->GetRecoilName(iRecoil), "pl");
            }
          } else if (theoryComparisonIndex == 6){
            // Theory comparison index 6, comparison to best k-value from Holguin and JEWEL

            // Add the perturbative calculation with k=0.3
            iKValue = holguinHistograms->FindKValueIndex(0.3);

            // There are some bins for which the prediction does not exist
            if(energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue] != NULL){

              // Give some nice styles for the predictions
              energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetLineColor(kMagenta);
              energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerColor(kMagenta);
              energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerStyle(kFullCircle);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->Draw("pl,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorHolguinPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue], "Holguin, k=0.3", "pl");
            } 

            // Add the JEWEL predictions with and without recoil to the same plot
            for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){

              // There might be some bins for which the prediction does not exist
              if(energyEnergyCorrelatorJewelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil] == NULL) continue;

              // Give some nice styles for the predictions
              energyEnergyCorrelatorJewelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->SetLineColor(color[iRecoil]);
              energyEnergyCorrelatorJewelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->SetMarkerColor(color[iRecoil]);
              energyEnergyCorrelatorJewelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->SetMarkerStyle(jewelMarkerStyle[iRecoil]);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorJewelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->Draw("pl,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorJewelPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil], jewelHistograms->GetRecoilName(iRecoil), "pl");
            }
          }

          // Draw the legends to the upper pad
          legend->Draw();
          anotherLegend->Draw();

          // Draw a line to one
          oneLine->Draw();

          // Draw latex messages to the plots
          mrLatexer->SetTextFont(62);
          mrLatexer->SetTextSize(0.09);

          // Need to move the CMS text if doing JEWEL comparisons, it will otherwise overlap with points
          if(theoryComparisonIndex >= 5 && !addPreliminaryTag){
            mrLatexer->DrawLatexNDC(0.16, 0.9, "CMS");
          } else if (theoryComparisonIndex == 2) {
            mrLatexer->DrawLatexNDC(0.2, 0.78, "CMS");
          } else {
            mrLatexer->DrawLatexNDC(0.19, 0.78, "CMS");
          }

          if(addPreliminaryTag){
            mrLatexer->SetTextFont(52);
            mrLatexer->SetTextSize(0.065);
            if(theoryComparisonIndex == 2){
              mrLatexer->DrawLatexNDC(0.31, 0.78, "Preliminary");
            } else {
              mrLatexer->DrawLatexNDC(0.3, 0.78, "Preliminary");
            }
          }

          // Luminosity
          mrLatexer->SetTextFont(42);
          mrLatexer->SetTextSize(0.06);
          mrLatexer->DrawLatexNDC(0.3, 0.9, "1.70 nb^{-1} PbPb (5.02 TeV) + 302 pb^{-1} pp (5.02 TeV)");
 
          // Binning
          mrLatexer->SetTextFont(62);
          mrLatexer->SetTextSize(0.07);
          if(theoryComparisonIndex == 2){
            mrLatexer->DrawLatexNDC(0.5, 0.79, jetPtString.Data());
            mrLatexer->DrawLatexNDC(0.632, 0.695, "anti-k_{T} R = 0.4");
            mrLatexer->DrawLatexNDC(0.715, 0.605, "|#eta_{jet}| < 1.6");
          } else if(theoryComparisonIndex == 5 && weightExponent == 2){
            mrLatexer->DrawLatexNDC(0.54, 0.38, jetPtString.Data());
            mrLatexer->DrawLatexNDC(0.22, 0.17, "anti-k_{T} R = 0.4");
            mrLatexer->DrawLatexNDC(0.22, 0.08, "|#eta_{jet}| < 1.6");
          } else if(theoryComparisonIndex >= 5){

            if(mainTextHighlightMode){
              mrLatexer->SetTextSize(0.11);
              mrLatexer->DrawLatexNDC(0.39, 0.08, jetPtString.Data());
              mrLatexer->SetTextSize(0.07);
              mrLatexer->DrawLatexNDC(0.2, 0.22, "anti-k_{T} R = 0.4");
              mrLatexer->DrawLatexNDC(0.2, 0.13, "|#eta_{jet}| < 1.6");
            } else {
              mrLatexer->DrawLatexNDC(0.54, 0.08, jetPtString.Data());
              mrLatexer->DrawLatexNDC(0.22, 0.17, "anti-k_{T} R = 0.4");
              mrLatexer->DrawLatexNDC(0.22, 0.08, "|#eta_{jet}| < 1.6");
            }
          } else {

            if(mainTextHighlightMode){
              mrLatexer->SetTextSize(0.11);
              mrLatexer->DrawLatexNDC(0.34, 0.75, jetPtString.Data());
              mrLatexer->SetTextSize(0.07);
              mrLatexer->DrawLatexNDC(0.642, 0.625, "anti-k_{T} R = 0.4");
              mrLatexer->DrawLatexNDC(0.642, 0.535, "|#eta_{jet}| < 1.6");
            } else {
              mrLatexer->DrawLatexNDC(0.51, 0.78, jetPtString.Data());
              mrLatexer->DrawLatexNDC(0.642, 0.685, "anti-k_{T} R = 0.4");
              mrLatexer->DrawLatexNDC(0.725, 0.595, "|#eta_{jet}| < 1.6");
            }
          }
          

          // Linear scale for the ratio
          drawer->SetLogY(false);

          // Create a new legend to show the different data uncertainty bands
          legendX1 = 0.25; legendX2 = 0.92;
          legendY1 = 0.86; legendY2 = 0.96;
          if(theoryComparisonIndex == 6){
            legendX1 = 0.18; legendX2 = 0.85;
            legendY1 = 0.3; legendY2 = 0.4;
          } else if(theoryComparisonIndex == 3 || theoryComparisonIndex == 4){
            legendX1 = 0.22; legendX2 = 0.89;
            legendY1 = 0.3; legendY2 = 0.4;
          } 
          legend = new TLegend(legendX1, legendY1, legendX2, legendY2);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(legendTextLowerCanvasSize); legend->SetTextFont(62);
          legend->SetNColumns(2);

          // Set the axis drawing ranges
          hRelativeUncertaintyPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintySystematic]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          hRelativeUncertaintyPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintySystematic]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

          // Adjust tick size for x-axis
          hRelativeUncertaintyPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintySystematic]->GetXaxis()->SetTickLength(tickSizeX*1.5);

          // Set the style for uncertainty bands for systematic and statistical uncertainties from data
          for(int iUncertainty = 0; iUncertainty < knRelativeUncertaintyTypes; iUncertainty++){
            hRelativeUncertaintyPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iUncertainty]->SetLineColor(kBlack);
            hRelativeUncertaintyPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iUncertainty]->SetLineStyle(9);
            hRelativeUncertaintyPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iUncertainty]->SetMarkerStyle(kFullCircle);
            hRelativeUncertaintyPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iUncertainty]->SetMarkerSize(0);
          }
          hRelativeUncertaintyPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintySystematic]->SetFillColorAlpha(kBlack, 0.4);
          hRelativeUncertaintyPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintySystematic]->SetLineWidth(0);


          // Draw the error bars from the data
          drawer->DrawHistogramToLowerPad(hRelativeUncertaintyPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintySystematic], "#Deltar", "#frac{Theory}{Data}", " ", "e3");
          hRelativeUncertaintyPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalUp]->Draw("same,HIST,C");
          hRelativeUncertaintyPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalDown]->Draw("same,HIST,C");

          // Add the different uncertainties to the legend
          legend->AddEntry(hRelativeUncertaintyPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintySystematic], "Data syst. unc.", "f");
          legend->AddEntry(hRelativeUncertaintyPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][kRelativeUncertaintyStatisticalUp], "Data stat. unc.", "l");

          // Choose a set of theory predictions to draw to the figure
          if(theoryComparisonIndex == 0){
            // Theory comparison index 0, Hybrid model with different wake settings
            for(int iWake = 0; iWake < HybridModelHistogramManager::kWakeConfigurations; iWake++){
              hybridModelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetMarkerColor(color[iWake]);
              hybridModelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetMarkerStyle(kFullCircle);
              hybridModelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetMarkerSize(0);
              hybridModelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->SetFillColorAlpha(color[iWake], jetWakeOpacity[iWake]);
              hybridModelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iWake]->Draw("same,e3");
            }
          } else if (theoryComparisonIndex == 1){
            // Theory comparison index 1, perturbative calculations from Holguin et. al. with different k-values
            styleIndex = 0;
            for(auto myKValue : holguinKValue){
              iKValue = holguinHistograms->FindKValueIndex(myKValue);
              if(holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue] == NULL) continue;
              holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerColor(color[iKValue]);
              holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetLineWidth(3);
              holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetLineColor(color[iKValue]);
              holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerStyle(markerStyleTheory[styleIndex++]);
              holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->Draw("same,HIST,C");
            }
          } else if(theoryComparisonIndex == 2){
            // Theory comparison index 2, CoLBT prediction
            for(auto myQValue : coLBTqValue){
              iQValue = coLBTHistograms->FindQValueIndex(myQValue);
              if(coLBTToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue] == NULL) continue;
              coLBTToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetMarkerColor(color[iQValue]);
              coLBTToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetMarkerStyle(kFullCircle);
              coLBTToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetMarkerSize(0);
              coLBTToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetFillColorAlpha(color[iQValue], colbtOpacity[iQValue]);
              coLBTToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->Draw("same,e3");
            }
          } else if(theoryComparisonIndex == 3){
            // Theory comparison index 3, comparison between Hybrid model, Holguin calculation and CoLBT

            // First, add the Hybrid model prediction
            if(hybridModelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][HybridModelHistogramManager::kFullWake] != NULL){
              hybridModelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][HybridModelHistogramManager::kFullWake]->SetMarkerColor(color[2]);
              hybridModelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][HybridModelHistogramManager::kFullWake]->SetMarkerStyle(kFullCircle);
              hybridModelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][HybridModelHistogramManager::kFullWake]->SetMarkerSize(0);
              hybridModelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][HybridModelHistogramManager::kFullWake]->SetFillColorAlpha(color[2], 0.4);
              hybridModelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][HybridModelHistogramManager::kFullWake]->Draw("same,e3");
            }

            // Then, add perturbative calculation by Holguin and friends
            iKValue = holguinHistograms->FindKValueIndex(0.3);
            if(holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue] != NULL){
              holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerColor(color[1]);
              holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetLineWidth(3);
              holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetLineColor(color[1]);
              holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerStyle(kFullCircle);
              holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->Draw("same,HIST,C");
            }

            // Finally, add the CoLBT prediction
            iQValue = coLBTHistograms->FindQValueIndex(1);
            if(coLBTToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue] != NULL){
              coLBTToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetMarkerColor(color[0]);
              coLBTToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetMarkerStyle(kFullCircle);
              coLBTToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetMarkerSize(0);
              coLBTToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetFillColorAlpha(color[0], 0.4);
              coLBTToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->Draw("same,e3");
            }

          } else if(theoryComparisonIndex == 4){
            // Theory comparison index 4, comparison between best k in Holguin calculation and CoLBT

            // Add the perturbative calculation by Holguin and friends
            iKValue = holguinHistograms->FindKValueIndex(0.3);
            if(holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue] != NULL){
              holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerColor(kMagenta);
              holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetLineWidth(3);
              holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetLineColor(kMagenta);
              holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerStyle(kFullCircle);
              holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->Draw("same,HIST,C");
            }

            // Then, add CoLBT prediction
            for(auto myQValue : coLBTqValue){
              iQValue = coLBTHistograms->FindQValueIndex(myQValue);
              if(coLBTToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue] == NULL) continue;
              coLBTToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetMarkerColor(color[iQValue]);
              coLBTToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetMarkerStyle(kFullCircle);
              coLBTToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetMarkerSize(0);
              coLBTToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->SetFillColorAlpha(color[iQValue], 0.4);
              coLBTToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iQValue]->Draw("same,e3");
            }

          } else if(theoryComparisonIndex == 5){
            // Theory comparison index 5, JEWEL prediction

            for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){
              if(jewelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil] == NULL) continue;
              jewelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->SetMarkerColor(color[iRecoil]);
              jewelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->SetMarkerStyle(jewelMarkerStyle[iRecoil]);
              jewelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->Draw("same,pl");
            }
          } else if(theoryComparisonIndex == 6){
            // Theory comparison index 6, best k-value from Holguin and JEWEL prediction

            // Add perturbative calculation by Holguin and friends
            iKValue = holguinHistograms->FindKValueIndex(0.3);
            if(holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue] != NULL){
              holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerColor(kMagenta);
              holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetLineWidth(3);
              holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetLineColor(kMagenta);
              holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->SetMarkerStyle(kFullCircle);
              holguinToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iKValue]->Draw("same,HIST,C");
            }

            for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){
              if(jewelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil] == NULL) continue;
              jewelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->SetMarkerColor(color[iRecoil]);
              jewelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->SetMarkerStyle(jewelMarkerStyle[iRecoil]);
              jewelToDataRatioPbPbToPpRatio[weightExponent-1][iCentrality][iJetPt][iTrackPt][iRecoil]->Draw("same,pl");
            }
          }

          // Draw the legend to the lower pad
          legend->Draw();

          // Draw a line to one
          oneLine->Draw();

          // If a plot name is given, save the plot in a file
          if(saveFigures) {
            gPad->GetCanvas()->SaveAs(Form("figures/energyEnergyCorrelator_%sRatio%s%s%s%s.%s", theorySaveName[theoryComparisonIndex].Data(), saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
          }

        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Comparing theory with the PbPb to pp rations

  // Draw the data/theory comparison for double ratios
  if(drawDoubleRatioDataToTheoryComparison){

    // We need to increase left margin for double ratio plots such that the title fits to screen
    drawer->SetTitleOffsetY(1.65);
    drawer->SetLeftMargin(0.22);

    for(auto centralityBin : drawnCentralityBin){
      iCentrality = card[kPbPb][weightExponent-1]->FindBinIndexCentrality(centralityBin);

      centralityString = Form("PbPb %.0f-%.0f%%", centralityBin.first, centralityBin.second);
      compactCentralityString = Form("_C=%.0f-%.0f", centralityBin.first, centralityBin.second);

      for(auto jetPtBin : drawnJetPtBin){
        iJetPt = card[kPbPb][weightExponent-1]->FindBinIndexJetPtEEC(jetPtBin);

        jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", jetPtBin.first, jetPtBin.second);
        compactJetPtString = Form("_J=%.0f-%.0f", jetPtBin.first, jetPtBin.second);

        // Create a new canvas for the plot
        drawer->CreateSplitCanvas();

        // No logarithmic drawing for ratio
        drawer->SetLogY(false);

        // Setup the legend for plots
        legendX1 = 0.6; legendX2 = 0.9;
        if(theoryComparisonIndex >= 5){
          legendX1 = 0.63; legendX2 = 0.93;
        }

        legend = new TLegend(legendX1, 0.05, legendX2, 0.2);
        legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(legendTextUpperCanvasSize); legend->SetTextFont(62);
        legend->AddEntry((TObject*)0, centralityString.Data(), "");
        legend->AddEntry((TObject*)0, energyWeightLegend[weightExponent-1].Data(), "");

        // Make another legend to which all the different histograms are collected.
        anotherLegend = new TLegend(0.25, 0.05, 0.55, 0.3);
        anotherLegend->SetFillStyle(0); anotherLegend->SetBorderSize(0); anotherLegend->SetTextSize(legendTextUpperCanvasSize); anotherLegend->SetTextFont(62);

        // Set the drawing style for double ratio histograms from data
        energyEnergyCorrelatorDoubleRatio[weightExponent-1][iCentrality][iJetPt]->SetMarkerStyle(kFullSquare);
        energyEnergyCorrelatorDoubleRatio[weightExponent-1][iCentrality][iJetPt]->SetMarkerColor(kBlack);
        energyEnergyCorrelatorDoubleRatio[weightExponent-1][iCentrality][iJetPt]->SetLineColor(kBlack);

        // Set drawing style for systematic uncertainties
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->SetFillColorAlpha(kBlack, 0.4);
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->SetLineColor(kBlack);
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->SetMarkerColor(kBlack);
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->SetMarkerStyle(kFullSquare);

        // Set the x-axis drawing range
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->GetYaxis()->SetRangeUser(doubleRatioZoom.first, doubleRatioZoom.second);

        // Adjust tick size for x-axis
        systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt]->GetXaxis()->SetTickLength(tickSizeX);

        // Draw first the systematic uncertainties to the upper canves
        drawer->DrawHistogramToUpperPad(systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt], "#Deltar", Form("#frac{PbPb/pp (p_{T}^{ch} > %.0f GeV)}{PbPb/pp (p_{T}^{ch} > %.0f GeV)}", trackPtCutsForDoubleDatio.second, trackPtCutsForDoubleDatio.first), " ", "e2");

        // Then draw the double ratio and add a legend for it
        energyEnergyCorrelatorDoubleRatio[weightExponent-1][iCentrality][iJetPt]->Draw("same,p");
        anotherLegend->AddEntry(systematicUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt], "Data", "lpf");

        if(theoryComparisonIndex == 0){
          // Theory comparison index 0: Hybrid model

          for(int iWake = 0; iWake < HybridModelHistogramManager::kWakeConfigurations; iWake++){

            // There are some bins for which the prediction does not exist
            if(histogrammifiedHybridModelDoubleRatio[weightExponent-1][iCentrality][iJetPt][iWake] == NULL) continue;

            // Give some nice styles for the predictions
            histogrammifiedHybridModelDoubleRatio[weightExponent-1][iCentrality][iJetPt][iWake]->SetLineColor(color[iWake]);
            histogrammifiedHybridModelDoubleRatio[weightExponent-1][iCentrality][iJetPt][iWake]->SetMarkerStyle(kFullCircle);
            histogrammifiedHybridModelDoubleRatio[weightExponent-1][iCentrality][iJetPt][iWake]->SetMarkerSize(0);
            histogrammifiedHybridModelDoubleRatio[weightExponent-1][iCentrality][iJetPt][iWake]->SetFillColorAlpha(color[iWake], jetWakeOpacity[iWake]);

            // Draw the prediction to the same canvas as the data
            histogrammifiedHybridModelDoubleRatio[weightExponent-1][iCentrality][iJetPt][iWake]->Draw("same,e3");

            // Add a legend for the theory prediction
            anotherLegend->AddEntry(histogrammifiedHybridModelDoubleRatio[weightExponent-1][iCentrality][iJetPt][iWake], hybridHistograms->GetWakeName(iWake), "f");
          }
        } else if (theoryComparisonIndex >= 5){
          // Theory comparison index >= 5: JEWEL

          for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){

            // There are some bins for which the prediction does not exist
            if(jewelDoubleRatio[weightExponent-1][iCentrality][iJetPt][iRecoil] == NULL) continue;

            // Give some nice styles for the predictions
            jewelDoubleRatio[weightExponent-1][iCentrality][iJetPt][iRecoil]->SetMarkerColor(color[iRecoil]);
            jewelDoubleRatio[weightExponent-1][iCentrality][iJetPt][iRecoil]->SetLineColor(color[iRecoil]);
            jewelDoubleRatio[weightExponent-1][iCentrality][iJetPt][iRecoil]->SetMarkerStyle(jewelMarkerStyle[iRecoil]);

            // Draw the prediction to the same canvas as the data
            jewelDoubleRatio[weightExponent-1][iCentrality][iJetPt][iRecoil]->Draw("same,pl");

            // Add a legend for the theory prediction
            anotherLegend->AddEntry(jewelDoubleRatio[weightExponent-1][iCentrality][iJetPt][iRecoil], jewelHistograms->GetRecoilName(iRecoil), "pl");
          }
        }

        // Draw the legends to the upper pad
        legend->Draw();
        anotherLegend->Draw();

        // Draw a line to one
        oneLine->Draw();

        // Draw latex messages to the plots
        mrLatexer->SetTextFont(62);
        mrLatexer->SetTextSize(0.09);
        mrLatexer->DrawLatexNDC(0.26, 0.78, "CMS");

        if(addPreliminaryTag){
          mrLatexer->SetTextFont(52);
          mrLatexer->SetTextSize(0.065);
          mrLatexer->DrawLatexNDC(0.37, 0.78, "Preliminary");
        }

        // Luminosity
        mrLatexer->SetTextFont(42);
        mrLatexer->SetTextSize(0.06);
        mrLatexer->DrawLatexNDC(0.3, 0.9, "1.70 nb^{-1} PbPb (5.02 TeV) + 302 pb^{-1} pp (5.02 TeV)");

        // Binning
        mrLatexer->SetTextFont(62);

        if(mainTextHighlightMode){
          mrLatexer->SetTextSize(0.11);
          mrLatexer->DrawLatexNDC(0.39, 0.76, jetPtString.Data());
          mrLatexer->SetTextSize(0.07);
          mrLatexer->DrawLatexNDC(0.39, 0.63, "anti-k_{T} R = 0.4, |#eta_{jet}| < 1.6");
        } else {
          mrLatexer->SetTextSize(0.07);
          mrLatexer->DrawLatexNDC(0.58, 0.78, jetPtString.Data());
          mrLatexer->DrawLatexNDC(0.712, 0.685, "anti-k_{T} R = 0.4");
          mrLatexer->DrawLatexNDC(0.795, 0.595, "|#eta_{jet}| < 1.6");
        }

        // Linear scale for the ratio
        drawer->SetLogY(false);

        // Create a new legend to show the different data uncertainty bands
        legend = new TLegend(0.25, 0.85, 0.95, 0.95);
        legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(legendTextLowerCanvasSize); legend->SetTextFont(62);
        legend->SetNColumns(2);

        // Set the axis drawing ranges
        hRelativeUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt][kRelativeUncertaintySystematic]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
        hRelativeUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt][kRelativeUncertaintySystematic]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

        // Adjust tick size for x-axis
        hRelativeUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt][kRelativeUncertaintySystematic]->GetXaxis()->SetTickLength(tickSizeX*1.5);

        // Set the style for uncertainty bands for systematic and statistical uncertainties from data
        for(int iUncertainty = 0; iUncertainty < knRelativeUncertaintyTypes; iUncertainty++){
          hRelativeUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt][iUncertainty]->SetLineColor(kBlack);
          hRelativeUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt][iUncertainty]->SetLineStyle(9);
          hRelativeUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt][iUncertainty]->SetMarkerStyle(kFullCircle);
          hRelativeUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt][iUncertainty]->SetMarkerSize(0);
        }
        hRelativeUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt][kRelativeUncertaintySystematic]->SetFillColorAlpha(kBlack, 0.4);
        hRelativeUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt][kRelativeUncertaintySystematic]->SetLineWidth(0);


        // Draw the error bars from the data
        drawer->DrawHistogramToLowerPad(hRelativeUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt][kRelativeUncertaintySystematic], "#Deltar", "#frac{Theory}{Data}", " ", "e3");
        hRelativeUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt][kRelativeUncertaintyStatisticalUp]->Draw("same,HIST,C");
        hRelativeUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt][kRelativeUncertaintyStatisticalDown]->Draw("same,HIST,C");

        // Add the different uncertainties to the legend
        legend->AddEntry(hRelativeUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt][kRelativeUncertaintySystematic], "Data syst. unc.", "f");
        legend->AddEntry(hRelativeUncertaintyDoubleRatio[weightExponent-1][iCentrality][iJetPt][kRelativeUncertaintyStatisticalUp], "Data stat. unc.", "l");

        if(theoryComparisonIndex == 0){
          // Theory comparison index 0: Hybrid model

          for(int iWake = 0; iWake < HybridModelHistogramManager::kWakeConfigurations; iWake++){
            hybridModelToDataRatioDoubleRatio[weightExponent-1][iCentrality][iJetPt][iWake]->SetMarkerStyle(kFullCircle);
            hybridModelToDataRatioDoubleRatio[weightExponent-1][iCentrality][iJetPt][iWake]->SetMarkerSize(0);
            hybridModelToDataRatioDoubleRatio[weightExponent-1][iCentrality][iJetPt][iWake]->SetFillColorAlpha(color[iWake], jetWakeOpacity[iWake]);
            hybridModelToDataRatioDoubleRatio[weightExponent-1][iCentrality][iJetPt][iWake]->Draw("same,e3");
          }
        } else if(theoryComparisonIndex >= 5){
          // Theory comparison index >= 5: JEWEL

          for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){

            // There are some bins for which the prediction does not exist
            if(jewelToDataRatioDoubleRatio[weightExponent-1][iCentrality][iJetPt][iRecoil] == NULL) continue;

            jewelToDataRatioDoubleRatio[weightExponent-1][iCentrality][iJetPt][iRecoil]->SetMarkerStyle(jewelMarkerStyle[iRecoil]);
            jewelToDataRatioDoubleRatio[weightExponent-1][iCentrality][iJetPt][iRecoil]->SetMarkerColor(color[iRecoil]);
            jewelToDataRatioDoubleRatio[weightExponent-1][iCentrality][iJetPt][iRecoil]->SetLineColor(color[iRecoil]);
            jewelToDataRatioDoubleRatio[weightExponent-1][iCentrality][iJetPt][iRecoil]->Draw("same,pl");
          }
        }

        // Draw the legend to the lower pad
        legend->Draw();

        // Draw a line to one
        oneLine->Draw();

        // If a plot name is given, save the plot in a file
        if(saveFigures) {
          gPad->GetCanvas()->SaveAs(Form("figures/energyEnergyCorrelator_%sDoubleRatio%s%s%s.%s", theorySaveName[theoryComparisonIndex].Data(), saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), figureFormat.Data()));
        }

      } // Jet pT loop
    } // Centrality loop
  } // Comparing theory with the data for double ratios
}
