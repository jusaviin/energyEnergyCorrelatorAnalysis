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
 *       7 = Full hybrid and JEWEL with recoils
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
  // 7 = Full hybrid and JEWEL with recoils
  TString theorySaveName[8] = {"hybridModel", "holguin", "colbt", "threeModels", "holguinAndColbt", "jewel", "holguinAndJewel", "hybridAndJewel"};


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

      // Create a new systematic uncertainty organizer
      uncertainties[iDataType][iWeightExponent] = new SystematicUncertaintyOrganizer(uncertaintyFile[iDataType][iWeightExponent]);

    } // Data type loop

    for(int iMCType = 0; iMCType < knPpMCTypes; iMCType++){

      // Create a new histogram manager
      ppMCHistograms[iMCType][iWeightExponent] = new EECHistogramManager(ppMCFile[iMCType][iWeightExponent], ppMCCard[iMCType][iWeightExponent]);

    }

  } // Weight exponent loop
 
  // Energy-energy correlators and PbPb to pp ratios
  std::map<std::tuple<int,int,int,int>, TH1D*> energyEnergyCorrelatorRawPbPb; // Dimensions: weight exponent, centrality, jet pT, track pT
  std::map<std::tuple<int,int,int>, TH1D*> energyEnergyCorrelatorRawPp; // Dimensions: weight exponent, jet pT, track pT
  std::map<std::tuple<int,int,int,int>, TH1D*> energyEnergyCorrelatorSignalPbPb; // Dimensions: weight exponent, centrality, jet pT, track pT
  std::map<std::tuple<int,int,int>, TH1D*> energyEnergyCorrelatorSignalPp; // Dimensions: weight exponent, jet pT, track pT
  std::map<std::tuple<int,int,int,int>, TH1D*> energyEnergyCorrelatorPbPbToPpRatio; // Dimensions: weight exponent, centrality, jet pT, track pT
  std::map<std::tuple<int,int,int,int,int>, TH1D*> systematicUncertaintyForPbPb; // Dimensions: weight exponent, systematic type, centrality, jet pT, track pT
  std::map<std::tuple<int,int,int,int>, TH1D*>  systematicUncertaintyForPp; // Dimensions: weight exponent, systematic type, jet pT, track pT
  std::map<std::tuple<int,int,int,int,int>, TH1D*>  systematicUncertaintyPbPbToPpRatio; // Dimensions: weight exponent, systematic type, centrality, jet pT, track pT

  // Histograms for relative uncertainties
  std::map<std::tuple<int,int,int,int>, TH1D*> hRelativeUncertaintyPp; // Dimensions: weight exponent, jet pT, track pT, relative uncertainty type
  std::map<std::tuple<int,int,int,int,int>, TH1D*> hRelativeUncertaintyPbPb; // Dimensions: weight exponent, centrality jet pT, track pT, relative uncertainty type
  std::map<std::tuple<int,int,int,int,int>, TH1D*> hRelativeUncertaintyPbPbToPpRatio; // Dimensions: weight exponent, centrality jet pT, track pT, relative uncertainty type

  // Double ratios from energy-energy correlator histograms
  std::map<std::tuple<int,int,int,int>, TH1D*> energyEnergyCorrelatorForDoubleRatioFromPbPb; // Dimensions: weight exponent, centrality, jet pT, track pT
  std::map<std::tuple<int,int,int>, TH1D*> energyEnergyCorrelatorForDoubleRatioFromPp; // Dimensions: weight exponent, jet pT, track pT
  std::map<std::tuple<int,int,int>, TH1D*> energyEnergyCorrelatorDoubleRatio; // Dimensions: weight exponent, centrality, jet pT
  std::map<std::tuple<int,int,int,int>, TH1D*>  systematicUncertaintyForDoubleRatioFromPbPb; // Dimensions: weight exponent, centrality, jet pT, track pT
  std::map<std::tuple<int,int,int>, TH1D*> systematicUncertaintyForDoubleRatioFromPp; // Dimensions: weight exponent, jet pT, track pT
  std::map<std::tuple<int,int,int>, TH1D*> systematicUncertaintyDoubleRatio; // Dimensions: weight exponent, centrality, jet pT

  // Relative uncertainties for double ratios
  std::map<std::tuple<int,int,int,int>, TH1D*> hRelativeUncertaintyDoubleRatio; // Dimensions: weight exponent, centrality, jet pT, relative uncertainty type

  // Energy-energy correlators from pp MC models
  std::map<std::tuple<int,int,int,int>, TH1D*> hEnergyEnergyCorrelatorPpMC; // Dimensions: weight exponent, jet pT, track pT, MC type index
  std::map<std::tuple<int,int,int,int>, TH1D*> hEnergyEnergyCorrelatorPpMCToDataRatio; // Dimensions: weight exponent, jet pT, track pT, MC type index

  // Energy-energy correlators for Hybrid model
  std::map<std::tuple<int, int, int, int>, TGraphErrors*> energyEnergyCorrelatorHybridModelPp; // Dimensions: weight exponent, jet pT, track pT, wake configuration
  std::map<std::tuple<int, int, int, int, int>, TGraphErrors*> energyEnergyCorrelatorHybridModelPbPb; // Dimensions: weight exponent, centrality, jet pT, track pT, wake configuration
  std::map<std::tuple<int, int, int, int, int>, TGraphErrors*> energyEnergyCorrelatorHybridModelPbPbToPpRatio; // Dimensions: weight exponent, centrality, jet pT, track pT, wake configuration
  std::map<std::tuple<int, int, int, int>, TH1D*> histogrammifiedHybridModelPp; // Dimensions: weight exponent, jet pT, track pT, wake configuration
  std::map<std::tuple<int, int, int, int, int>, TH1D*> histogrammifiedHybridModelPbPb; // Dimensions: weight exponent, centrality, jet pT, track pT, wake configuration
  std::map<std::tuple<int, int, int, int, int>, TH1D*> histogrammifiedHybridModelPbPbToPpRatio; // Dimensions: weight exponent, centrality, jet pT, track pT, wake configuration
  std::map<std::tuple<int, int, int, int>, TH1D*> hybridModelToDataRatioPp; // Dimensions: weight exponent, jet pT, track pT, wake configuration
  std::map<std::tuple<int, int, int, int, int>, TH1D*> hybridModelToDataRatioPbPb; // Dimensions: weight exponent, centrality, jet pT, track pT, wake configuration
  std::map<std::tuple<int, int, int, int, int>, TH1D*> hybridModelToDataRatioPbPbToPpRatio; // Dimensions: weight exponent, centrality, jet pT, track pT, wake configuration

   // Double ratio histograms for the Hybrid model
  std::map<std::tuple<int, int, int, int>, TH1D*> histogrammifiedHybridModelDoubleRatio; // Dimensions: weight exponent, centrality, jet pT, wake configuration
  std::map<std::tuple<int, int, int, int>, TH1D*> hybridModelToDataRatioDoubleRatio; // Dimensions: weight exponent, centrality, jet pT, wake configuration

  // Energy-energy correlators for perturbative calculations by Holguin et. al.
  std::map<std::tuple<int, int, int, int, int>, TGraph*>  energyEnergyCorrelatorHolguinPbPb; // Dimensions: weight exponent, centrality, jet pT, track pT, k-value
  std::map<std::tuple<int, int, int, int, int>, TGraph*> energyEnergyCorrelatorHolguinPbPbToPpRatio; // Dimensions: weight exponent, centrality, jet pT, track pT, k-value
  std::map<std::tuple<int, int, int, int, int>, TH1D*> histogrammifiedHolguinPbPb; // Dimensions: weight exponent, centrality, jet pT, track pT, k-value
  std::map<std::tuple<int, int, int, int, int>, TH1D*> histogrammifiedHolguinPbPbToPpRatio; // Dimensions: weight exponent, centrality, jet pT, track pT, k-value
  std::map<std::tuple<int, int, int, int, int>, TH1D*> holguinToDataRatioPbPb; // Dimensions: weight exponent, centrality, jet pT, track pT, k-value
  std::map<std::tuple<int, int, int, int, int>, TH1D*> holguinToDataRatioPbPbToPpRatio; // Dimensions: weight exponent, centrality, jet pT, track pT, k-value

  // Energy-energy correlators for CoLBT predictions
  std::map<std::tuple<int, int, int, int, int>, TGraphErrors*> energyEnergyCorrelatorCoLBTPbPbToPpRatio; // Dimensions: weight exponent, centrality, jet pT, track pT, q-value
  std::map<std::tuple<int, int, int, int, int>, TH1D*> histogrammifiedCoLBTPbPbToPpRatio; // Dimensions: weight exponent, centrality, jet pT, track pT, q-value
  std::map<std::tuple<int, int, int, int, int>, TH1D*> coLBTToDataRatioPbPbToPpRatio; // Dimensions: weight exponent, centrality, jet pT, track pT, q-value

  // Energy-energy correlators for JEWEL predictions
  std::map<std::tuple<int, int, int>, TH1D*> energyEnergyCorrelatorJewelPp; // Dimentions: weight exponent, jet pT, track pT
  std::map<std::tuple<int, int, int, int, int>, TH1D*> energyEnergyCorrelatorJewelPbPb; // Dimensions: weight exponent, centrality, jet pT, track pT, recoil setting
  std::map<std::tuple<int, int, int, int, int>, TH1D*> energyEnergyCorrelatorJewelPbPbToPpRatio; // Dimensions: weight exponent, centrality, jet pT, track pT, recoil setting

  std::map<std::tuple<int, int, int>, TH1D*> jewelToDataRatioPp; // Dimentions: weight exponent, jet pT, track pT
  std::map<std::tuple<int, int, int, int, int>, TH1D*> jewelToDataRatioPbPb; // Dimensions: weight exponent, centrality, jet pT, track pT, recoil setting
  std::map<std::tuple<int, int, int, int, int>, TH1D*> jewelToDataRatioPbPbToPpRatio; // Dimensions: weight exponent, centrality, jet pT, track pT, recoil setting
  std::map<std::tuple<int, int, int, int>, TH1D*> jewelDoubleRatio; // Dimensions: weight exponent, centrality, jet pT, recoil setting
  std::map<std::tuple<int, int, int, int>, TH1D*> jewelToDataRatioDoubleRatio; // Dimensions: weight exponent, centrality, jet pT, recoil setting
 
  // Helper variables for binning
  std::tuple<int, int, int, int> eecPbPbBin;
  std::tuple<int, int, int, int, int> eecPbPbCorrelatedSystematicsBin;
  std::tuple<int, int, int, int, int> eecPbPbCorrelatedSystematicsBinShapeUp;
  std::tuple<int, int, int, int, int> eecPbPbCorrelatedSystematicsBinShapeDown;
  std::tuple<int, int, int, int, int> eecPbPbUncorrelatedSystematicsBin;
  std::tuple<int, int, int, int, int> eecPbPbRelativeUncertaintyBin;
  std::tuple<int, int, int, int, int> eecPbPbRelativeUncertaintyBinSystematic;
  std::tuple<int, int, int, int, int> eecPbPbRelativeUncertaintyBinStatisticalUp;
  std::tuple<int, int, int, int, int> eecPbPbRelativeUncertaintyBinStatisticalDown;

  std::tuple<int, int, int> eecPpBin;
  std::tuple<int, int, int, int> eecPpCorrelatedSystematicsBin;
  std::tuple<int, int, int, int> eecPpCorrelatedSystematicsBinShapeUp;
  std::tuple<int, int, int, int> eecPpCorrelatedSystematicsBinShapeDown;
  std::tuple<int, int, int, int> eecPpUncorrelatedSystematicsBin;
  std::tuple<int, int, int, int> eecPpRelativeUncertaintyBin;
  std::tuple<int, int, int, int> eecPpRelativeUncertaintyBinSystematic;
  std::tuple<int, int, int, int> eecPpRelativeUncertaintyBinStatisticalUp;
  std::tuple<int, int, int, int> eecPpRelativeUncertaintyBinStatisticalDown;

  std::tuple<int, int, int> eecDoubleRatioBin;
  std::tuple<int, int, int, int> eecDoubleRatioRelativeUncertaintyBin;
  std::tuple<int, int, int, int> eecDoubleRatioRelativeUncertaintyBinSystematic;
  std::tuple<int, int, int, int> eecDoubleRatioRelativeUncertaintyBinStatisticalUp;
  std::tuple<int, int, int, int> eecDoubleRatioRelativeUncertaintyBinStatisticalDown;

  std::tuple<int, int, int, int> eecPpMCBin;
  std::tuple<int, int, int, int, int> eecPbPbMCBin;
  std::tuple<int, int, int, int> eecDoubleRatioMCBin;

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
        eecPpBin = std::make_tuple(iWeightExponent,iJetPt,iTrackPt);
        energyEnergyCorrelatorRawPp[eecPpBin] = histograms[kPp][iWeightExponent]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatchedPp, iTrackPtMatchedPp);

        for(auto centralityBin : drawnCentralityBin){

          iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);

          // Read the raw PbPb energy-energy correlator distributions
          eecPbPbBin = std::make_tuple(iWeightExponent,iCentrality,iJetPt,iTrackPt);
          energyEnergyCorrelatorRawPbPb[eecPbPbBin] = histograms[kPbPb][iWeightExponent]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt);
        }
      }

      // Then go to the main track pT loop
      for(auto trackPtBin : drawnTrackPtBin){
        iTrackPt = card[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);
        iTrackPtMatchedPp = card[kPp][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);
        iTrackPtMatchedPbPbUncertainty = uncertaintyCard[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);
        iTrackPtMatchedPpUncertainty = uncertaintyCard[kPp][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);
        eecPpBin = std::make_tuple(iWeightExponent, iJetPt, iTrackPt);
        eecPpUncorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kUncorrelatedUncertainty, iJetPt, iTrackPt);
        eecPpCorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kCorrelatedUncertainty, iJetPt, iTrackPt);
        eecPpRelativeUncertaintyBinSystematic = std::make_tuple(iWeightExponent, iJetPt, iTrackPt, kRelativeUncertaintySystematic);
        eecPpRelativeUncertaintyBinStatisticalUp = std::make_tuple(iWeightExponent, iJetPt, iTrackPt, kRelativeUncertaintyStatisticalUp);
        eecPpRelativeUncertaintyBinStatisticalDown = std::make_tuple(iWeightExponent, iJetPt, iTrackPt, kRelativeUncertaintyStatisticalDown);

        // Read the pp histograms that do not have centrality binning
        energyEnergyCorrelatorSignalPp[eecPpBin] = histograms[kPp][iWeightExponent]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatchedPp, iTrackPtMatchedPp, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);
        systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin] = uncertainties[kPp][iWeightExponent]->GetUncorrelatedSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty);
        systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin] = uncertainties[kPp][iWeightExponent]->GetCorrelatedSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty);

        // For relative uncertainties, read the total systematic uncertainties for pp collisions and clone the signal distribution
        hRelativeUncertaintyPp[eecPpRelativeUncertaintyBinSystematic] = uncertainties[kPp][iWeightExponent]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty);
        hRelativeUncertaintyPp[eecPpRelativeUncertaintyBinStatisticalUp] = (TH1D*) energyEnergyCorrelatorSignalPp[eecPpBin]->Clone(Form("relativeStatisticalUncertaintyUpPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));
        hRelativeUncertaintyPp[eecPpRelativeUncertaintyBinStatisticalDown] = (TH1D*) energyEnergyCorrelatorSignalPp[eecPpBin]->Clone(Form("relativeStatisticalUncertaintyDownPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));

        // Transform the uncertainties into relative uncertainties
        for(int iUncertainty = 0; iUncertainty < knRelativeUncertaintyTypes; iUncertainty++){
          eecPpRelativeUncertaintyBin = std::make_tuple(iWeightExponent, iJetPt, iTrackPt, iUncertainty);
          optimusPrimeTheTransformer->TransformToRelativeUncertainty(hRelativeUncertaintyPp[eecPpRelativeUncertaintyBin], true);
        }

        // Move the statistical uncertainty bin contents up and down by the uncertainties
        for(int iBin = 1; iBin <= hRelativeUncertaintyPp[eecPpRelativeUncertaintyBinStatisticalUp]->GetNbinsX(); iBin++){
          binError = hRelativeUncertaintyPp[eecPpRelativeUncertaintyBinStatisticalUp]->GetBinError(iBin);
          hRelativeUncertaintyPp[eecPpRelativeUncertaintyBinStatisticalUp]->SetBinContent(iBin, 1 + binError);
          hRelativeUncertaintyPp[eecPpRelativeUncertaintyBinStatisticalDown]->SetBinContent(iBin, 1 - binError);
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
          eecPpBin = std::make_tuple(iWeightExponent,iJetPt,iTrackPt);
          lowAnalysisBin = energyEnergyCorrelatorRawPp[eecPpBin]->GetXaxis()->FindBin(analysisDeltaR.first + epsilon);
          highAnalysisBin = energyEnergyCorrelatorRawPp[eecPpBin]->GetXaxis()->FindBin(analysisDeltaR.second - epsilon);
          lowPtIntegral = energyEnergyCorrelatorRawPp[std::make_tuple(iWeightExponent,iJetPt,trackPtBinsForDoubleRatio.first)]->Integral(lowAnalysisBin, highAnalysisBin, "width");
          highPtIntegral = energyEnergyCorrelatorRawPp[std::make_tuple(iWeightExponent,iJetPt,trackPtBinsForDoubleRatio.second)]->Integral(lowAnalysisBin, highAnalysisBin, "width");
          trackCorrelation = 1 - highPtIntegral / lowPtIntegral;

          // Now for the double ratio, we can combine the relevant uncertainties while scaling down the tracking related ones
          // by the expected overlap.
          systematicUncertaintyForDoubleRatioFromPp[eecPpBin] = (TH1D*) uncertaintyBackgroundSubtraction->Clone(Form("doubleRatioUncertaintyFromPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));
          for(int iBin = 1; iBin <= uncertaintyBackgroundSubtraction->GetNbinsX(); iBin++){
            systematicUncertaintyForDoubleRatioFromPp[eecPpBin]->SetBinError(iBin, TMath::Sqrt(TMath::Power(uncertaintyBackgroundSubtraction->GetBinError(iBin),2) + TMath::Power(uncertaintyTrackPairEfficiency->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintyTrackSelection->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintyMCnonclosure->GetBinError(iBin),2)));
          }

        } // Histograms relevant to double ratio


        for(auto centralityBin : drawnCentralityBin){

          iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);
          iCentralityMatched = uncertaintyCard[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);

          // Read the PbPb histograms
          eecPbPbBin = std::make_tuple(iWeightExponent,iCentrality,iJetPt,iTrackPt);
          eecPbPbUncorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kCorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
          eecPbPbRelativeUncertaintyBinSystematic = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt, kRelativeUncertaintySystematic);
          eecPbPbRelativeUncertaintyBinStatisticalUp = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt, kRelativeUncertaintyStatisticalUp);
          eecPbPbRelativeUncertaintyBinStatisticalDown = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt, kRelativeUncertaintyStatisticalDown);

          energyEnergyCorrelatorSignalPbPb[eecPbPbBin] = histograms[kPbPb][iWeightExponent]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin] = uncertainties[kPbPb][iWeightExponent]->GetUncorrelatedSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBin] = uncertainties[kPbPb][iWeightExponent]->GetCorrelatedSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty);

          // For relative uncertainties, read the total systematic uncertainties for PbPb collisions and clone the signal distribution
          hRelativeUncertaintyPbPb[eecPbPbRelativeUncertaintyBinSystematic] = uncertainties[kPbPb][iWeightExponent]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty);
          hRelativeUncertaintyPbPb[eecPbPbRelativeUncertaintyBinStatisticalUp] = (TH1D*) energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->Clone(Form("relativeStatisticalUncertaintyUpPbPb%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          hRelativeUncertaintyPbPb[eecPbPbRelativeUncertaintyBinStatisticalDown] = (TH1D*) energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->Clone(Form("relativeStatisticalUncertaintyDownPbPb%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));

          // Transform the uncertainties into relative uncertainties
          for(int iUncertainty = 0; iUncertainty < knRelativeUncertaintyTypes; iUncertainty++){
            eecPbPbRelativeUncertaintyBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt, iUncertainty);
            optimusPrimeTheTransformer->TransformToRelativeUncertainty(hRelativeUncertaintyPbPb[eecPbPbRelativeUncertaintyBin], true);
          }

          // Move the statistical uncertainty bin contents up and down by the uncertainties
          for(int iBin = 1; iBin <= hRelativeUncertaintyPbPb[eecPbPbRelativeUncertaintyBinStatisticalUp]->GetNbinsX(); iBin++){
            binError = hRelativeUncertaintyPbPb[eecPbPbRelativeUncertaintyBinStatisticalUp]->GetBinError(iBin);
            hRelativeUncertaintyPbPb[eecPbPbRelativeUncertaintyBinStatisticalUp]->SetBinContent(iBin, 1 + binError);
            hRelativeUncertaintyPbPb[eecPbPbRelativeUncertaintyBinStatisticalDown]->SetBinContent(iBin, 1 - binError);
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
            eecPbPbBin = std::make_tuple(iWeightExponent,iCentrality,iJetPt,iTrackPt);
            lowAnalysisBin = energyEnergyCorrelatorRawPbPb[eecPbPbBin]->GetXaxis()->FindBin(analysisDeltaR.first + epsilon);
            highAnalysisBin = energyEnergyCorrelatorRawPbPb[eecPbPbBin]->GetXaxis()->FindBin(analysisDeltaR.second - epsilon);
            lowPtIntegral = energyEnergyCorrelatorRawPbPb[std::make_tuple(iWeightExponent,iCentrality,iJetPt,trackPtBinsForDoubleRatio.first)]->Integral(lowAnalysisBin, highAnalysisBin, "width");
            highPtIntegral = energyEnergyCorrelatorRawPbPb[std::make_tuple(iWeightExponent,iCentrality,iJetPt,trackPtBinsForDoubleRatio.second)]->Integral(lowAnalysisBin, highAnalysisBin, "width");
            trackCorrelation = 1 - highPtIntegral / lowPtIntegral;

            // Now for the double ratio, we can combine the relevant uncertainties while scaling down the tracking related ones
            // by the expected overlap.
            systematicUncertaintyForDoubleRatioFromPbPb[eecPbPbBin] = (TH1D*) uncertaintyBackgroundSubtraction->Clone(Form("doubleRatioUncertaintyFromPbPb%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
            for(int iBin = 1; iBin <= uncertaintyBackgroundSubtraction->GetNbinsX(); iBin++){
              systematicUncertaintyForDoubleRatioFromPbPb[eecPbPbBin]->SetBinError(iBin, TMath::Sqrt(TMath::Power(uncertaintyBackgroundSubtraction->GetBinError(iBin),2) + TMath::Power(uncertaintyTrackPairEfficiency->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintyTrackSelection->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintySignalToBackgroundRatio->GetBinError(iBin)*trackCorrelation,2) + TMath::Power(uncertaintyMCnonclosure->GetBinError(iBin),2)));
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
    eecPbPbBin = std::make_tuple(iWeightExponent,iCentrality,iJetPt,iTrackPt);
    lowAnalysisBin = energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->GetXaxis()->FindBin(analysisDeltaR.first + epsilon);
    highAnalysisBin = energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->GetXaxis()->FindBin(analysisDeltaR.second - epsilon);
    for(auto jetPtBin : drawnJetPtBin){
      iJetPt = card[kPbPb][iWeightExponent]->FindBinIndexJetPtEEC(jetPtBin);
      for(auto trackPtBin : drawnTrackPtBin){
        iTrackPt = card[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);
        eecPpBin = std::make_tuple(iWeightExponent, iJetPt, iTrackPt);
        eecPpUncorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kUncorrelatedUncertainty, iJetPt, iTrackPt);
        eecPpCorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kCorrelatedUncertainty, iJetPt, iTrackPt);

        energyEnergyCorrelatorSignalPp[eecPpBin]->Scale(1.0 / energyEnergyCorrelatorSignalPp[eecPpBin]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
        systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->Scale(energyEnergyCorrelatorSignalPp[eecPpBin]->GetBinContent(10) / systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetBinContent(10));
        systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->Scale(energyEnergyCorrelatorSignalPp[eecPpBin]->GetBinContent(10) / systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->GetBinContent(10));

        if(drawDoubleRatioDataToTheoryComparison){
          systematicUncertaintyForDoubleRatioFromPp[eecPpBin]->Scale(energyEnergyCorrelatorSignalPp[eecPpBin]->GetBinContent(10) / systematicUncertaintyForDoubleRatioFromPp[eecPpBin]->GetBinContent(10));
        } // Stuff relevant for double ratios

        for(auto centralityBin : drawnCentralityBin){

          iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);
          eecPbPbBin = std::make_tuple(iWeightExponent,iCentrality,iJetPt,iTrackPt);
          eecPbPbUncorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kCorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);

          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->Scale(1.0 / energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->Scale(energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->GetBinContent(10) / systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetBinContent(10));
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBin]->Scale(energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->GetBinContent(10) / systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBin]->GetBinContent(10));
          
          if(drawDoubleRatioDataToTheoryComparison){
            systematicUncertaintyForDoubleRatioFromPbPb[eecPbPbBin]->Scale(energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->GetBinContent(10) / systematicUncertaintyForDoubleRatioFromPbPb[eecPbPbBin]->GetBinContent(10));
          } // Studd relevant for double ratios

          energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin] = (TH1D*) energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->Clone(Form("energyEnergyCorrelatorRatio%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->Divide(energyEnergyCorrelatorSignalPp[eecPpBin]);

          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin] = (TH1D*) systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->Clone(Form("systematicUncertaintyRatioUncorrelated%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->Divide(systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]);

          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBin] = (TH1D*) systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBin]->Clone(Form("systematicUncertaintyRatioCorrelated%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBin]->Divide(systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]);
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
          eecPpBin = std::make_tuple(iWeightExponent, iJetPt, iTrackPt);

          energyEnergyCorrelatorForDoubleRatioFromPp[eecPpBin] = (TH1D*) energyEnergyCorrelatorSignalPp[eecPpBin]->Clone(Form("energyEnergyCorrelatorForDoubleRatioFromPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));

          for(int iBin = 1; iBin <= energyEnergyCorrelatorForDoubleRatioFromPp[eecPpBin]->GetNbinsX(); iBin++){
            energyEnergyCorrelatorForDoubleRatioFromPp[eecPpBin]->SetBinError(iBin, energyEnergyCorrelatorForDoubleRatioFromPp[eecPpBin]->GetBinError(iBin) * trackCorrelation);
          }

          for(auto centralityBin : drawnCentralityBin){
            iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);

            eecPbPbBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt);
            energyEnergyCorrelatorForDoubleRatioFromPbPb[eecPbPbBin] = (TH1D*) energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->Clone(Form("energyEnergyCorrelatorForDoubleRatioFromPbPb%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));

            for(int iBin = 1; iBin <= energyEnergyCorrelatorForDoubleRatioFromPbPb[eecPbPbBin]->GetNbinsX(); iBin++){
              energyEnergyCorrelatorForDoubleRatioFromPbPb[eecPbPbBin]->SetBinError(iBin, energyEnergyCorrelatorForDoubleRatioFromPbPb[eecPbPbBin]->GetBinError(iBin) * trackCorrelation);
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
            eecPpBin = std::make_tuple(iWeightExponent, iJetPt, iTrackPt);
            eecPbPbBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt);

            // Calculate the single ratios with properly handled double ratio uncertainties
            energyEnergyCorrelatorForDoubleRatioFromPbPb[eecPbPbBin]->Divide(energyEnergyCorrelatorForDoubleRatioFromPp[eecPpBin]);
            systematicUncertaintyForDoubleRatioFromPbPb[eecPbPbBin]->Divide(systematicUncertaintyForDoubleRatioFromPp[eecPpBin]);
          } // Track pT loop

          eecDoubleRatioBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt);

          // Calculate the double ratios from the single ratios with properly handled uncertainties
          energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin] = (TH1D*) energyEnergyCorrelatorForDoubleRatioFromPbPb[std::make_tuple(iWeightExponent,iCentrality,iJetPt,trackPtBinsForDoubleRatio.second)]->Clone(Form("energyEnergyCorrelatorDoubleRatio%d%d%d", iWeightExponent, iCentrality, iJetPt));
          energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin]->Divide(energyEnergyCorrelatorForDoubleRatioFromPbPb[std::make_tuple(iWeightExponent,iCentrality,iJetPt,trackPtBinsForDoubleRatio.first)]);
          systematicUncertaintyDoubleRatio[eecDoubleRatioBin] = (TH1D*) systematicUncertaintyForDoubleRatioFromPbPb[std::make_tuple(iWeightExponent,iCentrality,iJetPt,trackPtBinsForDoubleRatio.second)]->Clone(Form("systematicUncertaintyDoubleRatio%d%d%d", iWeightExponent, iCentrality, iJetPt));
          systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->Divide(systematicUncertaintyForDoubleRatioFromPbPb[std::make_tuple(iWeightExponent,iCentrality,iJetPt,trackPtBinsForDoubleRatio.first)]);
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
        eecPpCorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kCorrelatedUncertainty, iJetPt, iTrackPt);
        eecPpCorrelatedSystematicsBinShapeUp = std::make_tuple(iWeightExponent, kCorrelatedUncertaintyShapeUp, iJetPt, iTrackPt);
        eecPpCorrelatedSystematicsBinShapeDown = std::make_tuple(iWeightExponent, kCorrelatedUncertaintyShapeDown, iJetPt, iTrackPt);

        systematicUncertaintyForPp[eecPpCorrelatedSystematicsBinShapeUp] = (TH1D*) systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->Clone(Form("correlatedUpBandPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));
        systematicUncertaintyForPp[eecPpCorrelatedSystematicsBinShapeDown] = (TH1D*) systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->Clone(Form("correlatedDownBandPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));
        calculateCorrelatedBands(systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin], systematicUncertaintyForPp[eecPpCorrelatedSystematicsBinShapeUp], systematicUncertaintyForPp[eecPpCorrelatedSystematicsBinShapeDown]);

        for(auto centralityBin : drawnCentralityBin){

          iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);
          eecPbPbCorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kCorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(iWeightExponent, kCorrelatedUncertaintyShapeUp, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(iWeightExponent, kCorrelatedUncertaintyShapeDown, iCentrality, iJetPt, iTrackPt);

          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp] = (TH1D*) systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBin]->Clone(Form("correlatedUpBand%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown] = (TH1D*) systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBin]->Clone(Form("correlatedDownBand%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          calculateCorrelatedBands(systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBin], systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp], systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]);

          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp] = (TH1D*) systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBin]->Clone(Form("correlatedUpBandForRatio%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown] = (TH1D*) systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBin]->Clone(Form("correlatedDownBandForRatio%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          calculateCorrelatedBands(systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBin], systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp], systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]);

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
          eecPbPbBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kCorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
          eecPbPbUncorrelatedSystematicsBin = std::make_tuple(iWeightExponent, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
          eecPbPbRelativeUncertaintyBinSystematic = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt, kRelativeUncertaintySystematic);
          eecPbPbRelativeUncertaintyBinStatisticalUp = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt, kRelativeUncertaintyStatisticalUp);
          eecPbPbRelativeUncertaintyBinStatisticalDown = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt, kRelativeUncertaintyStatisticalDown);

          // For relative systematic uncertainty, add the correlated and uncorrelated systematic uncertainties in quadrature
          hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBinSystematic] = (TH1D*) systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBin]->Clone(Form("relativeSystematicUncertaintyForPbPbToPpRatio%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));

          for(int iBin = 1; iBin <= hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBinSystematic]->GetNbinsX(); iBin++){
            binError = hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBinSystematic]->GetBinError(iBin);
            otherBinError = systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetBinError(iBin);
            hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBinSystematic]->SetBinError(iBin, TMath::Sqrt(binError*binError + otherBinError*otherBinError));
          }

          // Relative statistical uncertainties we can directly get from the ratio
          hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBinStatisticalUp] = (TH1D*) energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->Clone(Form("relativeStatisticalUncertaintyUpForPbPbToPpRatio%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));
          hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBinStatisticalDown] = (TH1D*) energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->Clone(Form("relativeStatisticalUncertaintyDownForPbPbToPpRatio%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt));

          // Transform the uncertainties into relative uncertainties
          for(int iUncertainty = 0; iUncertainty < knRelativeUncertaintyTypes; iUncertainty++){
            eecPbPbRelativeUncertaintyBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt, iUncertainty);
            optimusPrimeTheTransformer->TransformToRelativeUncertainty(hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBin], true);
          }

          // Move the statistical uncertainty bin contents up and down by the uncertainties
          for(int iBin = 1; iBin <= hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBinStatisticalUp]->GetNbinsX(); iBin++){
            binError = hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBinStatisticalUp]->GetBinError(iBin);
            hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBinStatisticalUp]->SetBinContent(iBin, 1 + binError);
            hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBinStatisticalDown]->SetBinContent(iBin, 1 - binError);
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
          eecDoubleRatioBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt);
          eecDoubleRatioRelativeUncertaintyBinSystematic = std::make_tuple(iWeightExponent, iCentrality, iJetPt, kRelativeUncertaintySystematic);
          eecDoubleRatioRelativeUncertaintyBinStatisticalUp = std::make_tuple(iWeightExponent, iCentrality, iJetPt, kRelativeUncertaintyStatisticalUp);
          eecDoubleRatioRelativeUncertaintyBinStatisticalDown = std::make_tuple(iWeightExponent, iCentrality, iJetPt, kRelativeUncertaintyStatisticalDown);

          // For relative uncertainties, read the total systematic uncertainties
          hRelativeUncertaintyDoubleRatio[eecDoubleRatioRelativeUncertaintyBinSystematic] = (TH1D*) systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->Clone(Form("relativeSystematicUncertaintyDoubleRatio%d%d%d", iWeightExponent, iCentrality, iJetPt));
          hRelativeUncertaintyDoubleRatio[eecDoubleRatioRelativeUncertaintyBinStatisticalUp] = (TH1D*) energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin]->Clone(Form("relativeStatisticalUncertaintyUpDoubleRatio%d%d%d", iWeightExponent, iCentrality, iJetPt));
          hRelativeUncertaintyDoubleRatio[eecDoubleRatioRelativeUncertaintyBinStatisticalDown] = (TH1D*) energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin]->Clone(Form("relativeStatisticalUncertaintyDownDoubleRatio%d%d%d", iWeightExponent, iCentrality, iJetPt));

          // Transform the uncertainties into relative uncertainties
          for(int iUncertainty = 0; iUncertainty < knRelativeUncertaintyTypes; iUncertainty++){
            eecDoubleRatioRelativeUncertaintyBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iUncertainty);
            optimusPrimeTheTransformer->TransformToRelativeUncertainty(hRelativeUncertaintyDoubleRatio[eecDoubleRatioRelativeUncertaintyBin], true);
          }

          // Move the statistical uncertainty bin contents up and down by the uncertainties
          for(int iBin = 1; iBin <= hRelativeUncertaintyDoubleRatio[eecDoubleRatioRelativeUncertaintyBinStatisticalUp]->GetNbinsX(); iBin++){
            binError = hRelativeUncertaintyDoubleRatio[eecDoubleRatioRelativeUncertaintyBinStatisticalUp]->GetBinError(iBin);
            hRelativeUncertaintyDoubleRatio[eecDoubleRatioRelativeUncertaintyBinStatisticalUp]->SetBinContent(iBin, 1 + binError);
            hRelativeUncertaintyDoubleRatio[eecDoubleRatioRelativeUncertaintyBinStatisticalDown]->SetBinContent(iBin, 1 - binError);
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
          eecPpBin = std::make_tuple(iWeightExponent, iJetPt, iTrackPt);
          eecPpMCBin = std::make_tuple(iWeightExponent, iJetPt, iTrackPt, iMCType);

          // Ask the histogram managers to deliver the proper histograms
          hEnergyEnergyCorrelatorPpMC[eecPpMCBin] = ppMCHistograms[iMCType][iWeightExponent]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatchedPp, iTrackPtMatchedPp, EECHistograms::kSameJetPair);

          // Normalize the histogram to the analysis region
          lowAnalysisBin = hEnergyEnergyCorrelatorPpMC[eecPpMCBin]->GetXaxis()->FindBin(analysisDeltaR.first + epsilon);
          highAnalysisBin = hEnergyEnergyCorrelatorPpMC[eecPpMCBin]->GetXaxis()->FindBin(analysisDeltaR.second - epsilon);
          hEnergyEnergyCorrelatorPpMC[eecPpMCBin]->Scale(1.0 / hEnergyEnergyCorrelatorPpMC[eecPpMCBin]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

          // Take a ratio between data distribution and MC distribution disregarding the uncertainty in the data
          hEnergyEnergyCorrelatorPpMCToDataRatio[eecPpMCBin] = (TH1D*) hEnergyEnergyCorrelatorPpMC[eecPpMCBin]->Clone(Form("monteCarloRatioToDataPp%d%d%d%d", iWeightExponent, iJetPt, iTrackPt, iMCType));
          errorlessData = (TH1D*) energyEnergyCorrelatorSignalPp[eecPpBin]->Clone(Form("errorlessComparisonToMC%d%d%d%d", iWeightExponent, iJetPt, iTrackPt, iMCType));
          optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
          hEnergyEnergyCorrelatorPpMCToDataRatio[eecPpMCBin]->Divide(errorlessData);

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
        eecPpBin = std::make_tuple(iWeightExponent, iJetPt, iTrackPt);
        for(int iWake = 0; iWake < HybridModelHistogramManager::kWakeConfigurations; iWake++){

          eecPpMCBin = std::make_tuple(iWeightExponent, iJetPt, iTrackPt, iWake);

          // Load the pp histograms
          energyEnergyCorrelatorHybridModelPp[eecPpMCBin] = hybridHistograms->GetEnergyEnergyCorrelatorPp(jetPtBin, trackPtBin, weightExponentHybrid, iWake);

          // Create a histogrammified version of the TGraphErrors
          histogrammifiedHybridModelPp[eecPpMCBin] = optimusPrimeTheTransformer->HistogrammifyWithErrors(energyEnergyCorrelatorHybridModelPp[eecPpMCBin], energyEnergyCorrelatorSignalPp[eecPpBin]);

          // Now we can take a ratio between hybrid model prediction and data
          hybridModelToDataRatioPp[eecPpMCBin] = (TH1D*) histogrammifiedHybridModelPp[eecPpMCBin]->Clone(Form("hybridModelRatioToDataPp%d%d%d%d", iWeightExponent, iJetPt, iTrackPt, iWake));
          errorlessData = (TH1D*) energyEnergyCorrelatorSignalPp[eecPpBin]->Clone(Form("errorlessHybrid%d%d%d%d", iWeightExponent, iJetPt, iTrackPt, iWake));
          optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
          hybridModelToDataRatioPp[eecPpMCBin]->Divide(errorlessData);

          for(auto centralityBin : drawnCentralityBin){
            iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);
            eecPbPbBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt);
            eecPbPbMCBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt, iWake);

            // Load the PbPb and ratio histograms
            energyEnergyCorrelatorHybridModelPbPb[eecPbPbMCBin] = hybridHistograms->GetEnergyEnergyCorrelatorPbPb(centralityBin, jetPtBin, trackPtBin, weightExponentHybrid, iWake);
            energyEnergyCorrelatorHybridModelPbPbToPpRatio[eecPbPbMCBin] = hybridHistograms->GetEnergyEnergyCorrelatorPbPbToPpRatio(centralityBin, jetPtBin, trackPtBin, weightExponentHybrid, iWake);

            // Create a histogrammified versions of the TGraphErrors
            histogrammifiedHybridModelPbPb[eecPbPbMCBin] = optimusPrimeTheTransformer->HistogrammifyWithErrors(energyEnergyCorrelatorHybridModelPbPb[eecPbPbMCBin], energyEnergyCorrelatorSignalPbPb[eecPbPbBin]);
            histogrammifiedHybridModelPbPbToPpRatio[eecPbPbMCBin] = optimusPrimeTheTransformer->HistogrammifyWithErrors(energyEnergyCorrelatorHybridModelPbPbToPpRatio[eecPbPbMCBin], energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]);

            // Now we can take a ratio between hybrid model prediction and data
            hybridModelToDataRatioPbPb[eecPbPbMCBin] = (TH1D*) histogrammifiedHybridModelPbPb[eecPbPbMCBin]->Clone(Form("hybridModelRatioToDataPbPb%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iWake));
            errorlessData = (TH1D*) energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->Clone(Form("errorlessHybridPbPb%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iWake));
            optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
            hybridModelToDataRatioPbPb[eecPbPbMCBin]->Divide(errorlessData);

            hybridModelToDataRatioPbPbToPpRatio[eecPbPbMCBin] = (TH1D*) histogrammifiedHybridModelPbPbToPpRatio[eecPbPbMCBin]->Clone(Form("hybridModelRatioToDataPbPbToPpRatio%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iWake));
            errorlessData = (TH1D*) energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->Clone(Form("errorlessHybridRatio%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iWake));
            optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
            hybridModelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->Divide(errorlessData);
        
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
          eecDoubleRatioBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt);
          for(int iWake = 0; iWake < HybridModelHistogramManager::kWakeConfigurations; iWake++){
            eecDoubleRatioMCBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iWake);
 
            // First, calculate the double ratio in Hybrid model
            histogrammifiedHybridModelDoubleRatio[eecDoubleRatioMCBin] = (TH1D*) histogrammifiedHybridModelPbPbToPpRatio[std::make_tuple(iWeightExponent,iCentrality,iJetPt,trackPtBinsForDoubleRatio.second,iWake)]->Clone(Form("energyEnergyCorrelatorDoubleRatioHybrid%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iWake));
            histogrammifiedHybridModelDoubleRatio[eecDoubleRatioMCBin]->Divide(histogrammifiedHybridModelPbPbToPpRatio[std::make_tuple(iWeightExponent,iCentrality,iJetPt,trackPtBinsForDoubleRatio.first,iWake)]);

            // Then, take a ratio between the double ratio in Hybrid model and double ratio in data
            hybridModelToDataRatioDoubleRatio[eecDoubleRatioMCBin] = (TH1D*) histogrammifiedHybridModelDoubleRatio[eecDoubleRatioMCBin]->Clone(Form("energyEnergyCorrelatorDoubleRatioDataHybridComparison%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iWake));
            errorlessData = (TH1D*) energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin]->Clone(Form("errorlessComparisonToHybridDoubleRatio%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iWake));
            optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
            hybridModelToDataRatioDoubleRatio[eecDoubleRatioMCBin]->Divide(errorlessData);

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
            eecPbPbBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt);
            eecPbPbMCBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt, iKValue);

            // Load the PbPb and ratio histograms
            energyEnergyCorrelatorHolguinPbPb[eecPbPbMCBin] = holguinHistograms->GetEnergyEnergyCorrelatorPbPb(centralityBin, jetPtBin, trackPtBin, iWeightExponent+1, myKValue);
            energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin] = holguinHistograms->GetEnergyEnergyCorrelatorPbPbToPpRatio(centralityBin, jetPtBin, trackPtBin, iWeightExponent+1, myKValue);

            // There are only a few bins for which predictions exist. Before continuing, check if the histograms are NULL
            if(energyEnergyCorrelatorHolguinPbPb[eecPbPbMCBin] == NULL) continue;

            // Create a histogrammified versions of the TGraphs
            histogrammifiedHolguinPbPb[eecPbPbMCBin] = optimusPrimeTheTransformer->Histogrammify(energyEnergyCorrelatorHolguinPbPb[eecPbPbMCBin], energyEnergyCorrelatorSignalPbPb[eecPbPbBin]);
            histogrammifiedHolguinPbPbToPpRatio[eecPbPbMCBin] = optimusPrimeTheTransformer->Histogrammify(energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin], energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]);

            // The absolute normalization for these theory predictions needs to be obtained from data
            // This is done by matching the integral in region 0.0419615 < DeltaR < 0.125882

            // Calculate the integrals in the data
            dataIntegralDistribution = energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->Integral(lowNormalizationBin, highNormalizationBin, "width");
            dataIntegralRatio = energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->Integral(lowNormalizationBin, highNormalizationBin, "width");

            // Calculate the integrals for the histogrammified distributions
            theoryIntegralDistribution = histogrammifiedHolguinPbPb[eecPbPbMCBin]->Integral(lowNormalizationBin, highNormalizationBin, "width");
            theoryIntegralRatio = histogrammifiedHolguinPbPbToPpRatio[eecPbPbMCBin]->Integral(lowNormalizationBin, highNormalizationBin, "width");

            // With these numbers, we can normalize the perturbative calculations
            energyEnergyCorrelatorHolguinPbPb[eecPbPbMCBin]->Scale(dataIntegralDistribution/theoryIntegralDistribution);
            histogrammifiedHolguinPbPb[eecPbPbMCBin]->Scale(dataIntegralDistribution/theoryIntegralDistribution);
            energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin]->Scale(dataIntegralRatio/theoryIntegralRatio);
            histogrammifiedHolguinPbPbToPpRatio[eecPbPbMCBin]->Scale(dataIntegralRatio/theoryIntegralRatio);

            // Now we can take a ratio between Holguin prediction and data
            holguinToDataRatioPbPb[eecPbPbMCBin] = (TH1D*) histogrammifiedHolguinPbPb[eecPbPbMCBin]->Clone(Form("holguinRatioToDataPbPb%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iKValue));
            errorlessData = (TH1D*) energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->Clone(Form("errorlessHolguinPbPb%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iKValue));
            optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
            holguinToDataRatioPbPb[eecPbPbMCBin]->Divide(errorlessData);

            holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin] = (TH1D*) histogrammifiedHolguinPbPbToPpRatio[eecPbPbMCBin]->Clone(Form("holguinRatioToDataPbPbToPpRatio%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iKValue));
            errorlessData = (TH1D*) energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->Clone(Form("errorlessHolguinRatio%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iKValue));
            optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
            holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin]->Divide(errorlessData);
        
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
          eecPbPbBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt);
          for(auto qValue : coLBTqValue){
            iQValue = coLBTHistograms->FindQValueIndex(qValue);
            if(iQValue < 0) continue;
            eecPbPbMCBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt, iQValue);

            // Load the PbPb to pp ratio histograms
            energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin] = coLBTHistograms->GetEnergyEnergyCorrelatorPbPbToPpRatio(centralityBin, jetPtBin, trackPtBin, iWeightExponent+1, qValue);

            // There are only a few bins for which predictions exist. Before continuing, check if the histograms are NULL
            if(energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin] == NULL) continue;

            // Create a histogrammified versions of the TGraphErrors
            histogrammifiedCoLBTPbPbToPpRatio[eecPbPbMCBin] = optimusPrimeTheTransformer->HistogrammifyWithErrors(energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin], energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]);

            // Now we can take a ratio between CoLBT prediction and data
            coLBTToDataRatioPbPbToPpRatio[eecPbPbMCBin] = (TH1D*) histogrammifiedCoLBTPbPbToPpRatio[eecPbPbMCBin]->Clone(Form("coLBTRatioToDataPbPbToPpRatio%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iQValue));
            errorlessData = (TH1D*) energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->Clone(Form("errorlessCoLBTRatio%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iQValue));
            optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
            coLBTToDataRatioPbPbToPpRatio[eecPbPbMCBin]->Divide(errorlessData);
        
          } // Q-vaule loop
        } // Centrality loop
      } // Track pT loop
    } // Jet pT loop
  } // Weight exponent loop

  // Read the histograms for JEWEL
  JewelHistogramManager* jewelHistograms = new JewelHistogramManager(jewelDataFolder);
  if(theoryComparisonIndex == 5 || theoryComparisonIndex == 6 || theoryComparisonIndex == 7 || drawnPlots == 0){
    for(int iWeightExponent = 0; iWeightExponent < nWeightExponents; iWeightExponent++){
      for(auto jetPtBin : drawnJetPtBin){
        iJetPt = card[kPbPb][iWeightExponent]->FindBinIndexJetPtEEC(jetPtBin);
        for(auto trackPtBin : drawnTrackPtBin){
          iTrackPt = card[kPbPb][iWeightExponent]->GetBinIndexTrackPtEEC(trackPtBin);
          eecPpBin = std::make_tuple(iWeightExponent, iJetPt, iTrackPt);

          // Load the pp histograms
          energyEnergyCorrelatorJewelPp[eecPpBin] = jewelHistograms->GetEnergyEnergyCorrelatorPp(jetPtBin, trackPtBin, iWeightExponent+1);

          // There are only a few bins for which predictions exist. Before continuing, check if the histograms are NULL
          if(energyEnergyCorrelatorJewelPp[eecPpBin] == NULL) continue;

          // Now we can take a ratio between JEWEL prediction and data
          jewelToDataRatioPp[eecPpBin] = (TH1D*) energyEnergyCorrelatorJewelPp[eecPpBin]->Clone(Form("jewelRatioToDataPp%d%d%d", iWeightExponent, iJetPt, iTrackPt));
          errorlessData = (TH1D*) energyEnergyCorrelatorSignalPp[eecPpBin]->Clone(Form("errorlessJewel%d%d%d", iWeightExponent, iJetPt, iTrackPt));
          optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
          jewelToDataRatioPp[eecPpBin]->Divide(errorlessData);

          for(auto centralityBin : drawnCentralityBin){
            iCentrality = card[kPbPb][iWeightExponent]->FindBinIndexCentrality(centralityBin);
            eecPbPbBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt);

            for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){
              eecPbPbMCBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iTrackPt, iRecoil);

              // Load the PbPb and ratio histograms
              energyEnergyCorrelatorJewelPbPb[eecPbPbMCBin] = jewelHistograms->GetEnergyEnergyCorrelatorPbPb(centralityBin, jetPtBin, trackPtBin, iWeightExponent+1, iRecoil);
              energyEnergyCorrelatorJewelPbPbToPpRatio[eecPbPbMCBin] = jewelHistograms->GetEnergyEnergyCorrelatorPbPbToPpRatio(centralityBin, jetPtBin, trackPtBin, iWeightExponent+1, iRecoil);

              // The prediction might not exist in all bins
              if(energyEnergyCorrelatorJewelPbPb[eecPbPbMCBin] == NULL) continue;

              // Take a ratio between JEWEL prediction and data
              jewelToDataRatioPbPb[eecPbPbMCBin] = (TH1D*) energyEnergyCorrelatorJewelPbPb[eecPbPbMCBin]->Clone(Form("jewelRatioToDataPbPb%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iRecoil));
              errorlessData = (TH1D*) energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->Clone(Form("errorlessJewelPbPb%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iRecoil));
              optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
              jewelToDataRatioPbPb[eecPbPbMCBin]->Divide(errorlessData);

              jewelToDataRatioPbPbToPpRatio[eecPbPbMCBin] = (TH1D*) energyEnergyCorrelatorJewelPbPbToPpRatio[eecPbPbMCBin]->Clone(Form("jewelRatioToDataPbPbToPpRatio%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iRecoil));
              errorlessData = (TH1D*) energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->Clone(Form("errorlessJewelRatio%d%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iTrackPt, iRecoil));
              optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
              jewelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->Divide(errorlessData);
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
            eecDoubleRatioBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt);
            for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){
              eecDoubleRatioMCBin = std::make_tuple(iWeightExponent, iCentrality, iJetPt, iRecoil);
 
              // First, calculate the double ratio in JEWEL
              jewelDoubleRatio[eecDoubleRatioMCBin] = (TH1D*) energyEnergyCorrelatorJewelPbPbToPpRatio[std::make_tuple(iWeightExponent,iCentrality,iJetPt,trackPtBinsForDoubleRatio.second,iRecoil)]->Clone(Form("energyEnergyCorrelatorDoubleRatioJewel%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iRecoil));
              jewelDoubleRatio[eecDoubleRatioMCBin]->Divide(energyEnergyCorrelatorJewelPbPbToPpRatio[std::make_tuple(iWeightExponent,iCentrality,iJetPt,trackPtBinsForDoubleRatio.first,iRecoil)]);

              // Then, take a ratio between the double ratio in JEWEL and double ratio in data
              jewelToDataRatioDoubleRatio[eecDoubleRatioMCBin] = (TH1D*) jewelDoubleRatio[eecDoubleRatioMCBin]->Clone(Form("energyEnergyCorrelatorDoubleRatioDataJewelComparison%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iRecoil));
              errorlessData = (TH1D*) energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin]->Clone(Form("errorlessComparisonToJewelDoubleRatio%d%d%d%d", iWeightExponent, iCentrality, iJetPt, iRecoil));
              optimusPrimeTheTransformer->RemoveUncertainties(errorlessData);
              jewelToDataRatioDoubleRatio[eecDoubleRatioMCBin]->Divide(errorlessData);

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
        eecPpBin = std::make_tuple(weightExponent-1, iJetPt, iTrackPt);
        eecPpUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iJetPt, iTrackPt);
        eecPpCorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kCorrelatedUncertainty, iJetPt, iTrackPt);
        eecPpRelativeUncertaintyBinSystematic = std::make_tuple(weightExponent-1, iJetPt, iTrackPt, kRelativeUncertaintySystematic);
        eecPpRelativeUncertaintyBinStatisticalUp = std::make_tuple(weightExponent-1, iJetPt, iTrackPt, kRelativeUncertaintyStatisticalUp);
        eecPpRelativeUncertaintyBinStatisticalDown = std::make_tuple(weightExponent-1, iJetPt, iTrackPt, kRelativeUncertaintyStatisticalDown);

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
        energyEnergyCorrelatorSignalPp[eecPpBin]->SetMarkerStyle(kFullSquare);
        energyEnergyCorrelatorSignalPp[eecPpBin]->SetMarkerColor(kBlack);

        systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetFillColorAlpha(kBlack, 0.4);
        systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetLineColor(kBlack);
        systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetMarkerColor(kBlack);
        systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->SetMarkerStyle(kFullSquare);

        systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->SetFillColorAlpha(kBlack, 0.4);

        // Set the x-axis drawing range
        systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
        systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetYaxis()->SetRangeUser(distributionZoom[weightExponent-1].first, distributionZoom[weightExponent-1].second);

        // Adjust tick size for x-axis
        systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin]->GetXaxis()->SetTickLength(tickSizeX);

        // Draw the data correlator to upper canves
        drawer->DrawHistogramToUpperPad(systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin], "#Deltar", "EEC", " ", "e2");
        systematicUncertaintyForPp[eecPpCorrelatedSystematicsBin]->Draw("same,e3");
        energyEnergyCorrelatorSignalPp[eecPpBin]->Draw("same,p");

        legend->AddEntry(systematicUncertaintyForPp[eecPpUncorrelatedSystematicsBin], "pp data", "lpf");

        // Pythia+Herwig+Hybrid

        // Add legend for Pythia8 and Herwig7, unless we are drawing only JEWEL
        if(theoryComparisonIndex != 5){
          legend->AddEntry(hEnergyEnergyCorrelatorPpMC[std::make_tuple(weightExponent-1,iJetPt,iTrackPt,kPythia)], "Pythia8 CP5", "pl");
          legend->AddEntry(hEnergyEnergyCorrelatorPpMC[std::make_tuple(weightExponent-1,iJetPt,iTrackPt,kHerwig)], "Herwig7 CH3", "pl");
    
          // There is no wake in pp, so do only comparison with Hybrid without wake
          eecPpMCBin = std::make_tuple(weightExponent-1, iJetPt, iTrackPt, 0);
          if(energyEnergyCorrelatorHybridModelPp[eecPpMCBin] != NULL){

            // Give some nice styles for the predictions
            energyEnergyCorrelatorHybridModelPp[eecPpMCBin]->SetLineColor(ppModelColor[kPpCompareHybrid]);
            energyEnergyCorrelatorHybridModelPp[eecPpMCBin]->SetLineWidth(0);
            energyEnergyCorrelatorHybridModelPp[eecPpMCBin]->SetMarkerColor(ppModelColor[kPpCompareHybrid]);
            energyEnergyCorrelatorHybridModelPp[eecPpMCBin]->SetMarkerStyle(kFullCircle);
            energyEnergyCorrelatorHybridModelPp[eecPpMCBin]->SetFillColorAlpha(ppModelColor[kPpCompareHybrid], 0.4);

            // Draw the prediction to the same canvas as the data
            energyEnergyCorrelatorHybridModelPp[eecPpMCBin]->Draw("3,same");

            // Add a legend for the theory prediction
            legend->AddEntry(energyEnergyCorrelatorHybridModelPp[eecPpMCBin], "Hybrid model", "f");
          }

          // After hybrid model, add also Pythia8 and Herwig7 predictions
          for(int iMCType = 0; iMCType < knPpMCTypes; iMCType++){

            eecPpMCBin = std::make_tuple(weightExponent-1, iJetPt, iTrackPt, iMCType);

            // Give some nice styles for the predictions
            hEnergyEnergyCorrelatorPpMC[eecPpMCBin]->SetLineColor(ppModelColor[ppMCIndexToModelColorMap[iMCType]]);
            hEnergyEnergyCorrelatorPpMC[eecPpMCBin]->SetMarkerColor(ppModelColor[ppMCIndexToModelColorMap[iMCType]]);
            hEnergyEnergyCorrelatorPpMC[eecPpMCBin]->SetMarkerStyle(ppMCMarkerStyle[iMCType]);
            hEnergyEnergyCorrelatorPpMC[eecPpMCBin]->SetFillColorAlpha(ppModelColor[ppMCIndexToModelColorMap[iMCType]], 0.4);

            // Draw the prediction to the same canvas as the data
            hEnergyEnergyCorrelatorPpMC[eecPpMCBin]->Draw("same,p");
          }
        }

        if(theoryComparisonIndex >= 5){
          // Option to also draw JEWEL in addition to other vacuum models

          // Set a nice style for the JEWEL histograms
          if(theoryComparisonIndex == 5){
            energyEnergyCorrelatorJewelPp[eecPpBin]->SetLineColor(kRed);
            energyEnergyCorrelatorJewelPp[eecPpBin]->SetMarkerColor(kRed);
            energyEnergyCorrelatorJewelPp[eecPpBin]->SetMarkerStyle(kFullCircle);
          } else {
            energyEnergyCorrelatorJewelPp[eecPpBin]->SetLineColor(kMagenta);
            energyEnergyCorrelatorJewelPp[eecPpBin]->SetMarkerColor(kMagenta);
            energyEnergyCorrelatorJewelPp[eecPpBin]->SetMarkerStyle(kFullCross);
          }

          // Draw the jewel histogram and add it to the legend
          energyEnergyCorrelatorJewelPp[eecPpBin]->Draw("same,pl");
          legend->AddEntry(energyEnergyCorrelatorJewelPp[eecPpBin], "JEWEL", "pl");
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
        hRelativeUncertaintyPp[eecPpRelativeUncertaintyBinSystematic]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
        hRelativeUncertaintyPp[eecPpRelativeUncertaintyBinSystematic]->GetYaxis()->SetRangeUser(ratioZoomPpDistribution.first, ratioZoomPpDistribution.second);

        // Adjust tick sizes for the x-axis
        hRelativeUncertaintyPp[eecPpRelativeUncertaintyBinSystematic]->GetXaxis()->SetTickLength(tickSizeX * 1.5);

        // Set the style for uncertainty bands for systematic and statistical uncertainties from data
        for(int iUncertainty = 0; iUncertainty < knRelativeUncertaintyTypes; iUncertainty++){
          eecPpRelativeUncertaintyBin = std::make_tuple(weightExponent-1, iJetPt, iTrackPt, iUncertainty);
          hRelativeUncertaintyPp[eecPpRelativeUncertaintyBin]->SetLineColor(kBlack);
          hRelativeUncertaintyPp[eecPpRelativeUncertaintyBin]->SetLineStyle(9);
          hRelativeUncertaintyPp[eecPpRelativeUncertaintyBin]->SetMarkerStyle(kFullCircle);
          hRelativeUncertaintyPp[eecPpRelativeUncertaintyBin]->SetMarkerSize(0);
        }
        hRelativeUncertaintyPp[eecPpRelativeUncertaintyBinSystematic]->SetFillColorAlpha(kBlack, 0.4);
        hRelativeUncertaintyPp[eecPpRelativeUncertaintyBinSystematic]->SetLineWidth(0);

        // Create a new legend to show the different data uncertainty bands
        legendX1 = 0.25; legendX2 = 0.95;
        legendY1 = 0.86; legendY2 = 0.96;
        legend = new TLegend(legendX1, legendY1, legendX2, legendY2);
        legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(legendTextLowerCanvasSize); legend->SetTextFont(62);
        legend->SetNColumns(2);

        // Draw the error bars from data
        drawer->DrawHistogramToLowerPad(hRelativeUncertaintyPp[eecPpRelativeUncertaintyBinSystematic], "#Deltar", "#frac{Theory}{Data}", " ", "e3");
        hRelativeUncertaintyPp[eecPpRelativeUncertaintyBinStatisticalUp]->Draw("same,HIST,C");
        hRelativeUncertaintyPp[eecPpRelativeUncertaintyBinStatisticalDown]->Draw("same,HIST,C");

        // Add the different uncertainties to the legend
        legend->AddEntry(hRelativeUncertaintyPp[eecPpRelativeUncertaintyBinSystematic], "Data syst. unc.", "f");
        legend->AddEntry(hRelativeUncertaintyPp[eecPpRelativeUncertaintyBinStatisticalUp], "Data stat. unc.", "l");

        // Draw Pythia, Herwig, and Hybrid unless we require only JEWEL to be drawn
        if(theoryComparisonIndex != 5){

          eecPpMCBin = std::make_tuple(weightExponent-1, iJetPt, iTrackPt, 0);

          // Set the style for MC to data comparison histograms
          hybridModelToDataRatioPp[eecPpMCBin]->SetMarkerStyle(kFullCircle);
          hybridModelToDataRatioPp[eecPpMCBin]->SetMarkerSize(0);
          hybridModelToDataRatioPp[eecPpMCBin]->SetFillColorAlpha(ppModelColor[kPpCompareHybrid], 0.4);

          // Set the style for pp MC simulation to data ratios
          for(int iMCType = 0; iMCType < knPpMCTypes; iMCType++){
            eecPpMCBin = std::make_tuple(weightExponent-1, iJetPt, iTrackPt, iMCType);

            hEnergyEnergyCorrelatorPpMCToDataRatio[eecPpMCBin]->SetMarkerStyle(ppMCMarkerStyle[iMCType]);
            hEnergyEnergyCorrelatorPpMCToDataRatio[eecPpMCBin]->SetMarkerColor(ppModelColor[ppMCIndexToModelColorMap[iMCType]]);
            hEnergyEnergyCorrelatorPpMCToDataRatio[eecPpMCBin]->SetLineColor(ppModelColor[ppMCIndexToModelColorMap[iMCType]]);
            hEnergyEnergyCorrelatorPpMCToDataRatio[eecPpMCBin]->SetFillColorAlpha(ppModelColor[ppMCIndexToModelColorMap[iMCType]], 0.4);
          }

          eecPpMCBin = std::make_tuple(weightExponent-1, iJetPt, iTrackPt, 0);

          // Draw the ratio with respect to hybrid model
          hybridModelToDataRatioPp[eecPpMCBin]->Draw("same,e3");

          // Draw the ratio with respect to Pythia and Herwig simulations
          hEnergyEnergyCorrelatorPpMCToDataRatio[std::make_tuple(weightExponent-1,iJetPt,iTrackPt,kPythia)]->Draw("same,p");
          hEnergyEnergyCorrelatorPpMCToDataRatio[std::make_tuple(weightExponent-1,iJetPt,iTrackPt,kHerwig)]->Draw("same,p");

        }

        if(theoryComparisonIndex >= 5){
          // Theory comparison index above 5: add JEWEL vacuum

          // Set the style for the histograms
          if(theoryComparisonIndex == 5){
            jewelToDataRatioPp[eecPpBin]->SetMarkerStyle(kFullCircle);
            jewelToDataRatioPp[eecPpBin]->SetMarkerColor(kRed);
            jewelToDataRatioPp[eecPpBin]->SetLineColor(kRed);
          } else {
            jewelToDataRatioPp[eecPpBin]->SetMarkerStyle(kFullCross);
            jewelToDataRatioPp[eecPpBin]->SetMarkerColor(kMagenta);
            jewelToDataRatioPp[eecPpBin]->SetLineColor(kMagenta);
          }

          // Draw the histogram to the same canvas as data
          jewelToDataRatioPp[eecPpBin]->Draw("same,pl");
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
          eecPbPbBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt);
          eecPbPbUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeUp, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeDown, iCentrality, iJetPt, iTrackPt);
          eecPbPbRelativeUncertaintyBinSystematic = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, kRelativeUncertaintySystematic);
          eecPbPbRelativeUncertaintyBinStatisticalUp = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, kRelativeUncertaintyStatisticalUp);
          eecPbPbRelativeUncertaintyBinStatisticalDown = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, kRelativeUncertaintyStatisticalDown);

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
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerStyle(kFullSquare);
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetMarkerColor(kBlack);
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->SetLineColor(kBlack);

          // Set drawing style for systematic uncertainties
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetFillColorAlpha(kBlack, 0.4);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetLineColor(kBlack);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetMarkerColor(kBlack);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->SetMarkerStyle(kFullSquare);

          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->SetFillColorAlpha(kGray+3, 0.4);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->SetFillColorAlpha(kBlue+4, 0.4);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerSize(0);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerStyle(9);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerSize(0);
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerStyle(9);

          // Set the axis drawing ranges
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetRangeUser(distributionZoom[weightExponent-1].first, distributionZoom[weightExponent-1].second);

          // Adjust tick size for x-axis
          systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTickLength(tickSizeX);

          // Draw the systematic uncertainties to the upper canvas
          drawer->DrawHistogramToUpperPad(systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin], "#Deltar", "EEC", " ", "e2");
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeUp]->Draw("same,e3");
          systematicUncertaintyForPbPb[eecPbPbCorrelatedSystematicsBinShapeDown]->Draw("same,e3");

          // Draw the data points to the same canvas and add the histogram to the legend
          energyEnergyCorrelatorSignalPbPb[eecPbPbBin]->Draw("same,p");

          legend->AddEntry(systematicUncertaintyForPbPb[eecPbPbUncorrelatedSystematicsBin], centralityString.Data(), "lpf");

          // Choose a set of theory predictions to draw to the figure
          if(theoryComparisonIndex == 0){
            // Theory comparison index 0, Hybrid model with different wake settings

            // Compare the prediction with and without wake
            for(int iWake = 0; iWake < HybridModelHistogramManager::kWakeConfigurations; iWake++){
              eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iWake);

              // There are some bins for which the prediction does not exist
              if(energyEnergyCorrelatorHybridModelPbPb[eecPbPbMCBin] == NULL) continue;

              // Give some nice styles for the predictions
              energyEnergyCorrelatorHybridModelPbPb[eecPbPbMCBin]->SetLineColor(color[iWake]);
              energyEnergyCorrelatorHybridModelPbPb[eecPbPbMCBin]->SetMarkerColor(color[iWake]);
              energyEnergyCorrelatorHybridModelPbPb[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);
              energyEnergyCorrelatorHybridModelPbPb[eecPbPbMCBin]->SetFillColorAlpha(color[iWake], jetWakeOpacity[iWake]);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorHybridModelPbPb[eecPbPbMCBin]->Draw("3,same");

              // Add a legend for the theory prediction
              legend->AddEntry(energyEnergyCorrelatorHybridModelPbPb[eecPbPbMCBin], hybridHistograms->GetWakeName(iWake), "f");
            }
          } else if(theoryComparisonIndex == 1){
            // Theory comparison index 1, perturbative calculations from Holguin et. al. with different k-values

            // Compare the prediction with different k-values
            styleIndex = 0;
            for(auto myKValue : holguinKValue){
              iKValue = holguinHistograms->FindKValueIndex(myKValue);
              eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iKValue);

              // There are some bins for which the prediction does not exist
              if(energyEnergyCorrelatorHolguinPbPb[eecPbPbMCBin] == NULL) continue;

              // Give some nice styles for the predictions
              energyEnergyCorrelatorHolguinPbPb[eecPbPbMCBin]->SetLineColor(color[iKValue]);
              energyEnergyCorrelatorHolguinPbPb[eecPbPbMCBin]->SetMarkerColor(color[iKValue]);
              energyEnergyCorrelatorHolguinPbPb[eecPbPbMCBin]->SetMarkerStyle(markerStyleTheory[styleIndex++]);
              energyEnergyCorrelatorHolguinPbPb[eecPbPbMCBin]->SetFillColorAlpha(color[iKValue], 0.4);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorHolguinPbPb[eecPbPbMCBin]->Draw("pl,same");

              // Add a legend for the theory prediction
              legend->AddEntry(energyEnergyCorrelatorHolguinPbPb[eecPbPbMCBin], Form("Holguin, k=%.1f", myKValue), "pl");
            }
          } else if(theoryComparisonIndex == 5){
            // Theory comparison index 5, comparison to JEWEL

            // Compare the prediction with and without recoils
            for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){
              eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iRecoil);

              // There might be some bins for which the prediction does not exist
              if(energyEnergyCorrelatorJewelPbPb[eecPbPbMCBin] == NULL) continue;

              // Give some nice styles for the predictions
              energyEnergyCorrelatorJewelPbPb[eecPbPbMCBin]->SetLineColor(color[iRecoil]);
              energyEnergyCorrelatorJewelPbPb[eecPbPbMCBin]->SetMarkerColor(color[iRecoil]);
              energyEnergyCorrelatorJewelPbPb[eecPbPbMCBin]->SetMarkerStyle(jewelMarkerStyle[iRecoil]);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorJewelPbPb[eecPbPbMCBin]->Draw("pl,same");

              // Add a legend for the theory prediction
              legend->AddEntry(energyEnergyCorrelatorJewelPbPb[eecPbPbMCBin], jewelHistograms->GetRecoilName(iRecoil), "pl");
            }
          } else if(theoryComparisonIndex == 6){
            // Theory comparison index 6, comparison to best k-value from Holguin and JEWEL

            // Add the perturbative calculation with k=0.3
            iKValue = holguinHistograms->FindKValueIndex(0.3);
            eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iKValue);

            // There are some bins for which the prediction does not exist
            if(energyEnergyCorrelatorHolguinPbPb[eecPbPbMCBin] != NULL){

              // Give some nice styles for the predictions
              energyEnergyCorrelatorHolguinPbPb[eecPbPbMCBin]->SetLineColor(color[2]);
              energyEnergyCorrelatorHolguinPbPb[eecPbPbMCBin]->SetMarkerColor(color[2]);
              energyEnergyCorrelatorHolguinPbPb[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorHolguinPbPb[eecPbPbMCBin]->Draw("pl,same");

              // Add a legend for the theory prediction
              legend->AddEntry(energyEnergyCorrelatorHolguinPbPb[eecPbPbMCBin], "Holguin, k=0.3", "pl");
            }

            // Add the JEWEL predictions with and without recoil to the same plot

            for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){
              eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iRecoil);

              // There might be some bins for which the prediction does not exist
              if(energyEnergyCorrelatorJewelPbPb[eecPbPbMCBin] == NULL) continue;

              // Give some nice styles for the predictions
              energyEnergyCorrelatorJewelPbPb[eecPbPbMCBin]->SetLineColor(color[iRecoil]);
              energyEnergyCorrelatorJewelPbPb[eecPbPbMCBin]->SetMarkerColor(color[iRecoil]);
              energyEnergyCorrelatorJewelPbPb[eecPbPbMCBin]->SetMarkerStyle(jewelMarkerStyle[iRecoil]);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorJewelPbPb[eecPbPbMCBin]->Draw("pl,same");

              // Add a legend for the theory prediction
              legend->AddEntry(energyEnergyCorrelatorJewelPbPb[eecPbPbMCBin], jewelHistograms->GetRecoilName(iRecoil), "pl");
            }
          } else if(theoryComparisonIndex == 7){

            // Draw the hybrid model with all bells and whistles
            eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, HybridModelHistogramManager::kFullWake);
            if(energyEnergyCorrelatorHybridModelPbPb[eecPbPbMCBin] != NULL){

              // Give some nice styles for the predictions
              energyEnergyCorrelatorHybridModelPbPb[eecPbPbMCBin]->SetLineColor(color[HybridModelHistogramManager::kFullWake]);
              energyEnergyCorrelatorHybridModelPbPb[eecPbPbMCBin]->SetMarkerColor(color[HybridModelHistogramManager::kFullWake]);
              energyEnergyCorrelatorHybridModelPbPb[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);
              energyEnergyCorrelatorHybridModelPbPb[eecPbPbMCBin]->SetFillColorAlpha(color[HybridModelHistogramManager::kFullWake], jetWakeOpacity[HybridModelHistogramManager::kFullWake]);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorHybridModelPbPb[eecPbPbMCBin]->Draw("3,same");

              // Add a legend for the theory prediction
              legend->AddEntry(energyEnergyCorrelatorHybridModelPbPb[eecPbPbMCBin], hybridHistograms->GetWakeName(HybridModelHistogramManager::kFullWake), "f");
            }

            // Draw the JEWEL with recoils
            eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, JewelHistogramManager::kRecoil);
            if(energyEnergyCorrelatorJewelPbPb[eecPbPbMCBin] != NULL){

              // Give some nice styles for the predictions
              energyEnergyCorrelatorJewelPbPb[eecPbPbMCBin]->SetLineColor(color[JewelHistogramManager::kRecoil]);
              energyEnergyCorrelatorJewelPbPb[eecPbPbMCBin]->SetMarkerColor(color[JewelHistogramManager::kRecoil]);
              energyEnergyCorrelatorJewelPbPb[eecPbPbMCBin]->SetMarkerStyle(jewelMarkerStyle[JewelHistogramManager::kRecoil]);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorJewelPbPb[eecPbPbMCBin]->Draw("pl,same");

              // Add a legend for the theory prediction
              legend->AddEntry(energyEnergyCorrelatorJewelPbPb[eecPbPbMCBin], jewelHistograms->GetRecoilName(JewelHistogramManager::kRecoil), "pl");
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
          hRelativeUncertaintyPbPb[eecPbPbRelativeUncertaintyBinSystematic]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          hRelativeUncertaintyPbPb[eecPbPbRelativeUncertaintyBinSystematic]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

          // Adjust tick size for x-axis
          hRelativeUncertaintyPbPb[eecPbPbRelativeUncertaintyBinSystematic]->GetXaxis()->SetTickLength(tickSizeX*1.5);

          // Set the style for uncertainty bands for systematic and statistical uncertainties from data
          for(int iUncertainty = 0; iUncertainty < knRelativeUncertaintyTypes; iUncertainty++){
            eecPbPbRelativeUncertaintyBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iUncertainty);
            hRelativeUncertaintyPbPb[eecPbPbRelativeUncertaintyBin]->SetLineColor(kBlack);
            hRelativeUncertaintyPbPb[eecPbPbRelativeUncertaintyBin]->SetLineStyle(9);
            hRelativeUncertaintyPbPb[eecPbPbRelativeUncertaintyBin]->SetMarkerStyle(kFullCircle);
            hRelativeUncertaintyPbPb[eecPbPbRelativeUncertaintyBin]->SetMarkerSize(0);
          }
          hRelativeUncertaintyPbPb[eecPbPbRelativeUncertaintyBinSystematic]->SetFillColorAlpha(kBlack, 0.4);
          hRelativeUncertaintyPbPb[eecPbPbRelativeUncertaintyBinSystematic]->SetLineWidth(0);


          // Create a new legend to show the different data uncertainty bands
          legendX1 = 0.25; legendX2 = 0.95;
          legendY1 = 0.86; legendY2 = 0.96;
          legend = new TLegend(legendX1, legendY1, legendX2, legendY2);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(legendTextLowerCanvasSize); legend->SetTextFont(62);
          legend->SetNColumns(2);

          // Draw the error bars from the data
          drawer->DrawHistogramToLowerPad(hRelativeUncertaintyPbPb[eecPbPbRelativeUncertaintyBinSystematic], "#Deltar", "#frac{Theory}{Data}", " ", "e3");
          hRelativeUncertaintyPbPb[eecPbPbRelativeUncertaintyBinStatisticalUp]->Draw("same,HIST,C");
          hRelativeUncertaintyPbPb[eecPbPbRelativeUncertaintyBinStatisticalDown]->Draw("same,HIST,C");

          // Add the different uncertainties to the legend
          legend->AddEntry(hRelativeUncertaintyPbPb[eecPbPbRelativeUncertaintyBinSystematic], "Data syst. unc.", "f");
          legend->AddEntry(hRelativeUncertaintyPbPb[eecPbPbRelativeUncertaintyBinStatisticalUp], "Data stat. unc.", "l");

          // Choose a set of theory predictions to draw to the figure
          if(theoryComparisonIndex == 0){
            // Theory comparison index 0, Hybrid model with different wake settings
            for(int iWake = 0; iWake < HybridModelHistogramManager::kWakeConfigurations; iWake++){
              eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iWake);
              hybridModelToDataRatioPbPb[eecPbPbMCBin]->SetMarkerColor(color[iWake]);
              hybridModelToDataRatioPbPb[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);
              hybridModelToDataRatioPbPb[eecPbPbMCBin]->SetMarkerSize(0);
              hybridModelToDataRatioPbPb[eecPbPbMCBin]->SetFillColorAlpha(color[iWake], jetWakeOpacity[iWake]);
              hybridModelToDataRatioPbPb[eecPbPbMCBin]->Draw("same,e3");
            }
          } else if(theoryComparisonIndex == 1){
            // Theory comparison index 1, perturbative calculations from Holguin et. al. with different k-values
            styleIndex = 0;
            for(auto myKValue : holguinKValue){
              iKValue = holguinHistograms->FindKValueIndex(myKValue);
              eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iKValue);
              if(holguinToDataRatioPbPb[eecPbPbMCBin] == NULL) continue;
              holguinToDataRatioPbPb[eecPbPbMCBin]->SetMarkerColor(color[iKValue]);
              holguinToDataRatioPbPb[eecPbPbMCBin]->SetLineWidth(3);
              holguinToDataRatioPbPb[eecPbPbMCBin]->SetLineColor(color[iKValue]);
              holguinToDataRatioPbPb[eecPbPbMCBin]->SetMarkerStyle(markerStyleTheory[styleIndex++]);
              holguinToDataRatioPbPb[eecPbPbMCBin]->Draw("same,HIST,C");
            }
          } else if(theoryComparisonIndex == 5){
            // Theory comparison index 5, JEWEL prediction

            for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){
              eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iRecoil);
              if(jewelToDataRatioPbPb[eecPbPbMCBin] == NULL) continue;
              jewelToDataRatioPbPb[eecPbPbMCBin]->SetMarkerColor(color[iRecoil]);
              jewelToDataRatioPbPb[eecPbPbMCBin]->SetMarkerStyle(jewelMarkerStyle[iRecoil]);
              jewelToDataRatioPbPb[eecPbPbMCBin]->Draw("same,pl");
            }
          } else if(theoryComparisonIndex == 6){
            // Theory comparison index 6, best Holguin k-value and JEWEL prediction

            // Draw the perturbative calculation
            iKValue = holguinHistograms->FindKValueIndex(0.3);
            eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iKValue);
            if(holguinToDataRatioPbPb[eecPbPbMCBin] != NULL){
              holguinToDataRatioPbPb[eecPbPbMCBin]->SetMarkerColor(color[2]);
              holguinToDataRatioPbPb[eecPbPbMCBin]->SetLineWidth(3);
              holguinToDataRatioPbPb[eecPbPbMCBin]->SetLineColor(color[2]);
              holguinToDataRatioPbPb[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);
              holguinToDataRatioPbPb[eecPbPbMCBin]->Draw("same,HIST,C");
            }

            // Add the JEWEL predictions to the same figure
            for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){
              eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iRecoil);
              if(jewelToDataRatioPbPb[eecPbPbMCBin] == NULL) continue;
              jewelToDataRatioPbPb[eecPbPbMCBin]->SetMarkerColor(color[iRecoil]);
              jewelToDataRatioPbPb[eecPbPbMCBin]->SetMarkerStyle(jewelMarkerStyle[iRecoil]);
              jewelToDataRatioPbPb[eecPbPbMCBin]->Draw("same,pl");
            }
          } else if(theoryComparisonIndex == 7){
            // Theory comparison index 7, full wake Hybrid and JEWEL with recoils

            eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, HybridModelHistogramManager::kFullWake);
            if(hybridModelToDataRatioPbPb[eecPbPbMCBin] != NULL){
              hybridModelToDataRatioPbPb[eecPbPbMCBin]->SetMarkerColor(color[HybridModelHistogramManager::kFullWake]);
              hybridModelToDataRatioPbPb[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);
              hybridModelToDataRatioPbPb[eecPbPbMCBin]->SetMarkerSize(0);
              hybridModelToDataRatioPbPb[eecPbPbMCBin]->SetFillColorAlpha(color[HybridModelHistogramManager::kFullWake], jetWakeOpacity[HybridModelHistogramManager::kFullWake]);
              hybridModelToDataRatioPbPb[eecPbPbMCBin]->Draw("same,e3");
            }

            eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, JewelHistogramManager::kRecoil);
            if(jewelToDataRatioPbPb[eecPbPbMCBin] != NULL){
              jewelToDataRatioPbPb[eecPbPbMCBin]->SetMarkerColor(color[JewelHistogramManager::kRecoil]);
              jewelToDataRatioPbPb[eecPbPbMCBin]->SetMarkerStyle(jewelMarkerStyle[JewelHistogramManager::kRecoil]);
              jewelToDataRatioPbPb[eecPbPbMCBin]->Draw("same,pl");
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
          eecPbPbBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt);
          eecPbPbUncorrelatedSystematicsBin = std::make_tuple(weightExponent-1, kUncorrelatedUncertainty, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeUp = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeUp, iCentrality, iJetPt, iTrackPt);
          eecPbPbCorrelatedSystematicsBinShapeDown = std::make_tuple(weightExponent-1, kCorrelatedUncertaintyShapeDown, iCentrality, iJetPt, iTrackPt);
          eecPbPbRelativeUncertaintyBinSystematic = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, kRelativeUncertaintySystematic);
          eecPbPbRelativeUncertaintyBinStatisticalUp = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, kRelativeUncertaintyStatisticalUp);
          eecPbPbRelativeUncertaintyBinStatisticalDown = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, kRelativeUncertaintyStatisticalDown);

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
          energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->SetMarkerStyle(kFullSquare);
          energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->SetMarkerColor(kBlack);
          energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->SetLineColor(kBlack);

          // Set drawing style for systematic uncertainties
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetFillColorAlpha(kBlack, 0.4);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetLineColor(kBlack);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetMarkerColor(kBlack);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->SetMarkerStyle(kFullSquare);

          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->SetFillColorAlpha(kGray+3, 0.4);
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->SetFillColorAlpha(kBlue+4, 0.4);
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerSize(0);
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->SetMarkerStyle(9);
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerSize(0);
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->SetMarkerStyle(9);

          // Set the x-axis drawing range
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetYaxis()->SetRangeUser(pbpbToPpRatioZoom.first, pbpbToPpRatioZoom.second);

          // Adjust tick size for x-axis
          systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin]->GetXaxis()->SetTickLength(tickSizeX);

          // Draw first the systematic uncertainties to the upper canves
          drawer->DrawHistogramToUpperPad(systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin], "#Deltar", "#frac{PbPb}{pp}", " ", "e2");
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeUp]->Draw("same,e3");
          systematicUncertaintyPbPbToPpRatio[eecPbPbCorrelatedSystematicsBinShapeDown]->Draw("same,e3");

          // Then draw the PbPb to pp ratio and add a legend for it
          energyEnergyCorrelatorPbPbToPpRatio[eecPbPbBin]->Draw("same,p");
          anotherLegend->AddEntry(systematicUncertaintyPbPbToPpRatio[eecPbPbUncorrelatedSystematicsBin], Form("%s / pp", centralityString.Data()), "lpf");

          // Choose a set of theory predictions to draw to the figure
          if(theoryComparisonIndex == 0){
            // Theory comparison index 0, Hybrid model with different wake settings

            for(int iWake = 0; iWake < HybridModelHistogramManager::kWakeConfigurations; iWake++){

              eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iWake);

              // There are some bins for which the prediction does not exist
              if(energyEnergyCorrelatorHybridModelPbPbToPpRatio[eecPbPbMCBin] == NULL) continue;

              // Give some nice styles for the predictions
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[eecPbPbMCBin]->SetLineColor(color[iWake]);
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(color[iWake]);
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[eecPbPbMCBin]->SetFillColorAlpha(color[iWake], jetWakeOpacity[iWake]);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[eecPbPbMCBin]->Draw("3,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorHybridModelPbPbToPpRatio[eecPbPbMCBin], hybridHistograms->GetWakeName(iWake), "f");
            }
          } else if (theoryComparisonIndex == 1){
            // Theory comparison index 1, perturbative calculations from Holguin et. al. with different k-values

            styleIndex = 0;
            for(auto myKValue : holguinKValue){
              iKValue = holguinHistograms->FindKValueIndex(myKValue);
              eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iKValue);

              // There are some bins for which the prediction does not exist
              if(energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin] == NULL) continue;

              // Give some nice styles for the predictions
              energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin]->SetLineColor(color[iKValue]);
              energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(color[iKValue]);
              energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(markerStyleTheory[styleIndex++]);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin]->Draw("pl,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin], Form("Holguin, k=%.1f", myKValue), "pl");
            }

          } else if(theoryComparisonIndex == 2){
            // Theory comparison index 2, CoLBT prediction

            for(auto myQValue : coLBTqValue){
              iQValue = coLBTHistograms->FindQValueIndex(myQValue);
              eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iQValue);

              // There are some bins for which the prediction does not exist
              if(energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin] == NULL) continue;

              // Give some nice styles for the predictions
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin]->SetLineColor(color[iQValue]);
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(color[iQValue]);
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin]->SetFillColorAlpha(color[iQValue], colbtOpacity[iQValue]);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin]->Draw("3,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin], Form("CoLBT, q=%.1f", myQValue), "f");
            }
          } else if(theoryComparisonIndex == 3){
            // Theory comparison index 3, comparison between Hybrid model, Holguin calculation and CoLBT

            // First, set the style and draw the Hybrid model prediction with full wake
            eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, HybridModelHistogramManager::kFullWake);

            // There are some bins for which the prediction does not exist
            if(energyEnergyCorrelatorHybridModelPbPbToPpRatio[eecPbPbMCBin] != NULL){

              // Give some nice styles for the predictions
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[eecPbPbMCBin]->SetLineColor(color[2]);
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(color[2]);
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[eecPbPbMCBin]->SetFillColorAlpha(color[2], 0.4);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[eecPbPbMCBin]->Draw("3,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorHybridModelPbPbToPpRatio[eecPbPbMCBin], hybridHistograms->GetWakeName(HybridModelHistogramManager::kFullWake), "f");
              }

            // Next, add the perturbative calculation with k=0.3
            iKValue = holguinHistograms->FindKValueIndex(0.3);
            eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iKValue);

            // There are some bins for which the prediction does not exist
            if(energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin] != NULL){

              // Give some nice styles for the predictions
              energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin]->SetLineColor(color[1]);
              energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(color[1]);
              energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin]->Draw("pl,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin], "Holguin, k=0.3", "pl");
            }

            // Finally, add CoLBT prediction with q=1
            iQValue = coLBTHistograms->FindQValueIndex(1);
            eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iQValue);

            // There are some bins for which the prediction does not exist
            if(energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin] != NULL){

              // Give some nice styles for the predictions
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin]->SetLineColor(color[0]);
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(color[0]);
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin]->SetFillColorAlpha(color[0], 0.4);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin]->Draw("3,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin], "CoLBT, q=1", "f");
            }

          } else if (theoryComparisonIndex == 4){
            // Theory comparison index 4, comparison between best k in Holguin calculation and CoLBT

            // Add the perturbative calculation with k=0.3
            iKValue = holguinHistograms->FindKValueIndex(0.3);
            eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iKValue);

            // There are some bins for which the prediction does not exist
            if(energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin] != NULL){

              // Give some nice styles for the predictions
              energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin]->SetLineColor(kMagenta);
              energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(kMagenta);
              energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin]->Draw("pl,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin], "Holguin, k=0.3", "pl");
            } 

            // Then add all different q-values for CoLBT
            for(auto myQValue : coLBTqValue){
              iQValue = coLBTHistograms->FindQValueIndex(myQValue);
              eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iQValue);

              // There are some bins for which the prediction does not exist
              if(energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin] == NULL) continue;

              // Give some nice styles for the predictions
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin]->SetLineColor(color[iQValue]);
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(color[iQValue]);
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin]->SetFillColorAlpha(color[iQValue], 0.4);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin]->Draw("3,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorCoLBTPbPbToPpRatio[eecPbPbMCBin], Form("CoLBT, q=%.1f", myQValue), "f");
            }

          } else if (theoryComparisonIndex == 5){
            // Theory comparison index 5, comparison to JEWEL

            // Compare distributions with and without recoil
            for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){
              eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iRecoil);

              // There might be some bins for which the prediction does not exist
              if(energyEnergyCorrelatorJewelPbPbToPpRatio[eecPbPbMCBin] == NULL) continue;

              // Give some nice styles for the predictions
              energyEnergyCorrelatorJewelPbPbToPpRatio[eecPbPbMCBin]->SetLineColor(color[iRecoil]);
              energyEnergyCorrelatorJewelPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(color[iRecoil]);
              energyEnergyCorrelatorJewelPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(jewelMarkerStyle[iRecoil]);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorJewelPbPbToPpRatio[eecPbPbMCBin]->Draw("pl,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorJewelPbPbToPpRatio[eecPbPbMCBin], jewelHistograms->GetRecoilName(iRecoil), "pl");
            }
          } else if (theoryComparisonIndex == 6){
            // Theory comparison index 6, comparison to best k-value from Holguin and JEWEL

            // Add the perturbative calculation with k=0.3
            iKValue = holguinHistograms->FindKValueIndex(0.3);
            eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iKValue);

            // There are some bins for which the prediction does not exist
            if(energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin] != NULL){

              // Give some nice styles for the predictions
              energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin]->SetLineColor(kMagenta);
              energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(kMagenta);
              energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin]->Draw("pl,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorHolguinPbPbToPpRatio[eecPbPbMCBin], "Holguin, k=0.3", "pl");
            } 

            // Add the JEWEL predictions with and without recoil to the same plot
            for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){
              eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iRecoil);

              // There might be some bins for which the prediction does not exist
              if(energyEnergyCorrelatorJewelPbPbToPpRatio[eecPbPbMCBin] == NULL) continue;

              // Give some nice styles for the predictions
              energyEnergyCorrelatorJewelPbPbToPpRatio[eecPbPbMCBin]->SetLineColor(color[iRecoil]);
              energyEnergyCorrelatorJewelPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(color[iRecoil]);
              energyEnergyCorrelatorJewelPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(jewelMarkerStyle[iRecoil]);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorJewelPbPbToPpRatio[eecPbPbMCBin]->Draw("pl,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorJewelPbPbToPpRatio[eecPbPbMCBin], jewelHistograms->GetRecoilName(iRecoil), "pl");
            }
          } else if(theoryComparisonIndex == 7){
            // Theory comparison index 7: Full wake hybrid and JEWEL with recoils

            // Hybrid full wake
            eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, HybridModelHistogramManager::kFullWake);
            if(energyEnergyCorrelatorHybridModelPbPbToPpRatio[eecPbPbMCBin] != NULL){

              // Give some nice styles for the predictions
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[eecPbPbMCBin]->SetLineColor(color[HybridModelHistogramManager::kFullWake]);
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(color[HybridModelHistogramManager::kFullWake]);
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[eecPbPbMCBin]->SetFillColorAlpha(color[HybridModelHistogramManager::kFullWake], jetWakeOpacity[HybridModelHistogramManager::kFullWake]);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorHybridModelPbPbToPpRatio[eecPbPbMCBin]->Draw("3,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorHybridModelPbPbToPpRatio[eecPbPbMCBin], hybridHistograms->GetWakeName(HybridModelHistogramManager::kFullWake), "f");
            }

            // JEWEL with recoil
            eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, JewelHistogramManager::kRecoil);
            if(energyEnergyCorrelatorJewelPbPbToPpRatio[eecPbPbMCBin] != NULL){

              // Give some nice styles for the predictions
              energyEnergyCorrelatorJewelPbPbToPpRatio[eecPbPbMCBin]->SetLineColor(color[JewelHistogramManager::kRecoil]);
              energyEnergyCorrelatorJewelPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(color[JewelHistogramManager::kRecoil]);
              energyEnergyCorrelatorJewelPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(jewelMarkerStyle[JewelHistogramManager::kRecoil]);

              // Draw the prediction to the same canvas as the data
              energyEnergyCorrelatorJewelPbPbToPpRatio[eecPbPbMCBin]->Draw("pl,same");

              // Add a legend for the theory prediction
              anotherLegend->AddEntry(energyEnergyCorrelatorJewelPbPbToPpRatio[eecPbPbMCBin], jewelHistograms->GetRecoilName(JewelHistogramManager::kRecoil), "pl");
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
          hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBinSystematic]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
          hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBinSystematic]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

          // Adjust tick size for x-axis
          hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBinSystematic]->GetXaxis()->SetTickLength(tickSizeX*1.5);

          // Set the style for uncertainty bands for systematic and statistical uncertainties from data
          for(int iUncertainty = 0; iUncertainty < knRelativeUncertaintyTypes; iUncertainty++){
            eecPbPbRelativeUncertaintyBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iUncertainty);
            hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBin]->SetLineColor(kBlack);
            hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBin]->SetLineStyle(9);
            hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBin]->SetMarkerStyle(kFullCircle);
            hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBin]->SetMarkerSize(0);
          }
          hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBinSystematic]->SetFillColorAlpha(kBlack, 0.4);
          hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBinSystematic]->SetLineWidth(0);


          // Draw the error bars from the data
          drawer->DrawHistogramToLowerPad(hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBinSystematic], "#Deltar", "#frac{Theory}{Data}", " ", "e3");
          hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBinStatisticalUp]->Draw("same,HIST,C");
          hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBinStatisticalDown]->Draw("same,HIST,C");

          // Add the different uncertainties to the legend
          legend->AddEntry(hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBinSystematic], "Data syst. unc.", "f");
          legend->AddEntry(hRelativeUncertaintyPbPbToPpRatio[eecPbPbRelativeUncertaintyBinStatisticalUp], "Data stat. unc.", "l");

          // Choose a set of theory predictions to draw to the figure
          if(theoryComparisonIndex == 0){
            // Theory comparison index 0, Hybrid model with different wake settings
            for(int iWake = 0; iWake < HybridModelHistogramManager::kWakeConfigurations; iWake++){
              eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iWake);
              hybridModelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(color[iWake]);
              hybridModelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);
              hybridModelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerSize(0);
              hybridModelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetFillColorAlpha(color[iWake], jetWakeOpacity[iWake]);
              hybridModelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->Draw("same,e3");
            }
          } else if (theoryComparisonIndex == 1){
            // Theory comparison index 1, perturbative calculations from Holguin et. al. with different k-values
            styleIndex = 0;
            for(auto myKValue : holguinKValue){
              iKValue = holguinHistograms->FindKValueIndex(myKValue);
              eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iKValue);
              if(holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin] == NULL) continue;
              holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(color[iKValue]);
              holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetLineWidth(3);
              holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetLineColor(color[iKValue]);
              holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(markerStyleTheory[styleIndex++]);
              holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin]->Draw("same,HIST,C");
            }
          } else if(theoryComparisonIndex == 2){
            // Theory comparison index 2, CoLBT prediction
            for(auto myQValue : coLBTqValue){
              iQValue = coLBTHistograms->FindQValueIndex(myQValue);
              eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iQValue);
              if(coLBTToDataRatioPbPbToPpRatio[eecPbPbMCBin] == NULL) continue;
              coLBTToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(color[iQValue]);
              coLBTToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);
              coLBTToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerSize(0);
              coLBTToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetFillColorAlpha(color[iQValue], colbtOpacity[iQValue]);
              coLBTToDataRatioPbPbToPpRatio[eecPbPbMCBin]->Draw("same,e3");
            }
          } else if(theoryComparisonIndex == 3){
            // Theory comparison index 3, comparison between Hybrid model, Holguin calculation and CoLBT

            // First, add the Hybrid model prediction
            eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt,HybridModelHistogramManager::kFullWake);
            if(hybridModelToDataRatioPbPbToPpRatio[eecPbPbMCBin] != NULL){
              hybridModelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(color[2]);
              hybridModelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);
              hybridModelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerSize(0);
              hybridModelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetFillColorAlpha(color[2], 0.4);
              hybridModelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->Draw("same,e3");
            }

            // Then, add perturbative calculation by Holguin and friends
            iKValue = holguinHistograms->FindKValueIndex(0.3);
            eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iKValue);
            if(holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin] != NULL){
              holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(color[1]);
              holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetLineWidth(3);
              holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetLineColor(color[1]);
              holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);
              holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin]->Draw("same,HIST,C");
            }

            // Finally, add the CoLBT prediction
            iQValue = coLBTHistograms->FindQValueIndex(1);
            eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iQValue);
            if(coLBTToDataRatioPbPbToPpRatio[eecPbPbMCBin] != NULL){
              coLBTToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(color[0]);
              coLBTToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);
              coLBTToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerSize(0);
              coLBTToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetFillColorAlpha(color[0], 0.4);
              coLBTToDataRatioPbPbToPpRatio[eecPbPbMCBin]->Draw("same,e3");
            }

          } else if(theoryComparisonIndex == 4){
            // Theory comparison index 4, comparison between best k in Holguin calculation and CoLBT

            // Add the perturbative calculation by Holguin and friends
            iKValue = holguinHistograms->FindKValueIndex(0.3);
            eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iKValue);
            if(holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin] != NULL){
              holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(kMagenta);
              holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetLineWidth(3);
              holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetLineColor(kMagenta);
              holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);
              holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin]->Draw("same,HIST,C");
            }

            // Then, add CoLBT prediction
            for(auto myQValue : coLBTqValue){
              iQValue = coLBTHistograms->FindQValueIndex(myQValue);
              eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iQValue);
              if(coLBTToDataRatioPbPbToPpRatio[eecPbPbMCBin] == NULL) continue;
              coLBTToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(color[iQValue]);
              coLBTToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);
              coLBTToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerSize(0);
              coLBTToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetFillColorAlpha(color[iQValue], 0.4);
              coLBTToDataRatioPbPbToPpRatio[eecPbPbMCBin]->Draw("same,e3");
            }

          } else if(theoryComparisonIndex == 5){
            // Theory comparison index 5, JEWEL prediction

            for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){
              eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iRecoil);
              if(jewelToDataRatioPbPbToPpRatio[eecPbPbMCBin] == NULL) continue;
              jewelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(color[iRecoil]);
              jewelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(jewelMarkerStyle[iRecoil]);
              jewelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->Draw("same,pl");
            }
          } else if(theoryComparisonIndex == 6){
            // Theory comparison index 6, best k-value from Holguin and JEWEL prediction

            // Add perturbative calculation by Holguin and friends
            iKValue = holguinHistograms->FindKValueIndex(0.3);
            eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iKValue);
            if(holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin] != NULL){
              holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(kMagenta);
              holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetLineWidth(3);
              holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetLineColor(kMagenta);
              holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);
              holguinToDataRatioPbPbToPpRatio[eecPbPbMCBin]->Draw("same,HIST,C");
            }

            for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){
              eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, iRecoil);
              if(jewelToDataRatioPbPbToPpRatio[eecPbPbMCBin] == NULL) continue;
              jewelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(color[iRecoil]);
              jewelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(jewelMarkerStyle[iRecoil]);
              jewelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->Draw("same,pl");
            }
          } else if(theoryComparisonIndex == 7){
            // Theory comparison index 7, Hybrid full wake and JEWEL with recoils

            // Hybrid full wake
            eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, HybridModelHistogramManager::kFullWake);
            if(hybridModelToDataRatioPbPbToPpRatio[eecPbPbMCBin] != NULL){
              hybridModelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(color[HybridModelHistogramManager::kFullWake]);
              hybridModelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(kFullCircle);
              hybridModelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerSize(0);
              hybridModelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetFillColorAlpha(color[HybridModelHistogramManager::kFullWake], jetWakeOpacity[HybridModelHistogramManager::kFullWake]);
              hybridModelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->Draw("same,e3");
            }

            // JEWEL with recoils
            eecPbPbMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iTrackPt, JewelHistogramManager::kRecoil);
            if(jewelToDataRatioPbPbToPpRatio[eecPbPbMCBin] != NULL){
              jewelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerColor(color[JewelHistogramManager::kRecoil]);
              jewelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->SetMarkerStyle(jewelMarkerStyle[JewelHistogramManager::kRecoil]);
              jewelToDataRatioPbPbToPpRatio[eecPbPbMCBin]->Draw("same,pl");
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

        eecDoubleRatioBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt);
        eecDoubleRatioRelativeUncertaintyBinSystematic = std::make_tuple(weightExponent-1, iCentrality, iJetPt, kRelativeUncertaintySystematic);
        eecDoubleRatioRelativeUncertaintyBinStatisticalUp = std::make_tuple(weightExponent-1, iCentrality, iJetPt, kRelativeUncertaintyStatisticalUp);
        eecDoubleRatioRelativeUncertaintyBinStatisticalDown = std::make_tuple(weightExponent-1, iCentrality, iJetPt, kRelativeUncertaintyStatisticalDown);

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
        energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin]->SetMarkerStyle(kFullSquare);
        energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin]->SetMarkerColor(kBlack);
        energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin]->SetLineColor(kBlack);

        // Set drawing style for systematic uncertainties
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->SetFillColorAlpha(kBlack, 0.4);
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->SetLineColor(kBlack);
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->SetMarkerColor(kBlack);
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->SetMarkerStyle(kFullSquare);

        // Set the x-axis drawing range
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetYaxis()->SetRangeUser(doubleRatioZoom.first, doubleRatioZoom.second);

        // Adjust tick size for x-axis
        systematicUncertaintyDoubleRatio[eecDoubleRatioBin]->GetXaxis()->SetTickLength(tickSizeX);

        // Draw first the systematic uncertainties to the upper canves
        drawer->DrawHistogramToUpperPad(systematicUncertaintyDoubleRatio[eecDoubleRatioBin], "#Deltar", Form("#frac{PbPb/pp (p_{T}^{ch} > %.0f GeV)}{PbPb/pp (p_{T}^{ch} > %.0f GeV)}", trackPtCutsForDoubleDatio.second, trackPtCutsForDoubleDatio.first), " ", "e2");

        // Then draw the double ratio and add a legend for it
        energyEnergyCorrelatorDoubleRatio[eecDoubleRatioBin]->Draw("same,p");
        anotherLegend->AddEntry(systematicUncertaintyDoubleRatio[eecDoubleRatioBin], "Data", "lpf");

        if(theoryComparisonIndex == 0){
          // Theory comparison index 0: Hybrid model

          for(int iWake = 0; iWake < HybridModelHistogramManager::kWakeConfigurations; iWake++){

            eecDoubleRatioMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iWake);

            // There are some bins for which the prediction does not exist
            if(histogrammifiedHybridModelDoubleRatio[eecDoubleRatioMCBin] == NULL) continue;

            // Give some nice styles for the predictions
            histogrammifiedHybridModelDoubleRatio[eecDoubleRatioMCBin]->SetLineColor(color[iWake]);
            histogrammifiedHybridModelDoubleRatio[eecDoubleRatioMCBin]->SetMarkerStyle(kFullCircle);
            histogrammifiedHybridModelDoubleRatio[eecDoubleRatioMCBin]->SetMarkerSize(0);
            histogrammifiedHybridModelDoubleRatio[eecDoubleRatioMCBin]->SetFillColorAlpha(color[iWake], jetWakeOpacity[iWake]);

            // Draw the prediction to the same canvas as the data
            histogrammifiedHybridModelDoubleRatio[eecDoubleRatioMCBin]->Draw("same,e3");

            // Add a legend for the theory prediction
            anotherLegend->AddEntry(histogrammifiedHybridModelDoubleRatio[eecDoubleRatioMCBin], hybridHistograms->GetWakeName(iWake), "f");
          }
        } else if (theoryComparisonIndex == 5 || theoryComparisonIndex == 6){
          // Theory comparison index 5 or 6: JEWEL

          for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){
            eecDoubleRatioMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iRecoil);

            // There are some bins for which the prediction does not exist
            if(jewelDoubleRatio[eecDoubleRatioMCBin] == NULL) continue;

            // Give some nice styles for the predictions
            jewelDoubleRatio[eecDoubleRatioMCBin]->SetMarkerColor(color[iRecoil]);
            jewelDoubleRatio[eecDoubleRatioMCBin]->SetLineColor(color[iRecoil]);
            jewelDoubleRatio[eecDoubleRatioMCBin]->SetMarkerStyle(jewelMarkerStyle[iRecoil]);

            // Draw the prediction to the same canvas as the data
            jewelDoubleRatio[eecDoubleRatioMCBin]->Draw("same,pl");

            // Add a legend for the theory prediction
            anotherLegend->AddEntry(jewelDoubleRatio[eecDoubleRatioMCBin], jewelHistograms->GetRecoilName(iRecoil), "pl");
          }
        } else if (theoryComparisonIndex == 7){
          // Theory comparison index 7: Hybrid full wake and JEWEL with recoils

          // Hybrid full wake
          eecDoubleRatioMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, HybridModelHistogramManager::kFullWake);
          if(histogrammifiedHybridModelDoubleRatio[eecDoubleRatioMCBin] != NULL){

            // Give some nice styles for the predictions
            histogrammifiedHybridModelDoubleRatio[eecDoubleRatioMCBin]->SetLineColor(color[HybridModelHistogramManager::kFullWake]);
            histogrammifiedHybridModelDoubleRatio[eecDoubleRatioMCBin]->SetMarkerStyle(kFullCircle);
            histogrammifiedHybridModelDoubleRatio[eecDoubleRatioMCBin]->SetMarkerSize(0);
            histogrammifiedHybridModelDoubleRatio[eecDoubleRatioMCBin]->SetFillColorAlpha(color[HybridModelHistogramManager::kFullWake], jetWakeOpacity[HybridModelHistogramManager::kFullWake]);

            // Draw the prediction to the same canvas as the data
            histogrammifiedHybridModelDoubleRatio[eecDoubleRatioMCBin]->Draw("same,e3");

            // Add a legend for the theory prediction
            anotherLegend->AddEntry(histogrammifiedHybridModelDoubleRatio[eecDoubleRatioMCBin], hybridHistograms->GetWakeName(HybridModelHistogramManager::kFullWake), "f");
          }

          // JEWEL with recoils
          eecDoubleRatioMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, JewelHistogramManager::kRecoil);
          if(jewelDoubleRatio[eecDoubleRatioMCBin] != NULL){

            // Give some nice styles for the predictions
            jewelDoubleRatio[eecDoubleRatioMCBin]->SetMarkerColor(color[JewelHistogramManager::kRecoil]);
            jewelDoubleRatio[eecDoubleRatioMCBin]->SetLineColor(color[JewelHistogramManager::kRecoil]);
            jewelDoubleRatio[eecDoubleRatioMCBin]->SetMarkerStyle(jewelMarkerStyle[JewelHistogramManager::kRecoil]);

            // Draw the prediction to the same canvas as the data
            jewelDoubleRatio[eecDoubleRatioMCBin]->Draw("same,pl");

            // Add a legend for the theory prediction
            anotherLegend->AddEntry(jewelDoubleRatio[eecDoubleRatioMCBin], jewelHistograms->GetRecoilName(JewelHistogramManager::kRecoil), "pl");
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
        hRelativeUncertaintyDoubleRatio[eecDoubleRatioRelativeUncertaintyBinSystematic]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
        hRelativeUncertaintyDoubleRatio[eecDoubleRatioRelativeUncertaintyBinSystematic]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

        // Adjust tick size for x-axis
        hRelativeUncertaintyDoubleRatio[eecDoubleRatioRelativeUncertaintyBinSystematic]->GetXaxis()->SetTickLength(tickSizeX*1.5);

        // Set the style for uncertainty bands for systematic and statistical uncertainties from data
        for(int iUncertainty = 0; iUncertainty < knRelativeUncertaintyTypes; iUncertainty++){
          eecDoubleRatioRelativeUncertaintyBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iUncertainty);
          hRelativeUncertaintyDoubleRatio[eecDoubleRatioRelativeUncertaintyBin]->SetLineColor(kBlack);
          hRelativeUncertaintyDoubleRatio[eecDoubleRatioRelativeUncertaintyBin]->SetLineStyle(9);
          hRelativeUncertaintyDoubleRatio[eecDoubleRatioRelativeUncertaintyBin]->SetMarkerStyle(kFullCircle);
          hRelativeUncertaintyDoubleRatio[eecDoubleRatioRelativeUncertaintyBin]->SetMarkerSize(0);
        }
        hRelativeUncertaintyDoubleRatio[eecDoubleRatioRelativeUncertaintyBinSystematic]->SetFillColorAlpha(kBlack, 0.4);
        hRelativeUncertaintyDoubleRatio[eecDoubleRatioRelativeUncertaintyBinSystematic]->SetLineWidth(0);


        // Draw the error bars from the data
        drawer->DrawHistogramToLowerPad(hRelativeUncertaintyDoubleRatio[eecDoubleRatioRelativeUncertaintyBinSystematic], "#Deltar", "#frac{Theory}{Data}", " ", "e3");
        hRelativeUncertaintyDoubleRatio[eecDoubleRatioRelativeUncertaintyBinStatisticalUp]->Draw("same,HIST,C");
        hRelativeUncertaintyDoubleRatio[eecDoubleRatioRelativeUncertaintyBinStatisticalDown]->Draw("same,HIST,C");

        // Add the different uncertainties to the legend
        legend->AddEntry(hRelativeUncertaintyDoubleRatio[eecDoubleRatioRelativeUncertaintyBinSystematic], "Data syst. unc.", "f");
        legend->AddEntry(hRelativeUncertaintyDoubleRatio[eecDoubleRatioRelativeUncertaintyBinStatisticalUp], "Data stat. unc.", "l");

        if(theoryComparisonIndex == 0){
          // Theory comparison index 0: Hybrid model

          for(int iWake = 0; iWake < HybridModelHistogramManager::kWakeConfigurations; iWake++){
            eecDoubleRatioMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iWake);
            hybridModelToDataRatioDoubleRatio[eecDoubleRatioMCBin]->SetMarkerStyle(kFullCircle);
            hybridModelToDataRatioDoubleRatio[eecDoubleRatioMCBin]->SetMarkerSize(0);
            hybridModelToDataRatioDoubleRatio[eecDoubleRatioMCBin]->SetFillColorAlpha(color[iWake], jetWakeOpacity[iWake]);
            hybridModelToDataRatioDoubleRatio[eecDoubleRatioMCBin]->Draw("same,e3");
          }
        } else if(theoryComparisonIndex == 5 || theoryComparisonIndex == 6){
          // Theory comparison index 5 or 6: JEWEL

          for(int iRecoil = 0; iRecoil < JewelHistogramManager::kRecoilSettings; iRecoil++){
            eecDoubleRatioMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, iRecoil);

            // There are some bins for which the prediction does not exist
            if(jewelToDataRatioDoubleRatio[eecDoubleRatioMCBin] == NULL) continue;

            jewelToDataRatioDoubleRatio[eecDoubleRatioMCBin]->SetMarkerStyle(jewelMarkerStyle[iRecoil]);
            jewelToDataRatioDoubleRatio[eecDoubleRatioMCBin]->SetMarkerColor(color[iRecoil]);
            jewelToDataRatioDoubleRatio[eecDoubleRatioMCBin]->SetLineColor(color[iRecoil]);
            jewelToDataRatioDoubleRatio[eecDoubleRatioMCBin]->Draw("same,pl");
          }
        } else if(theoryComparisonIndex == 7){
          // Theory comparison index 7: Hybrid full wake and JEWEL with recoils

          // Hybrid full wake
          eecDoubleRatioMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, HybridModelHistogramManager::kFullWake);
          if(hybridModelToDataRatioDoubleRatio[eecDoubleRatioMCBin] != NULL){
            hybridModelToDataRatioDoubleRatio[eecDoubleRatioMCBin]->SetMarkerStyle(kFullCircle);
            hybridModelToDataRatioDoubleRatio[eecDoubleRatioMCBin]->SetMarkerSize(0);
            hybridModelToDataRatioDoubleRatio[eecDoubleRatioMCBin]->SetFillColorAlpha(color[HybridModelHistogramManager::kFullWake], jetWakeOpacity[HybridModelHistogramManager::kFullWake]);
            hybridModelToDataRatioDoubleRatio[eecDoubleRatioMCBin]->Draw("same,e3");
          }

          // JEWEL with recoils
          eecDoubleRatioMCBin = std::make_tuple(weightExponent-1, iCentrality, iJetPt, JewelHistogramManager::kRecoil);
          if(jewelToDataRatioDoubleRatio[eecDoubleRatioMCBin] != NULL){

            jewelToDataRatioDoubleRatio[eecDoubleRatioMCBin]->SetMarkerStyle(jewelMarkerStyle[JewelHistogramManager::kRecoil]);
            jewelToDataRatioDoubleRatio[eecDoubleRatioMCBin]->SetMarkerColor(color[JewelHistogramManager::kRecoil]);
            jewelToDataRatioDoubleRatio[eecDoubleRatioMCBin]->SetLineColor(color[JewelHistogramManager::kRecoil]);
            jewelToDataRatioDoubleRatio[eecDoubleRatioMCBin]->Draw("same,pl");
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
