#include "SystematicUncertaintyOrganizer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "EECHistogramManager.h"
#include "JDrawer.h"

#include <algorithm>

const int kMaxPtBins = 6; // Maximum number of pT bins that is used to extract the pT variation uncertainty

// Function definitions. Implementations after the main macro.
TH1D* findTheDifference(TH1D* nominalResult, TH1D* variedResult[], const int nComparisonGraphs);
TH1D* findTheDifference(TH1D* nominalResult, TH1D* variedResult);
void drawIllustratingPlots(JDrawer* drawer, TH1D* nominalResult, TH1D* variedResult[], const int nVariations, const int iCentrality, const int iJetPt, const int iTrackPt, EECCard* card, TString comparisonLegend[], TString plotName, std::pair<double, double> drawingRange, std::pair<double, double> ratioZoom);
void drawIllustratingPlots(JDrawer* drawer, TH1D* nominalResult, TH1D* variedResult, const int iCentrality, const int iJetPt, const int iTrackPt, EECCard* card, TString comparisonLegend, TString plotName, std::pair<double, double> drawingRange, std::pair<double, double> ratioZoom);
void loadRelevantHistograms(EECHistogramManager* histograms);

/*
 * Macro for estimating systematic uncertainties for energy-energy correlators
 *
 * The following sources of systematic uncertainty are estimated:
 *  - Jet energy resolition uncertainty (modifies unfolding response matrix, correlated)
 *  - Jet energy scale uncertainty (modifies unfolding response matrix, correlated)
 *  - Truth reference for jet pT unfolding (modifies unfolding response matrix, correlated)
 *  - Tracking efficiency (uncorrelated)
 *  - Track pair efficiency (uncorrelated)
 *  - Background subtraction (uncorrelated)
 *
 *  The results will be saved to a file and also slides with different contributions can be printed
 */
void estimateSystematicUncertainties(){
  
  // ==================================================================
  // ============================= Input ==============================
  // ==================================================================

  // Nominal results
  TFile* nominalResultFile = TFile::Open("data/eecAnalysis_akFlowJet_wtaAxis_binningForUnfolding_processed_nominalResult_2023-05-23.root");
  EECCard* nominalResultCard = new EECCard(nominalResultFile);
  EECHistogramManager* nominalHistogramManager = new EECHistogramManager(nominalResultFile, nominalResultCard);
  loadRelevantHistograms(nominalHistogramManager);
  
  // Results unfolded with a response matrix smeared with jet energy resolution
  TFile* jetEnergyResolutionFile = TFile::Open("data/eecAnalysis_akFlowJet_wtaAxis_binningForUnfolding_processed_jetEnergyResolutionSystematicUncertainty_2023-05-23.root");
  EECCard* jetEnergyResolutionCard = new EECCard(jetEnergyResolutionFile);
  EECHistogramManager* jetEnergyResolutionHistogramManager = new EECHistogramManager(jetEnergyResolutionFile, jetEnergyResolutionCard);
  loadRelevantHistograms(jetEnergyResolutionHistogramManager);
  
  // Results unfolded with a response matrix smeared with jet energy scale
  TFile* jetEnergyScaleFile = TFile::Open("data/eecAnalysis_akFlowJet_wtaAxis_binningForUnfolding_processed_jetEnergyScaleSystematicUncertainty_2023-05-23.root");
  EECCard* jetEnergyScaleCard = new EECCard(jetEnergyScaleFile);
  EECHistogramManager* jetEnergyScaleHistogramManager = new EECHistogramManager(jetEnergyScaleFile, jetEnergyScaleCard);
  loadRelevantHistograms(jetEnergyScaleHistogramManager);

  // Results where background scaling factor is determined from 2% or 6% shifted simulation instead of nominal 4%
  TFile* backgroundSubtractionFile[2];
  backgroundSubtractionFile[0] = TFile::Open("data/eecAnalysis_akFlowJet_wtaAxis_binningForUnfolding_processed_backgroundScaleUncertainty2pShift_2023-05-23.root");
  backgroundSubtractionFile[1] = TFile::Open("data/eecAnalysis_akFlowJet_wtaAxis_binningForUnfolding_processed_backgroundScaleUncertainty6pShift_2023-05-23.root");
  EECCard* backgroundSubtractionCard[2];
  EECHistogramManager* backgroundSubtractionHistogramManager[2];
  for(int iBackgroundSubtractionFile = 0; iBackgroundSubtractionFile < 2; iBackgroundSubtractionFile++){
    backgroundSubtractionCard[iBackgroundSubtractionFile] = new EECCard(backgroundSubtractionFile[iBackgroundSubtractionFile]);
    backgroundSubtractionHistogramManager[iBackgroundSubtractionFile] = new EECHistogramManager(backgroundSubtractionFile[iBackgroundSubtractionFile], backgroundSubtractionCard[iBackgroundSubtractionFile]);
    loadRelevantHistograms(backgroundSubtractionHistogramManager[iBackgroundSubtractionFile]);
  }
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = nominalResultCard->GetNCentralityBins();
  const int nJetPtBinsEEC = nominalResultCard->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = nominalResultCard->GetNTrackPtBinsEEC();
  
  // Estimate the uncertainties for all bins that have been unfolded
  int firstStudiedCentralityBin = nominalResultCard->GetFirstUnfoldedCentralityBin();
  int lastStudiedCentralityBin = nominalResultCard->GetLastUnfoldedCentralityBin();
  
  int firstStudiedJetPtBinEEC = nominalResultCard->GetFirstUnfoldedJetPtBin();
  int lastStudiedJetPtBinEEC = nominalResultCard->GetLastUnfoldedJetPtBin();
  
  int firstStudiedTrackPtBinEEC = nominalResultCard->GetFirstUnfoldedTrackPtBin();
  int lastStudiedTrackPtBinEEC = nominalResultCard->GetLastUnfoldedTrackPtBin();

  // Only draw example plots from selected subset of bins
  vector<int> drawnCentralityBins = {0};
  vector<int> drawnJetPtBins = {6,7,8,9};
  vector<int> drawnTrackPtBins = {5};
  
  const bool printUncertainties = false;
  
  bool plotExample = true;                                                // Flag for drawing illustrative plots
  std::pair<double, double> analysisDeltaR = std::make_pair(0.006, 0.39); // DeltaR span in which tha analysis is done
  std::pair<double, double> ratioZoom = std::make_pair(0.6, 1.4);         // Y-axis zoom for rations
  
  TString outputFileName = "systematicUncertainties/systematicUncertainties_firstLook_2023-06-06.root";
  
  // Option to skip evaluating some of the sources defined in SystematicUncertaintyOrganizer
  bool skipUncertaintySource[SystematicUncertaintyOrganizer::knUncertaintySources];
  for(int iUncertainty = 0; iUncertainty < SystematicUncertaintyOrganizer::knUncertaintySources; iUncertainty++){
    skipUncertaintySource[iUncertainty] = false;
  }
  
  // ==================================================================
  // ====================== Configuration done ========================
  // ==================================================================
  
  // Systematic uncertainty organizer can easily provide names for all your naming purposes
  SystematicUncertaintyOrganizer* nameGiver = new SystematicUncertaintyOrganizer();
  
  // Histograms for uncertainty estimation
  TH1D* nominalEnergyEnergyCorrelators[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* jetEnergyResolutionUncertaintyCorrelators[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* jetEnergyScaleUncertaintyCorrelators[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* backgroundSubtractionUncertaintyCorrelators[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][2];

  // Histograms to hold the systematic uncertainty results
  TH1D* energyEnergyCorrelatorSystematicUncertainties[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][SystematicUncertaintyOrganizer::knUncertaintySources];
  
  // Initialize all the histograms to NULL
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt] = NULL;
        jetEnergyResolutionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt] = NULL;
        jetEnergyScaleUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt] = NULL;
        for(int iFile = 0; iFile < 2; iFile++){
          backgroundSubtractionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile] = NULL;
        }
        for(int iUncertainty = 0; iUncertainty < SystematicUncertaintyOrganizer::knUncertaintySources; iUncertainty++){
          energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][iUncertainty] = NULL;
        }
      }
    }
  }

  // Read the histograms from the defined range from the file
  int iCentralityMatched;
  int iTrackPtMatched;
  int iJetPtMatched;
  double epsilon = 0.0001;
  int lowAnalysisBin, highAnalysisBin;
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

        // Read the nominal histograms
        nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt] = nominalHistogramManager->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

        // Normalize nominal histograms to one in the drawing range
        lowAnalysisBin = nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(analysisDeltaR.first + epsilon);
        highAnalysisBin = nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(analysisDeltaR.second - epsilon);
        nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt]->Scale(1.0 / nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

        // Read the jet energy resolution histograms
        iCentralityMatched = jetEnergyResolutionCard->FindBinIndexCentrality(nominalResultCard->GetBinBordersCentrality(iCentrality));
        iTrackPtMatched = jetEnergyResolutionCard->FindBinIndexTrackPtEEC(nominalResultCard->GetBinBordersTrackPtEEC(iTrackPt));
        iJetPtMatched = jetEnergyResolutionCard->FindBinIndexJetPtEEC(nominalResultCard->GetBinBordersJetPtEEC(iJetPt));
        jetEnergyResolutionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt] = jetEnergyResolutionHistogramManager->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentralityMatched, iJetPtMatched, iTrackPtMatched, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

        // Normalize the jet energy resolution histograms to one
        jetEnergyResolutionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt]->Scale(1.0 / jetEnergyResolutionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

        // Read the jet energy scale histograms
        iCentralityMatched = jetEnergyScaleCard->FindBinIndexCentrality(nominalResultCard->GetBinBordersCentrality(iCentrality));
        iTrackPtMatched = jetEnergyScaleCard->FindBinIndexTrackPtEEC(nominalResultCard->GetBinBordersTrackPtEEC(iTrackPt));
        iJetPtMatched = jetEnergyScaleCard->FindBinIndexJetPtEEC(nominalResultCard->GetBinBordersJetPtEEC(iJetPt));
        jetEnergyScaleUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt] = jetEnergyScaleHistogramManager->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentralityMatched, iJetPtMatched, iTrackPtMatched, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

        // Normalize the jet energy scale histograms to one
        jetEnergyScaleUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt]->Scale(1.0 / jetEnergyScaleUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

        for(int iFile = 0; iFile < 2; iFile++){

          // Read the background subtraction uncertainty histograms
          iCentralityMatched = backgroundSubtractionCard[iFile]->FindBinIndexCentrality(nominalResultCard->GetBinBordersCentrality(iCentrality));
          iTrackPtMatched = backgroundSubtractionCard[iFile]->FindBinIndexTrackPtEEC(nominalResultCard->GetBinBordersTrackPtEEC(iTrackPt));
          iJetPtMatched = backgroundSubtractionCard[iFile]->FindBinIndexJetPtEEC(nominalResultCard->GetBinBordersJetPtEEC(iJetPt));
          backgroundSubtractionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile] = backgroundSubtractionHistogramManager[iFile]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentralityMatched, iJetPtMatched, iTrackPtMatched, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

          // Normalize the background subtraction histograms to one
          backgroundSubtractionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile]->Scale(1.0 / backgroundSubtractionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt][iFile]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

        }

      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
  
  // Drawer for illustrative plots
  JDrawer* drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  drawer->SetLogX(true);

  // Legends given to drawns graphs
  TString legendNames[2];
  
  // ================================================= //
  //   Uncertainty coming from jet energy resolution   //
  // ================================================= //
  
  if(!skipUncertaintySource[SystematicUncertaintyOrganizer::kJetEnergyResolution]){
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

          energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kJetEnergyResolution] = findTheDifference(nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt], jetEnergyResolutionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt]);

          // Draw example plots on how the uncertainty is obtained
          if(plotExample && std::binary_search(drawnCentralityBins.begin(), drawnCentralityBins.end(), iCentrality) && std::binary_search(drawnJetPtBins.begin(), drawnJetPtBins.end(), iJetPt) && std::binary_search(drawnTrackPtBins.begin(), drawnTrackPtBins.end(), iTrackPt)){
            legendNames[0] = "JER varied response";
            drawIllustratingPlots(drawer, nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt], jetEnergyResolutionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt], iCentrality, iJetPt, iTrackPt, nominalResultCard, legendNames[0], nameGiver->GetSystematicUncertaintyName(SystematicUncertaintyOrganizer::kJetEnergyResolution), analysisDeltaR, ratioZoom);
          }
      
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Jet energy resolution uncertainty
  
  // ============================================ //
  //   Uncertainty coming from jet energy scale   //
  // ============================================ //

  if(!skipUncertaintySource[SystematicUncertaintyOrganizer::kJetEnergyScale]){
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

          energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kJetEnergyScale] = findTheDifference(nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt], jetEnergyScaleUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt]);

          // Draw example plots on how the uncertainty is obtained
          if(plotExample && std::binary_search(drawnCentralityBins.begin(), drawnCentralityBins.end(), iCentrality) && std::binary_search(drawnJetPtBins.begin(), drawnJetPtBins.end(), iJetPt) && std::binary_search(drawnTrackPtBins.begin(), drawnTrackPtBins.end(), iTrackPt)){
            legendNames[0] = "JEC varied response";
            drawIllustratingPlots(drawer, nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt], jetEnergyScaleUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt], iCentrality, iJetPt, iTrackPt, nominalResultCard, legendNames[0], nameGiver->GetSystematicUncertaintyName(SystematicUncertaintyOrganizer::kJetEnergyScale), analysisDeltaR, ratioZoom);
          }
      
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Jet energy scale uncertainty
  
  // ================================================== //
  //   Uncertainty coming from background subtraction   //
  // ================================================== //
  
  if(!skipUncertaintySource[SystematicUncertaintyOrganizer::kBackgroundSubtraction]){
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

          energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kBackgroundSubtraction] = findTheDifference(nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt], backgroundSubtractionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt], 2);

          // Draw example plots on how the uncertainty is obtained
          if(plotExample && std::binary_search(drawnCentralityBins.begin(), drawnCentralityBins.end(), iCentrality) && std::binary_search(drawnJetPtBins.begin(), drawnJetPtBins.end(), iJetPt) && std::binary_search(drawnTrackPtBins.begin(), drawnTrackPtBins.end(), iTrackPt)){
            legendNames[0] = "Bg scale from 2% shift";
            legendNames[1] = "Bg scale from 6% shift";
            drawIllustratingPlots(drawer, nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt], backgroundSubtractionUncertaintyCorrelators[iCentrality][iJetPt][iTrackPt], 2, iCentrality, iJetPt, iTrackPt, nominalResultCard, legendNames, nameGiver->GetSystematicUncertaintyName(SystematicUncertaintyOrganizer::kBackgroundSubtraction), analysisDeltaR, ratioZoom);
          }
      
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Jet energy resolution uncertainty
  
  // ============================================= //
  //   Add all uncertainty sources in quadrature   //
  // ============================================= //
  
  double sumOfSquares;
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

        // Clone the uncertainty histogram from the nominal one
        energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kAll] = (TH1D*) nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt]->Clone(Form("summedUncertainty%d%d%d", iCentrality, iJetPt, iTrackPt));
      
        // Calculate a sum of squared from all evaluated uncertainties
        for(int iBin = 1; iBin <= nominalEnergyEnergyCorrelators[iCentrality][iJetPt][iTrackPt]->GetNbinsX(); iBin++){
          sumOfSquares = 0;
          for(int iUncertainty = 0; iUncertainty < SystematicUncertaintyOrganizer::kAll; iUncertainty++){
            if(energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][iUncertainty] != NULL){
              sumOfSquares = sumOfSquares + TMath::Power(energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][iUncertainty]->GetBinError(iBin),2);
            }
          } // Uncertainty type loop
          energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][SystematicUncertaintyOrganizer::kAll]->SetBinError(iBin, TMath::Sqrt(sumOfSquares));
        } // Bin loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
  
  // =========================================================== //
  //   Write the systematic uncertainty histograms into a file   //
  // =========================================================== //
  
  // If an output file name is given, put the total uncertainties to an output file
  if(outputFileName.EndsWith(".root")){
    TFile* outputFile = new TFile(outputFileName,"UPDATE");
    TString saveName;

    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
          for(int iUncertainty = 0; iUncertainty < SystematicUncertaintyOrganizer::knUncertaintySources; iUncertainty++){
            if(energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][iUncertainty] != NULL){
              saveName = Form("systematicUncertainty_%s_C%dJ%dT%d", nameGiver->GetSystematicUncertaintyName(iUncertainty).Data(), iCentrality, iJetPt, iTrackPt);
              energyEnergyCorrelatorSystematicUncertainties[iCentrality][iJetPt][iTrackPt][iUncertainty]->Write(saveName);
            }
          } // Uncertainty loop
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop

    outputFile->Close();
  }
  
  // =========================================================================== //
  // Print slides summarizing all the different sources of systemtic uncertainty //
  // =========================================================================== //
  
  // if(printUncertainties){
    
  //   // Print a slide with uncertainties for each source and each centrality for each V
  //   char namer[100];
  //   for(int iFlow = firstFlow-1; iFlow <= lastFlow-1; iFlow++){
      
  //     sprintf(namer,"\\frametitle{Uncertainty for jet $v_{%d}$}", iFlow+1);
      
  //     cout << endl;
  //     cout << "\\begin{frame}" << endl;
  //     cout << namer << endl;
  //     cout << "\\begin{center}" << endl;
  //     cout << "  \\begin{tabular}{cccccc}" << endl;
  //     cout << "    \\toprule" << endl;
  //     cout << "    Source & C: 0-10 & C: 10-30 & C: 30-50 \\\\" << endl;
  //     cout << "    \\midrule" << endl;
      
  //     // Set the correct precision for printing floating point numbers
  //     cout << fixed << setprecision(4);
      
  //     // First, print the actual values of vn:s
  //     cout << "$v_{" << iFlow+1 << "}$ ";
      
  //     for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
  //       nominalResultGraph[iFlow]->GetPoint(iCentrality, nominalX, nominalY);
  //       cout << " & $" << nominalY << "$";
  //     }
  //     cout << " \\\\" << endl;
  //     cout << "    \\midrule" << endl;
      
  //     // Print the relative uncertainties
      
  //     for(int iUncertainty = 0; iUncertainty < SystematicUncertaintyOrganizer::knUncertaintySources; iUncertainty++){
  //       cout << uncertaintyNamer->GetLongRangeUncertaintyName(iUncertainty).Data();
  //       for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
  //         cout << " & $" << relativeUncertaintyTable[iUncertainty][iFlow][iCentrality] << "$";
  //       }
  //       cout << " \\\\" << endl;
  //     }
      
  //     cout << "    \\midrule" << endl;
      
  //     // Print the absolute uncertainties
      
  //     for(int iUncertainty = 0; iUncertainty < SystematicUncertaintyOrganizer::knUncertaintySources; iUncertainty++){
  //       cout << uncertaintyNamer->GetLongRangeUncertaintyName(iUncertainty).Data();
  //       for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
  //         cout << " & $" << absoluteUncertaintyTable[iUncertainty][iFlow][iCentrality] << "$";
  //       }
  //       cout << " \\\\" << endl;
  //     }
      
      
  //     cout << "    \\bottomrule" << endl;
  //     cout << "  \\end{tabular}" << endl;
  //     cout << "\\end{center}" << endl;
  //     cout << "\\begin{itemize}" << endl;
  //     cout << "  \\item Top: value. Middle: relative uncertainty. Bottom: absolute uncertainty." << endl;
  //     cout << "\\end{itemize}" << endl;
  //     cout << "\\end{frame}" << endl;
  //     cout << endl;
      
  //   } // vn loop
    
  // } // Printing uncertainties
}

/*
 * Function for finding the relative and absolute difference of points in two graphs
 *
 *  TH1D* nominalResult = Histogram containing nominal results
 *  TH1D* variedResult = Histograms containing varied results for systematic uncertainty estimation
 *  const int nVariations = Number of variations used to estimate the systematic uncertainties
 *
 *  return: Histogram where the error bars in each point correspond to estimated systematic uncertainty
 *
 */
TH1D* findTheDifference(TH1D *nominalResult, TH1D *variedResult[], const int nVariations){
  
  // Create an uncertainty histogram from the nominal results
  TH1D* uncertaintyHistogram = (TH1D*) nominalResult->Clone(Form("uncertaintyFor%s", variedResult[0]->GetName()));
  
  // Loop over the bins in the histogram, and in each bin assign the greates variation from the nominal result as uncertainty
  double biggestDifference, currentDifference;
  for(int iBin = 1; iBin <= nominalResult->GetNbinsX(); iBin++){
    biggestDifference = 0;
    for(int iVariation = 0; iVariation < nVariations; iVariation++){
      currentDifference = TMath::Abs(nominalResult->GetBinContent(iBin) - variedResult[iVariation]->GetBinContent(iBin));
      if(currentDifference > biggestDifference) biggestDifference = currentDifference;
    }
    uncertaintyHistogram->SetBinError(iBin, currentDifference);
  } // Histogram bin loop
  
  // Return the histogram where the errors represent the systematic uncertainties
  return uncertaintyHistogram;
  
}

/*
 * Function for finding the realtive and absolute difference of points in two graphs
 *
 *  TH1D* nominalResult = Histogram containing nominal results
 *  TH1D* variedResult = Histogram containing varied results for systematic uncertainty estimation
 *
 *  return: Histogram where the error bars in each point correspond to estimated systematic uncertainty
 *
 */
TH1D* findTheDifference(TH1D* nominalResult, TH1D* variedResult){
  
  // Enclose the single graph to an array and use the difference finder for graph array
  TH1D* comparisonArray[1] = {variedResult};
  return findTheDifference(nominalResult, comparisonArray, 1);
  
}

/*
 * Draw plots illustrating how systematic uncertainties are estimated form a certain source
 *
 *  JDrawer *drawer = JDrawer doing the dirty work in drawing
 *  TH1D* nominalResult = Histogram containing nominal results
 *  TH1D* variedResult[] = Array of histograms containing variation of results around the nominal one
 *  const int nVariations = Number of variations around the nominal result
 *  const int iCentrality = Centrality index of the considered bin
 *  const int iJetPt = Jet pT index of the considered bin
 *  const int iTrackPt = Track pT index of the considered bin
 *  EECCard* card = Card used to interpret the bin index
 *  TString comparisonLegend[] = An array of strings describing each result variation
 *  TString plotName = String added to saved plots. If left empty, the plots are not saved into files.
 *  std::pair<double, double> drawingRange = Drawing range for x-axis
 *  std::pair<double, double> ratioZoom = Y-axis zoom for the ratio
 */
void drawIllustratingPlots(JDrawer* drawer, TH1D* nominalResult, TH1D* variedResult[], const int nVariations, const int iCentrality, const int iJetPt, const int iTrackPt, EECCard* card, TString comparisonLegend[], TString plotName, std::pair<double, double> drawingRange, std::pair<double, double> ratioZoom){
  
  const int markers[] = {kFullDiamond, kFullDoubleDiamond, kFullCross, kFullFourTrianglesPlus, kFullSquare, kFullStar};
  const int colors[] = {kBlue, kRed, kGreen+3, kMagenta, kCyan, kViolet, kBlack};
  
  // Interpret the given binning information
  TString centralityString = Form("Cent: %.0f-%.0f%%", card->GetLowBinBorderCentrality(iCentrality), card->GetHighBinBorderCentrality(iCentrality));
  TString compactCentralityString = Form("_C=%.0f-%.0f", card->GetLowBinBorderCentrality(iCentrality), card->GetHighBinBorderCentrality(iCentrality));
  TString jetPtString = Form("%.0f < jet p_{T} < %.0f", card->GetLowBinBorderJetPtEEC(iJetPt), card->GetHighBinBorderJetPtEEC(iJetPt));
  TString compactJetPtString = Form("_J=%.0f-%.0f", card->GetLowBinBorderJetPtEEC(iJetPt), card->GetHighBinBorderJetPtEEC(iJetPt));
  TString trackPtString = Form("%.1f < track p_{T}",card->GetLowBinBorderTrackPtEEC(iTrackPt));
  TString compactTrackPtString = Form("_T>%.1f",card->GetLowBinBorderTrackPtEEC(iTrackPt));
  compactTrackPtString.ReplaceAll(".","v");

  // Calculate ratios between nominal and comparison graphs
  TH1D *ratioHistogram[nVariations];
  for(int iVariation = 0; iVariation < nVariations; iVariation++){
    ratioHistogram[iVariation] = (TH1D*) variedResult[iVariation]->Clone(Form("ratio_%s", comparisonLegend[iVariation].Data()));
    ratioHistogram[iVariation]->Divide(nominalResult);
  } // Varied results loop

  // Create a new canvas for the plot
  drawer->CreateSplitCanvas();

  // Use logarithmic axis for EEC
  drawer->SetLogY(true);

  // Setup the legend for plots
  TLegend *legend = new TLegend(0.23,0.3-nVariations*0.06,0.53,0.6);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  legend->AddEntry((TObject*) 0, Form("%s 5.02 TeV",card->GetAlternativeDataType().Data()), "");
  legend->AddEntry((TObject*) 0, centralityString.Data(),"");
  legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
  legend->AddEntry((TObject*) 0, trackPtString.Data(),"");

  // Set the x-axis drawing range
  nominalResult->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);

  // Draw the nominal result to the upper canvas
  nominalResult->SetLineColor(kBlack);
  drawer->DrawHistogramToUpperPad(nominalResult, "#Deltar", "EEC Signal", " ");
  legend->AddEntry(nominalResult, "Nominal", "l");

  // Draw the variations to the same plot
  for(int iVariation = 0; iVariation < nVariations; iVariation++) {
    variedResult[iVariation]->SetLineColor(colors[iVariation]);
    variedResult[iVariation]->Draw("same");
    legend->AddEntry(variedResult[iVariation], comparisonLegend[iVariation], "l");
  }

  // Draw the legends to the upper pad
  legend->Draw();

  // Linear scale for the ratio
  drawer->SetLogY(false);

  // Set the x-axis drawing range
  ratioHistogram[0]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);

  ratioHistogram[0]->SetLineColor(colors[0]);
  ratioHistogram[0]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
  drawer->SetGridY(true);
  drawer->DrawHistogramToLowerPad(ratioHistogram[0], "#Deltar", "#frac{Variation}{Nominal}", " ");
  for(int iVariation = 1; iVariation < nVariations; iVariation++) {
    ratioHistogram[iVariation]->SetLineColor(colors[iVariation]);
    ratioHistogram[iVariation]->Draw("same");
  }
  drawer->SetGridY(false);
  
  // If a plot name is given, save the plot in a file
  if(plotName != ""){
    gPad->GetCanvas()->SaveAs(Form("figures/systematicUncertainty_%s%s%s%s.pdf", plotName.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data()));
  }
}

/*
 * Draw plots illustrating how systematic uncertainties are estimated form a certain source
 *
 *  JDrawer* drawer = JDrawer doing the dirty work in drawing
 *  TH1D* nominalResult = Histogram containing nominal results
 *  TH1D* variedResult = Histogram containing variation of results around the nominal one
 *  const int iCentrality = Centrality index of the considered bin
 *  const int iJetPt = Jet pT index of the considered bin
 *  const int iTrackPt = Track pT index of the considered bin
 *  EECCard* card = Card used to interpret the binning information
 *  TString comparisonLegend = String describing the variation
 *  TString plotName = String added to saved plots. If left empty, the plots are not saved into files.
 *  std::pair<double, double> drawingRange = Drawing range for x-axis
 *  std::pair<double, double> ratioZoom = Y-axis zoom for the ratio
 */
void drawIllustratingPlots(JDrawer* drawer, TH1D* nominalResult, TH1D* variedResult, const int iCentrality, const int iJetPt, const int iTrackPt, EECCard* card, TString comparisonLegend, TString plotName, std::pair<double, double> drawingRange, std::pair<double, double> ratioZoom){
  TH1D* variedResultArray[1] = {variedResult};
  TString comparisonLegendArray[1] = {comparisonLegend};
  drawIllustratingPlots(drawer, nominalResult, variedResultArray, 1, iCentrality, iJetPt, iTrackPt, card, comparisonLegendArray, plotName, drawingRange, ratioZoom);
}

/*
 * Load histograms relevant for systematic uncertainty analysis from the histogram manager 
 */
void loadRelevantHistograms(EECHistogramManager *histograms){

  // Load energy-energy correlators
  histograms->SetLoadEnergyEnergyCorrelators(true);

  // Load all available bins
  histograms->SetCentralityBinRange(0, histograms->GetCard()->GetNCentralityBins() - 1);
  histograms->SetJetPtBinRangeEEC(0, histograms->GetCard()->GetNJetPtBinsEEC() - 1);
  histograms->SetTrackPtBinRangeEEC(0, histograms->GetCard()->GetNTrackPtBinsEEC() - 1);

  // Load the histograms from the file
  histograms->LoadProcessedHistograms();

}
