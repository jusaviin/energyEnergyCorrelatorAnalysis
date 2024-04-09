#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Macro for validating the background subtraction algorithm. In Pythia+Hydjet simulation, we can directly compare the extracted
 * signal to the true signal. It is done in this macro. Note that the input file must contain both regular and processed histograms.
 */
void validateBackgroundSubtraction(){

  // Open the input file
  TString inputFileName = "data/PbPbMC2018_GenGen_eecAnalysis_4pCentShift_cutBadPhi_nominalEnergyWeight_mixedEventBackground_someJobsMissing_processed_2024-04-01.root";
  // data/PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_finalMcWeight_processed_2023-03-08.root
  // data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_finalMcWeight_matchJets_processed_2023-03-06.root
  // data/PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_moreLowPtBins_truthReferenceForUnfolding_part2_processed_2023-05-20.root
  // data/PbPbMC2018_GenGen_akFlowJets_4pCentShift_cutBadPhi_optimizedUnfoldingBins_nominalSmear_truthReference_processed_2024-01-16.root
  // data/PbPbMC2018_GenGen_akFlowJets_4pCentShift_cutBadPhi_optimizedUnfoldingBins_energyWeightSquared_nominalSmear_truthReference_processed_2024-01-18.root
  TFile* inputFile = TFile::Open(inputFileName);
  
  // Check that the files exist
  if(inputFile == NULL){
    cout << "Error! The file " << inputFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }

  EECCard* card = new EECCard(inputFile);
  const int weightExponent = card->GetWeightExponent();
  
  // ====================================================
  //               Binning configuration
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = card->GetNCentralityBins();
  const int nJetPtBinsEEC = card->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card->GetNTrackPtBinsEEC();
  
  // Default binning ranges for reference
  // centrality = {0,10,30,50,90}
  // jet pT in energy-energy correlators = {120,140,160,180,200,300,500,5030}
  // track pT in energy-energy correlators = {0.7,1,1.5,2,2.5,3,3.5,4,300}
  
  // Select which histograms are fitted
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = nCentralityBins-1;
  
  int firstDrawnJetPtBinEEC = 5;
  int lastDrawnJetPtBinEEC = 8; // Note: Jets integrated over all pT ranges are in nJetPtBinsEEC bin
  
  int firstDrawnTrackPtBinEEC = 1;
  int lastDrawnTrackPtBinEEC = 1;
  
  // Select the types of energy-energy correlators are studied
  bool studyEnergyEnergyCorrelator[EECHistogramManager::knEnergyEnergyCorrelatorTypes];
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator] = true;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationPlus] = false;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationMinus] = false;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationPlus] = false;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationMinus] = false;
  
  // Find index of one energy-energy correlator that is drawn
  int studiedEnergyEnergyCorrelatorType = -1;
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    if(studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]){
      studiedEnergyEnergyCorrelatorType = iEnergyEnergyCorrelator;
      break;
    }
  }
  
  // Logarithmic axes
  const bool logDeltaR = true;
  const bool logEEC = true;
  
  // Axis zooming
  std::pair<double,double> ratioZoom = std::make_pair(0.7, 1.3);

  // Extra tag for the plots
  const bool simulationTag = false;
  
  // Figure saving
  const bool saveFigures = false;  // Save figures
  const char* saveComment = "_energyWeightSquared_newBackground";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures
  
  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager *histograms;
  histograms = new EECHistogramManager(inputFile,card);
    
  // Choose the energy-energy correlator types to load
  histograms->SetLoadEnergyEnergyCorrelators(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator]);
  histograms->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationPlus(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationPlus]);
  histograms->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationMinus(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationMinus]);
  histograms->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationPlus(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationPlus]);
  histograms->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationMinus(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationMinus]);
    
  // Choose the bin ranges
  histograms->SetCentralityBinRange(0,nCentralityBins-1);
  histograms->SetJetPtBinRangeEEC(0,nJetPtBinsEEC);
  histograms->SetTrackPtBinRangeEEC(0,nTrackPtBinsEEC-1);
    
  // Load the histograms from the file
  histograms->LoadProcessedHistograms();
  
  // Energy-energy correlator histograms
  TH1D* hEnergyEnergyCorrelatorSignal[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC][2]; // Last bin, true signal/extracted signal
  TH1D* hEnergyEnergyCorrelatorSignalRatio[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC]; // Ratio between true and extracted signal
  TH1D* hEnergyEnergyCorrelatorBackground[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC][2]; // Last bin, true background/extracted background
  TH1D* hEnergyEnergyCorrelatorBackgroundRatio[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC]; // Ratio between true and extracted background
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        for(int iJetPt = 0; iJetPt < nJetPtBinsEEC+1; iJetPt++){
          for(int iSignalType = 0; iSignalType < 2; iSignalType++){
            hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSignalType] = NULL;
            hEnergyEnergyCorrelatorBackground[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSignalType] = NULL;
          } // Signal type loop
          hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = NULL;
          hEnergyEnergyCorrelatorBackgroundRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = NULL;
        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
  } // Energy-energy correlator type loop
  
  // Helper histograms
  TH1D* helperHistogram;
  double normalizationFactor;
  std::pair<double,double> drawingRange = std::make_pair(0.006, 0.39);
  double epsilon = 0.0000001;
  int lowSignalBin, highSignalBin;
  
  // Get the histograms from the histogram manager and normalize the signal histograms to one
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only load the selected energy-energy correlator types
    if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
    
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
        for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
          
          // Read the true energy-energy correlator signal
          hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0] = histograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJetPair, EECHistograms::kPythiaPythia);
          
          // The extracted energy-energy correlator signal
          hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][1] = histograms->GetHistogramEnergyEnergyCorrelatorProcessed(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistogramManager::kEnergyEnergyCorrelatorSignal);
          
          // Normalize the signal distributions to one in the region 0 to 0.4
          lowSignalBin = hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0]->GetXaxis()->FindBin(drawingRange.first+epsilon);
          highSignalBin = hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0]->GetXaxis()->FindBin(drawingRange.second-epsilon);
          for(int iSignalType = 0; iSignalType < 2; iSignalType++){
            hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSignalType]->Scale(1/hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSignalType]->Integral(lowSignalBin,highSignalBin,"width"));
          }
          
          // Calculate the extracted to true signal ratio
          hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][1]->Clone(Form("signalRatio%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iTrackPt, iJetPt));
          
          hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Divide(hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0]);
          
        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
  } // Energy-energy correlator type loop

  // ==========================================================================
  //        Draw the selected energy-energy correlator signal ratios
  // ==========================================================================
  
  JDrawer* drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  
  if(logDeltaR) drawer->SetLogX(true);

  TLegend* legend;
  TLegend* simulationLegend;
  TString centralityString, trackPtString, jetPtString;
  TString compactCentralityString, compactTrackPtString, compactJetPtString;
  TString legendString;
  double legendX1, legendX2;
  
  // Loop over all selected histograms
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only draw the selected energy-energy correlator types
    if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
    
    // Loop over centrality bins
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      
      // Set the centrality information for legends and figure saving
      centralityString = Form("Cent: %.0f-%.0f%%", histograms->GetCentralityBinBorder(iCentrality), histograms->GetCentralityBinBorder(iCentrality+1));
      //centralityString = "Cent: 0-10%";
      compactCentralityString = Form("_C=%.0f-%.0f", histograms->GetCentralityBinBorder(iCentrality), histograms->GetCentralityBinBorder(iCentrality+1));
      
      // Loop over jet pT bins
      for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
        
        // Set the jet pT information for legends and figure saving
        if(iJetPt == histograms->GetNJetPtBinsEEC()){
          jetPtString = Form("Jet p_{T} > %.0f", histograms->GetCard()->GetJetPtCut());
          compactJetPtString = "";
        } else {
          jetPtString = Form("%.0f < jet p_{T} < %.0f", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1));
          compactJetPtString = Form("_J=%.0f-%.0f", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1));
        }
        
        // Loop over track pT bins
        for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
          
          trackPtString = Form("%.1f < track p_{T}",histograms->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString = Form("_T>%.1f",histograms->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString.ReplaceAll(".","v");
          
          // Create a legend for the figure
          legend = new TLegend(0.2,0.04,0.45,0.44);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          legend->AddEntry((TObject*) 0, Form("%s 5.02 TeV",histograms->GetCard()->GetAlternativeDataType(false).Data()), "");
          //legend->AddEntry((TObject*) 0, "Energy weight squared","");
          legend->AddEntry((TObject*) 0, centralityString.Data(),"");
          legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
          legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
          
          // Create a new canvas for the plot
          drawer->CreateSplitCanvas();
          
          // Logarithmic EEC axis
          if(logEEC) drawer->SetLogY(true);
          
          // Set the x-axis drawing range
          hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
          
          // Draw the histograms to the upper canvas
          hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0]->SetLineColor(kBlack);
          drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0], "#Deltar", "EEC Signal", " ");
          legend->AddEntry(hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0], "True signal", "l");
          
          hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][1]->SetLineColor(kRed);
          hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][1]->Draw("same");
          legend->AddEntry(hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][1], "Extracted signal", "l");
          
          // Draw the legends to the upper pad
          legend->Draw();

          // Add simulation tag to the figures
          if(simulationTag){
            legendX1 = (weightExponent == 1) ? 0.23 : 0.18;
            legendX2 = (weightExponent == 1) ? 0.48 : 0.43;
            simulationLegend = new TLegend(legendX1,0.5,legendX2,0.65);
            simulationLegend->SetFillStyle(0);simulationLegend->SetBorderSize(0);simulationLegend->SetTextSize(0.05);simulationLegend->SetTextFont(62);
            simulationLegend->AddEntry((TObject*) 0, "CMS simulation +","");
            simulationLegend->AddEntry((TObject*) 0, "   Private work","");
            simulationLegend->Draw();
          }
          
          // Linear scale for the ratio
          drawer->SetLogY(false);
          
          // Set the x-axis drawing range
          hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
          
          hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetLineColor(kRed);
          hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
          drawer->SetGridY(true);
          drawer->DrawHistogramToLowerPad(hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt], "#Deltar", "#frac{Extracted}{True}", " ");
          drawer->SetGridY(false);
          
          // Save the figures to a file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/%sSignalValidityCheck%s%s%s%s.%s", histograms->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator), saveComment, compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
          }
          
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Energy-energy correlator type loop
}
