#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Macro for making closure plots for the analysis. It compares fully unfolded results to signal from MC truth.
 */
void fullAnalysisClosure(){

  // Enumeration for distribution type
  enum enumDistributionType{kMeasured, kTruth, kNDistributionTypes};
  const int nSplits = 2;

  // Open the input files
  TString fileName[kNDistributionTypes][nSplits];
  fileName[kMeasured][0] = "data/PbPbMC2018_RecoReco_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_noJetPtWeight_part2_processed_closureTest_2023-06-01.root";
  fileName[kMeasured][1] = "data/PbPbMC2018_RecoReco_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_noJetPtWeight_part1_processed_closureTest_2023-06-01.root";
  // data/PbPbMC2018_RecoReco_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_noJetPtWeight_part1_processed_closureTest_2023-06-01.root
  fileName[kTruth][0] = "data/PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_moreLowPtBins_truthReferenceForUnfolding_part2_processed_2023-05-20.root";
  fileName[kTruth][1] = "data/PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_moreLowPtBins_truthReferenceForUnfolding_part1_processed_2023-05-20.root";
  // data/PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_finalMcWeight_processed_2023-03-08.root
  // data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_finalMcWeight_matchJets_processed_2023-03-06.root
  // data/PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_moreLowPtBins_truthReferenceForUnfolding_part2_processed_2023-05-20.root

  TFile* inputFile[kNDistributionTypes][nSplits];
  EECCard* card[kNDistributionTypes][nSplits];

  for(int iFile = 0; iFile < kNDistributionTypes; iFile++){
    for(int iSplit = 0; iSplit < nSplits; iSplit++){
      inputFile[iFile][iSplit] = TFile::Open(fileName[iFile][iSplit]);
  
      // Check that the files exist
      if(inputFile[iFile][iSplit] == NULL){
        cout << "Error! The file " << fileName[iFile][iSplit].Data() << " does not exist!" << endl;
        cout << "Maybe you forgot the data/ folder path?" << endl;
        cout << "Will not execute the code" << endl;
        return;
      }

      card[iFile][iSplit] = new EECCard(inputFile[iFile][iSplit]);
    }
  }

  // It is assumed that the different splits have the same binning. It might be worth implementing a check here to avoid bugs producing scary closures.
  
  // ====================================================
  //               Binning configuration
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = card[kMeasured][0]->GetNCentralityBins();
  const int nJetPtBinsEEC = card[kMeasured][0]->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card[kMeasured][0]->GetNTrackPtBinsEEC();
  
  // Default binning ranges for reference
  // centrality = {0,10,30,50,90}
  // jet pT in energy-energy correlators = {120,140,160,180,200,300,500,5030}
  // track pT in energy-energy correlators = {0.7,1,1.5,2,2.5,3,3.5,4,300}
  
  // Draw the analysis closure for all bins that have been unfolded
  int firstDrawnCentralityBin = card[kMeasured][0]->GetFirstUnfoldedCentralityBin();
  int lastDrawnCentralityBin = card[kMeasured][0]->GetLastUnfoldedCentralityBin();
  
  int firstDrawnJetPtBinEEC = card[kMeasured][0]->GetFirstUnfoldedJetPtBin();
  int lastDrawnJetPtBinEEC = card[kMeasured][0]->GetLastUnfoldedJetPtBin();
  
  int firstDrawnTrackPtBinEEC = card[kMeasured][0]->GetFirstUnfoldedTrackPtBin();
  int lastDrawnTrackPtBinEEC = card[kMeasured][0]->GetLastUnfoldedTrackPtBin();
  
  // Select the types of energy-energy correlators are studied
  bool studyEnergyEnergyCorrelator[EECHistogramManager::knEnergyEnergyCorrelatorTypes];
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator] = true;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationPlus] = false;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationMinus] = false;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationPlus] = false;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationMinus] = false;
  
  // Logarithmic axes
  const bool logDeltaR = true;
  const bool logEEC = true;
  
  // Axis zooming
  std::pair<double,double> ratioZoom = std::make_pair(0.7, 1.3);
  
  // Figure saving
  const bool saveFigures = true;  // Save figures
  const char* saveComment = "_noSystematicUncertainties";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures
  
  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager *histograms[kNDistributionTypes][nSplits];
  for(int iFile = 0; iFile < kNDistributionTypes; iFile++){
    for(int iSplit = 0; iSplit < nSplits; iSplit++){
      histograms[iFile][iSplit] = new EECHistogramManager(inputFile[iFile][iSplit], card[iFile][iSplit]);

      // Choose the energy-energy correlator types to load
      histograms[iFile][iSplit]->SetLoadEnergyEnergyCorrelators(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator]);
      histograms[iFile][iSplit]->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationPlus(integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationPlus]);
      histograms[iFile][iSplit]->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationMinus(integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationMinus]);
      histograms[iFile][iSplit]->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationPlus(integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationPlus]);
      histograms[iFile][iSplit]->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationMinus(integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationMinus]);

      // Choose the bin ranges
      histograms[iFile][iSplit]->SetCentralityBinRange(0, card[iFile][iSplit]->GetNCentralityBins() - 1);
      histograms[iFile][iSplit]->SetJetPtBinRangeEEC(0, card[iFile][iSplit]->GetNJetPtBinsEEC() - 1);
      histograms[iFile][iSplit]->SetTrackPtBinRangeEEC(0, card[iFile][iSplit]->GetNTrackPtBinsEEC() - 1);

      // Load the histograms from the file
      histograms[iFile][iSplit]->LoadProcessedHistograms();
    }
  }

  // Energy-energy correlator histograms
  TH1D* hEnergyEnergyCorrelatorSignal[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][kNDistributionTypes][nSplits]; // Last bin, true signal/extracted signal
  TH1D* hEnergyEnergyCorrelatorSignalRatio[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][nSplits]; // Ratio between true and extracted signal
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
          for(int iSplit = 0; iSplit < nSplits; iSplit++){
            for(int iSignalType = 0; iSignalType < kNDistributionTypes; iSignalType++){
              hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSignalType][iSplit] = NULL;
            } // Signal type loop
            hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSplit] = NULL;
          } // Split loop
        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
  } // Energy-energy correlator type loop
  
  // Helper histograms
  TH1D *helperHistogram;
  double normalizationFactor;
  std::pair<double,double> drawingRange = std::make_pair(0.006, 0.39);
  double epsilon = 0.0001;
  int lowSignalBin, highSignalBin;
  int iCentralityTruth, iJetPtTruth, iTrackPtTruth;
  
  // Get the histograms from the histogram manager and normalize the signal histograms to one
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only load the selected energy-energy correlator types
    if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;

    for(int iSplit = 0; iSplit < nSplits; iSplit++){

      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        iCentralityTruth = card[kTruth][iSplit]->FindBinIndexCentrality(card[kMeasured][iSplit]->GetBinBordersCentrality(iCentrality));
        for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
          iTrackPtTruth = card[kTruth][iSplit]->FindBinIndexTrackPtEEC(card[kMeasured][iSplit]->GetBinBordersTrackPtEEC(iTrackPt));
          for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
            iJetPtTruth = card[kTruth][iSplit]->FindBinIndexJetPtEEC(card[kMeasured][iSplit]->GetBinBordersJetPtEEC(iJetPt));

            // Read the signal from the unfolded measurement
            hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kMeasured][iSplit] = histograms[kMeasured][iSplit]->GetHistogramEnergyEnergyCorrelatorProcessed(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);

            // Read the signal from MC truth
            hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTruth][iSplit] = histograms[kTruth][iSplit]->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentralityTruth, iJetPtTruth, iTrackPtTruth, EECHistograms::kSameJetPair, EECHistograms::kPythiaPythia);

            // Normalize the signal distributions to one in the drawingRange
            lowSignalBin = hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kMeasured][iSplit]->GetXaxis()->FindBin(drawingRange.first + epsilon);
            highSignalBin = hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kMeasured][iSplit]->GetXaxis()->FindBin(drawingRange.second - epsilon);
            for(int iSignalType = 0; iSignalType < kNDistributionTypes; iSignalType++){
              hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSignalType][iSplit]->Scale(1 / hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSignalType][iSplit]->Integral(lowSignalBin, highSignalBin, "width"));
            }

            // Calculate the extracted to true signal ratio
            hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSplit] = (TH1D*)hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kMeasured][iSplit]->Clone(Form("signalRatio%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iTrackPt, iJetPt, iSplit));

            hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSplit]->Divide(hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTruth][iSplit]);

          } // Jet pT loop
        } // Track pT loop
      } // Centrality loop
    } // Split loop
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
  TString centralityString, trackPtString, jetPtString;
  TString compactCentralityString, compactTrackPtString, compactJetPtString;
  TString legendString;
  
  // Loop over all selected histograms
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only draw the selected energy-energy correlator types
    if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
    
    // Loop over centrality bins
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      
      // Set the centrality information for legends and figure saving
      centralityString = Form("Cent: %.0f-%.0f%%", histograms[kMeasured][0]->GetCentralityBinBorder(iCentrality), histograms[kMeasured][0]->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f", histograms[kMeasured][0]->GetCentralityBinBorder(iCentrality), histograms[kMeasured][0]->GetCentralityBinBorder(iCentrality+1));
      
      // Loop over jet pT bins
      for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
        
        // Set the jet pT information for legends and figure saving
        if(iJetPt == histograms[kMeasured][0]->GetNJetPtBinsEEC()){
          jetPtString = Form("Jet p_{T} > %.0f", card[kMeasured][0]->GetJetPtCut());
          compactJetPtString = "";
        } else {
          jetPtString = Form("%.0f < jet p_{T} < %.0f", histograms[kMeasured][0]->GetJetPtBinBorderEEC(iJetPt), histograms[kMeasured][0]->GetJetPtBinBorderEEC(iJetPt+1));
          compactJetPtString = Form("_J=%.0f-%.0f", histograms[kMeasured][0]->GetJetPtBinBorderEEC(iJetPt), histograms[kMeasured][0]->GetJetPtBinBorderEEC(iJetPt+1));
        }
        
        // Loop over track pT bins
        for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
          
          trackPtString = Form("%.1f < track p_{T}",histograms[kMeasured][0]->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString = Form("_T>%.1f",histograms[kMeasured][0]->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString.ReplaceAll(".","v");
          
          // Create a legend for the figure
          legend = new TLegend(0.18,0.04,0.45,0.48);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          legend->AddEntry((TObject*) 0, Form("%s 5.02 TeV",card[kMeasured][0]->GetAlternativeDataType().Data()), "");
          legend->AddEntry((TObject*) 0, centralityString.Data(),"");
          legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
          legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
          
          // Create a new canvas for the plot
          drawer->CreateSplitCanvas();
          
          // Logarithmic EEC axis
          if(logEEC) drawer->SetLogY(true);
          
          // Set the x-axis drawing range
          hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTruth][0]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
          
          // Draw the histograms to the upper canvas
          hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTruth][0]->SetLineColor(kBlack);
          drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTruth][0], "#Deltar", "EEC Signal", " ");
          
          hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kMeasured][0]->SetLineColor(kRed);
          hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kMeasured][0]->Draw("same");

          for(int iSplit = 1; iSplit < nSplits; iSplit++){
            hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTruth][iSplit]->SetLineColor(kGreen+3);
            hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTruth][iSplit]->Draw("same");
            hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kMeasured][iSplit]->SetLineColor(kBlue);
            hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kMeasured][iSplit]->Draw("same");
          }

          for(int iSplit = 0; iSplit < nSplits; iSplit++){
            legend->AddEntry(hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kTruth][iSplit], Form("True signal, split %d", iSplit+1), "l");
            legend->AddEntry(hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][kMeasured][iSplit], Form("Unfolded signal, split %d", iSplit+1), "l");
          }
          
          // Draw the legends to the upper pad
          legend->Draw();
          
          // Linear scale for the ratio
          drawer->SetLogY(false);
          
          // Set the x-axis drawing range
          hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
          
          hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0]->SetLineColor(kRed);
          hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
          drawer->SetGridY(true);
          drawer->DrawHistogramToLowerPad(hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0], "#Deltar", "#frac{Unfolded}{True}", " ");
          for(int iSplit = 1; iSplit < nSplits; iSplit++){
            hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSplit]->SetLineColor(kBlue);
            hEnergyEnergyCorrelatorSignalRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSplit]->Draw("same");
          }
          drawer->SetGridY(false);
          
          // Save the figures to a file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/%sFullAnalysisClosure%s%s%s%s.%s", histograms[kMeasured][0]->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator), saveComment, compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
          }
          
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Energy-energy correlator type loop
}
