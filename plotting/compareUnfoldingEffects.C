#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Macro for comparing the unfolded-to-raw ratios in data and MC.
 * The prior assumption is that the effects to the distributions should be similar.
 */
void compareUnfoldingEffects(){
  
  // Enumeration for data and MC
  enum enumDataType{kData, kMCSplit1, kMCSplit2, knDataTypes};
  TString dataTypeString[knDataTypes] = {"data", "MC split 1", "MC split 2"};
  
  // Data and MC files for the comparison
  TString fileName[knDataTypes];
  fileName[kData] = "data/eecAnalysis_akFlowJet_wtaAxis_binningForUnfolding_processed_unfoldTest_2023-05-23.root";
  fileName[kMCSplit1] = "data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_moreLowPtBins_reconstructedReferenceForUnfolding_part2_processed_unfoldTest_2023-05-20.root";
  fileName[kMCSplit2] = "data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_moreLowPtBins_reconstructedReferenceForUnfolding_part1_processed_unfoldTest_2023-05-20.root";
  
  // Open the files and check that they exist
  TFile* inputFile[knDataTypes];
  EECCard* card[knDataTypes];
  for(int iFile = 0; iFile < knDataTypes; iFile++){
    
    inputFile[iFile] = TFile::Open(fileName[iFile]);
    
    if(inputFile[iFile] == NULL){
      cout << "Error! The file " << fileName[iFile].Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }

    card[iFile] = new EECCard(inputFile[iFile]);
  }
  
  // Check if we are using PbPb or pp data
  TString collisionSystem = card[kData]->GetDataType();
  bool isPbPbData = collisionSystem.Contains("PbPb");
  
  // ====================================================
  //               Binning configuration
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = card[kData]->GetNCentralityBins();
  const int nJetPtBinsEEC = card[kData]->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card[kData]->GetNTrackPtBinsEEC();
  
  // Default binning ranges for reference
  // centrality = {0,10,30,50,90}
  // jet pT in energy-energy correlators = {120,140,160,180,200,300,500,5030}
  // track pT in energy-energy correlators = {0.7,1,1.5,2,2.5,3,3.5,4,300}
  
  // Draw the analysis closure for all bins that have been unfolded
  int firstDrawnCentralityBin = card[kData]->GetFirstUnfoldedCentralityBin();
  int lastDrawnCentralityBin = card[kData]->GetLastUnfoldedCentralityBin();
  
  int firstDrawnJetPtBinEEC = card[kData]->GetFirstUnfoldedJetPtBin();
  int lastDrawnJetPtBinEEC = card[kData]->GetLastUnfoldedJetPtBin();
  
  int firstDrawnTrackPtBinEEC = card[kData]->GetFirstUnfoldedTrackPtBin();
  int lastDrawnTrackPtBinEEC = card[kData]->GetLastUnfoldedTrackPtBin();

  firstDrawnTrackPtBinEEC = 5;
  lastDrawnTrackPtBinEEC = 5;
  
  // Select the types of energy-energy correlators are studied
  bool studyEnergyEnergyCorrelator[EECHistogramManager::knEnergyEnergyCorrelatorTypes];
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator] = true;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationPlus] = false;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationMinus] = false;

  // ====================================================
  //                    Configuration
  // ====================================================
  
  // Figure saving
  const bool saveFigures = true;  // Save figures
  const char* saveComment = "_firstSanityCheck";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures

  // Drawing configuration
  std::pair<double, double> ratioZoom = std::make_pair(0.7, 1.3);
  
  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager* histograms[knDataTypes];
  for(int iFile = 0; iFile < knDataTypes; iFile++){
    histograms[iFile] = new EECHistogramManager(inputFile[iFile], card[iFile]);

    // Choose the energy-energy correlator types to load
    histograms[iFile]->SetLoadEnergyEnergyCorrelators(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator]);
    histograms[iFile]->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationPlus(integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationPlus]);
    histograms[iFile]->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationMinus(integrateEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationMinus]);

    // Choose the bin ranges
    histograms[iFile]->SetCentralityBinRange(0, card[iFile]->GetNCentralityBins() - 1);
    histograms[iFile]->SetJetPtBinRangeEEC(0, card[iFile]->GetNJetPtBinsEEC() - 1);
    histograms[iFile]->SetTrackPtBinRangeEEC(0, card[iFile]->GetNTrackPtBinsEEC() - 1);

    // Load the histograms from the file
    histograms[iFile]->LoadProcessedHistograms();
  }

  // Energy-energy correlator histograms
  TH1D* hEnergyEnergyCorrelatorRaw[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][knDataTypes];
  TH1D* hEnergyEnergyCorrelatorUnfolded[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][knDataTypes];
  TH1D* hEnergyEnergyCorrelatorRawToUnfoldedRatio[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC][knDataTypes]; // Ratio between raw and unfolded distributions
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
          for(int iSignalType = 0; iSignalType < knDataTypes; iSignalType++){
            hEnergyEnergyCorrelatorRaw[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSignalType] = NULL;
            hEnergyEnergyCorrelatorUnfolded[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSignalType] = NULL;
            hEnergyEnergyCorrelatorRawToUnfoldedRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iSignalType] = NULL;
          } // Signal type loop
        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
  } // Energy-energy correlator type loop
  
  // Helper histograms
  TH1D* helperHistogram;
  double normalizationFactor;
  std::pair<double, double> drawingRange = std::make_pair(0.006, 0.39);
  std::pair<double, double> centralityBinBorders;
  double epsilon = 0.0001;
  int lowNormalizationBin, highNormalizationBin;
  int iCentralityMatched[knDataTypes];
  int iTrackPtMatched[knDataTypes];
  int iJetPtMatched[knDataTypes];

  // Get the histograms from the histogram manager and normalize the signal histograms to one
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){

    // Only load the selected energy-energy correlator types
    if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;

    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      for(int iFile = 0; iFile < knDataTypes; iFile++){
        iCentralityMatched[iFile] = card[iFile]->FindBinIndexCentrality(card[kData]->GetBinBordersCentrality(iCentrality));

        // If we do not find a centrality bin match, try to shift the borders by 4 units to match between data and MC
        if(iCentralityMatched[iFile] < 0){
          centralityBinBorders = card[kData]->GetBinBordersCentrality(iCentrality);
          centralityBinBorders.first = centralityBinBorders.first + 4;
          centralityBinBorders.second = centralityBinBorders.second + 4;
          iCentralityMatched[iFile] = card[iFile]->FindBinIndexCentrality(centralityBinBorders);
        }
      }
      for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
        for(int iFile = 0; iFile < knDataTypes; iFile++){
          iTrackPtMatched[iFile] = card[iFile]->FindBinIndexTrackPtEEC(card[kData]->GetBinBordersTrackPtEEC(iTrackPt));
        }
        for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
          for(int iFile = 0; iFile < knDataTypes; iFile++){
            iJetPtMatched[iFile] = card[iFile]->FindBinIndexJetPtEEC(card[kData]->GetBinBordersJetPtEEC(iJetPt));

            // Read the raw and unfolded distributions from the file
            hEnergyEnergyCorrelatorRaw[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile] = histograms[iFile]->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentralityMatched[iFile], iJetPtMatched[iFile], iTrackPtMatched[iFile]);
            hEnergyEnergyCorrelatorUnfolded[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile] = histograms[iFile]->GetHistogramEnergyEnergyCorrelatorProcessed(iEnergyEnergyCorrelator, iCentralityMatched[iFile], iJetPtMatched[iFile], iTrackPtMatched[iFile], EECHistogramManager::kEnergyEnergyCorrelatorUnfolded);

            // Normalize the distributions to one in the drawingRange
            lowNormalizationBin = hEnergyEnergyCorrelatorRaw[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile]->GetXaxis()->FindBin(drawingRange.first + epsilon);
            highNormalizationBin = hEnergyEnergyCorrelatorRaw[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile]->GetXaxis()->FindBin(drawingRange.second - epsilon);

            hEnergyEnergyCorrelatorRaw[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile]->Scale(1 / hEnergyEnergyCorrelatorRaw[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile]->Integral(lowNormalizationBin, highNormalizationBin, "width"));
            hEnergyEnergyCorrelatorUnfolded[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile]->Scale(1 / hEnergyEnergyCorrelatorUnfolded[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile]->Integral(lowNormalizationBin, highNormalizationBin, "width"));

            // Calculate the raw to unfolded ratio
            hEnergyEnergyCorrelatorRawToUnfoldedRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile] = (TH1D*)hEnergyEnergyCorrelatorRaw[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile]->Clone(Form("rawToUnfoldedRatio%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iTrackPt, iJetPt, iFile));

            hEnergyEnergyCorrelatorRawToUnfoldedRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile]->Divide(hEnergyEnergyCorrelatorUnfolded[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile]);

          } // File type loop
        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
  } // Energy-energy correlator type loop

  // TODO: The code below needs to be updated
  
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

  TLegend* legend;
  TString centralityString, trackPtString, jetPtString;
  TString compactCentralityString, compactTrackPtString, compactJetPtString;
  TString legendString;
  int markerStyle[2] = {kOpenSquare, kOpenCircle};
  int color[knDataTypes] = {kBlack, kRed, kBlue};
  
  // Loop over all selected histograms
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only draw the selected energy-energy correlator types
    if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
    
    // Loop over centrality bins
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      
      // Set the centrality information for legends and figure saving
      centralityString = Form("Cent: %.0f-%.0f%%", histograms[kData]->GetCentralityBinBorder(iCentrality), histograms[kData]->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f", histograms[kData]->GetCentralityBinBorder(iCentrality), histograms[kData]->GetCentralityBinBorder(iCentrality+1));
      
      // Loop over jet pT bins
      for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
        
        // Set the jet pT information for legends and figure saving
        if(iJetPt == histograms[kData]->GetNJetPtBinsEEC()){
          jetPtString = Form("Jet p_{T} > %.0f", card[kData]->GetJetPtCut());
          compactJetPtString = "";
        } else {
          jetPtString = Form("%.0f < jet p_{T} < %.0f", histograms[kData]->GetJetPtBinBorderEEC(iJetPt), histograms[kData]->GetJetPtBinBorderEEC(iJetPt+1));
          compactJetPtString = Form("_J=%.0f-%.0f", histograms[kData]->GetJetPtBinBorderEEC(iJetPt), histograms[kData]->GetJetPtBinBorderEEC(iJetPt+1));
        }
        
        // Loop over track pT bins
        for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
          
          trackPtString = Form("%.1f < track p_{T}",histograms[kData]->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString = Form("_T>%.1f",histograms[kData]->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString.ReplaceAll(".","v");
          
          // Create a legend for the figure
          legend = new TLegend(0.18,0.04,0.45,0.48);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          legend->AddEntry((TObject*) 0, Form("%s 5.02 TeV",card[kData]->GetAlternativeDataType().Data()), "");
          legend->AddEntry((TObject*) 0, centralityString.Data(),"");
          legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
          legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
          
          // Create a new canvas for the plot
          drawer->CreateSplitCanvas();
          
          // Logarithmic EEC axis
          drawer->SetLogY(true);

          // Set drawing style for all histograms
          for(int iFile = 0; iFile < knDataTypes; iFile++){
            hEnergyEnergyCorrelatorRaw[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile]->SetMarkerStyle(markerStyle[0]);
            hEnergyEnergyCorrelatorRaw[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile]->SetMarkerColor(color[iFile]);
            hEnergyEnergyCorrelatorRaw[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile]->SetLineColor(color[iFile]);
            hEnergyEnergyCorrelatorUnfolded[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile]->SetMarkerStyle(markerStyle[1]);
            hEnergyEnergyCorrelatorUnfolded[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile]->SetMarkerColor(color[iFile]);
            hEnergyEnergyCorrelatorUnfolded[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile]->SetLineColor(color[iFile]);
          }
          
          // Set the x-axis drawing range
          hEnergyEnergyCorrelatorRaw[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
          
          // Draw the histograms to the upper canvas
          drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelatorRaw[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0], "#Deltar", "EEC", " ");
          hEnergyEnergyCorrelatorUnfolded[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0]->Draw("same");

          for(int iFile = 1; iFile < knDataTypes; iFile++){
            hEnergyEnergyCorrelatorRaw[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile]->Draw("same");
            hEnergyEnergyCorrelatorUnfolded[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile]->Draw("same");
          }

          for(int iFile = 0; iFile < knDataTypes; iFile++){
            legend->AddEntry(hEnergyEnergyCorrelatorRaw[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile], Form("Raw EEC, %s", dataTypeString[iFile].Data()), "p");
            legend->AddEntry(hEnergyEnergyCorrelatorUnfolded[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile], Form("Unfolded EEC, %s", dataTypeString[iFile].Data()), "p");
          }
          
          // Draw the legends to the upper pad
          legend->Draw();
          
          // Linear scale for the ratio
          drawer->SetLogY(false);

          // Set the drawing style for the histograms
          for(int iFile = 0; iFile < knDataTypes; iFile++){
            hEnergyEnergyCorrelatorRawToUnfoldedRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile]->SetMarkerStyle(markerStyle[0]);
            hEnergyEnergyCorrelatorRawToUnfoldedRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile]->SetMarkerColor(color[iFile]);
            hEnergyEnergyCorrelatorRawToUnfoldedRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile]->SetLineColor(color[iFile]);
          }
          
          // Set the drawing ranges
          hEnergyEnergyCorrelatorRawToUnfoldedRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);          
          hEnergyEnergyCorrelatorRawToUnfoldedRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

          // Draw the histograms
          drawer->SetGridY(true);
          drawer->DrawHistogramToLowerPad(hEnergyEnergyCorrelatorRawToUnfoldedRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][0], "#Deltar", "#frac{Raw}{Unfolded}", " ");
          for(int iFile = 1; iFile < knDataTypes; iFile++){
            hEnergyEnergyCorrelatorRawToUnfoldedRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iFile]->Draw("same");
          }
          drawer->SetGridY(false);
          
          // Save the figures to a file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/%sUnfoldingPerformanceComparison%s%s%s%s.%s", histograms[kData]->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator), saveComment, compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
          }
          
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Energy-energy correlator type loop
}
