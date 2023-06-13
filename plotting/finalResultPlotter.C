#include "SystematicUncertaintyOrganizer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECHistogramManager.h"
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Macro for making final result plots comparing energy-energy correlators between pp and PbPb
 */
void finalResultPlotter(){

  enum enumDataType{kPbPb, kPp, kNDataTypes};

  // ============= //
  // Configuration //
  // ============= //
  
  // Input files
  TString inputFileName[kNDataTypes];
  inputFileName[kPbPb] = "data/eecAnalysis_akFlowJet_wtaAxis_binningForUnfolding_processed_nominalResult_2023-05-23.root";
  inputFileName[kPp] = "data/ppData_pfJets_wtaAxis_nominalResults_processed_2023-06-12.root";
  TString uncertaintyFileName[kNDataTypes];
  uncertaintyFileName[kPbPb] = "systematicUncertainties/systematicUncertainties_firstLook_2023-06-13.root";
  uncertaintyFileName[kPp] = "systematicUncertainties/systematicUncertaintiesForPp_firstLook_2023-06-13.root";
  
  TFile* inputFile[kNDataTypes];
  TFile* uncertaintyFile[kNDataTypes];
  EECCard* card[kNDataTypes];
  EECCard* uncertaintyCard[kNDataTypes];
  for(int iFile = 0; iFile < kNDataTypes; iFile++){

    // Load the input file
    inputFile[iFile] = TFile::Open(inputFileName[iFile]);

    // Check that the input file exists
    if(inputFile[iFile] == NULL){
      cout << "Error! The file " << inputFileName[iFile].Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }
    
    // Load the card from the file and read the collision system
    card[iFile] = new EECCard(inputFile[iFile]);

    // Load the uncertainty file
    uncertaintyFile[iFile] = TFile::Open(uncertaintyFileName[iFile]);

    // Check that the uncertianty file exists
    if(uncertaintyFile[iFile] == NULL){
      cout << "Error! The file " << uncertaintyFileName[iFile].Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the systematicUncertainties/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }
    
    // Load the card from the file and read the collision system
    uncertaintyCard[iFile] = new EECCard(uncertaintyFile[iFile]);
  }

  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Find the number of bins from the PbPb card
  const int nCentralityBins = card[kPbPb]->GetNCentralityBins();
  const int nJetPtBinsEEC = card[kPbPb]->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card[kPbPb]->GetNTrackPtBinsEEC();
  
  // The final results are available for all the bins that are unfolded
  int firstDrawnCentralityBin = card[kPbPb]->GetFirstUnfoldedCentralityBin();
  int lastDrawnCentralityBin = card[kPbPb]->GetLastUnfoldedCentralityBin();
  
  int firstDrawnJetPtBinEEC = card[kPbPb]->GetFirstUnfoldedJetPtBin();
  int lastDrawnJetPtBinEEC = card[kPbPb]->GetLastUnfoldedJetPtBin();
  
  int firstDrawnTrackPtBinEEC = card[kPbPb]->GetFirstUnfoldedTrackPtBin();
  int lastDrawnTrackPtBinEEC = card[kPbPb]->GetLastUnfoldedTrackPtBin();

  firstDrawnTrackPtBinEEC = 3;
  lastDrawnTrackPtBinEEC = 5;
  
  // Save the final plots
  const bool saveFigures = false;
  TString saveComment = "_preliminary";

  // Ratio zoom settings
  std::pair<double, double> ratioZoom = std::make_pair(0.4, 1.6);
  
  // Marker colors and styles
  int markerStyle[] = {kFullSquare, kFullCircle, kFullCross, kFullFourTrianglesPlus};
  int markerColor[] = {kRed, kBlue, kMagenta, kGreen+3};

  // =============================================== //
  // Read the histograms from the histogram managers //
  // =============================================== //

  // Create histogram managers for result files and systematic uncertainty organizers for systematic uncertainty files
  EECHistogramManager* histograms[kNDataTypes];
  SystematicUncertaintyOrganizer* uncertainties[kNDataTypes];
  
  for(int iDataType = 0; iDataType < kNDataTypes; iDataType++){

    // Create a new histogram manager
    histograms[iDataType] = new EECHistogramManager(inputFile[iDataType], card[iDataType]);
  
    // Load all unfolded energy-energy correlators
    histograms[iDataType]->SetLoadEnergyEnergyCorrelators(true);
    histograms[iDataType]->SetCentralityBinRange(card[iDataType]->GetFirstUnfoldedCentralityBin(), card[iDataType]->GetLastUnfoldedCentralityBin());
    histograms[iDataType]->SetTrackPtBinRangeEEC(card[iDataType]->GetFirstUnfoldedTrackPtBin(), card[iDataType]->GetLastUnfoldedTrackPtBin());
    histograms[iDataType]->SetJetPtBinRangeEEC(card[iDataType]->GetFirstUnfoldedJetPtBin(), card[iDataType]->GetLastUnfoldedJetPtBin());

    // Load the histograms from the file
    histograms[iDataType]->LoadProcessedHistograms();

    // Create a new systematic uncertainty organizer
    uncertainties[iDataType] = new SystematicUncertaintyOrganizer(uncertaintyFile[iDataType]);

  }
 
  TH1D* energyEnergyCorrelatorSignalPbPb[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorSignalPp[nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorPbPbToPpRatio[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyForPbPb[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyForPp[nJetPtBinsEEC][nTrackPtBinsEEC];

  // Initialize histograms to NULL
  for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
      energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt] = NULL;
      systematicUncertaintyForPp[iJetPt][iTrackPt] = NULL;
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt] = NULL;
        energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt] = NULL;
        systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt] = NULL;
      } // Centrality loop
    } // Track pT loop
  } // Jet pT loop

  // Read the histograms from managers
  int iTrackPtMatchedPp, iTrackPtMatchedPbPbUncertainty, iTrackPtMatchedPpUncertainty;
  int iJetPtMatchedPp, iJetPtMatchedPbPbUncertainty, iJetPtMatchedPpUncertainty;
  int iCentralityMatched;

  for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
    iJetPtMatchedPp = card[kPp]->FindBinIndexJetPtEEC(card[kPbPb]->GetBinBordersJetPtEEC(iJetPt));
    iJetPtMatchedPbPbUncertainty = uncertaintyCard[kPbPb]->FindBinIndexJetPtEEC(card[kPbPb]->GetBinBordersJetPtEEC(iJetPt));
    iJetPtMatchedPpUncertainty = uncertaintyCard[kPp]->FindBinIndexJetPtEEC(card[kPbPb]->GetBinBordersJetPtEEC(iJetPt));
    for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
      iTrackPtMatchedPp = card[kPp]->FindBinIndexTrackPtEEC(card[kPbPb]->GetBinBordersTrackPtEEC(iTrackPt));
      iTrackPtMatchedPbPbUncertainty = uncertaintyCard[kPbPb]->FindBinIndexTrackPtEEC(card[kPbPb]->GetBinBordersTrackPtEEC(iTrackPt));
      iTrackPtMatchedPpUncertainty = uncertaintyCard[kPp]->FindBinIndexTrackPtEEC(card[kPbPb]->GetBinBordersTrackPtEEC(iTrackPt));

      // Read the pp histograms that do not have centrality binning
      energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt] = histograms[kPp]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPtMatchedPp, iTrackPtMatchedPp, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);
      systematicUncertaintyForPp[iJetPt][iTrackPt] = uncertainties[kPp]->GetSystematicUncertainty(0, iJetPtMatchedPpUncertainty, iTrackPtMatchedPpUncertainty);

      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        iCentralityMatched = uncertaintyCard[kPbPb]->FindBinIndexCentrality(card[kPbPb]->GetBinBordersCentrality(iCentrality));

        // Read the PbPb histograms
        energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt] = histograms[kPbPb]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal);
        systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt] = uncertainties[kPbPb]->GetSystematicUncertainty(iCentralityMatched, iJetPtMatchedPbPbUncertainty, iTrackPtMatchedPbPbUncertainty);
      } // Centrality loop
    } // Track pT loop
  } // Jet pT loop

  // ================================================ //
  //   Normalize the histograms and calculate ratios  //
  // ================================================ //

  std::pair<double, double> analysisDeltaR = std::make_pair(0.006, 0.39); // DeltaR span in which the analysis is done
  double epsilon = 0.0001;
  int lowAnalysisBin = energyEnergyCorrelatorSignalPbPb[firstDrawnCentralityBin][firstDrawnJetPtBinEEC][firstDrawnTrackPtBinEEC]->GetXaxis()->FindBin(analysisDeltaR.first + epsilon);
  int highAnalysisBin = energyEnergyCorrelatorSignalPbPb[firstDrawnCentralityBin][firstDrawnJetPtBinEEC][firstDrawnTrackPtBinEEC]->GetXaxis()->FindBin(analysisDeltaR.second - epsilon);
  for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
    for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
      energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->Scale(1.0 / energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
      systematicUncertaintyForPp[iJetPt][iTrackPt]->Scale(1.0 / systematicUncertaintyForPp[iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->Scale(1.0 / energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));
        systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->Scale(1.0 / systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->Integral(lowAnalysisBin, highAnalysisBin, "width"));

        energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt] = (TH1D*) energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->Clone(Form("energyEnergyCorrelatorRatio%d%d%d", iCentrality, iJetPt, iTrackPt));
        energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->Divide(energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]);
      } // Centrality loop
    } // Track pT loop
  } // Jet pT loop

  // ========================================= //
  //   Draw the histograms in separate plots   //
  // ========================================= //

  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  drawer->SetLogX(true); // Logarithmic deltaR axis

  //TString centralityString = Form("Cent: %.0f-%.0f%%", card->GetLowBinBorderCentrality(iCentrality), card->GetHighBinBorderCentrality(iCentrality));
  //TString compactCentralityString = Form("_C=%.0f-%.0f", card->GetLowBinBorderCentrality(iCentrality), card->GetHighBinBorderCentrality(iCentrality));
  TString jetPtString;
  TString compactJetPtString;
  TString trackPtString;
  TString compactTrackPtString;

  for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
    jetPtString = Form("%.0f < jet p_{T} < %.0f", card[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));
    compactJetPtString = Form("_J=%.0f-%.0f", card[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));
    for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
      trackPtString = Form("%.1f < track p_{T}",card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));
      compactTrackPtString = Form("_T>%.1f",card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));

      // Create a new canvas for the plot
      drawer->CreateSplitCanvas();

      // Use logarithmic axis for EEC
      drawer->SetLogY(true);

      // Setup the legend for plots
      TLegend *legend = new TLegend(0.23,0.05,0.53,0.6);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
      legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
      legend->AddEntry((TObject*) 0, trackPtString.Data(),"");

      // Set the drawing style for pp histogram
      energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->SetMarkerStyle(kFullDiamond);
      energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->SetMarkerColor(kBlack);

      // Set the x-axis drawing range
      energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);

      // Draw the pp correlator to upper canves
      drawer->DrawHistogramToUpperPad(energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt], "#Deltar", "EEC Signal", " ", "p");
      legend->AddEntry(energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt], "pp", "p");

      // Draw the different centrality bins to the same plot
      for(int iCentrality = lastDrawnCentralityBin; iCentrality >= firstDrawnCentralityBin; iCentrality--){
        energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColor[iCentrality - firstDrawnCentralityBin]);
        energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iCentrality - firstDrawnCentralityBin]);
        energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->Draw("same,p");
        legend->AddEntry(energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt], Form("PbPb %.0f-%.0f%%", card[kPbPb]->GetLowBinBorderCentrality(iCentrality), card[kPbPb]->GetHighBinBorderCentrality(iCentrality)), "p");
      }

      // Draw the legend to the upper pad
      legend->Draw();

      // Linear scale for the ratio
      drawer->SetLogY(false);

      // Set the axis drawing ranges
      energyEnergyCorrelatorPbPbToPpRatio[lastDrawnCentralityBin][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
      energyEnergyCorrelatorPbPbToPpRatio[lastDrawnCentralityBin][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

      // Set the style for histograms
      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColor[iCentrality - firstDrawnCentralityBin]);
        energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iCentrality - firstDrawnCentralityBin]);
      }

      drawer->SetGridY(true);
      drawer->DrawHistogramToLowerPad(energyEnergyCorrelatorPbPbToPpRatio[lastDrawnCentralityBin][iJetPt][iTrackPt], "#Deltar", "#frac{PbPb}{pp}", " ");
      for(int iCentrality = lastDrawnCentralityBin-1; iCentrality >= firstDrawnCentralityBin; iCentrality--) {
        energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->Draw("same,p");
      }
      drawer->SetGridY(false);
  
      // If a plot name is given, save the plot in a file
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/energyEnergyCorrelator_centralityComparison_%s%s%s.pdf", saveComment.Data(), compactJetPtString.Data(), compactTrackPtString.Data()));
      }


    } // Track pT loop
  } // Jet pT loop

  // // Open input file for reading and read graphs from it
  // TFile *inputFile = TFile::Open(directoryName+inputFileName);
  // for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
  //   jetVnGraph[iFlow] = (TGraphErrors*) inputFile->Get(Form("summaryV%d", iFlow+1));
    
  //   // If the graph we wanted to load does not exist, inform the user and end program
  //   if(jetVnGraph[iFlow] == NULL){
  //     cout << "Hey dude! The file: " << directoryName.Data() << inputFileName.Data() << " does not contain graph: " << Form("summaryV%d", iFlow+1) <<  "." << endl;
  //     cout << "Cannot do the plotting, mate! Be a good lad and make sure the graph is there next time, ok?" << endl;
  //     return;
  //   }
    
  //   // Set the style for the jet vn values
  //   jetVnGraph[iFlow]->SetMarkerStyle(bigCanvasMarker[iFlow]);
  //   jetVnGraph[iFlow]->SetMarkerColor(bigCanvasColor[iFlow]);
  //   jetVnGraph[iFlow]->SetMarkerSize(2);
    
  // }
  
  // // Read the systematic uncertainty graphs
  // TFile *uncertaintyFile = TFile::Open(directoryName+uncertaintyFileName);
  // LongRangeSystematicOrganizer *uncertaintyOrganizer = new LongRangeSystematicOrganizer(uncertaintyFile);
  // uncertaintyOrganizer->AdjustCentralPoints(jetVnGraph);
  
  // for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
  //   jetVnUncertainty[iFlow] = uncertaintyOrganizer->GetLongRangeSystematicUncertainty(iFlow);
    
  //   // If the graph we wanted to load does not exist, inform the user and end program
  //   if(jetVnUncertainty[iFlow] == NULL){
  //     cout << "Hey dude! The file: " << directoryName.Data() << uncertaintyFileName.Data() << " does not contain uncertainties for " << Form("v%d", iFlow+1) <<  "." << endl;
  //     cout << "No uncertainties can be plotted." << endl;
  //   }
    
  //   // Set the style for the jet vn uncertainties
  //   jetVnUncertainty[iFlow]->SetMarkerStyle(bigCanvasMarker[iFlow]);
  //   jetVnUncertainty[iFlow]->SetLineColor(bigCanvasColor[iFlow]);
  //   jetVnUncertainty[iFlow]->SetFillColorAlpha(bigCanvasColor[iFlow], 0.3);
  //   jetVnUncertainty[iFlow]->SetMarkerColor(bigCanvasColor[iFlow]);
  //   jetVnUncertainty[iFlow]->SetMarkerSize(2);
    
  // }
  
  // // Set the bin labels for x-axis
  // TString binLabels[] = {"0-10%"," ","10-30%"," ","30-50%"," ","50-90%"};
  // for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
  //   for(int iCentrality = 0; iCentrality < jetVnUncertainty[iFlow]->GetN()*2; iCentrality++){
  //     jetVnUncertainty[iFlow]->GetXaxis()->ChangeLabel(iCentrality+1,-1,-1,-1,-1,-1,binLabels[iCentrality]);
  //   } // Centrality loop
  // }
  
  // // ===================== //
  // // Previous measurements //
  // // ===================== //
  
  // // Previous results that can be plotted together with data from this analysis
  // double summaryXaxis[nCentralityBins];
  // double summaryXaxisError[nCentralityBins];
  // double summaryYaxisError[nCentralityBins];
  // for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
  //   summaryXaxis[iCentrality] = iCentrality+1;
  //   summaryXaxisError[iCentrality] = 0;
  //   summaryYaxisError[iCentrality] = 0;
  // }
  
  // // ATLAS results taken from those presented in Hard Probes 2020: ATLAS-CONF-2020-019 (https://cds.cern.ch/record/2720249?ln=en)
  // const double atlasV2Number[] = {0.018, 0.03, 0.035, 0.03};
  // TGraphErrors* atlasJetV2graph = new TGraphErrors(nCentralityBins, summaryXaxis, atlasV2Number, summaryXaxisError, summaryXaxisError);
  // atlasJetV2graph->SetMarkerStyle(kFullDiamond);
  // atlasJetV2graph->SetMarkerColor(kViolet-2);
  // atlasJetV2graph->SetMarkerSize(1.8);
  
  // // CMS high pT results extracted from analysis arXiv:1702.00630 (PLB 776 (2018) 195)
  // const double cmsHighPtV2Number[] = {0.0220, 0.0376, 0.0431, 0.04};
  // const double cmsHighPtV2Error[] = {0.0019, 0.0016, 0.0027, 0.04};
  // TGraphErrors* cmsHighPtV2 = new TGraphErrors(nCentralityBins, summaryXaxis, cmsHighPtV2Number, summaryXaxisError, cmsHighPtV2Error);
  // cmsHighPtV2->SetMarkerStyle(kFullStar);
  // cmsHighPtV2->SetMarkerColor(kAzure+9);
  // cmsHighPtV2->SetMarkerSize(1.8);
  
  // // ============== //
  // // Draw the plots //
  // // ============== //
  
  // // Setup the drawer for graphs
  // JDrawer *drawer = new JDrawer();
  // drawer->SetDefaultAppearanceGraph();
  // drawer->SetNDivisionsX(510);
  // drawer->SetNDivisionsY(510);
  // drawer->SetBottomMargin(0.16);
  // drawer->SetTitleOffsetX(1.3);
  // drawer->SetLabelOffsetX(0.02);
  // drawer->SetTitleOffsetY(1.77);
  // drawer->SetLabelOffsetY(0.01);
  
  // // Draw the graphs for selected flow components
  // TLegend *legend;
  // TLegend *anotherLegend;
  // double errorY;
  // double cmsYPosition, cmsXPosition;
  
  // TLine *zeroLine = new TLine(0.75,0,3.35,0);
  // zeroLine->SetLineStyle(2);
  
  // TLatex *preliminaryText = new TLatex();
  
  // // =======================================================================
  // // == Drawing style, where each histogram is drawn into it's own canvas ==
  // // =======================================================================
  
  // for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
  //   legend = new TLegend(0.2,0.7,0.5,0.9);
    
  //   // Create legends for the plot
  //   if(iFlow == 1){
  //     legend = new TLegend(0.2,0.83,0.7,0.92);
  //   } else {
  //     legend = new TLegend(0.21,0.82,0.71,0.91);
  //   }
  //   anotherLegend = new TLegend(0.2,0.71,0.7,0.83);
    
  //   legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  //   anotherLegend->SetFillStyle(0);anotherLegend->SetBorderSize(0);anotherLegend->SetTextSize(0.05);anotherLegend->SetTextFont(62);
    
  //   if(jetVnUncertainty[iFlow] != NULL){
      
  //     // Set the style for uncertainties
  //     for(int iCentrality; iCentrality < jetVnUncertainty[iFlow]->GetN(); iCentrality++){
  //       errorY = jetVnUncertainty[iFlow]->GetErrorY(iCentrality);
  //       jetVnUncertainty[iFlow]->SetPointError(iCentrality, 0.1, errorY);
  //     }
      
  //     //jetVnUncertainty[iFlow]->GetYaxis()->SetNdivisions(510);
  //     drawer->DrawGraphCustomAxes(jetVnUncertainty[iFlow], 0, 4, minZoom[iFlow], maxZoom[iFlow], "Centrality", "Dijet v_{n}", " ", "a,e2");
      
  //     jetVnGraph[iFlow]->Draw("p,same");
      
  //     legend->AddEntry(jetVnUncertainty[iFlow], Form("Dijet v_{%d}", iFlow+1), "pf");
      
  //   } else {
  //     drawer->DrawGraphCustomAxes(jetVnGraph[iFlow], 0, 4, minZoom[iFlow], maxZoom[iFlow], "Centrality", "Dijet v_{n}", " ", "ap");
  //     legend->AddEntry(jetVnGraph[iFlow], Form("Dijet v_{%d}", iFlow+1), "p");
  //   }
    
  //   zeroLine->Draw();
    
  //   if(iFlow == 1 && drawAtlasJetV2){
  //     atlasJetV2graph->Draw("p,same");
  //     legend->AddEntry(atlasJetV2graph, "ATLAS v_{2}", "p"); // (#scale[0.8]{HP 2020})
  //   }
    
  //   if(iFlow == 1 && drawCmsHigtPtV2){
  //     cmsHighPtV2->Draw("p,same");
  //     //legend->AddEntry(cmsHighPtV2, "CMS high p_{T} v_{2}", "p"); // (#scale[0.8]{PLB 776 (2018) 195})
  //     anotherLegend->AddEntry(cmsHighPtV2, "CMS charged hadron v_{2}", "p"); // (#scale[0.8]{PLB 776 (2018) 195})
  //     anotherLegend->AddEntry((TObject*)0, "p_{T} > 20 GeV, |#eta| < 1", "");
  //     anotherLegend->Draw();
  //   }
    
  //   legend->Draw();
    
  //   preliminaryText->SetTextFont(62);
  //   preliminaryText->SetTextSize(0.065);
  //   cmsYPosition = 0.32; cmsXPosition = 0.48;
  //   if(iFlow > 1){
  //     cmsYPosition = 0.83; cmsXPosition = 0.67;
  //   }
  //   preliminaryText->DrawLatexNDC(cmsXPosition, cmsYPosition, "CMS");
    
  //   preliminaryText->SetTextFont(42);
  //   preliminaryText->SetTextSize(0.055);
  //   preliminaryText->DrawLatexNDC(cmsXPosition-0.053, cmsYPosition-0.055, "Preliminary");
    
    
  //   // Save the figures to file
  //   if(saveFigures){
  //     gPad->GetCanvas()->SaveAs(Form("figures/finalJetV%d%s.pdf", iFlow+1, saveComment.Data()));
  //   }
  // }
  
  // // =============================================
  // // ===== Drawing style with one big canvas =====
  // // =============================================
  
  // if(drawBigCanvas){
    
  //   // Draw all the distributions to big canvases
  //   auxi_canvas *bigCanvas;
  //   TLatex *mainTitle;
    
  //   // Draw a big canvas and put all the plots in it
  //   bigCanvas = new auxi_canvas("bigCanvas", "", 1500, 550);
  //   bigCanvas->SetMargin(0.07, 0.01, 0.16, 0.01);
  //   bigCanvas->divide(1,3);
    
  //   mainTitle = new TLatex();
    
  //   for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
  //     bigCanvas->CD(iFlow);
      
      
  //     // Adjust the number of divisions
  //     jetVnUncertainty[iFlow]->GetYaxis()->SetNdivisions(510);
  //     jetVnUncertainty[iFlow]->GetYaxis()->SetLabelOffset(0.01);
  //     jetVnUncertainty[iFlow]->GetYaxis()->SetTitle("Dijet v_{n}");
  //     jetVnUncertainty[iFlow]->GetYaxis()->SetTitleOffset(1.5);
  //     jetVnUncertainty[iFlow]->GetYaxis()->SetTitleSize(0.06);
      
  //     // Adjust label sizes for all but the leftmost plot
  //     if(iFlow > firstDrawnVn-1){
  //       jetVnUncertainty[iFlow]->GetXaxis()->SetTitleSize(0.069);
  //       jetVnUncertainty[iFlow]->GetXaxis()->SetTitleOffset(1.09);
  //       jetVnUncertainty[iFlow]->GetXaxis()->SetLabelSize(0.058);
  //       jetVnUncertainty[iFlow]->GetXaxis()->SetLabelOffset(0.0135);
  //     } else {
  //       jetVnUncertainty[iFlow]->GetXaxis()->SetTitleOffset(1.25);
  //       jetVnUncertainty[iFlow]->GetXaxis()->SetLabelOffset(0.02);
  //     }
      
  //     // Draw the graphs
  //     jetVnUncertainty[iFlow]->Draw("a,e2");
  //     jetVnGraph[iFlow]->Draw("same,p");
      
  //     // Create a legend for the plot division
  //     if(iFlow > 1) {
  //       legend = new TLegend(0.06,0.89,0.7,0.98);
  //       anotherLegend = new TLegend(0.23,0.7,0.73,0.85);
  //     } else {
  //       legend = new TLegend(0.23,0.89,0.73,0.98);
  //       anotherLegend = new TLegend(0.23,0.77,0.73,0.89);
  //     }
      
  //     legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.054);legend->SetTextFont(62);
  //     anotherLegend->SetFillStyle(0);anotherLegend->SetBorderSize(0);anotherLegend->SetTextSize(0.054);anotherLegend->SetTextFont(62);
      
  //     if(iFlow > 1) legend->SetTextSize(0.06);  // Need to increase text size for smaller divisions
      
  //     legend->AddEntry(jetVnUncertainty[iFlow], Form("Dijet v_{%d}", iFlow+1), "pf");
      
  //     if(iFlow == 1 && drawAtlasJetV2){
  //       atlasJetV2graph->SetMarkerSize(2);
  //       atlasJetV2graph->Draw("p,same");
  //       legend->AddEntry(atlasJetV2graph, "ATLAS v_{2}", "p"); // (#scale[0.8]{HP 2020})
  //     }
      
  //     if(iFlow == 1 && drawCmsHigtPtV2){
  //       cmsHighPtV2->SetMarkerSize(2);
  //       cmsHighPtV2->Draw("p,same");
  //       anotherLegend->AddEntry(cmsHighPtV2, "CMS charged hadron v_{2}", "p"); // (#scale[0.8]{PLB 776 (2018) 195})
  //       anotherLegend->AddEntry((TObject*)0, "p_{T} > 20 GeV, |#eta| < 1", "");
  //     }
      
  //     legend->Draw();
      
  //     zeroLine->Draw();
      
  //     if(iFlow == 1) anotherLegend->Draw();
      
  //   } // Flow component loop
    
  //   // Draw all necessary CMS text to the plot
  //   bigCanvas->cd(0);
    
  //   mainTitle->SetTextFont(62);
  //   mainTitle->SetTextSize(0.065);
  //   cmsYPosition = 0.89;
  //   if(drawPreliminaryTag) cmsYPosition = 0.91;
  //   mainTitle->DrawLatexNDC(0.9, cmsYPosition, "CMS");
    
  //   mainTitle->SetTextFont(42);
  //   mainTitle->SetTextSize(0.055);
  //   if(drawPreliminaryTag)mainTitle->DrawLatexNDC(0.878, 0.86, "Preliminary");
  //   mainTitle->DrawLatexNDC(0.72, 0.76, "PbPb #sqrt{s_{NN}} = 5.02 TeV, 1.7 nb^{-1}");
  //   mainTitle->DrawLatexNDC(0.54, 0.92, "anti-k_{T} R = 0.4");
  //   mainTitle->DrawLatexNDC(0.54, 0.84, "|#eta_{jet}| < 1.3");
  //   if(leadSubTitles){
  //   mainTitle->DrawLatexNDC(0.54, 0.76, "p_{T}^{lead} > 120 GeV");
  //   mainTitle->DrawLatexNDC(0.54, 0.68, "p_{T}^{sub} > 50 GeV");
  //   mainTitle->DrawLatexNDC(0.54, 0.60, "#Delta#varphi > #frac{5#pi}{6}");
  //   } else {
  //     mainTitle->DrawLatexNDC(0.54, 0.76, "p_{T,1} > 120 GeV");
  //     mainTitle->DrawLatexNDC(0.54, 0.68, "p_{T,2} > 50 GeV");
  //     mainTitle->DrawLatexNDC(0.54, 0.60, "#Delta#varphi_{1,2} > #frac{5#pi}{6}");
  //   }
  //   mainTitle->DrawLatexNDC(0.13, 0.33, "Factorization region:");
  //   mainTitle->DrawLatexNDC(0.13, 0.25, "0.7 < Hadron p_{T} < 3 GeV");
    
  //   // Save the figures to file
  //   if(saveFigures){
  //     gPad->GetCanvas()->SaveAs(Form("figures/finalBigCanvasJetVn%s.pdf", saveComment.Data()));
  //   }
    
  // }
  
  // // Print the central values of the points to the console
  // if(printCentralValues){
    
  //   double vnValue;
  //   int centralityBinBorders[] = {0, 10, 30, 50, 90};
    
  //   for(int iFlow = firstDrawnVn-1; iFlow <= lastDrawnVn-1; iFlow++){
      
  //     for(int iPoint = 0; iPoint < 3; iPoint++){
  //       vnValue = jetVnGraph[iFlow]->GetPointY(iPoint);
  //       cout << "Dijet v" << iFlow+1 << " C: " << centralityBinBorders[iPoint] << "-" << centralityBinBorders[iPoint+1] << " = " << vnValue << endl;
  //     }
      
  //   }
  // }
}
