#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Macro for checking DeltaR and pT1*pT2 resolution between track and matched particle pairs
 */
void trackParticleMatchingPlotter(){

  // ============= //
  // Configuration //
  // ============= //
  
  // Open the input file
  TString inputFileName = "data/PbPbMC2018_RecoGen_akFlowJets_miniAOD_4pCentShift_noTrigger_nominalSmear_onlyTrackParticleMatching_processed_2023-07-03.root";
  // PbPbMC2018_RecoGen_akFlowJets_miniAOD_4pCentShift_noTrigger_nominalSmear_onlyTrackParticleMatching_processed_2023-07-03.root
  // ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_32deltaRBins_nominalSmear_onlyTrackParticleMatching_processed_2023-07-03.root
  TFile* inputFile = TFile::Open(inputFileName);

  // Check that the input file exists
  if(inputFile == NULL){
    cout << "Error! The file " << inputFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }

  // Load the card from the file
  EECCard* card = new EECCard(inputFile);
  TString collisionSystem = card->GetDataType();
  bool isPbPbData = collisionSystem.Contains("PbPb");

  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Find the number of bins from the Pythia+Hydjet card
  const int nCentralityBins = card->GetNCentralityBins();
  const int nJetPtBinsEEC = card->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card->GetNTrackPtBinsEEC();
  
  // Selection of drawn bins
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = nCentralityBins-1;
  
  int firstDrawnJetPtBinEEC = 6;
  int lastDrawnJetPtBinEEC = 6;
  
  int firstDrawnTrackPtBinEEC = 5;
  int lastDrawnTrackPtBinEEC = 5;
  
  // Save the final plots
  const bool saveFigures = false;
  TString saveComment = "_firstLook";

  // Constant: 0.771 pp -> 0.801 PbPb
  // Mean = 1
  // Sigma = 0.0244 pp -> 0.0237 PbPb

  // Ratio zoom settings
  std::pair<double, double> analysisDeltaR = std::make_pair(0.006, 0.39); // DeltaR span in which the analysis is done

  // ============================================== //
  // Read the histograms from the histogram manager //
  // ============================================== //

  EECHistogramManager* histograms = new EECHistogramManager(inputFile, card);

  // Load the histograms related to track/particle matching study
  histograms->SetLoadTrackParticleMatchingHistograms(true);
  histograms->SetCentralityBinRange(firstDrawnCentralityBin, lastDrawnCentralityBin);
  histograms->SetTrackPtBinRangeEEC(firstDrawnTrackPtBinEEC, lastDrawnTrackPtBinEEC);
  histograms->SetJetPtBinRangeEEC(firstDrawnJetPtBinEEC, lastDrawnJetPtBinEEC);

  // Load the histograms from the file
  histograms->LoadProcessedHistograms();
  
 
  // Define the histograms to be drawn
  TH2D* trackParticleDeltaRResponseMatrix[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH2D* trackParticlePtResponseMatrix[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* trackParticlePtClosure[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];

  // Initialize histograms to NULL
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        trackParticleDeltaRResponseMatrix[iCentrality][iJetPt][iTrackPt] = NULL;
        trackParticlePtResponseMatrix[iCentrality][iJetPt][iTrackPt] = NULL;
        trackParticlePtClosure[iCentrality][iJetPt][iTrackPt] = NULL;
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

  // Read the histograms from managers

  for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
    for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
      for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){

        trackParticleDeltaRResponseMatrix[iCentrality][iJetPt][iTrackPt] = histograms->GetHistogramTrackParticleDeltaRResponse(iCentrality, iJetPt, iTrackPt);
        trackParticlePtResponseMatrix[iCentrality][iJetPt][iTrackPt] = histograms->GetHistogramTrackParticlePtResponse(iCentrality, iJetPt, iTrackPt);
        trackParticlePtClosure[iCentrality][iJetPt][iTrackPt] = histograms->GetHistogramTrackParticlePtClosure(iCentrality, iJetPt, iTrackPt);

      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

  // ================================================ //
  //            Normalize the histograms              //
  // ================================================ //

  double rowContent;

  for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
    for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
      for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){

        // For closure histograms, just normalize the integral to one
        trackParticlePtClosure[iCentrality][iJetPt][iTrackPt]->Scale(1.0 / trackParticlePtClosure[iCentrality][iJetPt][iTrackPt]->Integral());

        // For the two-dimensional histograms, normalize each row to one
        for(int iRow = 1; iRow <= trackParticleDeltaRResponseMatrix[iCentrality][iJetPt][iTrackPt]->GetNbinsX(); iRow++){
          rowContent = 0;
          for(int iColumn = 1; iColumn <= trackParticleDeltaRResponseMatrix[iCentrality][iJetPt][iTrackPt]->GetNbinsY(); iColumn++){
            rowContent += trackParticleDeltaRResponseMatrix[iCentrality][iJetPt][iTrackPt]->GetBinContent(iRow, iColumn);
          } // Loop over columns in response matrix
          for(int iColumn = 1; iColumn <= trackParticleDeltaRResponseMatrix[iCentrality][iJetPt][iTrackPt]->GetNbinsY(); iColumn++){
            trackParticleDeltaRResponseMatrix[iCentrality][iJetPt][iTrackPt]->SetBinContent(iRow, iColumn, trackParticleDeltaRResponseMatrix[iCentrality][iJetPt][iTrackPt]->GetBinContent(iRow, iColumn) / rowContent);
            trackParticleDeltaRResponseMatrix[iCentrality][iJetPt][iTrackPt]->SetBinError(iRow, iColumn, trackParticleDeltaRResponseMatrix[iCentrality][iJetPt][iTrackPt]->GetBinError(iRow, iColumn) / rowContent);
          } // Loop over columns in DeltaR response matrix
        } // Loop over rows in DeltaR response matrix

        for(int iRow = 1; iRow <= trackParticlePtResponseMatrix[iCentrality][iJetPt][iTrackPt]->GetNbinsX(); iRow++){
          rowContent = 0;
          for(int iColumn = 1; iColumn <= trackParticlePtResponseMatrix[iCentrality][iJetPt][iTrackPt]->GetNbinsY(); iColumn++){
            rowContent += trackParticlePtResponseMatrix[iCentrality][iJetPt][iTrackPt]->GetBinContent(iRow, iColumn);
          } // Loop over columns in response matrix
          for(int iColumn = 1; iColumn <= trackParticlePtResponseMatrix[iCentrality][iJetPt][iTrackPt]->GetNbinsY(); iColumn++){
            trackParticlePtResponseMatrix[iCentrality][iJetPt][iTrackPt]->SetBinContent(iRow, iColumn, trackParticlePtResponseMatrix[iCentrality][iJetPt][iTrackPt]->GetBinContent(iRow, iColumn) / rowContent);
            trackParticlePtResponseMatrix[iCentrality][iJetPt][iTrackPt]->SetBinError(iRow, iColumn, trackParticlePtResponseMatrix[iCentrality][iJetPt][iTrackPt]->GetBinError(iRow, iColumn) / rowContent);
          } // Loop over columns in response matrix
        } // Loop over rows in response matrix

      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

  // ========================================= //
  //   Draw the histograms in separate plots   //
  // ========================================= //

  JDrawer* drawer = new JDrawer();
  drawer->SetRelativeCanvasSize(1.1, 1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);

  TString centralityString;
  TString compactCentralityString;
  TString jetPtString;
  TString compactJetPtString;
  TString trackPtString;
  TString compactTrackPtString;

  // Common variables for different plots
  TLegend* legend;

  // Draw plots for particle pair pT closure
  for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++) {
    if(isPbPbData){
      centralityString = Form("Cent: %.0f-%.0f%%", card->GetLowBinBorderCentrality(iCentrality), card->GetHighBinBorderCentrality(iCentrality));
      compactCentralityString = Form("_C=%.0f-%.0f", card->GetLowBinBorderCentrality(iCentrality), card->GetHighBinBorderCentrality(iCentrality));
    } else {
      centralityString = "pp";
      compactCentralityString = "_pp";
    }

    for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++) {
      jetPtString = Form("%.0f < jet p_{T} < %.0f", card->GetLowBinBorderJetPtEEC(iJetPt), card->GetHighBinBorderJetPtEEC(iJetPt));
      compactJetPtString = Form("_J=%.0f-%.0f", card->GetLowBinBorderJetPtEEC(iJetPt), card->GetHighBinBorderJetPtEEC(iJetPt));

      for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++) {
        trackPtString = Form("%.1f < track p_{T}", card->GetLowBinBorderTrackPtEEC(iTrackPt));
        compactTrackPtString = Form("_T>%.1f", card->GetLowBinBorderTrackPtEEC(iTrackPt));
        compactTrackPtString.ReplaceAll(".", "v");

        // Setup the legend for plots
        legend = new TLegend(0.53, 0.6, 0.83, 0.85);
        legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);
        legend->AddEntry((TObject*)0, centralityString.Data(), "");
        legend->AddEntry((TObject*)0, jetPtString.Data(), "");
        legend->AddEntry((TObject*)0, trackPtString.Data(), "");

        // Draw the histograms
        trackParticlePtClosure[iCentrality][iJetPt][iTrackPt]->Fit("gaus");
        trackParticlePtClosure[iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.7,1.3);
        drawer->DrawHistogram(trackParticlePtClosure[iCentrality][iJetPt][iTrackPt], "p_{T,1} p_{T,2} (Reco / Gen)", "Probability", " ");

        // Draw the legend
        legend->Draw();

        // If a plot name is given, save the plot in a file
        if(saveFigures) {
          gPad->GetCanvas()->SaveAs(Form("figures/trackParticlePairPtClosure%s%s%s%s.pdf", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data()));
        }

      } // Track pT loop
    } // Jet pT loop   
  } // Centrality loop

  drawer->SetLeftMargin(0.13);
  drawer->SetRightMargin(0.11);
  drawer->SetTopMargin(0.08);
  drawer->SetTitleOffsetX(1.17);
  drawer->SetTitleOffsetY(1);
  drawer->SetLogY(true);
  drawer->SetLogX(true);
  gStyle->SetPaintTextFormat("0.2f");

  // Draw the DeltaR response matrices
  for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++) {
    if(isPbPbData){
      centralityString = Form("Cent: %.0f-%.0f%%", card->GetLowBinBorderCentrality(iCentrality), card->GetHighBinBorderCentrality(iCentrality));
      compactCentralityString = Form("_C=%.0f-%.0f", card->GetLowBinBorderCentrality(iCentrality), card->GetHighBinBorderCentrality(iCentrality));
    } else {
      centralityString = "pp";
      compactCentralityString = "_pp";
    }

    for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++) {
      jetPtString = Form("%.0f < jet p_{T} < %.0f", card->GetLowBinBorderJetPtEEC(iJetPt), card->GetHighBinBorderJetPtEEC(iJetPt));
      compactJetPtString = Form("_J=%.0f-%.0f", card->GetLowBinBorderJetPtEEC(iJetPt), card->GetHighBinBorderJetPtEEC(iJetPt));

      for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++) {
        trackPtString = Form("%.1f < track p_{T}", card->GetLowBinBorderTrackPtEEC(iTrackPt));
        compactTrackPtString = Form("_T>%.1f", card->GetLowBinBorderTrackPtEEC(iTrackPt));
        compactTrackPtString.ReplaceAll(".", "v");

        // Setup the legend for plots
        legend = new TLegend(0.13, 0.65, 0.43, 0.9);
        legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);
        legend->AddEntry((TObject*)0, centralityString.Data(), "");
        legend->AddEntry((TObject*)0, jetPtString.Data(), "");
        legend->AddEntry((TObject*)0, trackPtString.Data(), "");

        // Draw the histograms
        trackParticleDeltaRResponseMatrix[iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.002,0.5);
        trackParticleDeltaRResponseMatrix[iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.002,0.5);
        drawer->DrawHistogram(trackParticleDeltaRResponseMatrix[iCentrality][iJetPt][iTrackPt], "#Deltar (Reco)", "#Deltar (Gen)", " ", "colz,text20");

        // Draw the legend
        legend->Draw();

        // If a plot name is given, save the plot in a file
        if(saveFigures) {
          gPad->GetCanvas()->SaveAs(Form("figures/trackParticleDeltaRResponseMatrix%s%s%s%s.pdf", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data()));
        }

      } // Track pT loop
    } // Jet pT loop   
  } // Centrality loop

  drawer->SetLogY(false);
  drawer->SetLogX(false);

  // Draw the pT1*pT2 response matrices
  for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++) {
    if(isPbPbData){
      centralityString = Form("Cent: %.0f-%.0f%%", card->GetLowBinBorderCentrality(iCentrality), card->GetHighBinBorderCentrality(iCentrality));
      compactCentralityString = Form("_C=%.0f-%.0f", card->GetLowBinBorderCentrality(iCentrality), card->GetHighBinBorderCentrality(iCentrality));
    } else {
      centralityString = "pp";
      compactCentralityString = "_pp";
    }

    for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++) {
      jetPtString = Form("%.0f < jet p_{T} < %.0f", card->GetLowBinBorderJetPtEEC(iJetPt), card->GetHighBinBorderJetPtEEC(iJetPt));
      compactJetPtString = Form("_J=%.0f-%.0f", card->GetLowBinBorderJetPtEEC(iJetPt), card->GetHighBinBorderJetPtEEC(iJetPt));

      for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++) {
        trackPtString = Form("%.1f < track p_{T}", card->GetLowBinBorderTrackPtEEC(iTrackPt));
        compactTrackPtString = Form("_T>%.1f", card->GetLowBinBorderTrackPtEEC(iTrackPt));
        compactTrackPtString.ReplaceAll(".", "v");

        // Setup the legend for plots
        legend = new TLegend(0.13, 0.65, 0.43, 0.9);
        legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);
        legend->AddEntry((TObject*)0, centralityString.Data(), "");
        legend->AddEntry((TObject*)0, jetPtString.Data(), "");
        legend->AddEntry((TObject*)0, trackPtString.Data(), "");

        // Draw the histograms
        trackParticlePtResponseMatrix[iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(4,100);
        trackParticlePtResponseMatrix[iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(4,100);
        drawer->DrawHistogram(trackParticlePtResponseMatrix[iCentrality][iJetPt][iTrackPt], "p_{T,1} p_{T,2} (Reco)", "p_{T,1} p_{T,2} (Gen)", " ", "colz,text20");

        // Draw the legend
        legend->Draw();

        // If a plot name is given, save the plot in a file
        if(saveFigures) {
          gPad->GetCanvas()->SaveAs(Form("figures/trackParticlePtResponseMatrix%s%s%s%s.pdf", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data()));
        }

      } // Track pT loop
    } // Jet pT loop   
  } // Centrality loop
}
