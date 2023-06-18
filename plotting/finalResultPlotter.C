#include "SystematicUncertaintyOrganizer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECHistogramManager.h"
#include "EECCard.h"
#include "JDrawer.h"
#include "SplitCanvas.h"

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

  // Choose which plots to draw
  bool drawIndividualPlotsAllCentralities = true;
  bool drawBigCanvasDistributions = false;
  bool drawCentralityRatios = false;
  
  // Save the final plots
  const bool saveFigures = false;
  TString saveComment = "_preliminary";

  // Ratio zoom settings
  std::pair<double, double> ratioZoom = std::make_pair(0.4, 1.6);
  
  // Marker colors and styles
  int markerStylePbPb[] = {kFullSquare, kFullCircle, kFullCross, kFullFourTrianglesPlus};
  int markerStylePp = kFullDiamond;
  int markerColorPbPb[] = {kRed, kBlue, kMagenta, kGreen+3};
  int markerColorPp = kBlack;

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
  TH1D* systematicUncertaintyPbPbToPpRatio[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];

  // Initialize histograms to NULL
  for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
      energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt] = NULL;
      systematicUncertaintyForPp[iJetPt][iTrackPt] = NULL;
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt] = NULL;
        energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt] = NULL;
        systematicUncertaintyForPbPb[iCentrality][][iTrackPt] = NULL;
        systematicUncertaintyPbPbToPpRatio[iCentrality][iJetPt][iTrackPt] = NULL;
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

        systematicUncertaintyPbPbToPpRatio[iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->Clone(Form("systematicUncertaintyRatio%d%d%d", iCentrality, iJetPt, iTrackPt));
        systematicUncertaintyPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->Divide(systematicUncertaintyForPp[iJetPt][iTrackPt]);
      } // Centrality loop
    } // Track pT loop
  } // Jet pT loop

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

  //TString centralityString = Form("Cent: %.0f-%.0f%%", card->GetLowBinBorderCentrality(iCentrality), card->GetHighBinBorderCentrality(iCentrality));
  //TString compactCentralityString = Form("_C=%.0f-%.0f", card->GetLowBinBorderCentrality(iCentrality), card->GetHighBinBorderCentrality(iCentrality));
  TString jetPtString;
  TString compactJetPtString;
  TString trackPtString;
  TString compactTrackPtString;

  // Common variables for different plots
  TLegend* legend;
  TLatex* mainTitle;
  int canvasIndex;
  double bottomRowScale, bottomPadMargin, leftPadMargin;
  double leftMarginAdder, bottomMarginAdder;
  double thisPadScale;

  // Draw individual plots with all centralities mixed together in a single figure
  if(drawIndividualPlotsAllCentralities){

    for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
      jetPtString = Form("%.0f < jet p_{T} < %.0f", card[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));
      compactJetPtString = Form("_J=%.0f-%.0f", card[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));

      for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
        trackPtString = Form("%.1f < track p_{T}", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));
        compactTrackPtString = Form("_T>%.1f", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));
        compactTrackPtString.ReplaceAll(".","v");

        // Create a new canvas for the plot
        drawer->CreateSplitCanvas();

        // Use logarithmic axis for EEC
        drawer->SetLogY(true);

        // Setup the legend for plots
        legend = new TLegend(0.23, 0.05, 0.53, 0.6);
        legend->SetFillStyle(0);
        legend->SetBorderSize(0);
        legend->SetTextSize(0.05);
        legend->SetTextFont(62);
        legend->AddEntry((TObject*)0, jetPtString.Data(), "");
        legend->AddEntry((TObject*)0, trackPtString.Data(), "");

        // Set the drawing style for pp histogram
        energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->SetMarkerStyle(kFullDiamond);
        energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->SetMarkerColor(kBlack);

        // Set the x-axis drawing range
        energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);

        // Draw the pp correlator to upper canves
        drawer->DrawHistogramToUpperPad(energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt], "#Deltar", "EEC Signal", " ", "p");
        legend->AddEntry(energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt], "pp", "p");

        // Draw the different centrality bins to the same plot
        for(int iCentrality = lastDrawnCentralityBin; iCentrality >= firstDrawnCentralityBin; iCentrality--) {
          energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
          energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin]);
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
        for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++) {
          energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
          energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin]);
        }

        drawer->SetGridY(true);
        drawer->DrawHistogramToLowerPad(energyEnergyCorrelatorPbPbToPpRatio[lastDrawnCentralityBin][iJetPt][iTrackPt], "#Deltar", "#frac{PbPb}{pp}", " ");
        for(int iCentrality = lastDrawnCentralityBin - 1; iCentrality >= firstDrawnCentralityBin; iCentrality--) {
          energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->Draw("same,p");
        }
        drawer->SetGridY(false);

        // If a plot name is given, save the plot in a file
        if(saveFigures) {
          gPad->GetCanvas()->SaveAs(Form("figures/energyEnergyCorrelator_centralityComparison%s%s%s.pdf", saveComment.Data(), compactJetPtString.Data(), compactTrackPtString.Data()));
        }

      } // Track pT loop
    } // Jet pT loop
  } // Drawing individual canvases with all centralities in one canvas

  // =============================================
  // ===== Drawing style with one big canvas =====
  // =============================================
  
  if(drawBigCanvasDistributions){
    
    // Draw all the distributions to big canvases
    SplitCanvas* bigCanvas[nTrackPtBinsEEC];
    mainTitle = new TLatex();

    for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
      compactTrackPtString = Form("_T>%.1f", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));
      compactTrackPtString.ReplaceAll(".","v");

      // Draw a big canvas and put all the plots in it
      bigCanvas[iTrackPt] = new SplitCanvas(Form("bigCanvas%d", iTrackPt), "", 1800, 1500);
      bigCanvas[iTrackPt]->SetMargin(0.08, 0.01, 0.08, 0.11);
      bigCanvas[iTrackPt]->DivideNeatly(4, 4);

      bottomRowScale = bigCanvas[iTrackPt]->GetBottomRowScale();
      bottomPadMargin = bigCanvas[iTrackPt]->GetBottomPadMargin();
      leftPadMargin = bigCanvas[iTrackPt]->GetLeftPadMargin();

      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++) {
          jetPtString = Form("%.0f < jet p_{T} < %.0f", card[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));

          canvasIndex = (iCentrality - firstDrawnCentralityBin) * (lastDrawnJetPtBinEEC - firstDrawnJetPtBinEEC + 1) + (iJetPt - firstDrawnJetPtBinEEC);
          bigCanvas[iTrackPt]->CD(canvasIndex);

          gPad->SetLogy();
          gPad->SetLogx();

          // The titles in the bottow row need to be scaled, since it has more margin than other rows
          thisPadScale = (iCentrality == lastDrawnCentralityBin) ? bottomRowScale : 1;

          // Set the axis titles and labels
          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleOffset(0.85);
          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitleSize(0.2 * thisPadScale);
          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelOffset(0.01);
          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetLabelSize(0.12 * thisPadScale);

          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleOffset(0.5 / thisPadScale);
          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitleSize(0.2 * thisPadScale);
          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelOffset(0.01 / thisPadScale);
          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetLabelSize(0.12 * thisPadScale);

          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetTitle("EEC");
          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->SetStats(0);
          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetNdivisions(505);
          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->GetYaxis()->CenterTitle();
          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetTitle("#Deltar");
          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetNdivisions(505);

          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->SetTitle("");

          // Set the drawing ranges
          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.39);
          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.15, 30);

          // Set the drawing style for histograms
          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin]);
          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);
          energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin]);
          energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);

          systematicUncertaintyForPp[iJetPt][iTrackPt]->SetMarkerStyle(markerStylePp);
          systematicUncertaintyForPp[iJetPt][iTrackPt]->SetMarkerSize(1.2);
          energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->SetMarkerStyle(markerStylePp);
          energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->SetMarkerSize(1.2);

          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin], 0.4);
          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
          energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
          energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);

          systematicUncertaintyForPp[iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPp, 0.4);
          systematicUncertaintyForPp[iJetPt][iTrackPt]->SetLineColor(markerColorPp);
          energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->SetLineColor(markerColorPp);
          systematicUncertaintyForPp[iJetPt][iTrackPt]->SetMarkerColor(markerColorPp);
          energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->SetMarkerColor(markerColorPp);

          // Draw the histograms
          systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt]->Draw("e2");
          systematicUncertaintyForPp[iJetPt][iTrackPt]->Draw("same,e2");
          energyEnergyCorrelatorSignalPp[iJetPt][iTrackPt]->Draw("same,p");
          energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][iTrackPt]->Draw("same,p");

          // Create a legend for each pad
          leftMarginAdder = (iJetPt == firstDrawnJetPtBinEEC) ? leftPadMargin : 0;
          bottomMarginAdder = (iCentrality == lastDrawnCentralityBin) ? bottomPadMargin : 0;
          legend = new TLegend(0.06 + leftMarginAdder, 0.05 + bottomMarginAdder, 0.5 / (1 - leftMarginAdder), 0.45 / thisPadScale);
          legend->SetFillStyle(0);
          legend->SetBorderSize(0);
          legend->SetTextSize(0.09 * thisPadScale);
          legend->SetTextFont(62);
          legend->AddEntry((TObject*)0, jetPtString.Data(), "");
          legend->AddEntry(systematicUncertaintyForPbPb[iCentrality][iJetPt][iTrackPt], Form("PbPb %.0f-%.0f%%", card[kPbPb]->GetLowBinBorderCentrality(iCentrality), card[kPbPb]->GetHighBinBorderCentrality(iCentrality)), "p");
          legend->AddEntry(systematicUncertaintyForPp[iJetPt][iTrackPt], "pp", "p");

          legend->Draw();

        }  // Jet pT loop
      }    // Centrality loop

      // Draw all necessary CMS text to the plot
      bigCanvas[iTrackPt]->cd(0);

      mainTitle->SetTextFont(62);
      mainTitle->SetTextSize(0.065);
      mainTitle->DrawLatexNDC(0.08, 0.93, "CMS");

      mainTitle->SetTextFont(42);
      mainTitle->SetTextSize(0.025);
      mainTitle->DrawLatexNDC(0.23, 0.96, "PbPb #sqrt{s_{NN}} = 5.02 TeV, 1.70 nb^{-1}");
      mainTitle->DrawLatexNDC(0.23, 0.92, "pp #sqrt{s} = 5.02 TeV, 302 pb^{-1}");
      mainTitle->DrawLatexNDC(0.53, 0.96, "anti-k_{T} R = 0.4");
      mainTitle->DrawLatexNDC(0.53, 0.92, "|#eta_{jet}| < 1.6");

      mainTitle->SetTextSize(0.035);
      if(iTrackPt == 4){
        mainTitle->DrawLatexNDC(0.75, 0.94, Form("p_{T}^{ch} > %.1f GeV", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt)));
      } else {
        mainTitle->DrawLatexNDC(0.75, 0.94, Form("p_{T}^{ch} > %.0f GeV", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt)));
      }

      // Save the figures to file
      if(saveFigures) {
        gPad->GetCanvas()->SaveAs(Form("figures/finalBigCanvas%s%s.pdf", saveComment.Data(), compactTrackPtString.Data()));
      }

    }  // Track pT loop
  } // If for drawing big canvases

  if(drawCentralityRatios){

    // Draw all dividable canvas to draw all the histograms
    SplitCanvas* centralityRatioCanvas[nTrackPtBinsEEC];
    mainTitle = new TLatex();

    for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
      compactTrackPtString = Form("_T>%.1f", card[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));
      compactTrackPtString.ReplaceAll(".", "v");

      // Draw a big canvas and put all the plots in it
      centralityRatioCanvas[iTrackPt] = new SplitCanvas(Form("centralityRatioCanvas%d", iTrackPt), "", 1300, 400);
      centralityRatioCanvas[iTrackPt]->SetMargin(0.07, 0.01, 0.2, 0.01);
      centralityRatioCanvas[iTrackPt]->DivideNeatly(1, 4);

      bottomRowScale = centralityRatioCanvas[iTrackPt]->GetBottomRowScale();
      bottomPadMargin = centralityRatioCanvas[iTrackPt]->GetBottomPadMargin();
      leftPadMargin = centralityRatioCanvas[iTrackPt]->GetLeftPadMargin();

      for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
        jetPtString = Form("%.0f < jet p_{T} < %.0f", card[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));

        centralityRatioCanvas[iTrackPt]->CD(iJetPt - firstDrawnJetPtBinEEC);

        gPad->SetLogy(false);
        gPad->SetLogx();

        thisPadScale = (iJetPt == firstDrawnJetPtBinEEC) ? 0.85 : 1;

        // Set the axis titles and labels
        systematicUncertaintyPbPbToPpRatio[firstDrawnCentralityBin][iJetPt][iTrackPt]->GetXaxis()->SetTitleOffset(1/thisPadScale);
        systematicUncertaintyPbPbToPpRatio[firstDrawnCentralityBin][iJetPt][iTrackPt]->GetXaxis()->SetTitleSize(0.1*thisPadScale);
        systematicUncertaintyPbPbToPpRatio[firstDrawnCentralityBin][iJetPt][iTrackPt]->GetXaxis()->SetLabelOffset(0.01/thisPadScale);
        systematicUncertaintyPbPbToPpRatio[firstDrawnCentralityBin][iJetPt][iTrackPt]->GetXaxis()->SetLabelSize(0.08*thisPadScale);

        systematicUncertaintyPbPbToPpRatio[firstDrawnCentralityBin][iJetPt][iTrackPt]->GetYaxis()->SetTitleOffset(0.8/thisPadScale);
        systematicUncertaintyPbPbToPpRatio[firstDrawnCentralityBin][iJetPt][iTrackPt]->GetYaxis()->SetTitleSize(0.1*thisPadScale);
        systematicUncertaintyPbPbToPpRatio[firstDrawnCentralityBin][iJetPt][iTrackPt]->GetYaxis()->SetLabelOffset(0.01/thisPadScale);
        systematicUncertaintyPbPbToPpRatio[firstDrawnCentralityBin][iJetPt][iTrackPt]->GetYaxis()->SetLabelSize(0.08*thisPadScale);

        systematicUncertaintyPbPbToPpRatio[firstDrawnCentralityBin][iJetPt][iTrackPt]->GetYaxis()->SetTitle("#frac{PbPb}{pp}");
        systematicUncertaintyPbPbToPpRatio[firstDrawnCentralityBin][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
        systematicUncertaintyPbPbToPpRatio[firstDrawnCentralityBin][iJetPt][iTrackPt]->SetStats(0);
        systematicUncertaintyPbPbToPpRatio[firstDrawnCentralityBin][iJetPt][iTrackPt]->GetYaxis()->SetNdivisions(505);
        systematicUncertaintyPbPbToPpRatio[firstDrawnCentralityBin][iJetPt][iTrackPt]->GetYaxis()->CenterTitle();
        systematicUncertaintyPbPbToPpRatio[firstDrawnCentralityBin][iJetPt][iTrackPt]->GetXaxis()->SetTitle("#Deltar");
        systematicUncertaintyPbPbToPpRatio[firstDrawnCentralityBin][iJetPt][iTrackPt]->GetXaxis()->CenterTitle();
        systematicUncertaintyPbPbToPpRatio[firstDrawnCentralityBin][iJetPt][iTrackPt]->GetXaxis()->SetNdivisions(505);

        systematicUncertaintyPbPbToPpRatio[firstDrawnCentralityBin][iJetPt][iTrackPt]->SetTitle("");

        // Set the drawing ranges
        systematicUncertaintyPbPbToPpRatio[firstDrawnCentralityBin][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.006, 0.39);
        systematicUncertaintyPbPbToPpRatio[firstDrawnCentralityBin][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.4, 1.6);

        // Set the drawing style for histograms
        for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++) {
          systematicUncertaintyPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin]);
          systematicUncertaintyPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);
          energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin]);
          energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetMarkerSize(1.2);

          systematicUncertaintyPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin], 0.4);
          systematicUncertaintyPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
          energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
          systematicUncertaintyPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
          energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
        }

        // Draw the histograms
        systematicUncertaintyPbPbToPpRatio[firstDrawnCentralityBin][iJetPt][iTrackPt]->Draw("e2");
        systematicUncertaintyPbPbToPpRatio[lastDrawnCentralityBin][iJetPt][iTrackPt]->Draw("same,e2");
        energyEnergyCorrelatorPbPbToPpRatio[lastDrawnCentralityBin][iJetPt][iTrackPt]->Draw("same,p");
        energyEnergyCorrelatorPbPbToPpRatio[firstDrawnCentralityBin][iJetPt][iTrackPt]->Draw("same,p");
        //for(int iCentrality = firstDrawnCentralityBin + 1; iCentrality <= lastDrawnCentralityBin; iCentrality++) {
        //  systematicUncertaintyPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->Draw("same,e2");
        //}
        //for(int iCentrality = lastDrawnCentralityBin; iCentrality >= firstDrawnCentralityBin; iCentrality--) {
        //  energyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt]->Draw("same,p");
        //}

      }  // Jet pT loop

      // Save the figures to file
      if(saveFigures) {
        gPad->GetCanvas()->SaveAs(Form("figures/finalCentralityRatio%s%s.pdf", saveComment.Data(), compactTrackPtString.Data()));
      }
    } // Track pT loop
  } // If for drawing centrality ratios
}
