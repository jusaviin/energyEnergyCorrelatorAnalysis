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

  firstDrawnTrackPtBinEEC = 5;
  lastDrawnTrackPtBinEEC = 5;

  // Choose which plots to draw
  bool drawIndividualPlotsAllCentralities = true;
  bool drawBigCanvasDistributions = false;
  
  // Save the final plots
  const bool saveFigures = true;
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
        TLegend* legend = new TLegend(0.23, 0.05, 0.53, 0.6);
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
    SplitCanvas* bigCanvas;
    TLatex* mainTitle;
    TLegend* legend;
    
    // Draw a big canvas and put all the plots in it
    bigCanvas = new SplitCanvas("bigCanvas", "", 1800, 1500);
    bigCanvas->SetMargin(0.08, 0.01, 0.08, 0.11);
    bigCanvas->DivideNeatly(4,4);
    
    mainTitle = new TLatex();

    int canvasIndex;
    double bottomRowScale = bigCanvas->GetBottomRowScale();
    double bottomPadMargin = bigCanvas->GetBottomPadMargin();
    double leftPadMargin = bigCanvas->GetLeftPadMargin();
    double thisScale = 1;
    double leftMarginAdder, bottomMarginAdder;
    
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){

        jetPtString = Form("%.0f < jet p_{T} < %.0f", card[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), card[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));

        canvasIndex = (iCentrality - firstDrawnCentralityBin) * (lastDrawnJetPtBinEEC - firstDrawnJetPtBinEEC + 1) + (iJetPt - firstDrawnJetPtBinEEC);
        bigCanvas->CD(canvasIndex);

        gPad->SetLogy();
        gPad->SetLogx();

        // The titles in the bottow row need to be scaled, since it has more margin than other rows
        thisScale = (iCentrality == lastDrawnCentralityBin) ? bottomRowScale : 1;
 
        // Set the axis titles and labels
        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->GetXaxis()->SetTitleOffset(0.85);
        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->GetXaxis()->SetTitleSize(0.2*thisScale);
        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->GetXaxis()->SetLabelOffset(0.01);
        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->GetXaxis()->SetLabelSize(0.12*thisScale);

        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->GetYaxis()->SetTitleOffset(0.5/thisScale);
        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->GetYaxis()->SetTitleSize(0.2*thisScale);
        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->GetYaxis()->SetLabelOffset(0.01/thisScale);
        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->GetYaxis()->SetLabelSize(0.12*thisScale);

        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->GetYaxis()->SetTitle("EEC");
        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->GetXaxis()->CenterTitle();
        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->SetStats(0);
        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->GetYaxis()->SetNdivisions(505);
        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->GetYaxis()->CenterTitle();
        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->GetXaxis()->SetTitle("#Deltar");
        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->GetXaxis()->CenterTitle();
        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->GetXaxis()->SetNdivisions(505);

        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->SetTitle("");

        // Set the drawing ranges
        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->GetXaxis()->SetRangeUser(0.006, 0.39);
        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->GetYaxis()->SetRangeUser(0.15, 30);

        // Set the drawing style for histograms
        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin]);
        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->SetMarkerSize(1.2);
        energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][5]->SetMarkerStyle(markerStylePbPb[iCentrality - firstDrawnCentralityBin]);
        energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][5]->SetMarkerSize(1.2);

        systematicUncertaintyForPp[iJetPt][5]->SetMarkerStyle(markerStylePp);
        systematicUncertaintyForPp[iJetPt][5]->SetMarkerSize(1.2);
        energyEnergyCorrelatorSignalPp[iJetPt][5]->SetMarkerStyle(markerStylePp);
        energyEnergyCorrelatorSignalPp[iJetPt][5]->SetMarkerSize(1.2);
    
        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->SetFillColorAlpha(markerColorPbPb[iCentrality - firstDrawnCentralityBin], 0.4);
        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
        energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][5]->SetLineColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);
        energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][5]->SetMarkerColor(markerColorPbPb[iCentrality - firstDrawnCentralityBin]);

        systematicUncertaintyForPp[iJetPt][5]->SetFillColorAlpha(markerColorPp, 0.4);
        systematicUncertaintyForPp[iJetPt][5]->SetLineColor(markerColorPp);
        energyEnergyCorrelatorSignalPp[iJetPt][5]->SetLineColor(markerColorPp);
        systematicUncertaintyForPp[iJetPt][5]->SetMarkerColor(markerColorPp);
        energyEnergyCorrelatorSignalPp[iJetPt][5]->SetMarkerColor(markerColorPp);
      
        // Draw the histograms
        systematicUncertaintyForPbPb[iCentrality][iJetPt][5]->Draw("e2");
        systematicUncertaintyForPp[iJetPt][5]->Draw("same,e2");
        energyEnergyCorrelatorSignalPp[iJetPt][5]->Draw("same,p");
        energyEnergyCorrelatorSignalPbPb[iCentrality][iJetPt][5]->Draw("same,p");
        
      
        // Create a legend for each pad
        leftMarginAdder = (iJetPt == firstDrawnJetPtBinEEC) ? leftPadMargin : 0;
        bottomMarginAdder = (iCentrality == lastDrawnCentralityBin) ? bottomPadMargin : 0;
        legend = new TLegend(0.06+leftMarginAdder,0.05+bottomMarginAdder,0.5 / (1 - leftMarginAdder),0.45/thisScale);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.09*thisScale);legend->SetTextFont(62);
        legend->AddEntry((TObject*)0, jetPtString.Data(), "");
        legend->AddEntry(systematicUncertaintyForPbPb[iCentrality][iJetPt][5], Form("PbPb %.0f-%.0f%%", card[kPbPb]->GetLowBinBorderCentrality(iCentrality), card[kPbPb]->GetHighBinBorderCentrality(iCentrality)), "p");
        legend->AddEntry(systematicUncertaintyForPp[iJetPt][5], "pp", "p");
      
        legend->Draw();
      
      } // Jet pT loop
    } // Centrality loop
    
    // Draw all necessary CMS text to the plot
    bigCanvas->cd(0);
    
    mainTitle->SetTextFont(62);
    mainTitle->SetTextSize(0.065);
    mainTitle->DrawLatexNDC(0.05, 0.93, "CMS");
    
    mainTitle->SetTextFont(42);
    mainTitle->SetTextSize(0.055);
    mainTitle->DrawLatexNDC(0.2, 0.93, "Add here luminosity, track p_{T} cut, etc.");
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
    
    // Save the figures to file
    if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/finalBigCanvas%s.pdf", saveComment.Data()));
    }
    
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
  }
}
