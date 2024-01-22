#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Macro for comparing energy-energy correlators in custom bins
 */
void compareEECinCustomBins(){
  
  // Files for comparison
  const int nComparisonFiles = 2;
  TString fileName[nComparisonFiles];
  fileName[0] = "data/eecAnalysis_akFlowJet_wtaAxis_energyWeightSquared_nominalReflectedCone_processed_2023-12-01.root";
  fileName[1] = "data/eecAnalysis_akFlowJet_energyWeightSquared_optimizedUnfoldingBins_noCovariance_processed_2024-01-17.root";
  
  
  // Open the files and check that they exist
  TFile* inputFile[nComparisonFiles];
  EECCard* card[nComparisonFiles];
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    
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
  TString collisionSystem = card[0]->GetDataType();
  bool isPbPbData = collisionSystem.Contains("PbPb");
  
  // ====================================================
  //               Binning configuration
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = card[0]->GetNCentralityBins();
  const int nJetPtBinsEEC = card[0]->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card[0]->GetNTrackPtBinsEEC();

  // For legends, collect common and not common selections
  std::vector<TString> commonLegend;
  std::vector<TString> individualLegend;
  std::vector<TString> legendComment;
  
  // Select explicitly which bins from the files are compared:
  std::vector<std::pair<double,double>> comparedCentralityBin;
  comparedCentralityBin.push_back(std::make_pair(0,100));
  comparedCentralityBin.push_back(std::make_pair(0,100));
  if(isPbPbData) commonLegend.push_back(Form("Cent: %.0f-%.0f%%", comparedCentralityBin.at(0).first, comparedCentralityBin.at(0).second));

  std::vector<std::pair<double,double>> comparedJetPtBin;
  comparedJetPtBin.push_back(std::make_pair(120,140));
  comparedJetPtBin.push_back(std::make_pair(135,155));
  for(auto jetPtBin : comparedJetPtBin){
    individualLegend.push_back(Form("%.0f < jet p_{T} < %.0f", jetPtBin.first, jetPtBin.second));
  }
  legendComment.push_back(" GeV (pp, no E-loss)");
  legendComment.push_back(" GeV (PbPb, before E-loss)");

  std::vector<double> comparedTrackPtBin;
  comparedTrackPtBin.push_back(3);
  comparedTrackPtBin.push_back(3);
  commonLegend.push_back(Form("%.1f < track p_{T}", comparedTrackPtBin.at(0)));

  // ====================================================
  //                    Configuration
  // ====================================================
  
  // Figure saving
  const bool saveFigures = true;  // Save figures
  const char* saveComment = "_pythia8";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures

  // Drawing configuration
  std::pair<double, double> ratioZoom = std::make_pair(0.7, 1.3);
  
  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager* histograms[nComparisonFiles];
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    histograms[iFile] = new EECHistogramManager(inputFile[iFile], card[iFile]);

    // Choose the energy-energy correlator types to load
    histograms[iFile]->SetLoadEnergyEnergyCorrelators(true);

    // Choose the bin ranges
    histograms[iFile]->SetCentralityBinRange(0, card[iFile]->GetNCentralityBins() - 1);
    histograms[iFile]->SetJetPtBinRangeEEC(0, card[iFile]->GetNJetPtBinsEEC() - 1);
    histograms[iFile]->SetTrackPtBinRangeEEC(0, card[iFile]->GetNTrackPtBinsEEC() - 1);

    // Load the histograms from the file
    histograms[iFile]->LoadProcessedHistograms();
  }

  // Energy-energy correlator histograms
  TH1D* hEnergyEnergyCorrelator[nComparisonFiles];
  TH1D* hEnergyEnergyCorrelatorRatio[nComparisonFiles]; // Ratio between raw and unfolded distributions
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    hEnergyEnergyCorrelator[iFile] = NULL;
    hEnergyEnergyCorrelatorRatio[iFile] = NULL;
  } // File loop
  
  // Helper histograms
  std::pair<double, double> drawingRange = std::make_pair(0.006, 0.39);
  double epsilon = 0.0001;
  int lowNormalizationBin, highNormalizationBin;
  int iCentrality;
  int iTrackPt;
  int iJetPt;

  // Get the histograms from the histogram manager and normalize the signal histograms to one
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){

    iCentrality = isPbPbData ? card[iFile]->FindBinIndexCentrality(comparedCentralityBin.at(iFile)) : 0;
    iJetPt = card[iFile]->FindBinIndexJetPtEEC(comparedJetPtBin.at(iFile));
    iTrackPt = card[iFile]->GetBinIndexTrackPtEEC(comparedTrackPtBin.at(iFile));

    cout << "iCentrality: " << iCentrality << " iJetPt: " << iJetPt << " iTrackPt: " << iTrackPt << endl;
    hEnergyEnergyCorrelator[iFile] = histograms[iFile]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt);
    if(hEnergyEnergyCorrelator[iFile] == NULL)  cout << "We Got NULL" << endl;

    // Normalize the distributions to one in the drawingRange
    lowNormalizationBin = hEnergyEnergyCorrelator[iFile]->GetXaxis()->FindBin(drawingRange.first + epsilon);
    highNormalizationBin = hEnergyEnergyCorrelator[iFile]->GetXaxis()->FindBin(drawingRange.second - epsilon);

    hEnergyEnergyCorrelator[iFile]->Scale(1 / hEnergyEnergyCorrelator[iFile]->Integral(lowNormalizationBin, highNormalizationBin, "width"));

    // Calculate the ratio to first index
    hEnergyEnergyCorrelatorRatio[iFile] = (TH1D*)hEnergyEnergyCorrelator[iFile]->Clone(Form("eecRatio%d", iFile));
    hEnergyEnergyCorrelatorRatio[iFile]->Divide(hEnergyEnergyCorrelator[0]);

  } // File loop

  
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

  TString centralityString, trackPtString, jetPtString;
  TString compactCentralityString, compactTrackPtString, compactJetPtString;
  TString legendString;
  int markerStyle[2] = {kOpenSquare, kOpenCircle};
  int color[] = {kBlack, kRed, kBlue};

  // Create a new canvas for the plot
  drawer->CreateSplitCanvas();
          
  // Logarithmic EEC axis
  drawer->SetLogY(true);

  TLegend* legend = new TLegend(0.18,0.04,0.45,0.48);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  legend->AddEntry((TObject*) 0, Form("%s 5.02 TeV",card[0]->GetAlternativeDataType(false).Data()), "");

  // Add common legend variables
  for(TString legendItem : commonLegend){
    legend->AddEntry((TObject*) 0, legendItem, "");
  }

  // Set drawing style for all histograms
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    hEnergyEnergyCorrelator[iFile]->SetMarkerStyle(markerStyle[0]);
    hEnergyEnergyCorrelator[iFile]->SetMarkerColor(color[iFile]);
    hEnergyEnergyCorrelator[iFile]->SetLineColor(color[iFile]);
    hEnergyEnergyCorrelatorRatio[iFile]->SetMarkerStyle(markerStyle[0]);
    hEnergyEnergyCorrelatorRatio[iFile]->SetMarkerColor(color[iFile]);
    hEnergyEnergyCorrelatorRatio[iFile]->SetLineColor(color[iFile]);
  }

  // Set the x-axis drawing range
  hEnergyEnergyCorrelator[0]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
          
  // Draw the histograms to the upper canvas
  drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelator[0], "#Deltar", "EEC", " ");

  for(int iFile = 1; iFile < nComparisonFiles; iFile++){
    hEnergyEnergyCorrelator[iFile]->Draw("same");
  }

  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    legend->AddEntry(hEnergyEnergyCorrelator[iFile], Form("%s%s", individualLegend.at(iFile).Data(), legendComment.at(iFile).Data()), "p");
  }
  
  // Draw the legends to the upper pad
  legend->Draw();
          
  // Linear scale for the ratio
  drawer->SetLogY(false);
          
  // Set the drawing ranges
  hEnergyEnergyCorrelatorRatio[1]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);          
  hEnergyEnergyCorrelatorRatio[1]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

  // Draw the histograms
  drawer->SetGridY(true);
  drawer->DrawHistogramToLowerPad(hEnergyEnergyCorrelatorRatio[1], "#Deltar", Form("#frac{%s}{%s}", individualLegend.at(1).Data(), individualLegend.at(0).Data()), " ");
  for(int iFile = 2; iFile < nComparisonFiles; iFile++){
    hEnergyEnergyCorrelatorRatio[iFile]->Draw("same");
  }
  drawer->SetGridY(false);
          
  // Save the figures to a file
  if(saveFigures){
  gPad->GetCanvas()->SaveAs(Form("figures/eecCustomBinComparison%s%s%s%s.%s", saveComment, compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
  }

}
