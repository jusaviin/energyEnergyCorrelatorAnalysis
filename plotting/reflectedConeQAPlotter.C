#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Macro plotting QA histograms for reflected cone
 */
void reflectedConeQAPlotter(){

  enum enumDataType{kPbPb, kPp, kNDataTypes};

  // ============= //
  // Configuration //
  // ============= //
  
  // Input files
  TString dataFileName = "data/eecAnalysis_akFlowJet_wtaAxis_reflectedConeQA_processed_2023-08-17.root";
  TFile* dataFile = TFile::Open(dataFileName);

  // Check that the input file exists
  if(dataFile == NULL){
    cout << "Error! The file " << dataFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }

  // Load the card from the file
  EECCard* dataCard = new EECCard(dataFile);

  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Find the number of centrality bins from the card
  const int nCentralityBins = dataCard->GetNCentralityBins();
  
  // Choose which bins are drawn
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = nCentralityBins-1;
  
  // Save the final plots
  const bool saveFigures = false;
  TString saveComment = "_firstLook";

  // =============================================== //
  // Read the histograms from the histogram managers //
  // =============================================== //

  // Create histogram manager for the reflected cone QA histograms
  EECHistogramManager* dataHistograms = new EECHistogramManager(dataFile, dataCard);
  dataHistograms->SetLoadReflectedConeQAHistograms(true);
  dataHistograms->SetCentralityBinRange(firstDrawnCentralityBin, lastDrawnCentralityBin);

  // Load the histograms from the file
  dataHistograms->LoadProcessedHistograms();
 
  TH1D* numberOfJetsWithinReflectedCone[nCentralityBins];
  TH1D* jetPtWithinReflectedCone[nCentralityBins];
  TH1D* jetPtWithinReflectedConeTrunkated[nCentralityBins];

  // Initialize histograms to NULL
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    numberOfJetsWithinReflectedCone[iCentrality] = NULL;
    jetPtWithinReflectedCone[iCentrality] = NULL;
    jetPtWithinReflectedConeTrunkated[iCentrality] = NULL;
  } // Centrality loop

  // Read the histograms from the histogram manager
  for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){

    // Read the PbPb histograms
    numberOfJetsWithinReflectedCone[iCentrality] = dataHistograms->GetHistogramNumberOfJetsWithinReflectedCone(iCentrality);
    jetPtWithinReflectedCone[iCentrality] = dataHistograms->GetHistogramJetPtWithinReflectedCone(iCentrality);

  } // Centrality loop

  // Trunkate the last bin in the jet pT histogram
  const int nJetPtBins = jetPtWithinReflectedCone[firstDrawnCentralityBin]->GetNbinsX();
  double trunkatedBins[nJetPtBins+1];
  for(int iBin = 1; iBin <= nJetPtBins; iBin++){
    trunkatedBins[iBin-1] = jetPtWithinReflectedCone[firstDrawnCentralityBin]->GetBinLowEdge(iBin);
  }
  trunkatedBins[nJetPtBins] = trunkatedBins[nJetPtBins-1] + jetPtWithinReflectedCone[firstDrawnCentralityBin]->GetBinWidth(nJetPtBins-1);

  // Create the histograms with trunkated bins
  for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
    jetPtWithinReflectedConeTrunkated[iCentrality] = new TH1D(Form("trunkatedBinJetPtWithinReflectedCone%d",iCentrality),Form("trunkatedBinJetPtWithinReflectedCone%d",iCentrality), nJetPtBins, trunkatedBins);
    jetPtWithinReflectedConeTrunkated[iCentrality]->Sumw2();
    for(int iBin = 1; iBin <= nJetPtBins; iBin++){
      jetPtWithinReflectedConeTrunkated[iCentrality]->SetBinContent(iBin,jetPtWithinReflectedCone[iCentrality]->GetBinContent(iBin));
      jetPtWithinReflectedConeTrunkated[iCentrality]->SetBinError(iBin,jetPtWithinReflectedCone[iCentrality]->GetBinError(iBin));
    }
  }

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

  // Common variables for different plots
  TLegend* legend;

  // Draw plots for all studied bins
  for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++) {
    centralityString = Form("Cent: %.0f-%.0f%%", dataCard->GetLowBinBorderCentrality(iCentrality), dataCard->GetHighBinBorderCentrality(iCentrality));
    compactCentralityString = Form("_C=%.0f-%.0f", dataCard->GetLowBinBorderCentrality(iCentrality), dataCard->GetHighBinBorderCentrality(iCentrality));

    // Setup the legend for plots
    legend = new TLegend(0.53, 0.6, 0.83, 0.85);
    legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);
    legend->AddEntry((TObject*)0, centralityString.Data(), "");

    // Draw the number of jets histogram
    drawer->DrawHistogram(numberOfJetsWithinReflectedCone[iCentrality], "N_{jets} within reflected cone", "Counts", " ");

    // Draw the legend
    legend->Draw();

    // Save the number of jets figure to a file
    if(saveFigures) {
      gPad->GetCanvas()->SaveAs(Form("figures/numberOfJetsWithinReflectedCone%s%s.pdf", saveComment.Data(), compactCentralityString.Data()));
    }

    // Draw the jets pT histogram
    drawer->DrawHistogram(jetPtWithinReflectedConeTrunkated[iCentrality], "Jet p_{T} within reflected cone", "Counts", " ");

    // Draw the legend
    legend->Draw();

    // Save the number of jets figure to a file
    if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/jetPtWithinReflectedCone%s%s.pdf", saveComment.Data(), compactCentralityString.Data()));
    }

  } // Centrality loop

  // Print a table of number of jets above 25 GeV within the reflected cone to console in a form of beamer slide
  cout << endl;
  cout << "\\begin{frame}" << endl;
  cout << "\\frametitle{Number of jets above 25 GeV within reflected cone}" << endl;
  cout << "\\begin{center}" << endl;
  cout << "  \\begin{tabular}{ccccc}" << endl;
  cout << "    \\toprule" << endl;
  cout << "    N_{jets} & 0-10\\% & 10-30\\% & 30-50\\% & 50-90\\% \\\\" << endl;
  cout << "    \\midrule" << endl;

  // Set the correct precision for printing floating point numbers
  cout << fixed << setprecision(3);

  for(int iBin = 1; iBin <= numberOfJetsWithinReflectedCone[firstDrawnCentralityBin]->GetNbinsX(); iBin++){
    cout << "    " << iBin-1;
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      cout << " & " << numberOfJetsWithinReflectedCone[iCentrality]->GetBinContent(iBin) / numberOfJetsWithinReflectedCone[iCentrality]->Integral();
    }
    cout << " \\\\" << endl;
  }

  cout << "    \\bottomrule" << endl;
  cout << "  \\end{tabular}" << endl;
  cout << "\\end{center}" << endl;
  cout << "\\end{frame}" << endl;
  cout << endl;

  // Print a table of interesting results to console in a form of beamer slide
  cout << endl;
  cout << "\\begin{frame}" << endl;
  cout << "\\frametitle{Jet pT within reflected cone}" << endl;
  cout << "\\begin{center}" << endl;
  cout << "  \\begin{tabular}{ccccc}" << endl;
  cout << "    \\toprule" << endl;
  cout << "    & 0-10\\% & 10-30\\% & 30-50\\% & 50-90\\% \\\\" << endl;
  cout << "    \\midrule" << endl;

  // Set the correct precision for printing floating point numbers
  cout << fixed << setprecision(3);

  for(int iBin = 1; iBin <= numberOfJetsWithinReflectedCone[firstDrawnCentralityBin]->GetNbinsX(); iBin++){
    if(iBin < numberOfJetsWithinReflectedCone[firstDrawnCentralityBin]->GetNbinsX()){
      cout << Form("   $%.0f < p_{\\mathrm{T}} < %.0f$", trunkatedBins[iBin-1], trunkatedBins[iBin]);
    } else {
      cout << "   $100 < p_{\\mathrm{T}}$";
    }
    for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
      cout << " & " << jetPtWithinReflectedCone[iCentrality]->GetBinContent(iBin) / jetPtWithinReflectedCone[iCentrality]->Integral();
    }
    cout << " \\\\" << endl;
  }

  cout << "    \\bottomrule" << endl;
  cout << "  \\end{tabular}" << endl;
  cout << "\\end{center}" << endl;
  cout << "\\end{frame}" << endl;
  cout << endl;
}
