#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Shift the jet pT spectrum in Pythia8 in order to reproduce the jet RAA in data
 */
void fakeJetRaa(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString inputFileName = "data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_regularCorrelators_processed_2023-10-02.root";
  
  // Drawing configuration
  std::pair<double, double> jetPtZoomRange = std::make_pair(100, 400);
  std::pair<double, double> ratioZoomRange = std::make_pair(0, 1.1);

  int color[] = {kBlack, kRed, kBlue, kGreen+3, kMagenta, kCyan};

  // Figure saving configuration
  const bool saveFigures = true;  // Save figures
  const char* saveComment = "";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures

  // ==================================================================
  // =================== Configuration ready ==========================
  // ==================================================================
  
  // Open the input files
  TFile* inputFile = TFile::Open(inputFileName);
  
  // Read the card
  EECCard* card = new EECCard(inputFile);
  
  // Create histogram managers to provide the histograms for the study
  EECHistogramManager *pythiaHistograms = new EECHistogramManager(inputFile);
  pythiaHistograms->SetLoadJetHistograms(true);
  pythiaHistograms->LoadProcessedHistograms();

  // The jet pT histogram has 5 GeV bin resolution, so shifts are done in 5 GeV steps
  const int nShifts = 5;
  TH1D* hJetPt[nShifts+1];
  TH1D* hJetRaa[nShifts+1];
  
  // Initialize the histograms to null
  for(int iShift = 0; iShift < nShifts+1; iShift++){
    hJetPt[iShift] = NULL;
    hJetRaa[iShift] = NULL;
  }

  // Read the jet pT distribution to bin 0
  hJetPt[0] = pythiaHistograms->GetHistogramJetPt(0);

  // Shift the bin contents by one bin
  const int nBins = hJetPt[0]->GetNbinsX();
  for(int iShift = 1; iShift <= nShifts; iShift++){
    hJetPt[iShift] = (TH1D*) hJetPt[0]->Clone(Form("shiftedJetPt%d",iShift));
    for(int iBin = 1; iBin <= nBins - iShift; iBin++){
      hJetPt[iShift]->SetBinContent(iBin, hJetPt[0]->GetBinContent(iBin+iShift));
      hJetPt[iShift]->SetBinError(iBin, hJetPt[0]->GetBinError(iBin+iShift));
    }
  }

  // Calculate jet RAA in all different cases
  for(int iShift = 0; iShift <= nShifts; iShift++){
    hJetPt[iShift]->SetLineColor(color[iShift]);
    hJetRaa[iShift] = (TH1D*) hJetPt[iShift]->Clone(Form("jetRaa%d",iShift));
    hJetRaa[iShift]->Divide(hJetPt[0]);
  }


  // Draw all the histograms to the canvas  
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);

  // Create a new canvas for the plot
  drawer->CreateSplitCanvas();

  TLegend* legend = new TLegend(0.58,0.4,0.85,0.84);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  legend->AddEntry((TObject*) 0, Form("%s 5.02 TeV",card->GetAlternativeDataType(false).Data()), "");

  drawer->SetLogY(true);
  hJetPt[0]->GetXaxis()->SetRangeUser(jetPtZoomRange.first, jetPtZoomRange.second);
  drawer->DrawHistogramToUpperPad(hJetPt[0], "Jet p_{T}", "dN/dp_{T}", " ");
  legend->AddEntry(hJetPt[0], "Reference spectrum", "l");
  for(int iShift = 1; iShift <= nShifts; iShift++){
    hJetPt[iShift]->Draw("same");
    legend->AddEntry(hJetPt[iShift], Form("E-loss = %d GeV", 5*iShift), "l");
  }
  legend->Draw();

  drawer->SetLogY(false);
  hJetRaa[0]->GetXaxis()->SetRangeUser(jetPtZoomRange.first, jetPtZoomRange.second);
  hJetRaa[0]->GetYaxis()->SetRangeUser(ratioZoomRange.first, ratioZoomRange.second);
  drawer->DrawHistogramToLowerPad(hJetRaa[0], "Jet p_{T}", "Jet R_{AA}", " ");
  for(int iShift = 1; iShift <= nShifts; iShift++){
    hJetRaa[iShift]->Draw("same");
  }

  // Save the figures to a file
  if(saveFigures){
  gPad->GetCanvas()->SaveAs(Form("figures/fakeJetRaa%s.%s", saveComment, figureFormat));
  }
  
}
