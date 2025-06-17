#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Macro for finding background normalization scale in MC. This is done by integrating the energy-energy correlator distributions for background within the signal cone, and in the reflected cone, and taking a ratio of these numbers.
 */
void hfEnergyPlotter(){

  // File from which the integrals are calculated
  TString inputFileName = "data/pPb/pPb_5TeV_pToMinusEta_pfJets_eschemeAxis_nominalEnergyWeight_minimumBias_jetEtaMCcut_HFEnergy_processed_2025-06-17.root";
  
  // Open the input file
  TFile* inputFile = TFile::Open(inputFileName);
  
  if(inputFile == NULL){
    cout << "Error! The file " << inputFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/pPb/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }
  
  // Load the card from the file and read the collision system
  EECCard* card = new EECCard(inputFile);
  
  // ====================================================
  //  Binning configuration for the integral calculation
  // ====================================================
  
  // Find the number of bins from the card
  const int nJetPtBins = card->GetNJetPtBinsEEC();
  
  // Bins jet pT bins to plot
  std::vector<std::pair<double,double>> comparedJetPtBin;
  comparedJetPtBin.push_back(std::make_pair(30,40));
  comparedJetPtBin.push_back(std::make_pair(40,50));
  comparedJetPtBin.push_back(std::make_pair(50,60));
  comparedJetPtBin.push_back(std::make_pair(60,80));
  comparedJetPtBin.push_back(std::make_pair(0,0));
  bool plotToSameFigure = false; // False = make different figure for each bin. True = plot all jet pT bins to the same figure.
    
  // Create and setup a new histogram manager to project and handle the histograms
  EECHistogramManager* histograms = new EECHistogramManager(inputFile,card);

  // Initialize the HF energy histograms to NULL
  TH1D* hHFPlus[nJetPtBins+1];
  TH1D* hHFMinus[nJetPtBins+1];
  TH1D* hHFSum[nJetPtBins+1];
  TH2D* hHFPlusVsHFMinus[nJetPtBins+1];

  for(int iJetPt = 0; iJetPt < nJetPtBins+1; iJetPt++){

    hHFPlus[iJetPt] = NULL;
    hHFMinus[iJetPt] = NULL;
    hHFSum[iJetPt] = NULL;
    hHFPlusVsHFMinus[iJetPt] = NULL;

  }  // Leading jet pT loop

  // Read the HF energy histograms from the file
  int iJetPt;
  for(auto jetPtBin : comparedJetPtBin){

    iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);
    if(iJetPt == -1) iJetPt = nJetPtBins; // If a given jet pT bin is not found, use the whole jet pT range

    // Read the histograms from the file
    hHFPlus[iJetPt] = histograms->GetHistogramHFPlus(iJetPt);
    hHFMinus[iJetPt] = histograms->GetHistogramHFMinus(iJetPt);
    hHFSum[iJetPt] = histograms->GetHistogramHFSum(iJetPt);
    hHFPlusVsHFMinus[iJetPt] = histograms->GetHistogramHFPlusVsHFMinus(iJetPt);

  }  // Centrality loop

  // **********************************
  //         Draw the figures
  // **********************************

  // Prepare a JDrawer for drawing purposes
  JDrawer* drawer = new JDrawer();
  drawer->SetRelativeCanvasSize(1.2, 1);
  drawer->SetLogZ(true);
  drawer->SetLeftMargin(0.15);
  drawer->SetRightMargin(0.14);
  drawer->SetTopMargin(0.09);
  drawer->SetBottomMargin(0.2);
  drawer->SetTitleOffsetX(1.17);
  drawer->SetTitleOffsetY(1.2);

  TLine* diagonal = new TLine(0,0,300,300);
  diagonal->SetLineColor(kRed);

  // Helper variables
  TLegend* legend;
  TString jetPtString;
  TString compactJetPtString;

  // Draw the two dimensional distributions to canvas
  for(auto jetPtBin : comparedJetPtBin){

    iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);

    if(iJetPt == -1){
      iJetPt = nJetPtBins;
      jetPtString = "No jet selection";
    } else {
      jetPtString = Form("%.0f < leading jet p_{T} < %.0f GeV", jetPtBin.first, jetPtBin.second);
    }


    drawer->DrawHistogram(hHFPlusVsHFMinus[iJetPt],"HF Pb going", "HF p going", jetPtString.Data(), "colz");
    diagonal->Draw();

  } // Jet pT loop

  // Reset the drawing style for 1-dimensional histograms
  drawer->Reset();

}
