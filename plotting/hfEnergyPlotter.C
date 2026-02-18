#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Macro for finding background normalization scale in MC. This is done by integrating the energy-energy correlator distributions for background within the signal cone, and in the reflected cone, and taking a ratio of these numbers.
 */
void hfEnergyPlotter(){

  // File from which the integrals are calculated
  TString inputFileName = "veryEposData_processed.root";

  // data/pPb/ppData_pfJets_eschemeAxis_nominalEnergyWeight_jetEtaCMcut_zeroBias_processed_2025-06-23.root
  // data/pPb/pPb_5TeV_pToMinusEta_pfJets_eschemeAxis_nominalEnergyWeight_minimumBias_jetEtaMCcut_HFEnergyNjet_processed_2025-06-20.root
  //"data/pPb/pPb_5TeV_pToMinusEta_pfJets_eschemeAxis_nominalEnergyWeight_minimumBias_jetEtaMCcut_HFEnergy_processed_2025-06-17.root";
  // data/pPb/monteCarlo/pythiaEpos_8TeV_GenGen_pToPlusEta_genJets_eschemeAxis_nominalEnergyWeight_jetEtaCMcut_megaSkimMix_9filesMissing_processed_2026-02-13.root
  
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
  double minJetPt = card->GetJetPtCut();

  // ====================================================
  //            Select which plots to draw
  // ====================================================
  const bool draw2DHFmap = false;
  const bool drawHFsum = true;
  const bool drawHFSumGraph = true;
  
  // ====================================================
  //  Binning configuration for the integral calculation
  // ====================================================
  
  // Find the number of bins from the card
  const int nJetPtBins = card->GetNJetPtBinsEEC();
  const int nNJetBins = 3;
  
  // Bins jet pT bins to plot
  std::vector<std::pair<double,double>> comparedJetPtBin;
  comparedJetPtBin.push_back(std::make_pair(30,40));
  comparedJetPtBin.push_back(std::make_pair(40,50));
  comparedJetPtBin.push_back(std::make_pair(50,60));
  comparedJetPtBin.push_back(std::make_pair(60,80));
  comparedJetPtBin.push_back(std::make_pair(80,100));
  comparedJetPtBin.push_back(std::make_pair(0,0));
  bool plotToSameFigure = false; // False = make different figure for each bin. True = plot all jet pT bins to the same figure.
    
  // Create and setup a new histogram manager to project and handle the histograms
  EECHistogramManager* histograms = new EECHistogramManager(inputFile,card);

  // Initialize the HF energy histograms to NULL
  TH1D* hHFPlus[nJetPtBins+1][nNJetBins+1];
  TH1D* hHFMinus[nJetPtBins+1][nNJetBins+1];
  TH1D* hHFSum[nJetPtBins+1][nNJetBins+1];
  TH2D* hHFPlusVsHFMinus[nJetPtBins+1][nNJetBins+1];

  for(int iJetPt = 0; iJetPt < nJetPtBins+1; iJetPt++){
    for(int iNJet = 0; iNJet < nNJetBins+1; iNJet++){

      hHFPlus[iJetPt][iNJet] = NULL;
      hHFMinus[iJetPt][iNJet] = NULL;
      hHFSum[iJetPt][iNJet] = NULL;
      hHFPlusVsHFMinus[iJetPt][iNJet] = NULL;

    } // Number of jets loop
  }  // Leading jet pT loop

  // Graphs with jet pT dependent mean HF sum
  TGraph* gAverageHFSum[nNJetBins];

  // Initialize the graphs to NULL
  for(int iNJet = 0; iNJet < nNJetBins; iNJet++){
    gAverageHFSum[iNJet] = NULL;
  }

  // Read the HF energy histograms from the file
  int iJetPt;
  int nJetIndex;
  for(auto jetPtBin : comparedJetPtBin){

    iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);
    if(iJetPt == -1) iJetPt = nJetPtBins; // If a given jet pT bin is not found, use the whole jet pT range

    for(int iNJet = 0; iNJet < nNJetBins+1; iNJet++){

      nJetIndex = iNJet;
      if(nJetIndex == nNJetBins) nJetIndex = -1;

      // Read the histograms from the file
      hHFPlus[iJetPt][iNJet] = histograms->GetHistogramHFPlus(iJetPt, nJetIndex);
      hHFMinus[iJetPt][iNJet] = histograms->GetHistogramHFMinus(iJetPt, nJetIndex);
      hHFSum[iJetPt][iNJet] = histograms->GetHistogramHFSum(iJetPt, nJetIndex);
      hHFPlusVsHFMinus[iJetPt][iNJet] = histograms->GetHistogramHFPlusVsHFMinus(iJetPt, nJetIndex);

    } // Number of jets loop
  }  // Jet pT loop

  // *******************************************************************************************
  //        Create graphs that contain jet pT dependent HF sum based on number of jets
  // *******************************************************************************************

  // Store the points that are used to create the graphs in vectors
  std::vector<double> grahpXpoints;
  std::vector<double> graphYpoints[nNJetBins];
  for(auto jetPtBin : comparedJetPtBin){
    
    iJetPt = card->FindBinIndexJetPtEEC(jetPtBin); 
    if(iJetPt == -1) continue; // Ignore the bin without jet pT selection

    // For now, place the x-point in the center of the bin
    // This can later the be improved to be the mean of the jet pT bin
    grahpXpoints.push_back((jetPtBin.first + jetPtBin.second) / 2.0);

    // The Y-axis value is the mean HF energy in this bin.
    // This can be later improved to contain the error of the mean
    for(int iNJet = 1; iNJet <= nNJetBins; iNJet++){
      nJetIndex = iNJet;
      if(iNJet == nNJetBins) nJetIndex = 0;
      graphYpoints[nJetIndex].push_back(hHFSum[iJetPt][iNJet]->GetMean());
    } // Number of jets loop
  } // Jet pT loop

  // Once we have the points needed to make the graphs, create the graphs
  for(int iNJet = 0; iNJet < nNJetBins; iNJet++){
    gAverageHFSum[iNJet] = new TGraph(grahpXpoints.size(), grahpXpoints.data(), graphYpoints[iNJet].data());
  }

  // Create also a line showing the HF energy when it is required the no jets are in the event
  TLine* zeroJetsLine = new TLine(comparedJetPtBin.at(0).first, hHFSum[nJetPtBins][0]->GetMean(), comparedJetPtBin.at(comparedJetPtBin.size()-2).second, hHFSum[nJetPtBins][0]->GetMean());
  zeroJetsLine->SetLineColor(kRed);
  zeroJetsLine->SetLineStyle(2);

  double hfShiftValue = 25;
  if(card->GetAlternativeDataType().Contains("pp")) hfShiftValue = 25;
  if(card->GetAlternativeDataType().Contains("pPb")) hfShiftValue = 28;
  TLine* estimatedShiftLine = new TLine(comparedJetPtBin.at(0).first, hHFSum[nJetPtBins][0]->GetMean()+hfShiftValue, comparedJetPtBin.at(comparedJetPtBin.size()-2).second, hHFSum[nJetPtBins][0]->GetMean()+hfShiftValue);
  estimatedShiftLine->SetLineColor(kBlack);
  estimatedShiftLine->SetLineStyle(2);
  TLine* estimatedShiftLineUp = new TLine(comparedJetPtBin.at(0).first, hHFSum[nJetPtBins][0]->GetMean()+hfShiftValue+3, comparedJetPtBin.at(comparedJetPtBin.size()-2).second, hHFSum[nJetPtBins][0]->GetMean()+hfShiftValue+3);
  estimatedShiftLineUp->SetLineColor(kBlack);
  estimatedShiftLineUp->SetLineStyle(2);
  TLine* estimatedShiftLineDown = new TLine(comparedJetPtBin.at(0).first, hHFSum[nJetPtBins][0]->GetMean()+hfShiftValue-3, comparedJetPtBin.at(comparedJetPtBin.size()-2).second, hHFSum[nJetPtBins][0]->GetMean()+hfShiftValue-3);
  estimatedShiftLineDown->SetLineColor(kBlack);
  estimatedShiftLineDown->SetLineStyle(2);

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
  TString nJetString;

  // Define nice markers and colors
  int markerStyle[] = {kFullCircle, kFullSquare, kFullCross, kFullDiamond, kFullStar, kFullCrossX, kFullDoubleDiamond};
  int color[] = {kBlack,kRed,kBlue,kGreen+3,kMagenta,kCyan,kOrange,kViolet+3,kPink-7,kSpring+3,kAzure-7};

  // Draw the two dimensional HF energy map
  if(draw2DHFmap){
    for(auto jetPtBin : comparedJetPtBin){

      iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);

      if(iJetPt == -1){
        iJetPt = nJetPtBins;
        jetPtString = "No jet selection";
      } else {
        jetPtString = Form("%.0f < leading jet p_{T} < %.0f GeV", jetPtBin.first, jetPtBin.second);
      }


      drawer->DrawHistogram(hHFPlusVsHFMinus[iJetPt][nNJetBins],"HF Pb going", "HF p going", jetPtString.Data(), "colz");
      diagonal->Draw();

    } // Jet pT loop
  } // Draw 2D HF energy map 

  // Reset the drawing style for 1-dimensional histograms
  drawer->Reset();

  // Draw the total energy in HF calorimeters
  if(drawHFsum){
    for(auto jetPtBin : comparedJetPtBin){
      for(int iNJet = 0; iNJet < nNJetBins+1; iNJet++){

        iJetPt = card->FindBinIndexJetPtEEC(jetPtBin);

        // Define legend
        legend = new TLegend(0.48,0.34,0.75,0.68);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*)0, Form("%s  5.02 TeV", card->GetAlternativeDataType().Data()), "");

        if(iJetPt == -1){
          iJetPt = nJetPtBins;
          jetPtString = "No jet selection";
        } else {
          jetPtString = Form("%.0f < leading jet p_{T} < %.0f GeV", jetPtBin.first, jetPtBin.second);
        }

        if(iNJet == nNJetBins){
          nJetString = "No number of jet requirement";
        } else {
          nJetString = Form("Require %d jets above %.0f GeV", iNJet, minJetPt);
        }


        drawer->DrawHistogram(hHFSum[iJetPt][iNJet], "HF energy", "Counts", " ");
        legend->AddEntry(hHFSum[iJetPt][iNJet], jetPtString.Data(), "l");
        legend->AddEntry((TObject*)0, nJetString.Data(), "");
        legend->AddEntry((TObject*)0, Form("Mean value: %.3f", hHFSum[iJetPt][iNJet]->GetMean()), "");
        legend->Draw();
      } // Number of jets loop
    } // Jet pT loop
  } // Draw 2D HF energy map

  // Draw the HF energy sum as a function of jet pT and number of jets
  if(drawHFSumGraph){

    // Define legend
    legend = new TLegend(0.48,0.34,0.75,0.68);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->AddEntry((TObject*)0, Form("%s  5.02 TeV", card->GetAlternativeDataType().Data()), "");

    // Loop over all number of jets bins
    for(int iNJet = 0; iNJet < nNJetBins; iNJet++){

      // Set a nice style for the graphs
      gAverageHFSum[iNJet]->SetMarkerStyle(markerStyle[iNJet]);
      gAverageHFSum[iNJet]->SetMarkerColor(color[iNJet]);
      gAverageHFSum[iNJet]->SetLineColor(color[iNJet]);

      // For the first graph create new canvas. Draw the others to the same canvas
      if(iNJet == 0){
        //drawer->DrawGraph(gAverageHFSum[0], comparedJetPtBin.at(0).first, comparedJetPtBin.at(comparedJetPtBin.size()-2).second, 35, 80, "Jet p_{T} (GeV)", "HF energy", "", "p");
        drawer->DrawGraph(gAverageHFSum[0], comparedJetPtBin.at(0).first, comparedJetPtBin.at(comparedJetPtBin.size()-2).second, 35, 105, "Jet p_{T} (GeV)", "HF energy", "", "p");
        legend->AddEntry(gAverageHFSum[0], Form("n_{jets > %.0f GeV} > 0", minJetPt), "p");
      } else {
        gAverageHFSum[iNJet]->Draw("p,same");
        legend->AddEntry(gAverageHFSum[iNJet], Form("n_{jets > %.0f GeV} = %d", minJetPt, iNJet), "p");
      }

    }

    // Draw a line showing the HF energy sum for the case when there are zero jets
    zeroJetsLine->Draw();
    estimatedShiftLine->Draw();
    estimatedShiftLineUp->Draw();
    estimatedShiftLineDown->Draw();
    legend->AddEntry(zeroJetsLine, Form("n_{jets > %.0f GeV} = 0", minJetPt), "l");

    // Draw the legend
    legend->Draw();

  } // Draw the HF energy sum as a function of jet pT and number of jets

}
