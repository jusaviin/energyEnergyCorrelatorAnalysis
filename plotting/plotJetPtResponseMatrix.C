#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Macro for finding background normalization scale in MC. This is done by integrating the energy-energy correlator distributions for background within the signal cone, and in the reflected cone, and taking a ratio of these numbers.
 */
void plotJetPtResponseMatrix(){

  // Enumeration for projection direction
  enum enumProjectionType{kConstantGen, kConstantReco, knProjectionTypes};

  // File from which the integrals are calculated
  TString inputFileName = "data/PbPbMC2018_RecoGen_akFlowJets_miniAOD_4pCentShift_noTrigger_jetPtClosure_newRatioBins_processed_2023-04-21.root";
  // ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_noCorrelations_jetPtResponseMatrixMoreBins_processed_2023-01-13.root
  // ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_noCorrelations_jetPtResponseMatrix_processed_2023-01-13.root
  // PbPbMC2018_RecoGen_akFlowJets_miniAOD_4pCentShift_noTrigger_jetPtClosure_newRatioBins_processed_2023-04-21.root
  // PbPbMC2018_GenGen_akFlowJets_miniAOD_4pCentShift_noTrigger_jetPtClosure_newRatioBins_processed_2023-04-21.root
  
  // Open the input file
  TFile* inputFile = TFile::Open(inputFileName);
  
  if(inputFile == NULL){
    cout << "Error! The file " << inputFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }
  
  // Load the card from the file and read the collision system
  EECCard* card = new EECCard(inputFile);
  TString collisionSystem = card->GetDataType();
  bool isPbPbData = collisionSystem.Contains("PbPb");
  
  // ====================================================
  //  Binning configuration for the integral calculation
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = (isPbPbData) ? card->GetNCentralityBins() : 1;
  const int nAnalysisJetPtBins = 4;
  
  // Default binning ranges for reference
  // centrality = {0,10,30,50,90}

  const double analysisJetPtBinBorders[nAnalysisJetPtBins+1] = {120,140,160,180,200};
    
  // Create and setup a new histogram manager to project and handle the histograms
  EECHistogramManager* histograms = new EECHistogramManager(inputFile,card);
  
  // Choose the energy-energy correlator types to load
  histograms->SetLoadJetPtResponseMatrix(true);
  
  // Choose the bin ranges
  histograms->SetCentralityBinRange(0,nCentralityBins);
  
  // Load the histograms from the file
  histograms->LoadProcessedHistograms();

  // Initialize the jet response matrices to NULL
  TH2D* hJetPtResponseMatrix[nCentralityBins];
  TH1D* hJetPtResponseMatrixProjection[knProjectionTypes][nCentralityBins][nAnalysisJetPtBins];

  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){

    hJetPtResponseMatrix[iCentrality] = NULL;

    for(int iJetPt = 0; iJetPt < nAnalysisJetPtBins; iJetPt++){
      for(int iProjectionType = 0; iProjectionType < knProjectionTypes; iProjectionType++){
        hJetPtResponseMatrixProjection[iProjectionType][iCentrality][iJetPt] = NULL;
      }
    }

  }  // Centrality loop

  // Read the jet response matrices from the file
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){

    // Read the two dimensional histograms from the file
    hJetPtResponseMatrix[iCentrality] = histograms->GetHistogramJetPtResponseMatrix(iCentrality);

    // Do the projections from the two-dimensional histogram
    for(int iJetPt = 0; iJetPt < nAnalysisJetPtBins; iJetPt++){

      // First project response for constant generator level jet pT
      hJetPtResponseMatrixProjection[kConstantGen][iCentrality][iJetPt] = hJetPtResponseMatrix[iCentrality]->ProjectionX(Form("projX%d%d", iCentrality, iJetPt), 4+iJetPt, 4+iJetPt);

      // Normilize the distribution to one to make it probability distribution
      hJetPtResponseMatrixProjection[kConstantGen][iCentrality][iJetPt]->Scale(1.0 / hJetPtResponseMatrixProjection[kConstantGen][iCentrality][iJetPt]->Integral());

      // Then project the response for constant reconstructed jet pT
      hJetPtResponseMatrixProjection[kConstantReco][iCentrality][iJetPt] = hJetPtResponseMatrix[iCentrality]->ProjectionY(Form("projY%d%d", iCentrality, iJetPt), 4+iJetPt, 4+iJetPt);

      // Normilize the distribution to one to make it probability distribution
      hJetPtResponseMatrixProjection[kConstantReco][iCentrality][iJetPt]->Scale(1.0 / hJetPtResponseMatrixProjection[kConstantReco][iCentrality][iJetPt]->Integral());
    } // Jet pT loop 

  }  // Centrality loop

  // **********************************
  //         Draw the figures
  // **********************************

  // Prepare a JDrawer for drawing purposes
  JDrawer* drawer = new JDrawer();
  drawer->SetLogZ(true);
  drawer->SetLeftMargin(0.13);
  drawer->SetRightMargin(0.11);
  drawer->SetTopMargin(0.08);
  drawer->SetTitleOffsetX(1.17);
  drawer->SetTitleOffsetY(1);

  // Helper variables
  TLegend* legend;
  TString centralityString;
  TString compactCentralityString;
  TString jetPtString;
  TString compactJetPtString;

  // Draw the response matrices to canvas
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){

    if(isPbPbData){
      centralityString = Form("Pythia+Hydjet: %.0f-%.0f", card->GetLowBinBorderCentrality(iCentrality), card->GetHighBinBorderCentrality(iCentrality));
      compactCentralityString = Form("_C%.0f-%.0f", card->GetLowBinBorderCentrality(iCentrality), card->GetHighBinBorderCentrality(iCentrality));
    } else {
      centralityString = "Pythia8";
      compactCentralityString = "_pythia8";
    }

    hJetPtResponseMatrix[iCentrality]->GetZaxis()->SetRangeUser(0.01,3000);
    drawer->DrawHistogram(hJetPtResponseMatrix[iCentrality],"Reconstructed p_{T}", "Generator level p_{T}", Form("Jet p_{T} response, %s", centralityString.Data()), "colz");

  } // Centrality loop

  // Reset the drawing style for 1-dimensional histograms
  drawer->Reset();

  // Draw the projections of the response matrices to the canvas
  TString jetTypeStringCompact[knProjectionTypes] = {"Gen", "Reco"};
  TString jetTypeString[knProjectionTypes] = {"Reconstructed", "Generator"};
  for(int iProjectionType = 0; iProjectionType < knProjectionTypes; iProjectionType++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      if(isPbPbData){
        centralityString = Form("Centrality: %.0f-%.0f", card->GetLowBinBorderCentrality(iCentrality), card->GetHighBinBorderCentrality(iCentrality));
        compactCentralityString = Form("_C%.0f-%.0f", card->GetLowBinBorderCentrality(iCentrality), card->GetHighBinBorderCentrality(iCentrality));
      } else {
        centralityString = "Pythia8";
        compactCentralityString = "_pythia8";
      }

      for(int iJetPt = 0; iJetPt < nAnalysisJetPtBins; iJetPt++){
        jetPtString = Form("%.0f < %s jet p_{T} < %.0f GeV", analysisJetPtBinBorders[iJetPt], jetTypeStringCompact[iProjectionType].Data(), analysisJetPtBinBorders[iJetPt + 1]);
        compactJetPtString = Form("_T%.0f-%.0f", analysisJetPtBinBorders[iJetPt], analysisJetPtBinBorders[iJetPt + 1]);

        drawer->DrawHistogram(hJetPtResponseMatrixProjection[iProjectionType][iCentrality][iJetPt], Form("%s p_{T}", jetTypeString[iProjectionType].Data()), "Probability", " ", "");

        // Create a legend for the figure
        legend = new TLegend(0.42, 0.7, 0.8, 0.88);
        legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);

        legend->AddEntry((TObject*)0, centralityString.Data(), "");
        legend->AddEntry((TObject*)0, jetPtString.Data(), "");

        legend->Draw();
      }
    }  // Centrality loop
  }    // Projection type loop
}
