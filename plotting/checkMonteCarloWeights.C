#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Macro for comparing weighted and un-weighted Monte Carlo simulations with data
 */
void checkMonteCarloWeights(){
  
  // Enumeration for data and MC
  enum enumDataType{kData, kMC, knDataTypes};
  
  // Data and MC files for the comparison
  TString fileName[knDataTypes];
  fileName[kData] = "data/pPb_5TeV_pToMinusEta_pfJets_eschemeAxis_nominalEnergyWeight_minimumBias_jetEtaMCcut_HFEnergyNjet_processed_2025-06-20.root";
  fileName[kMC] = "data/pythiaEpos_8TeV_RecoReco_pToMinusEta_pfJets_eschemeAxis_nominalEnergyWeight_jetEtaCMcut_processed_2025-07-16.root";
  
  // Some PbPb files
  // data/eecAnalysis_akFlowJets_wtaAxis_cutBadPhi_miniAODtesting_processed_2023-01-30.root
  // data/eecAnalysis_akFlowJets_removeBadAcceptance_wtaAxis_processed_2022-10-25.root
  // data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_mAOD_4pC_wtaAxis_jetTrig_cutBadPhi_processed_2023-02-10.root
  // data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJet_updatedMultiplicityAndDensity_wtaAxis_noTrigger_preprocessed_2022-10-17.root
  
  // Some pp files
  // data/ppData_pfJets_wtaAxis_processed_2022-12-16.root
  // data/ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_processed_2023-01-09.root
  
  // Open the files and check that they exist
  TFile* inputFile[knDataTypes];
  for(int iFile = 0; iFile < knDataTypes; iFile++){
    
    inputFile[iFile] = TFile::Open(fileName[iFile]);
    
    if(inputFile[iFile] == NULL){
      cout << "Error! The file " << fileName[iFile].Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }
  }
  
  // Check if we are using PbPb or pp data
  EECCard *systemCard = new EECCard(inputFile[0]);
  TString collisionSystem = systemCard->GetDataType();
  bool isPbPbData = collisionSystem.Contains("PbPb");
  
  // ====================================================
  //                    Configuration
  // ====================================================
  
  // Figure saving
  const bool saveFigures = false;  // Save figures
  const char* saveComment = "";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures
  
  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager *histograms[knDataTypes];
  
  for(int iFile = 0; iFile < knDataTypes; iFile++){
    histograms[iFile] = new EECHistogramManager(inputFile[iFile]);
  }
  
  // Loaded event information histograms used to check the MC weighted
  TH1D* hVz[knDataTypes+1]; // Vz distribution
  TH1D* hVzRatio[2];        // Ratios of data to MC vz distributions
  TH1D* hCentrality[knDataTypes+1]; // Ratio between true and extracted signal
  TH1D* hCentralityRatio[2]; // Ratios of data to MC centrality distributions
  TH1D* hJetPt[knDataTypes]; // Jet pT distribution
  TH1D* hJetPtRatio[knDataTypes-1]; // Jet pT ratop
  
  // Initialize the event information histograms to NULL
  for(int iDataType = 0; iDataType < knDataTypes+1; iDataType++){
    hVz[iDataType] = NULL;
    hCentrality[iDataType] = NULL;
  } // Data type loop
  for(int iRatioType = 0; iRatioType < 2; iRatioType++){
    hVzRatio[iRatioType] = NULL;
    hCentralityRatio[iRatioType] = NULL;
  }
  for(int iDataType = 0; iDataType < knDataTypes; iDataType++){
    hJetPt[iDataType] = NULL;
    if(iDataType == knDataTypes-1) continue;
    hJetPtRatio[iDataType] = NULL;
  }
  
  // Read the histograms from the file
  hVz[0] = histograms[0]->GetHistogramVertexZWeighted(); // Vertex-z distribution from data (weighing for trigger combination)
  hVz[1] = histograms[1]->GetHistogramVertexZ();         // Vertex-z distribution from MC
  hVz[2] = histograms[1]->GetHistogramVertexZWeighted(); // Weighted vertex-z distribution from MC
  
  hCentrality[0] = histograms[0]->GetHistogramCentralityWeighted(); // Centrality distribution from data (weighing for trigger combination)
  hCentrality[1] = histograms[1]->GetHistogramCentrality();         // Centrality distribution from MC
  hCentrality[2] = histograms[1]->GetHistogramCentralityWeighted(); // Weighted centrality distribution from MC

  hJetPt[kData] = histograms[kData]->GetHistogramJetPt(0); // Jet pT distribution without centrality selection in data
  hJetPt[kMC] = histograms[kMC]->GetHistogramJetPt(0);     // Jet pT distribution without centrality selection in MC

  // Normalize all distributions to one and set line colors for histograms
  int color[knDataTypes+1] = {kBlack, kRed, kBlue};
  int markerStyle[knDataTypes+1] = {kFullSquare, kFullCircle, kFullCircle};
  for(int iDataType = 0; iDataType < knDataTypes+1; iDataType++){
    hVz[iDataType]->Scale(1.0 / hVz[iDataType]->Integral("width"));
    hCentrality[iDataType]->Scale(1.0 / hCentrality[iDataType]->Integral("width"));
    
    hVz[iDataType]->SetLineColor(color[iDataType]);
    hVz[iDataType]->SetMarkerColor(color[iDataType]);
    hVz[iDataType]->SetMarkerStyle(markerStyle[iDataType]);
    hCentrality[iDataType]->SetLineColor(color[iDataType]);
    hCentrality[iDataType]->SetMarkerColor(color[iDataType]);
    hCentrality[iDataType]->SetMarkerStyle(markerStyle[iDataType]);
  }
  
  // Calculate the MC to data ratios
  for(int iRatioType = 0; iRatioType < 2; iRatioType++){
    hVzRatio[iRatioType] = (TH1D*) hVz[iRatioType+1]->Clone(Form("vzRatio%d", iRatioType));
    hVzRatio[iRatioType]->Divide(hVz[0]);
        
    hCentralityRatio[iRatioType] = (TH1D*) hCentrality[iRatioType+1]->Clone(Form("centralityRatio%d", iRatioType));
    hCentralityRatio[iRatioType]->Divide(hCentrality[0]);
  }

    // Normalize the jet pT distributions to some region and set drawing style
  for(int iDataType = 0; iDataType < knDataTypes; iDataType++){
    hJetPt[iDataType]->Scale(1.0 / histograms[iDataType]->GetJetPtIntegral(0,30,80));

    hJetPt[iDataType]->SetLineColor(color[iDataType]);
    hJetPt[iDataType]->SetMarkerColor(color[iDataType]);
    hJetPt[iDataType]->SetMarkerStyle(markerStyle[iDataType]);
  }

  // Calculate data to MC ratio also for the jet pT
  for(int iDataType = 1; iDataType < knDataTypes; iDataType++){
    hJetPtRatio[iDataType-1] = (TH1D*) hJetPt[iDataType]->Clone(Form("vzRatio%d", iRatioType));
    hJetPtRatio[iDataType-1]->Divide(hJetPt[kData]);
  }
  
  // ==========================================================================
  //        Draw the selected energy-energy correlator signal ratios
  // ==========================================================================
  
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);

  TLegend *legend;
  TString legendString[knDataTypes+1] = {"PbPb data", "Raw P+H", "Weighted P+H"};
  
  // Different naming for distributions for pp
  if(!isPbPbData){
    legendString[0] = "pp data";
    legendString[1] = "Raw Pythia8";
    legendString[2] = "Weighted Pythia8";
  }
  
  // =========================
  // Draw the vz distributions
  // =========================
  
  // Create a legend for the figure
  legend = new TLegend(0.45+0.01*isPbPbData,0.02,0.67+0.01*isPbPbData,0.31);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  
  // Create a new canvas for the plot
  drawer->CreateSplitCanvas();
  
  // Draw the histograms to the upper canvas
  hVz[0]->GetYaxis()->SetRangeUser(-0.0015, 0.086);
  drawer->DrawHistogramToUpperPad(hVz[0], "v_{z} (cm)", "A.U.", " ");
  legend->AddEntry(hVz[0], legendString[0].Data(), "p");
  
  for(int iDataType = 1; iDataType < knDataTypes+1; iDataType++){
    hVz[iDataType]->Draw("same");
    legend->AddEntry(hVz[iDataType], legendString[iDataType].Data(), "p");
  }
  
  // Draw the legends to the upper pad
  legend->Draw();
  
  // Draw the ratios to the lower pad
  drawer->SetGridY(true);
  hVzRatio[0]->GetYaxis()->SetRangeUser(0.5,1.5);
  drawer->DrawHistogramToLowerPad(hVzRatio[0], "v_{z} (cm)", "#frac{MC}{Data}", " ");
  drawer->SetGridY(false);
  hVzRatio[1]->Draw("same");
  
  // Save the figures to a file
  if(saveFigures){
    gPad->GetCanvas()->SaveAs(Form("figures/vzMonteCarloWeightCheck%s.%s", saveComment, figureFormat));
  }
  
  // =================================
  // Draw the centrality distributions
  // =================================
  
  if(isPbPbData){
    
    // Create a legend for the figure
    legend = new TLegend(0.61,0.55,0.83,0.84);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    
    // Create a new canvas for the plot
    drawer->CreateSplitCanvas();
    
    // Draw the histograms to the upper canvas
    drawer->DrawHistogramToUpperPad(hCentrality[0], "Centrality (%)", "A.U.", " ");
    hCentrality[0]->GetYaxis()->SetRangeUser(-0.001,0.054);
    legend->AddEntry(hCentrality[0], legendString[0].Data(), "p");
    
    for(int iDataType = 1; iDataType < knDataTypes+1; iDataType++){
      hCentrality[iDataType]->Draw("same");
      legend->AddEntry(hCentrality[iDataType], legendString[iDataType].Data(), "p");
    }
    
    // Draw the legends to the upper pad
    legend->Draw();
    
    // Draw the ratios to the lower pad
    drawer->SetGridY(true);
    hCentralityRatio[0]->GetYaxis()->SetRangeUser(0.5,1.5);
    drawer->DrawHistogramToLowerPad(hCentralityRatio[0], "Centrality (%)", "#frac{MC}{Data}", " ");
    drawer->SetGridY(false);
    hCentralityRatio[1]->Draw("same");
    
    // Save the figures to a file
    if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/centralityMonteCarloWeightCheck%s.%s", saveComment, figureFormat));
    }
  }
}
