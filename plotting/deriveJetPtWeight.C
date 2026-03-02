#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Macro for comparing weighted and un-weighted Monte Carlo simulations with data
 */
void deriveJetPtWeight(){
  
  // Enumeration for data and MC
  enum enumDataType{kData, kMC, knDataTypes};
  
  // Data and MC files for the comparison
  TString fileName[knDataTypes];
  fileName[kData] = "data/pPb/pPb_5TeV_pToMinusEta_pfJets_eschemeAxis_nominalEnergyWeight_minimumBias_jetEtaMCcut_mixedConeHFshift28_processed_2025-06-30.root";
  fileName[kMC] = "data/pPb/pythiaEpos_8TeV_RecoReco_pToMinusEta_pfJets_eschemeAxis_nominalEnergyWeight_jetEtaCMcut_processed_2025-07-16.root";
  
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
  
  // ====================================================
  //                    Configuration
  // ====================================================
  
  // Figure saving
  const bool saveFigures = false;  // Save figures
  const char* saveComment = "";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures
  
  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager* histograms[knDataTypes];
  
  for(int iFile = 0; iFile < knDataTypes; iFile++){
    histograms[iFile] = new EECHistogramManager(inputFile[iFile]);
  }
  
  // Histograms for jet pT
  TH1D* hJetPt[knDataTypes];        // Jet pT distribution
  TH1D* hJetPtRatio;               // Jet pT ratio
  
  // Initialize the jet histograms to NULL
  for(int iDataType = 0; iDataType < knDataTypes; iDataType++){
    hJetPt[iDataType] = NULL;
  }
  hJetPtRatio = NULL;
  
  // Read the histograms from the file
  hJetPt[kData] = histograms[kData]->GetHistogramJetPt(0); // Jet pT distribution without centrality selection in data
  hJetPt[kMC] = histograms[kMC]->GetHistogramJetPt(0);     // Jet pT distribution without centrality selection in MC

  // Normalize all distributions to one and set line colors for histograms
  int color[knDataTypes] = {kRed, kBlue};
  int markerStyle[knDataTypes] = {kFullSquare, kFullCircle};

    // Normalize the jet pT distributions to some region and set drawing style
  for(int iDataType = 0; iDataType < knDataTypes; iDataType++){
    hJetPt[iDataType]->Scale(1.0 / histograms[iDataType]->GetJetPtIntegral(0,30,80));

    hJetPt[iDataType]->SetLineColor(color[iDataType]);
    hJetPt[iDataType]->SetMarkerColor(color[iDataType]);
    hJetPt[iDataType]->SetMarkerStyle(markerStyle[iDataType]);
  }

  // Calculate data to MC ratio for the jet pT
  hJetPtRatio = (TH1D*) hJetPt[kMC]->Clone("jetPtRatio");
  hJetPtRatio->Divide(hJetPt[kData]);

  // Fit a pol3 from 40 to 150 to the ratio
  hJetPtRatio->Fit("pol3","","",40,150);
  TF1* jetPtFit = hJetPtRatio->GetFunction("pol3");
  
  // ==========================================================================
  //        Draw the selected energy-energy correlator signal ratios
  // ==========================================================================
  
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  drawer->SetLogY(true);

  TLegend *legend;
  TString legendString[knDataTypes] = {"pPb data 5 TeV", "Pythia+EPOS 8 TeV"};
  
  // =========================
  // Draw the jet pT distributions
  // =========================
  
  // Create a legend for the figure
  legend = new TLegend(0.45,0.02,0.67,0.31);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  
  // Create a new canvas for the plot
  drawer->CreateSplitCanvas();
  
  // Draw the histograms to the upper canvas
  drawer->DrawHistogramToUpperPad(hJetPt[kData], "jet p_{T} (GeV)", "A.U.", " ");
  legend->AddEntry(hJetPt[kData], legendString[kData].Data(), "p");
  
  hJetPt[kMC]->Draw("same");
  legend->AddEntry(hJetPt[kMC], legendString[kMC].Data(), "p");
  
  // Draw the legends to the upper pad
  legend->Draw();
  
  // Draw the ratios to the lower pad
  drawer->SetLogY(false);
  drawer->SetGridY(true);
  hJetPtRatio->GetYaxis()->SetRangeUser(0.5,1.5);
  drawer->DrawHistogramToLowerPad(hJetPtRatio, "jet p_{T} (GeV)", "#frac{MC}{Data}", " ");
  drawer->SetGridY(false);
  
  // Save the figures to a file
  if(saveFigures){
    gPad->GetCanvas()->SaveAs(Form("figures/vzMonteCarloWeightCheck%s.%s", saveComment, figureFormat));
  }
  
}
