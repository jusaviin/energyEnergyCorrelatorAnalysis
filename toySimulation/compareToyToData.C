#include "../plotting/JDrawer.h"

/*
 * Macro for comparing the results from toy simulation to only reflected cone part of the data
 */
void compareToyToData(){
  
  // Open the data and toy files
  TFile *dataFile = TFile::Open("../data/PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_finalMcWeight_processed_2023-03-08.root");
  // ../data/eecAnalysis_akFlowJets_updatedMultiplicityAndDensity_wtaAxis_preprocessed_2022-10-17.root
  // eecAnalysis_akFlowJets_removeBadAcceptance_wtaAxis_processed_2022-10-25.root
  // PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_finalMcWeight_processed_2023-03-08.root
  
  const int nToyFiles = 1;
  TString toyFileNames[] = {"toySimulation100kevents.root","toySimulation100kevents10pFlow.root"};
  TFile *toyFile[nToyFiles];
  for(int iFile = 0; iFile < nToyFiles; iFile++){
    toyFile[iFile] = TFile::Open(toyFileNames[iFile]);
  }
  
  TString toyLegendComment[] = {"Nominal", "10% v_{2}"};
  
  
  // Read the correcponding energy-energy correlator histograms from the files
  const int nDataCorrelators = 2;
  TH1D *hEnergyEnergyCorrelatorData[nDataCorrelators];
  hEnergyEnergyCorrelatorData[0] = (TH1D*) dataFile->Get("energyEnergyCorrelator/energyEnergyCorrelatorSignalReflectedConePair_C0T5S3");
  hEnergyEnergyCorrelatorData[1] = (TH1D*) dataFile->Get("energyEnergyCorrelator/energyEnergyCorrelatorReflectedConePair_C0T5S3");
  TH1D *hParticleDensityDataRaw[2];
  hParticleDensityDataRaw[0] = (TH1D*) dataFile->Get("particleDensity/particleDensityReflectedCone_C0T0");
  hParticleDensityDataRaw[1] = (TH1D*) dataFile->Get("particlePtDensity/particlePtDensityReflectedCone_C0T0");
  TH1D *hMultiplicityData = (TH1D*) dataFile->Get("multiplicityInReflectedCone/multiplicityInReflectedCone_C0T0");
  TH1D *hTrackPtData = (TH1D*) dataFile->Get("track/trackPt_C0");
  
  TH1D *hEnergyEnergyCorrelatorToy[nToyFiles];
  TH1D *hParticleDensityToy[nToyFiles][2];
  TH1D *hMultiplicityToy[nToyFiles];
  TH1D *hTrackPtToy[nToyFiles];
  for(int iFile = 0; iFile < nToyFiles; iFile++){
    hEnergyEnergyCorrelatorToy[iFile] = (TH1D*) toyFile[iFile]->Get("energyEnergyCorrelator");
    hParticleDensityToy[iFile][0] = (TH1D*) toyFile[iFile]->Get("particleDensity");
    hParticleDensityToy[iFile][1] = (TH1D*) toyFile[iFile]->Get("particlePtDensity");
    hMultiplicityToy[iFile] = (TH1D*) toyFile[iFile]->Get("multiplicity");
    hTrackPtToy[iFile] = (TH1D*) toyFile[iFile]->Get("trackPt");
  }
  
  
  
  // Particle density in data has more bins, but same bin borders. For comparison with toy, trunkate the data histogram
  TH1D *hParticleDensityData[2];
  hParticleDensityData[0] = (TH1D*) hParticleDensityToy[0][0]->Clone("particleDensityForData");
  hParticleDensityData[1] = (TH1D*) hParticleDensityToy[0][1]->Clone("particlePtDensityForData");
  for(int iParticleDensity = 0; iParticleDensity < 2; iParticleDensity++){
    for(int iBin = 1; iBin <= hParticleDensityData[iParticleDensity]->GetNbinsX(); iBin++){
      hParticleDensityData[iParticleDensity]->SetBinContent(iBin, hParticleDensityDataRaw[iParticleDensity]->GetBinContent(iBin));
      hParticleDensityData[iParticleDensity]->SetBinError(iBin, hParticleDensityDataRaw[iParticleDensity]->GetBinError(iBin));
    }
  }
  
  // Normalize the histograms to the number of particle-particle pairs
  for(int iEnergyEnergyCorrelatorData = 0; iEnergyEnergyCorrelatorData < nDataCorrelators; iEnergyEnergyCorrelatorData++){
    hEnergyEnergyCorrelatorData[iEnergyEnergyCorrelatorData]->Scale(1.0 / hEnergyEnergyCorrelatorData[iEnergyEnergyCorrelatorData]->Integral("width"));
  }
  for(int iParticleDensity = 0; iParticleDensity < 2; iParticleDensity++){
    hParticleDensityData[iParticleDensity]->Scale(1.0 / hParticleDensityData[iParticleDensity]->Integral("width"));
  }
  hMultiplicityData->Scale(1.0 / hMultiplicityData->Integral("width"));
  hTrackPtData->Scale(1.0 / hTrackPtData->Integral("width"));
  
  for(int iFile = 0; iFile < nToyFiles; iFile++){
    hEnergyEnergyCorrelatorToy[iFile]->Scale(1.0 / hEnergyEnergyCorrelatorToy[iFile]->Integral("width"));
    for(int iParticleDensity = 0; iParticleDensity < 2; iParticleDensity++){
      hParticleDensityToy[iFile][iParticleDensity]->Scale(1.0 / hParticleDensityToy[iFile][iParticleDensity]->Integral("width"));
    }
    
    hMultiplicityToy[iFile]->Scale(1.0 / hMultiplicityToy[iFile]->Integral("width"));
    hTrackPtToy[iFile]->Scale(1.0 / hTrackPtToy[iFile]->Integral("width"));
  }
  
  
  
  // Calculate the ratio between data and toy
  TH1D* hToyToDataEnergyEnergyCorrelatorRatio[nDataCorrelators];
  TH1D* hToyToDataParticleDensityRatio[nToyFiles][2];
  TH1D* hToyToDataMultiplicityRatio[nToyFiles];
  TH1D* hToyToDataTrackPtRatio[nToyFiles];
  
  for(int iDataCorrelator = 0; iDataCorrelator < nDataCorrelators; iDataCorrelator++){
    hToyToDataEnergyEnergyCorrelatorRatio[iDataCorrelator] = (TH1D*) hEnergyEnergyCorrelatorToy[0]->Clone(Form("energyEnergyCorrelatorRatio%d", iDataCorrelator));
    hToyToDataEnergyEnergyCorrelatorRatio[iDataCorrelator]->Divide(hEnergyEnergyCorrelatorData[iDataCorrelator]);
  }
  
  for(int iFile = 0; iFile < nToyFiles; iFile++){
    //hToyToDataEnergyEnergyCorrelatorRatio[iFile] = (TH1D*) hEnergyEnergyCorrelatorToy[iFile]->Clone(Form("energyEnergyCorrelatorRatio%d", iFile));
    //hToyToDataEnergyEnergyCorrelatorRatio[iFile]->Divide(hEnergyEnergyCorrelatorData[0]);
    
    for(int iParticleDensity = 0; iParticleDensity < 2; iParticleDensity++){
      hToyToDataParticleDensityRatio[iFile][iParticleDensity] = (TH1D*) hParticleDensityToy[iFile][iParticleDensity]->Clone(Form("particleDensityRatio%d%d", iFile, iParticleDensity));
      hToyToDataParticleDensityRatio[iFile][iParticleDensity]->Divide(hParticleDensityData[iParticleDensity]);
    }
    
    hToyToDataMultiplicityRatio[iFile] = (TH1D*) hMultiplicityToy[iFile]->Clone(Form("multiplicityRatio%d", iFile));
    hToyToDataMultiplicityRatio[iFile]->Divide(hMultiplicityData);
    
    hToyToDataTrackPtRatio[iFile] = (TH1D*) hTrackPtToy[iFile]->Clone(Form("trackPtRatio%d", iFile));
    hToyToDataTrackPtRatio[iFile]->Divide(hTrackPtData);
  }
  
  // Configure a JDrawer to do the drawing for us
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  
  int color[9] = {kRed, kBlue, kGreen+2, kCyan, kMagenta, kOrange-1, kAzure+5, kOrange-2, kGray};
  int dataColor[3] = {kBlack, kBlue, kRed};
  
  // ======================================
  // === Draw energy-energy correlators ===
  // ======================================
  
  // Create a legend for the figure
  TLegend *legend = new TLegend(0.58,0.08,0.83,0.44);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  legend->AddEntry((TObject*) 0, "Cent: 0-10%","");
  legend->AddEntry((TObject*) 0, "Jet p_{T} > 120 GeV","");
  legend->AddEntry((TObject*) 0, "Track p_{T} > 3 GeV","");
  //legend->AddEntry((TObject*) 0, "Pythia+Hydjet, GenGen","");
  
  // Create a canvas and draw the distributions and the ratio to the canvas
  drawer->CreateSplitCanvas();
  drawer->SetLogY(true);
  drawer->SetLogX(true);
  hEnergyEnergyCorrelatorData[0]->SetLineColor(kBlack);
  hEnergyEnergyCorrelatorData[0]->GetXaxis()->SetRangeUser(0.001,0.8);
  drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelatorData[0], "#Deltar", "EEC", " ");
  legend->AddEntry(hEnergyEnergyCorrelatorData[0], "Jet+ref, Hydjet", "l");
  
  hEnergyEnergyCorrelatorData[1]->SetLineColor(kBlue);
  hEnergyEnergyCorrelatorData[1]->Draw("same");
  legend->AddEntry(hEnergyEnergyCorrelatorData[1], "Ref+ref, Hydjet", "l");
  
  for(int iFile = 0; iFile < nToyFiles; iFile++){
    hEnergyEnergyCorrelatorToy[iFile]->SetLineColor(color[iFile]);
    hEnergyEnergyCorrelatorToy[iFile]->Draw("same");
    legend->AddEntry(hEnergyEnergyCorrelatorToy[iFile], Form("Toy, %s", toyLegendComment[iFile].Data()), "l");
  }
  
  // Draw the legend
  legend->Draw();
  
  // For the ratio, do linear y-axis
  drawer->SetLogY(false);
  drawer->SetGridY(true);
  for(int iFile = 0; iFile < nDataCorrelators; iFile++){
    hToyToDataEnergyEnergyCorrelatorRatio[iFile]->GetXaxis()->SetRangeUser(0.001,0.8);
    hToyToDataEnergyEnergyCorrelatorRatio[iFile]->GetYaxis()->SetRangeUser(0,2);
    hToyToDataEnergyEnergyCorrelatorRatio[iFile]->SetLineColor(dataColor[iFile]);
    if(iFile == 0){
      drawer->DrawHistogramToLowerPad(hToyToDataEnergyEnergyCorrelatorRatio[iFile], "#Deltar", "#frac{Toy result}{Reflected cone}", " ");
    } else {
      hToyToDataEnergyEnergyCorrelatorRatio[iFile]->Draw("same");
    }
  }
  drawer->SetGridY(false);

  // =====================================
  // ====== Draw particle densities ======
  // =====================================
  const char* densityAxisString[2] = {"Particle density", "p_{T} density"};
  
  for(int iParticleDensity = 0; iParticleDensity < 2; iParticleDensity++){
    legend = new TLegend(0.58,0.08,0.83,0.44);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->AddEntry((TObject*) 0, "Cent: 0-10%","");
    legend->AddEntry((TObject*) 0, "Jet p_{T} > 120 GeV","");
    legend->AddEntry((TObject*) 0, "Track p_{T} > 0.7 GeV","");
    
    // Create a canvas and draw the distributions and the ratio to the canvas
    drawer->CreateSplitCanvas();
    drawer->SetLogY(false);
    drawer->SetLogX(false);
    hParticleDensityData[iParticleDensity]->SetLineColor(kBlack);
    drawer->DrawHistogramToUpperPad(hParticleDensityData[iParticleDensity], "#Deltar", densityAxisString[iParticleDensity], " ");
    legend->AddEntry(hParticleDensityData[iParticleDensity], "PbPb, reflected cone", "l");
    
    for(int iFile = 0; iFile < nToyFiles; iFile++){
      hParticleDensityToy[iFile][iParticleDensity]->SetLineColor(color[iFile]);
      hParticleDensityToy[iFile][iParticleDensity]->Draw("same");
      legend->AddEntry(hParticleDensityToy[iFile][iParticleDensity], Form("Toy, %s", toyLegendComment[iFile].Data()), "l");
    }
    
    // Draw the legend
    legend->Draw();
    
    // For the ratio, do linear y-axis
    drawer->SetLogY(false);
    drawer->SetGridY(true);
    for(int iFile = 0; iFile < nToyFiles; iFile++){
      hToyToDataParticleDensityRatio[iFile][iParticleDensity]->GetYaxis()->SetRangeUser(0.8,1.2);
      hToyToDataParticleDensityRatio[iFile][iParticleDensity]->SetLineColor(color[iFile]);
      if(iFile == 0){
        drawer->DrawHistogramToLowerPad(hToyToDataParticleDensityRatio[iFile][iParticleDensity], "#Deltar", "#frac{Toy result}{Reflected cone}", " ");
      } else {
        hToyToDataParticleDensityRatio[iFile][iParticleDensity]->Draw("same");
      }
    }
    drawer->SetGridY(false);
  }
  
  // =======================================
  // === Draw multiplicities in jet cone ===
  // =======================================
  
  // Create a legend for the figure
  legend = new TLegend(0.58,0.08,0.83,0.44);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  legend->AddEntry((TObject*) 0, "Cent: 0-10%","");
  legend->AddEntry((TObject*) 0, "Jet p_{T} > 120 GeV","");
  legend->AddEntry((TObject*) 0, "Track p_{T} > 0.7 GeV","");
  
  // Create a canvas and draw the distributions and the ratio to the canvas
  drawer->CreateSplitCanvas();
  drawer->SetLogY(false);
  drawer->SetLogX(false);
  hMultiplicityData->SetLineColor(kBlack);
  drawer->DrawHistogramToUpperPad(hMultiplicityData, "Multiplicity", "Counts", " ");
  legend->AddEntry(hMultiplicityData, "PbPb, reflected cone", "l");
  
  for(int iFile = 0; iFile < nToyFiles; iFile++){
    hMultiplicityToy[iFile]->SetLineColor(color[iFile]);
    hMultiplicityToy[iFile]->Draw("same");
    legend->AddEntry(hMultiplicityToy[iFile], Form("Toy, %s", toyLegendComment[iFile].Data()), "l");
  }
  
  // Draw the legend
  legend->Draw();
  
  // For the ratio, do linear y-axis
  drawer->SetLogY(false);
  drawer->SetGridY(true);
  for(int iFile = 0; iFile < nToyFiles; iFile++){
    hToyToDataMultiplicityRatio[iFile]->GetYaxis()->SetRangeUser(0,2);
    hToyToDataMultiplicityRatio[iFile]->SetLineColor(color[iFile]);
    if(iFile == 0){
      drawer->DrawHistogramToLowerPad(hToyToDataMultiplicityRatio[iFile], "#Deltar", "#frac{Toy result}{Reflected cone}", " ");
    } else {
      hToyToDataMultiplicityRatio[iFile]->Draw("same");
    }
  }
  drawer->SetGridY(false);
}
