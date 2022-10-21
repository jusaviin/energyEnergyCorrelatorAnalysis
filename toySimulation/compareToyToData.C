#include "../plotting/JDrawer.h"

/*
 * Macro for comparing the results from toy simulation to only reflected cone part of the data
 */
void compareToyToData(){

  // Open the data and toy files
  TFile *dataFile = TFile::Open("../data/eecAnalysis_akFlowJets_updatedMultiplicityAndDensity_wtaAxis_preprocessed_2022-10-17.root");
  TFile *toyFile = TFile::Open("toySimulation100kevents2pSlope.root");
  
  // Read the correcponding energy-energy correlator histograms from the files
  TH1D *hEnergyEnergyCorrelatorData = (TH1D*) dataFile->Get("energyEnergyCorrelator/energyEnergyCorrelatorReflectedConePair_C0T0");
  TH1D *hEnergyEnergyCorrelatorToy = (TH1D*) toyFile->Get("energyEnergyCorrelator");
  TH1D *hParticleDensityDataRaw[2];
  hParticleDensityDataRaw[0] = (TH1D*) dataFile->Get("particleDensity/particleDensityReflectedCone_C0T0");
  hParticleDensityDataRaw[1] = (TH1D*) dataFile->Get("particlePtDensity/particlePtDensityReflectedCone_C0T0");
  TH1D *hParticleDensityToy[2];
  hParticleDensityToy[0] = (TH1D*) toyFile->Get("particleDensity");
  hParticleDensityToy[1] = (TH1D*) toyFile->Get("particlePtDensity");
  TH1D *hMultiplicityData = (TH1D*) dataFile->Get("multiplicityInReflectedCone/multiplicityInReflectedCone_C0T0");
  TH1D *hMultiplicityToy = (TH1D*) toyFile->Get("multiplicity");
  
  // Particle density in data has more bins, but same bin borders. For comparison with toy, trunkate the data histogram
  TH1D *hParticleDensityData[2];
  hParticleDensityData[0] = (TH1D*) hParticleDensityToy[0]->Clone("particleDensityForData");
  hParticleDensityData[1] = (TH1D*) hParticleDensityToy[1]->Clone("particlePtDensityForData");
  for(int iParticleDensity = 0; iParticleDensity < 2; iParticleDensity++){
    for(int iBin = 1; iBin <= hParticleDensityData[iParticleDensity]->GetNbinsX(); iBin++){
      hParticleDensityData[iParticleDensity]->SetBinContent(iBin, hParticleDensityDataRaw[iParticleDensity]->GetBinContent(iBin));
      hParticleDensityData[iParticleDensity]->SetBinError(iBin, hParticleDensityDataRaw[iParticleDensity]->GetBinError(iBin));
    }
  }
  
  // Normalize the histograms to the number of particle-particle pairs
  hEnergyEnergyCorrelatorData->Scale(1.0 / hEnergyEnergyCorrelatorData->Integral("width"));
  hEnergyEnergyCorrelatorToy->Scale(1.0 / hEnergyEnergyCorrelatorToy->Integral("width"));
  for(int iParticleDensity = 0; iParticleDensity < 2; iParticleDensity++){
    hParticleDensityData[iParticleDensity]->Scale(1.0 / hParticleDensityData[iParticleDensity]->Integral("width"));
    hParticleDensityToy[iParticleDensity]->Scale(1.0 / hParticleDensityToy[iParticleDensity]->Integral("width"));
  }
  hMultiplicityData->Scale(1.0 / hMultiplicityData->Integral("width"));
  hMultiplicityToy->Scale(1.0 / hMultiplicityToy->Integral("width"));
  
  // Calculate the ratio between data and toy
  TH1D* hToyToDataEnergyEnergyCorrelatorRatio = (TH1D*) hEnergyEnergyCorrelatorToy->Clone("energyEnergyCorrelatorRatio");
  hToyToDataEnergyEnergyCorrelatorRatio->Divide(hEnergyEnergyCorrelatorData);
  
  TH1D* hToyToDataParticleDensityRatio[2];
  for(int iParticleDensity = 0; iParticleDensity < 2; iParticleDensity++){
    hToyToDataParticleDensityRatio[iParticleDensity] = (TH1D*) hParticleDensityToy[iParticleDensity]->Clone(Form("particleDensityRatio%d", iParticleDensity));
    hToyToDataParticleDensityRatio[iParticleDensity]->Divide(hParticleDensityData[iParticleDensity]);
  }
  
  TH1D* hToyToDataMultiplicityRatio = (TH1D*) hMultiplicityToy->Clone("multiplicityrRatio");
  hToyToDataMultiplicityRatio->Divide(hMultiplicityData);
  
  // Configure a JDrawer to do the drawing for us
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  
  // ======================================
  // === Draw energy-energy correlators ===
  // ======================================
  
  // Create a legend for the figure
  TLegend *legend = new TLegend(0.58,0.08,0.83,0.44);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  legend->AddEntry((TObject*) 0, "Cent: 0-10%","");
  legend->AddEntry((TObject*) 0, "Jet p_{T} > 120 GeV","");
  legend->AddEntry((TObject*) 0, "Track p_{T} > 0.7 GeV","");
  
  // Create a canvas and draw the distributions and the ratio to the canvas
  drawer->CreateSplitCanvas();
  drawer->SetLogY(true);
  drawer->SetLogX(true);
  hEnergyEnergyCorrelatorData->SetLineColor(kBlack);
  hEnergyEnergyCorrelatorData->GetXaxis()->SetRangeUser(0.001,0.8);
  drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelatorData, "#Deltar", "EEC", " ");
  hEnergyEnergyCorrelatorToy->SetLineColor(kRed);
  hEnergyEnergyCorrelatorToy->Draw("same");
  
  // Draw the legend
  legend->AddEntry(hEnergyEnergyCorrelatorData, "PbPb, reflected cone", "l");
  legend->AddEntry(hEnergyEnergyCorrelatorToy, "Toy simulation", "l");
  legend->Draw();
  
  // For the ratio, do linear y-axis
  drawer->SetLogY(false);
  drawer->SetGridY(true);
  hToyToDataEnergyEnergyCorrelatorRatio->GetXaxis()->SetRangeUser(0.001,0.8);
  hToyToDataEnergyEnergyCorrelatorRatio->GetYaxis()->SetRangeUser(0.8,1.2);
  hToyToDataEnergyEnergyCorrelatorRatio->SetLineColor(kRed);
  drawer->DrawHistogramToLowerPad(hToyToDataEnergyEnergyCorrelatorRatio, "#Deltar", "#frac{Toy result}{Reflected cone}", " ");
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
    hParticleDensityToy[iParticleDensity]->SetLineColor(kRed);
    hParticleDensityToy[iParticleDensity]->Draw("same");
    
    // Draw the legend
    legend->AddEntry(hParticleDensityData[iParticleDensity], "PbPb, reflected cone", "l");
    legend->AddEntry(hParticleDensityToy[iParticleDensity], "Toy simulation", "l");
    legend->Draw();
    
    // For the ratio, do linear y-axis
    drawer->SetLogY(false);
    drawer->SetGridY(true);
    hToyToDataParticleDensityRatio[iParticleDensity]->GetYaxis()->SetRangeUser(0.8,1.2);
    hToyToDataParticleDensityRatio[iParticleDensity]->SetLineColor(kRed);
    drawer->DrawHistogramToLowerPad(hToyToDataParticleDensityRatio[iParticleDensity], "#Deltar", "#frac{Toy result}{Reflected cone}", " ");
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
  hMultiplicityToy->SetLineColor(kRed);
  hMultiplicityToy->Draw("same");
  
  // Draw the legend
  legend->AddEntry(hMultiplicityData, "PbPb, reflected cone", "l");
  legend->AddEntry(hMultiplicityToy, "Toy simulation", "l");
  legend->Draw();
  
  // For the ratio, do linear y-axis
  drawer->SetLogY(false);
  drawer->SetGridY(true);
  hToyToDataMultiplicityRatio->GetYaxis()->SetRangeUser(0,2);
  hToyToDataMultiplicityRatio->SetLineColor(kRed);
  drawer->DrawHistogramToLowerPad(hToyToDataMultiplicityRatio, "#Deltar", "#frac{Toy result}{Reflected cone}", " ");
  drawer->SetGridY(false);
}
