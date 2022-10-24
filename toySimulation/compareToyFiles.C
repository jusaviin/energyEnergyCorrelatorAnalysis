#include "../plotting/JDrawer.h"

/*
 * Macro for comparing the results from toy simulation to only reflected cone part of the data
 */
void compareToyFiles(){

  // Open the data and toy files
  const int nFiles = 4;
  TString fileNames[] = {"toySimulation100kevents.root", "toySimulation100kevents10pFlow.root", "toySimulation100kevents10pFlow30pBias.root", "toySimulation100kevents10pFlowMaxBiasSlope.root"};
  
  TString legendComment[] = {"Nominal", "10% v_{2}", "10% v_{2}, 30% bias", "10% v_{2}, max bias"};
  TString ratioName = "#frac{Flow}{Nominal}";
  
  int color[10] = {kBlack, kBlue, kRed, kGreen+2, kCyan, kMagenta, kOrange-1, kAzure+5, kOrange-2, kGray};
  
  // Open the toy simulation files
  TFile *inputFile[nFiles];
  for(int iFile = 0; iFile < nFiles; iFile++){
    inputFile[iFile] = TFile::Open(fileNames[iFile]);
  }
  
  // Read the energy-energy correlators and particle densities from the files
  TH1D *hEnergyEnergyCorrelator[nFiles];
  TH1D *hParticleDensity[nFiles][2];
  TH1D *hMultiplicity[nFiles];
  TH1D *hTrackPt[nFiles];
  TH1D *hEnergyEnergyCorrelatorRatio[nFiles-1];
  TH1D *hParticleDensityRatio[nFiles-1][2];
  TH1D *hMultiplicityRatio[nFiles-1];
  TH1D *hTrackPtRatio[nFiles-1];
  for(int iFile = 0; iFile < nFiles; iFile++){
    hEnergyEnergyCorrelator[iFile] = (TH1D*) inputFile[iFile]->Get("energyEnergyCorrelator");
    hParticleDensity[iFile][0] = (TH1D*) inputFile[iFile]->Get("particleDensity");
    hParticleDensity[iFile][1] = (TH1D*) inputFile[iFile]->Get("particlePtDensity");
    hMultiplicity[iFile] = (TH1D*) inputFile[iFile]->Get("multiplicity");
    hTrackPt[iFile] = (TH1D*) inputFile[iFile]->Get("trackPt");
    
    // Normalize the histograms to the number of particle-particle pairs
    hEnergyEnergyCorrelator[iFile]->Scale(1.0 / hEnergyEnergyCorrelator[iFile]->Integral("width"));
    hMultiplicity[iFile]->Scale(1.0 / hMultiplicity[iFile]->Integral("width"));
    hTrackPt[iFile]->Scale(1.0 / hTrackPt[iFile]->Integral("width"));
    for(int iParticleDensity = 0; iParticleDensity < 2; iParticleDensity++){
      hParticleDensity[iFile][iParticleDensity]->Scale(1.0 / hParticleDensity[iFile][iParticleDensity]->Integral("width"));
    }
    
    // Calculate the ratios to the first file
    if(iFile > 0){
      hEnergyEnergyCorrelatorRatio[iFile-1] = (TH1D*) hEnergyEnergyCorrelator[iFile]->Clone(Form("energyEnergyCorrelatorRatio%d",iFile));
      hEnergyEnergyCorrelatorRatio[iFile-1]->Divide(hEnergyEnergyCorrelator[0]);
      
      hMultiplicityRatio[iFile-1] = (TH1D*) hMultiplicity[iFile]->Clone(Form("multiplicityRatio%d",iFile));
      hMultiplicityRatio[iFile-1]->Divide(hMultiplicity[0]);
      
      hTrackPtRatio[iFile-1] = (TH1D*) hTrackPt[iFile]->Clone(Form("trackPtRatio%d",iFile));
      hTrackPtRatio[iFile-1]->Divide(hTrackPt[0]);
      
      for(int iParticleDensity = 0; iParticleDensity < 2; iParticleDensity++){
        hParticleDensityRatio[iFile-1][iParticleDensity] = (TH1D*) hParticleDensity[iFile][iParticleDensity]->Clone(Form("particleDensityRatio%d%d", iFile, iParticleDensity));
        hParticleDensityRatio[iFile-1][iParticleDensity]->Divide(hParticleDensity[0][iParticleDensity]);
      }
    }
    
  }
  
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
  legend->AddEntry((TObject*) 0, "Toy simulation","");
  legend->AddEntry((TObject*) 0, "Track p_{T} > 0.7 GeV","");
  
  // Create a canvas and draw the distributions and the ratio to the canvas
  drawer->CreateSplitCanvas();
  drawer->SetLogY(true);
  drawer->SetLogX(true);
  hEnergyEnergyCorrelator[0]->SetLineColor(color[0]);
  hEnergyEnergyCorrelator[0]->GetXaxis()->SetRangeUser(0.001,0.8);
  drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelator[0], "#Deltar", "EEC", " ");
  for(int iFile = 1; iFile < nFiles; iFile++){
    hEnergyEnergyCorrelator[iFile]->SetLineColor(color[iFile]);
    hEnergyEnergyCorrelator[iFile]->Draw("same");
  }
  
  // Draw the legend
  for(int iFile = 0; iFile < nFiles; iFile++){
    legend->AddEntry(hEnergyEnergyCorrelator[iFile], legendComment[iFile].Data(), "l");
  }
  legend->Draw();
  
  // For the ratio, do linear y-axis
  drawer->SetLogY(false);
  drawer->SetGridY(true);
  for(int iRatio = 0; iRatio < nFiles-1; iRatio++){
    hEnergyEnergyCorrelatorRatio[iRatio]->GetXaxis()->SetRangeUser(0.001,0.8);
    hEnergyEnergyCorrelatorRatio[iRatio]->GetYaxis()->SetRangeUser(0.9,1.1);
    hEnergyEnergyCorrelatorRatio[iRatio]->SetLineColor(color[iRatio+1]);
    if(iRatio == 0){
      drawer->DrawHistogramToLowerPad(hEnergyEnergyCorrelatorRatio[iRatio], "#Deltar", ratioName.Data(), " ");
    } else {
      hEnergyEnergyCorrelatorRatio[iRatio]->Draw("same");
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
    legend->AddEntry((TObject*) 0, "Track p_{T} > 0.7 GeV","");
    
    // Create a canvas and draw the distributions and the ratio to the canvas
    drawer->CreateSplitCanvas();
    drawer->SetLogY(false);
    drawer->SetLogX(false);
    hParticleDensity[0][iParticleDensity]->SetLineColor(color[0]);
    drawer->DrawHistogramToUpperPad(hParticleDensity[0][iParticleDensity], "#Deltar", densityAxisString[iParticleDensity], " ");
    for(int iFile = 1; iFile < nFiles; iFile++){
      hParticleDensity[iFile][iParticleDensity]->SetLineColor(color[iFile]);
      hParticleDensity[iFile][iParticleDensity]->Draw("same");
    }
    
    // Draw the legend
    for(int iFile = 0; iFile < nFiles; iFile++){
      legend->AddEntry(hParticleDensity[iFile][iParticleDensity], legendComment[iFile].Data(), "l");
    }
    legend->Draw();
    
    // For the ratio, do linear y-axis
    drawer->SetLogY(false);
    drawer->SetGridY(true);
    for(int iRatio = 0; iRatio < nFiles-1; iRatio++){
      hParticleDensityRatio[iRatio][iParticleDensity]->GetYaxis()->SetRangeUser(0.9,1.1);
      hParticleDensityRatio[iRatio][iParticleDensity]->SetLineColor(color[iRatio+1]);
      if(iRatio == 0){
        drawer->DrawHistogramToLowerPad(hParticleDensityRatio[iRatio][iParticleDensity], "#Deltar", ratioName.Data(), " ");
      } else {
        hParticleDensityRatio[iRatio][iParticleDensity]->Draw("same");
      }
    }
    drawer->SetGridY(false);
  }
  
  // =======================================
  // ======    Draw multiplicities    ======
  // =======================================
  
  // Create a legend for the figure
  legend = new TLegend(0.58,0.48,0.83,0.84);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  legend->AddEntry((TObject*) 0, "Toy simulation","");
  legend->AddEntry((TObject*) 0, "Track p_{T} > 0.7 GeV","");
  
  // Create a canvas and draw the distributions and the ratio to the canvas
  drawer->CreateSplitCanvas();
  hMultiplicity[0]->SetLineColor(color[0]);
  drawer->DrawHistogramToUpperPad(hMultiplicity[0], "Multiplicity", "Counts", " ");
  for(int iFile = 1; iFile < nFiles; iFile++){
    hMultiplicity[iFile]->SetLineColor(color[iFile]);
    hMultiplicity[iFile]->Draw("same");
  }
  
  // Draw the legend
  for(int iFile = 0; iFile < nFiles; iFile++){
    legend->AddEntry(hMultiplicity[iFile], legendComment[iFile].Data(), "l");
  }
  legend->Draw();
  
  // For the ratio, do linear y-axis
  drawer->SetLogY(false);
  drawer->SetGridY(true);
  for(int iRatio = 0; iRatio < nFiles-1; iRatio++){
    hMultiplicityRatio[iRatio]->GetYaxis()->SetRangeUser(0.8,1.2);
    hMultiplicityRatio[iRatio]->SetLineColor(color[iRatio+1]);
    if(iRatio == 0){
      drawer->DrawHistogramToLowerPad(hMultiplicityRatio[iRatio], "Multiplicity", ratioName.Data(), " ");
    } else {
      hMultiplicityRatio[iRatio]->Draw("same");
    }
  }
  drawer->SetGridY(false);
  
  // =================================
  // ======    Draw track pT    ======
  // =================================
  
  // Create a legend for the figure
  legend = new TLegend(0.58,0.48,0.83,0.84);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  legend->AddEntry((TObject*) 0, "Toy simulation","");
  legend->AddEntry((TObject*) 0, "Track p_{T} > 0.7 GeV","");
  
  // Create a canvas and draw the distributions and the ratio to the canvas
  drawer->CreateSplitCanvas();
  drawer->SetLogY(true);
  hTrackPt[0]->SetLineColor(color[0]);
  drawer->DrawHistogramToUpperPad(hTrackPt[0], "Particle p_{T} (GeV)", "Counts", " ");
  for(int iFile = 1; iFile < nFiles; iFile++){
    hTrackPt[iFile]->SetLineColor(color[iFile]);
    hTrackPt[iFile]->Draw("same");
  }
  
  // Draw the legend
  for(int iFile = 0; iFile < nFiles; iFile++){
    legend->AddEntry(hTrackPt[iFile], legendComment[iFile].Data(), "l");
  }
  legend->Draw();
  
  // For the ratio, do linear y-axis
  drawer->SetLogY(false);
  drawer->SetGridY(true);
  for(int iRatio = 0; iRatio < nFiles-1; iRatio++){
    hTrackPtRatio[iRatio]->GetYaxis()->SetRangeUser(0.8,1.2);
    hTrackPtRatio[iRatio]->SetLineColor(color[iRatio+1]);
    if(iRatio == 0){
      drawer->DrawHistogramToLowerPad(hTrackPtRatio[iRatio], "Particle p_{T} (GeV)", ratioName.Data(), " ");
    } else {
      hTrackPtRatio[iRatio]->Draw("same");
    }
  }
  drawer->SetGridY(false);
  
}
