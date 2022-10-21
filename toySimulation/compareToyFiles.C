#include "../plotting/JDrawer.h"

/*
 * Macro for comparing the results from toy simulation to only reflected cone part of the data
 */
void compareToyFiles(){

  // Open the data and toy files
  const int nFiles = 2;
  TString fileNames[] = {"toySimulation100kevents.root", "toySimulation100keventsSignalMultiplicity.root", "toySimulation100kevents10pSlope.root"};
  
  TString legendComment[] = {"Signal multi", "Reflected multi", "10% sloped density"};
  
  int color[10] = {kBlack, kBlue, kRed, kGreen+2, kCyan, kMagenta, kOrange-1, kAzure+5, kOrange-2, kGray};
  
  // Open the toy simulation files
  TFile *inputFile[nFiles];
  for(int iFile = 0; iFile < nFiles; iFile++){
    inputFile[iFile] = TFile::Open(fileNames[iFile]);
  }
  
  // Read the energy-energy correlators and particle densities from the files
  TH1D *hEnergyEnergyCorrelator[nFiles];
  TH1D *hParticleDensity[nFiles][2];
  TH1D* hEnergyEnergyCorrelatorRatio[nFiles-1];
  TH1D* hParticleDensityRatio[nFiles-1][2];
  for(int iFile = 0; iFile < nFiles; iFile++){
    hEnergyEnergyCorrelator[iFile] = (TH1D*) inputFile[iFile]->Get("energyEnergyCorrelator");
    hParticleDensity[iFile][0] = (TH1D*) inputFile[iFile]->Get("particleDensity");
    hParticleDensity[iFile][1] = (TH1D*) inputFile[iFile]->Get("particlePtDensity");
    
    // Normalize the histograms to the number of particle-particle pairs
    hEnergyEnergyCorrelator[iFile]->Scale(1.0 / hEnergyEnergyCorrelator[iFile]->Integral("width"));
    for(int iParticleDensity = 0; iParticleDensity < 2; iParticleDensity++){
      hParticleDensity[iFile][iParticleDensity]->Scale(1.0 / hParticleDensity[iFile][iParticleDensity]->Integral("width"));
    }
    
    // Calculate the ratios to the first file
    if(iFile > 0){
      hEnergyEnergyCorrelatorRatio[iFile-1] = (TH1D*) hEnergyEnergyCorrelator[iFile]->Clone(Form("energyEnergyCorrelatorRatio%d",iFile));
      hEnergyEnergyCorrelatorRatio[iFile-1]->Divide(hEnergyEnergyCorrelator[0]);
      
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
    hEnergyEnergyCorrelatorRatio[iRatio]->GetYaxis()->SetRangeUser(0.8,1.2);
    hEnergyEnergyCorrelatorRatio[iRatio]->SetLineColor(color[iRatio+1]);
    if(iRatio == 0){
      drawer->DrawHistogramToLowerPad(hEnergyEnergyCorrelatorRatio[iRatio], "#Deltar", "#frac{Color}{Black}", " ");
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
    legend->AddEntry((TObject*) 0, "Jet p_{T} > 120 GeV","");
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
      hParticleDensityRatio[iRatio][iParticleDensity]->GetYaxis()->SetRangeUser(0.8,1.2);
      hParticleDensityRatio[iRatio][iParticleDensity]->SetLineColor(color[iRatio+1]);
      if(iRatio == 0){
        drawer->DrawHistogramToLowerPad(hParticleDensityRatio[iRatio][iParticleDensity], "#Deltar", "#frac{Color}{Black}", " ");
      } else {
        hParticleDensityRatio[iRatio][iParticleDensity]->Draw("same");
      }
    }
    drawer->SetGridY(false);
  }
  
}
