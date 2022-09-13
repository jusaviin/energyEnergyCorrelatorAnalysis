/*
 * Implementation of EECDrawer
 */

// Root includes
#include <TPad.h>

// Own includes
#include "EECDrawer.h"

/*
 * Constructor
 */
EECDrawer::EECDrawer(EECHistogramManager *inputHistograms) :
  fHistograms(inputHistograms),
  fFigureSaveNameAppend(""),
  fDrawEventInformation(false),
  fDrawJets(false),
  fDrawIndividualEnergyEnergyCorrelators(true),
  fDrawEnergyEnergyCorrelatorsForConstantJetPt(false),
  fDrawEnergyEnergyCorrelatorsForConstantTrackPt(false),
  fSaveFigures(false),
  fFigureFormat("pdf"),
  fLogPt(true),
  fColorPalette(kRainBow),
  fStyle2D("colz"),
  fStyle3D("surf1")
{
  
  // Read card from inputfile and collision system from card
  TString collisionSystem = fHistograms->GetSystem();
  
  // Make a string for collision system based on information on the card
  fSystemAndEnergy = Form("%s 5.02 TeV",collisionSystem.Data());
  fCompactSystemAndEnergy = fSystemAndEnergy;
  fCompactSystemAndEnergy.ReplaceAll(" ","");
  fCompactSystemAndEnergy.ReplaceAll(".","v");
  
  // Create a new drawer
  fDrawer = new JDrawer();

  // Do not draw anything by default
  for(int iTrackType = 0; iTrackType < EECHistogramManager::knTrackCategories; iTrackType++){
    fDrawTracks[iTrackType] = false;
  }
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    fDrawEnergyEnergyCorrelators[iEnergyEnergyCorrelator] = false;
  }

  // Setup the centrality, track pT and energy-energy correlator jet and track pT bins to be drawn
  fFirstDrawnCentralityBin = fHistograms->GetFirstCentralityBin();
  fLastDrawnCentralityBin = fHistograms->GetLastCentralityBin();
  fFirstDrawnTrackPtBin = fHistograms->GetFirstTrackPtBin();
  fLastDrawnTrackPtBin = fHistograms->GetLastTrackPtBin();
  fFirstDrawnJetPtBinEEC = fHistograms->GetFirstJetPtBinEEC();
  fLastDrawnJetPtBinEEC = fHistograms->GetLastJetPtBinEEC();
  fFirstDrawnTrackPtBinEEC = fHistograms->GetFirstTrackPtBinEEC();
  fLastDrawnTrackPtBinEEC = fHistograms->GetLastTrackPtBinEEC();
  
}

/*
 * Destructor
 */
EECDrawer::~EECDrawer(){
  delete fDrawer;
}

/*
 * Draw all the selected histograms using JDrawer
 */
void EECDrawer::DrawHistograms(){
  
  // Draw the event information histograms
  DrawEventInformation();

  // Draw the jet histograms
  DrawJetHistograms();
  
  // Draw the track histograms
  DrawTrackHistograms();
  
  // Draw the energy-energy correlation histograms
  DrawEnergyEnergyCorrelationHistograms();
  
}

/*
 * Draw event information histograms
 */
void EECDrawer::DrawEventInformation(){
  
  if(!fDrawEventInformation) return;  // Only draw the event information histograms if they are selected for drawing
  
  // Helper variables for histogram drawing
  TH1D *drawnHistogram;
  TLegend *legend;
  
  // === Vertex z-position ===
  drawnHistogram = fHistograms->GetHistogramVertexZ();
  drawnHistogram->SetMarkerStyle(kFullDiamond);
  fDrawer->DrawHistogram(drawnHistogram,"v_{z} (cm)","#frac{dN}{dv_{z}}  (1/cm)"," ");
  legend = new TLegend(0.65,0.75,0.85,0.9);
  SetupLegend(legend);
  legend->Draw();
  
  // Save the figure to a file
  SaveFigure("vz");
  
  // === Event cuts ===
  drawnHistogram = fHistograms->GetHistogramEvents();
  fDrawer->DrawHistogram(drawnHistogram," ","Number of events", " ");
  legend = new TLegend(0.17,0.22,0.37,0.37);
  SetupLegend(legend);
  legend->Draw();
  
  // Save the figure to a file
  SaveFigure("eventCuts");
  
  // === Track cuts ===
  drawnHistogram = fHistograms->GetHistogramTrackCuts();
  fDrawer->DrawHistogram(drawnHistogram," ","Number of tracks", " ");
  legend = new TLegend(0.65,0.75,0.85,0.9);
  SetupLegend(legend);
  legend->Draw();
  
  // Save the figure to a file
  SaveFigure("trackCuts");
  
  // === Centrality ===
  drawnHistogram = fHistograms->GetHistogramCentrality();
  drawnHistogram->SetMarkerStyle(kFullDiamond);
  fDrawer->DrawHistogram(drawnHistogram,"Centrality percentile","N"," ");
  legend = new TLegend(0.63,0.75,0.83,0.9);
  SetupLegend(legend);
  legend->Draw();
  
  // Save the figure to a file
  SaveFigure("centrality");
  
}

/*
 * Draw jet histograms
 */
void EECDrawer::DrawJetHistograms(){
  
  // If jet histograms are not drawn, there is nothing to do here
  if(!fDrawJets) return;
  
  // Helper variables for histogram drawing
  TH1D *drawnHistogram;
  TH2D *drawnHistogram2D;
  TLegend *legend;
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  char namerX[100];
  char namerY[100];
  
  // Loop over centrality
  for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
    
    centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
    compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
    
    // Select logarithmic drawing for pT
    fDrawer->SetLogY(fLogPt);
    
    // === Jet pT ===
    drawnHistogram = fHistograms->GetHistogramJetPt(iCentrality);
    sprintf(namerX,"%s p_{T}  (GeV)",fHistograms->GetJetAxisName());
    fDrawer->DrawHistogram(drawnHistogram,namerX,"#frac{dN}{dp_{T}}  (1/GeV)"," ");
    legend = new TLegend(0.62,0.75,0.82,0.9);
    SetupLegend(legend,centralityString);
    legend->Draw();
    
    // Save the figure to a file
    sprintf(namerX,"%sPt",fHistograms->GetJetHistogramName());
    SaveFigure(namerX,compactCentralityString);
    
    // Set linear drawing
    fDrawer->SetLogY(false);
    
    // === Jet phi ===
    drawnHistogram = fHistograms->GetHistogramJetPhi(iCentrality);
    sprintf(namerX,"%s #varphi",fHistograms->GetJetAxisName());
    fDrawer->DrawHistogram(drawnHistogram,namerX,"#frac{dN}{d#varphi}"," ");
    legend = new TLegend(0.62,0.75,0.82,0.9);
    SetupLegend(legend,centralityString);
    legend->Draw();
    
    // Save the figure to a file
    sprintf(namerX,"%sPhi",fHistograms->GetJetHistogramName());
    SaveFigure(namerX,compactCentralityString);
    
    // === Jet eta ===
    drawnHistogram = fHistograms->GetHistogramJetEta(iCentrality);
    sprintf(namerX,"%s #eta",fHistograms->GetJetAxisName());
    fDrawer->DrawHistogram(drawnHistogram,namerX,"#frac{dN}{d#eta}"," ");
    legend = new TLegend(0.4,0.20,0.82,0.35);
    SetupLegend(legend,centralityString);
    legend->Draw();
    
    // Save the figure to a file
    sprintf(namerX,"%sEta",fHistograms->GetJetHistogramName());
    SaveFigure(namerX,compactCentralityString);
    
    // Change the right margin better suited for 2D-drawing
    fDrawer->SetRightMargin(0.1);
    
    // === Jet eta vs. phi ===
    drawnHistogram2D = fHistograms->GetHistogramJetEtaPhi(iCentrality);
    sprintf(namerX,"%s #varphi",fHistograms->GetJetAxisName());
    sprintf(namerY,"%s #eta",fHistograms->GetJetAxisName());
    fDrawer->DrawHistogram(drawnHistogram2D,namerX,namerY," ",fStyle2D);
    legend = new TLegend(0.17,0.78,0.37,0.93);
    SetupLegend(legend,centralityString);
    legend->Draw();
    
    // Save the figures to file
    sprintf(namerX,"%sEtaPhi",fHistograms->GetJetHistogramName());
    SaveFigure(namerX,compactCentralityString);
    
    // Change right margin back to 1D-drawing
    fDrawer->SetRightMargin(0.06);
    
  } // Centrality loop
}

/*
 * Draw the track histograms
 */
void EECDrawer::DrawTrackHistograms(){
  
  // Helper variables for histogram drawing
  TH1D *drawnHistogram;
  TH2D *drawnHistogram2D;
  TLegend *legend;
  double legendX1;
  double legendY1;
  double legendX2;
  double legendY2;
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  TString trackPtString;
  TString compactTrackPtString;
  char namerX[100];
  char namerY[100];
  
  // Number of events for normalization
  int numberOfEvents = fHistograms->GetNEvents();  // Normalize with the number of all events for inclusive histograms. TODO: Fix the normalization, currently not reliable
  
  // Loop over track types
  for(int iTrackType = 0; iTrackType < EECHistogramManager::knTrackCategories; iTrackType++){
    if(!fDrawTracks[iTrackType]) continue;  // Only draw selected tracks
    
    // Loop over centrality
    for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
      
      centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
      
      // Select logarithmic drawing for pT
      fDrawer->SetLogY(fLogPt);
      
      // === Track pT ===
      drawnHistogram = fHistograms->GetHistogramTrackPt(iTrackType,iCentrality);
      drawnHistogram->Scale(1.0/numberOfEvents);  // Normalize with the number of events
      sprintf(namerX,"%s p_{T}  (GeV)",fHistograms->GetTrackAxisName(iTrackType));
      fDrawer->DrawHistogram(drawnHistogram,namerX,"#frac{1}{N_{event}} #frac{dN}{dp_{T}}  (1/GeV)"," ");
      legend = new TLegend(0.62,0.75,0.82,0.9);
      SetupLegend(legend,centralityString);
      legend->Draw();
      
      // Save the figure to a file
      sprintf(namerX,"%sPt",fHistograms->GetTrackHistogramName(iTrackType));
      SaveFigure(namerX,compactCentralityString);
      
      // Select linear drawing
      fDrawer->SetLogY(false);
      
      // Loop over track pT bins
      for(int iTrackPt = fFirstDrawnTrackPtBin; iTrackPt <= fHistograms->GetNTrackPtBins(); iTrackPt++){
        
        // Draw the selected track pT bins and the special bin containing integrated distributions
        if(iTrackPt > fLastDrawnTrackPtBin && iTrackPt != fHistograms->GetNTrackPtBins()) continue;
        
        // Set the correct track pT bins
        trackPtString = Form("Track pT: %.1f-%.1f GeV",fHistograms->GetTrackPtBinBorder(iTrackPt),fHistograms->GetTrackPtBinBorder(iTrackPt+1));
        compactTrackPtString = Form("_pT=%.1f-%.1f",fHistograms->GetTrackPtBinBorder(iTrackPt),fHistograms->GetTrackPtBinBorder(iTrackPt+1));
        compactTrackPtString.ReplaceAll(".","v");
        
        // No pT selection for integrated distributions
        if(iTrackPt == fHistograms->GetNTrackPtBins()){
          trackPtString = "";
          compactTrackPtString = "";
        }
        
        // === Track phi ===
        drawnHistogram = fHistograms->GetHistogramTrackPhi(iTrackType,iCentrality,iTrackPt);
        drawnHistogram->Scale(1.0/numberOfEvents);  // Normalize with the number of events
        sprintf(namerX,"%s #varphi",fHistograms->GetTrackAxisName(iTrackType));
        fDrawer->DrawHistogram(drawnHistogram,namerX,"#frac{1}{N_{event}} #frac{dN}{d#varphi}"," ");
        legend = new TLegend(0.17,0.20,0.37,0.35);
        SetupLegend(legend,centralityString,trackPtString);
        legend->Draw();
        
        // Save the figure to a file
        sprintf(namerX,"%sPhi",fHistograms->GetTrackHistogramName(iTrackType));
        SaveFigure(namerX,compactCentralityString,compactTrackPtString);
        
        // === Track eta ===
        drawnHistogram = fHistograms->GetHistogramTrackEta(iTrackType,iCentrality,iTrackPt);
        drawnHistogram->Scale(1.0/numberOfEvents);  // Normalize with the number of events
        
        // Set nice position for the legend
        legendY1 = 0.20; legendY2 = legendY1+0.15;
        if(iTrackPt == fHistograms->GetNTrackPtBins()){
          legendX1 = 0.4; legendX2 = legendX1+0.2;
        } else if (iTrackPt == fHistograms->GetNTrackPtBins() - 1){
          legendX1 = 0.32; legendX2 = legendX1+0.2;
        } else {
          legendX1 = 0.34; legendX2 = legendX1+0.2;
        }
        
        sprintf(namerX,"%s #eta",fHistograms->GetTrackAxisName(iTrackType));
        fDrawer->DrawHistogram(drawnHistogram,namerX,"#frac{1}{N_{event}} #frac{dN}{d#eta}"," ");
        legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
        SetupLegend(legend,centralityString,trackPtString);
        legend->Draw();
        
        // Save the figure to a file
        sprintf(namerX,"%sEta",fHistograms->GetTrackHistogramName(iTrackType));
        SaveFigure(namerX,compactCentralityString,compactTrackPtString);
        
        // Change the right margin better suited for 2D-drawing
        fDrawer->SetRightMargin(0.1);
        
        // === Track eta-phi ===
        drawnHistogram2D = fHistograms->GetHistogramTrackEtaPhi(iTrackType,iCentrality,iTrackPt);
        drawnHistogram2D->Scale(1.0/numberOfEvents);  // Normalize with the number of events
        sprintf(namerX,"%s #varphi",fHistograms->GetTrackAxisName(iTrackType));
        sprintf(namerY,"%s #eta",fHistograms->GetTrackAxisName(iTrackType));
        fDrawer->DrawHistogram(drawnHistogram2D,namerX,namerY," ",fStyle2D);
        legend = new TLegend(0.17,0.78,0.37,0.93);
        SetupLegend(legend,centralityString,trackPtString);
        legend->Draw();
        
        // Save the figure to a file
        sprintf(namerX,"%sEtaPhi",fHistograms->GetTrackHistogramName(iTrackType));
        SaveFigure(namerX,compactCentralityString,compactTrackPtString);
        
        // Change right margin back to 1D-drawing
        fDrawer->SetRightMargin(0.06);
        
      } // Track pT loop
    } // Centrality loop
  } // Track type loop
}

/*
 * Draw energy-energy correlator histograms
 */
void EECDrawer::DrawEnergyEnergyCorrelationHistograms(){
  
  // Helper variables for histogram drawing
  TH1D *drawnHistogram;
  TLegend *legend;
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  TString trackPtString;
  TString compactTrackPtString;
  TString jetPtString;
  TString compactJetPtString;
  char namerY[100];
  
  int color[10] = {kBlack, kBlue, kRed, kGreen+3, kCyan, kMagenta, kOrange-1, kAzure-1, kOrange-1, kGray};
    
  // Loop over energy-energy correlator types
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
        
    // Only draw the selected histograms
    if(!fDrawEnergyEnergyCorrelators[iEnergyEnergyCorrelator]) continue;
        
    // Loop over centrality
    for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
      
      centralityString = Form("Cent: %.0f-%.0f%%",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f",fHistograms->GetCentralityBinBorder(iCentrality),fHistograms->GetCentralityBinBorder(iCentrality+1));
      
      // Draw individual energy-energy correlator histograms
      if(fDrawIndividualEnergyEnergyCorrelators){
        
        // Loop over track pT bins
        for(int iTrackPt = fFirstDrawnTrackPtBinEEC; iTrackPt <= fLastDrawnTrackPtBinEEC; iTrackPt++){
          
          trackPtString = Form("%.1f < track p_{T}",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString = Form("_T>%.1f",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString.ReplaceAll(".","v");
          
          // Loop over jet pT bins
          for(int iJetPt = fFirstDrawnJetPtBinEEC; iJetPt <= fLastDrawnJetPtBinEEC; iJetPt++){
            
            jetPtString = Form("%.0f < jet p_{T} < %.0f", fHistograms->GetJetPtBinBorderEEC(iJetPt), fHistograms->GetJetPtBinBorderEEC(iJetPt+1));
            compactJetPtString = Form("_J=%.0f-%.0f",fHistograms->GetJetPtBinBorderEEC(iJetPt), fHistograms->GetJetPtBinBorderEEC(iJetPt+1));
            
            // === Energy-energy correlator ===
            drawnHistogram = fHistograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt);
            sprintf(namerY,"%s", fHistograms->GetEnergyEnergyCorrelatorAxisName(iEnergyEnergyCorrelator));
            fDrawer->DrawHistogram(drawnHistogram,"#Deltar",namerY," ");
            legend = new TLegend(0.62,0.75,0.82,0.9);
            SetupLegend(legend,centralityString,trackPtString,jetPtString);
            legend->Draw();
            
            // Save the figure to a file
            SaveFigure(fHistograms->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator), compactCentralityString, compactTrackPtString, compactJetPtString);
            
          } // Jet pT loop
        } // Track pT loop
      } // Draw individual energy-energy correlator histograms
      
      // Draw all track pT cuts in the same plot with the jet pT integrated histogram TODO: Automatic scaling for y-axes, add jet pT bins
      if(fDrawEnergyEnergyCorrelatorsForConstantJetPt){
        
        // Only one legend for the plot
        legend = new TLegend(0.62,0.35,0.82,0.9);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*) 0, fSystemAndEnergy.Data(), "");
        legend->AddEntry((TObject*) 0, centralityString.Data(),"");
        legend->AddEntry((TObject*) 0,"Jet p_{T} > 120 GeV","");
        
        // Loop over track pT bins
        for(int iTrackPt = fFirstDrawnTrackPtBinEEC; iTrackPt <= fLastDrawnTrackPtBinEEC; iTrackPt++){
          
          trackPtString = Form("%.1f < track p_{T}",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
          
          drawnHistogram = fHistograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, fHistograms->GetNJetPtBinsEEC(), iTrackPt);
          drawnHistogram->Scale(1/drawnHistogram->Integral("width"));
          drawnHistogram->SetLineColor(color[iTrackPt]);
          
          if(iTrackPt == fFirstDrawnTrackPtBinEEC){
            sprintf(namerY,"%s", fHistograms->GetEnergyEnergyCorrelatorAxisName(iEnergyEnergyCorrelator));
            fDrawer->DrawHistogram(drawnHistogram,"#Deltar",namerY," ");
          } else {
            drawnHistogram->Draw("same");
          }
          
          legend->AddEntry(drawnHistogram, trackPtString.Data(), "l");
          
        } // Track pT loop
        
        legend->Draw();
        
        // Save the figure to a file
        sprintf(namerY,"%sConstantJetPt120",fHistograms->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator));
        SaveFigure(namerY, compactCentralityString);
        
      } // Draw all track pT cuts in the same plot with the jet pT integrated histogram
      
      // Draw all jet pT selections for a single track cut TODO: Automatic scaling for y-axes
      if(fDrawEnergyEnergyCorrelatorsForConstantTrackPt){
        
        // Loop over track pT bins
        for(int iTrackPt = fFirstDrawnTrackPtBinEEC; iTrackPt <= fLastDrawnTrackPtBinEEC; iTrackPt++){
          
          trackPtString = Form("%.1f < track p_{T}",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString = Form("%.1f",fHistograms->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString.ReplaceAll(".","v");
          
          // Only one legend for the plot
          legend = new TLegend(0.62,0.35,0.82,0.9);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          legend->AddEntry((TObject*) 0, fSystemAndEnergy.Data(), "");
          legend->AddEntry((TObject*) 0, centralityString.Data(),"");
          legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
          
          // Loop over jet pT bins
          for(int iJetPt = fFirstDrawnJetPtBinEEC; iJetPt <= fLastDrawnJetPtBinEEC; iJetPt++){
            
            jetPtString = Form("%.0f < jet p_{T} < %.0f", fHistograms->GetJetPtBinBorderEEC(iJetPt), fHistograms->GetJetPtBinBorderEEC(iJetPt+1));
            
            drawnHistogram = fHistograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt);
            drawnHistogram->Scale(1/drawnHistogram->Integral("width"));
            drawnHistogram->SetLineColor(color[iJetPt]);
            
            if(iJetPt == fFirstDrawnJetPtBinEEC){
              sprintf(namerY,"%s", fHistograms->GetEnergyEnergyCorrelatorAxisName(iEnergyEnergyCorrelator));
              fDrawer->DrawHistogram(drawnHistogram,"#Deltar",namerY," ");
            } else {
              drawnHistogram->Draw("same");
            }
            
            legend->AddEntry(drawnHistogram, jetPtString.Data(), "l");
            
          } // Jet pT loop
          
          legend->Draw();
          
          // Save the figure to a file
          sprintf(namerY,"%sConstantTrackPt%s",fHistograms->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator), compactTrackPtString.Data());
          SaveFigure(namerY, compactCentralityString);
          
        } // Track pT loop
      } // Draw all jet pT selections for a single track cut
      
    } // Centrality loop
  } // Energy-energy correlator type loop
}


/*
 * Common legend style setup for figures
 *
 *  TLegend *legend = Pointer to legend that needs setup
 *  TString centralityString = Collision centrality
 *  TString trackString = Track pT information
 *  TString extraString = Additional line to be put into the legend
 */
void EECDrawer::SetupLegend(TLegend *legend, TString centralityString, TString trackString, TString extraString){
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62); // Size: 0.05
  legend->AddEntry((TObject*) 0, fSystemAndEnergy.Data(), "");
  if(fSystemAndEnergy.Contains("PbPb")) legend->AddEntry((TObject*) 0,centralityString.Data(),"");
  if(trackString != "") legend->AddEntry((TObject*) 0,trackString.Data(),"");
  if(extraString != "") legend->AddEntry((TObject*) 0,extraString.Data(),"");
}

/*
 * Save the figure in current canvas to a file
 *
 *  TString figureName = Name for the saved figures
 *  TString centralityString = Information about collision centrality
 *  TString trackPtString = Information about track pT
 *  TString correlationTypeString = Information about correlation type (same/mixed event)
 *  TString deltaPhiString = Information about deltaPhi binning
 */
void EECDrawer::SaveFigure(TString figureName, TString centralityString, TString trackPtString, TString correlationTypeString, TString deltaPhiString){
  
  // Only save the figures if flag is set
  if(!fSaveFigures) return;
  
  // Write the figure to a file
  TString figName = Form("figures/%s_%s",figureName.Data(),fCompactSystemAndEnergy.Data());
  if(fCompactSystemAndEnergy.Contains("PbPb")) figName.Append(centralityString);
  figName.Append(trackPtString);
  figName.Append(correlationTypeString);
  figName.Append(deltaPhiString);
  figName.Append(fFigureSaveNameAppend);
  gPad->GetCanvas()->SaveAs(Form("%s.%s",figName.Data(),fFigureFormat));
  
}

// Setter for drawing event information
void EECDrawer::SetDrawEventInformation(const bool drawOrNot){
  fDrawEventInformation = drawOrNot;
}

// Setter for drawing all jet histograms
void EECDrawer::SetDrawJetHistograms(const bool drawOrNot){
  fDrawJets = drawOrNot;
}

// Setter for drawing tracks
void EECDrawer::SetDrawTracks(const bool drawOrNot){
  fDrawTracks[EECHistogramManager::kTrack] = drawOrNot;
}

// Setter for drawing uncorrected tracks
void EECDrawer::SetDrawTracksUncorrected(const bool drawOrNot){
  fDrawTracks[EECHistogramManager::kUncorrectedTrack] = drawOrNot;
}

// Setter for drawing track histograms
void EECDrawer::SetDrawAllTracks(const bool drawTracks, const bool drawUncorrected){
  SetDrawTracks(drawTracks);
  SetDrawTracksUncorrected(drawUncorrected);
}

// Setter for drawing energy-energy correlator
void EECDrawer::SetDrawEnergyEnergyCorrelor(const bool drawOrNot){
  fDrawEnergyEnergyCorrelators[EECHistogramManager::kEnergyEnergyCorrelator] = drawOrNot;
}

// Setter for drawing jet pT weighted energy-energy correlator
void EECDrawer::SetDrawEnergyEnergyCorrelorJetPt(const bool drawOrNot){
  fDrawEnergyEnergyCorrelators[EECHistogramManager::kEnergyEnergyCorrelatorJetPt] = drawOrNot;
}

// Setter for drawing uncorrected energy-energy correlator
void EECDrawer::SetDrawEnergyEnergyCorrelorUncorrected(const bool drawOrNot){
  fDrawEnergyEnergyCorrelators[EECHistogramManager::kEnergyEnergyCorrelatorUncorrected] = drawOrNot;
}

// Setter for drawing uncorrected jet pT weighted energy-energy correlator
void EECDrawer::SetDrawEnergyEnergyCorrelorJetPtUncorrected(const bool drawOrNot){
  fDrawEnergyEnergyCorrelators[EECHistogramManager::kEnergyEnergyCorrelatorJetPtUncorrected] = drawOrNot;
}

// Setter for drawing all energy-energy correlators
void EECDrawer::SetDrawAllEnergyEnergyCorrelors(const bool drawRegular, const bool drawJetPt, const bool drawUncorrected, const bool drawJetPtUncorrected){
  SetDrawEnergyEnergyCorrelor(drawRegular);
  SetDrawEnergyEnergyCorrelorJetPt(drawJetPt);
  SetDrawEnergyEnergyCorrelorUncorrected(drawUncorrected);
  SetDrawEnergyEnergyCorrelorJetPtUncorrected(drawJetPtUncorrected);
}

// Setter for drawing the individual energy-energy correlator histograms
void EECDrawer::SetDrawSingleEnergyEnergyCorrelators(const bool drawOrNot){
  fDrawIndividualEnergyEnergyCorrelators = drawOrNot;
}

// Setter for drawing all track pT cuts to the same figure for constant jet pT selection
void EECDrawer::SetDrawEnergyEnergyCorrelatorsForConstantJetPt(const bool drawOrNot){
  fDrawEnergyEnergyCorrelatorsForConstantJetPt = drawOrNot;
}

// Setter for drawing all jet pT selections to the same figure for constant track pT cut
void EECDrawer::SetDrawEnergyEnergyCorrelatorsForConstantTrackPt(const bool drawOrNot){
  fDrawEnergyEnergyCorrelatorsForConstantTrackPt = drawOrNot;
}

// Setter for saving the figures to a file
void EECDrawer::SetSaveFigures(const bool saveOrNot, const char *format, const TString suffix){
  fSaveFigures = saveOrNot;
  fFigureFormat = format;
  fFigureSaveNameAppend = suffix;
}

// Setter for logarithmic pT axis
void EECDrawer::SetLogPt(const bool isLog){
  fLogPt = isLog;
}

// Setter for color palette
void EECDrawer::SetColorPalette(const int color){
  fColorPalette = color;
  gStyle->SetPalette(color);
}

// Setter for 2D drawing style
void EECDrawer::SetDrawingStyle2D(const char* style){
  fStyle2D = style;
}

// Setter for 3D drawing style
void EECDrawer::SetDrawingStyle3D(const char* style){
  fStyle3D = style;
}

// Setter for 2D drawing style
void EECDrawer::SetDrawingStyles(const int color, const char* style2D, const char* style3D){
  SetColorPalette(color);
  SetDrawingStyle2D(style2D);
  SetDrawingStyle3D(style3D);
}
