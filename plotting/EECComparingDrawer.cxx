/*
 * Implementation of EECComparingDrawer
 */

// Root includes
#include <THnSparse.h>
#include <TPad.h>
#include <TF1.h>

// Own includes
#include "EECComparingDrawer.h"

/*
 * Constructor
 */
EECComparingDrawer::EECComparingDrawer(EECHistogramManager *fBaseHistograms) :
  fBaseHistograms(fBaseHistograms),
  fnAddedHistograms(0),
  fMainHistogram(0),
  fDrawEventInformation(false),
  fDrawJets(false),
  fTrackHistogramDrawn(false),
  fEnergyEnergyCorrelatorHistogramDrawn(false),
  fSaveFigures(false),
  fFigureFormat("pdf"),
  fFigureComment(""),
  fManualLegend(false),
  fApplyScaling(0),
  fLogPt(true),
  fLogDeltaR(true),
  fLogEEC(true),
  fUseDifferenceInsteadOfRatio(false),
  fRatioZoomMin(0.6),
  fRatioZoomMax(1.4),
  fRatioLabel(""),
  fColorPalette(kRainBow),
  fStyle2D("colz"),
  fStyle3D("surf1"),
  fRebinJetPt(false)
{
  
  // Create a new drawer
  fDrawer = new JDrawer();
  
  // Initialize histograms to NULL
  for(int iRatios = 0; iRatios < knMaxRatios; iRatios++){
    fAddedHistograms[iRatios] = NULL;
    fComparisonHistogram[iRatios] = NULL;
    fRatioHistogram[iRatios] = NULL;
  }
  
  // Do not draw anything by default
  for(int iTrackType = 0; iTrackType < EECHistogramManager::knTrackCategories; iTrackType++){
    fDrawTracks[iTrackType] = false;
  }
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    fDrawEnergyEnergyCorrelators[iEnergyEnergyCorrelator] = false;
  }
  for(int iComment = 0; iComment < knMaxRatios+1; iComment++){
    fLegendComment[iComment] = "";
    fScalingFactors[iComment] = 1;
  }
  
  // Find the bin ranges from the histogram manager
  fFirstDrawnCentralityBin = fBaseHistograms->GetFirstCentralityBin();
  fLastDrawnCentralityBin = fBaseHistograms->GetLastCentralityBin();
  fFirstDrawnTrackPtBin = fBaseHistograms->GetFirstTrackPtBin();
  fLastDrawnTrackPtBin = fBaseHistograms->GetLastTrackPtBin();
  fFirstDrawnJetPtBinEEC = fBaseHistograms->GetFirstJetPtBinEEC();
  fLastDrawnJetPtBinEEC = fBaseHistograms->GetLastJetPtBinEEC();
  fFirstDrawnTrackPtBinEEC = fBaseHistograms->GetFirstTrackPtBinEEC();
  fLastDrawnTrackPtBinEEC = fBaseHistograms->GetLastTrackPtBinEEC();
  
}

/*
 * Destructor
 */
EECComparingDrawer::~EECComparingDrawer(){
  delete fDrawer;
}

/*
 * Add histograms to draw together with base histograms
 *
 *  Arguments:
 *   EECHistogramManager *additionalHistogram = Histogram manager containing the set of histograms for this dataset
 */
void EECComparingDrawer::AddHistogramToDraw(EECHistogramManager *additionalHistogram){
  if(fnAddedHistograms == knMaxRatios){
    cout << "Already at maximum amount of histograms (" << knMaxRatios << "), cannot add more!" << endl;
    return;
  }
  
  fAddedHistograms[fnAddedHistograms++] = additionalHistogram;
}

/*
 * Add comment to legend
 *
 *  Arguments:
 *   TString comment = Comment given to the legend
 */
void EECComparingDrawer::AddLegendComment(TString comment){
  if(fnAddedHistograms == knMaxRatios){
    cout << "Already at maximum amount of histograms (" << knMaxRatios << "), cannot add more!" << endl;
    return;
  }
  
  fLegendComment[fnAddedHistograms] = comment;
}

/*
 * Draw all the selected histograms using JDrawer
 */
void EECComparingDrawer::DrawHistograms(){
  
  // Draw the event information histograms
  DrawEventInformation();
  
  // Draw the jet histograms
  DrawJetHistograms();
  
  // Draw the track histograms
  DrawTrackHistograms();
  
  // Draw the energy-energy correlator histograms
  DrawEnergyEnergyCorrelatorHistograms();

}

/*
 * Draw event information histograms
 */
void EECComparingDrawer::DrawEventInformation(){
  
  if(!fDrawEventInformation) return;
  
  // Legend helper variable
  TLegend *legend;
  
  // Namer helper
  //char namerX[100];
  TString centralityString;
  TString compactCentralityString;
  
  // === Centrality ===
  
  // Scale the histograms to the counts
  FindScalingFactors("centrality",-1,-1,-1);
  
  // Prepare the jet pT histograms and ratio to be drawn
  PrepareRatio("centrality", 1);
  
  // Draw the centrality histograms to the upper pad
  DrawToUpperPad("Centrality percentile", "Counts");
  
  // Add a legend to the plot
  legend = new TLegend(0.6,0.5,0.9,0.85);
  SetupLegend(legend, "All events");
  legend->Draw();
  
  // Draw the ratios to the lower portion of the split canvas
  DrawToLowerPad("Centrality percentile","Ratio",fRatioZoomMin,fRatioZoomMax);
  
  // Save the figure to a file
  SaveFigure("centralityInAllEvents");
  
  // === Weighted centrality ===
  
  // Scale the histograms to the counts
  FindScalingFactors("centralityWeighted",-1,-1,-1);
  
  // Prepare the jet pT histograms and ratio to be drawn
  PrepareRatio("centralityWeighted", 1);
  
  // Draw the centrality histograms to the upper pad
  DrawToUpperPad("Centrality percentile", "Counts");
  
  // Add a legend to the plot
  legend = new TLegend(0.6,0.5,0.9,0.85);
  SetupLegend(legend, "All events weighted");
  legend->Draw();
  
  // Draw the ratios to the lower portion of the split canvas
  DrawToLowerPad("Centrality percentile","Ratio",fRatioZoomMin,fRatioZoomMax);
  
  // Save the figure to a file
  SaveFigure("centralityInAllEventsWeighted");
  
  // === Multiplicity ===
  
  TH1D *testHistogram = fBaseHistograms->GetHistogramMultiplicity(0);
  
  // Loop over centrality
  if(testHistogram != NULL){
    for(int iCentrality = 0; iCentrality < fBaseHistograms->GetNCentralityBins(); iCentrality++){
      
      centralityString = Form("Cent: %.0f-%.0f%%",fBaseHistograms->GetCentralityBinBorder(iCentrality),fBaseHistograms->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f",fBaseHistograms->GetCentralityBinBorder(iCentrality),fBaseHistograms->GetCentralityBinBorder(iCentrality+1));
      
      // ==================================
      // === Multiplicity in all events ===
      // ==================================
      
      // Scale the histograms to the counts
      FindScalingFactors("multiplicity",iCentrality,-1,-1);
      
      // Prepare the jet pT histograms and ratio to be drawn
      PrepareRatio("multiplicity", 1, iCentrality);
      
      // Draw the centrality histograms to the upper pad
      DrawToUpperPad("Multiplicity", "Counts");
      
      // Add a legend to the plot
      legend = new TLegend(0.6,0.5,0.9,0.85);
      SetupLegend(legend, "All events", centralityString);
      legend->Draw();
      
      // Draw the ratios to the lower portion of the split canvas
      DrawToLowerPad("Multiplicity","Ratio",fRatioZoomMin,fRatioZoomMax);
      
      // Save the figure to a file
      SaveFigure("multiplicity", compactCentralityString);
      
    } // Centrality loop
  } // Multiplicity histograms
}

/*
 * Draw jet histograms
 */
void EECComparingDrawer::DrawJetHistograms(){
  
  if(!fDrawJets) return;
  
  // Legend helper variable
  TLegend *legend;
  
  // Helper variables for centrality naming in figures
  TString centralityString;
  TString compactCentralityString;
  char namerX[100];
  char namerY[100];
  
  // Rebin borders for jet pT
  const int nJetPtRebin = 25;
  double jetPtRebinBorders[] = {120,125,130,135,140,145,150,155,160,170,180,190,200,210,220,230,240,260,280,300,325,350,375,400,450,500};
  AlgorithmLibrary *rebinner = new AlgorithmLibrary();
  
  // Loop over centrality
  for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
    
    centralityString = Form("Cent: %.0f-%.0f%%",fBaseHistograms->GetCentralityBinBorder(iCentrality),fBaseHistograms->GetCentralityBinBorder(iCentrality+1));
    compactCentralityString = Form("_C=%.0f-%.0f",fBaseHistograms->GetCentralityBinBorder(iCentrality),fBaseHistograms->GetCentralityBinBorder(iCentrality+1));
    
    if(fApplyScaling) FindScalingFactors("jetPt", iCentrality);
    
    // === Jet pT ===
    
    // Prepare the jet pT histograms and ratio to be drawn
    PrepareRatio("jetPt", 1, iCentrality);
    
    if(fRebinJetPt){
      fMainHistogram = rebinner->RebinAsymmetric(fMainHistogram,nJetPtRebin,jetPtRebinBorders);
      for(int iAdditional = 0; iAdditional < fnAddedHistograms; iAdditional++){
        fComparisonHistogram[iAdditional] = rebinner->RebinAsymmetric(fComparisonHistogram[iAdditional],nJetPtRebin,jetPtRebinBorders);
        fRatioHistogram[iAdditional] = (TH1D*)fComparisonHistogram[iAdditional]->Clone(Form("thisRatio%d",iAdditional));
        fRatioHistogram[iAdditional]->Divide(fMainHistogram);
      }
    }
    
    // Draw the jet pT distributions to the upper panel of a split canvas plot
    sprintf(namerX,"%s p_{T}  (GeV)",fBaseHistograms->GetJetAxisName());
    if(fApplyScaling){
      sprintf(namerY,"#frac{1}{N_{jet}} #frac{dN}{dp_{T}}  (1/GeV)");
    } else {
      sprintf(namerY,"#frac{dN}{dp_{T}}  (1/GeV)");
    }
    DrawToUpperPad(namerX, namerY, false, fLogPt);
    
    // Add a legend to the plot
    legend = new TLegend(0.27,0.05,0.55,0.3);
    SetupLegend(legend,centralityString);
    legend->Draw();
    
    // Draw the ratios to the lower portion of the split canvas
    DrawToLowerPad(namerX,fRatioLabel.Data(),fRatioZoomMin,fRatioZoomMax);
    
    // Save the figure to a file
    sprintf(namerX,"%sPtRatio",fBaseHistograms->GetJetHistogramName());
    SaveFigure(namerX,compactCentralityString);
    
    // === Jet phi ===
    
    // Prepare the jet phi histograms and ratio to be drawn
    PrepareRatio("jetPhi", 1, iCentrality);
    
    // Draw the jet phi distributions to the upper panel of a split canvas plot
    sprintf(namerX,"%s #varphi",fBaseHistograms->GetJetAxisName());
    if(fApplyScaling){
      sprintf(namerY,"#frac{1}{N_{jet}} #frac{dN}{d#varphi}");
    } else {
      sprintf(namerY,"#frac{dN}{d#varphi}");
    }
    DrawToUpperPad(namerX,namerY);
    
    // Add a legend to the plot
    legend = new TLegend(0.62,0.65,0.82,0.9);
    SetupLegend(legend,centralityString);
    legend->Draw();
    
    // Draw the ratios to the lower portion of the split canvas
    DrawToLowerPad(namerX,fRatioLabel.Data(),fRatioZoomMin,fRatioZoomMax);
    
    // Save the figure to a file
    sprintf(namerX,"%sPhiRatio",fBaseHistograms->GetJetHistogramName());
    SaveFigure(namerX,compactCentralityString);
    
    // === Jet eta ===
    
    // Prepare the jet eta histograms and ratio to be drawn
    PrepareRatio("jetEta", 1, iCentrality);
    
    // Draw the jet eta distributions to the upper panel of a split canvas plot
    sprintf(namerX,"%s #eta",fBaseHistograms->GetJetAxisName());
    if(fApplyScaling){
      sprintf(namerY,"#frac{1}{N_{jet}} #frac{dN}{d#eta}");
    } else {
      sprintf(namerY,"#frac{dN}{d#eta}");
    }
    DrawToUpperPad(namerX,namerY);
    
    // Add a legend to the plot
    legend = new TLegend(0.42,0.07,0.62,0.32);
    SetupLegend(legend,centralityString);
    legend->Draw();
    
    // Draw the ratios to the lower portion of the split canvas
    DrawToLowerPad(namerX,fRatioLabel.Data(),fRatioZoomMin,fRatioZoomMax);
    
    // Save the figure to a file
    sprintf(namerX,"%sEtaRatio",fBaseHistograms->GetJetHistogramName());
    SaveFigure(namerX,compactCentralityString);
    
  } // Centrality loop
}

/*
 * Draw the track histograms
 */
void EECComparingDrawer::DrawTrackHistograms(){
  
  if(!fTrackHistogramDrawn) return;
  
  // Legend helper variable
  TLegend *legend;
  double legendX1;
  double legendY1;
  double legendX2;
  double legendY2;
  
  // Helper variables for figures
  TString centralityString;
  TString compactCentralityString;
  TString trackPtString;
  TString compactTrackPtString;
  char namerX[100];
  
  // Loop over track types
  for(int iTrackType = 0; iTrackType < EECHistogramManager::knTrackCategories; iTrackType++){
    if(!fDrawTracks[iTrackType]) continue;  // Only draw selected tracks
    
    // Loop over centrality
    for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
      
      centralityString = Form("Cent: %.0f-%.0f%%",fBaseHistograms->GetCentralityBinBorder(iCentrality),fBaseHistograms->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f",fBaseHistograms->GetCentralityBinBorder(iCentrality),fBaseHistograms->GetCentralityBinBorder(iCentrality+1));
      
      if(fApplyScaling) FindScalingFactors("trackPt", iTrackType, iCentrality);
      
      // Prepare the track pT histograms and ratio to be drawn
      PrepareRatio("trackPt", 4, iTrackType, iCentrality);
      
      // Draw the track pT distributions to the upper panel of a split canvas plot
      sprintf(namerX,"%s p_{T}  (GeV)",fBaseHistograms->GetTrackAxisName(iTrackType));
      DrawToUpperPad(namerX, "#frac{dN}{dp_{T}}  (1/GeV)", false, fLogPt);
      
      // Add a legend to the plot
      legend = new TLegend(0.56,0.66,0.76,0.81);
      SetupLegend(legend,centralityString);
      legend->Draw();
      
      // Draw the ratios to the lower portion of the split canvas
      fDrawer->SetGridY(true);
      DrawToLowerPad(namerX,fRatioLabel.Data(),fRatioZoomMin,fRatioZoomMax);
      fDrawer->SetGridY(false);
      
      // Save the figure to a file
      sprintf(namerX,"%sPtRatio",fBaseHistograms->GetTrackHistogramName(iTrackType));
      SaveFigure(namerX,compactCentralityString);
      
      // Loop over track pT bins
      for(int iTrackPt = fFirstDrawnTrackPtBin; iTrackPt <= fBaseHistograms->GetNTrackPtBins(); iTrackPt++){
        
        // Draw the selected track pT bins and the special bin containing integrated distributions
        if(iTrackPt > fLastDrawnTrackPtBin && iTrackPt != fBaseHistograms->GetNTrackPtBins()) continue;
        
        // Set the correct track pT bins
        trackPtString = Form("Track pT: %.1f-%.1f GeV",fBaseHistograms->GetTrackPtBinBorder(iTrackPt),fBaseHistograms->GetTrackPtBinBorder(iTrackPt+1));
        compactTrackPtString = Form("_pT=%.1f-%.1f",fBaseHistograms->GetTrackPtBinBorder(iTrackPt),fBaseHistograms->GetTrackPtBinBorder(iTrackPt+1));
        compactTrackPtString.ReplaceAll(".","v");
        
        // No pT selection for integrated distributions
        if(iTrackPt == fBaseHistograms->GetNTrackPtBins()){
          trackPtString = "";
          compactTrackPtString = "";
        }
        
        // === Track phi ===
        
        // Prepare the track phi histograms to be drawn
        PrepareRatio("trackPhi", 1, iTrackType, iCentrality, iTrackPt);
        
        // Draw the track phi distributions to the upper panel of a split canvas plot
        sprintf(namerX,"%s #varphi",fBaseHistograms->GetTrackAxisName(iTrackType));
        DrawToUpperPad(namerX, "#frac{dN}{d#varphi}");
        
        // Add a legend to the plot
        legend = new TLegend(0.24,0.11,0.44,0.26);
        SetupLegend(legend,centralityString,trackPtString);
        legend->Draw();
        
        // Draw the ratios to the lower portion of the split canvas
        fDrawer->SetGridY(true);
        DrawToLowerPad(namerX,fRatioLabel.Data(),fRatioZoomMin,fRatioZoomMax);
        fDrawer->SetGridY(false);
        
        // Save the figure to a file
        sprintf(namerX,"%sPhiRatio",fBaseHistograms->GetTrackHistogramName(iTrackType));
        SaveFigure(namerX,compactCentralityString,compactTrackPtString);
        
        // === Track eta ===
        
        // Set nice position for the legend
        legendY1 = 0.20; legendY2 = legendY1+0.15;
        if(iTrackPt == fBaseHistograms->GetNTrackPtBins()){
          legendX1 = 0.4; legendX2 = legendX1+0.2;
        } else if (iTrackPt == fBaseHistograms->GetNTrackPtBins() - 1){
          legendX1 = 0.32; legendX2 = legendX1+0.2;
        } else {
          legendX1 = 0.34; legendX2 = legendX1+0.2;
        }
        
        // Prepare the track eta histograms to be drawn
        PrepareRatio("trackEta", 1, iTrackType, iCentrality, iTrackPt);
        
        // Draw the track eta distributions to the upper panel of a split canvas plot
        sprintf(namerX,"%s #eta",fBaseHistograms->GetTrackAxisName(iTrackType));
        DrawToUpperPad(namerX, "#frac{dN}{d#eta}");
        
        // Add a legend to the plot
        legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
        SetupLegend(legend,centralityString,trackPtString);
        legend->Draw();
        
        // Draw the ratios to the lower portion of the split canvas
        fDrawer->SetGridY(true);
        DrawToLowerPad(namerX,fRatioLabel.Data(),fRatioZoomMin,fRatioZoomMax);
        fDrawer->SetGridY(false);
        
        // Save the figure to a file
        sprintf(namerX,"%sEtaRatio",fBaseHistograms->GetTrackHistogramName(iTrackType));
        SaveFigure(namerX,compactCentralityString,compactTrackPtString);
        
      } // Track pT loop
    } // Centrality loop
  } // Track type loop
}

// Draw energy-energy correlator histograms
void EECComparingDrawer::DrawEnergyEnergyCorrelatorHistograms(){
 
  if(!fEnergyEnergyCorrelatorHistogramDrawn) return;
  
  // Legend helper variable
  TLegend *legend;
  double legendX1 = 0.37; // Default x1 location for the legend
  double legendY1 = 0.1;  // Default y1 location for the legend
  double legendX2 = 0.57; // Default x2 location for the legend
  double legendY2 = 0.34 + 0.06*fnAddedHistograms; // Default y2 location for the legend
  
  // If we are not doing log-log drawing, move the legend to top right corner
  if(!fLogDeltaR || !fLogEEC){
    legendX1 = 0.52;
    legendY1 = 0.58 - 0.06*fnAddedHistograms;
    legendX2 = 0.72;
    legendY2 = 0.82;
  }
  
  // Helper variables for figures
  TString centralityString;
  TString compactCentralityString;
  TString trackPtString;
  TString compactTrackPtString;
  TString jetPtString;
  TString compactJetPtString;
  char namerY[100];
  
  // Change drawing ranges in case of logarithmic x-axis
  double minRangeX = 0;
  double maxRangeX = 0.8;
  if(fLogDeltaR) minRangeX = 0.001;
  
  // Loop over energy-energy correlator types
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only draw the selected histograms
    if(!fDrawEnergyEnergyCorrelators[iEnergyEnergyCorrelator]) continue;
    
    // Loop over centrality
    for(int iCentrality = fFirstDrawnCentralityBin; iCentrality <= fLastDrawnCentralityBin; iCentrality++){
      
      centralityString = Form("Cent: %.0f-%.0f%%",fBaseHistograms->GetCentralityBinBorder(iCentrality),fBaseHistograms->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f",fBaseHistograms->GetCentralityBinBorder(iCentrality),fBaseHistograms->GetCentralityBinBorder(iCentrality+1));
      
      // Loop over track pT bins
      for(int iTrackPt = fFirstDrawnTrackPtBinEEC; iTrackPt <= fLastDrawnTrackPtBinEEC; iTrackPt++){
        
        trackPtString = Form("%.1f < track p_{T}",fBaseHistograms->GetTrackPtBinBorderEEC(iTrackPt));
        compactTrackPtString = Form("_T>%.1f",fBaseHistograms->GetTrackPtBinBorderEEC(iTrackPt));
        compactTrackPtString.ReplaceAll(".","v");
        
        // Loop over jet pT bins
        for(int iJetPt = fFirstDrawnJetPtBinEEC; iJetPt <= fLastDrawnJetPtBinEEC; iJetPt++){
          
          jetPtString = Form("%.0f < jet p_{T} < %.0f", fBaseHistograms->GetJetPtBinBorderEEC(iJetPt), fBaseHistograms->GetJetPtBinBorderEEC(iJetPt+1));
          compactJetPtString = Form("_J=%.0f-%.0f",fBaseHistograms->GetJetPtBinBorderEEC(iJetPt), fBaseHistograms->GetJetPtBinBorderEEC(iJetPt+1));
          
          // === Energy-energy correlator ===
          
          // Determine scaling factors
          if(fApplyScaling) FindScalingFactors("energyEnergyCorrelator", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, 0, EECHistograms::knSubeventCombinations);
          
          // Prepare the energy-energy correlator histograms to be drawn TODO: Implementation for different pairing types
          PrepareRatio("energyEnergyCorrelator", 1, iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, 0, EECHistograms::knSubeventCombinations, minRangeX, maxRangeX);
          
          // Draw the track phi distributions to the upper panel of a split canvas plot
          sprintf(namerY,"%s", fBaseHistograms->GetEnergyEnergyCorrelatorAxisName(iEnergyEnergyCorrelator));
          DrawToUpperPad("#Deltar", namerY, fLogDeltaR, fLogEEC);
          
          // Add a legend to the plot
          legend = new TLegend(legendX1, legendY1, legendX2, legendY2);
          SetupLegend(legend,centralityString,trackPtString,jetPtString);
          legend->Draw();
          
          // Draw the ratios to the lower portion of the split canvas
          fDrawer->SetGridY(true);
          DrawToLowerPad("#Deltar",fRatioLabel.Data(),fRatioZoomMin,fRatioZoomMax,fLogDeltaR);
          fDrawer->SetGridY(false);
          
          // Save the figure to a file
          sprintf(namerY,"%sComparison", fBaseHistograms->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator));
          SaveFigure(namerY, compactCentralityString, compactTrackPtString, compactJetPtString);
          
        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
    
  } // Energy-energy correlator loop
}

/*
 * Prepare histograms for ratio plots
 *
 *  Arguments:
 *   TString name = Name for the histograms to be filled in arrays
 *   int rebin = Rebinning the histograms before taking the ratio
 *   int bin1 = First bin index for the loaded histograms
 *   int bin2 = Second bin index for the loaded histograms
 *   int bin3 = Third bin index for the loaded histograms
 *   int bin4 = Fourth bin index for the loaded histograms
 *   int bin5 = Fifth bin index for the loaded histograms
 *   int bin6 = Sixth bin index for the loaded histograms
 *   double minRange = Minimum range shown for the histograms
 *   double maxRange = Maximum range shown for the histograms
 */
void EECComparingDrawer::PrepareRatio(TString name, int rebin, int bin1, int bin2, int bin3, int bin4, int bin5, int bin6, double minRange, double maxRange){
  
  // Helper variables
  char namer[100];

  // Read the histograms, do rebinning or normalization if specified, and take the ratio
  fMainHistogram = (TH1D*)fBaseHistograms->GetOneDimensionalHistogram(name,bin1,bin2,bin3,bin4,bin5,bin6)->Clone();
  
  // Possibility to rebin the histograms
  if(rebin > 1) {
    fMainHistogram->Rebin(rebin);
    fMainHistogram->Scale(1.0/rebin);
  }
  
  // If range for the x-axis is specified, apply it
  if(maxRange > minRange) fMainHistogram->GetXaxis()->SetRangeUser(minRange,maxRange);
  
  // Apply normalization for the histograms
  if(fApplyScaling == 1) fMainHistogram->Scale(1.0/fScalingFactors[0]);
  if(fApplyScaling == 2){
    fMainHistogram->Fit("pol0","Q0");
    fMainHistogram->Scale(1/fMainHistogram->GetFunction("pol0")->GetParameter(0));
  }

  // Do the configuration also for the other histograms
  for(int iAdditional = 0; iAdditional < fnAddedHistograms; iAdditional++){
    fComparisonHistogram[iAdditional] = (TH1D*)fAddedHistograms[iAdditional]->GetOneDimensionalHistogram(name,bin1,bin2,bin3,bin4,bin5,bin6)->Clone();
    
    // Possibility to rebin the histograms
    if(rebin > 1) {
      fComparisonHistogram[iAdditional]->Rebin(rebin);
      fComparisonHistogram[iAdditional]->Scale(1.0/rebin);
    }
    
    // If range for the x-axis is specified, apply it
    if(maxRange > minRange) fComparisonHistogram[iAdditional]->GetXaxis()->SetRangeUser(minRange,maxRange);
    
    // Apply normalization for the histograms
    if(fApplyScaling == 1) fComparisonHistogram[iAdditional]->Scale(1.0/fScalingFactors[1+iAdditional]);
    if(fApplyScaling == 2){
      fComparisonHistogram[iAdditional]->Fit("pol0","Q0");
      fComparisonHistogram[iAdditional]->Scale(1/fComparisonHistogram[iAdditional]->GetFunction("pol0")->GetParameter(0));
    }

    sprintf(namer,"%sRatio%d",fMainHistogram->GetName(),iAdditional);
    

    fRatioHistogram[iAdditional] = (TH1D*)fMainHistogram->Clone(namer);
    if(fUseDifferenceInsteadOfRatio){
      fRatioHistogram[iAdditional]->Add(fComparisonHistogram[iAdditional],-1);
    } else {
      fRatioHistogram[iAdditional]->Divide(fComparisonHistogram[iAdditional]);
    }
    
  } // Loop over additional histograms
}

/*
 * Integrate over the drawn distributions, such that the distributions can be scaled to one.
 *
 *  const char* histogramName = Name of the histogram that is integrated
 *  const int bin1 = First bin index for the integrated histogram
 *  const int bin2 = Second bin index for the integrated histogram
 *  const int bin3 = Third bin index for the integrated histogram
 *  const int bin4 = Fourth bin index for the integrated histogram
 *  const int bin5 = Fifth bin index for the integrated histogram
 *  const int bin6 = Sixth bin index for the integrated histogram
 */
void EECComparingDrawer::FindScalingFactors(const char*  histogramName, const int bin1, const int bin2, const int bin3, const int bin4, const int bin5, const int bin6){

  // Helper variable for reading the normalization scales
  TH1D *scaleReader;
  scaleReader = (TH1D*)fBaseHistograms->GetOneDimensionalHistogram(histogramName, bin1, bin2, bin3, bin4, bin5, bin6)->Clone();
  fScalingFactors[0] = scaleReader->Integral("width");
  for(int iAdditional = 0; iAdditional < fnAddedHistograms; iAdditional++){
    scaleReader = (TH1D*)fAddedHistograms[iAdditional]->GetOneDimensionalHistogram(histogramName, bin1, bin2, bin3, bin4, bin5, bin6)->Clone();
    fScalingFactors[1+iAdditional] = scaleReader->Integral("width");
  }
  
}

/*
 * Draw to upper pad the histograms that are most recently prepared for drawing
 *
 *  Arguments:
 *   const char* xTitle = Title given to the x-axis
 *   const char* yTitle = Title given to the y-axis
 *   const bool logX = True: logarithmic x-axis, false = linear x-axis
 *   const bool logY = True: logarithmic y-axis, false = linear y-axis
 */
void EECComparingDrawer::DrawToUpperPad(const char* xTitle, const char* yTitle, const bool logX, const bool logY){

  // Define some nice colors for histograms
  fMainHistogram->SetLineColor(kBlack);

  // Create a split canvas and draw the histograms to the upped part of the canvas
  fDrawer->SetDefaultAppearanceSplitCanvas();
  fDrawer->CreateSplitCanvas();
  fDrawer->SetLogX(logX);
  fDrawer->SetLogY(logY);
  fDrawer->SetTitleOffsetY(1.5); // Y-axis has too much offset by default
  fDrawer->SetTitleOffsetX(1.1); // X-axis has too much offset by default
  fDrawer->DrawHistogramToUpperPad(fMainHistogram,xTitle,yTitle," ");
  
  for(int iAdditional = 0; iAdditional < fnAddedHistograms; iAdditional++){
    fComparisonHistogram[iAdditional]->SetLineColor(fColors[iAdditional]);
    fComparisonHistogram[iAdditional]->Draw("same");
  }
  
  // Reset back to linear for ratio
  fDrawer->SetLogX(false);
  fDrawer->SetLogY(false);
}

/*
 * Draw to the lower pad the ratio histograms that are most recently prepared
 *
 *  Arguments:
 *   const char* xTitle = Title given to the x-axis
 *   const char* yTitle = Title given to the y-axis
 *   const double zoomMin = Minimum value of the y-axis in the ratio plot
 *   const double zoomMax = Maximum value of the y-axis in the ratio plot
 *   const bool logX = Logarithmic x-axis
 */
void EECComparingDrawer::DrawToLowerPad(const char* xTitle, const char* yTitle, const double zoomMin, const double zoomMax, const bool logX){
  
  // Possibility to do logarithmic x-axis
  fDrawer->SetLogX(logX);
  
  // Draw theratio histograms
  if(fnAddedHistograms > 0){
    fRatioHistogram[0]->SetLineColor(fColors[0]);
    fRatioHistogram[0]->GetYaxis()->SetRangeUser(zoomMin,zoomMax);
    fDrawer->DrawHistogramToLowerPad(fRatioHistogram[0],xTitle,yTitle, " ");
  }
  for(int iAdditional = 1; iAdditional < fnAddedHistograms; iAdditional++){
    fRatioHistogram[iAdditional]->SetLineColor(fColors[iAdditional]);
    fRatioHistogram[iAdditional]->Draw("same");
  }
  
  // Reset x-axis to linear
  fDrawer->SetLogX(false);
}

/*
 * Zoom the y-axis scale to the specified region of the distribution
 *
 *  const double maxZoomValue = Maximum value in the distribution used to calculate the zoom
 *  const int nZoomBins = Number of bins around the maxZoomValue used to calculate the scale
 *  const double scaleFactor = How many times the variance is shown around the average value of the bin contents
 *  const bool bothSides = True: Use bins above -maxZoomValue together with bins below maxZoomValue, False = Use only bins below maxZoomValue
 *  const bool asymmetricZoom = Zoom more to the upper side than the lower. Good for zooming to tail below peak.
 */
void EECComparingDrawer::ZoomToRegion(const double maxZoomValue, const int nZoomBins, const double scaleFactor, const bool bothSides, const bool asymmetricZoom){
  
  // Create array for average values in the main histogram and added histograms
  double averageValues[fnAddedHistograms+1];
  double differences[fnAddedHistograms+1];
  
  // Fill the array with average values from the main histogram and added histograms
  std::tie(averageValues[0],differences[0]) = GetHistogramAverageAndDifferenceInRegion(fMainHistogram,maxZoomValue,nZoomBins,bothSides);
  for(int iAdditional = 0; iAdditional < fnAddedHistograms; iAdditional++){
    std::tie(averageValues[iAdditional+1],differences[iAdditional+1]) = GetHistogramAverageAndDifferenceInRegion(fComparisonHistogram[iAdditional],maxZoomValue,nZoomBins,bothSides);
  }
  
  // Find the maximum and minimum value from the array
  double maxAverage = -1e12;
  double minAverage = 1e12;
  double maxDifference = -1;
  for(int iBin = 0; iBin < fnAddedHistograms+1; iBin++){
    if(averageValues[iBin] < minAverage) minAverage = averageValues[iBin];
    if(averageValues[iBin] > maxAverage) maxAverage = averageValues[iBin];
    if(differences[iBin] > maxDifference) maxDifference = differences[iBin];
  }
  
  // Calculate the maximum and minimum scale for the histogram
  double asymmetryUp = asymmetricZoom ? 1.5 : 1;
  double asymmetryDown = asymmetricZoom ? 0.5 : 1;
  double maxScale = maxAverage + maxDifference*scaleFactor*asymmetryUp;
  double minScale = minAverage - maxDifference*scaleFactor*asymmetryDown;
  if(minScale < 0 || maxScale/minScale > 4) minScale = 0;
  
  // Set the new y-axis ranges for the histogram according to obtained values
  fMainHistogram->GetYaxis()->SetRangeUser(minScale,maxScale);
  for(int iAdditional = 0; iAdditional < fnAddedHistograms; iAdditional++){
    fComparisonHistogram[iAdditional]->GetYaxis()->SetRangeUser(minScale,maxScale);
  }
  
}

/*
 * Get the average and maximum difference in bin values from a histogram in the specific area
 *
 *  const double maxZoomValue = Maximum value in the distribution used to calculate the zoom
 *  const int nZoomBins = Number of bins around the maxZoomValue used to calculate the scale
 *  const double scaleFactor = How much around the average scale in the region is shown as a fraction of the average scale
 *  const bool bothSides = True: Use bins above -maxZoomValue together with bins below maxZoomValue, False = Use only bins below maxZoomValue
 *
 *  return: Average and difference of the bin contents in the given region
 */
std::tuple<double,double> EECComparingDrawer::GetHistogramAverageAndDifferenceInRegion(TH1D *histogram, const double maxZoomValue, const int nZoomBins, const bool bothSides){
  
  // Calculate the sum of bin contents from the positive side
  double binSum = 0;
  int nBinsSummedOver = 0;
  double maxBinValue = -1e12;
  double minBinValue = 1e12;
  double binContent = 0;
  int maxBin = histogram->FindBin(maxZoomValue-1e-4);
  for(int iBin = maxBin-nZoomBins+1; iBin <= maxBin; iBin++){
    binContent = histogram->GetBinContent(iBin);
    binSum += binContent;
    if(binContent < minBinValue) minBinValue = binContent;
    if(binContent > maxBinValue) maxBinValue = binContent;
    nBinsSummedOver++;
  }
  
  // If required, add the bin content also from the negative side
  if(bothSides){
    maxBin = histogram->FindBin(-maxZoomValue+1e-4);
    for(int iBin = maxBin; iBin < maxBin+nZoomBins; iBin++){
      binContent = histogram->GetBinContent(iBin);
      binSum += binContent;
      if(binContent < minBinValue) minBinValue = binContent;
      if(binContent > maxBinValue) maxBinValue = binContent;
      nBinsSummedOver++;
    }
  }
  
  // Calculate the average and max difference of the bin contents
  double binAverage = binSum/nBinsSummedOver;
  double binDifference = maxBinValue-minBinValue;
  
  // Return the average value and the difference in bin values of the considered region of the histogram
  return std::make_tuple(binAverage,binDifference);
  
}

/*
 * Common legend style setup for figures
 *
 *  TLegend *legend = Pointer to legend that needs setup
 *  TString centralityString = Collision centrality
 *  TString trackString = Track pT information
 *  TString asymmetryString = Asymmetry information
 */
void EECComparingDrawer::SetupLegend(TLegend *legend, TString centralityString, TString trackString, TString asymmetryString){
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
  if(fBaseHistograms->GetSystem().Contains("PbPb")) legend->AddEntry((TObject*) 0,centralityString.Data(),"");
  if(trackString != "") legend->AddEntry((TObject*) 0,trackString.Data(),"");
  if(asymmetryString != "") legend->AddEntry((TObject*) 0,asymmetryString.Data(),"");
  legend->AddEntry(fMainHistogram,/*fBaseHistograms->GetSystem() + " " +*/ fLegendComment[0],"l");
  for(int iAdditional = 0; iAdditional < fnAddedHistograms; iAdditional++){
    legend->AddEntry(fComparisonHistogram[iAdditional],/*fAddedHistograms[iAdditional]->GetSystem() + " " +*/ fLegendComment[iAdditional+1],"l");
  }
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
void EECComparingDrawer::SaveFigure(TString figureName, TString centralityString, TString trackPtString, TString correlationTypeString, TString deltaPhiString){
  
  // Only save the figures if flag is set
  if(!fSaveFigures) return;
  
  // Write the figure to a file
  TString figName = Form("figures/%s%s",figureName.Data(),fFigureComment);
  if(fBaseHistograms->GetSystem().Contains("PbPb")) figName.Append(centralityString);
  figName.Append(trackPtString);
  figName.Append(correlationTypeString);
  figName.Append(deltaPhiString);
  gPad->GetCanvas()->SaveAs(Form("%s.%s",figName.Data(),fFigureFormat));
  
}

/*
 * Set a nice drawing style for histograms drawn to big pad
 *
 *  TH1 * histogram = Histogram needing a nice style
 *  double rangeX = Maximum drawing range for x-axis.
 *  const char* xTitle = Title for the x-axis
 */
void EECComparingDrawer::SetHistogramStyle(TH1 *histogram, double rangeX, const char* xTitle){
  histogram->GetXaxis()->SetRangeUser(0,rangeX);
  histogram->GetYaxis()->SetRangeUser(fRatioZoomMin,fRatioZoomMax);
  histogram->SetStats(kFALSE);
  histogram->GetXaxis()->SetTitle(xTitle);
  histogram->SetLabelSize(0.09,"xy");
  histogram->SetTitleSize(0.09,"x");
}

// Setter for drawing event information histograms
void EECComparingDrawer::SetDrawEventInformation(const bool drawOrNot){
  fDrawEventInformation = drawOrNot;
}

// Setter for drawing all jet histograms
void EECComparingDrawer::SetDrawJets(const bool drawOrNot){
  fDrawJets = drawOrNot;
}

// Setter for drawing tracks
void EECComparingDrawer::SetDrawTracks(const bool drawOrNot){
  fDrawTracks[EECHistogramManager::kTrack] = drawOrNot;
  CheckFlagsTrack();
}

// Setter for drawing uncorrected tracks
void EECComparingDrawer::SetDrawTracksUncorrected(const bool drawOrNot){
  fDrawTracks[EECHistogramManager::kUncorrectedTrack] = drawOrNot;
  CheckFlagsTrack();
}

// Setter for drawing track histograms
void EECComparingDrawer::SetDrawAllTracks(const bool drawTracks, const bool drawUncorrected){
  SetDrawTracks(drawTracks);
  SetDrawTracksUncorrected(drawUncorrected);
}

// Setter for drawing energy-energy correlators
void EECComparingDrawer::SetDrawEnergyEnergyCorrelator(const bool drawOrNot){
  fDrawEnergyEnergyCorrelators[EECHistogramManager::kEnergyEnergyCorrelator] = drawOrNot;
  CheckFlagsEnergyEnergyCorrelator();
}

// Setter for drawing energy-energy correlators
void EECComparingDrawer::SetDrawEnergyEnergyCorrelatorJetPt(const bool drawOrNot){
  fDrawEnergyEnergyCorrelators[EECHistogramManager::kEnergyEnergyCorrelatorJetPt] = drawOrNot;
  CheckFlagsEnergyEnergyCorrelator();
}

// Setter for drawing energy-energy correlators
void EECComparingDrawer::SetDrawEnergyEnergyCorrelatorUncorrected(const bool drawOrNot){
  fDrawEnergyEnergyCorrelators[EECHistogramManager::kEnergyEnergyCorrelatorUncorrected] = drawOrNot;
  CheckFlagsEnergyEnergyCorrelator();
}

// Setter for drawing energy-energy correlators
void EECComparingDrawer::SetDrawEnergyEnergyCorrelatorJetPtUncorrected(const bool drawOrNot){
  fDrawEnergyEnergyCorrelators[EECHistogramManager::kEnergyEnergyCorrelatorJetPtUncorrected] = drawOrNot;
  CheckFlagsEnergyEnergyCorrelator(); 
}

// Setter for drawing all energy-energy correlators
void EECComparingDrawer::SetDrawAllEnergyEnergyCorrelators(const bool drawRegular, const bool drawJetPt, const bool drawUncorrected, const bool drawJetPtUncorrected){
  SetDrawEnergyEnergyCorrelator(drawRegular);
  SetDrawEnergyEnergyCorrelatorJetPt(drawJetPt);
  SetDrawEnergyEnergyCorrelatorUncorrected(drawUncorrected);
  SetDrawEnergyEnergyCorrelatorJetPtUncorrected(drawJetPtUncorrected);
}

// Setter for saving the figures to a file
void EECComparingDrawer::SetSaveFigures(const bool saveOrNot, const char *format, const char *comment){
  fSaveFigures = saveOrNot;
  fFigureFormat = format;
  fFigureComment = comment;
}

// Set if we should scale the histograms with their integral before comparing them
void EECComparingDrawer::SetApplyScaling(const int applyScaling){
  fApplyScaling = applyScaling;
}

// Tell if we want to rebin the jet pT histograms
void EECComparingDrawer::SetJetPtRebin(const bool doRebin){
  fRebinJetPt = doRebin;
}

// Setter for logarithmic pT axis
void EECComparingDrawer::SetLogPt(const bool isLog){
  fLogPt = isLog;
}

// Setter for logarithmic deltaR axis in energy-energy correlators
void EECComparingDrawer::SetLogDeltaR(const bool isLog){
  fLogDeltaR = isLog;
}

// Setter for logarithmic EEC axis in energy-energy correlators
void EECComparingDrawer::SetLogEEC(const bool isLog){
  fLogEEC = isLog;
}

// Setter for logarithmix axes
void EECComparingDrawer::SetLogAxes(const bool pt, const bool deltaR, const bool eec){
  SetLogPt(pt);
  SetLogDeltaR(deltaR);
  SetLogEEC(eec);
}

// Setter for manual legend setting
void EECComparingDrawer::SetManualLegend(const bool manualLegend){
  fManualLegend = manualLegend;
}

// Setter for plotting difference instead of ratio to lower pad
void EECComparingDrawer::SetUseDifferenceInRatioPlot(const bool useDifference){
  fUseDifferenceInsteadOfRatio = useDifference;
}

// Setter for minimum value of y-axis in ratio plots
void EECComparingDrawer::SetRatioZoomMin(const double minValue){
  fRatioZoomMin = minValue;
}

// Setter for maximum value of y-axis in ratio plots
void EECComparingDrawer::SetRatioZoomMax(const double maxValue){
  fRatioZoomMax = maxValue;
}

// Setter for y-axis values in ratio plots
void EECComparingDrawer::SetRatioZoom(const double minValue, const double maxValue){
  SetRatioZoomMin(minValue);
  SetRatioZoomMax(maxValue);
}

// Setter for the y-axis label in ratio plots
void EECComparingDrawer::SetRatioLabel(TString label){
  fRatioLabel = label;
}

// Setter for color palette
void EECComparingDrawer::SetColorPalette(const int color){
  fColorPalette = color;
  gStyle->SetPalette(color);
}

// Setter for 2D drawing style
void EECComparingDrawer::SetDrawingStyle2D(const char* style){
  fStyle2D = style;
}

// Setter for 3D drawing style
void EECComparingDrawer::SetDrawingStyle3D(const char* style){
  fStyle3D = style;
}

// Setter for 2D drawing style
void EECComparingDrawer::SetDrawingStyles(const int color, const char* style2D, const char* style3D){
  SetColorPalette(color);
  SetDrawingStyle2D(style2D);
  SetDrawingStyle3D(style3D);
}

// Check if any track histograms are drawn
void EECComparingDrawer::CheckFlagsTrack(){
  fTrackHistogramDrawn = false;
  
  for(int iTrackType = 0; iTrackType < EECHistogramManager::knTrackCategories; iTrackType++){
    if(!fDrawTracks[iTrackType]) continue;
    fTrackHistogramDrawn = true;
  }
}

// Check if any energy-energy correlators are drawn
void EECComparingDrawer::CheckFlagsEnergyEnergyCorrelator(){
  fEnergyEnergyCorrelatorHistogramDrawn = false;
  
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    if(!fDrawEnergyEnergyCorrelators[iEnergyEnergyCorrelator]) continue;
    fEnergyEnergyCorrelatorHistogramDrawn = true;
  }
}

// Check flags to see if we are drawing certain types of histograms
void EECComparingDrawer::CheckFlags(){
  
  CheckFlagsTrack();
  CheckFlagsEnergyEnergyCorrelator();
  
}
