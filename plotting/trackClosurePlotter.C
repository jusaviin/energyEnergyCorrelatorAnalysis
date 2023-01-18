#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "JDrawer.h"

//Define types for histogram arrays
enum enumCollisionSystem {kPp, kPbPb, knCollisionSystems};
enum enumDataType {kData, kMC, knDataTypes};
enum enumMonteCarloType{kRecoReco, kRecoGen, knMonteCarloTypes};

/*
 * Our stylist will give your histogram a fresh new look!
 *
 *  TH1* histogram = Histogram in need of a stylist
 *  int markerStyle = Marker style to be set for the histogram
 *  int markerColor = Marker color to be set for the histogram
 *  int rebin = Rebin to be applied for the histogram
 */
void stylistForHistogram(TH1* histogram, int markerStyle, int markerColor, int rebin){
  
  // Set the merker style and color
  histogram->SetMarkerStyle(markerStyle);
  histogram->SetMarkerColor(markerColor);
  histogram->SetLineColor(markerColor);
  
  // Do the rebin an scaling to retain the normalization
  if(rebin > 1){
    histogram->Rebin(rebin);
    histogram->Scale(1.0/rebin);
  }
}

/*
 * Plotter for track closure histograms
 *
 *  JDrawer *drawer = Drawer doing all the dirty work for drawing the histograms
 *  TH1 *genHistogram = Histogram with generator level tracks
 *  TH1 *recoHistogram = Histogram with corrected reconstructed tracks
 *  TH1 *uncorrectedHistogram = Histogram with uncorrected reconstructed tracks
 *  const char *xTitle = Title given for the x-axis of the histograms
 *  const char *yTitle = Title given for the y-axis of the histograms
 *  bool logScale = Use logarithmic scale for the main histogram
 *  int rebin = Rebin for the histograms
 *  double yZoom = Range for the y-axis in the main histogram
 *  bool bigZoom = Zoom to very close range to show exactly the deviation from one
 *  legendX = Left side x-position of the title
 *  legendY = Bottom side y-position of the title
 *  const char *header = Header for the legend
 *  const char *centralityString = String for the centrality
 *  const char *trackPtString = String for track pT
 */
void plotTrackClosure(JDrawer *drawer, TH1 *genHistogram, TH1 *recoHistogram, TH1 *uncorrectedHistogram, const char *xTitle, const char *yTitle, bool logScale, int rebin, double yZoom, bool bigZoom, double legendX, double legendY, const char *header, const char *centralityString, const char *trackPtString){

  // Helper variable to name the histograms
  const char* namer;
  
  // Style settings for the markers
  int markerStyleGen = 21;     // Full color square marker
  int markerColorGen = 38;     // Grey-blue
  int markerStyleReco = 20;    // Full color round masker
  int markerColorCorr = kRed;  // Red
  int markerColorReco = kBlue; // Blue
 
  stylistForHistogram(genHistogram,markerStyleGen,markerColorGen,rebin);
  genHistogram->SetMarkerSize(1.2);
  
  drawer->SetLogY(logScale);
  drawer->CreateSplitCanvas();
  if(yZoom > 0) genHistogram->GetYaxis()->SetRangeUser(0,yZoom);
  drawer->DrawHistogramToUpperPad(genHistogram,xTitle,yTitle," ");

  double ySpace = 0;
  if(strncmp(centralityString,"",2) != 0) ySpace += 0.025;
  if(strncmp(trackPtString,"",2) != 0) ySpace += 0.025;
  
  stylistForHistogram(recoHistogram,markerStyleReco,markerColorCorr,rebin);
  recoHistogram->Draw("same");
  stylistForHistogram(uncorrectedHistogram,markerStyleReco,markerColorReco,rebin);
  if(!bigZoom)  uncorrectedHistogram->Draw("same");

  TLegend *legend = new TLegend(legendX,legendY-ySpace,legendX+0.3,legendY+0.25+ySpace);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);

  legend->SetHeader(header);
  if(strncmp(centralityString,"",2) != 0) legend->AddEntry((TObject*) 0,centralityString,"");
  if(strncmp(trackPtString,"",2) != 0) legend->AddEntry((TObject*) 0,trackPtString,"");
  legend->AddEntry(genHistogram,"Gen. particles","pl");
  legend->AddEntry(recoHistogram,"Corr. tracks","pl");
  if(!bigZoom) legend->AddEntry(uncorrectedHistogram,"Reco. tracks","pl");
  legend->Draw();

  namer = Form("%sRatio",recoHistogram->GetName());
  TH1D *recoRatio = (TH1D*) recoHistogram->Clone(namer);
  recoRatio->Divide(genHistogram);
  
  namer = Form("%sRatio",uncorrectedHistogram->GetName());
  TH1D *uncorrectedRatio = (TH1D*) uncorrectedHistogram->Clone(namer);
  uncorrectedRatio->Divide(genHistogram);
  
  drawer->SetLogY(false);
  if(bigZoom){
    recoRatio->GetYaxis()->SetRangeUser(0.9,1.1);
  } else {
    recoRatio->GetYaxis()->SetRangeUser(0,2);
  }
  drawer->DrawHistogramToLowerPad(recoRatio,xTitle,"Reco/Gen"," ");
  if(!bigZoom) uncorrectedRatio->Draw("same");

}

/*
 * Plotter for track pT, eta and phi closure histograms
 */
void trackClosurePlotter(){
  
  // ============= //
  // Configuration //
  // ============= //
  
  bool saveFigures = true;          // Save the figures to a file
    
  int ptRebin = 10;                  // Rebin for track pT closure histograms (there are 500 bins)
  int trackAngleRebin = 2;          // Rebin for track eta and phi histograms
  
  // Read the number of bins from histogram manager
  const int nCentralityBins = 4; // Expected number of centrality bins
  const int nTrackPtBins = 7; // Expected number of ttack pT bins
  double centralityBinBorders[] = {0,10,30,50,90};       // Bin borders for centrality
  double trackPtBinBorders[] = {0.7,1,2,3,4,8,12,300};       // Bin borders for track pT
  
  // Select which bins to plot
  int firstCentralityBin = 0;
  int lastCentralityBin = nCentralityBins-1;
  
  int firstTrackPtBin = 0;
  int lastTrackPtBin = nTrackPtBins-1;
  
  // Zoom very close to one in closure plots
  bool bigClosureZoom = false;
  
  // ============= //
  //  Config done  //
  // ============= //
  
  const char *legendNames[knCollisionSystems][knDataTypes+1] = {{"pp","Raw Pythia", "Weighted Pythia"},{"PbPb","Raw P+H", "Weighted P+H"}};
  const char *monteCarloName[knCollisionSystems] = {"Pythia8","Pythia+Hydjet"};
  TString systemString[knCollisionSystems] = {"Pp","PbPb"};
  
  // Open files for the closure tests
  EECHistogramManager *closureManager[knCollisionSystems][knMonteCarloTypes];
  TFile *closureFile[knCollisionSystems][knMonteCarloTypes];
  closureFile[kPp][kRecoReco] = TFile::Open("data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_jetPtWeightForEEC_processed_2023-01-17.root");
  closureFile[kPp][kRecoGen] = TFile::Open("data/ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_processed_2023-01-09.root");
  closureFile[kPbPb][kRecoReco] = TFile::Open("data/PbPbMC2018_RecoReco_eecAnalysis_akFlowJet_MnD_wtaAxis_noTrigger_preprocessed_2022-10-19.root");
  closureFile[kPbPb][kRecoGen] = TFile::Open("data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJet_updatedMultiplicityAndDensity_wtaAxis_noTrigger_preprocessed_2022-10-17.root");
  
  // Load the necessary histograms to histogram managers
  for(int iSystem = 0; iSystem < knCollisionSystems; iSystem++){
    for(int iDataType = 0; iDataType < knMonteCarloTypes; iDataType++){
      closureManager[iSystem][iDataType] = new EECHistogramManager(closureFile[iSystem][iDataType]);
      closureManager[iSystem][iDataType]->SetLoadEventInformation(true);
      closureManager[iSystem][iDataType]->SetLoadAllTracks(true,true);
      closureManager[iSystem][iDataType]->LoadProcessedHistograms();
    }
  }
  
  // Define needed histograms
  TH1D *trackPt[EECHistogramManager::knTrackCategories][knCollisionSystems][knMonteCarloTypes][nCentralityBins];
  TH1D *trackEta[EECHistogramManager::knTrackCategories][knCollisionSystems][knMonteCarloTypes][nCentralityBins][nTrackPtBins+1];
  TH1D *trackPhi[EECHistogramManager::knTrackCategories][knCollisionSystems][knMonteCarloTypes][nCentralityBins][nTrackPtBins+1];
  TH1D *trackPtRatio[EECHistogramManager::knTrackCategories][knCollisionSystems][nCentralityBins];
  TH1D *trackEtaRatio[EECHistogramManager::knTrackCategories][knCollisionSystems][nCentralityBins][nTrackPtBins+1];
  TH1D *trackPhiRatio[EECHistogramManager::knTrackCategories][knCollisionSystems][nCentralityBins][nTrackPtBins+1];
  TH1D *trackMultiplicity[knCollisionSystems][knMonteCarloTypes][nCentralityBins];
  
  // String for finding inclusive histograms
  const char *correctionString[EECHistogramManager::knTrackCategories] = {"","Uncorrected"}; // 0 = Track correction included, 1 = No tracking corrections
  double normalizationFactor;
  
  // ******************************************
  // **    Read the histograms from files    **
  // ******************************************
  
  for(int iSystem = 0; iSystem < knCollisionSystems; iSystem++){
    for(int iMonteCarloType = 0; iMonteCarloType < knMonteCarloTypes; iMonteCarloType++){
      for(int iCentrality = firstCentralityBin; iCentrality <= lastCentralityBin; iCentrality++){
        
        // No centrality binning for pp
        if(iSystem == kPp && iCentrality > 0) continue;
        
        // Read track multiplicity histograms for normalizing the track histogram per event
        trackMultiplicity[iSystem][iMonteCarloType][iCentrality] =  closureManager[iSystem][iMonteCarloType]->GetHistogramMultiplicity(iCentrality);
        normalizationFactor = 1.0/trackMultiplicity[iSystem][iMonteCarloType][iCentrality]->Integral("width");
        
        for(int iCorrection = 0; iCorrection < EECHistogramManager::knTrackCategories; iCorrection++){
          
          // Read track pT histograms
          trackPt[iCorrection][iSystem][iMonteCarloType][iCentrality] = closureManager[iSystem][iMonteCarloType]->GetHistogramTrackPt(iCorrection, iCentrality);
          trackPt[iCorrection][iSystem][iMonteCarloType][iCentrality]->Scale(normalizationFactor);
          
          // Read track eta histograms without pT cut
          trackEta[iCorrection][iSystem][iMonteCarloType][iCentrality][nTrackPtBins] = closureManager[iSystem][iMonteCarloType]->GetHistogramTrackEta(iCorrection, iCentrality, closureManager[iSystem][iMonteCarloType]->GetNTrackPtBins());
          trackEta[iCorrection][iSystem][iMonteCarloType][iCentrality][nTrackPtBins]->Scale(normalizationFactor);
          
          // Read track phi histograms without pT cut
          trackPhi[iCorrection][iSystem][iMonteCarloType][iCentrality][nTrackPtBins] = closureManager[iSystem][iMonteCarloType]->GetHistogramTrackPhi(iCorrection, iCentrality, closureManager[iSystem][iMonteCarloType]->GetNTrackPtBins());
          trackPhi[iCorrection][iSystem][iMonteCarloType][iCentrality][nTrackPtBins]->Scale(normalizationFactor);
          
          for(int iTrackPt = firstTrackPtBin; iTrackPt <= lastTrackPtBin; iTrackPt++){
            
            // Read track eta histograms in pT bins
            trackEta[iCorrection][iSystem][iMonteCarloType][iCentrality][iTrackPt] = closureManager[iSystem][iMonteCarloType]->GetHistogramTrackEta(iCorrection, iCentrality, iTrackPt);
            trackEta[iCorrection][iSystem][iMonteCarloType][iCentrality][iTrackPt]->Scale(normalizationFactor);
            
            // Read track phi histograms in pT bins
            trackPhi[iCorrection][iSystem][iMonteCarloType][iCentrality][iTrackPt] = closureManager[iSystem][iMonteCarloType]->GetHistogramTrackPhi(iCorrection, iCentrality, iTrackPt);
            trackPhi[iCorrection][iSystem][iMonteCarloType][iCentrality][iTrackPt]->Scale(normalizationFactor);
            
          } // Track pT loop
        } // Corrected - Uncorrected loop
      } // Centrality loop
    } // Loop over Monte Carlo types
  } // Collision system loop
  
  
  // Drawing class for drawing your favorite root style
  JDrawer *drawer = new JDrawer();
  drawer->SetCanvasSize(700,700);
  TLegend *legend;
  
  // ***************************************
  // **    Drawing track closure plots    **
  // ***************************************
  
  TString header;
  TString trackPtString;
  TString trackPtName;
  TString centralityString;
  TString centralityName;
  TString figureName;
  TString zoomerName;
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetTitleOffsetY(1.5);
  drawer->SetTitleOffsetX(1.2);
  
  if(bigClosureZoom){
    zoomerName = "bigZoom_";
  } else {
    zoomerName = "";
  }
  
  for(int iSystem = 0; iSystem < knCollisionSystems; iSystem++){
    for(int iCentrality = firstCentralityBin; iCentrality <= lastCentralityBin; iCentrality++){
      
      if(iSystem == kPp && iCentrality > 0) continue;
      
      // Closure plots for track pT
      header = Form("%s, inclusive", monteCarloName[iSystem]);
      if(iSystem == kPp){
        centralityString = "";
        centralityName = "";
      } else {
        centralityString = Form("Cent: %.0f-%.0f%%",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
        centralityName = Form("_C=%.0f-%.0f",centralityBinBorders[iCentrality],centralityBinBorders[iCentrality+1]);
      }
      
      trackPtString = "p_{T} inclusive";
      
      /*
       * Parameters for plotter function
       *
       *  JDrawer *drawer = Drawer doing all the dirty work for drawing the histograms
       *  TH1 *genHistogram = Histogram with generator level tracks
       *  TH1 *recoHistogram = Histogram with corrected reconstructed tracks
       *  TH1 *uncorrectedHistogram = Histogram with uncorrected reconstructed tracks
       *  const char *xTitle = Title given for the x-axis of the histograms
       *  const char *yTitle = Title given for the y-axis of the histograms
       *  bool logScale = Use logarithmic scale for the main histogram
       *  int rebin = Rebin for the histograms
       *  double yZoom = Range for the y-axis in the main histogram
       *  bool bigZoom = Zoom to very close range to show exactly the deviation from one
       *  legendX = Left side x-position of the title
       *  legendY = Bottom side y-position of the title
       *  const char *header = Header for the legend
       *  const char *centralityString = String for the centrality
       *  const char *trackPtString = String for track pT
       */
      
      // Plot the closure for track pT
      plotTrackClosure(drawer, trackPt[0][iSystem][kRecoGen][iCentrality], trackPt[0][iSystem][kRecoReco][iCentrality], trackPt[1][iSystem][kRecoReco][iCentrality], "p_{T} (GeV)", "#frac{dN}{dp_{T}} (1/GeV)", true, ptRebin, -1, bigClosureZoom, 0.5, 0.5, header.Data(), centralityString.Data(), "");
      
      // Save the figures for track pT closure
      if(saveFigures){
        figureName = Form("figures/trackPtClosure_%s%s%s.pdf", zoomerName.Data(), systemString[iSystem].Data(), centralityName.Data());
        gPad->GetCanvas()->SaveAs(figureName);
      }
      
      // Plot the closure for pT inclusive track eta
      plotTrackClosure(drawer, trackEta[0][iSystem][kRecoGen][iCentrality][nTrackPtBins], trackEta[0][iSystem][kRecoReco][iCentrality][nTrackPtBins], trackEta[1][iSystem][kRecoReco][iCentrality][nTrackPtBins], "#eta", "#frac{dN}{d#eta}", false, trackAngleRebin, trackEta[0][iSystem][kRecoReco][iCentrality][nTrackPtBins]->GetMaximum()*2.6, bigClosureZoom, 0.5, 0.5, header.Data(), centralityString.Data(), trackPtString.Data());
      
      // Save the figures for track eta closure
      if(saveFigures){
        figureName = Form("figures/trackEtaClosure_%s%s%s.pdf",zoomerName.Data(),systemString[iSystem].Data(),centralityName.Data());
        gPad->GetCanvas()->SaveAs(figureName);
      }
      
      // Plot the closure for pT inclusive track phi
      plotTrackClosure(drawer, trackPhi[0][iSystem][kRecoGen][iCentrality][nTrackPtBins], trackPhi[0][iSystem][kRecoReco][iCentrality][nTrackPtBins], trackPhi[1][iSystem][kRecoReco][iCentrality][nTrackPtBins], "#phi", "#frac{dN}{d#phi}", false, trackAngleRebin, trackPhi[0][iSystem][kRecoReco][iCentrality][nTrackPtBins]->GetMaximum()*2.6, bigClosureZoom, 0.5, 0.5, header.Data(), centralityString.Data(), trackPtString.Data());
      
      // Save the figures for track phi closure
      if(saveFigures){
        figureName = Form("figures/trackPhiClosure_%s%s%s.pdf",zoomerName.Data(),systemString[iSystem].Data(),centralityName.Data());
        gPad->GetCanvas()->SaveAs(figureName);
      }
      
      for(int iTrackPt = firstTrackPtBin; iTrackPt <= lastTrackPtBin; iTrackPt++){
        
        trackPtString = Form("%.1f < p_{T} < %.1f",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
        trackPtName = Form("_pT=%.0f-%.0f",trackPtBinBorders[iTrackPt],trackPtBinBorders[iTrackPt+1]);
        
        // Closure plots for track eta
        plotTrackClosure(drawer, trackEta[0][iSystem][kRecoGen][iCentrality][iTrackPt], trackEta[0][iSystem][kRecoReco][iCentrality][iTrackPt], trackEta[1][iSystem][kRecoReco][iCentrality][iTrackPt], "#eta", "#frac{dN}{d#eta}", false, trackAngleRebin, trackEta[0][iSystem][kRecoReco][iCentrality][iTrackPt]->GetMaximum()*2.6, bigClosureZoom, 0.5, 0.5, header.Data(), centralityString.Data(), trackPtString.Data());
        
        // Save the figures for track eta closure
        if(saveFigures){
          figureName = Form("figures/trackEtaClosure_%s%s%s%s.pdf",zoomerName.Data(),systemString[iSystem].Data(),centralityName.Data(),trackPtName.Data());
          gPad->GetCanvas()->SaveAs(figureName);
        }
        
        // Closure plots for track phi
        plotTrackClosure(drawer, trackPhi[0][iSystem][kRecoGen][iCentrality][iTrackPt], trackPhi[0][iSystem][kRecoReco][iCentrality][iTrackPt], trackPhi[1][iSystem][kRecoReco][iCentrality][iTrackPt], "#phi", "#frac{dN}{d#phi}", false, trackAngleRebin,  trackPhi[0][iSystem][kRecoReco][iCentrality][iTrackPt]->GetMaximum()*2.6, bigClosureZoom, 0.5, 0.5, header.Data(), centralityString.Data(), trackPtString.Data());
        
        // Save the figures for track phi closure
        if(saveFigures){
          figureName = Form("figures/trackPhiClosure_%s%s%s%s.pdf",zoomerName.Data(),systemString[iSystem].Data(),centralityName.Data(),trackPtName.Data());
          gPad->GetCanvas()->SaveAs(figureName);
        }
      } // Track pT loop
      
    } // Centrality loop
  } // Collision system loop
}
