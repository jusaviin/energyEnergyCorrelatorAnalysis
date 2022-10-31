#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Macro studying particle densities in different situations
 */
void particleDensityFitter(){

  // Open the input file
  TString inputFileName = "data/eecAnalysis_akFlowJets_removeBadAcceptance_wtaAxis_processed_2022-10-25.root";
  TFile* inputFile = TFile::Open(inputFileName);
  
  // Check that the files exist
  if(inputFile == NULL){
    cout << "Error! The file " << inputFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }

  EECCard* card = new EECCard(inputFile);
  
  // ====================================================
  //    Binning configuration for the background study
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = card->GetNCentralityBins();
  const int nJetPtBinsEEC = card->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card->GetNTrackPtBinsEEC();
  
  // Default binning ranges for reference
  // centrality = {0,10,30,50,90}
  // track pT =  {0.7,1,2,3,4,8,12,300}
  // jet pT in energy-energy correlators = {120,140,160,180,200,300,500,5030}
  // track pT in energy-energy correlators = {0.7,1,2,3,4,6,300}
  
  // Select which histograms are fitted
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = 0;
  
  int firstDrawnJetPtBinEEC = 0;
  int lastDrawnJetPtBinEEC = 0; // Note: Jets integrated over all pT ranges are in nJetPtBinsEEC bin
  
  int firstDrawnTrackPtBinEEC = 0;
  int lastDrawnTrackPtBinEEC = nTrackPtBinsEEC-1;
  
  // Select the types of energy-energy correlators are studied
  bool drawParticleDensityType[EECHistogramManager::knParticleDensityAroundJetAxisTypes];
  drawParticleDensityType[EECHistogramManager::kParticleDensityAroundJetAxis] = false;
  drawParticleDensityType[EECHistogramManager::kParticlePtDensityAroundJetAxis] = false;
  drawParticleDensityType[EECHistogramManager::kParticleDensityAroundJetAxisPtBinned] = true;
  drawParticleDensityType[EECHistogramManager::kParticlePtDensityAroundJetAxisPtBinned] = false;
  
  // Logarithmic axes
  const bool logY = true;
  
  // Axis zooming
  std::pair<double,double> ratioZoom = std::make_pair(0.9, 1.1);
  
  // Figure saving
  const bool saveFigures = false;  // Save figures
  const char* saveComment = "";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures
  
  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager *histograms;
  histograms = new EECHistogramManager(inputFile,card);
    
  // Choose the particle density types to load
  histograms->SetLoadParticleDensityAroundJets(true);
  histograms->SetLoadParticlePtDensityAroundJets(true);
  histograms->SetLoadParticleDensityAroundJetsPtBinned(true);
  histograms->SetLoadParticlePtDensityAroundJetsPtBinned(true);
    
  // Choose the bin ranges
  histograms->SetCentralityBinRange(0,nCentralityBins-1);
  histograms->SetJetPtBinRangeEEC(0,nJetPtBinsEEC);
  histograms->SetTrackPtBinRangeEEC(0,nTrackPtBinsEEC-1);
    
  // Load the histograms from the file
  histograms->LoadProcessedHistograms();
  
  // Particle density histograms from the data
  TH1D* hParticleDensity[EECHistogramManager::knParticleDensityAroundJetAxisTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC][EECHistograms::knJetConeTypes];
  TH1D* hParticleDensityToFitRatio[EECHistogramManager::knParticleDensityAroundJetAxisTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC][EECHistograms::knJetConeTypes];
  
  // Initialize the energy-energy correlator histogram array to NULL
  for(int iParticleDensity = 0; iParticleDensity < EECHistogramManager::knParticleDensityAroundJetAxisTypes; iParticleDensity++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        for(int iJetCone = 0; iJetCone < EECHistograms::knJetConeTypes; iJetCone++){
          for(int iJetPt = 0; iJetPt < nJetPtBinsEEC+1; iJetPt++){
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone] = NULL;
            hParticleDensityToFitRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone] = NULL;
          } // Jet pT loop
        } // Jet cone type loop
      } // Track pT loop
    } // Centrality loop
  } // Particle density type loop
  
  double normalizationFactor;
  double lowFitRegion[2] = {0.341,0.05};
  double highFitRegion[2] = {0.449,0.5};
  
  // Get the histograms from the histogram manager and fit them with a constant
  for(int iParticleDensity = 0; iParticleDensity < EECHistogramManager::knParticleDensityAroundJetAxisTypes; iParticleDensity++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        for(int iJetCone = 0; iJetCone < EECHistograms::knJetConeTypes; iJetCone++){
          for(int iJetPt = 0; iJetPt <= nJetPtBinsEEC; iJetPt++){
            
            // Find the jet pT normalization factor
            if(iJetPt == histograms->GetNJetPtBinsEEC()){
              normalizationFactor = histograms->GetJetPtIntegral(iCentrality);
            } else {
              normalizationFactor = histograms->GetJetPtIntegral(iCentrality, histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1));
            }
            
            // Read the particle density histograms and normalize integrals to one
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone] = histograms->GetHistogramParticleDensityAroundJetAxis(iCentrality, iJetPt, iTrackPt, iJetCone, iParticleDensity);
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->Scale(1/normalizationFactor);
            hParticleDensityToFitRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone] = (TH1D*) hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->Clone(Form("distributionToFitRatio%d%d%d%d%d", iParticleDensity, iCentrality, iTrackPt, iJetCone, iJetPt));
            
            // Fit a constant to the tail of the distribution, use region 0.3 < DeltaR < 0.5 for now. This can be optimized later if the approach works
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->Fit("pol1","Q","",lowFitRegion[iJetCone],highFitRegion[iJetCone]);
            
            // Divide the distribution with a with a fit the get the ratio
            hParticleDensityToFitRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->Divide(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetFunction("pol1"));
            
          } // Jet pT loop
        } // Jet cone type loop
      } // Track pT loop
    } // Centrality loop
  } // Particle density type loop
  
  // =======================================================
  //    Collect the slope parameters from the linear fits
  // =======================================================
  for(int iParticleDensity = 0; iParticleDensity < EECHistogramManager::knParticleDensityAroundJetAxisTypes; iParticleDensity++){
    for(int iJetCone = 0; iJetCone < EECHistograms::knJetConeTypes; iJetCone++){
      cout << " double " << histograms->GetParticleDensityAroundJetAxisSaveName(iParticleDensity) << histograms->GetJetConeTypeSaveName(iJetCone) << "[" << nCentralityBins << "][" << nJetPtBinsEEC << "][" << nTrackPtBinsEEC << "] = {" << endl;
      for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
        cout << " // " << Form("Centrality %.0f-%.0f%%", histograms->GetCentralityBinBorder(iCentrality), histograms->GetCentralityBinBorder(iCentrality+1)) << endl;
        cout << "{";
        for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
          cout << "{";
          for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
            if(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetFunction("pol1") == NULL){
              cout << "-1";
            } else {
              cout << hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetFunction("pol1")->Eval(0.34) / hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetFunction("pol1")->Eval(0.4);
            }
            if(iTrackPt < nTrackPtBinsEEC-1){
              cout << ",";
            }
          } // Track pT loop
          cout << "}";
          if(iJetPt < nJetPtBinsEEC-1){
            cout << ", // " << Form("%.0f < jet pT < %.0f GeV", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1)) << endl;
          } else {
            cout << "}";
            if(iCentrality < nCentralityBins-1){
              cout << ", // " << Form("%.0f < jet pT < %.0f GeV", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1)) << endl;
            } else {
              cout << " // " << Form("%.0f < jet pT < %.0f GeV", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1)) << endl;
            }
          }
        } // Jet pT loop
      } // Centrality loop
      cout << "};" << endl;
    } // Jet cone type loop
  } // Particle density type loop
  
  // ==========================================================================
  //    Draw the particle density histograms to check the quality of the fit
  // ==========================================================================
  
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);

  TLegend *legend;
  TString centralityString, trackPtString, jetPtString;
  TString compactCentralityString, compactTrackPtString, compactJetPtString;
  TString legendString;
  
  // Loop over all selected histograms
  for(int iParticleDensity = 0; iParticleDensity < EECHistogramManager::knParticleDensityAroundJetAxisTypes; iParticleDensity++){
    
    // Only read the selected energy-energy correlator types
    if(!drawParticleDensityType[iParticleDensity]) continue;
    
    for(int iJetCone = 0; iJetCone < EECHistograms::knJetConeTypes; iJetCone++){
      
      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        
        // Set the centrality information for legends and figure saving
        centralityString = Form("Cent: %.0f-%.0f%%", histograms->GetCentralityBinBorder(iCentrality), histograms->GetCentralityBinBorder(iCentrality+1));
        compactCentralityString = Form("_C=%.0f-%.0f", histograms->GetCentralityBinBorder(iCentrality), histograms->GetCentralityBinBorder(iCentrality+1));
        
        // Loop over jet pT bins
        for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
          
          // Set the jet pT information for legends and figure saving
          if(iJetPt == histograms->GetNJetPtBinsEEC()){
            jetPtString = Form("Jet p_{T} > %.0f", histograms->GetCard()->GetJetPtCut());
            compactJetPtString = "";
          } else {
            jetPtString = Form("%.0f < jet p_{T} < %.0f", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1));
            compactJetPtString = Form("_J=%.0f-%.0f", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1));
          }
          
          // Loop over track pT bins
          for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
            
            // Track pT binning is different depending on the particle density type
            if(iParticleDensity == EECHistogramManager::kParticleDensityAroundJetAxisPtBinned || iParticleDensity == EECHistogramManager::kParticlePtDensityAroundJetAxisPtBinned){
              trackPtString = Form("%.1f < track p_{T} < %.1f",histograms->GetTrackPtBinBorderEEC(iTrackPt), histograms->GetTrackPtBinBorderEEC(iTrackPt+1));
              compactTrackPtString = Form("_T%.1f-%.1f",histograms->GetTrackPtBinBorderEEC(iTrackPt), histograms->GetTrackPtBinBorderEEC(iTrackPt+1));
              compactTrackPtString.ReplaceAll(".","v");
            } else {
              trackPtString = Form("%.1f < track p_{T}",histograms->GetTrackPtBinBorderEEC(iTrackPt));
              compactTrackPtString = Form("_T>%.1f",histograms->GetTrackPtBinBorderEEC(iTrackPt));
              compactTrackPtString.ReplaceAll(".","v");
            }
            
            // Create a legend for the figure
            legend = new TLegend(0.58,0.48,0.83,0.84);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->AddEntry((TObject*) 0, histograms->GetCard()->GetAlternativeDataType().Data(), "");
            legend->AddEntry((TObject*) 0, centralityString.Data(),"");
            legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
            legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
            if(iJetCone == EECHistograms::kReflectedCone) legend->AddEntry((TObject*) 0,"Reflected cone","");
            
            // Create a new canvas for the plot
            drawer->CreateSplitCanvas();
            
            // Logarithmic EEC axis
            if(logY) drawer->SetLogY(true);
            
            // Set the x-axis range
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetXaxis()->SetRangeUser(0.0, 0.6);
            
            // Draw the histograms to the upper canvas
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetLineColor(kBlack);
            drawer->DrawHistogramToUpperPad(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone], "#Deltar", histograms->GetParticleDensityAroundJetAxisAxisName(iParticleDensity), " ");
            
            // Draw the legends to the upper pad
            legend->Draw();
            
            // Linear scale for the ratio
            drawer->SetLogY(false);
            
            // Set the x-axis range for the ratio histogram
            hParticleDensityToFitRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetXaxis()->SetRangeUser(0.0, 0.6);
            
            hParticleDensityToFitRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetLineColor(kBlack);
            drawer->SetGridY(true);
            hParticleDensityToFitRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
            drawer->DrawHistogramToLowerPad(hParticleDensityToFitRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone], "#Deltar", "#frac{Distribution}{Fit}", " ");
            drawer->SetGridY(false);
            
            
          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Jet cone loop
  } // Particle density type loop


}
