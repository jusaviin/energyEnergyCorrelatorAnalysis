#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"
#include "AlgorithmLibrary.h"

/*
 * Macro studying the processed energy-energy correlator distributions
 */
void eecSignalExaminer(){

  // Open the PbPb input file
  TString inputFileName = "data/eecAnalysis_akFlowJets_removeBadAcceptance_wtaAxis_processed_2022-10-25.root";
  // data/eecAnalysis_akFlowJets_removeBadAcceptance_wtaAxis_processed_2022-10-25.root
  // data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_4pC_newSRC_cutBadPhiAndComb_wtaAxis_jetTrigger_preprocessed_2022-11-29.root
  // data/PbPbMC2018_GenGen_eecAnalysis_akFlowJet_MnD_wtaAxis_noTrigger_preprocessed_2022-10-21.root
  TFile* inputFile = TFile::Open(inputFileName);
  
  // Check that the files exist
  if(inputFile == NULL){
    cout << "Error! The file " << inputFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }

  EECCard* card = new EECCard(inputFile);
  
  // Open the pp input file
  TString ppFileName = "data/ppData_pfJets_wtaAxis_subtractReflectedCone_processed_2022-12-16.root";
  // data/ppData_pfJets_wtaAxis_noBackgroundSubtraction_processed_2022-12-16.root
  // data/ppData_pfJets_wtaAxis_subtractReflectedCone_processed_2022-12-16.root
  TFile *ppFile = TFile::Open(ppFileName);
  
  // Check that the files exist
  if(ppFile == NULL){
    cout << "Error! The file " << ppFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }

  EECCard* ppCard = new EECCard(ppFile);
  
  // ====================================================
  //               Binning configuration
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = card->GetNCentralityBins();
  const int nJetPtBinsEEC = card->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card->GetNTrackPtBinsEEC();
  
  // Default binning ranges for reference
  // centrality = {0,10,30,50,90}
  // jet pT in energy-energy correlators = {120,140,160,180,200,300,500,5030}
  // track pT in energy-energy correlators = {0.7,1,1.5,2,2.5,3,3.5,4,300}
  
  // Select which histograms are fitted
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = nCentralityBins;
  
  int firstDrawnJetPtBinEEC = 0;
  int lastDrawnJetPtBinEEC = nJetPtBinsEEC-1; // Note: Jets integrated over all pT ranges are in nJetPtBinsEEC bin
  
  int firstDrawnTrackPtBinEEC = 5;
  int lastDrawnTrackPtBinEEC = 5;
  
  // Select the types of energy-energy correlators are studied
  bool studyEnergyEnergyCorrelator[EECHistogramManager::knEnergyEnergyCorrelatorTypes];
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator] = true;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPt] = false;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorUncorrected] = false;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPtUncorrected] = false;
  
  // Select the processing levels the are studied
  bool studyEnergyEnergyCorrelatorProcessingLevel[EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels];
  studyEnergyEnergyCorrelatorProcessingLevel[EECHistogramManager::kEnergyEnergyCorrelatorNormalized] = false;
  studyEnergyEnergyCorrelatorProcessingLevel[EECHistogramManager::kEnergyEnergyCorrelatorBackground] = false;
  studyEnergyEnergyCorrelatorProcessingLevel[EECHistogramManager::kEnergyEnergyCorrelatorSignal] = true;
  
  // Find index of one energy-energy correlator that is drawn
  int studiedEnergyEnergyCorrelatorType = -1;
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    if(studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]){
      studiedEnergyEnergyCorrelatorType = iEnergyEnergyCorrelator;
      break;
    }
  }
  
  // Find index of one energy-energy correlator processing level that is drawn
  int studiedEnergyEnergyCorrelatorProcessLevel = -1;
  for(int iProcessLevel = 0; iProcessLevel < EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels; iProcessLevel++){
    if(studyEnergyEnergyCorrelatorProcessingLevel[iProcessLevel]){
      studiedEnergyEnergyCorrelatorProcessLevel = iProcessLevel;
      break;
    }
  }
  
  // Types of comparisons drawn
  const bool drawCentralityComparison = false; // Centrality comparison for constant jet and track pT selection
  const bool drawTrackPtComparison = false;    // Track pT comparison for constant centrality and jet pT selection
  const bool drawJetPtComparison = true;      // Jet pT comparison for constant centrality and track pT selection
  
  // Logarithmic axes
  const bool logDeltaR = true;
  const bool logEEC = true;
  
  // Scaling options
  const bool applyAdHocScale = false;  // Instead of scaling signals to one, match the yield at the peak for different centrality bins
  double adHocScale[nJetPtBinsEEC+1][nTrackPtBinsEEC];
  double adHocLowBinBorder[] = {0.019, 0.017, 0.015, 0.012, 0.01, 0.006, 0.005, 0.017};
  double adHocLowHighBorder[] = {0.026, 0.024, 0.021, 0.019, 0.015, 0.01, 0.009, 0.024};
  
  // Drawing range for x-axis
  std::pair<double,double> drawingRange = std::make_pair(0.001, 0.4); // Cut the histograms at 0.4, since background subtraction does not work beyond that
  
  // Axis zooming
  std::pair<double,double> ratioZoomCentrality = std::make_pair(0.7, 1.3);
  std::pair<double,double> ratioZoomTrackPt = std::make_pair(0.4, 1.6);
  std::pair<double,double> ratioZoomJetPt = std::make_pair(0, 2);
  
  // Figure saving
  const bool saveFigures = false;  // Save figures
  const char* saveComment = "AdHocScale";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures
  
  // Create and setup a new histogram managers for PbPb collisions to project and handle the histograms
  EECHistogramManager *histograms;
  histograms = new EECHistogramManager(inputFile,card);
    
  // Choose the energy-energy correlator types to load
  histograms->SetLoadEnergyEnergyCorrelators(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator]);
  histograms->SetLoadEnergyEnergyCorrelatorsJetPt(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPt]);
  histograms->SetLoadEnergyEnergyCorrelatorsUncorrected(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorUncorrected]);
  histograms->SetLoadEnergyEnergyCorrelatorsJetPtUncorrected(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPtUncorrected]);
    
  // Choose the bin ranges
  histograms->SetCentralityBinRange(0,nCentralityBins-1);
  histograms->SetJetPtBinRangeEEC(0,nJetPtBinsEEC);
  histograms->SetTrackPtBinRangeEEC(0,nTrackPtBinsEEC-1);
    
  // Load the histograms from the file
  histograms->LoadProcessedHistograms();
  
  // Do the same for pp collisions
  EECHistogramManager *ppHistograms;
  ppHistograms = new EECHistogramManager(ppFile,ppCard);
    
  // Choose the energy-energy correlator types to load
  ppHistograms->SetLoadEnergyEnergyCorrelators(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator]);
  ppHistograms->SetLoadEnergyEnergyCorrelatorsJetPt(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPt]);
  ppHistograms->SetLoadEnergyEnergyCorrelatorsUncorrected(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorUncorrected]);
  ppHistograms->SetLoadEnergyEnergyCorrelatorsJetPtUncorrected(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPtUncorrected]);
    
  // Choose the bin ranges
  ppHistograms->SetCentralityBinRange(0,0);
  ppHistograms->SetJetPtBinRangeEEC(0,nJetPtBinsEEC);
  ppHistograms->SetTrackPtBinRangeEEC(0,nTrackPtBinsEEC-1);
    
  // Load the histograms from the file
  ppHistograms->LoadProcessedHistograms();
  
  // Energy-energy correlator histograms
  TH1D* hEnergyEnergyCorrelator[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins+1][nJetPtBinsEEC+1][nTrackPtBinsEEC][EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels];
  TH1D* hEnergyEnergyCorrelatorCentralityRatio[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins+1][nJetPtBinsEEC+1][nTrackPtBinsEEC][EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels];
  TH1D* hEnergyEnergyCorrelatorTrackPtRatio[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins+1][nJetPtBinsEEC+1][nTrackPtBinsEEC][EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels];
  TH1D* hEnergyEnergyCorrelatorJetPtRatio[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins+1][nJetPtBinsEEC+1][nTrackPtBinsEEC][EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels];
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    for(int iCentrality = 0; iCentrality < nCentralityBins+1; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        for(int iJetPt = 0; iJetPt < nJetPtBinsEEC+1; iJetPt++){
          for(int iProcessLevel = 0; iProcessLevel < EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels; iProcessLevel++){
            hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel] = NULL;
            hEnergyEnergyCorrelatorCentralityRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel] = NULL;
            hEnergyEnergyCorrelatorTrackPtRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel] = NULL;
            hEnergyEnergyCorrelatorJetPtRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel] = NULL;
          } // Processing level loop
        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
  } // Energy-energy correlator type loop
  
  double epsilon = 0.0000001;
  int lowSignalBin, highSignalBin;
  int centralityBin;
  double peakIntegral;
  EECHistogramManager *histogramLoader;
  
  // Get the histograms from the histogram manager and normalize the signal histograms to one
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only load the selected energy-energy correlator types
    if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
    
    for(int iCentrality = lastDrawnCentralityBin; iCentrality >= firstDrawnCentralityBin; iCentrality--){
      
      // Select the histogram manager based on the centrality bin
      if(iCentrality == nCentralityBins){
        histogramLoader = ppHistograms;
        centralityBin = 0;
      } else {
        histogramLoader = histograms;
        centralityBin = iCentrality;
      }
      
      for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
        for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
          for(int iProcessLevel = 0; iProcessLevel < EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels; iProcessLevel++){
            
            // Only load the selected processing levels
            if(!studyEnergyEnergyCorrelatorProcessingLevel[iProcessLevel]) continue;
            
            // Read the processed energy-energy correlator histograms
            hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel] = histogramLoader->GetHistogramEnergyEnergyCorrelatorProcessed(iEnergyEnergyCorrelator, centralityBin, iJetPt, iTrackPt, iProcessLevel);
            
            // Normalize the signal distributions to one over the drawn range
            if(iProcessLevel == EECHistogramManager::kEnergyEnergyCorrelatorSignal){
              
              // Scale the signal distribution to one
              lowSignalBin = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->GetXaxis()->FindBin(drawingRange.first+epsilon);
              highSignalBin = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->GetXaxis()->FindBin(drawingRange.second-epsilon);
              hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->Scale(1/hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->Integral(lowSignalBin, highSignalBin, "width"));
              
              // Instead of normalization to one, match the different centrality bins at
              if(applyAdHocScale){
                lowSignalBin = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->GetXaxis()->FindBin(adHocLowBinBorder[iJetPt]+epsilon);
                highSignalBin = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->GetXaxis()->FindBin(adHocLowHighBorder[iJetPt]-epsilon);
                if(iCentrality == lastDrawnCentralityBin){
                  adHocScale[iJetPt][iTrackPt] = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->Integral(lowSignalBin, highSignalBin, "width");
                }
                peakIntegral = hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->Integral(lowSignalBin, highSignalBin, "width");
                hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->Scale(adHocScale[iJetPt][iTrackPt]/peakIntegral);
              }
            }
            
            // Calculate the ratio with respect to the first centrality bin
            hEnergyEnergyCorrelatorCentralityRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel] = (TH1D*) hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->Clone(Form("centralityRatio%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iProcessLevel));
            
            hEnergyEnergyCorrelatorCentralityRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->Divide(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][lastDrawnCentralityBin][iJetPt][iTrackPt][iProcessLevel]);
            
            // Calculate the ratio with respect to the first track pT bin
            hEnergyEnergyCorrelatorTrackPtRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel] = (TH1D*) hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->Clone(Form("trackPtRatio%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iProcessLevel));
            
            hEnergyEnergyCorrelatorTrackPtRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->Divide(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][firstDrawnTrackPtBinEEC][iProcessLevel]);
            
            // Calculate the ratio with respect to the first jet pT bin
            hEnergyEnergyCorrelatorJetPtRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel] = (TH1D*) hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->Clone(Form("jetPtRatio%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iProcessLevel));
            
            hEnergyEnergyCorrelatorJetPtRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->Divide(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][firstDrawnJetPtBinEEC][iTrackPt][iProcessLevel]);
            
          } // Processing level loop
        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
  } // Energy-energy correlator type loop
  
  // Find the binning information for DeltaR bins from a histogram
  const int nDeltaRBins = hEnergyEnergyCorrelator[studiedEnergyEnergyCorrelatorType][firstDrawnCentralityBin][firstDrawnJetPtBinEEC][firstDrawnTrackPtBinEEC][studiedEnergyEnergyCorrelatorProcessLevel]->GetNbinsX();
  double deltaRBinBorders[nDeltaRBins+1];
  deltaRBinBorders[0] = hEnergyEnergyCorrelator[studiedEnergyEnergyCorrelatorType][firstDrawnCentralityBin][firstDrawnJetPtBinEEC][firstDrawnTrackPtBinEEC][studiedEnergyEnergyCorrelatorProcessLevel]->GetXaxis()->GetBinLowEdge(1);
  for(int iBin = 1; iBin <= nDeltaRBins; iBin++){
    deltaRBinBorders[iBin] = hEnergyEnergyCorrelator[studiedEnergyEnergyCorrelatorType][firstDrawnCentralityBin][firstDrawnJetPtBinEEC][firstDrawnTrackPtBinEEC][studiedEnergyEnergyCorrelatorProcessLevel]->GetXaxis()->GetBinUpEdge(iBin);
  }
  
  // ==========================================================================
  //        Draw the selected processed energy-energy correlator styles
  // ==========================================================================
  
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  
  if(logDeltaR) drawer->SetLogX(true);
  int color[10] = {kBlack, kRed, kBlue, kGreen+2, kCyan, kMagenta, kOrange-1, kAzure+5, kOrange-2, kGray};

  TLegend *legend;
  TString centralityString, trackPtString, jetPtString;
  TString compactCentralityString, compactTrackPtString, compactJetPtString;
  TString legendString;
  
  std::pair<double,double> histogramYrange;
  
  // Create an algorithm library object to help find y-axis scaling ranges for the histograms
  AlgorithmLibrary *algorithms = new AlgorithmLibrary();
  
  // ============================================== //
  // Energy-energy correlator centrality comparison //
  // ============================================== //
  
  if(drawCentralityComparison){
    
    // Loop over all selected histograms
    for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
      
      // Only draw the selected energy-energy correlator types
      if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
      
      for(int iProcessLevel = 0; iProcessLevel < EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels; iProcessLevel++){
        
        // Only draw the selected processing levels
        if(!studyEnergyEnergyCorrelatorProcessingLevel[iProcessLevel]) continue;
        
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
            
            trackPtString = Form("%.1f < track p_{T}",histograms->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString = Form("_T>%.1f",histograms->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString.ReplaceAll(".","v");
            
            // Create a legend for the figure
            legend = new TLegend(0.44,0.08,0.69,0.44);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->AddEntry((TObject*) 0, Form("%s 5.02 TeV",histograms->GetCard()->GetAlternativeDataType().Data()), "");
            legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
            legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
            
            // Create a new canvas for the plot
            drawer->CreateSplitCanvas();
            
            // Logarithmic EEC axis
            if(logEEC) drawer->SetLogY(true);
            
            for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
              
              // Set the centrality information for legends and figure saving
              centralityString = Form("Cent: %.0f-%.0f%%", histograms->GetCentralityBinBorder(iCentrality), histograms->GetCentralityBinBorder(iCentrality+1));
              compactCentralityString = Form("_C=%.0f-%.0f", histograms->GetCentralityBinBorder(iCentrality), histograms->GetCentralityBinBorder(iCentrality+1));
              
              // The pp distributions are stored in the nCentralityBins bin
              if(iCentrality == nCentralityBins){
                centralityString = "pp";
                compactCentralityString = "_pp";
              }
              
              // For logarithmic drawing, cannot go down to zero in x-axis
              if(logDeltaR){
                hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
              }
              
              // Draw the histograms to the upper canvas
              hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->SetLineColor(color[iCentrality-firstDrawnCentralityBin]);
              if(iCentrality == firstDrawnCentralityBin){
                drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel], "#Deltar", Form("EEC %s", histograms->GetEnergyEnergyCorrelatorProcessSaveName(iProcessLevel)), " ");
              } else {
                hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->Draw("same");
              }
              
              legend->AddEntry(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel], centralityString.Data(), "l");
              
            } // Centrality loop
            
            // Draw the legends to the upper pad
            legend->Draw();
            
            // Linear scale for the ratio
            drawer->SetLogY(false);
            
            if(lastDrawnCentralityBin == nCentralityBins){
              centralityString = "#frac{Centrality}{pp}";
            } else {
              centralityString = Form("#frac{Centrality}{%.0f-%.0f%%}", histograms->GetCentralityBinBorder(lastDrawnCentralityBin), histograms->GetCentralityBinBorder(lastDrawnCentralityBin+1));
            }
            
            for(int iCentrality = firstDrawnCentralityBin; iCentrality < lastDrawnCentralityBin; iCentrality++){
              
              // For logarithmic drawing, cannot go down to zero in x-axis
              if(logDeltaR){
                hEnergyEnergyCorrelatorCentralityRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
              }
              
              hEnergyEnergyCorrelatorCentralityRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->SetLineColor(color[iCentrality-firstDrawnCentralityBin]);
              hEnergyEnergyCorrelatorCentralityRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->GetYaxis()->SetRangeUser(ratioZoomCentrality.first, ratioZoomCentrality.second);
              if(iCentrality == firstDrawnCentralityBin){
                drawer->SetGridY(true);
                drawer->DrawHistogramToLowerPad(hEnergyEnergyCorrelatorCentralityRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel], "#Deltar", centralityString.Data(), " ");
                drawer->SetGridY(false);
              } else {
                hEnergyEnergyCorrelatorCentralityRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->Draw("same");
              }
              
            } // Centrality loop
            
            // Save the figures to a file
            if(saveFigures){
              gPad->GetCanvas()->SaveAs(Form("figures/%s%sCentralityComparison%s%s%s.%s", histograms->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator), histograms->GetEnergyEnergyCorrelatorProcessSaveName(iProcessLevel), saveComment, compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
            }
            
          } // Track pT loop
        } // Jet pT loop
      } // Energy-energy correlator processing level loop
    } // Energy-energy correlator type loop
  }

  // ============================================== //
  //  Energy-energy correlator track pT comparison  //
  // ============================================== //
  
  if(drawTrackPtComparison){
    
    // Loop over all selected histograms
    for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
      
      // Only draw the selected energy-energy correlator types
      if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
      
      for(int iProcessLevel = 0; iProcessLevel < EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels; iProcessLevel++){
        
        // Only draw the selected processing levels
        if(!studyEnergyEnergyCorrelatorProcessingLevel[iProcessLevel]) continue;
        
        // Loop over centrality bins
        for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
          
          // Select the correct system loader based on the centrality bin
          if(iCentrality == nCentralityBins){
            histogramLoader = ppHistograms;
            compactCentralityString = "_pp";
          } else {
            histogramLoader = histograms;
            centralityString = Form("Cent: %.0f-%.0f%%", histograms->GetCentralityBinBorder(iCentrality), histograms->GetCentralityBinBorder(iCentrality+1));
            compactCentralityString = Form("_C=%.0f-%.0f", histograms->GetCentralityBinBorder(iCentrality), histograms->GetCentralityBinBorder(iCentrality+1));
          }
          
          // Loop over jet pT bins
          for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
            
            // Set the jet pT information for legends and figure saving
            if(iJetPt == histograms->GetNJetPtBinsEEC()){
              jetPtString = Form("Jet p_{T} > %.0f", histogramLoader->GetCard()->GetJetPtCut());
              compactJetPtString = "";
            } else {
              jetPtString = Form("%.0f < jet p_{T} < %.0f", histogramLoader->GetJetPtBinBorderEEC(iJetPt), histogramLoader->GetJetPtBinBorderEEC(iJetPt+1));
              compactJetPtString = Form("_J=%.0f-%.0f", histogramLoader->GetJetPtBinBorderEEC(iJetPt), histogramLoader->GetJetPtBinBorderEEC(iJetPt+1));
            }
            
            // Create a legend for the figure
            legend = new TLegend(0.44,0.08,0.69,0.44);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->AddEntry((TObject*) 0, Form("%s 5.02 TeV",histogramLoader->GetCard()->GetAlternativeDataType().Data()), "");
            if(iCentrality < nCentralityBins) legend->AddEntry((TObject*) 0, centralityString.Data(),"");
            legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
            
            // Create a new canvas for the plot
            drawer->CreateSplitCanvas();
            
            // Logarithmic EEC axis
            if(logEEC) drawer->SetLogY(true);
            
            // Loop over track pT bins
            for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
              
              trackPtString = Form("%.1f < track p_{T}",histogramLoader->GetTrackPtBinBorderEEC(iTrackPt));
              compactTrackPtString = Form("_T>%.1f",histogramLoader->GetTrackPtBinBorderEEC(iTrackPt));
              compactTrackPtString.ReplaceAll(".","v");
              
              // For logarithmic drawing, cannot go down to zero in x-axis
              if(logDeltaR){
                hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
              }
              
              // Draw the histograms to the upper canvas
              hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->SetLineColor(color[iTrackPt-firstDrawnTrackPtBinEEC]);
              if(iTrackPt == firstDrawnTrackPtBinEEC){
                drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel], "#Deltar", Form("EEC %s", histogramLoader->GetEnergyEnergyCorrelatorProcessSaveName(iProcessLevel)), " ");
              } else {
                hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->Draw("same");
              }
              
              legend->AddEntry(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel], trackPtString.Data(), "l");
              
            } // Track pT loop
            
            // Draw the legends to the upper pad
            legend->Draw();
            
            // Linear scale for the ratio
            drawer->SetLogY(false);
            
            trackPtString = Form("#frac{Track p_{T} cut}{%.1f < track p_{T}}", histogramLoader->GetTrackPtBinBorderEEC(firstDrawnTrackPtBinEEC));
            
            for(int iTrackPt = firstDrawnTrackPtBinEEC+1; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
              
              // For logarithmic drawing, cannot go down to zero in x-axis
              if(logDeltaR){
                hEnergyEnergyCorrelatorTrackPtRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
              }
              
              hEnergyEnergyCorrelatorTrackPtRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->SetLineColor(color[iTrackPt-firstDrawnTrackPtBinEEC]);
              hEnergyEnergyCorrelatorTrackPtRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->GetYaxis()->SetRangeUser(ratioZoomTrackPt.first, ratioZoomTrackPt.second);
              if(iTrackPt == firstDrawnTrackPtBinEEC+1){
                drawer->SetGridY(true);
                drawer->DrawHistogramToLowerPad(hEnergyEnergyCorrelatorTrackPtRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel], "#Deltar", trackPtString.Data(), " ");
                drawer->SetGridY(false);
              } else {
                hEnergyEnergyCorrelatorTrackPtRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->Draw("same");
              }
              
            } // Track pT loop
            
            // Save the figures to a file
            if(saveFigures){
              gPad->GetCanvas()->SaveAs(Form("figures/%s%sTrackPtComparison%s%s%s.%s", histogramLoader->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator), histogramLoader->GetEnergyEnergyCorrelatorProcessSaveName(iProcessLevel), saveComment, compactCentralityString.Data(), compactJetPtString.Data(), figureFormat));
            }
            
          } // Jet pT loop
        } // Centrality loop
      } // Energy-energy correlator processing level loop
    } // Energy-energy correlator type loop
  }
  
  // ============================================== //
  //   Energy-energy correlator jet pT comparison   //
  // ============================================== //
  
  if(drawJetPtComparison){
    
    // Loop over all selected histograms
    for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
      
      // Only draw the selected energy-energy correlator types
      if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
      
      for(int iProcessLevel = 0; iProcessLevel < EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels; iProcessLevel++){
        
        // Only draw the selected processing levels
        if(!studyEnergyEnergyCorrelatorProcessingLevel[iProcessLevel]) continue;
        
        // Loop over centrality bins
        for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
          
          // Select the correct system loader based on the centrality bin
          if(iCentrality == nCentralityBins){
            histogramLoader = ppHistograms;
            compactCentralityString = "_pp";
          } else {
            histogramLoader = histograms;
            centralityString = Form("Cent: %.0f-%.0f%%", histograms->GetCentralityBinBorder(iCentrality), histograms->GetCentralityBinBorder(iCentrality+1));
            compactCentralityString = Form("_C=%.0f-%.0f", histograms->GetCentralityBinBorder(iCentrality), histograms->GetCentralityBinBorder(iCentrality+1));
          }
          
          // Loop over track pT bins
          for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
            
            trackPtString = Form("%.1f < track p_{T}",histogramLoader->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString = Form("_T>%.1f",histogramLoader->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString.ReplaceAll(".","v");
            
            // Create a legend for the figure
            legend = new TLegend(0.41,0.04,0.66,0.56);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->AddEntry((TObject*) 0, Form("%s 5.02 TeV",histogramLoader->GetCard()->GetAlternativeDataType().Data()), "");
            if(iCentrality < nCentralityBins) legend->AddEntry((TObject*) 0, centralityString.Data(),"");
            legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
            
            // Create a new canvas for the plot
            drawer->CreateSplitCanvas();
            
            // Logarithmic EEC axis
            if(logEEC) drawer->SetLogY(true);
            
            // First, loop over jet pT bins and determine a good drawing range from that
            histogramYrange = std::make_pair(1e10, -1e10);
            for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
              histogramYrange = algorithms->FindHistogramMinMax(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel], histogramYrange, drawingRange);
            }
            
            // Add some empty space to the top and the bottom of the histogram
            histogramYrange.first = histogramYrange.first - histogramYrange.first*0.15;
            histogramYrange.second = histogramYrange.second + histogramYrange.second*0.15;
            
            // Loop over jet pT bins
            for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
              
              // Set the jet pT information for legends and figure saving
              if(iJetPt == histogramLoader->GetNJetPtBinsEEC()){
                jetPtString = Form("Jet p_{T} > %.0f", histogramLoader->GetCard()->GetJetPtCut());
                compactJetPtString = "";
              } else {
                jetPtString = Form("%.0f < jet p_{T} < %.0f", histogramLoader->GetJetPtBinBorderEEC(iJetPt), histogramLoader->GetJetPtBinBorderEEC(iJetPt+1));
                compactJetPtString = Form("_J=%.0f-%.0f", histogramLoader->GetJetPtBinBorderEEC(iJetPt), histogramLoader->GetJetPtBinBorderEEC(iJetPt+1));
              }
              
              // For logarithmic drawing, cannot go down to zero in x-axis
              if(logDeltaR){
                hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
              }
              
              // Set a good y-axis drawing range for the histograms
              hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->GetYaxis()->SetRangeUser(histogramYrange.first, histogramYrange.second);
              
              // Draw the histograms to the upper canvas
              hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->SetLineColor(color[iJetPt-firstDrawnJetPtBinEEC]);
              if(iJetPt == firstDrawnJetPtBinEEC){
                drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel], "#Deltar", Form("EEC %s", histogramLoader->GetEnergyEnergyCorrelatorProcessSaveName(iProcessLevel)), " ");
              } else {
                hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->Draw("same");
              }
              
              legend->AddEntry(hEnergyEnergyCorrelator[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel], jetPtString.Data(), "l");
              
            } // Jet pT loop
            
            // Draw the legends to the upper pad
            legend->Draw();
            
            // Linear scale for the ratio
            drawer->SetLogY(false);
            
            jetPtString = Form("#frac{Jet p_{T} selection}{%.0f < jet p_{T} < %.0f}", histograms->GetJetPtBinBorderEEC(firstDrawnJetPtBinEEC), histograms->GetJetPtBinBorderEEC(firstDrawnJetPtBinEEC+1));
            
            for(int iJetPt = firstDrawnJetPtBinEEC+1; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
              
              // For logarithmic drawing, cannot go down to zero in x-axis
              if(logDeltaR){
                hEnergyEnergyCorrelatorJetPtRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
              }
              
              hEnergyEnergyCorrelatorJetPtRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->SetLineColor(color[iJetPt-firstDrawnJetPtBinEEC]);
              hEnergyEnergyCorrelatorJetPtRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->GetYaxis()->SetRangeUser(ratioZoomJetPt.first, ratioZoomJetPt.second);
              if(iJetPt == firstDrawnJetPtBinEEC+1){
                drawer->SetGridY(true);
                drawer->DrawHistogramToLowerPad(hEnergyEnergyCorrelatorJetPtRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel], "#Deltar", jetPtString.Data(), " ");
                drawer->SetGridY(false);
              } else {
                hEnergyEnergyCorrelatorJetPtRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iProcessLevel]->Draw("same");
              }
              
            } // Jet pT loop
            
            // Save the figures to a file
            if(saveFigures){
              gPad->GetCanvas()->SaveAs(Form("figures/%s%sJetPtComparison%s%s%s.%s", histogramLoader->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator), histogramLoader->GetEnergyEnergyCorrelatorProcessSaveName(iProcessLevel), saveComment, compactCentralityString.Data(), compactTrackPtString.Data(), figureFormat));
            }
            
          } // Track pT loop
        } // Centrality loop
      } // Energy-energy correlator processing level loop
    } // Energy-energy correlator type loop
  }
  
}
