#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Draw comparison plots from
 */

/*
 * Macro for studying the reflected cone background estimation using MC
 */
void studyReflectedConeBackground(){

  // File from which the integrals are calculated
  TString inputFileName = "data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJet_includeReflectedCone_noTrigger_2022-09-22_fiveMissing_preprocessed.root";
  
  // Open the input file
  TFile *inputFile = TFile::Open(inputFileName);
  
  if(inputFile == NULL){
    cout << "Error! The file " << inputFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }
  
  // Load the card from the file and read the collision system
  EECCard *card = new EECCard(inputFile);
  TString collisionSystem = card->GetDataType();
  
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
  
  // Bin range to be integrated
  int firstStudiedCentralityBin = 0;
  int lastStudiedCentralityBin = 0;
  
  int firstStudiedJetPtBinEEC = 0;
  int lastStudiedJetPtBinEEC = nJetPtBinsEEC-1; // Note: Jets integrated over all pT ranges are in nJetPtBinsEEC bin
  
  int firstStudiedTrackPtBinEEC = 5;
  int lastStudiedTrackPtBinEEC = 5;
  
  // Select the types of energy-energy correlators are studied
  bool studyEnergyEnergyCorrelator[EECHistogramManager::knEnergyEnergyCorrelatorTypes];
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator] = true;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPt] = false;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorUncorrected] = false;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPtUncorrected] = false;
    
  // Select the types of plots to look at
  const bool drawBackgroundToReflectedConeRatio = false;
  const bool drawSignalFakeToReflectedConeRatio = true;
  const bool drawSignalToSubtrctedRatio = false;
  
  // Logarithmic axes
  const bool logDeltaR = true;
  const bool logEEC = true;
  
  // Zoom for ratio
  const double ratioZoomLow = 0;
  const double ratioZoomHigh = 2;
  
  // Create and setup a new histogram manager to project and handle the histograms
  EECHistogramManager *histograms = new EECHistogramManager(inputFile,card);
  
  // Choose the energy-energy correlator types to load
  histograms->SetLoadEnergyEnergyCorrelators(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator]);
  histograms->SetLoadEnergyEnergyCorrelatorsJetPt(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPt]);
  histograms->SetLoadEnergyEnergyCorrelatorsUncorrected(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorUncorrected]);
  histograms->SetLoadEnergyEnergyCorrelatorsJetPtUncorrected(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPtUncorrected]);
  
  // Choose the bin ranges
  histograms->SetCentralityBinRange(firstStudiedCentralityBin,lastStudiedCentralityBin);
  histograms->SetJetPtBinRangeEEC(firstStudiedJetPtBinEEC,lastStudiedJetPtBinEEC);
  histograms->SetTrackPtBinRangeEEC(firstStudiedTrackPtBinEEC,lastStudiedTrackPtBinEEC);
  
  // Load the histograms from the file
  histograms->LoadProcessedHistograms();
  
  // Initialize the energy-energy correlator histogram array to NULL
  const int nBackgroundNormalizationBins = 9;
  TH1D* hEnergyEnergyCorrelatorTotal[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC];
  TH1D* hEnergyEnergyCorrelatorSignal[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC];
  TH1D* hEnergyEnergyCorrelatorSignalFake[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC];
  TH1D* hEnergyEnergyCorrelatorBackground[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC];
  TH1D* hEnergyEnergyCorrelatorReflectedCone[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC][nBackgroundNormalizationBins];
  TH1D* hBackgroundToReflectedConeRatio[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC][nBackgroundNormalizationBins];
  TH1D* hSignalFakeToReflectedConeRatio[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC][nBackgroundNormalizationBins];
  TH1D* hBackgroundSubtracted[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC][nBackgroundNormalizationBins];
  TH1D* hSignalToBackgroundSubtractedRatio[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC][nBackgroundNormalizationBins];
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC+1; iJetPt++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
          hEnergyEnergyCorrelatorTotal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = NULL;
          hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = NULL;
          hEnergyEnergyCorrelatorSignalFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = NULL;
          hEnergyEnergyCorrelatorBackground[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = NULL;
          for(int iNormalization = 0; iNormalization < nBackgroundNormalizationBins; iNormalization++){
            hEnergyEnergyCorrelatorReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization] = NULL;
            hBackgroundToReflectedConeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization] = NULL;
            hSignalFakeToReflectedConeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization] = NULL;
            hBackgroundSubtracted[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization] = NULL;
            hSignalToBackgroundSubtractedRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization] = NULL;
          }
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Energy-energy correlator loop

  // Helper variables
  TH1D *helperHistogram;
  double normalizationFactor;
  int lowIntegralBin, highIntegralBin;
  int studiedEnergyEnergyCorrelatorIndex = -1;
  
  // Get the histograms from the histogram manager
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only read the selected energy-energy correlator types
    if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
    studiedEnergyEnergyCorrelatorIndex = iEnergyEnergyCorrelator;
    
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
          
          // Histogram with all pair combinations. Normalize it to one
          hEnergyEnergyCorrelatorTotal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = histograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJet, EECHistogramManager::knSubeventTypes);
          normalizationFactor = hEnergyEnergyCorrelatorTotal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Integral("width");
          hEnergyEnergyCorrelatorTotal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Scale(1/normalizationFactor);
          
          // Histogram with only signal. Normalize the sum of signal and background to total.
          hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = histograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJet, EECHistogramManager::kPythiaPythia);
          hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Scale(1/normalizationFactor);
          
          // Histogram with all background contributions. Normalize the sum of signal and background to total.
          hEnergyEnergyCorrelatorSignalFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = histograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJet, EECHistogramManager::kPythiaHydjet);
          hEnergyEnergyCorrelatorBackground[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt] = (TH1D *) hEnergyEnergyCorrelatorSignalFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Clone(Form("energyEnergyCorrelatorBackground%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt));
          helperHistogram = histograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kSameJet, EECHistogramManager::kHydjetHydjet);
          hEnergyEnergyCorrelatorBackground[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Add(helperHistogram);
          hEnergyEnergyCorrelatorSignalFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Scale(1/normalizationFactor);
          hEnergyEnergyCorrelatorBackground[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Scale(1/normalizationFactor);
          
          // The reflected cone histogram.
          helperHistogram = histograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistograms::kReflectedCone, EECHistogramManager::knSubeventTypes);
          
          // Normalize the reflected cone histogram to different DeltaR regions of the tail of the total distribution
          highIntegralBin = hEnergyEnergyCorrelatorTotal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetNbinsX();
          for(int iNormalization = 0; iNormalization < nBackgroundNormalizationBins; iNormalization++){
            lowIntegralBin = highIntegralBin - (nBackgroundNormalizationBins - iNormalization - 1);
            
            hEnergyEnergyCorrelatorReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization] = (TH1D*) helperHistogram->Clone(Form("reflectedCone%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iNormalization));
            normalizationFactor = hEnergyEnergyCorrelatorTotal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Integral(lowIntegralBin, highIntegralBin, "width") / hEnergyEnergyCorrelatorReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->Integral(lowIntegralBin, highIntegralBin, "width");
            hEnergyEnergyCorrelatorReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->Scale(normalizationFactor);
            
            // Calculate the ratio between reflected cone histograms and background histograms
            hBackgroundToReflectedConeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization] = (TH1D*) hEnergyEnergyCorrelatorReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->Clone(Form("backgroundToReflectedConeRatio%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iNormalization));
            hBackgroundToReflectedConeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->Divide(hEnergyEnergyCorrelatorBackground[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]);
            
            // Calculate the ratio between reflected cone histograms and signal-fake histograms
            hSignalFakeToReflectedConeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization] = (TH1D*) hEnergyEnergyCorrelatorReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->Clone(Form("signalFakeToReflectedConeRatio%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iNormalization));
            hSignalFakeToReflectedConeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->Divide(hEnergyEnergyCorrelatorSignalFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]);
            
            // Subtract the reflected cone background from the all pairs distribution
            hBackgroundSubtracted[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization] = (TH1D*) hEnergyEnergyCorrelatorTotal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Clone(Form("backgroundSubtracted%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iNormalization));
            hBackgroundSubtracted[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->Add(hEnergyEnergyCorrelatorReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization], -1);
            
            // Calculate the signal to background subtracted ratio
            hSignalToBackgroundSubtractedRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization] = (TH1D*) hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Clone(Form("signalToBackgroundSubtractedRatio%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, iNormalization));
            hSignalToBackgroundSubtractedRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->Divide(hBackgroundSubtracted[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]);
          }
            
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Energy-energy correlator loop
  
  // Find the binning information for DeltaR bins from a histogram
  const int nDeltaRBins = hEnergyEnergyCorrelatorTotal[studiedEnergyEnergyCorrelatorIndex][firstStudiedCentralityBin][firstStudiedJetPtBinEEC][firstStudiedTrackPtBinEEC]->GetNbinsX();
  double deltaRBinBorders[nDeltaRBins+1];
  deltaRBinBorders[0] = hEnergyEnergyCorrelatorTotal[studiedEnergyEnergyCorrelatorIndex][firstStudiedCentralityBin][firstStudiedJetPtBinEEC][firstStudiedTrackPtBinEEC]->GetXaxis()->GetBinLowEdge(1);
  for(int iBin = 1; iBin <= nDeltaRBins; iBin++){
    deltaRBinBorders[iBin] = hEnergyEnergyCorrelatorTotal[studiedEnergyEnergyCorrelatorIndex][firstStudiedCentralityBin][firstStudiedJetPtBinEEC][firstStudiedTrackPtBinEEC]->GetXaxis()->GetBinUpEdge(iBin);
  }
  
  
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  
  if(logDeltaR) drawer->SetLogX(true);
  
  TLegend *legend;
  TString centralityString, trackPtString, jetPtString;
  TString compactCentralityString, compactTrackPtString, compactJetPtString;
  int color[9] = {kBlue, kRed, kGreen+2, kCyan, kMagenta, kOrange-1, kAzure+5, kOrange-2, kGray};
  double maximumY, currentY;
  
  // First set of plots: compare background histograms to reflected cone histograms with different normalization region
  if(drawBackgroundToReflectedConeRatio){
    
    // Loop over all selevted histograms
    for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
      
      // Only read the selected energy-energy correlator types
      if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
      
      for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
        
        // Set the centrality information for legends and figure saving
        centralityString = Form("Cent: %.0f-%.0f%%",histograms->GetCentralityBinBorder(iCentrality),histograms->GetCentralityBinBorder(iCentrality+1));
        compactCentralityString = Form("_C=%.0f-%.0f",histograms->GetCentralityBinBorder(iCentrality),histograms->GetCentralityBinBorder(iCentrality+1));
        
        for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
          
          // Set the jet pT information for legends and figure saving
          if(iJetPt == histograms->GetNJetPtBinsEEC()){
            jetPtString = Form("Jet p_{T} > %.0f", histograms->GetCard()->GetJetPtCut());
            compactJetPtString = "";
          } else {
            jetPtString = Form("%.0f < jet p_{T} < %.0f", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1));
            compactJetPtString = Form("_J=%.0f-%.0f", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1));
          }
          
          for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
            
            // Set the track pT information for legends and figure saving
            trackPtString = Form("%.1f < track p_{T}",histograms->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString = Form("%.1f",histograms->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString.ReplaceAll(".","v");
            
            // Create a legend for the figure
            legend = new TLegend(0.19,0.09,0.39,0.94);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->AddEntry((TObject*) 0, histograms->GetCard()->GetAlternativeDataType().Data(), "");
            legend->AddEntry((TObject*) 0, centralityString.Data(),"");
            legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
            legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
            
            // Create a new canvas for the plot
            drawer->CreateSplitCanvas();
            
            // Logarithmic EEC axis
            if(logEEC) drawer->SetLogY(true);
            
            // For logarithmic drawing, cannot go down to zero
            if(logDeltaR){
              hEnergyEnergyCorrelatorBackground[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
            }
            
            // Set good y-ranges for plotting
            hEnergyEnergyCorrelatorBackground[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.0005, 2);
            
            // Draw the background histogram to the upper canvas
            hEnergyEnergyCorrelatorBackground[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetLineColor(kBlack);
            drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelatorBackground[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt], "#Deltar", histograms->GetEnergyEnergyCorrelatorAxisName(iEnergyEnergyCorrelator), " ");
            legend->AddEntry(hEnergyEnergyCorrelatorBackground[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt], "All background pairs", "l");
            
            // Draw different normalization regions for the reflected cone histogram to the same plot
            for(int iNormalization = 0; iNormalization < nBackgroundNormalizationBins; iNormalization++){
              
              // For logarithmic drawing, cannot go down to zero
              if(logDeltaR){
                hEnergyEnergyCorrelatorReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
              }
              
              hEnergyEnergyCorrelatorReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->SetLineColor(color[iNormalization]);
              hEnergyEnergyCorrelatorReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->Draw("same");
              legend->AddEntry(hEnergyEnergyCorrelatorReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization], Form("Reflected, #Deltar > %.3f", deltaRBinBorders[nDeltaRBins-(nBackgroundNormalizationBins-iNormalization)]), "l");
              
            } // Normalization region loop
            
            // Draw the legend to the upper pag
            legend->Draw();
            
            // Linear scale for the ratio
            drawer->SetLogY(false);
            
            // Draw the ratios to the lower pad
            for(int iNormalization = 0; iNormalization < nBackgroundNormalizationBins; iNormalization++){
              
              // For logarithmic drawing, cannot go down to zero
              if(logDeltaR){
                hBackgroundToReflectedConeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
              }
              
              hBackgroundToReflectedConeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->SetLineColor(color[iNormalization]);
              if(iNormalization == 0){
                drawer->SetGridY(true);
                hBackgroundToReflectedConeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->GetYaxis()->SetRangeUser(ratioZoomLow, ratioZoomHigh);
                drawer->DrawHistogramToLowerPad(hBackgroundToReflectedConeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization], "#Deltar", "Reflected / BG", " ");
                drawer->SetGridY(false);
              } else {
                hBackgroundToReflectedConeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->Draw("same");
              }
              
            } // Normalization region loop
            
          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Energy-energy correlator loop
    
  } // Background to reflected cone ratio
  
  // Second set of plots: compare signal-fake histograms to reflected cone histograms with different normalization region
  if(drawSignalFakeToReflectedConeRatio){
    
    // Loop over all selevted histograms
    for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
      
      // Only read the selected energy-energy correlator types
      if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
      
      for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
        
        // Set the centrality information for legends and figure saving
        centralityString = Form("Cent: %.0f-%.0f%%",histograms->GetCentralityBinBorder(iCentrality),histograms->GetCentralityBinBorder(iCentrality+1));
        compactCentralityString = Form("_C=%.0f-%.0f",histograms->GetCentralityBinBorder(iCentrality),histograms->GetCentralityBinBorder(iCentrality+1));
        
        for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
          
          // Set the jet pT information for legends and figure saving
          if(iJetPt == histograms->GetNJetPtBinsEEC()){
            jetPtString = Form("Jet p_{T} > %.0f", histograms->GetCard()->GetJetPtCut());
            compactJetPtString = "";
          } else {
            jetPtString = Form("%.0f < jet p_{T} < %.0f", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1));
            compactJetPtString = Form("_J=%.0f-%.0f", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1));
          }
          
          for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
            
            // Set the track pT information for legends and figure saving
            trackPtString = Form("%.1f < track p_{T}",histograms->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString = Form("%.1f",histograms->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString.ReplaceAll(".","v");
            
            // Create a legend for the figure
            legend = new TLegend(0.19,0.09,0.39,0.94);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->AddEntry((TObject*) 0, histograms->GetCard()->GetAlternativeDataType().Data(), "");
            legend->AddEntry((TObject*) 0, centralityString.Data(),"");
            legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
            legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
            
            // Create a new canvas for the plot
            drawer->CreateSplitCanvas();
            
            // Logarithmic EEC axis
            if(logEEC) drawer->SetLogY(true);
            
            // For logarithmic drawing, cannot go down to zero
            if(logDeltaR){
              hEnergyEnergyCorrelatorSignalFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
            }
            
            // Set good y-ranges for plotting
            hEnergyEnergyCorrelatorSignalFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.0005, 2);
            
            // Draw the background histogram to the upper canvas
            hEnergyEnergyCorrelatorSignalFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetLineColor(kBlack);
            drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelatorSignalFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt], "#Deltar", histograms->GetEnergyEnergyCorrelatorAxisName(iEnergyEnergyCorrelator), " ");
            legend->AddEntry(hEnergyEnergyCorrelatorSignalFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt], "Signal-fake pairs", "l");
            
            // Draw different normalization regions for the reflected cone histogram to the same plot
            for(int iNormalization = 0; iNormalization < nBackgroundNormalizationBins; iNormalization++){
              
              // For logarithmic drawing, cannot go down to zero
              if(logDeltaR){
                hEnergyEnergyCorrelatorReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
              }
              
              hEnergyEnergyCorrelatorReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->SetLineColor(color[iNormalization]);
              hEnergyEnergyCorrelatorReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->Draw("same");
              legend->AddEntry(hEnergyEnergyCorrelatorReflectedCone[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization], Form("Reflected, #Deltar > %.3f", deltaRBinBorders[nDeltaRBins-(nBackgroundNormalizationBins-iNormalization)]), "l");
              
            } // Normalization region loop
            
            // Draw the legend to the upper pag
            legend->Draw();
            
            // Linear scale for the ratio
            drawer->SetLogY(false);
            
            // Draw the ratios to the lower pad
            for(int iNormalization = 0; iNormalization < nBackgroundNormalizationBins; iNormalization++){
              
              // For logarithmic drawing, cannot go down to zero
              if(logDeltaR){
                hSignalFakeToReflectedConeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
              }
              
              hSignalFakeToReflectedConeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->SetLineColor(color[iNormalization]);
              if(iNormalization == 0){
                drawer->SetGridY(true);
                hSignalFakeToReflectedConeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->GetYaxis()->SetRangeUser(ratioZoomLow, ratioZoomHigh);
                drawer->DrawHistogramToLowerPad(hSignalFakeToReflectedConeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization], "#Deltar", "Reflected / SF", " ");
                drawer->SetGridY(false);
              } else {
                hSignalFakeToReflectedConeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->Draw("same");
              }
              
            } // Normalization region loop
            
          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Energy-energy correlator loop
    
  } // Background to reflected cone ratio
  
  // Third set of plots: compare the signal distribution to total distribution from which the reflected cone background has been subtracted
  if(drawSignalToSubtrctedRatio){
    
    // Loop over all selevted histograms
    for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
      
      // Only read the selected energy-energy correlator types
      if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
      
      for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
        
        // Set the centrality information for legends and figure saving
        centralityString = Form("Cent: %.0f-%.0f%%",histograms->GetCentralityBinBorder(iCentrality),histograms->GetCentralityBinBorder(iCentrality+1));
        compactCentralityString = Form("_C=%.0f-%.0f",histograms->GetCentralityBinBorder(iCentrality),histograms->GetCentralityBinBorder(iCentrality+1));
        
        for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
          
          // Set the jet pT information for legends and figure saving
          if(iJetPt == histograms->GetNJetPtBinsEEC()){
            jetPtString = Form("Jet p_{T} > %.0f", histograms->GetCard()->GetJetPtCut());
            compactJetPtString = "";
          } else {
            jetPtString = Form("%.0f < jet p_{T} < %.0f", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1));
            compactJetPtString = Form("_J=%.0f-%.0f", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1));
          }
          
          for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
            
            // Set the track pT information for legends and figure saving
            trackPtString = Form("%.1f < track p_{T}",histograms->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString = Form("%.1f",histograms->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString.ReplaceAll(".","v");
            
            // Create a legend for the figure
            legend = new TLegend(0.19,0.09,0.39,0.94);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->AddEntry((TObject*) 0, histograms->GetCard()->GetAlternativeDataType().Data(), "");
            legend->AddEntry((TObject*) 0, centralityString.Data(),"");
            legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
            legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
            
            // Create a new canvas for the plot
            drawer->CreateSplitCanvas();
            
            // Logarithmic EEC axis
            if(logEEC) drawer->SetLogY(true);
            
            // For logarithmic drawing, cannot go down to zero
            if(logDeltaR){
              hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
            }
            
            // Set good y-ranges for plotting
            hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.001, 50);
            
            // Draw the background histogram to the upper canvas
            hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetLineColor(kBlack);
            drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt], "#Deltar", histograms->GetEnergyEnergyCorrelatorAxisName(iEnergyEnergyCorrelator), " ");
            legend->AddEntry(hEnergyEnergyCorrelatorSignal[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt], "Signal pairs", "l");
            
            // Draw different normalization regions for the reflected cone histogram to the same plot
            for(int iNormalization = 0; iNormalization < nBackgroundNormalizationBins; iNormalization++){
              
              // For logarithmic drawing, cannot go down to zero
              if(logDeltaR){
                hBackgroundSubtracted[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
              }
              
              hBackgroundSubtracted[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->SetLineColor(color[iNormalization]);
              hBackgroundSubtracted[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->Draw("same");
              legend->AddEntry(hBackgroundSubtracted[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization], Form("BG subtracted, #Deltar > %.3f", deltaRBinBorders[nDeltaRBins-(nBackgroundNormalizationBins-iNormalization)]), "l");
              
            } // Normalization region loop
            
            // Draw the legend to the upper pag
            legend->Draw();
            
            drawer->SetLogY(false);  // Linear scale for the ratio
            
            // Draw the ratios to the lower pad
            for(int iNormalization = 0; iNormalization < nBackgroundNormalizationBins; iNormalization++){
              
              // For logarithmic drawing, cannot go down to zero
              if(logDeltaR){
                hSignalToBackgroundSubtractedRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
              }
              
              hSignalToBackgroundSubtractedRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->SetLineColor(color[iNormalization]);
              if(iNormalization == 0){
                drawer->SetGridY(true);
                hSignalToBackgroundSubtractedRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->GetYaxis()->SetRangeUser(ratioZoomLow, ratioZoomHigh);
                drawer->DrawHistogramToLowerPad(hSignalToBackgroundSubtractedRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization], "#Deltar", "Signal / BG sub", " ");
                drawer->SetGridY(false);
              } else {
                hSignalToBackgroundSubtractedRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt][iNormalization]->Draw("same");
              }
              
            } // Normalization region loop
            
          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Energy-energy correlator loop
    
  } // Background to reflected cone ratio
  
}
