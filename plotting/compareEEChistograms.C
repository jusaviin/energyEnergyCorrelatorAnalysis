#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECComparingDrawer.h"
#include "EECCard.h"

/*
 * Macro for configuring the EECComparingDrawer and defining which histograms are drawn
 */
void compareEEChistograms(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Define the used data files, and a comment describing the data in each file
  const int nDatasets = 2;
  TString inputFileName[] = { "data/eecAnalysis_akFlowJets_moreLowTrackPtCuts_processed_2022-09-15.root", "data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJet_firstLook_2022-09-15_someMissing_processed.root", "data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJet_firstLook_withTrigger_2022-09-15_someMissing_processed.root"};
  
  TString legendComment[] = {"PbPb data", "Pythia+Hydjet", "Pythia+Hydjet, jet trigger"};
  
  // Try to open the files
  TFile *inputFile[nDatasets];
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    // Open the file for the given dataset
    inputFile[iDataset] = TFile::Open(inputFileName[iDataset]);
    if(inputFile[iDataset] == NULL){
      cout << "Error! The file " << inputFileName[iDataset].Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Please give a file that exists. Will not exacute the code" << endl;
      return;
    }
  }
   
  // Read the EECCard from the first file. This is used to gain binning information, and assumed to be the same in other files too
  EECCard *card = new EECCard(inputFile[0]);
  
  // Choose which figure sets to draw
  bool drawEventInformation = false;
  bool drawJets = false;
  bool drawTracks = false;
  bool drawUncorrectedTracks = false;
  bool drawEnergyEnergyCorrelators = true;
  bool drawEnergyEnergyCorrelatorsJetPt = false;
  bool drawEnergyEnergyCorrelatorsUncorrected = false;
  bool drawEnergyEnergyCorrelatorsJetPtUncorrected = false;
  
  // Choose if you want to write the figures to pdf file
  bool saveFigures = false;
  const char* figureFormat = "pdf";
  const char* figureComment = "_dataMCfirstLook";
  
  // Logarithmic scales for figures
  bool logPt = true;       // pT axis for jet
  bool logDeltaR = true;   // deltaR axis for energy-energy correlators
  bool logEEC = true;      // EEC axis for energy-energy correlators
  
  // Plotting style for 2D and 3D plots
  int colorPalette = kRainBow;
  const char* style2D = "colz";
  const char* style3D = "surf1";
  
  // Settings for ratios
  bool useDifferenceInsteadOfRatio = false;
  double minZoom = 0.1;
  double maxZoom = 1.9;
  TString ratioLabel = "Data / MC";
  bool manualLegend = false; // Set this true if you want to set legend manually in EECComparingDrawer.cxx code instead of using automatic legend generation
  
  // Scaling for histograms
  int scaleHistograms = 1; // 0 = Do not scale histograms. 1 = Scale integral to one. 2 = Scale average to one
  bool rebinJetPt = false;
  
  // ====================================================
  //  Binning configuration for the drawn histograms
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = card->GetNCentralityBins();
  const int nTrackPtBins = card->GetNTrackPtBins();
  const int nJetPtBinsEEC = card->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card->GetNTrackPtBinsEEC();
  
  // Default binning ranges for reference
  // centrality = {0,10,30,50,90}
  // track pT =  {0.7,1,2,3,4,8,12,300}
  // jet pT in energy-energy correlators = {120,140,160,180,200,300,500,5030}
  // track pT in energy-energy correlators = {0.7,1,2,3,4,6,300}
  
  // Bin range to be drawn
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = 0;
  
  int firstDrawnTrackPtBin = 0;
  int lastDrawnTrackPtBin = nTrackPtBins-1;
  
  int firstDrawnJetPtBinEEC = 0;
  int lastDrawnJetPtBinEEC = 0; // Note: Jets integrated over all pT ranges are in nJetPtBinsEEC bin
  
  int firstDrawnTrackPtBinEEC = 7;
  int lastDrawnTrackPtBinEEC = 7;
  
  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // ==================================================================
  // ================= Setup EECHistogramManagers =====================
  // ==================================================================
  EECHistogramManager *histograms[nDatasets];
  
  for(int iDataset = 0; iDataset < nDatasets; iDataset++){
    
    // Create a new histogram manager
    histograms[iDataset] = new EECHistogramManager(inputFile[iDataset]);

    // Set which histograms to draw from the input file
    histograms[iDataset]->SetLoadEventInformation(drawEventInformation);
    histograms[iDataset]->SetLoadJetHistograms(drawJets);
    histograms[iDataset]->SetLoadTracks(drawTracks);
    histograms[iDataset]->SetLoadTracksUncorrected(drawUncorrectedTracks);
    histograms[iDataset]->SetLoadEnergyEnergyCorrelators(drawEnergyEnergyCorrelators);
    histograms[iDataset]->SetLoadEnergyEnergyCorrelatorsJetPt(drawEnergyEnergyCorrelatorsJetPt);
    histograms[iDataset]->SetLoadEnergyEnergyCorrelatorsUncorrected(drawEnergyEnergyCorrelatorsUncorrected);
    histograms[iDataset]->SetLoadEnergyEnergyCorrelatorsJetPtUncorrected(drawEnergyEnergyCorrelatorsJetPtUncorrected);
    
    histograms[iDataset]->SetCentralityBinRange(firstDrawnCentralityBin,lastDrawnCentralityBin);
    histograms[iDataset]->SetTrackPtBinRange(firstDrawnTrackPtBin,lastDrawnTrackPtBin);
    histograms[iDataset]->SetJetPtBinRangeEEC(firstDrawnJetPtBinEEC,lastDrawnJetPtBinEEC);
    histograms[iDataset]->SetTrackPtBinRangeEEC(firstDrawnTrackPtBinEEC,lastDrawnTrackPtBinEEC);
    
    // Load the histograms from the file
    histograms[iDataset]->LoadProcessedHistograms();

  } // Loop over datasets
  
  // ==================================================================
  // ==== Setup EECComparingDrawer and draw the comparison plots ======
  // ==================================================================
  
  EECComparingDrawer *drawer = new EECComparingDrawer(histograms[0]);
  drawer->AddLegendComment(legendComment[0]);
  for(int i = 1; i < nDatasets; i++){
    drawer->AddHistogramToDraw(histograms[i]);
    drawer->AddLegendComment(legendComment[i]);
  }
  
  drawer->SetDrawEventInformation(drawEventInformation);
  drawer->SetDrawJets(drawJets);
  drawer->SetDrawAllTracks(drawTracks,drawUncorrectedTracks);
  drawer->SetDrawEnergyEnergyCorrelator(drawEnergyEnergyCorrelators);
  drawer->SetDrawEnergyEnergyCorrelatorJetPt(drawEnergyEnergyCorrelatorsJetPt);
  drawer->SetDrawEnergyEnergyCorrelatorUncorrected(drawEnergyEnergyCorrelatorsUncorrected);
  drawer->SetDrawEnergyEnergyCorrelatorJetPtUncorrected(drawEnergyEnergyCorrelatorsJetPtUncorrected);
  
  drawer->SetSaveFigures(saveFigures,figureFormat,figureComment);
  drawer->SetLogAxes(logPt, logDeltaR, logEEC);
  drawer->SetDrawingStyles(colorPalette,style2D,style3D);
  drawer->SetUseDifferenceInRatioPlot(useDifferenceInsteadOfRatio);
  drawer->SetRatioZoom(minZoom,maxZoom);
  drawer->SetRatioLabel(ratioLabel);
  drawer->SetApplyScaling(scaleHistograms);
  drawer->SetJetPtRebin(rebinJetPt);
  drawer->SetManualLegend(manualLegend);

  // Draw the selected histograms
  drawer->DrawHistograms();
  
}
