#include "EECCard.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECHistogramManager.h"
#include "EECDrawer.h"

#include <bitset>

/*
 * Macro for projecting the histograms needed in the energy-energy correlator analysis from the THnSparses
 *
 *  Arguments:
 *   TString inputFileName = File from which the histograms are read
 *   int histogramSelection = If > 0, select bitwise histograms. Intended to be used for easier production of output files.
 */
void plotEEChistograms(TString inputFileName = "veryCoolData_processed.root", int histogramSelection = 16){

  // Print the file name to console
  cout << "Projecting histograms histograms from " << inputFileName.Data() << endl;
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // If we write a file, define the output name and write mode
  const char* fileWriteMode = "UPDATE";
  
  // Choose which figure sets to draw
  bool drawEventInformation = false;
  bool drawJets = false;
  bool drawTracks = false;
  bool drawUncorrectedTracks = false;
  bool drawEnergyEnergyCorrelators = true;
  bool drawEnergyEnergyCorrelatorsJetPt = false;
  bool drawEnergyEnergyCorrelatorsUncorrected = false;
  bool drawEnergyEnergyCorrelatorsJetPtUncorrected = false;
  
  /*
   * Drawing only selected histograms. Done with bitwise check of an integer
   *
   *  Bit 0 = Load event information histograms (to set: 1)
   *  Bit 1 = Load jet histograms (to set: 2)
   *  Bit 2 = Load track histograms (to set: 4)
   *  Bit 3 = Load uncorrected track histograms (to set: 8)
   *  Bit 4 = Load energy-energy correlator histograms (to set: 16)
   *  Bit 5 = Load jet pT weighted energy-energy correlator histograms (to set: 32)
   *  Bit 6 = Load uncorrected energy-energy correlator histograms (to set: 64)
   *  Bit 7 = Load uncorrected jet pT weighted energy-energy correlator histograms (to set: 128)
   */
  if(histogramSelection > 0){
    std::bitset<8> bitChecker(histogramSelection);
    drawEventInformation = bitChecker.test(0);
    drawJets = bitChecker.test(1);
    drawTracks = bitChecker.test(2);
    drawUncorrectedTracks = bitChecker.test(3);
    drawEnergyEnergyCorrelators = bitChecker.test(4);
    drawEnergyEnergyCorrelatorsJetPt = bitChecker.test(5);
    drawEnergyEnergyCorrelatorsUncorrected = bitChecker.test(6);
    drawEnergyEnergyCorrelatorsJetPtUncorrected = bitChecker.test(7);
  }
  
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
  int lastDrawnJetPtBinEEC = nJetPtBinsEEC-1; // Note: Jets integrated over all pT ranges are in nJetPtBinsEEC bin
  
  int firstDrawnTrackPtBinEEC = 2;
  int lastDrawnTrackPtBinEEC = 2;
  
  // Remove centrality selection from pp data
  if(collisionSystem.Contains("pp")){
    lastDrawnCentralityBin = 0;
  }
  
  // Figure saving
  const bool saveFigures = false;
  const char* figureFormat = "png";
  TString figureNameSuffix = "_test";
  
  // Logarithmic scales for figures
  const bool logPt = true;          // pT distributions
  const bool logDeltaR = true;      // DeltaR axis for energy-energy correlators
  const bool logEEC = true;         // EEC axis for energy-energy correlators
  
  // Plotting style for 2D and 3D plots
  const int colorPalette = kLightTemperature;  // kRainBow kTemperatureMap kLightTemperature
  const char* style2D = "colz";
  const char* style3D = "surf1";
  
  // Select the style of histograms drawn for energy-energy correlators
  const bool drawIndividualEnergyEnergyCorrelators = false;
  const bool drawEnergyEnergyCorrelatorsForConstantJetPt = false;
  const bool drawEnergyEnergyCorrelatorsForConstantTrackPt = false;
  bool drawEnergyEnergyCorrelatorsSubevent = true;
  
  // If the collision system in not PbPb MC, we cannot draw subevent decomposition
  if(!collisionSystem.Contains("PbPb MC")){
    drawEnergyEnergyCorrelatorsSubevent = false;
  }

  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
  // ============================ //
  //     EECHistogramManager    //
  // ============================ //
    
  // Create and setup a new histogram manager to project and handle the histograms
  EECHistogramManager *histograms = new EECHistogramManager(inputFile,card);
  
  // Set which histograms to draw from the input file
  histograms->SetLoadEventInformation(drawEventInformation);
  histograms->SetLoadJetHistograms(drawJets);
  histograms->SetLoadTracks(drawTracks);
  histograms->SetLoadTracksUncorrected(drawUncorrectedTracks);
  histograms->SetLoadEnergyEnergyCorrelators(drawEnergyEnergyCorrelators);
  histograms->SetLoadEnergyEnergyCorrelatorsJetPt(drawEnergyEnergyCorrelatorsJetPt);
  histograms->SetLoadEnergyEnergyCorrelatorsUncorrected(drawEnergyEnergyCorrelatorsUncorrected);
  histograms->SetLoadEnergyEnergyCorrelatorsJetPtUncorrected(drawEnergyEnergyCorrelatorsJetPtUncorrected);
  histograms->SetLoad2DHistograms(true);
  
  histograms->SetCentralityBinRange(firstDrawnCentralityBin,lastDrawnCentralityBin);
  histograms->SetTrackPtBinRange(firstDrawnTrackPtBin,lastDrawnTrackPtBin);
  histograms->SetJetPtBinRangeEEC(firstDrawnJetPtBinEEC,lastDrawnJetPtBinEEC);
  histograms->SetTrackPtBinRangeEEC(firstDrawnTrackPtBinEEC,lastDrawnTrackPtBinEEC);
  
  // Load the histograms from the file
  histograms->LoadProcessedHistograms();
  
  
  // ============================ //
  //           EECDrawer          //
  // ============================ //
  
  // Create a new EECDrawer
  EECDrawer *resultDrawer = new EECDrawer(histograms);

  // Set which histograms to draw and the drawing style to use
  resultDrawer->SetDrawEventInformation(drawEventInformation);
  resultDrawer->SetDrawJetHistograms(drawJets);
  resultDrawer->SetDrawTracks(drawTracks);
  resultDrawer->SetDrawTracksUncorrected(drawUncorrectedTracks);
  resultDrawer->SetDrawEnergyEnergyCorrelor(drawEnergyEnergyCorrelators);
  resultDrawer->SetDrawEnergyEnergyCorrelorJetPt(drawEnergyEnergyCorrelatorsJetPt);
  resultDrawer->SetDrawEnergyEnergyCorrelorUncorrected(drawEnergyEnergyCorrelatorsUncorrected);
  resultDrawer->SetDrawEnergyEnergyCorrelorJetPtUncorrected(drawEnergyEnergyCorrelatorsJetPtUncorrected);
  
  resultDrawer->SetDrawSingleEnergyEnergyCorrelators(drawIndividualEnergyEnergyCorrelators);
  resultDrawer->SetDrawEnergyEnergyCorrelatorsForConstantJetPt(drawEnergyEnergyCorrelatorsForConstantJetPt);
  resultDrawer->SetDrawEnergyEnergyCorrelatorsForConstantTrackPt(drawEnergyEnergyCorrelatorsForConstantTrackPt);
  resultDrawer->SetDrawEnergyEnergyCorrelatorsSubevent(drawEnergyEnergyCorrelatorsSubevent);
  
  resultDrawer->SetSaveFigures(saveFigures,figureFormat,figureNameSuffix);
  resultDrawer->SetLogPt(logPt);
  resultDrawer->SetLogDeltaR(logDeltaR);
  resultDrawer->SetLogEEC(logEEC);
  resultDrawer->SetDrawingStyles(colorPalette,style2D,style3D);
  
  // Draw the selected histograms
  resultDrawer->DrawHistograms();
  
}
