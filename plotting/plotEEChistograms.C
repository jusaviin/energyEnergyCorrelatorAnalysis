#include "EECCard.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECHistogramManager.h"
#include "EECDrawer.h"

#include <bitset>

/*
 * Macro for projecting the histograms needed in the energy-energy correlator analysis from the THnSparses
 *
 *  Arguments:
 *   TString inputFileName = File from which the histograms are read
 */
void plotEEChistograms(TString inputFileName = "veryCoolData_processed.root"){

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
  bool drawMultiplicityHistograms = false;
  bool drawParticleDensityAroundJets = false;
  bool drawParticlePtDensityAroundJets = false;
  bool drawParticleDensityAroundJetsPtBinned = false;
  bool drawParticlePtDensityAroundJetsPtBinned = false;
  bool drawMaxParticlePtWithinJetCone = false;
  bool drawMaxBackgroundParticlePtWithinJetCone = false;
  
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
  
  int firstDrawnJetPtBinEEC = nJetPtBinsEEC-1;
  int lastDrawnJetPtBinEEC = nJetPtBinsEEC-1; // Note: Jets integrated over all pT ranges are in nJetPtBinsEEC bin
  
  int firstDrawnTrackPtBinEEC = 0;
  int lastDrawnTrackPtBinEEC = nTrackPtBinsEEC-1;
  
  // Remove centrality selection from pp data
  if(collisionSystem.Contains("pp")){
    lastDrawnCentralityBin = 0;
  }
  
  // Select with multiplicity histograms to draw from all possible ones
  bool drawMultiplicityInJetCone = true && drawMultiplicityHistograms;
  bool drawMultiplicityInReflectedCone = false && drawMultiplicityHistograms;
  bool drawMultiplicityInJetConeUncorrected = false && drawMultiplicityHistograms;
  bool drawMultiplicityInReflectedConeUncorrected = false && drawMultiplicityHistograms;
  
  // Select track pairing type to be draw
  const bool drawSameJetEnergyEnergyCorrelator = true;       // Draw energy-energy correlator where tracks from the same jet are paired
  const bool drawSignalReflectedConeEnergyEnergyCorrelator = false; // Draw energy-energy correlator where tracks from jet cone are paired with tracks from reflected jet cone
  const bool drawReflectedConeOnlyEnergyEnergyCorrelator = false; // Draw energy-energy correlator where tracks from reflected jet cone are paired with tracks from reflected jet cone
  
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
  
  // Select the style of histograms drawn for particle density histograms
  const bool drawIndividualParticleDensities = true;
  const bool drawParticleDensitiesForConstantJetPt = false;
  
  // Select the style of histograms drawn for energy-energy correlators
  const bool drawIndividualEnergyEnergyCorrelators = false;
  const bool drawEnergyEnergyCorrelatorsForConstantJetPt = true;
  const bool drawEnergyEnergyCorrelatorsForConstantTrackPt = false;
  bool drawEnergyEnergyCorrelatorsSubevent = false;
  
  // Select which subevents to draw
  bool drawAllSubevents = true;   // Draw histograms without subevent selection
  bool drawPythiaOnly = false;    // Draw only Pythia histograms in Pythia+Hydjet simulation
  bool drawHydjetOnly = false;    // Draw only Hydjet histograms in Pythia+Hydjet simulation
  
  bool drawAllSubeventPairs = true;  // Draw energy-energy correlators without subevent selection
  bool drawSignalOnly = false;        // Draw Pythia+Pythia correlations from MC
  bool drawSignalFake = false;        // Draw Pythia+Hydjet correlations from MC
  bool drawFakeFake = false;          // Draw Hydjet+Hydjet correlations from MC
  
  // If the collision system in not PbPb MC, we cannot draw subevent decomposition
  if(!collisionSystem.Contains("PbPb MC")){
    drawEnergyEnergyCorrelatorsSubevent = false;
    
    drawPythiaOnly = false;
    drawHydjetOnly = false;
    drawAllSubevents = true;
    
    drawSignalOnly = false;
    drawSignalFake = false;
    drawFakeFake = false;
    drawAllSubeventPairs = true;
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
  histograms->SetLoadMultiplicityInJets(drawMultiplicityHistograms);
  histograms->SetLoadParticleDensityAroundJets(drawParticleDensityAroundJets);
  histograms->SetLoadParticlePtDensityAroundJets(drawParticlePtDensityAroundJets);
  histograms->SetLoadParticleDensityAroundJetsPtBinned(drawParticleDensityAroundJetsPtBinned);
  histograms->SetLoadParticlePtDensityAroundJetsPtBinned(drawParticlePtDensityAroundJetsPtBinned);
  histograms->SetLoadMaxParticlePtWithinJetCone(drawMaxParticlePtWithinJetCone || drawMaxBackgroundParticlePtWithinJetCone);
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
  
  resultDrawer->SetDrawMultiplicityInJetCone(drawMultiplicityInJetCone);
  resultDrawer->SetDrawMultiplicityInReflectedCone(drawMultiplicityInReflectedCone);
  resultDrawer->SetDrawMultiplicityInJetConeUncorrected(drawMultiplicityInJetConeUncorrected);
  resultDrawer->SetDrawMultiplicityInReflectedConeUncorrected(drawMultiplicityInReflectedConeUncorrected);
  
  resultDrawer->SetDrawParticleDensityAroundJetAxis(drawParticleDensityAroundJets);
  resultDrawer->SetDrawParticlePtDensityAroundJetAxis(drawParticlePtDensityAroundJets);
  resultDrawer->SetDrawParticleDensityAroundJetAxisPtBinned(drawParticleDensityAroundJetsPtBinned);
  resultDrawer->SetDrawParticlePtDensityAroundJetAxisPtBinned(drawParticlePtDensityAroundJetsPtBinned);
  
  resultDrawer->SetDrawSingleParticleDensityHistograms(drawIndividualParticleDensities);
  resultDrawer->SetDrawParticleDensityForConstantJetPt(drawParticleDensitiesForConstantJetPt);
  
  resultDrawer->SetDrawMaxParticlePtWithinJetCone(drawMaxParticlePtWithinJetCone);
  resultDrawer->SetDrawMaxBackgroundParticlePtWithinJetCone(drawMaxBackgroundParticlePtWithinJetCone);
  
  resultDrawer->SetDrawSingleEnergyEnergyCorrelators(drawIndividualEnergyEnergyCorrelators);
  resultDrawer->SetDrawEnergyEnergyCorrelatorsForConstantJetPt(drawEnergyEnergyCorrelatorsForConstantJetPt);
  resultDrawer->SetDrawEnergyEnergyCorrelatorsForConstantTrackPt(drawEnergyEnergyCorrelatorsForConstantTrackPt);
  resultDrawer->SetDrawEnergyEnergyCorrelatorsSubevent(drawEnergyEnergyCorrelatorsSubevent);
  
  resultDrawer->SetDrawSameJetEnergyEnergyCorrelators(drawSameJetEnergyEnergyCorrelator);
  resultDrawer->SetDrawSignalReflectedConeEnergyEnergyCorrelators(drawSignalReflectedConeEnergyEnergyCorrelator);
  resultDrawer->SetDrawReflectedConeOnlyEnergyEnergyCorrelators(drawReflectedConeOnlyEnergyEnergyCorrelator);
  
  resultDrawer->SetDrawAllSubeventTypes(drawAllSubevents, drawPythiaOnly, drawHydjetOnly);
  resultDrawer->SetDrawAllSubeventCombinations(drawAllSubeventPairs, drawSignalOnly, drawSignalFake, drawFakeFake);
  
  resultDrawer->SetSaveFigures(saveFigures,figureFormat,figureNameSuffix);
  resultDrawer->SetLogPt(logPt);
  resultDrawer->SetLogDeltaR(logDeltaR);
  resultDrawer->SetLogEEC(logEEC);
  resultDrawer->SetDrawingStyles(colorPalette,style2D,style3D);
  
  // Draw the selected histograms
  resultDrawer->DrawHistograms();
  
}
