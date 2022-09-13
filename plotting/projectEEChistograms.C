#include "EECCard.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECHistogramManager.h"

#include <bitset>

/*
 * Macro for projecting the histograms needed in the energy-energy correlator analysis from the THnSparses
 *
 *  Arguments:
 *   TString inputFileName = File from which the histograms are read
 *   const char* outputFileName = If we are producing output file, name of the output file
 *   int histogramSelection = If > 0, select a preset group of histograms. Intended to be used for easier production of output files.
 */
void projectEEChistograms(TString inputFileName = "veryCoolData.root", const char* outputFileName = "veryCoolData_processed.root", int histogramSelection = 255){

  // Print the file name to console
  cout << "Projecting histograms histograms from " << inputFileName.Data() << endl;
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // If we write a file, define the output name and write mode
  const char* fileWriteMode = "UPDATE";
  
  // Choose which figure sets to draw
  bool loadEventInformation = false;
  bool loadJets = false;
  bool loadTracks = false;
  bool loadUncorrectedTracks = false;
  bool loadEnergyEnergyCorrelators = true;
  bool loadEnergyEnergyCorrelatorsJetPt = false;
  bool loadEnergyEnergyCorrelatorsUncorrected = false;
  bool loadEnergyEnergyCorrelatorsJetPtUncorrected = false;
  bool loadJetPtClosure = false;
  
  /*
   * Loading only selected histograms. Done with bitwise check of an integer
   *
   *  Bit 0 = Load event information histograms (to set: 1)
   *  Bit 1 = Load jet histograms (to set: 2)
   *  Bit 2 = Load track histograms (to set: 4)
   *  Bit 3 = Load uncorrected track histograms (to set: 8)
   *  Bit 4 = Load energy-energy correlator histograms (to set: 16)
   *  Bit 5 = Load jet pT weighted energy-energy correlator histograms (to set: 32)
   *  Bit 6 = Load uncorrected energy-energy correlator histograms (to set: 64)
   *  Bit 7 = Load uncorrected jet pT weighted energy-energy correlator histograms (to set: 128)
   *  Bit 8 = Load jet pT closure histograms (to set: 256)
   */
  if(histogramSelection > 0){
    std::bitset<9> bitChecker(histogramSelection);
    loadEventInformation = bitChecker.test(0);
    loadJets = bitChecker.test(1);
    loadTracks = bitChecker.test(2);
    loadUncorrectedTracks = bitChecker.test(3);
    loadEnergyEnergyCorrelators = bitChecker.test(4);
    loadEnergyEnergyCorrelatorsJetPt = bitChecker.test(5);
    loadEnergyEnergyCorrelatorsUncorrected = bitChecker.test(6);
    loadEnergyEnergyCorrelatorsJetPtUncorrected = bitChecker.test(7);
    loadJetPtClosure = bitChecker.test(8);
  }
  
  // ====================================================
  //  Binning configuration for the projected histograms
  // ====================================================
  
  // Option to read all the binning information from EECCard used to create the file
  const bool readCentralityBinsFromFile = false;
  const bool readTrackPtBinsFromFile = true;
  const bool readEECJetPtBinsFromFile = true;
  const bool readEECTrackPtBinsFromFile = true;
  
  // If not reading the bins from the file, manually define new bin borders
  const int nCentralityBins = 4;
  const int nTrackPtBins = 7;
  const int nJetPtBinsEEC = 7;
  const int nTrackPtBinsEEC = 6;
  double centralityBinBorders[nCentralityBins+1] = {0,10,30,50,90};   // Bin borders for centrality
  double trackPtBinBorders[nTrackPtBins+1] = {0.7,1,2,3,4,8,12,300};  // Bin borders for track pT
  double jetPtBinBordersEEC[nJetPtBinsEEC+1] = {120,140,160,180,200,300,500,5030}; // Bin borders for jet pT in energy-energy correlator histograms
  double trackPtBinBordersEEC[nTrackPtBinsEEC+1] = {0.7,1,2,3,4,6,300}; // Bin borders for track pT in energy-energy correlator histograms
  
  // Projected bin range
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = nCentralityBins-1;
  
  int firstDrawnTrackPtBin = 0;
  int lastDrawnTrackPtBin = nTrackPtBins-1;
  
  int firstDrawnJetPtBinEEC = 0;
  int lastDrawnJetPtBinEEC = nJetPtBinsEEC-1;
  
  int firstDrawnTrackPtBinEEC = 0;
  int lastDrawnTrackPtBinEEC = nTrackPtBinsEEC-1;
  
  
  // Jet flavor selection
  int jetFlavor = 0;   // Select jet flavor for MC: 0 = Any, 1 = Quark, 2 = Gluon. For data 0 = All vz, 1 = Negative vz, 2 = Positive vz
  
  
  // ==================================================================
  // ===================== Configuration ready ========================
  // ==================================================================
  
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
  
  // Remove centrality selection from pp data
  if(collisionSystem.Contains("pp")){
    lastDrawnCentralityBin = 0;
    centralityBinBorders[0] = -0.5;
  }
  
  // If we change the binning, save the new binning to the card
  if(!readCentralityBinsFromFile) card->AddVector(EECCard::kCentralityBinEdges,nCentralityBins+1,centralityBinBorders);
  if(!readTrackPtBinsFromFile) card->AddVector(EECCard::kTrackPtBinEdges,nTrackPtBins+1,trackPtBinBorders);
  if(!readEECJetPtBinsFromFile) card->AddVector(EECCard::kJetPtBinEdgesEEC,nJetPtBinsEEC+1,jetPtBinBordersEEC);
  if(!readEECTrackPtBinsFromFile) card->AddVector(EECCard::kTrackPtBinEdgesEEC,nTrackPtBinsEEC+1,trackPtBinBordersEEC);
  
  // Add information about the used input files to the card
  card->AddFileName(EECCard::kInputFileName,inputFileName);
  
  // ============================ //
  //     EECHistogramManager    //
  // ============================ //
    
  // Create and setup a new histogram manager to project and handle the histograms
  EECHistogramManager *histograms = new EECHistogramManager(inputFile,card);
  
  // Set which histograms to project from the input file
  histograms->SetLoadEventInformation(loadEventInformation);
  histograms->SetLoadJetHistograms(loadJets);
  histograms->SetLoadTracks(loadTracks);
  histograms->SetLoadTracksUncorrected(loadUncorrectedTracks);
  histograms->SetLoadEnergyEnergyCorrelators(loadEnergyEnergyCorrelators);
  histograms->SetLoadEnergyEnergyCorrelatorsJetPt(loadEnergyEnergyCorrelatorsJetPt);
  histograms->SetLoadEnergyEnergyCorrelatorsUncorrected(loadEnergyEnergyCorrelatorsUncorrected);
  histograms->SetLoadEnergyEnergyCorrelatorsJetPtUncorrected(loadEnergyEnergyCorrelatorsJetPtUncorrected);
  histograms->SetLoad2DHistograms(true);
  histograms->SetLoadJetPtClosureHistograms(loadJetPtClosure);
  histograms->SetJetFlavor(jetFlavor);

  // Set the binning information
  histograms->SetCentralityBins(readCentralityBinsFromFile,nCentralityBins,centralityBinBorders,true);
  if(!readCentralityBinsFromFile) histograms->SetCentralityBinRange(firstDrawnCentralityBin,lastDrawnCentralityBin);
  histograms->SetTrackPtBins(readTrackPtBinsFromFile,nTrackPtBins,trackPtBinBorders,true);
  if(!readTrackPtBinsFromFile) histograms->SetTrackPtBinRange(firstDrawnTrackPtBin,lastDrawnTrackPtBin);
  histograms->SetJetPtBinsEEC(readEECJetPtBinsFromFile,nJetPtBinsEEC,jetPtBinBordersEEC,true);
  if(!readEECJetPtBinsFromFile) histograms->SetJetPtBinRangeEEC(firstDrawnJetPtBinEEC,lastDrawnJetPtBinEEC);
  histograms->SetTrackPtBinsEEC(readEECTrackPtBinsFromFile,nTrackPtBinsEEC,trackPtBinBordersEEC,true);
  if(!readEECTrackPtBinsFromFile) histograms->SetTrackPtBinRangeEEC(firstDrawnTrackPtBinEEC,lastDrawnTrackPtBinEEC);
  
  // Project the one dimensional histograms from the THnSparses
  histograms->LoadHistograms();
  
  // Save the histograms to an output file
  histograms->Write(outputFileName,fileWriteMode);
  
}
