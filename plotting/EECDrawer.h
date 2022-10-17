#ifndef EECDRAWER_H
#define EECDRAWER_H

// C++ includes
#include <bitset>
#include <vector>

// Root includes
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TString.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>

// Own includes
#include "JDrawer.h"
#include "EECHistogramManager.h"

/*
 * Class for drawing the histograms produced in the energy-energy correlator analysis
 */
class EECDrawer {
  
public:
  
  EECDrawer(EECHistogramManager *inputHistograms);  // Constructor
  ~EECDrawer();                                     // Destructor
  
  void DrawHistograms();          // Draw the histograms
  
  // Setter for event information
  void SetDrawEventInformation(const bool drawOrNot); // Setter for drawing event information
  
  // Setter for jets
  void SetDrawJetHistograms(const bool drawOrNot);    // Setter for drawing jet histograms
  
  // Setters for tracks
  void SetDrawTracks(const bool drawOrNot);            // Setter for drawing tracks
  void SetDrawTracksUncorrected(const bool drawOrNot); // Setter for drawing uncorrected tracks
  void SetDrawAllTracks(const bool drawTracks, const bool drawUncorrected); // Setter for drawing all track histograms
  
  // Setters for multiplicity within the jet cone
  void SetDrawMultiplicityInJetCone(const bool drawOrNot);                  // Setter for drawing multiplicity within the jet cone
  void SetDrawMultiplicityInReflectedCone(const bool drawOrNot);            // Setter for drawing multiplicity within the reflected cone
  void SetDrawMultiplicityInJetConeUncorrected(const bool drawOrNot);       // Setter for drawing uncorrected multiplicity within the jet cone
  void SetDrawMultiplicityInReflectedConeUncorrected(const bool drawOrNot); // Setter for drawing uncorrected multiplicity within the reflected cone
  void SetDrawAllMultiplicitiesInJetCone(const bool regular, const bool reflectedCone, const bool regularUncorrected, const bool reflectedConeUncorrected); // Setter for drawing all multiplicity histograms within the jet cone
  
  // Setters for particle density around the jet axis
  void SetDrawParticleDensityAroundJetAxis(const bool drawOrNot);   // Setter for drawing particle density around jet axis
  void SetDrawParticlePtDensityAroundJetAxis(const bool drawOrNot); // Setter for drawing particle pT density around jet axis
  void SetDrawAllParticleDensitiesAroundJetAxis(const bool drawRegular, const bool drawPt); // Setter for drawing all particle densities around jet axis
  
  // Setters for energy-energy correlators
  void SetDrawEnergyEnergyCorrelor(const bool drawOrNot);                  // Setter for drawing energy-energy correlator
  void SetDrawEnergyEnergyCorrelorJetPt(const bool drawOrNot);             // Setter for drawing jet pT weighted energy-energy correlator
  void SetDrawEnergyEnergyCorrelorUncorrected(const bool drawOrNot);       // Setter for drawing uncorrected energy-energy correlator
  void SetDrawEnergyEnergyCorrelorJetPtUncorrected(const bool drawOrNot);  // Setter for drawing uncorrected jet pT weighted energy-energy correlator
  void SetDrawAllEnergyEnergyCorrelors(const bool drawRegular, const bool drawJetPt, const bool drawUncorrected, const bool drawJetPtUncorrected); // Setter for drawing all energy-energy correlators
  
  // Setters for drawing different figures for energy-energy correlators
  void SetDrawSingleEnergyEnergyCorrelators(const bool drawOrNot);             // Setter for drawing the individual energy-energy correlator histograms
  void SetDrawEnergyEnergyCorrelatorsForConstantJetPt(const bool drawOrNot);   // Setter for drawing all track pT cuts to the same figure for constant jet pT selection
  void SetDrawEnergyEnergyCorrelatorsForConstantTrackPt(const bool drawOrNot); // Setter for drawing all jet pT selections to the same figure for constant track pT cut
  void SetDrawEnergyEnergyCorrelatorsSubevent(const bool drawOrNot);           // Setter for drawing subevent decomposition for energy-energy correlators
  
  // Setters for drawing different pairing types
  void SetDrawSameJetEnergyEnergyCorrelators(const bool drawOrNot);              // Setter for drawing same jet energy-energy correlators
  void SetDrawSignalReflectedConeEnergyEnergyCorrelators(const bool drawOrNot);  // Setter for drawing signal-reflected cone energy-energy correlators
  void SetDrawReflectedConeOnlyEnergyEnergyCorrelators(const bool drawOrNot);    // Setter for drawing reflected cone-reflected cone energy-energy correlators
  void SetDrawAllEnergyEnergyCorrelatorPairingTypes(const bool drawSameJet, const bool drawSignalReflectedCone, const bool drawReflectedConeOnly); // Setter for drawing all different energy-energy correlator pairing types
  
  // Setters for drawing different subevent types
  void SetDrawAllSubevents(const bool drawOrNot); // Setter for drawing histograms without subevent selection
  void SetDrawPythiaOnly(const bool drawOrNot);   // Setter for drawing only Pythia histograms
  void SetDrawHydjetOnly(const bool drawOrNot);   // Setter for drawing only Hydjet histograms
  void SetDrawAllSubeventTypes(const bool drawAll, const bool drawPythia, const bool drawHydjet); // Setter for drawing all subevent types
  
  void SetDrawAllCombinations(const bool drawOrNot); // Setter for drawing all pairing combinations
  void SetDrawSignalOnly(const bool drawOrNot);      // Setter for drawing Pythia+Pythia correlations for simulation
  void SetDrawSignalFake(const bool drawOrNot);      // Setter for drawing Pythia+Hydjet correlations for simulation
  void SetDrawFakeFake(const bool drawOrNot);        // Setter for drawing Hydjet+Hydjet correlations for simulation
  void SetDrawAllSubeventCombinations(const bool drawAll, const bool drawSignal, const bool drawSignalFake, const bool drawFakeFake); // Setter for drawing all subevent combinations for energy-energy correlators
  
  // Setters for figure saving and logarithmic axes
  void SetSaveFigures(const bool saveOrNot, const char *format, const TString suffix);  // Setter for saving the figures to a file
  void SetLogPt(const bool isLog);          // Setter for logarithmic pT axis
  void SetLogDeltaR(const bool isLog);      // Setter for logarithmic deltaR axis in energy-energy correlators
  void SetLogEEC(const bool isLog);         // Setter for logarithmic EEC axis in energy-energy correlators
  
  // Setters for drawing style and colors
  void SetColorPalette(const int color);     // Setter for color palette
  void SetDrawingStyle2D(const char* style); // Setter for 2D drawing style
  void SetDrawingStyle3D(const char* style); // Setter for 3D drawing style
  void SetDrawingStyles(const int color, const char* style2D, const char* style3D); // Setter for drawing styles
  
private:
  
  // Data members
  EECHistogramManager *fHistograms; // Manager for all the drawn histograms
  TString fSystemAndEnergy;           // Collision system (pp,PbPb,pp MC,PbPb MC,localTest) and energy
  TString fCompactSystemAndEnergy;    // Same a before but without white spaces and dots
  TString fFigureSaveNameAppend;      // Text that can be appended to standard figure naming scheme
  JDrawer *fDrawer;                   // JDrawer for drawing the histograms
  
  // ==============================================
  // ======== Flags for histograms to draw ========
  // ==============================================
  
  bool fDrawEventInformation;                                                                     // Draw the event information histograms
  bool fDrawJets;                                                                                 // Draw the jet histograms
  bool fDrawTracks[EECHistogramManager::knTrackCategories];                                       // Draw the track histograms
  bool fDrawMultiplicityInJetCone[EECHistogramManager::knMultiplicityInJetConeTypes];             // Draw the multiplicity in jet cone histograms
  bool fDrawParticleDensityAroundJets[EECHistogramManager::knParticleDensityAroundJetAxisTypes];  // Draw the particle densities around jet axis
  bool fDrawEnergyEnergyCorrelators[EECHistogramManager::knEnergyEnergyCorrelatorTypes];          // Draw the energy-energy correlator histograms
  
  bool fDrawIndividualEnergyEnergyCorrelators;          // Draw the individual energy-energy correlator histograms
  bool fDrawEnergyEnergyCorrelatorsForConstantJetPt;    // Draw all track pT cuts to the same figure for constant jet pT selection
  bool fDrawEnergyEnergyCorrelatorsForConstantTrackPt;  // Draw all jet pT selections to the same figure for constant track pT cut
  bool fDrawEnergyEnergyCorrelatorsSubevents;           // Draw subevent decomposition of energy-energy correlators
  
  bool fDrawPairingType[EECHistograms::knPairingTypes];                   // Select which pairing types to draw
  bool fDrawSubeventType[EECHistograms::knSubeventTypes+1];               // Select which subevent types to draw
  bool fDrawSubeventCombination[EECHistograms::knSubeventCombinations+1]; // Select which subevent combinations to draw
  
  // ==============================================
  // ============== Drawing settings ==============
  // ==============================================
  
  // Choose if you want to write the figures to pdf file
  bool fSaveFigures;           // Flag for saving the figures to file
  const char* fFigureFormat;   // Format in which the figures are saved
  
  // Logarithmic scales for figures
  bool fLogPt;          // pT distributions
  bool fLogDeltaR;      // DeltaR axis in energy-energy correlators
  bool fLogEEC;         // EEC axis in energy-energy correlators
  
  // Plotting style for 2D and 3D plots
  int fColorPalette;      // Used color palatte for drawing
  const char* fStyle2D;   // Used option for two-dimensional drawing style
  const char* fStyle3D;   // Used option for three-dimensional drawing style
  
  // Drawn bins
  int fFirstDrawnCentralityBin;  // First centrality bin that is drawn
  int fLastDrawnCentralityBin;   // Last centrality bin that is drawn
  int fFirstDrawnTrackPtBin;     // First track pT bin that is drawn
  int fLastDrawnTrackPtBin;      // Last track pT bin that is drawn
  int fFirstDrawnJetPtBinEEC;    // First energy-energy correlator jet pT bin that is drawn
  int fLastDrawnJetPtBinEEC;     // Last energy-energy correlator jet pT bin that is drawn
  int fFirstDrawnTrackPtBinEEC;  // First energy-energy correlator track pT bin that is drawn
  int fLastDrawnTrackPtBinEEC;   // Last energy-energy correlator track pT bin that is drawn
  
  // Methods for drawing
  void DrawEventInformation();              // Draw the event information histograms
  void DrawJetHistograms();                 // Draw jet histograms
  void DrawMultiplicityInJetCone();         // Draw the multiplicities within the jet cone
  void DrawParticleDensityAroundJetAxis();  // Draw particle densities around the jet axis
  void DrawTrackHistograms();               // Draw track histograms
  void DrawEnergyEnergyCorrelationHistograms();     // Draw the energy-energy correlation histograms
  void SetupLegend(TLegend *legend, TString centralityString = "", TString jetString = "", TString trackString = "", TString extraString = "", TString anotherString = ""); // Common legend style setup for figures
  void SaveFigure(TString figureName, TString centralityString = "", TString trackPtString = "", TString correlationTypeString = "", TString deltaPhiString = ""); // Save the figure from current canvas to file
  
};

#endif
