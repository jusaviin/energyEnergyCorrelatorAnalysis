#ifndef EECCOMPARINGDRAWER_H
#define EECCOMPARINGDRAWER_H

// C++ includes
#include <iostream>
#include <tuple>      // For returning several arguments in a transparent manner

// Root includes
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>

// Own includes
#include "EECCard.h"
#include "JDrawer.h"
#include "EECHistogramManager.h"
#include "AlgorithmLibrary.h"

/*
 * Class for comparing histograms from different data files for the energy-energy correlator analysis
 */
class EECComparingDrawer {
  
private:
  static const int knMaxRatios = 10;  // Maximum number or ratio plots per one canvas
  int fColors[knMaxRatios] = {kRed,kBlue,kMagenta,kCyan,kGreen+4,kOrange,kViolet+3,kPink-7,kSpring+3,kAzure-7};
  
  void SetHistogramStyle(TH1 *histogram, double rangeX, const char* xTitle); // Style setting for histograms in big canvas
  
public:
  
  EECComparingDrawer(EECHistogramManager *baseHistograms);  // Constructor
  ~EECComparingDrawer();                 // Destructor
  
  void DrawHistograms();          // Draw the histograms
  
  // Add histograms to draw together with base histograms
  void AddHistogramToDraw(EECHistogramManager *additionalHistogram);
  void AddLegendComment(TString comment);
  
  // Setter for event information
  void SetDrawEventInformation(const bool drawOrNot); // Setter for drawing event information histograms
  
  // Setters jets
  void SetDrawJets(const bool drawOrNot);        // Setter for drawing all jet histograms
  
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
  void SetDrawParticleDensityAroundJetAxisPtBinned(const bool drawOrNot);   // Setter for drawing pT binned particle density around jet axis
  void SetDrawParticlePtDensityAroundJetAxisPtBinned(const bool drawOrNot); // Setter for drawing pT binned particle pT density around jet axis
  void SetDrawAllParticleDensitiesAroundJetAxis(const bool drawRegular, const bool drawPt, const bool drawPtBinned, const bool drawPtWeightedPtBinned); // Setter for drawing all particle densities around jet axis
  
  // Setters for energy-energy correaltors
  void SetDrawEnergyEnergyCorrelator(const bool drawOrNot);                 // Setter for drawing energy-energy correlators
  void SetDrawEnergyEnergyCorrelatorJetPt(const bool drawOrNot);            // Setter for drawing energy-energy correlators
  void SetDrawEnergyEnergyCorrelatorUncorrected(const bool drawOrNot);      // Setter for drawing energy-energy correlators
  void SetDrawEnergyEnergyCorrelatorJetPtUncorrected(const bool drawOrNot); // Setter for drawing energy-energy correlators
  void SetDrawAllEnergyEnergyCorrelators(const bool drawRegular, const bool drawJetPt, const bool drawUncorrected, const bool drawJetPtUncorrected);  // Setter for drawing all energy-energy correlators
  
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
  void SetSaveFigures(const bool saveOrNot, const char *format, const char *comment = "");  // Setter for saving the figures to a file
  void SetLogPt(const bool isLog);      // Setter for logarithmic jet pT axis
  void SetLogDeltaR(const bool isLog);  // Setter for logarithmic deltaR axis in energy-energy correlators
  void SetLogEEC(const bool isLog);     // Setter for logarithmic EEC axis in energy-energy correlators
  void SetLogAxes(const bool pt, const bool deltaR, const bool eec);  // Setter for logarithmic axes
  void SetLogParticleDensity(const bool isLog); // Setter for logarithmic y-axis in particle density histograms
  void SetManualLegend(const bool manualLegend);  // Setter for manual legend setting
  void SetAddSystemToLegend(const bool addSystem); // Setter for adding collision system to the legend
  void SetAddEnergyToLegend(const bool addSystem); // Setter for adding collision energy to the legend
  
  // Setters for drawing style
  void SetLineWidth(const int lineWidth); // Setter for line width in histograms
  
  // Setter for scaling and rebinning the histograms
  void SetApplyScaling(const int applyScaling); // Set if we should scale the histograms with their integral before comparing them
  void SetJetPtRebin(const bool doRebin);  // Tell if we want to rebin the jet pT histograms
  
  // Setters for ratio plots
  void SetUseDifferenceInRatioPlot(const bool useDifference);  // Setter for plotting difference instead of ratio to lower pad
  void SetRatioZoomMin(const double minValue);  // Setter for minimum value of y-axis in ratio plots
  void SetRatioZoomMax(const double maxValue);  // Setter for maximum value of y-axis in ratio plots
  void SetRatioZoom(const double minValue, const double maxValue);  // Setter for y-axis values in ratio plots
  void SetRatioLabel(TString label);        // Setter for y-axis label in ratio plots
  
  // Setters for drawing style and colors
  void SetColorPalette(const int color);     // Setter for color palette
  void SetDrawingStyle2D(const char* style); // Setter for 2D drawing style
  void SetDrawingStyle3D(const char* style); // Setter for 3D drawing style
  void SetDrawingStyles(const int color, const char* style2D, const char* style3D); // Setter for drawing styles
  
private:
  
  // Data members
  JDrawer *fDrawer;                       // JDrawer for drawing the histograms
  EECHistogramManager *fBaseHistograms; // Histograms with respect to which ratios are takes
  EECHistogramManager *fAddedHistograms[knMaxRatios];  // Histograms drawn together with the base histogram
  int fnAddedHistograms;                  // Number of histograms added for drawing
  
  // ==============================================================
  // == Helper variables to store histograms before drawing them ==
  // ==============================================================
  TH1D *fMainHistogram;
  TH1D *fComparisonHistogram[knMaxRatios];
  TH1D *fRatioHistogram[knMaxRatios];
  
  // ==============================================
  // ======== Flags for histograms to draw ========
  // ==============================================
  
  bool fDrawEventInformation;      // Draw the event information histograms
  bool fDrawJets;                  // Draw the single jet histograms
  bool fDrawTracks[EECHistogramManager::knTrackCategories];                                       // Draw the track histograms
  bool fDrawMultiplicityInJetCone[EECHistogramManager::knMultiplicityInJetConeTypes];             // Draw the multiplicity in jet cone histograms
  bool fDrawParticleDensityAroundJets[EECHistogramManager::knParticleDensityAroundJetAxisTypes];  // Draw the particle densities around jet axis
  bool fDrawEnergyEnergyCorrelators[EECHistogramManager::knEnergyEnergyCorrelatorTypes];          // Draw the energy-energy correlators
  bool fTrackHistogramDrawn;                      // Flag telling that at least one track histogram will be drawn
  bool fEnergyEnergyCorrelatorHistogramDrawn;     // Flag telling that at least one energy-energy correlator histogram will be drawn
  
  bool fDrawPairingType[EECHistograms::knPairingTypes];                   // Select which pairing types (signal+signal, signal+reflected cone, ref+ref) to draw
  bool fDrawSubeventType[EECHistograms::knSubeventTypes+1];               // Select which subevent types (Pythia, Hydjet, all) to draw
  bool fDrawSubeventCombination[EECHistograms::knSubeventCombinations+1]; // Select which subevent combinations to draw
  
  // ==============================================
  // ============== Drawing settings ==============
  // ==============================================
  
  // Choose if you want to write the figures to pdf file
  bool fSaveFigures;           // Flag for saving the figures to file
  const char* fFigureFormat;   // Format in which the figures are saved
  const char* fFigureComment;  // Comment added to figure name
  bool fManualLegend;          // Flag for using manual instead of automatic legend
  
  // Choose if we should scale the histograms before comparing them
  int fApplyScaling;
  double fScalingFactors[knMaxRatios+1];
  
  // Logarithmic scales for figures
  bool fLogPt;          // jet pT distributions
  bool fLogDeltaR;      // logarithmic deltaR axis for the energy-energy correlators
  bool fLogEEC;         // logarithmic EEC axis for the energy-energy correaltors
  bool fLogParticleDensity; // logarithmic y-axis for particle densities
  
  // Zooming for ratio plots
  bool fUseDifferenceInsteadOfRatio;  // Instead of ratio, draw the difference of the two distributions
  double fRatioZoomMin;               // Lower y-axis boundary in ratio plots
  double fRatioZoomMax;               // Upper y-axis boundary in ratio plots
  TString fRatioLabel;                // Label given to ratio plots y-axes
  
  // Plotting style for 2D and 3D plots
  int fColorPalette;      // Used color palatte for drawing
  const char* fStyle2D;   // Used option for two-dimensional drawing style
  const char* fStyle3D;   // Used option for three-dimensional drawing style
  
  // Rebinning
  bool fRebinJetPt;       // Rebin the single jet pT distributions
  
  // Drawing style
  bool fAddSystemToLegend; // Add the collision system to the legend
  bool fAddEnergyToLegend; // Add the collision energy to the legend
  int fLineWidth;          // Line width for the drawn histgrams
  
  // Drawn bins
  int fFirstDrawnCentralityBin;  // First centrality bin that is drawn
  int fLastDrawnCentralityBin;   // Last centrality bin that is drawn
  int fFirstDrawnTrackPtBin;     // First track pT bin that is drawn
  int fLastDrawnTrackPtBin;      // Last track pT bin that is drawn
  int fFirstDrawnJetPtBinEEC;    // First energy-energy correlator jet pT bin that is drawn
  int fLastDrawnJetPtBinEEC;     // Last energy-energy correlator jet pT bin that is drawn
  int fFirstDrawnTrackPtBinEEC;  // First energy-energy correlator track pT bin that is drawn
  int fLastDrawnTrackPtBinEEC;   // Last energy-energy correlator track pT bin that is drawn
  
  // Comments given for legends
  TString fLegendComment[knMaxRatios+1];
  
  // Methods for drawing
  void DrawEventInformation();                  // Draw event information histograms
  void DrawJetHistograms();                     // Draw jet histograms
  void DrawTrackHistograms();                   // Draw track histograms
  void DrawMultiplicityInJetCone();             // Draw the multiplicities within the jet cone
  void DrawParticleDensityAroundJetAxis();      // Draw particle densities around the jet axis
  void DrawEnergyEnergyCorrelatorHistograms();  // Draw energy-energy correlator histograms
  void SetupLegend(TLegend *legend, TString centralityString = "", TString trackString = "", TString asymmetryString = "", TString extraString = "", TString additionalString = ""); // Common legend style setup for figures
  void SaveFigure(TString figureName, TString centralityString = "", TString trackPtString = "", TString correlationTypeString = "", TString deltaPhiString = ""); // Save the figure from current canvas to file
  void PrepareRatio(TString name, int rebin, int bin1 = 0, int bin2 = 0, int bin3 = 0, int bin4 = 0, int bin5 = 0, int bin6 = 0, double minRange = 1000, double maxRange = -1000); // Prepare the ratio histograms out of input histograms
  void DrawToUpperPad(const char* xTitle, const char* yTitle, const bool logX = false, const bool logY = false); // Draw the histograms to the same figure in the upper pad of JDrawer
  void DrawToLowerPad(const char* xTitle, const char* yTitle, const double zoomMin, const double zoomMax, const bool logX = false); // Draw the ratios to the lower pad of the JDrawer
  void ZoomToRegion(const double maxZoomValue, const int nZoomBins, const double scaleFactor, const bool bothSides, const bool asymmetricZoom);  // Zoom the y-axis scale to the specified region of the distribution
  std::tuple<double,double> GetHistogramAverageAndDifferenceInRegion(TH1D *histogram, const double maxZoomValue, const int nZoomBins, const bool bothSides); // Get the average and absolute difference from a histogram in the specific area
  
  // Find the scaling factors that scale the integral of the distributions to one
  void FindScalingFactors(const char* histogramName, const int bin1, const int bin2 = 0, const int bin3 = 0, const int bin4 = 0, const int bin5 = 0, const int bin6 = 0);
  
  // Check flags to see if we are drawing certain types of histograms
  void CheckFlags();                       // Check flags for all histogram categories
  void CheckFlagsTrack();                  // Check flags for track histograms
  void CheckFlagsEnergyEnergyCorrelator(); // Check flags for energy-energy correlators
  
};

#endif
