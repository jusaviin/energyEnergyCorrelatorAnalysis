#ifndef EECHISTOGRAMMANAGER_H
#define EECHISTOGRAMMANAGER_H

// Root includes
#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TLegend.h>
#include <TStyle.h>

// Own includes
#include "EECCard.h"
#include "../src/EECHistograms.h"

class SeagullConfiguration; // Also SeagullConfiguration depends on EECHistogramManager

/*
 * Class for drawing the histograms produced in the dijet analysis
 */
class EECHistogramManager {
 
public:
  
  // Indices for different track histogram categories
  enum enumTrackHistograms{kTrack, kUncorrectedTrack, knTrackCategories};
  
  // Indices for different energy-energy correlator categories
  enum enumEnergyEnergyCorrelators{kEnergyEnergyCorrelator, kEnergyEnergyCorrelatorJetPt, kEnergyEnergyCorrelatorUncorrected, kEnergyEnergyCorrelatorJetPtUncorrected, knEnergyEnergyCorrelatorTypes};
  
  // Indices for different subevent types
  enum enumSubeventTypes{kPythiaPythia, kPythiaHydjet, kHydjetHydjet, knSubeventTypes};
  
  // Dimensions for histogram arrays
  static const int kMaxCentralityBins = 5;       // Maximum allowed number of centrality bins
  static const int kMaxTrackPtBins = 10;         // Maximum allowed number of track pT bins
  static const int knGenJetPtBins = 45;          // Number of generator level jet pT bins for jet pT closures
  static const int knJetEtaBins = 50;            // Number of jet eta bins for jet pT closures
  static const int kMaxJetPtBinsEEC = 12;       // Maximum allowed number of jet pT bins for energy-energy correlators
  static const int kMaxTrackPtBinsEEC = 10;     // Maximum allowed number of track pT bins for energy-energy correlators
  
private:
  
  // Naming for track histograms
  const char* fTrackHistogramNames[knTrackCategories] = {"track","trackUncorrected"}; // Names that different track histograms have in the input file
  const char* fTrackAxisNames[knTrackCategories] = {"Track","Uncorrected track"}; // Names attached to the figure axes
  
  // Naming for jet histograms
  const char* fJetHistogramName = "inclusiveJet";
  const char* fJetAxisName = "Jet";
  
  // Naming for energy-energy correlator histograms
  const char* fEnergyEnergyCorrelatorHistogramNames[knEnergyEnergyCorrelatorTypes] = {"energyEnergyCorrelator", "energyEnergyCorrelatorJetPt", "energyEnergyCorrelatorUncorrected", "energyEnergyCorrelatorJetPtUncorrected"};
  const char* fEnergyEnergyCorrelatorAxisNames[knEnergyEnergyCorrelatorTypes] = {"EEC", "EEC/jet pT", "Uncorrected EEC", "Uncorrected EEC/jet pT"};
  
  // Naming for closure particle
  const char* fClosureParticleName[EECHistograms::knClosureParticleTypes+1] = {"_quark","_gluon",""};
  
  // Naming for subevent types
  const char* fSubeventTypeName[knSubeventTypes] = {"Pythia-Pythia", "Pythia-Hydjet", "Hydjet-Hydjet"};
  
public:
  
  EECHistogramManager();                                    // Default constructor
  EECHistogramManager(TFile *inputFile);                    // Constructor
  EECHistogramManager(TFile *inputFile, EECCard *card);   // Constructor with card
  EECHistogramManager(const EECHistogramManager& in);     // Copy constructor
  ~EECHistogramManager();                                   // Destructor
  
  void LoadHistograms();          // Load the histograms from the inputfile
  void Write(const char* fileName, const char* fileOption);  // Write all the loaded histograms into a file
  void LoadProcessedHistograms(); // Load processed histograms from the inputfile
  
  // Setters for binning information
  void SetCentralityBins(const bool readBinsFromFile, const int nBins, const double *binBorders, bool setIndices = true); // Set up centrality bin indices according to provided bin borders
  void SetTrackPtBins(const bool readBinsFromFile, const int nBins, const double *binBorders, bool setIndices = true);    // Set up track pT bin indices according to provided bin borders
  void SetJetPtBinsEEC(const bool readBinsFromFile, const int nBins, const double *binBorders, bool setIndices = true);   // Set up jet pT bin indices for energy-energy correlator according to provided bin borders
  void SetTrackPtBinsEEC(const bool readBinsFromFile, const int nBins, const double *binBorders, bool setIndices = true); // Set up track pT bin indices for energy-energy correlator according to provided bin borders
  
  // Setters for event information and dijets
  void SetLoadEventInformation(const bool loadOrNot); // Setter for loading event information
  
  // Setters for jets
  void SetLoadJetHistograms(const bool loadOrNot);        // Setter for loading all jet histograms
  
  // Setters for tracks
  void SetLoadTracks(const bool loadOrNot);            // Setter for loading tracks
  void SetLoadTracksUncorrected(const bool loadOrNot); // Setter for loading uncorrected tracks
  void SetLoadAllTracks(const bool loadTracks, const bool loadUncorrected); // Setter for loading all track histograms
  
  // Setters for energy-energy correlator histograms
  void SetLoadEnergyEnergyCorrelators(const bool loadOrNot);                 // Setter for loading energy-energy correlators
  void SetLoadEnergyEnergyCorrelatorsJetPt(const bool loadOrNot);            // Setter for loading jet pT weighted energy-energy correlators
  void SetLoadEnergyEnergyCorrelatorsUncorrected(const bool loadOrNot);      // Setter for loading energy-energy correlators without tracking efficiency
  void SetLoadEnergyEnergyCorrelatorsJetPtUncorrected(const bool loadOrNot); // Setter for loading jet pT weighted energy-energy correlators without tracking efficiency
  void SetLoadAllEnergyEnergyCorrelators(const bool loadRegular, const bool loadJetPt, const bool loadUncorrected, const bool loadJetPtUncorrected); // Setter for loading all energy-energy correlators
  
  // Setter for jet flavor
  void SetJetFlavor(const int iFlavor);  // For Monte Carlo, can select if we are looking for quark or gluon initiated jets
  
  // Setter for loading additional histograms
  void SetLoad2DHistograms(const bool loadOrNot);           // Setter for loading two-dimensional histograms
  void SetLoadJetPtClosureHistograms(const bool loadOrNot); // Setter for loading jet pT closure histograms
  
  // Setters for ranges for different bins
  void SetCentralityBinRange(const int first, const int last); // Setter for centrality bin range
  void SetTrackPtBinRange(const int first, const int last);    // Setter for track pT bin range
  void SetJetPtBinRangeEEC(const int first, const int last);   // Setter for jet pT bin range in energy-energy correlator histograms
  void SetTrackPtBinRangeEEC(const int first, const int last); // Setter for track pT bin range in energy-energy correlator histograms
  
  // Getters for number of bins in histograms
  int GetNCentralityBins() const;  // Getter for the number of centrality bins
  int GetNTrackPtBins() const;     // Getter for the number of track pT bins
  int GetNJetPtBinsEEC() const;    // Getter for the number of jet pT bins in energy-energy correlator histograms
  int GetNTrackPtBinsEEC() const;  // Getter for the number of track pT bins in energy-energy correlator histograms
  double GetCentralityBinBorder(const int iCentrality) const;  // Getter for i:th centrality bin border
  double GetTrackPtBinBorder(const int iTrackPt) const;        // Getter for i:th track pT bin border
  double GetJetPtBinBorderEEC(const int iJetPt) const;         // Getter for i:th jet pT bin border in energy-energy correlator histograms
  double GetTrackPtBinBorderEEC(const int iTrackPt) const;     // Getter for i:th track pT bin border in energy-energy correlator histograms
  
  // Getters for histogram and axis naming
  const char* GetTrackHistogramName(int iTrackType) const; // Getter for track histogram name
  const char* GetTrackAxisName(int iTrackType) const;      // Getter for name suitable for x-axis in a given track histogram
  
  const char* GetJetHistogramName() const; // Getter for the jet histogram name
  const char* GetJetAxisName() const;      // Getter for name suitable for x-axis in a jet histogram
  
  const char* GetEnergyEnergyCorrelatorHistogramName(int iEnergyEnergyCorrelatorType) const; // Getter for energy-energy correlator histogram name
  const char* GetEnergyEnergyCorrelatorAxisName(int iEnergyEnergyCorrelatorType) const;      // Getter for energy-energy correlator axis name
  
  const char* GetSubeventType(const int iSubeventType) const; // Getter for subevent types
  
  TString GetSystem() const;  // Getter for collision system
  
  // Getters for event information histograms
  TH1D* GetHistogramVertexZ() const;            // Getter for z-vertex histogram
  TH1D* GetHistogramVertexZWeighted() const;    // Getter for weighted z-vertex histogram
  TH1D* GetHistogramEvents() const;             // Getter for histogram for number of events surviving different event cuts
  TH1D* GetHistogramTrackCuts() const;          // Getter for histogram for number of tracks surviving different track cuts
  TH1D* GetHistogramCentrality() const;         // Getter for centrality histogram in all events
  TH1D* GetHistogramCentralityWeighted() const; // Getter for weighted centrality histogram in all events
  
  TH1D* GetHistogramMultiplicity(int iCentrality) const;               // Getter for multiplicity from all events
  TH1D* GetHistogramMultiplicityWeighted(int iCentrality) const;       // Getter for efficiency weighted multiplicity from all events
  TH2D* GetHistogramMultiplicityMap() const;                           // Getter for multiplicity vs. centrality map
  TH2D* GetHistogramWeightedMultiplicityMap() const;                   // Getter for efficiency weighted multiplicity vs. centrality map
  
  // Getters for jet histograms
  TH1D* GetHistogramJetPt(int iCentrality) const;     // Jet pT histograms
  TH1D* GetHistogramJetPhi(int iCentrality) const;    // Jet phi histograms
  TH1D* GetHistogramJetEta(int iCentrality) const;    // Jet eta histograms
  TH2D* GetHistogramJetEtaPhi(int iCentrality) const; // 2D eta-phi histogram for jets
  
  // Getters for histograms for tracks
  TH1D* GetHistogramTrackPt(const int iTrackType, const int iCentrality) const;                      // Track pT histograms
  TH1D* GetHistogramTrackPhi(const int iTrackType, const int iCentrality, const int iTrackPt) const;    // Track phi histograms
  TH1D* GetHistogramTrackEta(const int iTrackType, const int iCentrality, const int iTrackPt) const;    // Track eta histograms
  TH2D* GetHistogramTrackEtaPhi(const int iTrackType, const int iCentrality, const int iTrackPt) const; // 2D eta-phi histogram for track
  
  // Getters for energy-energy correlator histograms
  TH1D* GetHistogramEnergyEnergyCorrelator(const int iEnergyEnergyCorrelatorType, const int iCentrality, const int iJetPt, const int iTrackPt, const int iSubevent = knSubeventTypes) const;  // Energy-energy correlator histograms
  
  // Getter for jet pT closure histograms
  TH1D* GetHistogramJetPtClosure(const int iGenPtBin, const int iEtaBin, const int iCentrality, const int iClosureParticle) const; // Jet pT closure
  
  TH1D* GetOneDimensionalHistogram(TString name, int bin1 = 0, int bin2 = 0, int bin3 = 0, int bin4 = 0, int bin5 = 0, int bin6 = 0) const; // Getter for any one-dimensional histogram based on input string
  TH2D* GetTwoDimensionalHistogram(TString name, int bin1 = 0, int bin2 = 0, int bin3 = 0, int bin4 = 0, int bin5 = 0) const; // Getter for any two-dimensional histogram based on input string
  
  
  // Getters for the loaded centrality and track pT bins
  int GetFirstCentralityBin() const;  // Get the first loaded centrality bin
  int GetLastCentralityBin() const;   // Get the last loaded centrality bin
  int GetFirstTrackPtBin() const;     // Get the first loaded track pT bin
  int GetLastTrackPtBin() const;      // Get the last loaded track pT bin
  int GetFirstJetPtBinEEC() const;    // Get the first loaded energy-energy correlator jet pT bin
  int GetLastJetPtBinEEC() const;     // Get the last loaded energy-energy correlator jet pT bin
  int GetFirstTrackPtBinEEC() const;  // Get the first loaded energy-energy correlator track pT bin
  int GetLastTrackPtBinEEC() const;   // Get the last loaded energy-energy correlator track pT bin
  
  // Getters for normalization information
  int GetNEvents() const;                      // Getter for the number of events passing the cuts
  double GetJetPtIntegral(const int iCentrality) const; // Getter for integral over inclusive jet pT in a given centrality
  double GetJetPtIntegral(int iCentrality, const double minPt, const double maxPt) const; // Getter for integral over inclusive jet pT in a given pT range within a given centrality bin
  
  // Getter for the card
  EECCard* GetCard() const;  // Getter for the JCard
  
private:
  
  // Data members
  TFile *fInputFile;                  // File from which the histograms are read
  EECCard *fCard;                   // Card inside the data file for binning, cut collision system etc. information
  TString fSystemAndEnergy;           // Collision system (pp,PbPb,pp MC,PbPb MC,localTest) and energy
  TString fCompactSystemAndEnergy;    // Same a before but without white spaces and dots
  
  // ==============================================
  // ======== Flags for histograms to load ========
  // ==============================================
  
  bool fLoadEventInformation;                              // Load the event information histograms
  bool fLoadJets;                                          // Load the jet histograms
  bool fLoadTracks[knTrackCategories];                     // Load the track histograms
  bool fLoad2DHistograms;                                  // Load also two-dimensional (eta,phi) histograms
  bool fLoadJetPtClosureHistograms;                        // Load the jet pT closure histograms
  bool fLoadEnergyEnergyCorrelatorHistograms[knEnergyEnergyCorrelatorTypes]; // Load the energy-energy correlator histograms
  int  fJetFlavor;                                         // Select the flavor for loaded jets (1 = Quark, 2 = Gluon)
  
  // ==============================================
  // ======== Ranges of histograms to load ========
  // ==============================================
  
  int fFirstLoadedCentralityBin;  // First centrality bin that is loaded
  int fLastLoadedCentralityBin;   // Last centrality bin that is loaded
  int fFirstLoadedTrackPtBin;     // First track pT bin that is loaded
  int fLastLoadedTrackPtBin;      // Last track pT bin that is loaded
  int fFirstLoadedJetPtBinEEC;    // First loaded jet pT bin for the energy-energy correlator histograms
  int fLastLoadedJetPtBinEEC;     // Last loaded jet pT bin for the energy-energy correlator histograms
  int fFirstLoadedTrackPtBinEEC;  // First loaded track pT bin for the energy-energy correlator histograms
  int fLastLoadedTrackPtBinEEC;   // Last loaded track pT bin for the energy-energy correlator histograms
  
  // =============================================
  // ============ Binning information ============
  // =============================================
  int fCentralityBinIndices[kMaxCentralityBins+1];    // Indices for centrality bins in centrality binned histograms
  double fCentralityBinBorders[kMaxCentralityBins+1]; // Centrality bin borders, from which bin indices are obtained
  int fTrackPtBinIndices[kMaxTrackPtBins+1];          // Indices for track pT bins in track eta and phi histograms
  double fTrackPtBinBorders[kMaxTrackPtBins+1];       // Track pT bin borders, from which bin indices are obtained
  int fJetPtIndicesEEC[kMaxJetPtBinsEEC+1];           // Indices for jet pT bins in energy-energy correlator histograms
  double fJetPtBinBordersEEC[kMaxJetPtBinsEEC+1];     // Jet pT bin borders in energy-energy correlator histograms
  int fTrackPtIndicesEEC[kMaxJetPtBinsEEC+1];         // Indices for track pT bins in energy-energy correlator histograms
  double fTrackPtBinBordersEEC[kMaxJetPtBinsEEC+1];   // Track pT bin borders in energy-energy correlator histograms
  int fnCentralityBins;                               // Number of centrality bins in the JCard of the data file
  int fnTrackPtBins;                                  // Number of track pT bins in the JCard of the data file
  int fnJetPtBinsEEC;                                 // Number of jet pT bins for the energy-energy correlator histograms
  int fnTrackPtBinsEEC;                               // Number of track pT bins for the energy-energy correlator histograms
  
  // =============================================
  // ===== Histograms for the dijet analysis =====
  // =============================================
  
  // Event information histograms
  TH1D *fhVertexZ;            // Vertex z position
  TH1D *fhVertexZWeighted;    // Weighted vertex z-position (only meaningfull for MC)
  TH1D *fhEvents;             // Number of events surviving different event cuts
  TH1D *fhTrackCuts;          // Number of tracks surviving different track cuts
  TH1D *fhCentrality;         // Centrality of all events
  TH1D *fhCentralityWeighted; // Weighted centrality distribution in all events (only meaningful for MC)
  TH1D *fhPtHat;              // pT hat for MC events (only meaningful for MC)
  TH1D *fhPtHatWeighted;      // Weighted pT hat distribution (only meaningful for MC)
  
  TH1D *fhMultiplicity[kMaxCentralityBins+1];              // Multiplicity form all events
  TH1D *fhMultiplicityWeighted[kMaxCentralityBins+1];      // Efficiency weighted multiplicity form all events
  TH2D *fhMultiplicityMap;                                 // Multiplicity vs. centrality map
  TH2D *fhMultiplicityMapWeighted;                         // Efficiency weighted multiplicity vs. centrality map
  
  // Histograms for jets
  TH1D *fhJetPt[kMaxCentralityBins];          // Jet pT histograms
  TH1D *fhJetPhi[kMaxCentralityBins];         // Jet phi histograms
  TH1D *fhJetEta[kMaxCentralityBins];         // Jet eta histograms
  TH2D *fhJetEtaPhi[kMaxCentralityBins];      // 2D eta-phi histogram for jets
  
  // Histograms for tracks
  TH1D *fhTrackPt[knTrackCategories][kMaxCentralityBins];                        // Track pT histograms
  TH1D *fhTrackPhi[knTrackCategories][kMaxCentralityBins][kMaxTrackPtBins+1];    // Track phi histograms
  TH1D *fhTrackEta[knTrackCategories][kMaxCentralityBins][kMaxTrackPtBins+1];    // Track eta histograms
  TH2D *fhTrackEtaPhi[knTrackCategories][kMaxCentralityBins][kMaxTrackPtBins+1]; // 2D eta-phi histogram for track
  
  // Histograms for energy-energy correlators
  TH1D *fhEnergyEnergyCorrelator[knEnergyEnergyCorrelatorTypes][kMaxCentralityBins][kMaxJetPtBinsEEC][kMaxTrackPtBinsEEC][knSubeventTypes+1];
  
  // Histograms for jet pT closure
  TH1D *fhJetPtClosure[knGenJetPtBins+1][knJetEtaBins+1][kMaxCentralityBins][EECHistograms::knClosureParticleTypes+1]; // Jet pT closure
  
  // Private methods
  void InitializeFromCard(); // Initialize several member variables from EECCard
  
  // Binning related methods
  void SetBinIndices(const char* histogramName, const int nBins, int *binIndices, const double *binBorders, const int iAxis); // Read the bin indices for given bin borders
  void SetBinBordersAndIndices(const char* histogramName, const int nBins, double *copyBinBorders, int *binIndices, const double *binBorders, const int iAxis, const bool setIndices); // Read the bin indices for given bin borders
  
  // Finders for histograms with different amount of restrictions
  TH2D* FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex); // Extract a 2D histogram using given axis restrictions from THnSparseD
  TH2D* FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2 = 0, int lowBinIndex2 = 0, int highBinIndex2 = 0); // Extract a 2D histogram using given axis restrictions from THnSparseD
  TH1D* FindHistogram(TFile *inputFile, const char *name, int xAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex); // Extract a histogram using given axis restrictions from THnSparseD
  TH1D* FindHistogram(TFile *inputFile, const char *name, int xAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2 = 0, int lowBinIndex2 = 0, int highBinIndex2 = 0); // Extract a histogram using given axis restrictions from THnSparseD
  
  // Loaders for different groups of histograms
  void LoadMultiplicityHistograms(); // Loader for multiplicity histograms
  void LoadJetHistograms(); // Loader for jet histograms
  void LoadTrackHistograms(); // Loader for track histograms
  void LoadEnergyEnergyCorrelatorHistograms(); // Loader for energy-energy correlator histograms
  void LoadJetPtClosureHistograms(); // Loader for jet pT closure histograms
  
  // Generic setter for bin indice and borders
  void SetGenericBins(const bool readBinsFromFile, const char* histogramName, const int iAxis, int nSetBins, double* setBinBorders, int* setBinIndices, const int nBins, const double *binBorders, const char* errorMessage, const int maxBins, const bool setIndices); // Generic bin setter
  
  // Methods for binning
  void BinSanityCheck(const int nBins, int& first, int& last); // Sanity check for given binning
  int BinIndexCheck(const int nBins, const int binIndex) const; // Check that given index is in defined range
  
  // Methods for histogram writing
  void WriteJetHistograms();                    // Write the jet histograms to the file that is currently open
  void WriteTrackHistograms();                  // Write the track histograms to the file that is currently open
  void WriteEnergyEnergyCorrelatorHistograms(); // Write the energy-energy correlator histograms to the file that is currently open
  void WriteClosureHistograms();                // Write the closure histograms to the file that is currently open
  
};

#endif
