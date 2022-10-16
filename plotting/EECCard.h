#ifndef EECCARD_H
#define EECCARD_H

// C++ includes
#include <iostream>

// Root includes
#include <TFile.h>
#include <TDirectory.h>
#include <TString.h>
#include <TObjString.h>
#include <TVectorT.h>

/*
 * EECCard class
 *
 * This class reads the ConfigurationCard from the input root file and decodes the
 * necessary information for the analysis.
 */
class EECCard {
  
public:
 
  // Indices for card entries
  enum enumCardEntries{
    kDataType,                  // Data type in the data file (pp, PbPb, pp MC, PbPb MC)
    kMcCorrelationType,         // Monte Carlo correlation type (RecoReco, RecoGen, GenReco, GenGen)
    kMatchJets,                 // 0 = Reco and Gen jets are not matched, 1 = They are matched
    kUseJetTrigger,             // 0 = Do not use any triggers, 1 = Require jet trigger
    kJetType,                   // 0 = Calorimeter jets, 1 = PF jets
    kJetAxis,                   // 0 = Anti-kt axis, 1 = Leading PF candidate axis, 2 = WTA axis
    kJetEtaCut,                 // Eta cut for jets
    kMinPtCut,                  // Minimum allowed pT for the inclusive jets
    kMaxPtCut,                  // Maximum allowed pT for the inclusive jets
    kMinMaxTrackPtFraction,     // Minimum fraction of jet pT taken by the highest pT track in jet
    kMaxMaxTrackPtFraction,     // Maximum fraction of jet pT taken by the highest pT track in jet
    kJetUncertaintyMode,        // 0 = Nominal jet pT, 1 = pT minus uncertainty, 2 = pT plus uncertainty
    kTrackEtaCut,               // Eta cut for tracks
    kMinTrackPtCut,             // Minimum accepted track pT
    kMaxTrackPtRelativeError,   // Maximum relative error allowed for track pT
    kVertexMaxDistance,         // Maximum allowed distance of tracks from reconstructed vertex
    kCalorimeterSignalLimitPt,  // Limit for track pT above which a signal in calorimeters is required
    kHighPtEtFraction,          // Minimum fraction between pT and Et for high pT tracks
    kChi2QualityCut,            // Maximum accepted chi2 for reconstructed tracks
    kMinimumTrackHits,          // Minimum number of hits in tracking for a track
    kSubeventCut,               // 0 = Subevent 0 (Pythia), 1 = Subevent > 0 (Hydjet), 2 = No subevent selection
    kZVertexCut,                // Maximum accepted vz in the event
    kLowPtHatCut,               // Minimum accepted pT hat
    kHighPtHatCut,              // Maximum accepted pT hat
    kMultiplicityMode,          // Select between multiplicity weighting and centrality weighting
    kJetRadius,                 // Radius around the jet axis which is used to select the tracks for energy-energy correlation analysis
    kJetPtBinEdgesEEC,          // Jet pT binning for energy-energy correlation analysis
    kTrackPtBinEdgesEEC,        // Track pT binning for energy-energy correlation analysis
    kCentralityBinEdges,        // Centrality bin edges
    kTrackPtBinEdges,           // Track pT bin edges
    kPtHatBinEdges,             // pT hat bin edges
    kDoReflectedCone,           // 0 = No background estimation, 1 = Estimate background using reflected cone
    knEntries};                 // Number of entries in the card
  
  // Enumeration for input files used in postprocessing
  enum enumFileNames{kInputFileName,knFileNames};
  
private:
  
  // Names for each entry read from the configuration card
  const char *fCardEntryNames[knEntries] = {"DataType","McCorrelationType","MatchJets","UseJetTrigger","JetType","JetAxis","JetEtaCut","MinJetPtCut","MaxJetPtCut","MinMaxTrackPtFraction","MaxMaxTrackPtFraction","JetUncertainty","TrackEtaCut","MinTrackPtCut","MaxTrackPtRelativeError","VertexMaxDistance","CalorimeterSignalLimitPt","HighPtEtFraction","Chi2QualityCut","MinimumTrackHits","SubeventCut","ZVertexCut","LowPtHatCut","HighPtHatCut","MultiplicityMode","JetRadius","JetPtBinEdgesEEC","TrackPtBinEdgesEEC","CentralityBinEdges","TrackPtBinEdges","PtHatBinEdges","DoReflectedCone"};
  const char *fFileNameType[knFileNames] = {"input"};
  const char *fFileNameSaveName[knFileNames] = {"InputFile"};
  
  TFile *fInputFile;         // Input file from which all the data is read
  TString fCardDirectory;    // Path to the ConfigurationCard directory
  int fDataType;             // Total number of centrality bins in the analysis
  int fMonteCarloType;       // Type of Monte Carlo used for jet-track correlations
  TString fDataTypeString;   // Total number of eta gaps in the analysis
  TString fAlternativeDataTypeString; // Alternative data type string
  
  void FindDataTypeString(); // Construct a data type string based on information on the card
  void ReadVectors();        // Read the vectors from the file
  
  // String for git hash
  TObjString *fGitHash;
  
  // Vectors for all the lines inside the card
  TVectorT<float> *fCardEntries[knEntries];   // Array of all the vectors in the card
  TObjString *fFileNames[knFileNames];        // Array for filenames used in postprocessing
  
  // Private methods
  int GetNBins(const int index) const;                            // Get the number of bins for internal index
  double GetLowBinBorder(const int index, const int iBin) const;  // Get the low border of i:th bin from internal index
  double GetHighBinBorder(const int index, const int iBin) const; // Get the high border of i:th bin from internal index
  int GetBinIndex(const int index, const double value) const;     // Get the bin index in the i:th bin from internal index based on given value
   
public:
  
  EECCard(TFile *inFile); // Contructor with input file
  ~EECCard();             // Destructor
  
  TString GetDataType() const;            // Getter for data type string
  TString GetAlternativeDataType() const; // Getter for alternative data type string
  void Write(TDirectory *file);           // Write the contents of the card to a file
  void Print() const;                     // Print the contents of the card to the console
  
  int GetNCentralityBins() const; // Get the number of centrality bins
  int GetNTrackPtBins() const;    // Get the number of track pT bins
  int GetNJetPtBinsEEC() const;   // Get the number of jet pT bins in energy-energy correlator analysis
  int GetNTrackPtBinsEEC() const; // Get the number of track pT bins in energy-energy correlator analysis
  double GetLowBinBorderCentrality(const int iBin) const;  // Get the low border of i:th centrality bin
  double GetLowBinBorderTrackPt(const int iBin) const;     // Get the low border of i:th track pT bin
  double GetLowBinBorderJetPtEEC(const int iBin) const;    // Get the low border of i:th jet pT bin in energy-energy correlator analysis
  double GetLowBinBorderTrackPtEEC(const int iBin) const;  // Get the low border of i:th track pT bin in energy-energy correlator analysis
  double GetHighBinBorderCentrality(const int iBin) const; // Get the high border of i:th centrality bin
  double GetHighBinBorderTrackPt(const int iBin) const;    // Get the high border of i:th track pT bin
  double GetHighBinBorderJetPtEEC(const int iBin) const;   // Get the high border of i:th jet pT bin in energy-energy correlator analysis
  double GetHighBinBorderTrackPtEEC(const int iBin) const; // Get the high border of i:th track pT bin in energy-energy correlator analysis
  int GetBinIndexCentrality(const double value) const;     // Get the bin index for a given centrality value
  int GetBinIndexTrackPt(const double value) const;        // Get the bin index for a given track pT value
  int GetBinIndexJetPtEEC(const double value) const;       // Get the bin index for a given jet pT value in energy-energy correlator analysis
  int GetBinIndexTrackPtEEC(const double value) const;     // Get the bin index for a given track pT value in energy-energy correlator analysis
  int GetSubeventCut() const;      // Get the index for used subevent cut
  int GetJetType() const;          // Get the jet type index
  double GetJetPtCut() const;      // Get the minimum jet pT cut
  bool GetDoReflectedCone() const; // Get the information if reflected cone histograms are filled
  
  void AddOneDimensionalVector(int entryIndex, float entryContent); // Add one dimensional vector to the card
  void AddVector(int entryIndex, int dimension, double *contents); // Add a vector to the card
  void AddFileName(int entryIndex, TString fileName); // Add a file name to the card
  
};

#endif
