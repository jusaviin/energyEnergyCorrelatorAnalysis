#ifndef TRACKPAIREFFICIENCYCARD_H
#define TRACKPAIREFFICIENCYCARD_H

// C++ includes
#include <iostream>

// Root includes
#include <TFile.h>
#include <TDirectory.h>
#include <TString.h>
#include <TObjString.h>
#include <TVectorT.h>

/*
 * TrackPairEfficiencyCard class
 *
 * This class reads the ConfigurationCard from the input root file and decodes the
 * necessary information for the analysis.
 */
class TrackPairEfficiencyCard {
  
public:
 
  // Indices for card entries
  enum enumCardEntries{
    kDataType,                  // Data type in the data file (pp, PbPb, pp MC, PbPb MC)
    kUseTrigger,                // Flag to tell if a jet trigger is used in the analysis
    kJetType,                   // 0 = Calorimeter jets, 1 = PF jets
    kJetAxis,                   // 0 = Anti-kt axis, 1 = WTA axis
    kJetEtaCut,                 // Eta cut for jets
    kMinPtCut,                  // Minimum allowed pT for the inclusive jets
    kMaxPtCut,                  // Maximum allowed pT for the inclusive jets
    kCutBadPhiRegion,           // Cut the phi region with bad tracker performance from the analysis
    kMinMaxTrackPtFraction,     // Minimum fraction of jet pT taken by the highest pT track in jet
    kMaxMaxTrackPtFraction,     // Maximum fraction of jet pT taken by the highest pT track in jet
    kTrackEtaCut,               // Eta cut for tracks
    kTriggerEtaCut,             // Stricter eta cut for the trigger particle
    kCutBadPhiRegionTrigger,    // Do not let the trigger particle to be in the phi region with bad tracker performance
    kMinTrackPtCut,             // Minimum accepted track pT
    kMaxTrackPtCut,             // Maximum accepted track pT
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
    kCentralityBinEdges,        // Centrality bin edges
    kTrackPtBinEdges,           // Track pT bin edges
    kTrackPairPtBinEdges,       // Track pT bin edges for track pair histogram
    kPtHatBinEdges,             // pT hat bin edges
    knEntries};                 // Number of entries in the card
  
  // Enumeration for input files used in postprocessing
  enum enumFileNames{kInputFileName,knFileNames};
  
private:
  
  // Names for each entry read from the configuration card
  const char *fCardEntryNames[knEntries] = {"DataType","UseTrigger","JetType","JetAxis","JetEtaCut","MinJetPtCut","MaxJetPtCut","CutBadPhi","MinMaxTrackPtFraction","MaxMaxTrackPtFraction","TrackEtaCut","TriggerEtaCut","CutBadPhiTrigger","MinTrackPtCut","MaxTrackPtCut","MaxTrackPtRelativeError","VertexMaxDistance","CalorimeterSignalLimitPt","HighPtEtFraction","Chi2QualityCut","MinimumTrackHits","SubeventCut","ZVertexCut","LowPtHatCut","HighPtHatCut","CentralityBinEdges","TrackPtBinEdges","TrackPairPtBinEdges","PtHatBinEdges"};
  const char *fFileNameType[knFileNames] = {"input"};
  const char *fFileNameSaveName[knFileNames] = {"InputFile"};
  
  TFile *fInputFile;         // Input file from which all the data is read
  TString fCardDirectory;    // Path to the ConfigurationCard directory
  int fDataType;             // Total number of centrality bins in the analysis
  TString fDataTypeString;   // Total number of eta gaps in the analysis
  TString fAlternativeDataTypeString; // Alternative data type string
  
  void FindDataTypeString(); // Construct a data type string based on information on the card
  void ReadVectors();        // Read the vectors from the file
  
  // Strings for git hash
  TObjString *fGitHash;
  TObjString *fProjectionGitHash;
  
  // Vectors for all the lines inside the card
  TVectorT<float> *fCardEntries[knEntries];   // Array of all the vectors in the card
  TObjString *fFileNames[knFileNames];        // Array for filenames used in postprocessing
  
  // Private methods
  int GetNBins(const int index) const;                            // Get the number of bins for internal index
  double GetLowBinBorder(const int index, const int iBin) const;  // Get the low border of i:th bin from internal index
  double GetHighBinBorder(const int index, const int iBin) const; // Get the high border of i:th bin from internal index
  int GetBinIndex(const int index, const double value) const;     // Get the bin index in the i:th bin from internal index based on given value
   
public:
  
  TrackPairEfficiencyCard(TFile *inFile); // Contructor with input file
  ~TrackPairEfficiencyCard();             // Destructor
  
  TString GetDataType() const;             // Getter for data type string
  TString GetAlternativeDataType() const;  // Getter for alternative data type string
  void Write(TDirectory *file);            // Write the contents of the card to a file
  void Print() const;                      // Print the contents of the card to the console
  
  int GetNCentralityBins() const;   // Get the number of centrality bins
  int GetNTrackPtBins() const;      // Get the number of track pT bins
  int GetNTrackPairPtBins() const;  // Get the number of track pT bins in track pair histograms
  double GetLowBinBorderCentrality(const int iBin) const;    // Get the low border of i:th centrality bin
  double GetHighBinBorderCentrality(const int iBin) const;   // Get the high border of i:th centrality bin
  double GetLowBinBorderTrackPt(const int iBin) const;       // Get the low border of i:th track pT bin
  double GetHighBinBorderTrackPt(const int iBin) const;      // Get the high border of i:th track pT bin
  double GetLowBinBorderTrackPairPt(const int iBin) const;   // Get the low border of i:th track pT bin in track pair histograms
  double GetHighBinBorderTrackPairPt(const int iBin) const;  // Get the high border of i:th track pT bin in track pair histograms
  int GetBinIndexCentrality(const double value) const;       // Get the bin index for a given centrality value
  int GetBinIndexTrackPt(const double value) const;          // Get the bin index for a given track pT value
  int GetBinIndexTrackPairPt(const double value) const;      // Get the bin index for a given track pT value in track pair histograms
  int GetJetType() const;          // Get the jet type index
  double GetJetPtCut() const;      // Get the minimum jet pT cut
  
  void AddOneDimensionalVector(int entryIndex, float entryContent); // Add one dimensional vector to the card
  void AddVector(int entryIndex, int dimension, double *contents); // Add a vector to the card
  void AddFileName(int entryIndex, TString fileName); // Add a file name to the card
  void AddProjectionGitHash(const char* gitHash); // Add a git hash used to project the histograms to the file
  
};

#endif
