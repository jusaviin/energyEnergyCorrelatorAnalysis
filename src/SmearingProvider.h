#ifndef SMEARINGPROVIDER_H
#define SMEARINGPROVIDER_H

// Root includes
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TRandom3.h>

// Own Includes
#include "BinningCard.h"

/*
 * SmearingProvider class
 *
 * Class whose purpose is to provide weight for reflected cone histograms to take into account jet reconstructions effects to background particles
 */
class SmearingProvider {
  
public:
 
  // Maximum dimensions for the arrays
  static const int kMaxCentralityBins = 5;      // Maximum allowed number of centrality bins
  static const int kMaxTrackPtBins = 12;        // Maximum allowed number of track pT bins
  static const int kMaxJetPtBins = 20;          // Maximum allowed number of jet pT bins
  
  SmearingProvider();                   // Contructor
  SmearingProvider(TString inputFileName, TString histogramName, bool useSmearingFactors, bool isPbPbData);   // Custom constructor
  ~SmearingProvider() = default;        // Destructor
  
  double GetSmearedValue(const double valueInNeedOfSmearing, const double centrality, const double jetPt, const double trackPt) const;

  void ReadSmearingTables(); // Initialize the smearing tables from file
  
  void SetDisableSmearing(const bool disable); // Setter for disabling the correction
  
private:
  
  // Card with binning information
  BinningCard* fCard;
  
  // File from which the histograms are read
  TFile* fInputFile;

  // Random number generator for smearing factors
  TRandom3* fRng;
  
  // Binning information read to be read from the input file
  double fCentralityBinBorders[kMaxCentralityBins+1] = {0};
  double fTrackPtBinBorders[kMaxTrackPtBins+1] = {0};
  double fJetPtBinBorders[kMaxJetPtBins+1] = {0};
  int fnCentralityBins;
  int fnTrackPtBins;
  int fnJetPtBins;
  
  // Flags for selecting correct corrections
  bool fDisableSmearing;        // Disable the smearing and just return the input value
  TString fHistogramName;       // Name of the histogram that is read from the file
  bool fUseSmearingFactors;     // True: smearing is done using smearing factors. False: smearing is done using response matrix
  bool fIsPbPbData;             // PbPb and pp have different fluctuation reduction settings. Thus we need a flag for data type
  
  // Histograms for smearing factors and response matrices
  TH1D* fSmearingFactors[kMaxCentralityBins][kMaxJetPtBins][kMaxTrackPtBins];
  TH2D* fResponseMatrix[kMaxCentralityBins][kMaxJetPtBins][kMaxTrackPtBins];

  
  // Methods
  int FindBinIndex(const double* array, const int nBins, const double value) const; // Find a bin index from a given array
  int FindCentralityBin(const double centrality) const; // Find centrality bin index
  int FindTrackPtBin(const double trackPt) const; // Find track pT bin index
  int FindJetPtBin(const double jetPt) const; // Find jet pT bin index
  
};

#endif
