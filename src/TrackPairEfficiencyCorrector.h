#ifndef TRACKPAIREFFICIENCYCORRECTOR_H
#define TRACKPAIREFFICIENCYCORRECTOR_H

// Root includes
#include <TH1.h>
#include <TFile.h>

// Own Includes
#include "TrackPairEfficiencyCard.h"

/*
 * TrackPairEfficiencyCorrector class
 *
 * Class whose purpose is to provide weight for reflected cone histograms to take into account jet reconstructions effects to background particles
 */
class TrackPairEfficiencyCorrector {
  
public:
 
  // Maximum dimensions for the arrays
  static const int kMaxCentralityBins = 5;      // Maximum allowed number of centrality bins
  static const int kMaxTrackPtBins = 12;        // Maximum allowed number of track pT bins
  
  TrackPairEfficiencyCorrector();                   // Contructor
  TrackPairEfficiencyCorrector(TString inputFileName);   // Custom constructor
  ~TrackPairEfficiencyCorrector() = default;        // Destructor
  
  double GetTrackPairEfficiencyCorrection(const double deltaR, const double centrality, const double triggerPt, const double associatedPt) const;
  
  void ReadCorrectionTables(); // Initialize the correction tables from file
  
  void SetDisableCorrection(const bool disable); // Setter for disabling the correction
  
private:
  
  // Card with binning information
  TrackPairEfficiencyCard *fCard;
  
  // File from which the histograms are read
  TFile *fInputFile;
  
  // Binning information read to be read from the input file
  double fCentralityBinBorders[kMaxCentralityBins+1] = {0};
  double fTrackPtBinBorders[kMaxTrackPtBins+1] = {0};
  int fnCentralityBins;
  int fnTrackPtBins;
  
  // Flag for disabling the correction
  bool fDisableCorrection;
  
  // Histograms from which the corrections are read
  TH1D *fCorrectionTable[kMaxCentralityBins][kMaxTrackPtBins][kMaxTrackPtBins];
  
  // Methods
  int FindBinIndex(const double* array, const int nBins, const double value) const; // Find a bin index from a given array
  int FindCentralityBin(const double centrality) const; // Find centrality bin index
  int FindTrackPtBin(const double trackPt) const; // Find track pT bin index
  
};

#endif
