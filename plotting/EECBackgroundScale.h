#ifndef EECBACKGROUNDSCALE_H
#define EECBACKGROUNDSCALE_H

// Own includes
#include "EECCard.h"

/*
 * EECBackgroundScale class
 *
 * Class whose purpose is to provide weight for reflected cone histograms to take into account jet reconstructions effects to background particles
 */
class EECBackgroundScale {
  
public:
 
  // Dimensions for the arrays are defined by the files used the obtain the weights. They should be copied from the file to here
  static const int kNCentralityBins = 4;      // Number of centrality bins for which background scale is determined
  static const int kNTrackPtBins = 8;         // Number of track pT bins for which the background scale is determined
  static const int kNJetPtBins = 20;          // Number of jet pT bins for which the background scale is determined
  
  EECBackgroundScale();                      // Constructor
  EECBackgroundScale(EECCard* card);         // Custom constructor
  EECBackgroundScale(const bool useGenJets); // Custom constructor
  ~EECBackgroundScale();                     // Destructor

  // Getter for the background scale
  double GetEECBackgroundScale(const std::pair<double,double> centralityBinBorders, const std::pair<double,double> jetPtBinBorders, const double trackPtLowBorder) const;
  
private:
  
  // Binning information for the scaling tables
  const double fCentralityBinBorderLow[kNCentralityBins] = {0, 10, 30, 50};
  const double fCentralityBinBorderHigh[kNCentralityBins] = {10, 30, 50, 90};
  const double fCentralityBinBorderLowShifted[kNCentralityBins] = {4, 14, 34, 54};
  const double fCentralityBinBorderHighShifted[kNCentralityBins] = {14, 34, 54, 94};
  const double fJetPtBinBorderLow[kNJetPtBins] = {60, 70, 80, 90, 100, 110, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 500, 60, 200, 120};
  const double fJetPtBinBorderHigh[kNJetPtBins] = {70, 80, 90, 100, 110, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 500, 5020, 5020, 300, 5020};
  const double fTrackPtBinBorderLow[kNTrackPtBins] = {0.7, 1, 1.5, 2, 2.5, 3, 3.5, 4};

  // Arrays that hold the information about linear fits to the particle density distributions
  double fBackgroundScale[kNCentralityBins][kNJetPtBins][kNTrackPtBins];
  
  // Methods
  void InitializeArrays(const bool useGenJets = false);
  
};

#endif
