#ifndef EECBACKGROUNDSCALE_H
#define EECBACKGROUNDSCALE_H

// Own includes
#include "../src/EECHistograms.h"
#include "EECCard.h"

/*
 * EECBackgroundScale class
 *
 * Class whose purpose is to provide weight for reflected cone histograms to take into account jet reconstructions effects to background particles
 */
class EECBackgroundScale {
  
public:
 
  // Dimensions for the arrays are defined by the files used the obtain the weights. They should be copied from the file to here
  static const int kNCentralityBins = 4;      // Maximum allowed number of centrality bins
  static const int kNTrackPtBins = 8;         // Maximum allowed number of track pT bins
  static const int kNJetPtBins = 7;           // Number of generator level jet pT bins for jet pT closures
  
  EECBackgroundScale();   // Contructor
  ~EECBackgroundScale();  // Destructor

  // Getter for the background scale
  double GetEECBackgroundScale(const int iCentrality, const int iJetPt, const int iTrackPt) const;
  
  // Check that the binning is the same in the file we are correcting, and in the file used to produce the correction
  bool CheckBinBorders(EECCard* card) const;
  
private:
  
  // Binning information in the study used to fit the particle density distributions
  const double kCentralityBinBorders[kNCentralityBins+1] = {0, 10, 30, 50, 90};
  const double kJetPtBinBorders[kNJetPtBins+1] = {120, 140, 160, 180, 200, 300, 500, 5020};
  const double kTrackPtBinBorders[kNTrackPtBins+1] = {0.7, 1, 1.5, 2, 2.5, 3, 3.5, 4, 300};

  // Arrays that hold the information about linear fits to the particle density distributions
  double fBackgroundScale[kNCentralityBins][kNJetPtBins+1][kNTrackPtBins];
  
  // Methods
  void InitializeArrays();
  
};

#endif
