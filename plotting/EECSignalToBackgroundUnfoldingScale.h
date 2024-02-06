#ifndef EECSIGNALTOBACKGROUNDUNFOLDINGSCALE_H
#define EECSIGNALTOBACKGROUNDUNFOLDINGSCALE_H

// Root includes
#include <TF1.h>

// Own includes
#include "EECCard.h"

/*
 * EECSignalToBackgroundUnfoldingScale class
 *
 * Class whose purpose is to provide weight for reflected cone histograms to take into account jet reconstructions effects to background particles
 */
class EECSignalToBackgroundUnfoldingScale {
  
public:

  // Dimensions for the arrays are defined by the files used the obtain the weights. They should be copied from the file to here
  static const int kNCentralityBins = 4;      // Number of centrality bins for which background scale is determined
  static const int kNTrackPtBins = 5;         // Number of track pT bins for which the background scale is determined
  static const int kNJetPtBins = 6;         // Number of track pT bins for which the background scale is determined
  
  EECSignalToBackgroundUnfoldingScale();                          // Constructor
  EECSignalToBackgroundUnfoldingScale(EECCard* card);             // Custom constructor
  EECSignalToBackgroundUnfoldingScale(const int iWeightExponent); // Custom constructor
  ~EECSignalToBackgroundUnfoldingScale() = default;      // Destructor

  // Getter for the background scale
  double GetEECSignalToBackgroundUnfoldingScale(const std::pair<double,double> centralityBinBorders, const std::pair<double,double> jetPtBinBorders, const double trackPtLowBorder, const bool isPbPbData) const;
  
private:
  
  // Binning information for the scaling tables
  const double fCentralityBinBorderLow[kNCentralityBins] = {0, 10, 30, 50};
  const double fCentralityBinBorderHigh[kNCentralityBins] = {10, 30, 50, 90};
  const double fCentralityBinBorderLowShifted[kNCentralityBins] = {4, 14, 34, 54};
  const double fCentralityBinBorderHighShifted[kNCentralityBins] = {14, 34, 54, 94};
  const double fJetPtBinBorderLow[kNJetPtBins] = {100, 120, 140, 160, 180, 200};
  const double fJetPtBinBorderHigh[kNJetPtBins] = {120, 140, 160, 180, 200, 220};
  const double fTrackPtBinBorderLow[kNTrackPtBins] = {1, 1.5, 2, 2.5, 3};
  const double fUpshift[kNCentralityBins] = {31.7, 20.6, 12.9, 9.6};

  // Functions from which signal to background ratio scales are read
  TF1* fSignalToBackgroundFunction;

  // Arrays that holds the information about parameters for the functions
  double fMeanJetPt[kNCentralityBins][kNJetPtBins];
  double fFunctionParameters[kNCentralityBins][kNTrackPtBins][3];
  
  // Methods
  void InitializeFunctions(const int iWeightExponent = 1);
  
};

#endif
