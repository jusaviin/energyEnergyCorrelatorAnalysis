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
  double GetEECSignalToBackgroundUnfoldingScale(const std::pair<double,double> centralityBinBorders, const std::pair<double,double> jetPtBinBorders, const double trackPtLowBorder, const bool isPbPbData, int iSystematic = 0) const;
  
private:
  
  // Binning information for the scaling tables
  const double fCentralityBinBorderLow[kNCentralityBins] = {0, 10, 30, 50};
  const double fCentralityBinBorderHigh[kNCentralityBins] = {10, 30, 50, 90};
  const double fCentralityBinBorderLowShifted[kNCentralityBins] = {4, 14, 34, 54};
  const double fCentralityBinBorderHighShifted[kNCentralityBins] = {14, 34, 54, 94};
  const double fJetPtBinBorderLow[kNJetPtBins] = {100, 120, 140, 160, 180, 200};
  const double fJetPtBinBorderHigh[kNJetPtBins] = {120, 140, 160, 180, 200, 220};
  const double fTrackPtBinBorderLow[kNTrackPtBins] = {1, 1.5, 2, 2.5, 3};

  // Nominal values
  const double fRelativeUpshift[3][kNCentralityBins] = {{0.076, 0.050, 0.033, 0.029},
                                                        {0.052, 0.042, 0.027, 0.021},
                                                        {0.101, 0.058, 0.039, 0.037}};

  // Alternative values to be only used with MC closures
  //const double fRelativeUpshift[3][kNCentralityBins] = {{0.100, 0.081, 0.064, 0.029},
  //                                                      {0.067, 0.046, 0.030, 0.025},
  //                                                      {0.092, 0.053, 0.036, 0.033}};

  // Functions from which signal to background ratio scales are read
  TF1* fSignalToBackgroundFunction;

  // Arrays that holds the information about parameters for the functions
  double fMeanJetPt[kNCentralityBins][kNJetPtBins];
  double fFunctionParameters[kNCentralityBins][kNTrackPtBins][3];
  
  // Methods
  void InitializeFunctions(const int iWeightExponent = 1);
  
};

#endif
