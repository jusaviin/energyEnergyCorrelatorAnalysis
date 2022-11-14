#ifndef REFLECTEDCONEWEIGHT_H
#define REFLECTEDCONEWEIGHT_H

// Root includes
#include <TF1.h>

/*
 * ReflectedConeWeight class
 *
 * Class whose purpose is to provide weight for reflected cone histograms to take into account jet reconstructions effects to background particles
 */
class ReflectedConeWeight {
  
public:
 
  // Dimensions for the arrays are defined by the files used the obtain the weights. They should be copied from the file to here
  static const int kNCentralityBins = 4;      // Maximum allowed number of centrality bins
  static const int kNTrackPtBins = 8;         // Maximum allowed number of track pT bins
  static const int kNTrackPtBinsMC = 16;         // Maximum allowed number of track pT bins
  static const int kNJetPtBins = 7;           // Number of generator level jet pT bins for jet pT closures
  static const int kMaxParameters = 8;        // Maximum number of parameters in the weight functions
  
  ReflectedConeWeight();   // Contructor
  ~ReflectedConeWeight();  // Destructor
  
  double GetReflectedConeWeight(const bool isData, const double deltaR, const double centrality, const double jetPt, const double trackPt) const;
  double GetReflectedConeWeightData(const double deltaR, const double centrality, const double jetPt, const double trackPt) const;
  double GetReflectedConeWeightMC(const double deltaR, const double centrality, const double jetPt, const double trackPt) const;
  double GetReflectedConeWeightFakeInSignalCone(const double deltaR, const double centrality, const double jetPt, const double trackPt) const;
  
  void SetUseFakeInSignalConeWeight(const bool useFakeInSignalCone);
  void SetDisableWeights(const bool disableWeight);
  
private:
  
  // Binning information in the study used to fit the particle density distributions
  const double kCentralityBinBorders[kNCentralityBins+1] = {0, 10, 30, 50, 90};
  const double kJetPtBinBorders[kNJetPtBins+1] = {120, 140, 160, 180, 200, 300, 500, 5020};
  const double kTrackPtBinBorders[kNTrackPtBins+1] = {0.7, 1, 1.5, 2, 2.5, 3, 3.5, 4, 300};
  const double kTrackPtBinBordersMC[kNTrackPtBinsMC+1] = {0.7, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 10, 12, 16, 20, 300};

  // Arrays that hold the information about linear fits to the particle density distributions
  double fSignalConeMainShape[kNCentralityBins][kNJetPtBins][kNTrackPtBins];
  double fSignalConeKink[kNCentralityBins][kNJetPtBins][kNTrackPtBins];
  double fReflectedConeShape[kNCentralityBins][kNJetPtBins][kNTrackPtBins];
  double fMonteCarloParameter[kMaxParameters][kNCentralityBins][kNJetPtBins][kNTrackPtBinsMC]; // Monte Carlo weight obtained from fit to only Hydjet particles
  double fMonteCarloParameterDataLike[kMaxParameters][kNCentralityBins][kNJetPtBins][kNTrackPtBinsMC]; // Monte Carlo weight obtained similarly as for data

  // Functions that can be initialized using the parameters from the arrays
  TF1* fThreePieceLinear;
  TF1* fTwoPieceLinear;
  TF1* fLinear;
  TF1* fExpoLinear;
  
  // Flags for selected and disabled weights
  bool fUseFakeInSignalConeMC;
  bool fDisableWeight;
  
  // Methods
  void InitializeArrays();
  int FindBinIndex(const double* array, const int nBins, const double value) const; // Find a bin index from a given array
  int FindCentralityBin(const double centrality) const; // Find centrality bin index
  int FindJetPtBin(const double jetPt) const; // Find jet pT bin index
  int FindTrackPtBin(const double trackPt, const double isMC) const; // Find track pT bin index
  
};

#endif
