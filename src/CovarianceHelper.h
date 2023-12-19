#ifndef COVARIANCEHELPER_H
#define COVARIANCEHELPER_H

// Root includes
#include <TH1.h>
#include <TFile.h>

// Own Includes
#include "BinningCard.h"

/*
 * CoverianceHelper class
 *
 * To calculate covariances, we need to know the per jet average values of EECs. This class is there to help and provide the said values!
 */
class CovarianceHelper {
  
public:
 
  // Maximum dimensions for the arrays
  static const int kMaxCentralityBins = 5;      // Maximum allowed number of centrality bins
  static const int kMaxTrackPtBins = 12;        // Maximum allowed number of track pT bins
  static const int kMaxJetPtBins = 20;          // Maximum allowed number of jet pT bins

  CovarianceHelper();                       // Contructor
  CovarianceHelper(TString inputFileName);  // Custom constructor with input file
  ~CovarianceHelper() = default;            // Destructor
  
  double GetAverageValue(const double deltaR, const double centrality, const double jetPt, const double trackPt) const; // Getter for average EEC value for a given deltaR value
  double GetAverageValue(const int iDeltaR, const double centrality, const double jetPt, const double trackPt) const; // Getter for average EEC value for a given deltaR bin

  void ReadAverageValueTables(); // Initialize the average EEC values from a file
    
private:
  
  // Card with binning information
  BinningCard* fCard;
  
  // File from which the histograms are read
  TFile* fInputFile;
  
  // Binning information read to be read from the input file
  double fCentralityBinBorders[kMaxCentralityBins+1] = {0};
  double fTrackPtBinBorders[kMaxTrackPtBins+1] = {0};
  double fJetPtBinBorders[kMaxJetPtBins+1] = {0};
  int fnCentralityBins;
  int fnTrackPtBins;
  int fnJetPtBins;

  // Other class vairables
  TString fHistogramName = "eecAverageForCovariance"; // Naming for jet histograms
  bool fIsPbPbData; // Flag for PbPb data
  
  // Histograms for smearing factors and response matrices
  TH1D* fAverageEECvalue[kMaxCentralityBins][kMaxJetPtBins][kMaxTrackPtBins];
  
  // Methods
  int FindBinIndex(const double* array, const int nBins, const double value) const; // Find a bin index from a given array
  int FindCentralityBin(const double centrality) const; // Find centrality bin index
  int FindTrackPtBin(const double trackPt) const; // Find track pT bin index
  int FindJetPtBin(const double jetPt) const; // Find jet pT bin index
  
};

#endif
