#ifndef JEWELHISTOGRAMMANAGER_H
#define JEWELHISTOGRAMMANAGER_H

// C++ includes
#include <iostream>
#include <fstream>    // File stream for input/output to/from files
#include <filesystem> // Checking that files exists

// Root includes
#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>

/*
 * Class for drawing the histograms produced in the dijet analysis
 */
class JewelHistogramManager {
 
public:
  
  // Dimensions for histogram arrays
  static const int kCentralityBins = 4;       // Number of centrality bins for which the predictions are available
  static const int kTrackPtBins = 2;          // Number of track pT bins for which the predictions are available
  static const int kJetPtBins = 4;            // Number of jet pT bins for which the predictions are available
  static const int kEnergyWeights = 2;        // Number of energy weights for which the predictions are available
  
private:
  
  // File name convention to find predictions for different bins
  const char* fCentralityName[kCentralityBins] = {"010", "1030", "3050", "5090"}; // String telling which centrality bin a file corresponds to
  const char* fTrackPtName[kTrackPtBins] = {"ptgt1", "ptgt2"}; // String telling which track pT cut a histograms corresponds to
  const char* fJetPtName[kJetPtBins] = {"120140", "140160", "160180", "180200"}; // String telling which jet pT cut a file corresponds to
  const char* fEnergyWeightName[kEnergyWeights] = {"hEEC", "hE2E2C"}; // String telling which energy weight the histogram has
  
  // Bin borders for different binning variables
  std::pair<int,int> fCentralityBinBorders[kCentralityBins] = {std::make_pair(0,10), std::make_pair(10,30), std::make_pair(30,50), std::make_pair(50,90)};
  double fTrackPtCuts[kTrackPtBins] = {1, 2};
  std::pair<int,int> fJetPtBinBorders[kJetPtBins] = {std::make_pair(120,140), std::make_pair(140,160), std::make_pair(160,180), std::make_pair(180,200)};
  double fEnergyWeights[kEnergyWeights] = {1, 2};
  
public:
  
  JewelHistogramManager();                                      // Default constructor
  JewelHistogramManager(TString inputDirectory);                // Constructor with input directory
  JewelHistogramManager(const JewelHistogramManager& in); // Copy constructor
  ~JewelHistogramManager() = default;                           // Destructor
  
  void LoadHistograms(TString inputDirectory);  // Load the graphs from the data files
  
  // Getters for energy-energy correlator graphs
  TH1D* GetEnergyEnergyCorrelatorPp(const int iJetPt, const int iTrackPt, const int iEnergyWeight) const;
  TH1D* GetEnergyEnergyCorrelatorPp(std::pair<int,int> jetPtBin, double trackPtCut, double energyWeight) const;
  TH1D* GetEnergyEnergyCorrelatorPbPb(const int iCentrality, const int iJetPt, const int iTrackPt, const int iEnergyWeight) const;
  TH1D* GetEnergyEnergyCorrelatorPbPb(std::pair<int,int> centralityBin, std::pair<int,int> jetPtBin, double trackPtCut, double energyWeight) const;
  TH1D* GetEnergyEnergyCorrelatorPbPbToPpRatio(const int iCentrality, const int iJetPt, const int iTrackPt, const int iEnergyWeight) const;
  TH1D* GetEnergyEnergyCorrelatorPbPbToPpRatio(std::pair<int,int> centralityBin, std::pair<int,int> jetPtBin, double trackPtCut, double energyWeight) const;

  // Set the normalization region for the analysis
  void SetNormalizationRegion(double lowDeltaR, double highDeltaR);
  
private:

  // Normalization region for the histograms
  double fLowNormalizationDeltaR;
  double fHighNormalizationDeltaR;

  // Histograms for energy-energy correlators
  TH1D* fEnergyEnergyCorrelatorPp[kJetPtBins][kTrackPtBins][kEnergyWeights];  // Energy-energy correlators for pp
  TH1D* fEnergyEnergyCorrelatorPbPb[kCentralityBins][kJetPtBins][kTrackPtBins][kEnergyWeights];  // Energy-energy correlators for PbPb
  TH1D* fEnergyEnergyCorrelatorPbPbToPpRatio[kCentralityBins][kJetPtBins][kTrackPtBins][kEnergyWeights];  // Energy-energy correlators for PbPb to pp ratio

  // Normalize all histograms
  void NormalizeHistograms();
  
  // Find bin indices from bin borders
  int FindCentralityBinIndex(std::pair<int,int> centralityBin) const; // Get a centrality bin index from a given centrality bin borders
  int FindJetPtBinIndex(std::pair<int,int> jetPtBin) const; // Get a jet pT bin index from a given jet pT bin borders
  int FindTrackPtBinIndex(double trackPtBin) const; // Get a track pT bin index from a given track pT cut
  int FindEnergyWeightIndex(double energyWeight) const; // Get an energy weight bin index from a given energy weight value
  
};

#endif
