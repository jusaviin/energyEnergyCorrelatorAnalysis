#ifndef COLBTHISTOGRAMMANAGER_H
#define COLBTHISTOGRAMMANAGER_H

// C++ includes
#include <iostream>
#include <fstream>    // File stream for input/output to/from files

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
class CoLBTHistogramManager {
 
public:
  
  // Dimensions for histogram arrays
  static const int kCentralityBins = 1;     // Number of centrality bins for which the predictions are available
  static const int kTrackPtBins = 1;        // Number of track pT bins for which the predictions are available
  static const int kJetPtBins = 1;          // Number of jet pT bins for which the predictions are available
  static const int kEnergyWeights = 1;      // Number of energy weights for which the predictions are available
  static const int kQValues = 2;            // Number of q-values available in the file directory
  
private:
  
  // File name convention to find predictions for different bins
  const char* fJetPtName[kJetPtBins] = {"120_140"}; // String telling which jet pT cut a file corresponds to
  const char* fTrackPtName[kTrackPtBins] = {"1"}; // String telling which jet pT cut a file corresponds to
  const char* fQValueName[kQValues] = {"q0p5", "q1"};
  
  // Bin borders for different binning variables
  std::pair<int,int> fCentralityBinBorders[kCentralityBins] = {std::make_pair(0,10)};
  double fTrackPtCuts[kTrackPtBins] = {1};
  std::pair<int,int> fJetPtBinBorders[kJetPtBins] = {std::make_pair(120,140)};
  double fEnergyWeights[kEnergyWeights] = {1};
  double fQValues[kQValues] = {0.5, 1};
  
public:
  
  CoLBTHistogramManager();                                      // Default constructor
  CoLBTHistogramManager(TString inputDirectory);                // Constructor with input directory
  CoLBTHistogramManager(const CoLBTHistogramManager& in); // Copy constructor
  ~CoLBTHistogramManager() = default;                           // Destructor
  
  void LoadGraphs(TString inputDirectory);  // Load the graphs from the data files
  
  // Getters for energy-energy correlator graphs
  TGraphErrors* GetEnergyEnergyCorrelatorPbPbToPpRatio(const int iCentrality, const int iJetPt, const int iTrackPt, const int iEnergyWeight, const int iQValue) const;
  TGraphErrors* GetEnergyEnergyCorrelatorPbPbToPpRatio(std::pair<int,int> centralityBin, std::pair<double,double> jetPtBin, double trackPtCut, double energyWeight, double qValue) const;

  // Find a bin index for k-value
  int FindQValueIndex(double qValue) const; // Get an energy weight bin index from a given energy weight value
  
private:

  // Histograms for energy-energy correlators
  TGraphErrors* fEnergyEnergyCorrelatorPbPbToPpRatio[kCentralityBins][kJetPtBins][kTrackPtBins][kEnergyWeights][kQValues];  // Energy-energy correlators for pp
  
  // Find bin indices from bin borders
  int FindCentralityBinIndex(std::pair<int,int> centralityBin) const; // Get a centrality bin index from a given centrality bin borders
  int FindJetPtBinIndex(std::pair<double,double> jetPtBin) const; // Get a jet pT bin index from a given jet pT bin borders
  int FindTrackPtBinIndex(double trackPtBin) const; // Get a track pT bin index from a given track pT cut
  int FindEnergyWeightIndex(double energyWeight) const; // Get an energy weight bin index from a given energy weight value
  
};

#endif
