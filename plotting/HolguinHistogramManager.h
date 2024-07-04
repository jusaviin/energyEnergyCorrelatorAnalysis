#ifndef HOLGUINHISTOGRAMMANAGER_H
#define HOLGUINHISTOGRAMMANAGER_H

// C++ includes
#include <iostream>
#include <fstream>    // File stream for input/output to/from files

// Root includes
#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TGraph.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>

/*
 * Class for drawing the histograms produced in the dijet analysis
 */
class HolguinHistogramManager {
 
public:
  
  // Dimensions for histogram arrays
  static const int kCentralityBins = 1;     // Number of centrality bins for which the predictions are available
  static const int kTrackPtBins = 1;        // Number of track pT bins for which the predictions are available
  static const int kJetPtBins = 4;          // Number of jet pT bins for which the predictions are available
  static const int kEnergyWeights = 1;      // Number of energy weights for which the predictions are available
  static const int kMaxKValues = 7;         // Maximum number of k-values that arrays can hold
  
private:
  
  // File name convention to find predictions for different bins
  const char* fJetPtName[kJetPtBins] = {"pt120-140", "pt140-160", "pt160-180", "pt180-200"}; // String telling which jet pT cut a file corresponds to

  // Vactor for k-value names that are read from the file
  std::vector<double> fKValueVector;
  
  // Bin borders for different binning variables
  std::pair<int,int> fCentralityBinBorders[kCentralityBins] = {std::make_pair(0,10)};
  double fTrackPtCuts[kTrackPtBins] = {1};
  std::pair<int,int> fJetPtBinBorders[kJetPtBins] = {std::make_pair(120,140), std::make_pair(140,160), std::make_pair(160,180), std::make_pair(180,200)};
  double fEnergyWeights[kEnergyWeights] = {1};
  
public:
  
  HolguinHistogramManager();                                      // Default constructor
  HolguinHistogramManager(TString inputDirectory);                // Constructor with input directory
  HolguinHistogramManager(const HolguinHistogramManager& in); // Copy constructor
  ~HolguinHistogramManager() = default;                           // Destructor
  
  void LoadGraphs(TString inputDirectory);  // Load the graphs from the data files
  
  // Getters for energy-energy correlator graphs
  TGraph* GetEnergyEnergyCorrelatorPbPb(const int iCentrality, const int iJetPt, const int iTrackPt, const int iEnergyWeight, const int iKValue) const;
  TGraph* GetEnergyEnergyCorrelatorPbPb(std::pair<int,int> centralityBin, std::pair<double,double> jetPtBin, double trackPtCut, double energyWeight, double valueOfK) const;
  TGraph* GetEnergyEnergyCorrelatorPbPbToPpRatio(const int iCentrality, const int iJetPt, const int iTrackPt, const int iEnergyWeight, const int iKValue) const;
  TGraph* GetEnergyEnergyCorrelatorPbPbToPpRatio(std::pair<int,int> centralityBin, std::pair<double,double> jetPtBin, double trackPtCut, double energyWeight, double valueOfK) const;

  // Find a bin index for k-value
  int FindKValueIndex(double valueOfK) const; // Get an energy weight bin index from a given energy weight value
  
private:

  // Histograms for energy-energy correlators
  TGraph* fEnergyEnergyCorrelatorPbPb[kCentralityBins][kJetPtBins][kTrackPtBins][kEnergyWeights][kMaxKValues];  // Energy-energy correlators for pp
  TGraph* fEnergyEnergyCorrelatorPbPbToPpRatio[kCentralityBins][kJetPtBins][kTrackPtBins][kEnergyWeights][kMaxKValues];  // Energy-energy correlators for pp

  // Number of k values in the histograms loaded from the directory
  int fnKValues;
  bool fKValuesFound;

  // Find graphs from a .dat file formatted in the Holguin way
  std::vector<std::pair<double, TGraph*>> GetGraphsFromDatFile(TString fileName);

  // Getter for the number of k-values in for the histograms
  int GetNKValues() const;
  
  // Find bin indices from bin borders
  int FindCentralityBinIndex(std::pair<int,int> centralityBin) const; // Get a centrality bin index from a given centrality bin borders
  int FindJetPtBinIndex(std::pair<double,double> jetPtBin) const; // Get a jet pT bin index from a given jet pT bin borders
  int FindTrackPtBinIndex(double trackPtBin) const; // Get a track pT bin index from a given track pT cut
  int FindEnergyWeightIndex(double energyWeight) const; // Get an energy weight bin index from a given energy weight value
  
};

#endif
