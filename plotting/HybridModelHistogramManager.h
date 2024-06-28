#ifndef HYBRIDMODELHISTOGRAMMANAGER_H
#define HYBRIDMODELHISTOGRAMMANAGER_H

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
class HybridModelHistogramManager {
 
public:

  // Enumeration for different wake configurations
  enum enumTrackHistograms{kNoWake, kPositiveWake, kFullWake, kWakeConfigurations};
  
  // Dimensions for histogram arrays
  static const int kCentralityBins = 4;       // Number of centrality bins for which the predictions are available
  static const int kTrackPtBins = 3;          // Number of track pT bins for which the predictions are available
  static const int kJetPtBins = 4;            // Number of jet pT bins for which the predictions are available
  static const int kEnergyWeights = 2;        // Number of energy weights for which the predictions are available
  
private:
  
  // File name convention to find predictions for different bins
  const char* fCentralityName[kCentralityBins] = {"010", "1030", "3050", "5070"}; // String telling which centrality bin a file corresponds to
  const char* fTrackPtName[kTrackPtBins] = {"iC_0", "iC_1", "iC_2"}; // String telling which track pT cut a file corresponds to
  const char* fJetPtName[kJetPtBins] = {"pTbin_0", "pTbin_1", "pTbin_2", "pTbin_3"}; // String telling which jet pT cut a file corresponds to
  const char* fWakeName[kWakeConfigurations] = {"WantWake_0_IgnoreNeg_1", "WantWake_1_IgnoreNeg_1", "WantWake_1_IgnoreNeg_0"};
  const char* fEnergyWeightName[kEnergyWeights] = {"iN_0", "iN_1"};
  const char* fSubtractionName[kWakeConfigurations] = {"NoElastic", "NoElastic", "NegaSub_Rsub0p2_NoElastic"};
  const char* fWakeLegendName[kWakeConfigurations] = {"Hybrid, no wake", "Hybrid, pos. wake", "Hybrid, full wake"};
  
  // Bin borders for different binning variables
  std::pair<int,int> fCentralityBinBorders[kCentralityBins] = {std::make_pair(0,10), std::make_pair(10,30), std::make_pair(30,50), std::make_pair(50,70)};
  double fTrackPtCuts[kTrackPtBins] = {1, 2, 0};
  std::pair<int,int> fJetPtBinBorders[kJetPtBins] = {std::make_pair(120,140), std::make_pair(140,160), std::make_pair(160,180), std::make_pair(180,200)};
  double fEnergyWeights[kEnergyWeights] = {1, 2};
  
public:
  
  HybridModelHistogramManager();                                      // Default constructor
  HybridModelHistogramManager(TString inputDirectory);                // Constructor with input directory
  HybridModelHistogramManager(const HybridModelHistogramManager& in); // Copy constructor
  ~HybridModelHistogramManager() = default;                           // Destructor
  
  void LoadGraphs(TString inputDirectory);  // Load the graphs from the data files
  
  // Getters for energy-energy correlator graphs
  TGraphErrors* GetEnergyEnergyCorrelatorPp(const int iJetPt, const int iTrackPt, const int iEnergyWeight, const int iWake) const;
  TGraphErrors* GetEnergyEnergyCorrelatorPp(std::pair<int,int> jetPtBin, double trackPtCut, double energyWeight, int iWake) const;
  TGraphErrors* GetEnergyEnergyCorrelatorPbPb(const int iCentrality, const int iJetPt, const int iTrackPt, const int iEnergyWeight, const int iWake) const;
  TGraphErrors* GetEnergyEnergyCorrelatorPbPb(std::pair<int,int> centralityBin, std::pair<int,int> jetPtBin, double trackPtCut, double energyWeight, int iWake) const;
  TGraphErrors* GetEnergyEnergyCorrelatorPbPbToPpRatio(const int iCentrality, const int iJetPt, const int iTrackPt, const int iEnergyWeight, const int iWake) const;
  TGraphErrors* GetEnergyEnergyCorrelatorPbPbToPpRatio(std::pair<int,int> centralityBin, std::pair<int,int> jetPtBin, double trackPtCut, double energyWeight, int iWake) const;

  // Getter for a nice legend name for the wake configuration
  const char* GetWakeName(const int iWake) const;
  
private:

  // Histograms for energy-energy correlators
  TGraphErrors* fEnergyEnergyCorrelatorPp[kJetPtBins][kTrackPtBins][kEnergyWeights][kWakeConfigurations];  // Energy-energy correlators for pp
  TGraphErrors* fEnergyEnergyCorrelatorPbPb[kCentralityBins][kJetPtBins][kTrackPtBins][kEnergyWeights][kWakeConfigurations];  // Energy-energy correlators for pp
  TGraphErrors* fEnergyEnergyCorrelatorPbPbToPpRatio[kCentralityBins][kJetPtBins][kTrackPtBins][kEnergyWeights][kWakeConfigurations];  // Energy-energy correlators for pp
  
  // Find bin indices from bin borders
  int FindCentralityBinIndex(std::pair<int,int> centralityBin) const; // Get a centrality bin index from a given centrality bin borders
  int FindJetPtBinIndex(std::pair<int,int> jetPtBin) const; // Get a jet pT bin index from a given jet pT bin borders
  int FindTrackPtBinIndex(double trackPtBin) const; // Get a track pT bin index from a given track pT cut
  int FindEnergyWeightIndex(double energyWeight) const; // Get an energy weight bin index from a given energy weight value
  
};

#endif
