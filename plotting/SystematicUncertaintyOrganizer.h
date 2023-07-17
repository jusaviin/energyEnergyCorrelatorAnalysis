#ifndef SYSTEMATICUNCERTAINTYORGANIZER_H
#define SYSTEMATICUNCERTAINTYORGANIZER_H

// C++ includes
#include <iostream>

// Root includes
#include <TFile.h>
#include <TH1.h>
#include <TMath.h>
#include <TString.h>

// Own includes
#include "EECCard.h"

/*
 * Class for drawing the histograms produced in the dijet analysis
 */
class SystematicUncertaintyOrganizer {

private:
  static const int kMaxCentralityBins = 5;      // Maximum allowed number of centrality bins
  static const int kMaxJetPtBinsEEC = 20;       // Maximum allowed number of jet pT bins for energy-energy correlators
  static const int kMaxTrackPtBinsEEC = 20;     // Maximum allowed number of track pT bins for energy-energy correlators
  
public:
  
  enum enumUncertaintySources{kJetEnergyResolution, kJetEnergyScale, kUnfoldingTruth, kTrackSelection, kSingleTrackEfficiency, kTrackPairEfficiency, kBackgroundSubtraction, kCentralityShift, kMonteCarloNonClosure, kAll, knUncertaintySources};
  enum enumGroupFlagExplanation{kCorrelatedInDeltaR, kUncorrelatedInDeltaR, kSkipped, kGroupForAll, knUncertaintyGroups};
  
  SystematicUncertaintyOrganizer();                                          // Default constructor
  SystematicUncertaintyOrganizer(TFile* inputFile);                          // Constructor
  SystematicUncertaintyOrganizer(const SystematicUncertaintyOrganizer& in);  // Copy constructor
  ~SystematicUncertaintyOrganizer() = default;                               // Default destructor
  
  // Setter for input file
  void ReadInputFile(TFile* inputFile);               // Read the long range systematic uncertainty file
  
  // Getters for the long range systematic uncertainties
  TString GetSystematicUncertaintyName(const int iUncertainty) const;
  TString GetSystematicUncertaintyLegendName(const int iUncertainty) const;
  TString GetUncertaintyAxisName(const int iUncertainty) const;
  TH1D* GetSystematicUncertainty(const int iCentrality, const int iJetPt, const int iTrackPt, const int iUncertainty = kAll) const;
  TH1D* GetCorrelatedSystematicUncertainty(const int iCentrality, const int iJetPt, const int iTrackPt) const;
  TH1D* GetUncorrelatedSystematicUncertainty(const int iCentrality, const int iJetPt, const int iTrackPt) const;
  
  int GetNUncertaintySources() const;
  int GetUncertaintyColor(const int iUncertainty) const;
  bool GetSystematicUncertaintyRelevancy(const int iUncertainty, const bool isPbPb) const;
  bool GetSystematicUncertaintyRelevancyForPp(const int iUncertainty) const;
  bool GetSystematicUncertaintyRelevancyForPbPb(const int iUncertainty) const;
  
  
private:
    
  TString fSystematicUncertaintyName[knUncertaintySources] = {"jetEnergyResolution", "jetEnergyCorrection", "unfoldingTruth", "trackSelection", "trackEfficiency", "trackPairEfficiency", "backgroundSubtraction", "centralityShift", "monteCarloNonClosure", "all"};
  TString fSystematicUncertaintyLegendName[knUncertaintySources] = {"Jet energy resolution", "Jet energy scale", "Unfolding prior shape", "Track selection", "Single track efficiency", "Track pair efficiency", "Background subtraction", "Centrality shift", "MC non-closure", "all"};
  TString fUncertaintyAxisName[knUncertaintySources] = {"JER", "JEC", "unfold truth", "track selection", "track eff", "pair eff", "bg sub", "cent shift", "MC non-closure", "all"};

  // Systematic uncertainty for energy-energy correlators
  TH1D* fhEnergyEnergyCorrelatorUncertainty[kMaxCentralityBins][kMaxJetPtBinsEEC][kMaxTrackPtBinsEEC][knUncertaintySources];

  // Combine a predefined group of systematic uncertainty sources
  TH1D* CombineUncertaintySources(const int iCentrality, const int iJetPt, const int iTrackPt, const int iGroup, const char* newName) const;

  // Systematic uncertainty group flag
  int fSystematicsGroupFlag[knUncertaintySources];
  int fUncertaintyColor[knUncertaintySources];      // Standardized color scheme for different uncertainty sources
  bool fIsRelevant[2][knUncertaintySources];        // Define if this uncertainty is included in total uncertainty estimate (0 = pp, 1 = PbPb)
};

#endif
