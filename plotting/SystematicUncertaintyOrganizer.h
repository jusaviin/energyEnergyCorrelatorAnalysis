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
  
  enum enumUncertaintySources{kJetEnergyResolution, kJetEnergyScale, kUnfoldingTruth, kTrackSelection, kSingleTrackEfficiency, kTrackPairEfficiency, kBackgroundSubtraction, kAll, knUncertaintySources};
  
  SystematicUncertaintyOrganizer();                                          // Default constructor
  SystematicUncertaintyOrganizer(TFile* inputFile);                          // Constructor
  SystematicUncertaintyOrganizer(const SystematicUncertaintyOrganizer& in);  // Copy constructor
  ~SystematicUncertaintyOrganizer() = default;                               // Default destructor
  
  // Setter for input file
  void ReadInputFile(TFile* inputFile);               // Read the long range systematic uncertainty file
  
  // Getters for the long range systematic uncertainties
  TString GetSystematicUncertaintyName(const int iUncertainty) const;
  TString GetUncertaintyAxisName(const int iUncertainty) const;
  TH1D* GetSystematicUncertainty(const int iCentrality, const int iJetPt, const int iTrackPt, const int iUncertainty = kAll) const;
  
  int GetNUncertaintySources() const;
  
  
private:
    
  TString fSystematicUncertaintyName[knUncertaintySources] = {"jetEnergyResolution", "jetEnergyCorrection", "unfoldingTruth", "trackSelection", "trackEfficiency", "trackPairEfficiency", "backgroundSubtraction", "all"};
  TString fUncertaintyAxisName[knUncertaintySources] = {"JER", "JEC", "unfold truth", "track selection", "track eff", "pair eff", "bg sub" ,"all"};

  // Systematic uncertainty for energy-energy correlators
  TH1D* fhEnergyEnergyCorrelatorUncertainty[kMaxCentralityBins][kMaxJetPtBinsEEC][kMaxTrackPtBinsEEC][knUncertaintySources];  
};

#endif
