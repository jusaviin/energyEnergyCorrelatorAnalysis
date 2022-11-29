// Class for the main analysis algorithms for leading-subleading jet analysis

#ifndef EECANALYZER_H
#define EECANALYZER_H

// C++ includes
#include <vector>
#include <bitset>
#include <assert.h>   // Standard c++ debugging tool. Terminates the program if expression given evaluates to 0.
#include <tuple>      // For returning several arguments in a transparent manner
#include <fstream>
#include <string>

// Root includes
#include <TString.h>
#include <TRandom3.h>
#include <TMath.h>

// Own includes
#include "ConfigurationCard.h"
#include "EECHistograms.h"
#include "HighForestReader.h"
#include "GeneratorLevelForestReader.h"
#include "JetCorrector.h"
#include "JetUncertainty.h"
#include "trackingEfficiency2018PbPb.h"
#include "trackingEfficiency2017pp.h"
#include "TrackingEfficiencyInterface.h"
#include "ReflectedConeWeight.h"

class EECAnalyzer{
  
private:
  
  enum enumFilledHistograms{kFillEventInformation, kFillJets, kFillTracks, kFillJetConeHistograms, kFillEnergyEnergyCorrelators, kFillEnergyEnergyCorrelatorsUncorrected, kFillJetPtClosure, knFillTypes}; // Histograms to fill
  enum enumSubeventCuts{kSubeventZero,kSubeventNonZero,kSubeventAny,knSubeventCuts}; // Cuts for subevent index
  enum enumMcCorrelationType{kRecoReco,kRecoGen,kGenReco,kGenGen,knMcCorrelationTypes}; // How to correlate jets and tracks in MC
  
public:
  
  // Constructors and destructor
  EECAnalyzer(); // Default constructor
  EECAnalyzer(std::vector<TString> fileNameVector, ConfigurationCard *newCard, bool runLocal); // Custom constructor
  EECAnalyzer(const EECAnalyzer& in); // Copy constructor
  virtual ~EECAnalyzer(); // Destructor
  EECAnalyzer& operator=(const EECAnalyzer& obj); // Equal sign operator
  
  // Methods
  void RunAnalysis();                     // Run the dijet analysis
  EECHistograms* GetHistograms() const;   // Getter for histograms
  
private:
  
  // Private methods
  void CalculateEnergyEnergyCorrelator(const vector<double> selectedTrackPt[2], const vector<double> relativeTrackEta[2], const vector<double> relativeTrackPhi[2], const vector<int> selectedTrackSubevent[2], const double jetPt);  // Calculate energy-energy correlators
  void FillJetPtClosureHistograms(const Int_t jetIndex); // Fill jet pT closure histograms
  void ReadConfigurationFromCard(); // Read all the configuration from the input card
  
  Bool_t PassSubeventCut(const Int_t subeventIndex) const;  // Check if the track passes the set subevent cut
  Bool_t PassTrackCuts(const Int_t iTrack, TH1F *trackCutHistogram, const Bool_t bypassFill = false); // Check if a track passes all the track cuts
  Bool_t PassEventCuts(ForestReader *eventReader, const Bool_t fillHistograms); // Check if the event passes the event cuts
  Double_t GetTrackEfficiencyCorrection(const Int_t iTrack); // Get the track efficiency correction for a given track
  Double_t  GetTrackEfficiencyCorrection(const Float_t trackPt, const Float_t trackEta, const Int_t hiBin); // Get the track efficiency correction for given track and event information
  Double_t GetVzWeight(const Double_t vz) const;  // Get the proper vz weighting depending on analyzed system
  Double_t GetCentralityWeight(const Int_t hiBin) const; // Get the proper centrality weighting depending on analyzed system
  Double_t GetMultiplicityWeight(const Double_t multiplicity) const; // Get the proper multiplicity weight for MC
  Double_t GetJetPtWeight(const Double_t jetPt) const; // Get the proper jet pT weighting for 2017 and 2018 MC
  Double_t GetSmearingFactor(Double_t jetPt, const Double_t centrality); // Getter for jet pT smearing factor
  Int_t GetCentralityBin(const Double_t centrality) const; // Getter for centrality bin
  Double_t GetMultiplicity(); // Get the track multiplicity in the current event
  Double_t GetCentralityFromMultiplicity(const Double_t multiplicity) const; // Get the analysis centrality bin corresponding to the given multiplicity
  Double_t GetDeltaR(const Double_t eta1, const Double_t phi1, const Double_t eta2, const Double_t phi2) const; // Get deltaR between two objects
  Int_t GetSubeventType(const Int_t subevent1, const Int_t subevent2) const; // Get the subevent type from two track subevents
  Int_t GetSubeventIndex(const Int_t subevent) const; // Get the subevent index for a track
  Double_t GetReflectedEta(const Double_t eta) const; // Get jet eta reflected around zero, avoiding overlapping jet cones
  
  // Private data members
  ForestReader *fJetReader;            // Reader for jets in the event
  ForestReader *fTrackReader;          // Readers for tracks in the event
  std::vector<TString> fFileNames;     // Vector for all the files to loop over
  ConfigurationCard *fCard;            // Configuration card for the analysis
  EECHistograms *fHistograms;        // Filled histograms
  TF1 *fVzWeightFunction;              // Weighting function for vz. Needed for MC.
  TF1 *fCentralityWeightFunction;      // Weighting function for centrality. Needed for MC.
  TF1 *fMultiplicityWeightFunction;    // Track multiplicity based weighting function. Can be done instead of centrality weight.
  TF1 *fPtWeightFunction;              // Weighting function for jet pT. Needed for MC.
  TF1 *fSmearingFunction;              // Additional smearing for jets. Needed in systematic uncertainty study.
  TrackingEfficiencyInterface *fTrackEfficiencyCorrector2018;  // Tracking efficiency corrector for 2018 PbPb and 2017 pp data.
  JetCorrector *fJetCorrector2018;     // Class for making jet energy correction for 2018 data
  JetUncertainty *fJetUncertainty2018; // Class for finding uncertainty for jet pT for 2018 data
  ReflectedConeWeight *fReflectedConeWeighter;  // Class for weighting the tracks in the reflected cone
  TRandom3 *fRng;                      // Random number generator
  
  // Analyzed data and forest types
  Int_t fDataType;                   // Analyzed data type
  Int_t fUseJetTrigger;              // 0 = Do not use any triggers, 1 = Require jet trigger
  Int_t fJetType;                    // Type of jets used for analysis. 0 = Calo jets, 1 = PF jets
  Bool_t fMatchJets;                 // Match generator and reconstruction level jets
  Int_t fDebugLevel;                 // Amount of debug messages printed to console
  Int_t fLocalRun;                   // Flag for running locally or on crab
  
  // Weights for filling the MC histograms
  Double_t fVzWeight;                // Weight for vz in MC
  Double_t fCentralityWeight;        // Weight for centrality in MC
  Double_t fPtHatWeight;             // Weight for pT hat in MC
  Double_t fTotalEventWeight;        // Combined weight factor for MC
  
  // Jet and track selection cuts
  Int_t fJetAxis;                      // Used jet axis type. 0 = Anti-kT jet axis, 1 = Axis from leading PF candidate
  Double_t fVzCut;                     // Cut for vertez z-position in an event
  Double_t fMinimumPtHat;              // Minimum accepted pT hat value
  Double_t fMaximumPtHat;              // Maximum accepted pT hat value
  Double_t fJetEtaCut;                 // Eta cut around midrapidity
  Double_t fJetMinimumPtCut;           // Minimum pT cut for jets
  Double_t fJetMaximumPtCut;           // Maximum pT accepted for jets (and tracks)
  Bool_t fCutBadPhiRegion;             // Cut the phi region with bad tracker performance from the analysis
  Double_t fMinimumMaxTrackPtFraction; // Cut for jets consisting only from soft particles
  Double_t fMaximumMaxTrackPtFraction; // Cut for jets consisting only from one high pT
  Int_t fJetUncertaintyMode;           // Use uncertainty for jet pT. 0 = Nominal, 1 = Minus uncertainty, 2 = Plus uncertainty
  Int_t fCombinatorialJetCut;          // Use a cut to remove combinatorial jets in MC
  Double_t fCombinatoriaJetMargin;     // Margin used with the combinatorial jet cut
  Double_t fTrackEtaCut;               // Eta cut around midrapidity
  Double_t fTrackMinPtCut;             // Minimum pT cut
  Double_t fMaxTrackPtRelativeError;   // Maximum relative error for pT
  Double_t fMaxTrackDistanceToVertex;  // Maximum distance to primary vetrex
  Double_t fCalorimeterSignalLimitPt;  // Require signal in calorimeters for track above this pT
  Double_t fHighPtEtFraction;          // For high pT tracks, minimum required Et as a fraction of track pT
  Double_t fChi2QualityCut;            // Quality cut for track reconstruction
  Double_t fMinimumTrackHits;          // Quality cut for track hits
  Int_t fSubeventCut;                  // Cut for the subevent index
  
  // Correlation type for Monte Carlo
  Int_t fMcCorrelationType;            // Correlation type for Monte Carlo. See enumeration enumMcCorrelationType
  
  // Configuration for energy-energy correlators
  Double_t fJetRadius;       // Jet radius parameter
  Bool_t fDoReflectedCone;   // Estimate background from eta-reflected cones
  
  // Which histograms are filled. Do not fill all in order to save memory and not to crash jobs.
  Bool_t fFillEventInformation;                   // Fill event information histograms
  Bool_t fFillJetHistograms;                      // Fill single and dijet histograms
  Bool_t fFillTrackHistograms;                    // Fill inclusive tracks and tracks in dijet events
  Bool_t fFillJetConeHistograms;                  // Fill particle multiplicity and density histograms around the jet cone
  Bool_t fFillEnergyEnergyCorrelators;            // Fill energy-energy correlator histograms
  Bool_t fFillEnergyEnergyCorrelatorsUncorrected; // Fill uncorrected energy-energy correlator histograms
  Bool_t fFillJetPtClosure;                       // Fill jet pT closure histograms
  
  // Weighting mode flag for MC
  Bool_t fMultiplicityMode;       // True: Weight multiplicity to match the data. False: Weight centrality to match the data

};

#endif
