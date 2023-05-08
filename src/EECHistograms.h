// Class for histograms needed in the dijet analysis

#ifndef EECHISTOGRAMS_H
#define EECHISTOGRAMS_H

// Root includes
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>

// Own includes
#include "ConfigurationCard.h"

class EECHistograms{
  
public:
  
  // Enumeration for event types to event histogram and track cuts for track cut histogram
  enum enumEventTypes {kAll, kPrimaryVertex, kHfCoincidence, kClusterCompatibility, kHBHENoise, kBeamScraping, kVzCut, knEventTypes};
  enum enumTriggerBits {kNoTrigger, kOnlyCaloJet60, kOnlyCaloJet80, kOnlyCaloJet100, kCaloJet60And80, kCaloJet60And100, kCaloJet80And100, kCaloJet60And80And100, knTriggerCombinations};
  enum enumTrackCuts {kAllTracks, kMcCharge, kMcSube, kMcStatus, kPtCuts, kEtaCut, kTrackAlgorithm, kHighPurity, kPtError, kVertexDistance, kCaloSignal, kReconstructionQuality, knTrackCuts};
  enum enumClosureParticleType {kQuark,kGluon,knClosureParticleTypes};
  enum enumJetConeTypes{kSignalCone, kReflectedCone, knJetConeTypes};
  enum enumPairingType{kSameJetPair, kSignalReflectedConePair, kReflectedConePair, knPairingTypes}; // Pair tracks from same jet or from reflected jet cone
  enum enumSubeventTypes{kPythia, kHydjet, knSubeventTypes};
  enum enumSubeventCombinations{kPythiaPythia, kPythiaHydjet, kHydjetPythia, kHydjetHydjet, knSubeventCombinations};
    
  // Constructors and destructor
  EECHistograms(); // Default constructor
  EECHistograms(ConfigurationCard* newCard); // Custom constructor
  EECHistograms(const EECHistograms& in); // Copy constructor
  virtual ~EECHistograms(); // Destructor
  EECHistograms& operator=(const EECHistograms& obj); // Equal sign operator
  
  // Methods
  void CreateHistograms();                   // Create all histograms
  void Write() const;                        // Write the histograms to a file that is opened somewhere else
  void Write(TString outputFileName) const;  // Write the histograms to a file
  void SetCard(ConfigurationCard* newCard);  // Set a new configuration card for the histogram class
  
  // Histograms defined public to allow easier access to them. Should not be abused
  // Notation in comments: l = leading jet, s = subleading jet, inc - inclusive jet, uc = uncorrected, ptw = pT weighted
  TH1F* fhVertexZ;                 // Vertex z-position
  TH1F* fhVertexZWeighted;         // Weighted vertex z-position (only meaningfull for MC)
  TH1F* fhEvents;                  // Number of events. For binning see enumEventTypes.
  TH1F* fhTriggers;                // Number of events selected for each trigger. For binning, see enumTriggerBits
  TH1F* fhTriggersAfterSelection;  // Same as before, but after trigger selection. For binning, see enumTriggerBits
  TH1F* fhTrackCuts;               // Number of tracks. For binning see enumTrackCuts.
  TH1F* fhCentrality;              // Centrality information. -0.5 for pp or PYTHIA.
  TH1F* fhCentralityWeighted;      // Weighted centrality distribution (only meaningful for MC)
  TH1F* fhPtHat;                   // pT hat for MC events (only meaningful for MC)
  TH1F* fhPtHatWeighted;           // Weighted pT hat distribution
  THnSparseF* fhMultiplicity;      // Track multiplicity from all events [multiplicity][centrality]
  THnSparseF* fhInclusiveJet;      // Inclusive jet information. Axes: [jet pT][jet phi][jet eta][cent]
  THnSparseF* fhTrack;             // Track histogram. Axes: [pT][phi][eta][cent]
  THnSparseF* fhTrackUncorrected;  // Track histogram for uncorrected tracks. Axes: [uc pT][uc phi][uc eta][cent]
  THnSparseF* fhParticleDensityAroundJet;   // Particle density around the studied jets
  THnSparseF* fhParticlePtDensityAroundJet; // pT weighted particle density around the studied jets
  THnSparseF* fhParticleMultiplicityInJet;                      // Multiplicity of particles within the studied jet cones
  THnSparseF* fhParticleMultiplicityInReflectedCone;            // Multiplicity of particles within the reflected jet cones
  THnSparseF* fhParticleMultiplicityInJetUncorrected;           // Uncorrected multiplicity of particles within the studied jet cones
  THnSparseF* fhParticleMultiplicityInReflectedConeUncorrected; // Uncorrected multiplicity of particles within the reflected jet cones
  THnSparseF* fhMaxPtParticleInJet;                // Histogram for maximum pT particle found in the jet
  THnSparseF* fhEnergyEnergyCorrelator;            // Histogram for energy-energy correlator. Axes: [dR][jetPt][trackPt][centrality]
  THnSparseF* fhEnergyEnergyCorrelatorUncorrected; // Histogram for energy-energy correlator with uncorrected tracks. Axes: [dR][jetPt][trackPt][centrality]
  THnSparseF* fhEnergyEnergyCorrelatorJetPt;            // Histogram for energy-energy correlator with jet pT normalization. Axes: [dR][jetPt][trackPt][centrality]
  THnSparseF* fhEnergyEnergyCorrelatorJetPtUncorrected; // Histogram for energy-energy correlator with jet pT normalization with uncorrected tracks. Axes: [dR][jetPt][trackPt][centrality]
  THnSparseF* fhJetPtClosure; // Jet pT closure histograms. Also information for response matrix. [gen pT][reco pT][centrality][q/g][reco/gen]
  
  // Extra histogram for unfolding study
  THnSparseF* fhUnfoldingMeasured;  // Measured disribution for unfolding
  THnSparseF* fhUnfoldingTruth;     // Truth distribution for unfolding
  THnSparseF* fhUnfoldingResponse;  // Unfolding response matrix

private:
  
  ConfigurationCard* fCard;    // Card for binning info
  const TString kEventTypeStrings[knEventTypes] = {"All", "PrimVertex", "HfCoin2Th4", "ClustCompt", "HBHENoise", "BeamScrape", "v_{z} cut"}; // Strings corresponding to event types
  const TString kTriggerStrings[knTriggerCombinations] = {"No trigger", "60 && !80 && !100", "!60 && 80 && !100", "!60 && !80 && 100", "60 && 80 && !100", "60 && !80 && 100", "!60 && 80 && 100", "60 && 80 && 100"};
  const TString kTrackCutStrings[knTrackCuts] = {"All", "MC Charge", "MC sube", "MC status", "p_{T} cut", "#eta cut", "Track algo", "HighPurity", "p_{T} error", "vertexDist", "caloSignal", "RecoQuality"}; // String corresponding to track cuts
  
};

#endif
