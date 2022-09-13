// Reader for jet trees from CMS data
//
//===========================================================
// GeneratorLevelForestReader.h
//
// Author: Jussi Viinikainen
//===========================================================

#ifndef GENERATORLEVELFORESTREADER_H
#define GENERATORLEVELFORESTREADER_H

// Own includes
#include "ForestReader.h"

using namespace std;

class GeneratorLevelForestReader : public ForestReader{
  
private:
  static const Int_t fnMaxJet = 250;        // Maximum number of jets in an event
  
public:
  
  // Constructors and destructors
  GeneratorLevelForestReader();                                                 // Default constructor
  GeneratorLevelForestReader(Int_t dataType, Int_t readMode, Int_t jetType, Int_t jetAxis, Bool_t matchJets, Bool_t readTrackTree = true);    // Custom constructor
  GeneratorLevelForestReader(const GeneratorLevelForestReader& in);             // Copy constructor
  virtual ~GeneratorLevelForestReader();                                        // Destructor
  GeneratorLevelForestReader& operator=(const GeneratorLevelForestReader& obj); // Equal sign operator
  
  // Methods
  void ReadForestFromFile(TFile *inputFile);   // Read the forest from a file
  void ReadForestFromFileList(std::vector<TString> fileList);  // Read the forest from a file list
  void BurnForest();                           // Burn the forest  
  void GetEvent(Int_t nEvent);                 // Get the nEventh event from the file
  
  // Getters for leaves in jet tree
  Float_t GetJetPt(Int_t iJet) const;         // Getter for jet pT
  Float_t GetJetPhi(Int_t iJet) const;        // Getter for jet phi
  Float_t GetJetEta(Int_t iJet) const;        // Getter for jet eta
  Float_t GetJetRawPt(Int_t iJet) const;      // Getter for jet raw pT (not relevant for generator jets, just return value that passes cuts)
  Float_t GetJetMaxTrackPt(Int_t iJet) const; // Getter for maximum track pT inside a jet (not relevant for generator jets, just return value that passes cuts)
  
  // Getters for leaves in the track tree relevant for generator level tracks
  Float_t GetTrackPt(Int_t iTrack) const;          // Getter for track pT
  Float_t GetTrackPhi(Int_t iTrack) const;         // Getter for track phi
  Float_t GetTrackEta(Int_t iTrack) const;         // Getter for track eta
  Int_t GetTrackCharge(Int_t iTrack) const;        // Getter for track charge
  Int_t GetTrackSubevent(Int_t iTrack) const;      // Getter for track subevent index
  Int_t GetTrackMCStatus(Int_t iTrack) const;      // Getter for track MC status (1 = final state particle, other number = not final state particle)
  
  // Getters for leaves in the track tree that are just there to provide values that pass cuts
  Float_t GetTrackPtError(Int_t iTrack) const;               // Getter for track pT error
  Bool_t GetTrackHighPurity(Int_t iTrack) const;             // Getter for the high purity of the track
  Float_t GetTrackVertexDistanceZ(Int_t iTrack) const;       // Getter for track distance from primary vertex in z-direction
  Float_t GetTrackVertexDistanceZError(Int_t iTrack) const;  // Getter for error of track distance from primary vertex in z-direction
  Float_t GetTrackVertexDistanceXY(Int_t iTrack) const;      // Getter for track distance from primary vertex in xy-direction
  Float_t GetTrackVertexDistanceXYError(Int_t iTrack) const; // Getter for error of track distance from primary vertex in xy-direction
  Float_t GetTrackChi2(Int_t iTrack) const;                  // Getter for track chi2 value from reconstruction fit
  Int_t GetNTrackDegreesOfFreedom(Int_t iTrack) const;       // Getter for number of degrees of freedom in reconstruction fit
  Int_t GetNHitsTrackerLayer(Int_t iTrack) const;            // Getter for number of hits in tracker layers
  Int_t GetNHitsTrack(Int_t iTrack) const;                   // Getter for number of hits for the track
  Float_t GetTrackEnergyEcal(Int_t iTrack) const;            // Getter for track energy in ECal
  Float_t GetTrackEnergyHcal(Int_t iTrack) const;            // Getter for track energy in HCal
  
  // Getters for leaves in the particle flow candidate tree
  Int_t GetParticleFlowCandidateId(Int_t iCandidate) const;      // Getter for particle flow candidate ID
  Float_t GetParticleFlowCandidatePt(Int_t iCandidate) const;    // Getter for particle flow candidate pT
  Float_t GetParticleFlowCandidatePhi(Int_t iCandidate) const;   // Getter for particle flow candidate phi
  Float_t GetParticleFlowCandidateEta(Int_t iCandidate) const;   // Getter for particle flow candidate eta
  Int_t GetNParticleFlowCandidates() const;                      // Getter for number of particle flow candidates in an event
  
  // Check if generator level jet has a matching reconstructed jet
  Bool_t HasMatchingJet(Int_t iJet) const;  // Check if generator level jet has a matching reconstructed jet
  Float_t GetMatchedPt(Int_t iJet) const;   // Getter for matched reconstructed jet pT
  Float_t GetMatchedEta(Int_t iJet) const;  // Getter for matched reconstructed jet eta
  Float_t GetMatchedPhi(Int_t iJet) const;  // Getter for matched reconstructed jet phi
  Int_t GetPartonFlavor(Int_t iJet) const;  // Get the flavor of the parton initiating the jet
  
private:
  
  // Methods
  void Initialize();       // Connect the branches to the tree
  Int_t GetMatchingIndex(Int_t iJet) const; // Get the mathing reconstructed jet index for the given generator level jet

  // Trees in the forest
  TTree *fHeavyIonTree;    // Tree for heavy ion event information
  TTree *fJetTree;         // Tree for jet information
  TTree *fHltTree;         // Tree for HLT information
  TTree *fSkimTree;        // Tree for event selection information
  TTree *fTrackTree;       // Tree for tracks  PbPb: anaTrack/trackTree pp: ppTrack/trackTree GenParticles: HiGenParticleAna/hi
  
  // Leaves for jet tree
  Float_t fJetPtArray[fnMaxJet] = {0};         // pT:s of all the jets in an event
  Float_t fJetPhiArray[fnMaxJet] = {0};        // phis of all the jets in an event
  Float_t fJetEtaArray[fnMaxJet] = {0};        // etas of all the jets in an event
  Float_t fJetRefPtArray[fnMaxJet] = {0};      // reference generator level pT for a reconstructed jet
  Int_t fJetRefFlavorArray[fnMaxJet] = {0};    // flavor for initiating parton for the reference gen jet
  Float_t fRecoJetPtArray[fnMaxJet] = {0};     // Array for matched reconstructed jet pT
  Float_t fRecoJetEtaArray[fnMaxJet] = {0};    // Array for matched reconstructed jet pta
  Float_t fRecoJetPhiArray[fnMaxJet] = {0};    // Array for matched reconstructed jet phi
  Float_t fRecoJetRawPtArray[fnMaxJet] = {0};  // Array for matched reconstructed jet raw pT
  
  // Leaves for the track tree
  vector<float> *fTrackPtArray;       // Array for track pT:s
  vector<float> *fTrackPhiArray;      // Array for track phis
  vector<float> *fTrackEtaArray;      // Array for track etas
  vector<int> *fTrackChargeArray;     // Array for track charges
  vector<int> *fTrackSubeventArray;   // Array for track subevent indices (0 = PYTHIA, (>0) = HYDJET)
  vector<int> *fTrackStatusArray;     // Array for Monte Carlo status (1 = final state, others = not final state)
  
};

#endif
