// Reader for jet trees from CMS data
//
//===========================================================
// ForestReader.h
//
// Author: Jussi Viinikainen
//===========================================================

#ifndef FORESTREADER_H
#define FORESTREADER_H

// C++ includes
#include <iostream>
#include <assert.h>
#include <vector>

// Root includes
#include <TString.h>
#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>
#include <TFile.h>

using namespace std;

class ForestReader{
  
public:
  
  // Possible data types to be read with the reader class
  enum enumDataTypes{kPp, kPbPb, kPpMC, kPbPbMC, kLocalTest, knDataTypes};
  
  // Constructors and destructors
  ForestReader();                                          // Default constructor
  ForestReader(Int_t dataType, Int_t readMode, Int_t jetType, Int_t jetAxis, Bool_t matchJets, Bool_t readTrackTree = true); // Custom constructor
  ForestReader(const ForestReader& in);                    // Copy constructor
  virtual ~ForestReader();                                 // Destructor
  ForestReader& operator=(const ForestReader& obj);        // Equal sign operator
  
  // Methods
  virtual void GetEvent(Int_t nEvent) = 0;                 // Get the nth event in tree
  virtual Int_t GetNEvents() const;                        // Get the number of events
  virtual void ReadForestFromFile(TFile *inputFile) = 0;   // Read the forest from a file
  virtual void ReadForestFromFileList(std::vector<TString> fileList) = 0;   // Read the forest from a file list
  virtual void BurnForest() = 0;                           // Burn the forest
  
  // Getters for leaves in heavy ion tree
  Float_t GetVz() const;              // Getter for vertex z position
  Float_t GetCentrality() const;      // Getter for centrality
  Int_t GetHiBin() const;             // Getter for CMS hiBin
  Float_t GetPtHat() const;           // Getter for pT hat
  Float_t GetEventWeight() const;     // Getter for jet weight in MC
  
  // Getters for leaves in jet tree
  virtual Float_t GetJetPt(Int_t iJet) const = 0;         // Getter for jet pT
  virtual Float_t GetJetPhi(Int_t iJet) const = 0;        // Getter for jet phi
  virtual Float_t GetJetEta(Int_t iJet) const = 0;        // Getter for jet eta
  virtual Int_t GetNJets() const;                         // Getter for number of jets
  virtual Float_t GetJetRawPt(Int_t iJet) const = 0;      // Getter for jet raw pT
  virtual Float_t GetJetMaxTrackPt(Int_t iJet) const = 0; // Getter for maximum track pT inside a jet
  
  // Getters for leaves in HLT tree
  Int_t GetCaloJetFilterBit() const;           // Getter for calorimeter jet filter bit
  
  // Getters for leaves in skim tree
  Int_t GetPrimaryVertexFilterBit() const;           // Getter for primary vertex filter bit
  Int_t GetBeamScrapingFilterBit() const;            // Getter got beam scraping filter bit
  Int_t GetHBHENoiseFilterBit() const;               // Getter for HB/HE noise filter bit
  Int_t GetHfCoincidenceFilterBit() const;           // Getter for hadronic forward coincidence filter bit
  Int_t GetClusterCompatibilityFilterBit() const;    // Getter for cluster compatibility filter bit
  
  // Specific functions for jet closure plots
  virtual Bool_t HasMatchingJet(Int_t iJet) const = 0; // Check if generator level jet has a matching reconstructed jet
  virtual Float_t GetMatchedPt(Int_t iJet) const = 0;  // Getter for matched jet pT (reco for gen and vice versa)
  virtual Float_t GetMatchedEta(Int_t iJet) const = 0; // Getter for matched jet eta (reco for gen and vice versa)
  virtual Float_t GetMatchedPhi(Int_t iJet) const = 0; // Getter for matched jet phi (reco for gen and vice versa)
  virtual Int_t GetPartonFlavor(Int_t iJet) const = 0; // Parton flavor for the parton initiating the jet
  
  // Getters for leaves in the track tree
  virtual Float_t GetTrackPt(Int_t iTrack) const = 0;                    // Getter for track pT
  virtual Float_t GetTrackPtError(Int_t iTrack) const = 0;               // Getter for track pT error
  virtual Float_t GetTrackPhi(Int_t iTrack) const = 0;                   // Getter for track phi
  virtual Float_t GetTrackEta(Int_t iTrack) const = 0;                   // Getter for track eta
  Int_t GetNTracks() const;                                              // Getter for number of tracks
  virtual Bool_t GetTrackHighPurity(Int_t iTrack) const = 0;             // Getter for the high purity of the track
  virtual Float_t GetTrackVertexDistanceZ(Int_t iTrack) const = 0;       // Getter for track distance from primary vertex in z-direction
  virtual Float_t GetTrackVertexDistanceZError(Int_t iTrack) const = 0;  // Getter for error of track distance from primary vertex in z-direction
  virtual Float_t GetTrackVertexDistanceXY(Int_t iTrack) const = 0;      // Getter for track distance from primary vertex in xy-direction
  virtual Float_t GetTrackVertexDistanceXYError(Int_t iTrack) const = 0; // Getter for error of track distance from primary vertex in xy-direction
  virtual Float_t GetTrackChi2(Int_t iTrack) const = 0;                  // Getter for track chi2 value from reconstruction fit
  virtual Int_t GetNTrackDegreesOfFreedom(Int_t iTrack) const = 0;       // Getter for number of degrees of freedom in reconstruction fit
  virtual Int_t GetNHitsTrackerLayer(Int_t iTrack) const = 0;            // Getter for number of hits in tracker layers
  virtual Int_t GetNHitsTrack(Int_t iTrack) const = 0;                   // Getter for number of hits for the track
  virtual Float_t GetTrackEnergyEcal(Int_t iTrack) const = 0;            // Getter for track energy in ECal
  virtual Float_t GetTrackEnergyHcal(Int_t iTrack) const = 0;            // Getter for track energy in HCal
  virtual Int_t GetTrackCharge(Int_t iTrack) const = 0;                  // Getter for track charge (only for generator level tracks)
  virtual Int_t GetTrackSubevent(Int_t iTrack) const = 0;                // Getter for track subevent index (only for generator level tracks)
  virtual Int_t GetTrackMCStatus(Int_t iTrack) const = 0;                // Getter for track MC status (only for generator level tracks)
  
  virtual Int_t GetTrackAlgorithm(Int_t iTrack) const;                   // Getter for track algorithm
  virtual Int_t GetTrackOriginalAlgorithm(Int_t iTrack) const;           // Getter for track original algorithm
  virtual Float_t GetTrackMVA(Int_t iTrack) const;                       // Getter for track MVA
  
  // Setter for data type
  void SetDataType(Int_t dataType); // Setter for data type
  
protected:
  
  // Methods
  virtual void Initialize() = 0;  // Connect the branches to the tree
  
  Int_t fDataType;        // Type of data read with the tree. 0 = pp, 1 = PbPb, 2 = ppMC, 3 = PbPbMC, 4 = LocalTest
  Int_t fUseJetTrigger;   // 0 = Do not use any triggers, 1 = Require jet trigger
  Int_t fJetType;         // Choose the type of jets usedfor analysis. 0 = Calo jets, 1 = PF jets
  Int_t fJetAxis;         // Jet axis used for the jets. 0 = Anti-kT, 1 = Leading particle flow candidate, 2 = WTA
  Bool_t fMatchJets;      // Match generator and reconstructed level jets
  Bool_t fReadTrackTree;  // Read the track trees from the forest
  
  // Branches for heavy ion tree
  TBranch *fHiVzBranch;                   // Branch for vertex z-position
  TBranch *fHiBinBranch;                  // Branch for centrality
  TBranch *fPtHatBranch;                  // Branch for pT hat
  
  // Branches for jet tree
  TBranch *fJetPtBranch;         // Branch for jet pT
  TBranch *fJetPhiBranch;        // Branch for jet phi
  TBranch *fJetEtaBranch;        // Branch for jet eta
  TBranch *fJetRawPtBranch;      // Branch for raw jet pT
  TBranch *fJetMaxTrackPtBranch; // Maximum pT for a track inside a jet
  TBranch *fJetRefPtBranch;      // Branch for reference generator level pT for a reconstructed jet
  TBranch *fJetRefFlavorBranch;  // Branch for flavor for the parton initiating the jet
  TBranch *fJetMatchedPtBranch;  // Branch for the matched jet pT (reco to gen or vice versa)
  TBranch *fJetMatchedEtaBranch; // Branch for the matched jet eta (reco to gen or vice versa)
  TBranch *fJetMatchedPhiBranch; // Branch for the matched jet phi (reco to gen or vice versa)
  TBranch *fEventWeightBranch;     // Branch for jet weight in MC
  
  // Branches for HLT tree
  TBranch *fCaloJetFilterBranch;         // Branch for calo jet filter bit
  
  // Branches for skim tree
  TBranch *fPrimaryVertexBranch;           // Branch for primary vertex filter bit
  TBranch *fBeamScrapingBranch;            // Branch for beam scraping filter bit
  TBranch *fHBHENoiseBranch;               // Branch for HB/HE noise filter bit
  TBranch *fHfCoincidenceBranch;           // Branch for energy recorded in at least 3 HF calorimeter towers
  TBranch *fClusterCompatibilityBranch;    // Branch for cluster compatibility
  
  // Branches for track tree
  TBranch *fTrackPtBranch;                    // Branch for track pT
  TBranch *fTrackPtErrorBranch;               // Branch for track pT error
  TBranch *fTrackPhiBranch;                   // Branch for track phi
  TBranch *fTrackEtaBranch;                   // Branch for track eta
  TBranch *fHighPurityTrackBranch;            // Branch for high purity of the track
  TBranch *fTrackVertexDistanceZBranch;       // Branch for track distance from primary vertex in z-direction
  TBranch *fTrackVertexDistanceZErrorBranch;  // Branch for error for track distance from primary vertex in z-direction
  TBranch *fTrackVertexDistanceXYBranch;      // Branch for track distance from primary vertex in xy-direction
  TBranch *fTrackVertexDistanceXYErrorBranch; // Branch for error for track distance from primary vertex in xy-direction
  TBranch *fTrackChi2Branch;                  // Branch for track chi2 value from reconstruction fit
  TBranch *fnTrackDegreesOfFreedomBranch;     // Branch for number of degrees of freedom in reconstruction fit
  TBranch *fnHitsTrackerLayerBranch;          // Branch for number of hits in tracker layers
  TBranch *fnHitsTrackBranch;                 // Branch for number of hits for the track
  TBranch *fTrackEnergyEcalBranch;            // Branch for track energy in ECal
  TBranch *fTrackEnergyHcalBranch;            // Branch for track energy in HCal
    
  // Leaves for heavy ion tree
  Float_t fVertexZ;    // Vertex z-position
  Int_t fHiBin;        // HiBin = Centrality percentile * 2
  Float_t fPtHat;      // pT hat
  
  // Leaves for jet tree
  Int_t fnJets;          // number of jets in an event
  Int_t fnMatchedJets;   // number of matched jets in an event
  Float_t fEventWeight;    // jet weight in the MC tree
  
  // Leaves for the HLT tree
  Int_t fCaloJetFilterBit;         // Filter bit for calorimeter jets
  
  // Leaves for the skim tree
  Int_t fPrimaryVertexFilterBit;           // Filter bit for primary vertex
  Int_t fBeamScrapingFilterBit;            // Filter bit for beam scraping
  Int_t fHBHENoiseFilterBit;               // Filter bit for HB/HE noise
  Int_t fHfCoincidenceFilterBit;           // Filter bit for energy recorded in at least 3 HF calorimeter towers
  Int_t fClusterCompatibilityFilterBit;    // Filter bit for cluster compatibility
  
  // Leaves for the track tree
  Int_t fnTracks;  // Number of tracks
  
};

#endif
