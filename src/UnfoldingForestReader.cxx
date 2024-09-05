// Implementation for UnfoldingForestReader

// Own includes
#include "UnfoldingForestReader.h"

/*
 * Default constructor
 */
UnfoldingForestReader::UnfoldingForestReader() :
  fDataType(0),
  fJetType(0),
  fJetAxis(0),
  fIsMiniAOD(false),
  fHeavyIonTree(0),
  fJetTree(0),
  fTrackTree(0),
  fGenParticleTree(0),
  fHiVzBranch(0),
  fHiBinBranch(0),
  fPtHatBranch(0),
  fEventWeightBranch(0),
  fnJetsBranch(0),
  fJetPtBranch(0),
  fJetPhiBranch(0),
  fJetWTAPhiBranch(0),
  fJetEtaBranch(0),
  fJetWTAEtaBranch(0),
  fJetRawPtBranch(0),
  fJetMaxTrackPtBranch(0),
  fJetRefPtBranch(0),
  fJetRefEtaBranch(0),
  fJetRefPhiBranch(0),
  fnGenJetsBranch(0),
  fGenJetPtBranch(0),
  fGenJetPhiBranch(0),
  fGenJetWTAPhiBranch(0),
  fGenJetEtaBranch(0),
  fGenJetWTAEtaBranch(0),
  fnTracksBranch(0),
  fTrackPtBranch(0),
  fTrackPtErrorBranch(0),
  fTrackPhiBranch(0),
  fTrackEtaBranch(0),
  fHighPurityTrackBranch(0),
  fTrackVertexDistanceZBranch(0),
  fTrackVertexDistanceZErrorBranch(0),
  fTrackVertexDistanceXYBranch(0),
  fTrackVertexDistanceXYErrorBranch(0),
  fTrackChi2Branch(0),
  fnTrackDegreesOfFreedomBranch(0),
  fnHitsTrackerLayerBranch(0),
  fnHitsTrackBranch(0),
  fTrackEnergyEcalBranch(0),
  fTrackEnergyHcalBranch(0),
  fTrackChargeBranch(0),
  fGenParticlePtBranch(0),
  fGenParticlePhiBranch(0),
  fGenParticleEtaBranch(0),
  fGenParticleChargeBranch(0),
  fGenParticleSubeventBranch(0),
  fVertexZ(-100),
  fHiBin(-1),
  fPtHat(0),
  fnJets(0),
  fnGenJets(0),
  fEventWeight(1),
  fJetPtArray(),
  fJetPhiArray(),
  fJetWTAPhiArray(),
  fJetEtaArray(),
  fJetWTAEtaArray(),
  fJetRawPtArray(),
  fJetRefPtArray(),
  fJetRefEtaArray(),
  fJetRefPhiArray(),
  fGenJetPtArray(),
  fGenJetPhiArray(),
  fGenJetWTAPhiArray(),
  fGenJetEtaArray(),
  fGenJetWTAEtaArray(),
  fnTracks(0),
  fTrackPtArray(),
  fTrackPtErrorArray(),
  fTrackPhiArray(),
  fTrackEtaArray(),
  fHighPurityTrackArray(),
  fTrackVertexDistanceZArray(),
  fTrackVertexDistanceZErrorArray(),
  fTrackVertexDistanceXYArray(),
  fTrackVertexDistanceXYErrorArray(),
  fTrackChi2Array(),
  fnTrackDegreesOfFreedomArray(),
  fnHitsTrackerLayerArray(),
  fnHitsTrackArray(),
  fTrackEnergyEcalArray(),
  fTrackEnergyHcalArray(),
  fTrackChargeArray(),
  fTrackPtVector(0),
  fTrackPtErrorVector(0),
  fTrackPhiVector(0),
  fTrackEtaVector(0),
  fHighPurityTrackVector(0),
  fTrackVertexDistanceZVector(0),
  fTrackVertexDistanceZErrorVector(0),
  fTrackVertexDistanceXYVector(0),
  fTrackVertexDistanceXYErrorVector(0),
  fTrackNormalizedChi2Vector(0),
  fnHitsTrackerLayerVector(0),
  fnHitsTrackVector(0),
  fTrackEnergyEcalVector(0),
  fTrackEnergyHcalVector(0),
  fTrackChargeVector(0),
  fnGenParticles(0),
  fGenParticlePtArray(0),
  fGenParticlePhiArray(0),
  fGenParticleEtaArray(0),
  fGenParticleChargeArray(0),
  fGenParticleSubeventArray(0)
{
  // Default constructor
  
  // Initialize fJetMaxTrackPtArray to -1
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetMaxTrackPtArray[i] = -1;
  }
  
}

/*
 * Custom constructor
 *
 *  Arguments:
 *   Int_t dataType: 0 = pp, 1 = PbPb, 2 = pp MC, 3 = PbPb MC
 *   Int_t jetType: 0 = Calo jets, 1 = CSPF jets, 2 = PuPF jets, 3 = Flow jets
 *   Int_t jetAxis: 0 = Anti-kT axis, 1 = WTA axis
 */
UnfoldingForestReader::UnfoldingForestReader(Int_t dataType, Int_t jetType, Int_t jetAxis) :
  fDataType(0),
  fJetType(jetType),
  fJetAxis(jetAxis),
  fIsMiniAOD(false),
  fHeavyIonTree(0),
  fJetTree(0),
  fTrackTree(0),
  fGenParticleTree(0),
  fHiVzBranch(0),
  fHiBinBranch(0),
  fPtHatBranch(0),
  fEventWeightBranch(0),
  fnJetsBranch(0),
  fJetPtBranch(0),
  fJetPhiBranch(0),
  fJetWTAPhiBranch(0),
  fJetEtaBranch(0),
  fJetWTAEtaBranch(0),
  fJetRawPtBranch(0),
  fJetMaxTrackPtBranch(0),
  fJetRefPtBranch(0),
  fJetRefEtaBranch(0),
  fJetRefPhiBranch(0),
  fnGenJetsBranch(0),
  fGenJetPtBranch(0),
  fGenJetPhiBranch(0),
  fGenJetWTAPhiBranch(0),
  fGenJetEtaBranch(0),
  fGenJetWTAEtaBranch(0),
  fnTracksBranch(0),
  fTrackPtBranch(0),
  fTrackPtErrorBranch(0),
  fTrackPhiBranch(0),
  fTrackEtaBranch(0),
  fHighPurityTrackBranch(0),
  fTrackVertexDistanceZBranch(0),
  fTrackVertexDistanceZErrorBranch(0),
  fTrackVertexDistanceXYBranch(0),
  fTrackVertexDistanceXYErrorBranch(0),
  fTrackChi2Branch(0),
  fnTrackDegreesOfFreedomBranch(0),
  fnHitsTrackerLayerBranch(0),
  fnHitsTrackBranch(0),
  fTrackEnergyEcalBranch(0),
  fTrackEnergyHcalBranch(0),
  fTrackChargeBranch(0),
  fGenParticlePtBranch(0),
  fGenParticlePhiBranch(0),
  fGenParticleEtaBranch(0),
  fGenParticleChargeBranch(0),
  fGenParticleSubeventBranch(0),
  fVertexZ(-100),
  fHiBin(-1),
  fPtHat(0),
  fnJets(0),
  fnGenJets(0),
  fEventWeight(1),
  fJetPtArray(),
  fJetPhiArray(),
  fJetWTAPhiArray(),
  fJetEtaArray(),
  fJetWTAEtaArray(),
  fJetRawPtArray(),
  fJetRefPtArray(),
  fJetRefEtaArray(),
  fJetRefPhiArray(),
  fGenJetPtArray(),
  fGenJetPhiArray(),
  fGenJetWTAPhiArray(),
  fGenJetEtaArray(),
  fGenJetWTAEtaArray(),
  fnTracks(0),
  fTrackPtArray(),
  fTrackPtErrorArray(),
  fTrackPhiArray(),
  fTrackEtaArray(),
  fHighPurityTrackArray(),
  fTrackVertexDistanceZArray(),
  fTrackVertexDistanceZErrorArray(),
  fTrackVertexDistanceXYArray(),
  fTrackVertexDistanceXYErrorArray(),
  fTrackChi2Array(),
  fnTrackDegreesOfFreedomArray(),
  fnHitsTrackerLayerArray(),
  fnHitsTrackArray(),
  fTrackEnergyEcalArray(),
  fTrackEnergyHcalArray(),
  fTrackChargeArray(),
  fTrackPtVector(0),
  fTrackPtErrorVector(0),
  fTrackPhiVector(0),
  fTrackEtaVector(0),
  fHighPurityTrackVector(0),
  fTrackVertexDistanceZVector(0),
  fTrackVertexDistanceZErrorVector(0),
  fTrackVertexDistanceXYVector(0),
  fTrackVertexDistanceXYErrorVector(0),
  fTrackNormalizedChi2Vector(0),
  fnHitsTrackerLayerVector(0),
  fnHitsTrackVector(0),
  fTrackEnergyEcalVector(0),
  fTrackEnergyHcalVector(0),
  fTrackChargeVector(0),
  fnGenParticles(0),
  fGenParticlePtArray(0),
  fGenParticlePhiArray(0),
  fGenParticleEtaArray(0),
  fGenParticleChargeArray(0),
  fGenParticleSubeventArray(0)
{
  // Custom constructor
  
  SetDataType(dataType);
  
  // Initialize fJetMaxTrackPtArray to -1
  for(int i = 0; i < fnMaxJet; i++){
    fJetMaxTrackPtArray[i] = -1;
  }
  
}

/*
 * Copy constructor
 */
UnfoldingForestReader::UnfoldingForestReader(const UnfoldingForestReader& in) :
  fDataType(in.fDataType),
  fJetType(in.fJetType),
  fJetAxis(in.fJetAxis),
  fIsMiniAOD(in.fIsMiniAOD),
  fHeavyIonTree(in.fHeavyIonTree),
  fJetTree(in.fJetTree),
  fTrackTree(in.fTrackTree),
  fGenParticleTree(in.fGenParticleTree),
  fHiVzBranch(in.fHiVzBranch),
  fHiBinBranch(in.fHiBinBranch),
  fPtHatBranch(in.fPtHatBranch),
  fEventWeightBranch(in.fEventWeightBranch),
  fnJetsBranch(in.fnJetsBranch),
  fJetPtBranch(in.fJetPtBranch),
  fJetPhiBranch(in.fJetPhiBranch),
  fJetWTAPhiBranch(in.fJetWTAPhiBranch),
  fJetEtaBranch(in.fJetEtaBranch),
  fJetWTAEtaBranch(in.fJetWTAEtaBranch),
  fJetRawPtBranch(in.fJetRawPtBranch),
  fJetMaxTrackPtBranch(in.fJetMaxTrackPtBranch),
  fJetRefPtBranch(in.fJetRefPtBranch),
  fJetRefEtaBranch(in.fJetRefEtaBranch),
  fJetRefPhiBranch(in.fJetRefPhiBranch),
  fnGenJetsBranch(in.fnGenJetsBranch),
  fGenJetPtBranch(in.fGenJetPtBranch),
  fGenJetPhiBranch(in.fGenJetPhiBranch),
  fGenJetWTAPhiBranch(in.fGenJetWTAPhiBranch),
  fGenJetEtaBranch(in.fGenJetEtaBranch),
  fGenJetWTAEtaBranch(in.fGenJetWTAEtaBranch),
  fnTracksBranch(in.fnTracksBranch),
  fTrackPtBranch(in.fTrackPtBranch),
  fTrackPtErrorBranch(in.fTrackPtErrorBranch),
  fTrackPhiBranch(in.fTrackPhiBranch),
  fTrackEtaBranch(in.fTrackEtaBranch),
  fHighPurityTrackBranch(in.fHighPurityTrackBranch),
  fTrackVertexDistanceZBranch(in.fTrackVertexDistanceZBranch),
  fTrackVertexDistanceZErrorBranch(in.fTrackVertexDistanceZErrorBranch),
  fTrackVertexDistanceXYBranch(in.fTrackVertexDistanceXYBranch),
  fTrackVertexDistanceXYErrorBranch(in.fTrackVertexDistanceXYErrorBranch),
  fTrackChi2Branch(in.fTrackChi2Branch),
  fnTrackDegreesOfFreedomBranch(in.fnTrackDegreesOfFreedomBranch),
  fnHitsTrackerLayerBranch(in.fnHitsTrackerLayerBranch),
  fnHitsTrackBranch(in.fnHitsTrackBranch),
  fTrackEnergyEcalBranch(in.fTrackEnergyEcalBranch),
  fTrackEnergyHcalBranch(in.fTrackEnergyHcalBranch),
  fTrackChargeBranch(in.fTrackChargeBranch),
  fGenParticlePtBranch(in.fGenParticlePtBranch),
  fGenParticlePhiBranch(in.fGenParticlePhiBranch),
  fGenParticleEtaBranch(in.fGenParticleEtaBranch),
  fGenParticleChargeBranch(in.fGenParticleChargeBranch),
  fGenParticleSubeventBranch(in.fGenParticleSubeventBranch),
  fVertexZ(in.fVertexZ),
  fHiBin(in.fHiBin),
  fPtHat(in.fPtHat),
  fnJets(in.fnJets),
  fnGenJets(in.fnGenJets),
  fEventWeight(in.fEventWeight),
  fnTracks(in.fnTracks),
  fTrackPtVector(in.fTrackPtVector),
  fTrackPtErrorVector(in.fTrackPtErrorVector),
  fTrackPhiVector(in.fTrackPhiVector),
  fTrackEtaVector(in.fTrackEtaVector),
  fHighPurityTrackVector(in.fHighPurityTrackVector),
  fTrackVertexDistanceZVector(in.fTrackVertexDistanceZVector),
  fTrackVertexDistanceZErrorVector(in.fTrackVertexDistanceZErrorVector),
  fTrackVertexDistanceXYVector(in.fTrackVertexDistanceXYVector),
  fTrackVertexDistanceXYErrorVector(in.fTrackVertexDistanceXYErrorVector),
  fTrackNormalizedChi2Vector(in.fTrackNormalizedChi2Vector),
  fnHitsTrackerLayerVector(in.fnHitsTrackerLayerVector),
  fnHitsTrackVector(in.fnHitsTrackVector),
  fTrackEnergyEcalVector(in.fTrackEnergyEcalVector),
  fTrackEnergyHcalVector(in.fTrackEnergyHcalVector),
  fTrackChargeVector(in.fTrackChargeVector),
  fnGenParticles(in.fnGenParticles),
  fGenParticlePtArray(in.fGenParticlePtArray),
  fGenParticlePhiArray(in.fGenParticlePhiArray),
  fGenParticleEtaArray(in.fGenParticleEtaArray),
  fGenParticleChargeArray(in.fGenParticleChargeArray),
  fGenParticleSubeventArray(in.fGenParticleSubeventArray)
{
  // Copy constructor
  
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetPtArray[i] = in.fJetPtArray[i];
    fJetPhiArray[i] = in.fJetPhiArray[i];
    fJetWTAPhiArray[i] = in.fJetWTAPhiArray[i];
    fJetEtaArray[i] = in.fJetEtaArray[i];
    fJetWTAEtaArray[i] = in.fJetWTAEtaArray[i];
    fJetRawPtArray[i] = in.fJetRawPtArray[i];
    fJetMaxTrackPtArray[i] = in.fJetMaxTrackPtArray[i];
    fJetRefPtArray[i] = in.fJetRefPtArray[i];
    fJetRefEtaArray[i] = in.fJetRefEtaArray[i];
    fJetRefPhiArray[i] = in.fJetRefPhiArray[i];
    fGenJetPtArray[i] = in.fGenJetPtArray[i];
    fGenJetPhiArray[i] = in.fGenJetPhiArray[i];
    fGenJetWTAPhiArray[i] = in.fGenJetWTAPhiArray[i];
    fGenJetEtaArray[i] = in.fGenJetEtaArray[i];
    fGenJetWTAEtaArray[i] = in.fGenJetWTAEtaArray[i];
  }
  
  // Copy the track arrays
  for(Int_t i = 0; i < fnMaxTrack; i++){
    fTrackPtArray[i] = in.fTrackPtArray[i];
    fTrackPtErrorArray[i] = in.fTrackPtErrorArray[i];
    fTrackPhiArray[i] = in.fTrackPhiArray[i];
    fTrackEtaArray[i] = in.fTrackEtaArray[i];
    fHighPurityTrackArray[i] = in.fHighPurityTrackArray[i];
    fTrackVertexDistanceZArray[i] = in.fTrackVertexDistanceZArray[i];
    fTrackVertexDistanceZErrorArray[i] = in.fTrackVertexDistanceZErrorArray[i];
    fTrackVertexDistanceXYArray[i] = in.fTrackVertexDistanceXYArray[i];
    fTrackVertexDistanceXYErrorArray[i] = in.fTrackVertexDistanceXYErrorArray[i];
    fTrackChi2Array[i] = in.fTrackChi2Array[i];
    fnTrackDegreesOfFreedomArray[i] = in.fnTrackDegreesOfFreedomArray[i];
    fnHitsTrackerLayerArray[i] = in.fnHitsTrackerLayerArray[i];
    fnHitsTrackArray[i] = in.fnHitsTrackArray[i];
    fTrackEnergyEcalArray[i] = in.fTrackEnergyEcalArray[i];
    fTrackEnergyHcalArray[i] = in.fTrackEnergyHcalArray[i];
    fTrackChargeArray[i] = in.fTrackChargeArray[i];
  }
}

/*
 * Assignment operator
 */
UnfoldingForestReader& UnfoldingForestReader::operator=(const UnfoldingForestReader& in){
  // Assignment operator
  
  if (&in==this) return *this;
  
  fDataType = in.fDataType;
  fJetType = in.fJetType;
  fJetAxis = in.fJetAxis;
  fIsMiniAOD = in.fIsMiniAOD;
  fHeavyIonTree = in.fHeavyIonTree;
  fJetTree = in.fJetTree;
  fTrackTree = in.fTrackTree;
  fGenParticleTree = in.fGenParticleTree;
  fHiVzBranch = in.fHiVzBranch;
  fHiBinBranch = in.fHiBinBranch;
  fPtHatBranch = in.fPtHatBranch;
  fEventWeightBranch = in.fEventWeightBranch;
  fnJetsBranch = in.fnJetsBranch;
  fJetPtBranch = in.fJetPtBranch;
  fJetPhiBranch = in.fJetPhiBranch;
  fJetWTAPhiBranch = in.fJetWTAPhiBranch;
  fJetEtaBranch = in.fJetEtaBranch;
  fJetWTAEtaBranch = in.fJetWTAEtaBranch;
  fJetRawPtBranch = in.fJetRawPtBranch;
  fJetMaxTrackPtBranch = in.fJetMaxTrackPtBranch;
  fJetRefPtBranch = in.fJetRefPtBranch;
  fJetRefEtaBranch = in.fJetRefEtaBranch;
  fJetRefPhiBranch = in.fJetRefPhiBranch;
  fnGenJetsBranch = in.fnGenJetsBranch;
  fGenJetPtBranch = in.fGenJetPtBranch;
  fGenJetPhiBranch = in.fGenJetPhiBranch;
  fGenJetWTAPhiBranch = in.fGenJetWTAPhiBranch;
  fGenJetEtaBranch = in.fGenJetEtaBranch;
  fGenJetWTAEtaBranch = in.fGenJetWTAEtaBranch;
  fnTracksBranch = in.fnTracksBranch;
  fTrackPtBranch = in.fTrackPtBranch;
  fTrackPtErrorBranch = in.fTrackPtErrorBranch;
  fTrackPhiBranch = in.fTrackPhiBranch;
  fTrackEtaBranch = in.fTrackEtaBranch;
  fHighPurityTrackBranch = in.fHighPurityTrackBranch;
  fTrackVertexDistanceZBranch = in.fTrackVertexDistanceZBranch;
  fTrackVertexDistanceZErrorBranch = in.fTrackVertexDistanceZErrorBranch;
  fTrackVertexDistanceXYBranch = in.fTrackVertexDistanceXYBranch;
  fTrackVertexDistanceXYErrorBranch = in.fTrackVertexDistanceXYErrorBranch;
  fTrackChi2Branch = in.fTrackChi2Branch;
  fnTrackDegreesOfFreedomBranch = in.fnTrackDegreesOfFreedomBranch;
  fnHitsTrackerLayerBranch = in.fnHitsTrackerLayerBranch;
  fnHitsTrackBranch = in.fnHitsTrackBranch;
  fTrackEnergyEcalBranch = in.fTrackEnergyEcalBranch;
  fTrackEnergyHcalBranch = in.fTrackEnergyHcalBranch;
  fTrackChargeBranch = in.fTrackChargeBranch;
  fGenParticlePtBranch = in.fGenParticlePtBranch;
  fGenParticlePhiBranch = in.fGenParticlePhiBranch;
  fGenParticleEtaBranch = in.fGenParticleEtaBranch;
  fGenParticleChargeBranch = in.fGenParticleChargeBranch;
  fGenParticleSubeventBranch = in.fGenParticleSubeventBranch;
  fVertexZ = in.fVertexZ;
  fHiBin = in.fHiBin;
  fPtHat = in.fPtHat;
  fnJets = in.fnJets;
  fnGenJets = in.fnGenJets;
  fEventWeight = in.fEventWeight;
  fnTracks = in.fnTracks;
  
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetPtArray[i] = in.fJetPtArray[i];
    fJetPhiArray[i] = in.fJetPhiArray[i];
    fJetWTAPhiArray[i] = in.fJetWTAPhiArray[i];
    fJetEtaArray[i] = in.fJetEtaArray[i];
    fJetWTAEtaArray[i] = in.fJetWTAEtaArray[i];
    fJetRawPtArray[i] = in.fJetRawPtArray[i];
    fJetMaxTrackPtArray[i] = in.fJetMaxTrackPtArray[i];
    fJetRefPtArray[i] = in.fJetRefPtArray[i];
    fJetRefEtaArray[i] = in.fJetRefEtaArray[i];
    fJetRefPhiArray[i] = in.fJetRefPhiArray[i];
    fGenJetPtArray[i] = in.fGenJetPtArray[i];
    fGenJetPhiArray[i] = in.fGenJetPhiArray[i];
    fGenJetWTAPhiArray[i] = in.fGenJetWTAPhiArray[i];
    fGenJetEtaArray[i] = in.fGenJetEtaArray[i];
    fGenJetWTAEtaArray[i] = in.fGenJetWTAEtaArray[i];
  }
  
  // Copy the track arrays
  for(Int_t i = 0; i < fnMaxTrack; i++){
    fTrackPtArray[i] = in.fTrackPtArray[i];
    fTrackPtErrorArray[i] = in.fTrackPtErrorArray[i];
    fTrackPhiArray[i] = in.fTrackPhiArray[i];
    fTrackEtaArray[i] = in.fTrackEtaArray[i];
    fHighPurityTrackArray[i] = in.fHighPurityTrackArray[i];
    fTrackVertexDistanceZArray[i] = in.fTrackVertexDistanceZArray[i];
    fTrackVertexDistanceZErrorArray[i] = in.fTrackVertexDistanceZErrorArray[i];
    fTrackVertexDistanceXYArray[i] = in.fTrackVertexDistanceXYArray[i];
    fTrackVertexDistanceXYErrorArray[i] = in.fTrackVertexDistanceXYErrorArray[i];
    fTrackChi2Array[i] = in.fTrackChi2Array[i];
    fnTrackDegreesOfFreedomArray[i] = in.fnTrackDegreesOfFreedomArray[i];
    fnHitsTrackerLayerArray[i] = in.fnHitsTrackerLayerArray[i];
    fnHitsTrackArray[i] = in.fnHitsTrackArray[i];
    fTrackEnergyEcalArray[i] = in.fTrackEnergyEcalArray[i];
    fTrackEnergyHcalArray[i] = in.fTrackEnergyHcalArray[i];
    fTrackChargeArray[i] = in.fTrackChargeArray[i];
  }
  
  // Copy the track vectors
  fTrackPtVector = in.fTrackPtVector;
  fTrackPtErrorVector = in.fTrackPtErrorVector;
  fTrackPhiVector = in.fTrackPhiVector;
  fTrackEtaVector = in.fTrackEtaVector;
  fHighPurityTrackVector = in.fHighPurityTrackVector;
  fTrackVertexDistanceZVector = in.fTrackVertexDistanceZVector;
  fTrackVertexDistanceZErrorVector = in.fTrackVertexDistanceZErrorVector;
  fTrackVertexDistanceXYVector = in.fTrackVertexDistanceXYVector;
  fTrackVertexDistanceXYErrorVector = in.fTrackVertexDistanceXYErrorVector;
  fTrackNormalizedChi2Vector = in.fTrackNormalizedChi2Vector;
  fnHitsTrackerLayerVector = in.fnHitsTrackerLayerVector;
  fnHitsTrackVector = in.fnHitsTrackVector;
  fTrackEnergyEcalVector = in.fTrackEnergyEcalVector;
  fTrackEnergyHcalVector = in.fTrackEnergyHcalVector;
  fTrackChargeVector = in.fTrackChargeVector;
  
  // Copy the generator level particle vectors
  fnGenParticles = in.fnGenParticles;
  fGenParticlePtArray = in.fGenParticlePtArray;
  fGenParticlePhiArray = in.fGenParticlePhiArray;
  fGenParticleEtaArray = in.fGenParticleEtaArray;
  fGenParticleChargeArray = in.fGenParticleChargeArray;
  fGenParticleSubeventArray = in.fGenParticleSubeventArray;
  
  return *this;
}

/*
 * Destructor
 */
UnfoldingForestReader::~UnfoldingForestReader(){
  // destructor
}

/*
 * Initialization, meaning that the branches are connected to the tree
 */
void UnfoldingForestReader::Initialize(){
  
  // Connect the branches of the heavy ion tree
  fHeavyIonTree->SetBranchStatus("*",0);
  fHeavyIonTree->SetBranchStatus("vz",1);
  fHeavyIonTree->SetBranchAddress("vz",&fVertexZ,&fHiVzBranch);
  fHeavyIonTree->SetBranchStatus("hiBin",1);
  fHeavyIonTree->SetBranchAddress("hiBin",&fHiBin,&fHiBinBranch);
  if(fDataType == kPpMC || fDataType == kPbPbMC){
    fHeavyIonTree->SetBranchStatus("pthat",1);
    fHeavyIonTree->SetBranchAddress("pthat",&fPtHat,&fPtHatBranch); // pT hat only for MC
    fHeavyIonTree->SetBranchStatus("weight",1);
    fHeavyIonTree->SetBranchAddress("weight",&fEventWeight,&fEventWeightBranch); // event weight only for MC
  } else {
    fPtHat = 0; // We do not have pT hat information for real data
    fEventWeight = 1;
  }
  
  // Connect the branches to the jet tree  
  fJetTree->SetBranchStatus("*",0);
  fJetTree->SetBranchStatus("jtpt",1);
  fJetTree->SetBranchAddress("jtpt",&fJetPtArray,&fJetPtBranch);
  
  // Load jet phi with E-scheme and WTA axes
  fJetTree->SetBranchStatus("jtphi",1);
  fJetTree->SetBranchAddress("jtphi",&fJetPhiArray,&fJetPhiBranch);
  fJetTree->SetBranchStatus("WTAphi",1);
  fJetTree->SetBranchAddress("WTAphi",&fJetWTAPhiArray,&fJetWTAPhiBranch);
  
  // Load jet eta with E-scheme and WTA axes
  fJetTree->SetBranchStatus("jteta",1);
  fJetTree->SetBranchAddress("jteta",&fJetEtaArray,&fJetEtaBranch);
  fJetTree->SetBranchStatus("WTAeta",1);
  fJetTree->SetBranchAddress("WTAeta",&fJetWTAEtaArray,&fJetWTAEtaBranch);
  
  fJetTree->SetBranchStatus("nref",1);
  fJetTree->SetBranchAddress("nref",&fnJets,&fnJetsBranch);
  fJetTree->SetBranchStatus("rawpt",1);
  fJetTree->SetBranchAddress("rawpt",&fJetRawPtArray,&fJetRawPtBranch);
  fJetTree->SetBranchStatus("trackMax",1);
  fJetTree->SetBranchAddress("trackMax",&fJetMaxTrackPtArray,&fJetMaxTrackPtBranch);

  // If we are looking at Monte Carlo, connect the reference jet and gen jet arrays
  if(fDataType > kPbPb){

    fJetTree->SetBranchStatus("refpt",1);
    fJetTree->SetBranchAddress("refpt",&fJetRefPtArray,&fJetRefPtBranch);
    fJetTree->SetBranchStatus("refeta",1);
    fJetTree->SetBranchAddress("refeta",&fJetRefEtaArray,&fJetRefEtaBranch);
    fJetTree->SetBranchStatus("refphi",1);
    fJetTree->SetBranchAddress("refphi",&fJetRefPhiArray,&fJetRefPhiBranch);

    fJetTree->SetBranchStatus("genpt",1);
    fJetTree->SetBranchAddress("genpt",&fGenJetPtArray,&fGenJetPtBranch);
    
    // Load jet phi with E-scheme and WTA axes
    fJetTree->SetBranchStatus("genphi",1);
    fJetTree->SetBranchAddress("genphi",&fGenJetPhiArray,&fGenJetPhiBranch);
    fJetTree->SetBranchStatus("WTAgenphi",1);
    fJetTree->SetBranchAddress("WTAgenphi",&fGenJetWTAPhiArray,&fGenJetWTAPhiBranch);
    
    // Load jet eta with E-scheme and WTA axes
    fJetTree->SetBranchStatus("geneta",1);
    fJetTree->SetBranchAddress("geneta",&fGenJetEtaArray,&fGenJetEtaBranch);
    fJetTree->SetBranchStatus("WTAgeneta",1);
    fJetTree->SetBranchAddress("WTAgeneta",&fGenJetWTAEtaArray,&fGenJetWTAEtaBranch);
    
    fJetTree->SetBranchStatus("ngen",1);
    fJetTree->SetBranchAddress("ngen",&fnGenJets,&fnGenJetsBranch);
  }
  
  // Connect the branches to the track tree
  
  fTrackTree->SetBranchStatus("*",0);
  
  // We need to read the forest to vectors for MiniAODs and to arrays for AODs
  if(fIsMiniAOD){
    
    fTrackTree->SetBranchStatus("trkPt",1);
    fTrackTree->SetBranchAddress("trkPt",&fTrackPtVector,&fTrackPtBranch);
    fTrackTree->SetBranchStatus("trkPtError",1);
    fTrackTree->SetBranchAddress("trkPtError",&fTrackPtErrorVector,&fTrackPtErrorBranch);
    fTrackTree->SetBranchStatus("trkPhi",1);
    fTrackTree->SetBranchAddress("trkPhi",&fTrackPhiVector,&fTrackPhiBranch);
    fTrackTree->SetBranchStatus("trkEta",1);
    fTrackTree->SetBranchAddress("trkEta",&fTrackEtaVector,&fTrackEtaBranch);
    fTrackTree->SetBranchStatus("nTrk",1);
    fTrackTree->SetBranchAddress("nTrk",&fnTracks,&fnTracksBranch);
    fTrackTree->SetBranchStatus("highPurity",1);
    fTrackTree->SetBranchAddress("highPurity",&fHighPurityTrackVector,&fHighPurityTrackBranch);
    fTrackTree->SetBranchStatus("trkDzFirstVtx",1);
    fTrackTree->SetBranchAddress("trkDzFirstVtx",&fTrackVertexDistanceZVector,&fTrackVertexDistanceZBranch);
    fTrackTree->SetBranchStatus("trkDzErrFirstVtx",1);
    fTrackTree->SetBranchAddress("trkDzErrFirstVtx",&fTrackVertexDistanceZErrorVector,&fTrackVertexDistanceZErrorBranch);
    fTrackTree->SetBranchStatus("trkDxyFirstVtx",1);
    fTrackTree->SetBranchAddress("trkDxyFirstVtx",&fTrackVertexDistanceXYVector,&fTrackVertexDistanceXYBranch);
    fTrackTree->SetBranchStatus("trkDxyErrFirstVtx",1);
    fTrackTree->SetBranchAddress("trkDxyErrFirstVtx",&fTrackVertexDistanceXYErrorVector,&fTrackVertexDistanceXYErrorBranch);
    fTrackTree->SetBranchStatus("trkNormChi2",1);
    fTrackTree->SetBranchAddress("trkNormChi2",&fTrackNormalizedChi2Vector,&fTrackChi2Branch);
    fTrackTree->SetBranchStatus("trkNLayers",1);
    fTrackTree->SetBranchAddress("trkNLayers",&fnHitsTrackerLayerVector,&fnHitsTrackerLayerBranch);
    fTrackTree->SetBranchStatus("trkNHits",1);
    fTrackTree->SetBranchAddress("trkNHits",&fnHitsTrackVector,&fnHitsTrackBranch);
    fTrackTree->SetBranchStatus("pfEcal",1);
    fTrackTree->SetBranchAddress("pfEcal",&fTrackEnergyEcalVector,&fTrackEnergyEcalBranch);
    fTrackTree->SetBranchStatus("pfHcal",1);
    fTrackTree->SetBranchAddress("pfHcal",&fTrackEnergyHcalVector,&fTrackEnergyHcalBranch);
    fTrackTree->SetBranchStatus("trkCharge",1);
    fTrackTree->SetBranchAddress("trkCharge",&fTrackChargeVector,&fTrackChargeBranch);
    
  } else { // Read the tree from AOD files
    
    fTrackTree->SetBranchStatus("trkPt",1);
    fTrackTree->SetBranchAddress("trkPt",&fTrackPtArray,&fTrackPtBranch);
    fTrackTree->SetBranchStatus("trkPtError",1);
    fTrackTree->SetBranchAddress("trkPtError",&fTrackPtErrorArray,&fTrackPtErrorBranch);
    fTrackTree->SetBranchStatus("trkPhi",1);
    fTrackTree->SetBranchAddress("trkPhi",&fTrackPhiArray,&fTrackPhiBranch);
    fTrackTree->SetBranchStatus("trkEta",1);
    fTrackTree->SetBranchAddress("trkEta",&fTrackEtaArray,&fTrackEtaBranch);
    fTrackTree->SetBranchStatus("nTrk",1);
    fTrackTree->SetBranchAddress("nTrk",&fnTracks,&fnTracksBranch);
    fTrackTree->SetBranchStatus("highPurity",1);
    fTrackTree->SetBranchAddress("highPurity",&fHighPurityTrackArray,&fHighPurityTrackBranch);
    fTrackTree->SetBranchStatus("trkDz1",1);
    fTrackTree->SetBranchAddress("trkDz1",&fTrackVertexDistanceZArray,&fTrackVertexDistanceZBranch);
    fTrackTree->SetBranchStatus("trkDzError1",1);
    fTrackTree->SetBranchAddress("trkDzError1",&fTrackVertexDistanceZErrorArray,&fTrackVertexDistanceZErrorBranch);
    fTrackTree->SetBranchStatus("trkDxy1",1);
    fTrackTree->SetBranchAddress("trkDxy1",&fTrackVertexDistanceXYArray,&fTrackVertexDistanceXYBranch);
    fTrackTree->SetBranchStatus("trkDxyError1",1);
    fTrackTree->SetBranchAddress("trkDxyError1",&fTrackVertexDistanceXYErrorArray,&fTrackVertexDistanceXYErrorBranch);
    fTrackTree->SetBranchStatus("trkChi2",1);
    fTrackTree->SetBranchAddress("trkChi2",&fTrackChi2Array,&fTrackChi2Branch);
    fTrackTree->SetBranchStatus("trkNdof",1);
    fTrackTree->SetBranchAddress("trkNdof",&fnTrackDegreesOfFreedomArray,&fnTrackDegreesOfFreedomBranch);
    fTrackTree->SetBranchStatus("trkNlayer",1);
    fTrackTree->SetBranchAddress("trkNlayer",&fnHitsTrackerLayerArray,&fnHitsTrackerLayerBranch);
    fTrackTree->SetBranchStatus("trkNHit",1);
    fTrackTree->SetBranchAddress("trkNHit",&fnHitsTrackArray,&fnHitsTrackBranch);
    fTrackTree->SetBranchStatus("pfEcal",1);
    fTrackTree->SetBranchAddress("pfEcal",&fTrackEnergyEcalArray,&fTrackEnergyEcalBranch);
    fTrackTree->SetBranchStatus("pfHcal",1);
    fTrackTree->SetBranchAddress("pfHcal",&fTrackEnergyHcalArray,&fTrackEnergyHcalBranch);
    fTrackTree->SetBranchStatus("trkCharge",1);
    fTrackTree->SetBranchAddress("trkCharge",&fTrackChargeArray,&fTrackChargeBranch);
  }
  
  // Connect the branches to the generator level particle tree
  if(fDataType == kPpMC || fDataType == kPbPbMC){
    fGenParticleTree->SetBranchStatus("*",0);
    fGenParticleTree->SetBranchStatus("pt",1);
    fGenParticleTree->SetBranchAddress("pt",&fGenParticlePtArray,&fGenParticlePtBranch);
    fGenParticleTree->SetBranchStatus("phi",1);
    fGenParticleTree->SetBranchAddress("phi",&fGenParticlePhiArray,&fGenParticlePhiBranch);
    fGenParticleTree->SetBranchStatus("eta",1);
    fGenParticleTree->SetBranchAddress("eta",&fGenParticleEtaArray,&fGenParticleEtaBranch);
    fGenParticleTree->SetBranchStatus("chg",1);
    fGenParticleTree->SetBranchAddress("chg",&fGenParticleChargeArray,&fGenParticleChargeBranch);
    fGenParticleTree->SetBranchStatus("sube",1);
    fGenParticleTree->SetBranchAddress("sube",&fGenParticleSubeventArray,&fGenParticleSubeventBranch);
  } // Reading track trees
  
}


/*
 * Setter for fDataType
 */
void UnfoldingForestReader::SetDataType(Int_t dataType){
  
  //Sanity check for given data type
  if(dataType < 0 || dataType > knDataTypes-1){
    cout << "ERROR: Data type input " << dataType << " is invalid in UnfoldingForestReader.cxx!" << endl;
    cout << "Please give integer between 0 and " << knDataTypes-1 << "." << endl;
    cout << "Setting data type to 0 (pp)." << endl;
    fDataType = 0;
  } else {
    
    // If the sanity check passes, set the given data type
    fDataType = dataType;
  }
}

/*
 * Connect a new tree to the reader
 */
void UnfoldingForestReader::ReadForestFromFile(TFile* inputFile){
  
  // When reading a forest, we need to check if it is AOD or MiniAOD forest as there are some differences
  // The HiForest tree is renamed to HiForestInfo in MiniAODs, so we can determine the forest type from this.
  TTree* miniAODcheck = (TTree*)inputFile->Get("HiForestInfo/HiForest");
  fIsMiniAOD = !(miniAODcheck == NULL);
  
  // Helper variable for finding the correct tree
  const char *treeName[4] = {"none","none","none","none"};
  
  // Connect a trees from the file to the reader
  fHeavyIonTree = (TTree*)inputFile->Get("hiEvtAnalyzer/HiTree");
  
  // The jet tree has different name in different datasets
  if(fDataType == kPp || fDataType == kPpMC){
    treeName[0] = "ak4CaloJetAnalyzer/t"; // Tree for calo jets
    treeName[1] = "ak4PFJetAnalyzer/t";   // Tree for PF jets
  } else if (fDataType == kPbPb || fDataType == kPbPbMC){
    treeName[0] = "akPu4CaloJetAnalyzer/t";     // Tree for calo jets
    treeName[1] = "akCs4PFJetAnalyzer/t";       // Tree for csPF jets
    treeName[2] = "akPu4PFJetAnalyzer/t";       // Tree for puPF jets
    treeName[3] = "akFlowPuCs4PFJetAnalyzer/t"; // Tree for flow subtracted csPF jets
  }
  
  fJetTree = (TTree*)inputFile->Get(treeName[fJetType]);
  
  // The track tree has different name for pp and PbPb
  if(fIsMiniAOD && (fDataType == kPbPb || fDataType == kPbPbMC)){
    fTrackTree = (TTree*)inputFile->Get("PbPbTracks/trackTree");
  } else {
    fTrackTree = (TTree*)inputFile->Get("ppTrack/trackTree");
  }
  
  // Read the generator level particle tree if we are analyzing simulation
  if(fDataType == kPpMC || fDataType == kPbPbMC){
    fGenParticleTree = (TTree*)inputFile->Get("HiGenParticleAna/hi");
  }
  
  Initialize();
}

/*
 * Connect a new tree to the reader
 */
void UnfoldingForestReader::ReadForestFromFileList(std::vector<TString> fileList){
  TFile *inputFile = TFile::Open(fileList.at(0));
  ReadForestFromFile(inputFile);
}

/*
 * Burn the current forest.
 */
void UnfoldingForestReader::BurnForest(){
  fHeavyIonTree->Delete();
  fJetTree->Delete();
  fTrackTree->Delete();
  if(fDataType == kPpMC || fDataType == kPbPbMC) fGenParticleTree->Delete();
}

/*
 * Load an event to memory
 */
void UnfoldingForestReader::GetEvent(Int_t nEvent){
  fHeavyIonTree->GetEntry(nEvent);
  fJetTree->GetEntry(nEvent);
  fTrackTree->GetEntry(nEvent);
  if(fDataType == kPpMC || fDataType == kPbPbMC) {
    fGenParticleTree->GetEntry(nEvent);
   
    // Read the numbers of generator level particles for this event
    fnGenParticles = fGenParticlePtArray->size();
  }
}

// Getter for number of events in the tree
Int_t UnfoldingForestReader::GetNEvents() const{
  return fJetPtBranch->GetEntries();
}

// Getter for number of jets in an event
Int_t UnfoldingForestReader::GetNJets() const{
  return fnJets;
}

// Getter for number of jets in an event
Int_t UnfoldingForestReader::GetNGeneratorJets() const{
  return fnGenJets;
}

// Getter for jet pT
Float_t UnfoldingForestReader::GetJetPt(Int_t iJet) const{
  return fJetPtArray[iJet];
}

// Getter for jet phi
Float_t UnfoldingForestReader::GetJetPhi(Int_t iJet) const{
  if(fJetAxis == 0) return fJetPhiArray[iJet];
  return fJetWTAPhiArray[iJet];
}

// Getter for jet eta
Float_t UnfoldingForestReader::GetJetEta(Int_t iJet) const{
  if(fJetAxis == 0) return fJetEtaArray[iJet];
  return fJetWTAEtaArray[iJet];
}

// Getter for jet raw pT
Float_t UnfoldingForestReader::GetJetRawPt(Int_t iJet) const{
  return fJetRawPtArray[iJet];
}

// Getter for maximum track pT inside a jet
Float_t UnfoldingForestReader::GetJetMaxTrackPt(Int_t iJet) const{
  return fJetMaxTrackPtArray[iJet];
}

// Check if reconstructed jet has a matching generator level jet
Bool_t UnfoldingForestReader::HasMatchingGenJet(Int_t iJet) const{
  
  // This is meaningless for data, just return true
  if(fDataType <= kPbPb) return true;
  
  // For each reconstructed jet there is a reference pT, which tells the the pT of a matched generator level jet
  // If this number is -999, it means that there are no generator level jets matching the reconstructed jet
  if(fJetRefPtArray[iJet] < 0) return false;
  return true;
}

// Get the matching generator level jet index for the given reconstructed jet
Int_t UnfoldingForestReader::GetMatchingGenIndex(Int_t iJet) const{
  
  // Ref pT array has pT for all the generator level jets that are matched with reconstructed jets
  // If our generator level pT is found from this array, it is matched to the given reconstructed jet
  for(Int_t iGenJet = 0; iGenJet < fnGenJets; iGenJet++){
    if(TMath::Abs(fJetRefEtaArray[iJet] - fGenJetEtaArray[iGenJet]) < 0.015){
      if(TMath::Abs(fJetRefPhiArray[iJet] - fGenJetPhiArray[iGenJet]) < 0.015){
        if(TMath::Abs(fJetRefPtArray[iJet] - fGenJetPtArray[iGenJet]) < 0.03*fGenJetPtArray[iGenJet]){
          return iGenJet;
        }
      }
    }
  }
  // If a matching index is not found, return -1 to show that
  return -1;  
}

// Getter for matched generator level jet pT
Float_t UnfoldingForestReader::GetMatchedGenPt(Int_t iJet) const{
  
  // This value has no meaning in real data
  if(fDataType <= kPbPb) return 0;
  
  // Find the index of the matching generator level jet
  Int_t matchedIndex = GetMatchingGenIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchedIndex == -1) return -999;

  // Return matched gen pT
  return fGenJetPtArray[matchedIndex];
}

// Get the eta of the matched generator level jet
Float_t UnfoldingForestReader::GetMatchedGenEta(Int_t iJet) const{
  
  // This value has no meaning in real data
  if(fDataType <= kPbPb) return 0;
  
  // Find the index of the matching generator level jet
  Int_t matchedIndex = GetMatchingGenIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchedIndex == -1) return -999;
  
  // Return the pT of the matching reconstructed jet
  if(fJetAxis == 0) return fGenJetEtaArray[matchedIndex];
  return fGenJetWTAEtaArray[matchedIndex];
  
}

// Get the phi of the matched reconstructed jet
Float_t UnfoldingForestReader::GetMatchedGenPhi(Int_t iJet) const{
  
  // This value has no meaning in real data
  if(fDataType <= kPbPb) return 0;
  
  // Find the index of the matching generator level jet
  Int_t matchedIndex = GetMatchingGenIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchedIndex == -1) return -999;
  
  // Return the pT of the matching reconstructed jet
  if(fJetAxis == 0) return fGenJetPhiArray[matchedIndex];
  return fGenJetWTAPhiArray[matchedIndex];
  
}

// Getter for generator level jet pT
Float_t UnfoldingForestReader::GetGeneratorJetPt(Int_t iJet) const{
  return fGenJetPtArray[iJet];
}

// Getter for generator level jet phi
Float_t UnfoldingForestReader::GetGeneratorJetPhi(Int_t iJet) const{
  if(fJetAxis == 0) return fGenJetPhiArray[iJet];
  return fGenJetWTAPhiArray[iJet];
}

// Getter for generator level jet eta
Float_t UnfoldingForestReader::GetGeneratorJetEta(Int_t iJet) const{
  if(fJetAxis == 0) return fGenJetEtaArray[iJet];
  return fGenJetWTAEtaArray[iJet];
}

// Check if generator level jet has a matching reconstructed jet
Bool_t UnfoldingForestReader::HasMatchingRecoJet(Int_t iJet) const{
  
  // This is meaningless for data, just return true
  if(fDataType <= kPbPb) return true;
  
  // Ref pT array has pT for all the generator level jets that are matched with reconstructed jets
  // If our generator level pT is found from this array, check also eta and phi
  // If all values are close by, we must have a matching jet
  Double_t jetPt = GetGeneratorJetPt(iJet);
  for(Int_t iRef = 0; iRef < fnJets; iRef++){
    if(TMath::Abs(fGenJetEtaArray[iJet] - fJetRefEtaArray[iRef]) < 0.015){
      if(TMath::Abs(fGenJetPhiArray[iJet] - fJetRefPhiArray[iRef]) < 0.015){
        if(TMath::Abs(jetPt - fJetRefPtArray[iRef]) < 0.03*jetPt) {
          return true;
        }
      }
    }
  }
  
  return false;
}

// Get the index of the matched reconstructed jet
Int_t UnfoldingForestReader::GetMatchingRecoIndex(Int_t iJet) const{
  
  // Ref pT array has pT for all the generator level jets that are matched with reconstructed jets
  // If our generator level pT is found from this array, check also eta and phi
  // If all values are close by, we must have a matching jet
  Double_t jetPt = GetGeneratorJetPt(iJet);
  for(Int_t iRef = 0; iRef < fnJets; iRef++){
    if(TMath::Abs(fGenJetEtaArray[iJet] - fJetRefEtaArray[iRef]) < 0.015){
      if(TMath::Abs(fGenJetPhiArray[iJet] - fJetRefPhiArray[iRef]) < 0.015){
        if(TMath::Abs(jetPt - fJetRefPtArray[iRef]) < 0.03*jetPt) {
          return iRef;
        }
      }
    }
  }
  
  // Return -1 to show a match is not found
  return -1;
}

// Get the pT of the matched reconstructed jet
Float_t UnfoldingForestReader::GetMatchedRecoPt(Int_t iJet) const{
  
  // This value has no meaning in real data
  if(fDataType <= kPbPb) return 0;
  
  // Find the index of the matching reconstructed jet
  Int_t matchingIndex = GetMatchingRecoIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchingIndex == -1) return -999;
  
  // Return the matching jet pT. Need raw pT for reco jets before jet corrections are in the forest
  return fJetRawPtArray[matchingIndex]; // Use this if doing jet pT correction manually
  //return fJetPtArray[matchingIndex]; // Use this if jet pT corrected in the forest
}

// Get the phi of the matched reconstructed jet
Float_t UnfoldingForestReader::GetMatchedRecoPhi(Int_t iJet) const{
  
  // This value has no meaning in real data
  if(fDataType <= kPbPb) return 0;
  
  // Find the index of the matching reconstructed jet
  Int_t matchingIndex = GetMatchingRecoIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchingIndex == -1) return -999;
  
  // Return the matching jet phi
  if(fJetAxis == 0) return fJetPhiArray[matchingIndex];
  return fJetWTAPhiArray[matchingIndex];
}

// Get the eta of the matched reconstructed jet
Float_t UnfoldingForestReader::GetMatchedRecoEta(Int_t iJet) const{
  
  // This value has no meaning in real data
  if(fDataType <= kPbPb) return 0;
  
  // Find the index of the matching reconstructed jet
  Int_t matchingIndex = GetMatchingRecoIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchingIndex == -1) return -999;
  
  // Return the matching jet eta
  if(fJetAxis == 0) return fJetEtaArray[matchingIndex];
  return fJetWTAEtaArray[matchingIndex];
}

// Getter for vertex z position
Float_t UnfoldingForestReader::GetVz() const{
  return fVertexZ;
}

// Getter for centrality. CMS has integer centrality bins from 0 to 200, thus division by 2.
Float_t UnfoldingForestReader::GetCentrality() const{
  return fHiBin/2.0;
}

// Getter for hiBin. Return 1 for negative values (for easier handling of tracking efficiency correction)
Int_t UnfoldingForestReader::GetHiBin() const{
  if(fHiBin < 0) return 1;
  return fHiBin;
}

// Getter for pT hat
Float_t UnfoldingForestReader::GetPtHat() const{
  return fPtHat;
}

// Getter for pT hat
Float_t UnfoldingForestReader::GetEventWeight() const{
  return fEventWeight;
}

// Getter for number of tracks in an event
Int_t UnfoldingForestReader::GetNTracks() const{
  return fnTracks;
}

// Getter for track pT
Float_t UnfoldingForestReader::GetTrackPt(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackPtVector->at(iTrack);
  return fTrackPtArray[iTrack];
}

// Getter for track pT error
Float_t UnfoldingForestReader::GetTrackPtError(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackPtErrorVector->at(iTrack);
  return fTrackPtErrorArray[iTrack];
}

// Getter for track phi
Float_t UnfoldingForestReader::GetTrackPhi(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackPhiVector->at(iTrack);
  return fTrackPhiArray[iTrack];
}

// Getter for track eta
Float_t UnfoldingForestReader::GetTrackEta(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackEtaVector->at(iTrack);
  return fTrackEtaArray[iTrack];
}

// Getter for high purity of the track
Bool_t UnfoldingForestReader::GetTrackHighPurity(Int_t iTrack) const{
  if(fIsMiniAOD) return fHighPurityTrackVector->at(iTrack);
  return fHighPurityTrackArray[iTrack];
}

// Getter for track distance from primary vertex in z-direction
Float_t UnfoldingForestReader::GetTrackVertexDistanceZ(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackVertexDistanceZVector->at(iTrack);
  return fTrackVertexDistanceZArray[iTrack];
}

// Getter for error of track distance from primary vertex in z-direction
Float_t UnfoldingForestReader::GetTrackVertexDistanceZError(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackVertexDistanceZErrorVector->at(iTrack);
  return fTrackVertexDistanceZErrorArray[iTrack];
}

// Getter for track distance from primary vertex in xy-direction
Float_t UnfoldingForestReader::GetTrackVertexDistanceXY(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackVertexDistanceXYVector->at(iTrack);
  return fTrackVertexDistanceXYArray[iTrack];
}

// Getter for error of track distance from primary vertex in xy-direction
Float_t UnfoldingForestReader::GetTrackVertexDistanceXYError(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackVertexDistanceXYErrorVector->at(iTrack);
  return fTrackVertexDistanceXYErrorArray[iTrack];
}

// Getter for normalized track chi2 value from reconstruction fit
Float_t UnfoldingForestReader::GetTrackNormalizedChi2(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackNormalizedChi2Vector->at(iTrack);
  return GetTrackChi2(iTrack) / (1.0*GetNTrackDegreesOfFreedom(iTrack));
}

// Getter for track chi2 value from reconstruction fit
Float_t UnfoldingForestReader::GetTrackChi2(Int_t iTrack) const{
  if(fIsMiniAOD) return -1; // Does not exist in MiniAOD forest
  return fTrackChi2Array[iTrack];
}

// Getter for number of degrees of freedom in reconstruction fit
Int_t UnfoldingForestReader::GetNTrackDegreesOfFreedom(Int_t iTrack) const{
  if(fIsMiniAOD) return -1; // Does not exist in MiniAOD forest
  return fnTrackDegreesOfFreedomArray[iTrack];
}

// Getter for number of hits in tracker layers
Int_t UnfoldingForestReader::GetNHitsTrackerLayer(Int_t iTrack) const{
  if(fIsMiniAOD) return fnHitsTrackerLayerVector->at(iTrack);
  return fnHitsTrackerLayerArray[iTrack];
}

// Getter for number of hits for the track
Int_t UnfoldingForestReader::GetNHitsTrack(Int_t iTrack) const{
  if(fIsMiniAOD) return fnHitsTrackVector->at(iTrack);
  return fnHitsTrackArray[iTrack];
}

// Getter for track energy in ECal
Float_t UnfoldingForestReader::GetTrackEnergyEcal(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackEnergyEcalVector->at(iTrack);
  return fTrackEnergyEcalArray[iTrack];
}

// Getter for track energy in HCal
Float_t UnfoldingForestReader::GetTrackEnergyHcal(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackEnergyHcalVector->at(iTrack);
  return fTrackEnergyHcalArray[iTrack];
}

// Getter for track charge
Int_t UnfoldingForestReader::GetTrackCharge(Int_t iTrack) const{
  if(fIsMiniAOD) return fTrackChargeVector->at(iTrack);
  return fTrackChargeArray[iTrack];
}

// Getter for number of generator level particles
Int_t UnfoldingForestReader::GetNGenParticles() const{
  return fnGenParticles;
}

// Getter for generator level particle pT
Float_t UnfoldingForestReader::GetGenParticlePt(Int_t iTrack) const{
  return fGenParticlePtArray->at(iTrack);
}

// Getter for generator level particle phi
Float_t UnfoldingForestReader::GetGenParticlePhi(Int_t iTrack) const{
  return fGenParticlePhiArray->at(iTrack);
}

// Getter for generator level particle eta
Float_t UnfoldingForestReader::GetGenParticleEta(Int_t iTrack) const{
  return fGenParticleEtaArray->at(iTrack);
}

// Getter for generator level particle charge
Int_t UnfoldingForestReader::GetGenParticleCharge(Int_t iTrack) const{
  return fGenParticleChargeArray->at(iTrack);
}

// Getter for generator level particle subevent index
Int_t UnfoldingForestReader::GetGenParticleSubevent(Int_t iTrack) const{
  return fGenParticleSubeventArray->at(iTrack);
}
