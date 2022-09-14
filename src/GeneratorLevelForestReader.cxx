// Implementation for GeneratorLevelForestReader

// Own includes
#include "GeneratorLevelForestReader.h"

/*
 * Default constructor
 */
GeneratorLevelForestReader::GeneratorLevelForestReader() :
  ForestReader(),
  fHeavyIonTree(0),
  fJetTree(0),
  fHltTree(0),
  fSkimTree(0),
  fTrackTree(0),
  fJetPtArray(),
  fJetPhiArray(),
  fJetEtaArray(),
  fTrackPtArray(0),
  fTrackPhiArray(0),
  fTrackEtaArray(0),
  fTrackChargeArray(0),
  fTrackSubeventArray(0),
  fTrackStatusArray(0)
{
  // Default constructor
}

/*
 * Custom constructor
 *
 *  Arguments:
 *   Int_t dataType: 0 = pp, 1 = PbPb, 2 = pp MC, 3 = PbPb MC, 4 = Local Test
 *   Int_t readMode: 0 = Regular forests, 1 = Official PYTHIA8 forest
 *   Int_t jetType: 0 = Calo jets, 1 = PF jets
 *   Int_t jetAxis: 0 = Anti-kT axis, 1 = Leading particle flow candidate axis, 2 = WTA axis
 *   Bool_t matchJets: True = Do matching for reco and gen jets. False = Do not require matching
 *   Bool_t readTrackTree: Read the track trees from the forest. Optimizes speed if tracks are not needed
 */
GeneratorLevelForestReader::GeneratorLevelForestReader(Int_t dataType, Int_t readMode, Int_t jetType, Int_t jetAxis, Bool_t matchJets, Bool_t readTrackTree) :
  ForestReader(dataType,readMode,jetType,jetAxis,matchJets,readTrackTree),
  fHeavyIonTree(0),
  fJetTree(0),
  fHltTree(0),
  fSkimTree(0),
  fTrackTree(0),
  fJetPtArray(),
  fJetPhiArray(),
  fJetEtaArray(),
  fJetRefFlavorArray(),
  fRecoJetPtArray(),
  fRecoJetEtaArray(),
  fRecoJetPhiArray(),
  fRecoJetRawPtArray(),
  fTrackPtArray(0),
  fTrackPhiArray(0),
  fTrackEtaArray(0),
  fTrackChargeArray(0),
  fTrackSubeventArray(0),
  fTrackStatusArray(0)
{
  // Custom constructor
  
}

/*
 * Copy constructor
 */
GeneratorLevelForestReader::GeneratorLevelForestReader(const GeneratorLevelForestReader& in) :
  ForestReader(in),
  fHeavyIonTree(in.fHeavyIonTree),
  fJetTree(in.fJetTree),
  fHltTree(in.fHltTree),
  fSkimTree(in.fSkimTree),
  fTrackTree(in.fTrackTree),
  fTrackPtArray(in.fTrackPtArray),
  fTrackPhiArray(in.fTrackPhiArray),
  fTrackEtaArray(in.fTrackEtaArray),
  fTrackChargeArray(in.fTrackChargeArray),
  fTrackSubeventArray(in.fTrackSubeventArray),
  fTrackStatusArray(in.fTrackStatusArray)
{
  // Copy constructor
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetPtArray[i] = in.fJetPtArray[i];
    fJetPhiArray[i] = in.fJetPhiArray[i];
    fJetEtaArray[i] = in.fJetEtaArray[i];
    fJetRefFlavorArray[i] = in.fJetRefFlavorArray[i];
    fRecoJetPtArray[i] = in.fRecoJetPtArray[i];
    fRecoJetEtaArray[i] = in.fRecoJetEtaArray[i];
    fRecoJetPhiArray[i] = in.fRecoJetPhiArray[i];
    fRecoJetRawPtArray[i] = in.fRecoJetRawPtArray[i];
  }
}

/*
 * Assignment operator
 */
GeneratorLevelForestReader& GeneratorLevelForestReader::operator=(const GeneratorLevelForestReader& in){
  // Assignment operator
  
  if (&in==this) return *this;
  
  ForestReader::operator=(in);
  
  fHeavyIonTree = in.fHeavyIonTree;
  fJetTree = in.fJetTree;
  fHltTree = in.fHltTree;
  fSkimTree = in.fSkimTree;
  fTrackTree = in.fTrackTree;
  fTrackPtArray = in.fTrackPtArray;
  fTrackPhiArray = in.fTrackPhiArray;
  fTrackEtaArray = in.fTrackEtaArray;
  fTrackChargeArray = in.fTrackChargeArray;
  fTrackSubeventArray = in.fTrackSubeventArray;
  fTrackStatusArray = in.fTrackStatusArray;
  
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetPtArray[i] = in.fJetPtArray[i];
    fJetPhiArray[i] = in.fJetPhiArray[i];
    fJetEtaArray[i] = in.fJetEtaArray[i];
    fJetRefFlavorArray[i] = in.fJetRefFlavorArray[i];
    fRecoJetPtArray[i] = in.fRecoJetPtArray[i];
    fRecoJetEtaArray[i] = in.fRecoJetEtaArray[i];
    fRecoJetPhiArray[i] = in.fRecoJetPhiArray[i];
    fRecoJetRawPtArray[i] = in.fRecoJetRawPtArray[i];
  }
  
  return *this;
}

/*
 * Destructor
 */
GeneratorLevelForestReader::~GeneratorLevelForestReader(){
  // destructor
}

/*
 * Initialization, meaning that the branches are connected to the tree
 */
void GeneratorLevelForestReader::Initialize(){
  
  // Connect the branches of the heavy ion tree
  fHeavyIonTree->SetBranchStatus("*",0);
  fHeavyIonTree->SetBranchStatus("vz",1);
  fHeavyIonTree->SetBranchAddress("vz",&fVertexZ,&fHiVzBranch);
  fHeavyIonTree->SetBranchStatus("hiBin",1);
  fHeavyIonTree->SetBranchAddress("hiBin",&fHiBin,&fHiBinBranch);
  fHeavyIonTree->SetBranchStatus("pthat",1);
  fHeavyIonTree->SetBranchAddress("pthat",&fPtHat,&fPtHatBranch); // pT hat only for MC
  
  // Read event weight for 2018 simulation
  fHeavyIonTree->SetBranchStatus("weight",1);
  fHeavyIonTree->SetBranchAddress("weight",&fEventWeight,&fEventWeightBranch);
  
  // Connect the branches to the jet tree
  const char *jetAxis[3] = {"jt","jt","WTA"};
  const char *genJetAxis[3] = {"","","WTA"};
  char branchName[30];
  
  fJetTree->SetBranchStatus("*",0);
  fJetTree->SetBranchStatus("genpt",1);
  fJetTree->SetBranchAddress("genpt",&fJetPtArray,&fJetPtBranch);
  
  // If specified, select WTA axis for jet phi
  sprintf(branchName,"%sgenphi",genJetAxis[fJetAxis]);
  fJetTree->SetBranchStatus(branchName,1);
  fJetTree->SetBranchAddress(branchName,&fJetPhiArray,&fJetPhiBranch);
  
  // If specified, select WTA axis for jet eta
  sprintf(branchName,"%sgeneta",genJetAxis[fJetAxis]);
  fJetTree->SetBranchStatus(branchName,1);
  fJetTree->SetBranchAddress(branchName,&fJetEtaArray,&fJetEtaBranch);
  
  fJetTree->SetBranchStatus("ngen",1);
  fJetTree->SetBranchAddress("ngen",&fnJets,&fJetRawPtBranch); // Reuse a branch from ForestReader that is not otherwise needed here
  
  // If we are matching jets, connect the ref pT array
  if(fMatchJets){
    fJetTree->SetBranchStatus("refpt",1);
    fJetTree->SetBranchAddress("refpt",&fJetRefPtArray,&fJetRefPtBranch);
    fJetTree->SetBranchStatus("refparton_flavor",1);
    fJetTree->SetBranchAddress("refparton_flavor",&fJetRefFlavorArray,&fJetRefFlavorBranch);
    fJetTree->SetBranchStatus("jtpt",1);
    fJetTree->SetBranchAddress("jtpt",&fRecoJetPtArray,&fJetMatchedPtBranch);
    fJetTree->SetBranchStatus("rawpt",1);
    fJetTree->SetBranchAddress("rawpt",&fRecoJetRawPtArray,&fHighPurityTrackBranch); // Reuse a branch from ForestReader that is not otherwise needed here
    
    // If specified, select WTA axis for jet phi
    sprintf(branchName,"%sphi",jetAxis[fJetAxis]);
    fJetTree->SetBranchStatus(branchName,1);
    fJetTree->SetBranchAddress(branchName,&fRecoJetPhiArray,&fJetMatchedPhiBranch);
    
    // If specified, select WTA axis for jet eta
    sprintf(branchName,"%seta",jetAxis[fJetAxis]);
    fJetTree->SetBranchStatus(branchName,1);
    fJetTree->SetBranchAddress(branchName,&fRecoJetEtaArray,&fJetMatchedEtaBranch);
    
    fJetTree->SetBranchStatus("nref",1);
    fJetTree->SetBranchAddress("nref",&fnMatchedJets,&fnTrackDegreesOfFreedomBranch); // Reuse a branch from ForestReader that is not otherwise needed here
  }
  
  // Helper variable for choosing correct branches in HLT tree
  const char *branchNameHlt[2] = {"none","none"};
  
  // Connect the branches to the HLT tree
  fHltTree->SetBranchStatus("*",0);
  if(fDataType == kPp){ // pp data
    branchNameHlt[0] = "HLT_AK4CaloJet80_Eta5p1_v1";
    branchNameHlt[1] = "HLT_AK4PFJet80_Eta5p1_v1";
    fHltTree->SetBranchStatus(branchNameHlt[fJetType],1);
    fHltTree->SetBranchAddress(branchNameHlt[fJetType],&fCaloJetFilterBit,&fCaloJetFilterBranch);
  } else if (fDataType == kPpMC){
    fCaloJetFilterBit = 1; // No filtering for Monte Carlo
  } else if (fDataType == kPbPb || (fDataType == kPbPbMC && fReadMode == 2019)){ // PbPb or MC is specifically required
    fHltTree->SetBranchStatus("HLT_HIPuAK4CaloJet100Eta5p1_v1",1);
    fHltTree->SetBranchAddress("HLT_HIPuAK4CaloJet100Eta5p1_v1",&fCaloJetFilterBit,&fCaloJetFilterBranch);
  } else {
    fCaloJetFilterBit = 1;  // No filter for local test of MC if not specifically required
  }
  fCaloJetFilterBitPrescale = 1; // Set the prescaled filter bit to 1. Only relevant for minimum bias PbPb (data skim)
  
  // Connect the branches to the skim tree (different for pp and PbPb Monte Carlo)
  fSkimTree->SetBranchStatus("*",0);
  if(fDataType == kPpMC){ // pp MC
    fSkimTree->SetBranchStatus("pPAprimaryVertexFilter",1);
    fSkimTree->SetBranchAddress("pPAprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    fSkimTree->SetBranchStatus("pBeamScrapingFilter",1);
    fSkimTree->SetBranchAddress("pBeamScrapingFilter",&fBeamScrapingFilterBit,&fBeamScrapingBranch);
    fSkimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
    fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
    fCollisionEventSelectionFilterBit = 1;  // No collision event selection filter for pp
    fHfCoincidenceFilterBit = 1; // No HF energy coincidence requirement for pp
    fClusterCompatibilityFilterBit = 1; // No cluster compatibility requirement for pp
  } else if (fDataType == kPbPbMC){ // PbPb MC
    
    // Primary vertex has at least two tracks, is within 25 cm in z-rirection and within 2 cm in xy-direction
    fSkimTree->SetBranchStatus("pprimaryVertexFilter",1);
    fSkimTree->SetBranchAddress("pprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    
    // Cut on noise on HCAL
    fSkimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
    fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
    
    // Combination of hfCoincFilter2Th4, pprimaryVertexFilter and pclusterCompatibilityFilter
    fSkimTree->SetBranchStatus("collisionEventSelectionAODv2",1);
    fSkimTree->SetBranchAddress("collisionEventSelectionAODv2", &fCollisionEventSelectionFilterBit, &fCollisionEventSelectionBranch);
    
    // Have at least two HF towers on each side of the detector with an energy deposit of 4 GeV
    fSkimTree->SetBranchStatus("phfCoincFilter2Th4",1);
    fSkimTree->SetBranchAddress("phfCoincFilter2Th4", &fHfCoincidenceFilterBit, &fHfCoincidenceBranch);

    // Calculated from pixel clusters. Ensures that measured and predicted primary vertices are compatible
    fSkimTree->SetBranchStatus("pclusterCompatibilityFilter",1);
    fSkimTree->SetBranchAddress("pclusterCompatibilityFilter", &fClusterCompatibilityFilterBit, &fClusterCompatibilityBranch);
    
    fBeamScrapingFilterBit = 1;  // No beam scraping filter for PbPb
  } else { // Local test
    fPrimaryVertexFilterBit = 1;
    fBeamScrapingFilterBit = 1;
    fCollisionEventSelectionFilterBit = 1;
    fHBHENoiseFilterBit = 1;
    fHfCoincidenceFilterBit = 1;
    fClusterCompatibilityFilterBit = 1;
  }
  
  // Connect the branches to the track tree
  if(fReadTrackTree){
    fTrackTree->SetBranchStatus("*",0);
    fTrackTree->SetBranchStatus("pt",1);
    fTrackTree->SetBranchAddress("pt",&fTrackPtArray,&fTrackPtBranch);
    fTrackTree->SetBranchStatus("phi",1);
    fTrackTree->SetBranchAddress("phi",&fTrackPhiArray,&fTrackPhiBranch);
    fTrackTree->SetBranchStatus("eta",1);
    fTrackTree->SetBranchAddress("eta",&fTrackEtaArray,&fTrackEtaBranch);
    fTrackTree->SetBranchStatus("chg",1);
    fTrackTree->SetBranchAddress("chg",&fTrackChargeArray,&fTrackPtErrorBranch);  // Reuse a branch from ForestReader that is not otherwise needed here
    fTrackTree->SetBranchStatus("sube",1);
    fTrackTree->SetBranchAddress("sube",&fTrackSubeventArray,&fTrackChi2Branch);  // Reuse a branch from ForestReader that is not otherwise needed here
    /*if(fDataType != kLocalTest && fReadMode == 0) {
     fTrackTree->SetBranchStatus("sta",1);
     fTrackTree->SetBranchAddress("sta",&fTrackStatusArray,&fTrackEnergyEcalBranch); // Reuse a branch from ForestReader that is not otherwise needed here. Not available for local test or PYTHIA8 forest
     } */ // Quick test. New PbPb MC files do not have this branch
  } // Reading track trees

}

/*
 * Connect a new tree to the reader
 */
void GeneratorLevelForestReader::ReadForestFromFile(TFile *inputFile){
  
  // Helper variable for finding the correct tree
  const char *treeName[4] = {"none","none","none","none"};
  
  // Connect a trees from the file to the reader
  fHeavyIonTree = (TTree*)inputFile->Get("hiEvtAnalyzer/HiTree");
  fHltTree = (TTree*)inputFile->Get("hltanalysis/HltTree");
  fSkimTree = (TTree*)inputFile->Get("skimanalysis/HltTree");
  
  // The jet tree has different name in different datasets
  if(fDataType == kPp || fDataType == kPpMC){
    treeName[0] = "ak4CaloJetAnalyzer/t"; // Tree for calo jets
    treeName[1] = "ak4PFJetAnalyzer/t";   // Tree for PF jets
  } else if (fDataType == kPbPb || fDataType == kPbPbMC){
    treeName[0] = "akPu4CaloJetAnalyzer/t";     // Tree for calo jets
    treeName[1] = "akCs4PFJetAnalyzer/t";       // Tree for PF jets
    treeName[2] = "akPu4PFJetAnalyzer/t";       // Tree for PF jets
    treeName[3] = "akFlowPuCs4PFJetAnalyzer/t"; // Tree for flow subtracted csPF jets
  } else if (fDataType == kLocalTest){
    treeName[0] = "ak4PFJetAnalyzer/t";  // Only PF jets in local test file
    treeName[1] = "ak4PFJetAnalyzer/t";  // Only PF jets in local test file
  }
  fJetTree = (TTree*)inputFile->Get(treeName[fJetType]);
  
  // The track tree is different than in other types of forests
  if(fReadTrackTree) fTrackTree = (TTree*)inputFile->Get("HiGenParticleAna/hi");
  
  // Connect branches to trees
  Initialize();
}

/*
 * Connect a new tree to the reader
 */
void GeneratorLevelForestReader::ReadForestFromFileList(std::vector<TString> fileList){
  TFile *inputFile = TFile::Open(fileList.at(0));
  ReadForestFromFile(inputFile);
}

/*
 * Burn the current forest.
 */
void GeneratorLevelForestReader::BurnForest(){
  fHeavyIonTree->Delete();
  fJetTree->Delete();
  fHltTree->Delete();
  fSkimTree->Delete();
  fTrackTree->Delete();
}

/*
 * Load an event to memory
 */
void GeneratorLevelForestReader::GetEvent(Int_t nEvent){
  fHeavyIonTree->GetEntry(nEvent);
  fJetTree->GetEntry(nEvent);
  fHltTree->GetEntry(nEvent);
  fSkimTree->GetEntry(nEvent);
  if(fReadTrackTree) fTrackTree->GetEntry(nEvent);
  
  // Read the numbers of tracks for this event
  if(fReadTrackTree) fnTracks = fTrackPtArray->size();
}

// Getter for jet pT
Float_t GeneratorLevelForestReader::GetJetPt(Int_t iJet) const{
  return fJetPtArray[iJet];
}

// Getter for jet phi
Float_t GeneratorLevelForestReader::GetJetPhi(Int_t iJet) const{
  return fJetPhiArray[iJet];
}

// Getter for jet eta
Float_t GeneratorLevelForestReader::GetJetEta(Int_t iJet) const{
  return fJetEtaArray[iJet];
}

// Getter for jet raw pT (not relevant for generator jets, just return value that passes cuts)
Float_t GeneratorLevelForestReader::GetJetRawPt(Int_t iJet) const{
  return 2; // The cut is on maxTrackPt/rawPt. Giving 1 and 2 here passes the analysis cuts
}

// Getter for maximum track pT inside a jet (not relevant for generator jets, just return value that passes cuts)
Float_t GeneratorLevelForestReader::GetJetMaxTrackPt(Int_t iJet) const{
  return 1; // The cut is on maxTrackPt/rawPt. Giving 1 and 2 here passes the analysis cuts
}

// Getter for track pT
Float_t GeneratorLevelForestReader::GetTrackPt(Int_t iTrack) const{
  return fTrackPtArray->at(iTrack);
}

// Getter for track phi
Float_t GeneratorLevelForestReader::GetTrackPhi(Int_t iTrack) const{
  return fTrackPhiArray->at(iTrack);
}

// Getter for track eta
Float_t GeneratorLevelForestReader::GetTrackEta(Int_t iTrack) const{
  return fTrackEtaArray->at(iTrack);
}

// Getter for track charge
Int_t GeneratorLevelForestReader::GetTrackCharge(Int_t iTrack) const{
  return fTrackChargeArray->at(iTrack);
}

// Getter for track subevent index
Int_t GeneratorLevelForestReader::GetTrackSubevent(Int_t iTrack) const{
  return fTrackSubeventArray->at(iTrack);
}

// Getter for track MC status
Int_t GeneratorLevelForestReader::GetTrackMCStatus(Int_t iTrack) const{
  /*if(fDataType == kLocalTest || fReadMode == 1 || fReadMode == 2) return 1;
  return fTrackStatusArray->at(iTrack);*/ // There is no status for MC tracks in new PbPb MC trees.
  return 1;
}

// Getter for track pT error (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackPtError(Int_t iTrack) const{
  return 0; // Setting all errors to 0 always passes the track quality cut
}

// Getter for high purity of the track (not relevant for generator tracks)
Bool_t GeneratorLevelForestReader::GetTrackHighPurity(Int_t iTrack) const{
  return true; // All the generator tracks are of high purity
}

// Getter for track distance from primary vertex in z-direction (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackVertexDistanceZ(Int_t iTrack) const{
  return 0; // The cut is on distance/error. Setting this to 0 and error to 1 always passes the cut
}

// Getter for error of track distance from primary vertex in z-direction (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackVertexDistanceZError(Int_t iTrack) const{
  return 1; // The cut is on distance/error. Setting this to 1 and distance to 0 always passes the cut
}

// Getter for track distance from primary vertex in xy-direction (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackVertexDistanceXY(Int_t iTrack) const{
  return 0; // The cut is on distance/error. Setting this to 0 and error to 1 always passes the cut
}

// Getter for error of track distance from primary vertex in xy-direction (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackVertexDistanceXYError(Int_t iTrack) const{
  return 1; // The cut is on distance/error. Setting this to 1 and distance to 0 always passes the cut
}

// Getter for track chi2 value from reconstruction fit (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackChi2(Int_t iTrack) const{
  return 1; // Note: The chi2 quality cut is disabled for generator level tracks in the main analysis code
}

// Getter for number of degrees of freedom in reconstruction fit (not relevant for generator tracks)
Int_t GeneratorLevelForestReader::GetNTrackDegreesOfFreedom(Int_t iTrack) const{
  return 1; // Note: The chi2 quality cut is disabled for generator level tracks in the main analysis code
}

// Getter for number of hits in tracker layers (not relevant for generator tracks)
Int_t GeneratorLevelForestReader::GetNHitsTrackerLayer(Int_t iTrack) const{
  return 1; // Note: The chi2 quality cut is disabled for generator level tracks in the main analysis code
}

// Getter for number of hits for the track (not relevant for generator tracks)
Int_t GeneratorLevelForestReader::GetNHitsTrack(Int_t iTrack) const{
  return 1; // Note: The cut on NHits is disabled for generator level tracks in the main analysis code
}

// Getter for track energy in ECal (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackEnergyEcal(Int_t iTrack) const{
  return 1; // Note: The cut on Et is disabled for generator level tracks in the main analysis code
}

// Getter for track energy in HCal (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackEnergyHcal(Int_t iTrack) const{
  return 1; // Note: The cut on Et is disabled for generator level tracks in the main analysis code
}

// Getter for particle flow candidate ID (not relevant for generator tracks)
Int_t GeneratorLevelForestReader::GetParticleFlowCandidateId(Int_t iCandidate) const{
  return 1; // Return 1 as we use regular generator lavel tracks to determine leading particle flow candidate
}

// Getter for particle flow candidate pT (just regular tracks in generator level)
Float_t GeneratorLevelForestReader::GetParticleFlowCandidatePt(Int_t iCandidate) const{
  return fTrackPtArray->at(iCandidate); // Use regular generated particle pT for particle flow candidates
}

// Getter for particle flow candidate phi (just regular tracks in generator level)
Float_t GeneratorLevelForestReader::GetParticleFlowCandidatePhi(Int_t iCandidate) const{
  return fTrackPhiArray->at(iCandidate); // Use regular generated particle phi for particle flow candidates
}

// Getter for particle flow candidate eta (just regular tracks in generator level)
Float_t GeneratorLevelForestReader::GetParticleFlowCandidateEta(Int_t iCandidate) const{
  return fTrackEtaArray->at(iCandidate); // Use regular generated particle eta for particle flow candidates
}

// Getter number of particle flow candidates in an event (just regular tracks in generator level)
Int_t GeneratorLevelForestReader::GetNParticleFlowCandidates() const{
  return fnTracks; // Use regular tracks in generator level
}

// Check if generator level jet has a matching reconstructed jet
Bool_t GeneratorLevelForestReader::HasMatchingJet(Int_t iJet) const{
  
  // If not matching jets, just tell that everything is fine
  if(!fMatchJets) return true;
  
  // Ref pT array has pT for all the generator level jets that are matched with reconstructed jets
  // If our generator level pT is found from this array, it has a matching reconstructed jet
  Double_t jetPt = GetJetPt(iJet);
  for(Int_t iRef = 0; iRef < fnMatchedJets; iRef++){
    if(TMath::Abs(jetPt - fJetRefPtArray[iRef]) < 0.001) return true;
  }
  
  return false;
}

// Get the pT of the matched reconstructed jet
Int_t GeneratorLevelForestReader::GetMatchingIndex(Int_t iJet) const{
  
  // Ref pT array has pT for all the generator level jets that are matched with reconstructed jets
  // If our generator level pT is found from this array, it has a matching reconstructed jet
  Double_t jetPt = GetJetPt(iJet);
  Int_t matchingIndex = -1;
  for(Int_t iRef = 0; iRef < fnMatchedJets; iRef++){
    if(TMath::Abs(jetPt - fJetRefPtArray[iRef]) < 0.001){
      matchingIndex = iRef;
      break;
    }
  }
  
  // Return the matching index
  return matchingIndex;
}

// Get the pT of the matched reconstructed jet
Float_t GeneratorLevelForestReader::GetMatchedPt(Int_t iJet) const{
  
  // If not matching jets, just return 0 as this value has no meaning
  if(!fMatchJets) return 0;
  
  // Find the index of the matching reconstructed jet
  Int_t matchingIndex = GetMatchingIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchingIndex == -1) return -999;
  
  // Return the matching jet pT. Need raw pT for reco jets before jet corrections are in the forest
  if(fReadMode > 2000) return fRecoJetRawPtArray[matchingIndex];
  return fRecoJetPtArray[matchingIndex];
}

// Get the phi of the matched reconstructed jet
Float_t GeneratorLevelForestReader::GetMatchedPhi(Int_t iJet) const{
  
  // If not matching jets, just return 0 as this value has no meaning
  if(!fMatchJets) return 0;
  
  // Find the index of the matching reconstructed jet
  Int_t matchingIndex = GetMatchingIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchingIndex == -1) return -999;
  
  // Return the matching jet phi
  return fRecoJetPhiArray[matchingIndex];
}

// Get the eta of the matched reconstructed jet
Float_t GeneratorLevelForestReader::GetMatchedEta(Int_t iJet) const{
  
  // If not matching jets, just return 0 as this value has no meaning
  if(!fMatchJets) return 0;
  
  // Find the index of the matching reconstructed jet
  Int_t matchingIndex = GetMatchingIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchingIndex == -1) return -999;
  
  // Return the matching jet eta
  return fRecoJetEtaArray[matchingIndex];
}

// Check if generator level jet has a matching reconstructed jet
Int_t GeneratorLevelForestReader::GetPartonFlavor(Int_t iJet) const{
  
  // If not matching jets, just return 0 as this value has no meaning
  if(!fMatchJets) return 0;
  
  // Find the index of the matching reconstructed jet
  Int_t matchingIndex = GetMatchingIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchingIndex == -1) return -999;
  
  // Return the matching parton flavor
  return fJetRefFlavorArray[matchingIndex];
}
