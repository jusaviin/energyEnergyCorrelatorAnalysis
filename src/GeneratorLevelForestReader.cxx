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
 *   Int_t dataType: 0 = pp, 1 = PbPb, 2 = pp MC, 3 = PbPb MC
 *   Int_t useJetTrigger: 0 = Do not use any triggers, > 0 = Require jet trigger
 *   Int_t jetType: 0 = Calo jets, 1 = PF jets
 *   Int_t jetAxis: 0 = Anti-kT axis, 1 = WTA axis
 *   Int_t matchJets: non-0 = Do matching for reco and gen jets. 0 = Do not require matching
 *   Bool_t readTrackTree: Read the track trees from the forest. Optimizes speed if tracks are not needed
 */
GeneratorLevelForestReader::GeneratorLevelForestReader(Int_t dataType, Int_t useJetTrigger, Int_t jetType, Int_t jetAxis, Int_t matchJets, Bool_t readTrackTree) :
  ForestReader(dataType,useJetTrigger,jetType,jetAxis,matchJets,readTrackTree),
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
  
  // Read event weight for simulation
  fHeavyIonTree->SetBranchStatus("weight",1);
  fHeavyIonTree->SetBranchAddress("weight",&fEventWeight,&fEventWeightBranch);
  
  // Connect the branches to the jet tree
  const char* jetAxis[2] = {"jt","WTA"};
  const char* branchName;
  
  fJetTree->SetBranchStatus("*",0);
  fJetTree->SetBranchStatus("genpt",1);
  fJetTree->SetBranchAddress("genpt",&fJetPtArray,&fJetPtBranch);
  
  // If specified, select WTA axis for jet phi
  fJetTree->SetBranchStatus("genphi",1);
  fJetTree->SetBranchAddress("genphi",&fJetPhiArray,&fJetPhiBranch);
  fJetTree->SetBranchStatus("WTAgenphi",1);
  fJetTree->SetBranchAddress("WTAgenphi",&fJetWTAPhiArray,&fJetWTAPhiBranch);
  
  // If specified, select WTA axis for jet eta
  fJetTree->SetBranchStatus("geneta",1);
  fJetTree->SetBranchAddress("geneta",&fJetEtaArray,&fJetEtaBranch);
  fJetTree->SetBranchStatus("WTAgeneta",1);
  fJetTree->SetBranchAddress("WTAgeneta",&fJetWTAEtaArray,&fJetWTAEtaBranch);
  
  fJetTree->SetBranchStatus("ngen",1);
  fJetTree->SetBranchAddress("ngen",&fnJets,&fJetRawPtBranch); // Reuse a branch from ForestReader that is not otherwise needed here
  
  // If we are matching jets, connect the ref pT array
  if(fMatchJets){
    fJetTree->SetBranchStatus("refpt",1);
    fJetTree->SetBranchAddress("refpt",&fJetRefPtArray,&fJetRefPtBranch);
    fJetTree->SetBranchStatus("refeta",1);
    fJetTree->SetBranchAddress("refeta",&fJetRefEtaArray,&fJetRefEtaBranch);
    fJetTree->SetBranchStatus("refphi",1);
    fJetTree->SetBranchAddress("refphi",&fJetRefPhiArray,&fJetRefPhiBranch);
    fJetTree->SetBranchStatus("refparton_flavorForB",1);
    fJetTree->SetBranchAddress("refparton_flavorForB",&fJetRefFlavorArray,&fJetRefFlavorBranch);
    fJetTree->SetBranchStatus("jtpt",1);
    fJetTree->SetBranchAddress("jtpt",&fRecoJetPtArray,&fJetMatchedPtBranch);
    fJetTree->SetBranchStatus("rawpt",1);
    fJetTree->SetBranchAddress("rawpt",&fRecoJetRawPtArray,&fHighPurityTrackBranch); // Reuse a branch from ForestReader that is not otherwise needed here
    
    // If specified, select WTA axis for jet phi
    branchName = Form("%sphi",jetAxis[fJetAxis]);
    fJetTree->SetBranchStatus(branchName,1);
    fJetTree->SetBranchAddress(branchName,&fRecoJetPhiArray,&fJetMatchedPhiBranch);
    
    // If specified, select WTA axis for jet eta
    branchName = Form("%seta",jetAxis[fJetAxis]);
    fJetTree->SetBranchStatus(branchName,1);
    fJetTree->SetBranchAddress(branchName,&fRecoJetEtaArray,&fJetMatchedEtaBranch);
    
    fJetTree->SetBranchStatus("nref",1);
    fJetTree->SetBranchAddress("nref",&fnMatchedJets,&fnTrackDegreesOfFreedomBranch); // Reuse a branch from ForestReader that is not otherwise needed here
  }
  
  // Connect the branches to the HLT tree
  if(fUseJetTrigger){
    
    fHltTree->SetBranchStatus("*",0);
    
    if(fDataType == kPp || fDataType == kPpMC){ // pp data or MC

      // Calo jet 80 trigger
      fHltTree->SetBranchStatus("HLT_HIAK4CaloJet60_v1", 1);
      fHltTree->SetBranchAddress("HLT_HIAK4CaloJet60_v1",&fCaloJet60FilterBit,&fCaloJet60FilterBranch);

      // Calo jet 80 trigger
      fHltTree->SetBranchStatus("HLT_HIAK4CaloJet80_v1",1);
      fHltTree->SetBranchAddress("HLT_HIAK4CaloJet80_v1",&fCaloJet80FilterBit,&fCaloJet80FilterBranch);
      
      // Calo jet 100 trigger
      fHltTree->SetBranchStatus("HLT_HIAK4CaloJet100_v1",1);
      fHltTree->SetBranchAddress("HLT_HIAK4CaloJet100_v1",&fCaloJet100FilterBit,&fCaloJet100FilterBranch);
      
    } else { // PbPb data or MC

      // Calo jet 80 trigger
      fHltTree->SetBranchStatus("HLT_HIPuAK4CaloJet60Eta5p1_v1",1);
      fHltTree->SetBranchAddress("HLT_HIPuAK4CaloJet60Eta5p1_v1",&fCaloJet60FilterBit,&fCaloJet60FilterBranch);

      // Calo jet 80 trigger
      fHltTree->SetBranchStatus("HLT_HIPuAK4CaloJet80Eta5p1_v1",1);
      fHltTree->SetBranchAddress("HLT_HIPuAK4CaloJet80Eta5p1_v1",&fCaloJet80FilterBit,&fCaloJet80FilterBranch);
      
      // Calo jet 100 trigger
      fHltTree->SetBranchStatus("HLT_HIPuAK4CaloJet100Eta5p1_v1",1);
      fHltTree->SetBranchAddress("HLT_HIPuAK4CaloJet100Eta5p1_v1",&fCaloJet100FilterBit,&fCaloJet100FilterBranch);
      
    }
  } else {
    fCaloJet60FilterBit = 1;
    fCaloJet80FilterBit = 1;
    fCaloJet100FilterBit = 1;
  }
  
  // Connect the branches to the skim tree (different for pp and PbPb Monte Carlo)
  fSkimTree->SetBranchStatus("*",0);
  if(fDataType == kPpMC){ // pp MC
    fSkimTree->SetBranchStatus("pPAprimaryVertexFilter",1);
    fSkimTree->SetBranchAddress("pPAprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    fSkimTree->SetBranchStatus("pBeamScrapingFilter",1);
    fSkimTree->SetBranchAddress("pBeamScrapingFilter",&fBeamScrapingFilterBit,&fBeamScrapingBranch);
    
    if(fIsMiniAOD){
      fHBHENoiseFilterBit = 1; // This filter bit is not available in MiniAODs
    } else {
      fSkimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
      fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
    }
    
    fHfCoincidenceFilterBit = 1; // No HF energy coincidence requirement for pp
    fClusterCompatibilityFilterBit = 1; // No cluster compatibility requirement for pp
  } else { // PbPb MC
    
    // Primary vertex has at least two tracks, is within 25 cm in z-rirection and within 2 cm in xy-direction
    fSkimTree->SetBranchStatus("pprimaryVertexFilter",1);
    fSkimTree->SetBranchAddress("pprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    
    if(fIsMiniAOD){
      
      fHBHENoiseFilterBit = 1; // This filter bit is not available in MiniAODs
      
      // Have at least two HF towers on each side of the detector with an energy deposit of 4 GeV
      fSkimTree->SetBranchStatus("pphfCoincFilter2Th4",1);
      fSkimTree->SetBranchAddress("pphfCoincFilter2Th4", &fHfCoincidenceFilterBit, &fHfCoincidenceBranch);
      
    } else {
      
      // Cut on noise on HCAL
      fSkimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
      fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
      
      // Have at least two HF towers on each side of the detector with an energy deposit of 4 GeV
      fSkimTree->SetBranchStatus("phfCoincFilter2Th4",1);
      fSkimTree->SetBranchAddress("phfCoincFilter2Th4", &fHfCoincidenceFilterBit, &fHfCoincidenceBranch);
      
    }

    // Calculated from pixel clusters. Ensures that measured and predicted primary vertices are compatible
    fSkimTree->SetBranchStatus("pclusterCompatibilityFilter",1);
    fSkimTree->SetBranchAddress("pclusterCompatibilityFilter", &fClusterCompatibilityFilterBit, &fClusterCompatibilityBranch);
    
    fBeamScrapingFilterBit = 1;  // No beam scraping filter for PbPb
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
  } // Reading track trees

}

/*
 * Connect a new tree to the reader
 */
void GeneratorLevelForestReader::ReadForestFromFile(TFile *inputFile){
  
  // When reading a forest, we need to check if it is AOD or MiniAOD forest as there are some differences
  // The HiForest tree is renamed to HiForestInfo in MiniAODs, so we can determine the forest type from this.
  TTree* miniAODcheck = (TTree*)inputFile->Get("HiForestInfo/HiForest");
  fIsMiniAOD = !(miniAODcheck == NULL);
  
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
  } else { // PbPb data or MC
    treeName[0] = "akPu4CaloJetAnalyzer/t";     // Tree for calo jets
    treeName[1] = "akCs4PFJetAnalyzer/t";       // Tree for PF jets
    treeName[2] = "akPu4PFJetAnalyzer/t";       // Tree for PF jets
    treeName[3] = "akFlowPuCs4PFJetAnalyzer/t"; // Tree for flow subtracted csPF jets
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
  if(fReadTrackTree) {
    fTrackTree->GetEntry(nEvent);
    
    // Read the numbers of tracks for this event
    fnTracks = fTrackPtArray->size();
  }
}

// Getter for jet pT
Float_t GeneratorLevelForestReader::GetJetPt(Int_t iJet) const{
  return fJetPtArray[iJet];
}

// Getter for jet phi
Float_t GeneratorLevelForestReader::GetJetPhi(Int_t iJet) const{
  if(fJetAxis == 0) return fJetPhiArray[iJet];
  return fJetWTAPhiArray[iJet];
}

// Getter for jet eta
Float_t GeneratorLevelForestReader::GetJetEta(Int_t iJet) const{
  if(fJetAxis == 0) return fJetEtaArray[iJet];
  return fJetWTAEtaArray[iJet];
}

// Getter for jet raw pT (not relevant for generator jets, just return regular jet pT)
Float_t GeneratorLevelForestReader::GetJetRawPt(Int_t iJet) const{
  return fJetPtArray[iJet];;
}

// Getter for maximum track pT inside a jet (not relevant for generator jets, return negative value to show this)
Float_t GeneratorLevelForestReader::GetJetMaxTrackPt(Int_t iJet) const{
  return -1;
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

// Getter for normalized track chi2 value from reconstruction fit (not relevant for generator tracks)
Float_t GeneratorLevelForestReader::GetTrackNormalizedChi2(Int_t iTrack) const{
  return 1; // Note: The chi2 quality cut is disabled for generator level tracks in the main analysis code
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

// Check if generator level jet has a matching reconstructed jet
Bool_t GeneratorLevelForestReader::HasMatchingJet(Int_t iJet) const{
  
  // If not matching jets, just tell that everything is fine
  if(!fMatchJets) return true;
  
  // Ref pT array has pT for all the generator level jets that are matched with reconstructed jets
  // If our generator level pT is found from this array, check also eta and phi
  // If all values are close by, we must have a matching jet
  Double_t jetPt = GetJetPt(iJet);
  for(Int_t iRef = 0; iRef < fnMatchedJets; iRef++){
    if(TMath::Abs(fJetEtaArray[iJet] - fJetRefEtaArray[iRef]) < 0.015){
      if(TMath::Abs(fJetPhiArray[iJet] - fJetRefPhiArray[iRef]) < 0.015){
        if(TMath::Abs(jetPt - fJetRefPtArray[iRef]) < 0.03*jetPt) {
          return true;
        }
      }
    }
  }
  
  return false;
}

// Get the pT of the matched reconstructed jet
Int_t GeneratorLevelForestReader::GetMatchingIndex(Int_t iJet) const{
  
  // Ref pT array has pT for all the generator level jets that are matched with reconstructed jets
  // If our generator level pT is found from this array, check also eta and phi
  // If all values are close by, we must have a matching jet
  Double_t jetPt = GetJetPt(iJet);
  Int_t matchingIndex = -1;
  for(Int_t iRef = 0; iRef < fnMatchedJets; iRef++){
    if(TMath::Abs(fJetEtaArray[iJet] - fJetRefEtaArray[iRef]) < 0.015){
      if(TMath::Abs(fJetPhiArray[iJet] - fJetRefPhiArray[iRef]) < 0.015){
        if(TMath::Abs(jetPt - fJetRefPtArray[iRef]) < 0.03*jetPt) {
          matchingIndex = iRef;
          break;
        }
      }
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
  //return fRecoJetRawPtArray[matchingIndex]; // Use this if doing jet pT correction manually
  return fRecoJetPtArray[matchingIndex]; // Use this if jet pT corrected in the forest
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
