// Implementation for HighForestReader

// Own includes
#include "HighForestReader.h"

/*
 * Default constructor
 */
HighForestReader::HighForestReader() :
  ForestReader(),
  fHeavyIonTree(0),
  fJetTree(0),
  fHltTree(0),
  fSkimTree(0),
  fTrackTree(0),
  fnJetsBranch(0),
  fnTracksBranch(0),
  fTrackAlgorithmBranch(0),
  fTrackOriginalAlgorithmBranch(0),
  fTrackMVABranch(0),
  fJetPtArray(),
  fJetPhiArray(),
  fJetEtaArray(),
  fJetRawPtArray(),
  fJetMaxTrackPtArray(),
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
  fTrackAlgorithmArray(),
  fTrackOriginalAlgorithmArray(),
  fTrackMVAArray()
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
HighForestReader::HighForestReader(Int_t dataType, Int_t readMode, Int_t jetType, Int_t jetAxis, Bool_t matchJets, Bool_t readTrackTree) :
  ForestReader(dataType,readMode,jetType,jetAxis,matchJets,readTrackTree),
  fHeavyIonTree(0),
  fJetTree(0),
  fHltTree(0),
  fSkimTree(0),
  fTrackTree(0),
  fnJetsBranch(0),
  fnTracksBranch(0),
  fTrackAlgorithmBranch(0),
  fTrackOriginalAlgorithmBranch(0),
  fTrackMVABranch(0),
  fJetPtArray(),
  fJetPhiArray(),
  fJetEtaArray(),
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
  fTrackAlgorithmArray(),
  fTrackOriginalAlgorithmArray(),
  fTrackMVAArray()
{
  // Custom constructor
  
}

/*
 * Copy constructor
 */
HighForestReader::HighForestReader(const HighForestReader& in) :
  ForestReader(in),
  fHeavyIonTree(in.fHeavyIonTree),
  fJetTree(in.fJetTree),
  fHltTree(in.fHltTree),
  fSkimTree(in.fSkimTree),
  fTrackTree(in.fTrackTree),
  fnJetsBranch(in.fnJetsBranch),
  fnTracksBranch(in.fnTracksBranch),
  fTrackAlgorithmBranch(in.fTrackAlgorithmBranch),
  fTrackOriginalAlgorithmBranch(in.fTrackOriginalAlgorithmBranch),
  fTrackMVABranch(in.fTrackMVABranch)
{
  // Copy constructor
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetPtArray[i] = in.fJetPtArray[i];
    fJetPhiArray[i] = in.fJetPhiArray[i];
    fJetEtaArray[i] = in.fJetEtaArray[i];
    fJetRawPtArray[i] = in.fJetRawPtArray[i];
    fJetMaxTrackPtArray[i] = in.fJetMaxTrackPtArray[i];
  }
  
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
    fTrackAlgorithmArray[i] = in.fTrackAlgorithmArray[i];
    fTrackOriginalAlgorithmArray[i] = in.fTrackOriginalAlgorithmArray[i];
    fTrackMVAArray[i] = in.fTrackMVAArray[i];
  }
}

/*
 * Assignment operator
 */
HighForestReader& HighForestReader::operator=(const HighForestReader& in){
  // Assignment operator
  
  if (&in==this) return *this;
  
  ForestReader::operator=(in);
  
  fHeavyIonTree = in.fHeavyIonTree;
  fJetTree = in.fJetTree;
  fHltTree = in.fHltTree;
  fSkimTree = in.fSkimTree;
  fTrackTree = in.fTrackTree;
  fnJetsBranch = in.fnJetsBranch;
  fnTracksBranch = in.fnTracksBranch;
  fTrackAlgorithmBranch = in.fTrackAlgorithmBranch;
  fTrackOriginalAlgorithmBranch = in.fTrackOriginalAlgorithmBranch;
  fTrackMVABranch = in.fTrackMVABranch;
  
  for(Int_t i = 0; i < fnMaxJet; i++){
    fJetPtArray[i] = in.fJetPtArray[i];
    fJetPhiArray[i] = in.fJetPhiArray[i];
    fJetEtaArray[i] = in.fJetEtaArray[i];
    fJetRawPtArray[i] = in.fJetRawPtArray[i];
    fJetMaxTrackPtArray[i] = in.fJetMaxTrackPtArray[i];
  }
  
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
    fTrackAlgorithmArray[i] = in.fTrackAlgorithmArray[i];
    fTrackOriginalAlgorithmArray[i] = in.fTrackOriginalAlgorithmArray[i];
    fTrackMVAArray[i] = in.fTrackMVAArray[i];
  }
  
  return *this;
}

/*
 * Destructor
 */
HighForestReader::~HighForestReader(){
  // destructor
}

/*
 * Initialization, meaning that the branches are connected to the tree
 */
void HighForestReader::Initialize(){
  
  // Connect the branches of the heavy ion tree
  fHeavyIonTree->SetBranchStatus("*",0);
  fHeavyIonTree->SetBranchStatus("vz",1);
  fHeavyIonTree->SetBranchAddress("vz",&fVertexZ,&fHiVzBranch);
  fHeavyIonTree->SetBranchStatus("hiBin",1);
  fHeavyIonTree->SetBranchAddress("hiBin",&fHiBin,&fHiBinBranch);
  if(fDataType == kPpMC || fDataType == kPbPbMC || fDataType == kLocalTest){
    fHeavyIonTree->SetBranchStatus("pthat",1);
    fHeavyIonTree->SetBranchAddress("pthat",&fPtHat,&fPtHatBranch); // pT hat only for MC
  } else {
    fPtHat = 0; // We do not have pT hat information for real data
  }
  // Event weight for 2018 MC
  if((fDataType == kPbPbMC || fDataType == kPpMC) && fReadMode > 2000){
    fHeavyIonTree->SetBranchStatus("weight",1);
    fHeavyIonTree->SetBranchAddress("weight",&fEventWeight,&fEventWeightBranch);
  } else {
    fEventWeight = 1;
  }
  
  // Connect the branches to the jet tree
  const char *jetAxis[3] = {"jt","jt","WTA"};
  const char *genJetAxis[3] = {"","","WTA"};
  char branchName[30];
  
  fJetTree->SetBranchStatus("*",0);
  fJetTree->SetBranchStatus("jtpt",1);
  fJetTree->SetBranchAddress("jtpt",&fJetPtArray,&fJetPtBranch);
  
  // If specified, select WTA axis for jet phi
  sprintf(branchName,"%sphi",jetAxis[fJetAxis]);
  fJetTree->SetBranchStatus(branchName,1);
  fJetTree->SetBranchAddress(branchName,&fJetPhiArray,&fJetPhiBranch);
  
  // If specified, select WTA axis for jet eta
  sprintf(branchName,"%seta",jetAxis[fJetAxis]);
  fJetTree->SetBranchStatus(branchName,1);
  fJetTree->SetBranchAddress(branchName,&fJetEtaArray,&fJetEtaBranch);
  
  fJetTree->SetBranchStatus("nref",1);
  fJetTree->SetBranchAddress("nref",&fnJets,&fnJetsBranch);
  fJetTree->SetBranchStatus("rawpt",1);
  fJetTree->SetBranchAddress("rawpt",&fJetRawPtArray,&fJetRawPtBranch);
  fJetTree->SetBranchStatus("trackMax",1);
  fJetTree->SetBranchAddress("trackMax",&fJetMaxTrackPtArray,&fJetMaxTrackPtBranch);
  
  // If we are looking at Monte Carlo, connect the reference pT and parton arrays
  if(fDataType > kPbPb){
    fJetTree->SetBranchStatus("refpt",1);
    fJetTree->SetBranchAddress("refpt",&fJetRefPtArray,&fJetRefPtBranch);
    fJetTree->SetBranchStatus("refparton_flavor",1);
    fJetTree->SetBranchAddress("refparton_flavor",&fJetRefFlavorArray,&fJetRefFlavorBranch);
    fJetTree->SetBranchStatus("genpt",1);
    fJetTree->SetBranchAddress("genpt",&fMatchedJetPtArray,&fJetMatchedPtBranch);
    
    // If specified, select WTA axis for jet phi
    sprintf(branchName,"%sgenphi",genJetAxis[fJetAxis]);
    fJetTree->SetBranchStatus(branchName,1);
    fJetTree->SetBranchAddress(branchName,&fMatchedJetPhiArray,&fJetMatchedPhiBranch);
    
    // If specified, select WTA axis for jet eta
    sprintf(branchName,"%sgeneta",genJetAxis[fJetAxis]);
    fJetTree->SetBranchStatus(branchName,1);
    fJetTree->SetBranchAddress(branchName,&fMatchedJetEtaArray,&fJetMatchedEtaBranch);
    
    fJetTree->SetBranchStatus("ngen",1);
    fJetTree->SetBranchAddress("ngen",&fnMatchedJets,&fnMatchedJetsBranch);
  }
  
  // Event selection summary
  //
  //         tree                      branch                         What it is
  //  hltanalysis/HltTree   HLT_HIPuAK4CaloJet100_Eta5p1_v1      Event selection for PbPb
  //  hltanalysis/HltTree      HLT_AK4CaloJet80_Eta5p1_v1         Event selection for pp
  // skimanalysis/HltTree         pprimaryVertexFilter           Event selection for PbPb
  // skimanalysis/HltTree    HBHENoiseFilterResultRun2Loose   Event selection for pp and PbPb
  // skimanalysis/HltTree         pPAprimaryVertexFilter          Event selection for pp
  // skimanalysis/HltTree           pBeamScrapingFilter           Event selection for pp
  
  // Connect the branches to the HLT tree
  fHltTree->SetBranchStatus("*",0);
  if(fDataType == kPp || (fDataType == kPpMC && fReadMode == 2019)){ // pp data
    
    if(fReadMode > 2000){
      fHltTree->SetBranchStatus("HLT_HIAK4CaloJet80_v1",1); // 2017 syntax
      fHltTree->SetBranchAddress("HLT_HIAK4CaloJet80_v1",&fCaloJetFilterBit,&fCaloJetFilterBranch);
    } else {
      fHltTree->SetBranchStatus("HLT_AK4CaloJet80_Eta5p1_v1",1); // 2015 syntax
      fHltTree->SetBranchAddress("HLT_AK4CaloJet80_Eta5p1_v1",&fCaloJetFilterBit,&fCaloJetFilterBranch);
      
    }
    
  }else if ((fDataType == kPbPb ||fDataType == kPbPbMC) && fReadMode == 2019){ // PbPb data
    
    // Trigger branch naming is different for 2018 and 2015 forests
    if(fReadMode > 2000){
      
      // Calo jet 80 trigger
      //fHltTree->SetBranchStatus("HLT_HIPuAK4CaloJet80Eta5p1_v1",1); // 2018 syntax
      //fHltTree->SetBranchAddress("HLT_HIPuAK4CaloJet80Eta5p1_v1",&fCaloJetFilterBit,&fCaloJetFilterBranch); // 2018 syntax
      
      // Calo jet 100 trigger
      fHltTree->SetBranchStatus("HLT_HIPuAK4CaloJet100Eta5p1_v1",1); // 2018 syntax
      fHltTree->SetBranchAddress("HLT_HIPuAK4CaloJet100Eta5p1_v1",&fCaloJetFilterBit,&fCaloJetFilterBranch); // 2018 syntax
      
      // PF jet 80 trigger
      //fHltTree->SetBranchStatus("HLT_HICsAK4PFJet80Eta1p5_v1",1); // 2018 syntax
      //fHltTree->SetBranchAddress("HLT_HICsAK4PFJet80Eta1p5_v1",&fCaloJetFilterBit,&fCaloJetFilterBranch); // 2018 syntax
      
    } else {
      fHltTree->SetBranchStatus("HLT_HIPuAK4CaloJet100_Eta5p1_v1",1); // 2015 syntax
      fHltTree->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_v1",&fCaloJetFilterBit,&fCaloJetFilterBranch); // 2015 syntax
    }
  } else { // Local test or MC
    fCaloJetFilterBit = 1;  // No filter for local test or MC forests
  }
  fCaloJetFilterBitPrescale = 1; // Set the prescaled filter bit to 1. Only relevant for minimum bias PbPb (data skim)
  
  // Connect the branches to the skim tree (different for pp and PbPb data and Monte Carlo)
  fSkimTree->SetBranchStatus("*",0);
  // pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter
  if(fDataType == kPp || fDataType == kPpMC){ // pp data or MC
    fSkimTree->SetBranchStatus("pPAprimaryVertexFilter",1);
    fSkimTree->SetBranchAddress("pPAprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    fSkimTree->SetBranchStatus("pBeamScrapingFilter",1);
    fSkimTree->SetBranchAddress("pBeamScrapingFilter",&fBeamScrapingFilterBit,&fBeamScrapingBranch);
    fSkimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
    fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
    fHfCoincidenceFilterBit = 1; // No HF energy coincidence requirement for pp
    fClusterCompatibilityFilterBit = 1; // No cluster compatibility requirement for pp
  } else if (fDataType == kPbPb || fDataType == kPbPbMC){ // PbPb data or MC
    
    // Primary vertex has at least two tracks, is within 25 cm in z-rirection and within 2 cm in xy-direction
    fSkimTree->SetBranchStatus("pprimaryVertexFilter",1);
    fSkimTree->SetBranchAddress("pprimaryVertexFilter",&fPrimaryVertexFilterBit,&fPrimaryVertexBranch);
    
    // Cut on noise on HCAL
    fSkimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
    fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&fHBHENoiseFilterBit,&fHBHENoiseBranch);
    
    // Have at least two HF towers on each side of the detector with an energy deposit of 4 GeV
    fSkimTree->SetBranchStatus("phfCoincFilter2Th4",1);
    fSkimTree->SetBranchAddress("phfCoincFilter2Th4", &fHfCoincidenceFilterBit, &fHfCoincidenceBranch); // 2018 syntax
    
    // Calculated from pixel clusters. Ensures that measured and predicted primary vertices are compatible
    fSkimTree->SetBranchStatus("pclusterCompatibilityFilter",1);
    fSkimTree->SetBranchAddress("pclusterCompatibilityFilter",&fClusterCompatibilityFilterBit,&fClusterCompatibilityBranch);
    
    fBeamScrapingFilterBit = 1;  // No beam scraping filter for PbPb
  } else { // Local test
    fPrimaryVertexFilterBit = 1;
    fBeamScrapingFilterBit = 1;
    fHBHENoiseFilterBit = 1;
    fHfCoincidenceFilterBit = 1;
    fClusterCompatibilityFilterBit = 1;     
  }
  
  // Connect the branches to the track tree
  if(fReadTrackTree){
    fTrackTree->SetBranchStatus("*",0);
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
    
    // Additional information needed for 2018 track cuts
    fTrackTree->SetBranchStatus("trkAlgo",1);
    fTrackTree->SetBranchAddress("trkAlgo",&fTrackAlgorithmArray,&fTrackAlgorithmBranch);
    fTrackTree->SetBranchStatus("trkOriginalAlgo",1);
    fTrackTree->SetBranchAddress("trkOriginalAlgo",&fTrackOriginalAlgorithmArray,&fTrackOriginalAlgorithmBranch);
    
    // Track MVA only in 2018 PbPb trees
    if(fReadMode > 2000 && (fDataType == kPbPb || fDataType == kPbPbMC)){
      fTrackTree->SetBranchStatus("trkMVA",1);
      fTrackTree->SetBranchAddress("trkMVA",&fTrackMVAArray,&fTrackMVABranch);
    }
  } // Reading track trees
  
}

/*
 * Connect a new tree to the reader
 */
void HighForestReader::ReadForestFromFile(TFile *inputFile){
  
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
    treeName[1] = "akCs4PFJetAnalyzer/t";       // Tree for csPF jets
    treeName[2] = "akPu4PFJetAnalyzer/t";       // Tree for puPF jets
    treeName[3] = "akFlowPuCs4PFJetAnalyzer/t"; // Tree for flow subtracted csPF jets
  } else if (fDataType == kLocalTest){
    treeName[0] = "ak4PFJetAnalyzer/t";  // Only PF jets in local test file
    treeName[1] = "ak4PFJetAnalyzer/t";  // Only PF jets in local test file
  }
  fJetTree = (TTree*)inputFile->Get(treeName[fJetType]);
  if(fReadTrackTree) fTrackTree = (TTree*)inputFile->Get("ppTrack/trackTree");
  
  Initialize();
}

/*
 * Connect a new tree to the reader
 */
void HighForestReader::ReadForestFromFileList(std::vector<TString> fileList){
  TFile *inputFile = TFile::Open(fileList.at(0));
  ReadForestFromFile(inputFile);
}

/*
 * Burn the current forest.
 */
void HighForestReader::BurnForest(){
  fHeavyIonTree->Delete();
  fHltTree->Delete();
  fSkimTree->Delete();
  fJetTree->Delete();
  fTrackTree->Delete();
}

/*
 * Load an event to memory
 */
void HighForestReader::GetEvent(Int_t nEvent){
  fHeavyIonTree->GetEntry(nEvent);
  fJetTree->GetEntry(nEvent);
  fHltTree->GetEntry(nEvent);
  fSkimTree->GetEntry(nEvent);
  if(fReadTrackTree) fTrackTree->GetEntry(nEvent);
}

// Getter for jet pT
Float_t HighForestReader::GetJetPt(Int_t iJet) const{
  
  // TODO: For 2018, JEC is applied on the fly, so raw pT should be returned here
  // If this is fixed in the forest, can return the regular pT instead
  if(fReadMode > 2000) return fJetRawPtArray[iJet];
  
  return fJetPtArray[iJet];
}

// Getter for jet phi
Float_t HighForestReader::GetJetPhi(Int_t iJet) const{
  return fJetPhiArray[iJet];
}

// Getter for jet eta
Float_t HighForestReader::GetJetEta(Int_t iJet) const{
  return fJetEtaArray[iJet];
}

// Getter for jet raw pT
Float_t HighForestReader::GetJetRawPt(Int_t iJet) const{
  return fJetRawPtArray[iJet];
}

// Getter for maximum track pT inside a jet
Float_t HighForestReader::GetJetMaxTrackPt(Int_t iJet) const{
  return fJetMaxTrackPtArray[iJet];
}

// Getter for track pT
Float_t HighForestReader::GetTrackPt(Int_t iTrack) const{
  return fTrackPtArray[iTrack];
}

// Getter for track pT error
Float_t HighForestReader::GetTrackPtError(Int_t iTrack) const{
  return fTrackPtErrorArray[iTrack];
}

// Getter for track phi
Float_t HighForestReader::GetTrackPhi(Int_t iTrack) const{
  return fTrackPhiArray[iTrack];
}

// Getter for track eta
Float_t HighForestReader::GetTrackEta(Int_t iTrack) const{
  return fTrackEtaArray[iTrack];
}

// Getter for high purity of the track
Bool_t HighForestReader::GetTrackHighPurity(Int_t iTrack) const{
  return fHighPurityTrackArray[iTrack];
}

// Getter for track distance from primary vertex in z-direction
Float_t HighForestReader::GetTrackVertexDistanceZ(Int_t iTrack) const{
  return fTrackVertexDistanceZArray[iTrack];
}

// Getter for error of track distance from primary vertex in z-direction
Float_t HighForestReader::GetTrackVertexDistanceZError(Int_t iTrack) const{
  return fTrackVertexDistanceZErrorArray[iTrack];
}

// Getter for track distance from primary vertex in xy-direction
Float_t HighForestReader::GetTrackVertexDistanceXY(Int_t iTrack) const{
  return fTrackVertexDistanceXYArray[iTrack];
}

// Getter for error of track distance from primary vertex in xy-direction
Float_t HighForestReader::GetTrackVertexDistanceXYError(Int_t iTrack) const{
  return fTrackVertexDistanceXYErrorArray[iTrack];
}

// Getter for track chi2 value from reconstruction fit
Float_t HighForestReader::GetTrackChi2(Int_t iTrack) const{
  return fTrackChi2Array[iTrack];
}

// Getter for number of degrees of freedom in reconstruction fit
Int_t HighForestReader::GetNTrackDegreesOfFreedom(Int_t iTrack) const{
  return fnTrackDegreesOfFreedomArray[iTrack];
}

// Getter for number of hits in tracker layers
Int_t HighForestReader::GetNHitsTrackerLayer(Int_t iTrack) const{
  return fnHitsTrackerLayerArray[iTrack];
}

// Getter for number of hits for the track
Int_t HighForestReader::GetNHitsTrack(Int_t iTrack) const{
  return fnHitsTrackArray[iTrack];
}

// Getter for track energy in ECal
Float_t HighForestReader::GetTrackEnergyEcal(Int_t iTrack) const{
  return fTrackEnergyEcalArray[iTrack];
}

// Getter for track energy in HCal
Float_t HighForestReader::GetTrackEnergyHcal(Int_t iTrack) const{
  return fTrackEnergyHcalArray[iTrack];
}

// Getter for track charge. Charge is only relevant for generator level tracks.
Int_t HighForestReader::GetTrackCharge(Int_t iTrack) const{
  return 1;
}

// Getter for track subevent index. Relevant only for generator level tracks.
Int_t HighForestReader::GetTrackSubevent(Int_t iTrack) const{
  return -1;
}

// Getter for track MC status. Relevant only for generator level tracks.
Int_t HighForestReader::GetTrackMCStatus(Int_t iTrack) const{
  return 1;
}

// Getter for track algorithm
Int_t HighForestReader::GetTrackAlgorithm(Int_t iTrack) const{
  return fTrackAlgorithmArray[iTrack];
}

// Getter for track original algorithm
Int_t HighForestReader::GetTrackOriginalAlgorithm(Int_t iTrack) const{
  return fTrackOriginalAlgorithmArray[iTrack];
}

// Getter for track MVA
Float_t HighForestReader::GetTrackMVA(Int_t iTrack) const{
  return fTrackMVAArray[iTrack];
}

// Check if generator level jet has a matching reconstructed jet
Bool_t HighForestReader::HasMatchingJet(Int_t iJet) const{
  
  // If we are not matching jets or are considering real data, everything is fine
  if(!fMatchJets || fDataType <= kPbPb) return true;
  
  // For each reconstructed jet there is a reference pT, which tells the the pT of a matched generator level jet
  // If this number is -999, it means that there are no generator level jets matching the reconstructed jet
  if(fJetRefPtArray[iJet] < 0) return false;
  return true;
}

// Getter for matched generator level jet pT
Float_t HighForestReader::GetMatchedPt(Int_t iJet) const{
  
  // If we are not matching jets or are considering real data, this value has no meaning
  if(!fMatchJets || fDataType <= kPbPb) return 0;
  
  // Return matched gen pT
  return fJetRefPtArray[iJet];
}

// Parton flavor for the parton initiating the jet
Int_t HighForestReader::GetPartonFlavor(Int_t iJet) const{
  
  // If we are considering real data, this value has no meaning
  if(fDataType <= kPbPb) return 0;
  
  // Return the matching parton flavor
  return fJetRefFlavorArray[iJet];
}

// Get the matching generator level jet index for the given reconstructed jet
Int_t HighForestReader::GetMatchingIndex(Int_t iJet) const{
  
  // Ref pT array has pT for all the generator level jets that are matched with reconstructed jets
  // If our generator level pT is found from this array, it has a matching reconstructed jet
  Double_t jetPt = GetMatchedPt(iJet);
  Int_t matchedIndex = -1;
  for(Int_t iRef = 0; iRef < fnMatchedJets; iRef++){
    if(TMath::Abs(jetPt - fMatchedJetPtArray[iRef]) < 0.001){
      matchedIndex = iRef;
      break;
    }
  }
  
  return matchedIndex;
  
}

// Get the eta of the matched generator level jet
Float_t HighForestReader::GetMatchedEta(Int_t iJet) const{
  
  // If not matching jets, just return something because this has no meaning
  if(!fMatchJets) return 0;
  
  // Find the index of the matching reconstructed jet
  Int_t matchedIndex = GetMatchingIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchedIndex == -1) return -999;
  
  // Return the pT of the matching reconstructed jet
  return fMatchedJetEtaArray[matchedIndex];
  
}

// Get the pT of the matched reconstructed jet
Float_t HighForestReader::GetMatchedPhi(Int_t iJet) const{
  
  // If not matching jets, just return something because this has no meaning
  if(!fMatchJets) return 0;
  
  // Find the index of the matching reconstructed jet
  Int_t matchedIndex = GetMatchingIndex(iJet);
  
  // If we did not find macth, something went wrong. Return -999
  if(matchedIndex == -1) return -999;
  
  // Return the pT of the matching reconstructed jet
  return fMatchedJetPhiArray[matchedIndex];
  
}
