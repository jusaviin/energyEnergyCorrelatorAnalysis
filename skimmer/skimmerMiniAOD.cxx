// C++ includes
#include <iostream>   // Input/output stream. Needed for cout.
#include <fstream>    // File stream for intup/output to/from files
#include <stdlib.h>   // Standard utility libraries
#include <assert.h>   // Standard c++ debugging tool. Terminates the program if expression given evaluates to 0.
#include <vector>     // C++ vector class
#include <sstream>    // Libraries for checking boolean input
#include <string>     // Libraries for checking boolean input
#include <iomanip>    // Libraries for checking boolean input
#include <algorithm>  // Libraries for checking boolean input
#include <cctype>     // Libraries for checking boolean input

// Root includes
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TChain.h"
#include "TString.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"

using namespace std;

/*
 * File list reader
 *
 *  Arguments:
 *    std::vector<TString> &fileNameVector = Vector filled with filenames found in the file
 *    TString fileNameFile = Text file containing one analysis file name in each line
 *    int debug = Level of debug messages shown
 *    int locationIndex = Where to find analysis files: 0 = Purdue EOS, 1 = CERN EOS, 2 = Vanderbilt T2, 3 = Use xrootd to find the data
 *    bool runLocal = True: Local run mode. False: Crab run mode
 */
void ReadFileList(std::vector<TString> &fileNameVector, TString fileNameFile, int debug, int locationIndex, bool runLocal)
{
  
  // Possible location for the input files
  const char *fileLocation[] = {"root://xrootd.rcac.purdue.edu/", "root://eoscms.cern.ch/", "root://xrootd-vanderbilt.sites.opensciencegrid.org/", "root://cmsxrootd.fnal.gov/"};
  
  // Set up the file names file for reading
  ifstream file_stream(fileNameFile);
  std::string line;
  fileNameVector.clear();
  if( debug > 0 ) std::cout << "Open file " << fileNameFile.Data() << " to extract files to run over" << std::endl;
  
  // Open the file names file for reading
  if( file_stream.is_open() ) {
    if( debug > 0) std::cout << "Opened " << fileNameFile.Data() << " for reading" << std::endl;
    int lineNumber = 0;
    
    // Loop over the lines in the file
    while( !file_stream.eof() ) {
      getline(file_stream, line);
      if( debug > 0) std::cout << lineNumber << ": " << line << std::endl;
      TString lineString(line);
      
      // Put all non-empty lines to file names vector
      if( lineString.CompareTo("", TString::kExact) != 0 ) {
        
        if(runLocal){
          // For local running, it is assumed that the file name is directly the centents of the line
          fileNameVector.push_back(lineString);
          
        } else {
          // For crab running, the line will have format ["file1", "file2", ... , "fileN"]
          TObjArray *fileNameArray = lineString.Tokenize(" ");  // Tokenize the string from every ' ' character
          int numberOfFiles = fileNameArray->GetEntries();
          TObjString *currentFileNameObject;
          TString currentFileName;
          for(int i = 0; i < numberOfFiles; i++){   // Loop over all the files in the array
            currentFileNameObject = (TObjString *)fileNameArray->At(i);
            currentFileName = currentFileNameObject->String();
            
            // Strip unwanted characters
            currentFileName.Remove(TString::kBoth,'['); // Remove possible parantheses
            currentFileName.Remove(TString::kBoth,']'); // Remove possible parantheses
            currentFileName.Remove(TString::kBoth,','); // Remove commas
            currentFileName.Remove(TString::kBoth,'"'); // Remove quotation marks
            
            // After stripping characters not belonging to the file name, we can add the file to list
            currentFileName.Prepend(fileLocation[locationIndex]);  // If not running locally, we need to give xrootd path
            fileNameVector.push_back(currentFileName);
          }
        }
        
      } // Empty line if
      
      
      lineNumber++;
    } // Loop over lines in the file
    
    // If cannot read the file, give error and end program
  } else {
    std::cout << "Error, could not open " << fileNameFile.Data() << " for reading" << std::endl;
    assert(0);
  }
}

/*
 *  Convert string to boolean value
 */
bool checkBool(string str) {
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);
  std::istringstream is(str);
  bool b;
  is >> std::boolalpha >> b;
  return b;
}

/*
 * Main function of the skimmer
 */
int main(int argc,char *argv[]){
  
  //==== Read arguments =====
  if ( argc<5 ) {
    cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
    cout<<"+ Usage of the macro: " << endl;
    cout<<"+  "<<argv[0]<<" [fileNameFile] [isMC] [outputFileName] [fileLocation] <runLocal>"<<endl;
    cout<<"+  fileNameFile: Text file containing the list of files used in the analysis. For crab analysis a job id should be given here." <<endl;
    cout<<"+  isMC: 1 for MC, 0 for data." <<endl;
    cout<<"+  outputFileName: .root file to which the histograms are written." <<endl;
    cout<<"+  fileLocation: Where to find analysis files: 0 = Purdue EOS, 1 = CERN EOS, 2 = Vanderbilt T2, 3 = Use xrootd to find the data." << endl;
    cout<<"+  runLocal: True: Search input files from local machine. False (default): Search input files from grid with xrootd." << endl;
    cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
    cout << endl << endl;
    exit(1);
  }
  
  // First, check if we are supposed to run locally or on crab
  bool runLocal = false;
  if(argc >= 6) runLocal = checkBool(argv[5]);
  
  // Find the file list name depending on if we run locally or on crab
  TString fileNameFile;
  if(runLocal){
    fileNameFile = argv[1];
  } else {
    fileNameFile = Form("job_input_file_list_%d.txt",atoi(argv[1]));
  }
  
  // Read the other command line arguments
  int isMC = atoi(argv[2]);
  TString outputFileName = argv[3];
  const int fileSearchIndex = atoi(argv[4]);
  
  // Read the file names used for skimming to a vector
  std::vector<TString> fileNameVector;
  fileNameVector.clear();
  ReadFileList(fileNameVector,fileNameFile,1,fileSearchIndex,runLocal);
  
  // Maximum size of arrays
  const Int_t nMaxJet = 250;        // Maximum number of jets in an event
  
  bool includePfCandidates = false;
  
  // Define trees to be read from the files
  const int nJetTrees = 2; // 2 For PbPb, 1 for pp
  TChain *hiForestInfoTree = new TChain("HiForestInfo/HiForest");
  TChain *heavyIonTree = new TChain("hiEvtAnalyzer/HiTree");
  TChain *hltTree = new TChain("hltanalysis/HltTree");
  TChain *skimTree = new TChain("skimanalysis/HltTree");
  TChain *jetTree[nJetTrees];
  jetTree[0] = new TChain("akCs4PFJetAnalyzer/t");       // For pp: ak4PFJetAnalyzer/t
  jetTree[1] = new TChain("akFlowPuCs4PFJetAnalyzer/t"); // For pp: remove
  TChain *trackTree = new TChain("PbPbTracks/trackTree");   // For pp: ppTracks/trackTree
  TChain *genTrackTree = new TChain("HiGenParticleAna/hi");
  TChain *particleFlowCandidateTree = new TChain("particleFlowAnalyser/pftree");
  
  // All the branches and leaves come in arrat of two, one for input and one for output
  
  // Branches for the HiForestInfo tree
  TBranch *infoBranch;
  TBranch *versionBranch;
  TBranch *globalTagBranch;
  
  // Leaves for the HiForestInfo tree
  std::string infoString;
  std::string versionString;
  std::string globalTagString;
  
  // Branches for heavy ion tree
  TBranch *runBranch;             // Branch for run
  TBranch *eventBranch;           // Branch for event
  TBranch *lumiBranch;            // Branch for lumi
  TBranch *hiVzBranch;            // Branch for vertex z-position
  TBranch *hiBinBranch;           // Branch for centrality
  TBranch *ptHatBranch;           // Branch for pT hat
  TBranch *eventWeightBranch;     // Branch for jet weight for 2018 MC
  
  // Leaves for heavy ion tree
  UInt_t run;           // Run number
  ULong64_t event;      // Event number
  UInt_t lumi;          // Luminosity block
  Float_t vertexZ;      // Vertex z-position
  Int_t hiBin;          // HiBin = Centrality percentile * 2
  Float_t ptHat;        // pT hat
  Float_t eventWeight;  // jet weight in the 2018 MC tree
  
  // Branches for HLT tree
  TBranch *caloJetFilterBranch80;         // Branch for calo jet 80 filter bit
  TBranch *caloJetFilterBranch100;        // Branch for calo jet 100 filter bit
  
  // Leaves for the HLT tree
  Int_t caloJetFilterBit80;          // Filter bit for calorimeter jets 80
  Int_t caloJetFilterBit100;         // Filter bit for calorimeter jets 100
  
  // Branches for skim tree
  TBranch *primaryVertexBranch;             // Branch for primary vertex filter bit
  // TBranch *beamScrapingBranch;              // Branch for beam scraping filter bit (only in pp)
  TBranch *hfCoincidenceBranch2Th4;         // Branch for energy recorded in at least 2 HF calorimeter towers
  TBranch *clusterCompatibilityBranch;      // Branch for cluster compatibility
  
  // Leaves for the skim tree
  Int_t primaryVertexFilterBit;             // Filter bit for primary vertex
  // Int_t beamScrapingFilterBit;              // Filter bit for beam scraping (only in pp)
  Int_t hfCoincidenceFilterBit2Th4;         // Filter bit for energy recorded in at least 2 HF calorimeter towers
  Int_t clusterCompatibilityFilterBit;      // Filter bit for cluster compatibility
  
  // Branches for jet tree
  TBranch *nJetsBranch[nJetTrees];            // Branch for number of jets in an event
  TBranch *nGenJetsBranch[nJetTrees];               // Branch for the number of generator level jets in an event
  TBranch *jetPtBranch[nJetTrees];                  // Branch for jet pT
  TBranch *jetPhiBranch[nJetTrees];                 // Branch for jet phi
  TBranch *jetPhiBranchWTA[nJetTrees];              // Branch for jet phi with WTA axis
  TBranch *jetEtaBranch[nJetTrees];                 // Branch for jet eta
  TBranch *jetEtaBranchWTA[nJetTrees];              // Branch for jet eta with WTA axis
  TBranch *jetRawPtBranch[nJetTrees];               // Branch for raw jet pT
  TBranch *jetMaxTrackPtBranch[nJetTrees];          // Maximum pT for a track inside a jet
  TBranch *jetRefPtBranch[nJetTrees];               // Branch for reference generator level pT for a reconstructed jet
  TBranch *jetRefFlavorBranch[nJetTrees];           // Branch for flavor for the parton initiating the jet
  TBranch *jetRefFlavorForBBranch[nJetTrees];       // Branch for flavor for the parton initiating the jet
  TBranch *genJetPtBranch[nJetTrees];               // Branch for the generator level jet pT
  TBranch *genJetEtaBranch[nJetTrees];              // Branch for the generetor level jet eta
  TBranch *genJetEtaBranchWTA[nJetTrees];           // Branch for the generetor level jet eta with WTA axis
  TBranch *genJetPhiBranch[nJetTrees];              // Branch for the generator level jet phi
  TBranch *genJetPhiBranchWTA[nJetTrees];           // Branch for the generator level jet phi with WTA axis
  
  // Leaves for jet tree
  Int_t nJets[nJetTrees];                                        // number of jets in an event
  Int_t nGenJets[nJetTrees];                                     // number of generator level jets in an event
  Float_t jetPtArray[nJetTrees][nMaxJet] = {{0}};                // pT:s of all the jets in an event
  Float_t jetPhiArray[nJetTrees][nMaxJet] = {{0}};               // phis of all the jets in an event
  Float_t jetPhiArrayWTA[nJetTrees][nMaxJet] = {{0}};            // phis of all the jets in an event  with WTA axis
  Float_t jetEtaArray[nJetTrees][nMaxJet] = {{0}};               // etas of all the jets in an event
  Float_t jetEtaArrayWTA[nJetTrees][nMaxJet] = {{0}};            // etas of all the jets in an event  with WTA axis
  Float_t jetRawPtArray[nJetTrees][nMaxJet] = {{0}};             // raw jet pT for all the jets in an event
  Float_t jetMaxTrackPtArray[nJetTrees][nMaxJet] = {{0}};        // maximum track pT inside a jet for all the jets in an event
  Float_t jetRefPtArray[nJetTrees][nMaxJet] = {{0}};             // reference generator level pT for a reconstructed jet
  Int_t jetRefFlavorArray[nJetTrees][nMaxJet] = {{0}};           // flavor for initiating parton for the reference gen jet
  Int_t jetRefFlavorForBArray[nJetTrees][nMaxJet] = {{0}};       // heavy flavor for initiating parton for the reference gen jet
  Float_t genJetPtArray[nJetTrees][nMaxJet] = {{0}};             // pT:s of all the generator level jets in an event
  Float_t genJetPhiArray[nJetTrees][nMaxJet] = {{0}};            // phis of all the generator level jets in an event
  Float_t genJetPhiArrayWTA[nJetTrees][nMaxJet] = {{0}};         // phis of all the generator level jets in an event with WTA axis
  Float_t genJetEtaArray[nJetTrees][nMaxJet] = {{0}};            // etas of all the generator level jets in an event
  Float_t genJetEtaArrayWTA[nJetTrees][nMaxJet] = {{0}};         // etas of all the generator level jets in an event with WTA axis
  
  // Branches for track tree
  TBranch *nTracksBranch;                    // Branch for number of tracks
  TBranch *trackPtBranch;                    // Branch for track pT
  TBranch *trackPtErrorBranch;               // Branch for track pT error
  TBranch *trackPhiBranch;                   // Branch for track phi
  TBranch *trackEtaBranch;                   // Branch for track eta
  TBranch *trackHighPurityBranch;            // Branch for high purity of the track
  TBranch *trackVertexDistanceZBranch;       // Branch for track distance from primary vertex in z-direction
  TBranch *trackVertexDistanceZErrorBranch;  // Branch for error for track distance from primary vertex in z-direction
  TBranch *trackVertexDistanceXYBranch;      // Branch for track distance from primary vertex in xy-direction
  TBranch *trackVertexDistanceXYErrorBranch; // Branch for error for track distance from primary vertex in xy-direction
  TBranch *trackNormalizedChi2Branch;        // Branch for track chi2 value from reconstruction fit
  TBranch *nHitsTrackerLayerBranch;          // Branch for number of hits in tracker layers
  TBranch *nHitsTrackBranch;                 // Branch for number of hits for the track
  TBranch *trackEnergyEcalBranch;            // Branch for track energy in ECal
  TBranch *trackEnergyHcalBranch;            // Branch for track energy in HCal
  TBranch *trackChargeBranch;                // Branch for track charge
  
  // Leaves for the track tree
  Int_t nTracks;                                  // Number of tracks
  vector<float> *trackPtArray;                    // Array for track pT:s
  vector<float> *trackPtErrorArray;               // Array for track pT errors
  vector<float> *trackPhiArray;                   // Array for track phis
  vector<float> *trackEtaArray;                   // Array for track etas
  vector<bool> *trackHighPurityArray;             // Array for the high purity of tracks
  vector<float> *trackVertexDistanceZArray;       // Array for track distance from primary vertex in z-direction
  vector<float> *trackVertexDistanceZErrorArray;   // Array for error for track distance from primary vertex in z-direction
  vector<float> *trackVertexDistanceXYArray;      // Array for track distance from primary vertex in xy-direction
  vector<float> *trackVertexDistanceXYErrorArray; // Array for error for track distance from primary vertex in xy-direction
  vector<float> *trackNormalizedChi2Array;        // Array for track chi2 value from reconstruction fit
  vector<char> *nHitsTrackerLayerArray;           // Array for number of hits in tracker layers
  vector<char> *nHitsTrackArray;                  // Array for number of hits for the track
  vector<float> *trackEnergyEcalArray;            // Array for track energy in ECal
  vector<float> *trackEnergyHcalArray;            // Array for track energy in HCal
  vector<char> *trackChargeArray;                 // Array for track charge
  
  // Branches for generator level track tree
  TBranch *genTrackPtBranch;         // Branch for generator level track pT:s
  TBranch *genTrackPhiBranch;        // Branch for generator level track phis
  TBranch *genTrackEtaBranch;        // Branch for generator level track etas
  TBranch *genTrackPdgBranch;        // Branch for generator level track PDG code
  TBranch *genTrackChargeBranch;     // Branch for generator level track charges
  TBranch *genTrackSubeventBranch;   // Branch for generator level track subevent indices (0 = PYTHIA, (>0) = HYDJET)
  
  // Leaves for generator level track tree
  vector<float> *genTrackPtArray;       // Array for generator level track pT:s
  vector<float> *genTrackPhiArray;      // Array for generator level track phis
  vector<float> *genTrackEtaArray;      // Array for generator level track etas
  vector<int> *genTrackPdgArray;        // Array for generator level track PDG code
  vector<int> *genTrackChargeArray;     // Array for generator level track charges
  vector<int> *genTrackSubeventArray;   // Array for generator level track subevent indices (0 = PYTHIA, (>0) = HYDJET)
  
  // Branches for particle flow candidate ID tree
  TBranch *particleFlowCandidateIdBranch;    // Branch for particle flow candidate ID
  TBranch *particleFlowCandidatePtBranch;    // Branch for particle flow candidate pT
  TBranch *particleFlowCandidatePhiBranch;   // Branch for particle flow candidate phi
  TBranch *particleFlowCandidateEtaBranch;   // Branch for particle flow candidate eta
  TBranch *particleFlowCandidateMassBranch;  // Branch for particle flow candidate mass
  
  // Leaves for particle flow candidate tree
  vector<int> *particleFlowCandidateIdVector;       // Vector for particle flow candidate ID:s
  vector<float> *particleFlowCandidatePtVector;     // Vector for particle flow candidate pT:s
  vector<float> *particleFlowCandidatePhiVector;    // Vector for particle flow candidate phis
  vector<float> *particleFlowCandidateEtaVector;    // Vector for particle flow candidate etas
  vector<float> *particleFlowCandidateMassVector;   // Vector for particle flow candidate masses
  
  // Add all the files to the chain
  for(std::vector<TString>::iterator listIterator = fileNameVector.begin(); listIterator != fileNameVector.end(); listIterator++){
    
    cout << "Adding file " << *listIterator << " to the chains" << endl;
    
    hiForestInfoTree->Add(*listIterator);
    heavyIonTree->Add(*listIterator);
    hltTree->Add(*listIterator);
    skimTree->Add(*listIterator);
    trackTree->Add(*listIterator);
    genTrackTree->Add(*listIterator);
    if(includePfCandidates) particleFlowCandidateTree->Add(*listIterator);
    
    for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
      jetTree[iJetType]->Add(*listIterator);
    }
    
  }

  
  // ========================================== //
  // Read all the branches from the input trees //
  // ========================================== //
  
  // Connect the branches of the HiForestInfo tree
  hiForestInfoTree->SetBranchStatus("*",0);
  hiForestInfoTree->SetBranchStatus("info_0",1);
  hiForestInfoTree->SetBranchAddress("info_0",&infoString,&infoBranch);
  hiForestInfoTree->SetBranchStatus("HiForestVersion",1);
  hiForestInfoTree->SetBranchAddress("HiForestVersion",&versionString,&versionBranch);
  hiForestInfoTree->SetBranchStatus("GlobalTag",1);
  hiForestInfoTree->SetBranchAddress("GlobalTag",&globalTagString,&globalTagBranch);
  
  // Connect the branches of the heavy ion tree
  heavyIonTree->SetBranchStatus("*",0);
  heavyIonTree->SetBranchStatus("run",1);
  heavyIonTree->SetBranchAddress("run",&run,&runBranch);
  heavyIonTree->SetBranchStatus("evt",1);
  heavyIonTree->SetBranchAddress("evt",&event,&eventBranch);
  heavyIonTree->SetBranchStatus("lumi",1);
  heavyIonTree->SetBranchAddress("lumi",&lumi,&lumiBranch);
  heavyIonTree->SetBranchStatus("vz",1);
  heavyIonTree->SetBranchAddress("vz",&vertexZ,&hiVzBranch);
  heavyIonTree->SetBranchStatus("hiBin",1);
  heavyIonTree->SetBranchAddress("hiBin",&hiBin,&hiBinBranch);
  
  // ptHat and event weight only for MC
  if(isMC){
    heavyIonTree->SetBranchStatus("pthat",1);
    heavyIonTree->SetBranchAddress("pthat",&ptHat,&ptHatBranch);
    heavyIonTree->SetBranchStatus("weight",1);
    heavyIonTree->SetBranchAddress("weight",&eventWeight,&eventWeightBranch);
  }
  
  // Connect the branches to the HLT tree
  hltTree->SetBranchStatus("*",0);
  
  // PbPb syntax
  hltTree->SetBranchStatus("HLT_HIPuAK4CaloJet80Eta5p1_v1",1); // 2018 syntax
  hltTree->SetBranchAddress("HLT_HIPuAK4CaloJet80Eta5p1_v1",&caloJetFilterBit80,&caloJetFilterBranch80); // 2018 syntax
  hltTree->SetBranchStatus("HLT_HIPuAK4CaloJet100Eta5p1_v1",1); // 2018 syntax
  hltTree->SetBranchAddress("HLT_HIPuAK4CaloJet100Eta5p1_v1",&caloJetFilterBit100,&caloJetFilterBranch100); // 2018 syntax
  
  // pp syntax
  // hltTree->SetBranchStatus("HLT_HIAK4CaloJet80_v1",1); // 2017 syntax
  // hltTree->SetBranchAddress("HLT_HIAK4CaloJet80_v1",&fCaloJetFilterBit,&fCaloJetFilterBranch);
  
  // Connect the branches to the skim tree
  skimTree->SetBranchStatus("*",0);

  skimTree->SetBranchStatus("pprimaryVertexFilter",1);
  skimTree->SetBranchAddress("pprimaryVertexFilter",&primaryVertexFilterBit,&primaryVertexBranch);
  
  skimTree->SetBranchStatus("pphfCoincFilter2Th4",1);
  skimTree->SetBranchAddress("pphfCoincFilter2Th4", &hfCoincidenceFilterBit2Th4, &hfCoincidenceBranch2Th4);
  
  skimTree->SetBranchStatus("pclusterCompatibilityFilter",1);
  skimTree->SetBranchAddress("pclusterCompatibilityFilter",&clusterCompatibilityFilterBit,&clusterCompatibilityBranch);
  
  // Only in pp
  // skimTree->SetBranchStatus("pBeamScrapingFilter",1);
  // skimTree->SetBranchAddress("pBeamScrapingFilter",&beamScrapingFilterBit,&beamScrapingBranch);
  
  // pp syntax
  // skimTree->SetBranchStatus("pPAprimaryVertexFilter",1);
  // skimTree->SetBranchAddress("pPAprimaryVertexFilter",&primaryVertexFilterBit,&primaryVertexBranch);
  
  // Same branch names for all jet collections
  for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
    
    // Connect the branches to the jet tree
    jetTree[iJetType]->SetBranchStatus("*",0);
    jetTree[iJetType]->SetBranchStatus("jtpt",1);
    jetTree[iJetType]->SetBranchAddress("jtpt",&jetPtArray[iJetType],&jetPtBranch[iJetType]);
    
    // Jet eta with E-scheme and WTA axes
    jetTree[iJetType]->SetBranchStatus("jtphi",1);
    jetTree[iJetType]->SetBranchAddress("jtphi",&jetPhiArray[iJetType],&jetPhiBranch[iJetType]);
    jetTree[iJetType]->SetBranchStatus("WTAphi",1);
    jetTree[iJetType]->SetBranchAddress("WTAphi",&jetPhiArrayWTA[iJetType],&jetPhiBranchWTA[iJetType]);
    
    // Jet phi with E-scheme and WTA axes
    jetTree[iJetType]->SetBranchStatus("jteta",1);
    jetTree[iJetType]->SetBranchAddress("jteta",&jetEtaArray[iJetType],&jetEtaBranch[iJetType]);
    jetTree[iJetType]->SetBranchStatus("WTAeta",1);
    jetTree[iJetType]->SetBranchAddress("WTAeta",&jetEtaArrayWTA[iJetType],&jetEtaBranchWTA[iJetType]);
    
    jetTree[iJetType]->SetBranchStatus("nref",1);
    jetTree[iJetType]->SetBranchAddress("nref",&nJets[iJetType],&nJetsBranch[iJetType]);
    jetTree[iJetType]->SetBranchStatus("rawpt",1);
    jetTree[iJetType]->SetBranchAddress("rawpt",&jetRawPtArray[iJetType],&jetRawPtBranch[iJetType]);
    jetTree[iJetType]->SetBranchStatus("trackMax",1);
    jetTree[iJetType]->SetBranchAddress("trackMax",&jetMaxTrackPtArray[iJetType],&jetMaxTrackPtBranch[iJetType]);
    
    // If we are looking at Monte Carlo, connect the reference pT and parton arrays
    if(isMC){
      jetTree[iJetType]->SetBranchStatus("refpt",1);
      jetTree[iJetType]->SetBranchAddress("refpt",&jetRefPtArray[iJetType],&jetRefPtBranch[iJetType]);
      jetTree[iJetType]->SetBranchStatus("refparton_flavor",1);
      jetTree[iJetType]->SetBranchAddress("refparton_flavor",&jetRefFlavorArray[iJetType],&jetRefFlavorBranch[iJetType]);
      jetTree[iJetType]->SetBranchStatus("refparton_flavorForB",1);
      jetTree[iJetType]->SetBranchAddress("refparton_flavorForB", &jetRefFlavorForBArray[iJetType], &jetRefFlavorForBBranch[iJetType]);
      
      // Gen jet variables
      jetTree[iJetType]->SetBranchStatus("genpt",1);
      jetTree[iJetType]->SetBranchAddress("genpt",&genJetPtArray[iJetType],&genJetPtBranch[iJetType]);
      
      // Gen jet phi for e-scheme and WTA axes
      jetTree[iJetType]->SetBranchStatus("genphi",1);
      jetTree[iJetType]->SetBranchAddress("genphi",&genJetPhiArray[iJetType],&genJetPhiBranch[iJetType]);
      jetTree[iJetType]->SetBranchStatus("WTAgenphi",1);
      jetTree[iJetType]->SetBranchAddress("WTAgenphi",&genJetPhiArrayWTA[iJetType],&genJetPhiBranchWTA[iJetType]);
      
      // Gen jet eta for e-scheme and WTA axes
      jetTree[iJetType]->SetBranchStatus("geneta",1);
      jetTree[iJetType]->SetBranchAddress("geneta",&genJetEtaArray[iJetType],&genJetEtaBranch[iJetType]);
      jetTree[iJetType]->SetBranchStatus("WTAgeneta",1);
      jetTree[iJetType]->SetBranchAddress("WTAgeneta",&genJetEtaArrayWTA[iJetType],&genJetEtaBranchWTA[iJetType]);
      
      jetTree[iJetType]->SetBranchStatus("ngen",1);
      jetTree[iJetType]->SetBranchAddress("ngen",&nGenJets[iJetType],&nGenJetsBranch[iJetType]);
    }
    
  } // Loop over different jet collections
  
  // Connect the branches to the track tree
  trackTree->SetBranchStatus("*",0);
  
  trackTree->SetBranchStatus("trkPt",1);
  trackTree->SetBranchAddress("trkPt",&trackPtArray,&trackPtBranch);
  trackTree->SetBranchStatus("trkPtError",1);
  trackTree->SetBranchAddress("trkPtError",&trackPtErrorArray,&trackPtErrorBranch);
  trackTree->SetBranchStatus("trkPhi",1);
  trackTree->SetBranchAddress("trkPhi",&trackPhiArray,&trackPhiBranch);
  trackTree->SetBranchStatus("trkEta",1);
  trackTree->SetBranchAddress("trkEta",&trackEtaArray,&trackEtaBranch);
  trackTree->SetBranchStatus("nTrk",1);
  trackTree->SetBranchAddress("nTrk",&nTracks,&nTracksBranch);
  trackTree->SetBranchStatus("highPurity",1);
  trackTree->SetBranchAddress("highPurity",&trackHighPurityArray,&trackHighPurityBranch);
  trackTree->SetBranchStatus("trkDzFirstVtx",1);
  trackTree->SetBranchAddress("trkDzFirstVtx",&trackVertexDistanceZArray,&trackVertexDistanceZBranch);
  trackTree->SetBranchStatus("trkDzErrFirstVtx",1);
  trackTree->SetBranchAddress("trkDzErrFirstVtx",&trackVertexDistanceZErrorArray,&trackVertexDistanceZErrorBranch);
  trackTree->SetBranchStatus("trkDxyFirstVtx",1);
  trackTree->SetBranchAddress("trkDxyFirstVtx",&trackVertexDistanceXYArray,&trackVertexDistanceXYBranch);
  trackTree->SetBranchStatus("trkDxyErrFirstVtx",1);
  trackTree->SetBranchAddress("trkDxyErrFirstVtx",&trackVertexDistanceXYErrorArray,&trackVertexDistanceXYErrorBranch);
  trackTree->SetBranchStatus("trkNormChi2",1);
  trackTree->SetBranchAddress("trkNormChi2",&trackNormalizedChi2Array,&trackNormalizedChi2Branch);
  trackTree->SetBranchStatus("trkNLayers",1);
  trackTree->SetBranchAddress("trkNLayers",&nHitsTrackerLayerArray,&nHitsTrackerLayerBranch);
  trackTree->SetBranchStatus("trkNHits",1);
  trackTree->SetBranchAddress("trkNHits",&nHitsTrackArray,&nHitsTrackBranch);
  trackTree->SetBranchStatus("pfEcal",1);
  trackTree->SetBranchAddress("pfEcal",&trackEnergyEcalArray,&trackEnergyEcalBranch);
  trackTree->SetBranchStatus("pfHcal",1);
  trackTree->SetBranchAddress("pfHcal",&trackEnergyHcalArray,&trackEnergyHcalBranch);
  
  // Additional saved branches
  trackTree->SetBranchStatus("trkCharge",1);
  trackTree->SetBranchAddress("trkCharge",&trackChargeArray,&trackChargeBranch);
  
  // Generator level tracks only in Monte Carlo
  if(isMC){
    
    // Connect the branches to generator level track tree
    genTrackTree->SetBranchStatus("*",0);
    
    genTrackTree->SetBranchStatus("pt",1);
    genTrackTree->SetBranchAddress("pt",&genTrackPtArray,&genTrackPtBranch);
    genTrackTree->SetBranchStatus("phi",1);
    genTrackTree->SetBranchAddress("phi",&genTrackPhiArray,&genTrackPhiBranch);
    genTrackTree->SetBranchStatus("eta",1);
    genTrackTree->SetBranchAddress("eta",&genTrackEtaArray,&genTrackEtaBranch);
    genTrackTree->SetBranchStatus("pdg",1);
    genTrackTree->SetBranchAddress("pdg",&genTrackPdgArray,&genTrackPdgBranch);
    genTrackTree->SetBranchStatus("chg",1);
    genTrackTree->SetBranchAddress("chg",&genTrackChargeArray,&genTrackChargeBranch);
    genTrackTree->SetBranchStatus("sube",1);
    genTrackTree->SetBranchAddress("sube",&genTrackSubeventArray,&genTrackSubeventBranch);
    
  }
  
  // Connect the branches to the particle flow candidate tree if requested
  if(includePfCandidates){
    particleFlowCandidateTree->SetBranchStatus("*",0);
    
    particleFlowCandidateTree->SetBranchStatus("pfId",1);
    particleFlowCandidateTree->SetBranchAddress("pfId",&particleFlowCandidateIdVector,&particleFlowCandidateIdBranch);
    particleFlowCandidateTree->SetBranchStatus("pfPt",1);
    particleFlowCandidateTree->SetBranchAddress("pfPt",&particleFlowCandidatePtVector,&particleFlowCandidatePtBranch);
    particleFlowCandidateTree->SetBranchStatus("pfPhi",1);
    particleFlowCandidateTree->SetBranchAddress("pfPhi",&particleFlowCandidatePhiVector,&particleFlowCandidatePhiBranch);
    particleFlowCandidateTree->SetBranchStatus("pfEta",1);
    particleFlowCandidateTree->SetBranchAddress("pfEta",&particleFlowCandidateEtaVector,&particleFlowCandidateEtaBranch);
    particleFlowCandidateTree->SetBranchStatus("pfM",1);
    particleFlowCandidateTree->SetBranchAddress("pfM",&particleFlowCandidateMassVector,&particleFlowCandidateMassBranch);
  }
  
  
  // ========================================== //
  //           Define output trees
  // ========================================== //
  
  // Copy the HiForestInfo tree to the output
  TTree *hiForestInfoTreeOutput = new TTree("HiForest","");
  
  // Connect the branches of the HiForestInfo tree
  hiForestInfoTreeOutput->Branch("info_0",infoString.data(),"info/C");
  hiForestInfoTreeOutput->Branch("HiForestVersion",versionString.data(),"HiForestVersion/C");
  hiForestInfoTreeOutput->Branch("GlobalTag",globalTagString.data(),"GlobalTag/C");
  
  // Copy the heavy ion tree to the output
  TTree *heavyIonTreeOutput = new TTree("HiTree","");
  
  // Connect the branches of the heavy ion tree
  heavyIonTreeOutput->Branch("run",&run,"run/i");
  heavyIonTreeOutput->Branch("evt",&event,"evt/l");
  heavyIonTreeOutput->Branch("lumi",&lumi,"lumi/i");
  heavyIonTreeOutput->Branch("vz",&vertexZ,"vz/F");
  heavyIonTreeOutput->Branch("hiBin",&hiBin,"hiBin/I");
  
  // ptHat and event weight only for MC
  if(isMC){
    heavyIonTreeOutput->Branch("pthat",&ptHat,"pthat/F");
    heavyIonTreeOutput->Branch("weight",&eventWeight,"weight/F");
  }
  
  // Copy the HLT tree to the output
  TTree *hltTreeOutput = new TTree("HltTree","");
  
  // Connect the branches of the HLT tree
  hltTreeOutput->Branch("HLT_HIPuAK4CaloJet80Eta5p1_v1",&caloJetFilterBit80,"HLT_HIPuAK4CaloJet80Eta5p1_v1/I");
  hltTreeOutput->Branch("HLT_HIPuAK4CaloJet100Eta5p1_v1",&caloJetFilterBit100,"HLT_HIPuAK4CaloJet100Eta5p1_v1/I");
  
  // Copy the skim tree to the output
  TTree *skimTreeOutput = new TTree("HltTree","");
  
  skimTreeOutput->Branch("pprimaryVertexFilter",&primaryVertexFilterBit,"pprimaryVertexFilter/I");
  skimTreeOutput->Branch("pphfCoincFilter2Th4", &hfCoincidenceFilterBit2Th4, "pphfCoincFilter2Th4/I");
  skimTreeOutput->Branch("pclusterCompatibilityFilter",&clusterCompatibilityFilterBit,"pclusterCompatibilityFilter/I");
  // skimTreeOutput->Branch("pBeamScrapingFilter",&beamScrapingFilterBit,"pBeamScrapingFilter/I"); // Only in pp
       
  // Copy the jet trees to the output
  TTree *jetTreeOutput[nJetTrees];
  
  // Leaves for jet tree
  Int_t nJetsOutput[nJetTrees];                                        // number of jets in an event
  Int_t nGenJetsOutput[nJetTrees];                                     // number of generator level jets in an event
  Float_t jetPtArrayOutput[nJetTrees][nMaxJet] = {{0}};                // pT:s of all the jets in an event
  Float_t jetPhiArrayOutput[nJetTrees][nMaxJet] = {{0}};               // phis of all the jets in an event
  Float_t jetPhiArrayWTAOutput[nJetTrees][nMaxJet] = {{0}};            // phis of all the jets in an event  with WTA axis
  Float_t jetEtaArrayOutput[nJetTrees][nMaxJet] = {{0}};               // etas of all the jets in an event
  Float_t jetEtaArrayWTAOutput[nJetTrees][nMaxJet] = {{0}};            // etas of all the jets in an event  with WTA axis
  Float_t jetRawPtArrayOutput[nJetTrees][nMaxJet] = {{0}};             // raw jet pT for all the jets in an event
  Float_t jetMaxTrackPtArrayOutput[nJetTrees][nMaxJet] = {{0}};        // maximum track pT inside a jet for all the jets in an event
  Float_t jetRefPtArrayOutput[nJetTrees][nMaxJet] = {{0}};             // reference generator level pT for a reconstructed jet
  Int_t jetRefFlavorArrayOutput[nJetTrees][nMaxJet] = {{0}};           // flavor for initiating parton for the reference gen jet
  Int_t jetRefFlavorForBArrayOutput[nJetTrees][nMaxJet] = {{0}};       // heavy flavor for initiating parton for the reference gen jet
  Float_t genJetPtArrayOutput[nJetTrees][nMaxJet] = {{0}};             // pT:s of all the generator level jets in an event
  Float_t genJetPhiArrayOutput[nJetTrees][nMaxJet] = {{0}};            // phis of all the generator level jets in an event
  Float_t genJetPhiArrayWTAOutput[nJetTrees][nMaxJet] = {{0}};         // phis of all the generator level jets in an event with WTA axis
  Float_t genJetEtaArrayOutput[nJetTrees][nMaxJet] = {{0}};            // etas of all the generator level jets in an event
  Float_t genJetEtaArrayWTAOutput[nJetTrees][nMaxJet] = {{0}};         // etas of all the generator level jets in an event with WTA axis
  
  for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
    
    jetTreeOutput[iJetType] = new TTree("t","");
    
    jetTreeOutput[iJetType]->Branch("nref",&nJetsOutput[iJetType],"nref/I");
    jetTreeOutput[iJetType]->Branch("jtpt",&jetPtArrayOutput[iJetType],"jtpt[nref]/F");
    
    // Jet eta with E-scheme and WTA axes
    jetTreeOutput[iJetType]->Branch("jtphi",&jetPhiArrayOutput[iJetType],"jtphi[nref]/F");
    jetTreeOutput[iJetType]->Branch("WTAphi",&jetPhiArrayWTAOutput[iJetType],"WTAphi[nref]/F");
    
    // Jet phi with E-scheme and WTA axes
    jetTreeOutput[iJetType]->Branch("jteta",&jetEtaArrayOutput[iJetType],"jteta[nref]/F");
    jetTreeOutput[iJetType]->Branch("WTAeta",&jetEtaArrayWTAOutput[iJetType],"WTAeta[nref]/F");
    
    jetTreeOutput[iJetType]->Branch("rawpt",&jetRawPtArrayOutput[iJetType],"rawpt[nref]/F");
    jetTreeOutput[iJetType]->Branch("trackMax",&jetMaxTrackPtArrayOutput[iJetType],"trackMax[nref]/F");
    
    // If we are looking at Monte Carlo, connect the reference pT and parton arrays
    if(isMC){
      jetTreeOutput[iJetType]->Branch("refpt",&jetRefPtArrayOutput[iJetType],"refpt[nref]/F");
      jetTreeOutput[iJetType]->Branch("refparton_flavor", &jetRefFlavorArrayOutput[iJetType], "refparton_flavor[nref]/I");
      jetTreeOutput[iJetType]->Branch("refparton_flavorForB", &jetRefFlavorForBArrayOutput[iJetType], "refparton_flavorForB[nref]/I");
      
      jetTreeOutput[iJetType]->Branch("ngen",&nGenJetsOutput[iJetType],"ngen/I");
      
      jetTreeOutput[iJetType]->Branch("genpt",&genJetPtArrayOutput[iJetType],"genpt[ngen]/F");
      
      // Gen jet phi for e-scheme and WTA axes
      jetTreeOutput[iJetType]->Branch("genphi",&genJetPhiArrayOutput[iJetType],"genphi[ngen]/F");
      jetTreeOutput[iJetType]->Branch("WTAgenphi",&genJetPhiArrayWTAOutput[iJetType],"WTAgenphi[ngen]/F");
      
      // Gen jet eta for e-scheme and WTA axes
      jetTreeOutput[iJetType]->Branch("geneta",&genJetEtaArrayOutput[iJetType],"geneta[ngen]/F");
      jetTreeOutput[iJetType]->Branch("WTAgeneta",&genJetEtaArrayWTAOutput[iJetType],"WTAgeneta[ngen]/F");
      
    } // Branches only for MC

  } // Jet type loop
  
  // Copy the track trees to the output
  TTree *trackTreeOutput = new TTree("trackTree","");
  
  Int_t nTracksOutput;                                       // Number of tracks
  vector<float> *trackPtOutput = new std::vector<float>(); trackPtOutput->clear();
  vector<float> *trackPtErrorOutput = new std::vector<float>(); trackPtErrorOutput->clear();
  vector<float> *trackPhiOutput = new std::vector<float>(); trackPhiOutput->clear();
  vector<float> *trackEtaOutput = new std::vector<float>(); trackEtaOutput->clear();
  vector<bool> *trackHighPurityOutput = new std::vector<bool>(); trackHighPurityOutput->clear();
  vector<float> *trackVertexDistanceZOutput = new std::vector<float>(); trackVertexDistanceZOutput->clear();
  vector<float> *trackVertexDistanceZErrorOutput = new std::vector<float>(); trackVertexDistanceZErrorOutput->clear();
  vector<float> *trackVertexDistanceXYOutput = new std::vector<float>(); trackVertexDistanceXYOutput->clear();
  vector<float> *trackVertexDistanceXYErrorOutput = new std::vector<float>(); trackVertexDistanceXYErrorOutput->clear();
  vector<float> *trackNormalizedChi2Output = new std::vector<float>(); trackNormalizedChi2Output->clear();
  vector<char> *nHitsTrackerLayerOutput = new std::vector<char>(); nHitsTrackerLayerOutput->clear();
  vector<char> *nHitsTrackOutput = new std::vector<char>(); nHitsTrackOutput->clear();
  vector<float> *trackEnergyEcalOutput = new std::vector<float>(); trackEnergyEcalOutput->clear();
  vector<float> *trackEnergyHcalOutput = new std::vector<float>(); trackEnergyHcalOutput->clear();
  vector<char> *trackChargeOutput = new std::vector<char>(); trackChargeOutput->clear();
  
  trackTreeOutput->Branch("nTrk",&nTracksOutput,"nTrk/I");
  trackTreeOutput->Branch("trkPt","vector<float>",&trackPtOutput);
  trackTreeOutput->Branch("trkPtError","vector<float>",&trackPtErrorOutput);
  trackTreeOutput->Branch("trkPhi","vector<float>",&trackPhiOutput);
  trackTreeOutput->Branch("trkEta","vector<float>",&trackEtaOutput);
  
  trackTreeOutput->Branch("highPurity","vector<bool>",&trackHighPurityOutput);
  trackTreeOutput->Branch("trkDzFirstVtx","vector<float>",&trackVertexDistanceZOutput);
  trackTreeOutput->Branch("trkDzErrFirstVtx","vector<float>",&trackVertexDistanceZErrorOutput);
  trackTreeOutput->Branch("trkDxyFirstVtx","vector<float>",&trackVertexDistanceXYOutput);
  trackTreeOutput->Branch("trkDxyErrFirstVtx","vector<float>",&trackVertexDistanceXYErrorOutput);
  trackTreeOutput->Branch("trkNormChi2","vector<float>",&trackNormalizedChi2Output);
  trackTreeOutput->Branch("trkNLayers","vector<char>",&nHitsTrackerLayerOutput);
  trackTreeOutput->Branch("trkNHits","vector<char>",&nHitsTrackOutput);
  trackTreeOutput->Branch("pfEcal","vector<float>",&trackEnergyEcalOutput);
  trackTreeOutput->Branch("pfHcal","vector<float>",&trackEnergyHcalOutput);
  
  // Additional branches
  trackTreeOutput->Branch("trkCharge","vector<char>",&trackChargeOutput);
  
  // Generator level tracks only in Monte Carlo
  TTree *genTrackTreeOutput = new TTree("hi","");
  
  std::vector<float> *genTrackPtVector = new std::vector<float>(); genTrackPtVector->clear();
  std::vector<float> *genTrackPhiVector = new std::vector<float>(); genTrackPhiVector->clear();
  std::vector<float> *genTrackEtaVector = new std::vector<float>(); genTrackEtaVector->clear();
  std::vector<int> *genTrackPdgVector = new std::vector<int>(); genTrackPdgVector->clear();
  std::vector<int> *genTrackChargeVector = new std::vector<int>(); genTrackChargeVector->clear();
  std::vector<int> *genTrackSubeventVector = new std::vector<int>(); genTrackSubeventVector->clear();
  
  // Connect the branches to generator level track tree
  if(isMC){
    genTrackTreeOutput->Branch("pt","vector<float>", &genTrackPtVector);
    genTrackTreeOutput->Branch("phi","vector<float>", &genTrackPhiVector);
    genTrackTreeOutput->Branch("eta","vector<float>", &genTrackEtaVector);
    genTrackTreeOutput->Branch("pdg","vector<int>", &genTrackPdgVector);
    genTrackTreeOutput->Branch("chg","vector<int>", &genTrackChargeVector);
    genTrackTreeOutput->Branch("sube","vector<int>", &genTrackSubeventVector);
  }
  
  // Copy the particle flow candidate tree to the output
  TTree *particleFlowCandidateTreeOutput = new TTree("pftree","");
  
  std::vector<int> *particleFlowCandidateIdOutputVector = new std::vector<int>(); particleFlowCandidateIdOutputVector->clear();
  std::vector<float> *particleFlowCandidatePtOutputVector = new std::vector<float>(); particleFlowCandidatePtOutputVector->clear();
  std::vector<float> *particleFlowCandidatePhiOutputVector = new std::vector<float>(); particleFlowCandidatePhiOutputVector->clear();
  std::vector<float> *particleFlowCandidateEtaOutputVector = new std::vector<float>(); particleFlowCandidateEtaOutputVector->clear();
  std::vector<float> *particleFlowCandidateMassOutputVector = new std::vector<float>(); particleFlowCandidateMassOutputVector->clear();
  
  // Only define branched is requested
  if(includePfCandidates){
    particleFlowCandidateTreeOutput->Branch("pfId","vector<int>",&particleFlowCandidateIdOutputVector);
    particleFlowCandidateTreeOutput->Branch("pfPt","vector<float>",&particleFlowCandidatePtOutputVector);
    particleFlowCandidateTreeOutput->Branch("pfPhi","vector<float>",&particleFlowCandidatePhiOutputVector);
    particleFlowCandidateTreeOutput->Branch("pfEta","vector<float>",&particleFlowCandidateEtaOutputVector);
    particleFlowCandidateTreeOutput->Branch("pfM","vector<float>",&particleFlowCandidateMassOutputVector);
  }
  
  // ========================================== //
  //          Loop over all events              //
  // ========================================== //
  
  int nEvents = heavyIonTree->GetEntries();
  cout << "There are " << nEvents << " events" << endl;
  
  bool passTrackCuts;
  bool passJetCuts;
  int iJetOutput;
  
  for(int iEvent = 0; iEvent < nEvents; iEvent++) {
    
    if( iEvent % 1000 == 0 )  std::cout << "iEvent: " << iEvent <<  " of " << nEvents << std::endl;
    
    // ========================================== //
    //        Read the event to input trees
    // ========================================== //
    
    heavyIonTree->GetEntry(iEvent);
    hltTree->GetEntry(iEvent);
    skimTree->GetEntry(iEvent);
    trackTree->GetEntry(iEvent);
    if(isMC) genTrackTree->GetEntry(iEvent);
    particleFlowCandidateTree->GetEntry(iEvent);
    
    for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
      jetTree[iJetType]->GetEntry(iEvent);
    }
    
    heavyIonTreeOutput->Fill();
    hltTreeOutput->Fill();
    skimTreeOutput->Fill();
    
    // Fill jet histograms using basic jet cuts
    for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
      
      iJetOutput = 0;
      nJetsOutput[iJetType] = nJets[iJetType];
      
      for(int iJet = 0; iJet < nJets[iJetType]; iJet++){
        
        passJetCuts = true;
        
        // Apply very basic jet cuts
        if(jetPtArray[iJetType][iJet] < 25) passJetCuts = false;    // Minumum pT cut of 25 GeV
        if(jetEtaArray[iJetType][iJet] > 2 && jetEtaArrayWTA[iJetType][iJet] > 2) passJetCuts = false;    // Maximum eta cut of 2
        
        // Fill the jet arrays with reconstructed jets
        if(passJetCuts){
          jetPtArrayOutput[iJetType][iJetOutput] = jetPtArray[iJetType][iJet];
          jetPhiArrayOutput[iJetType][iJetOutput] = jetPhiArray[iJetType][iJet];
          jetPhiArrayWTAOutput[iJetType][iJetOutput] = jetPhiArrayWTA[iJetType][iJet];
          jetEtaArrayOutput[iJetType][iJetOutput] = jetEtaArray[iJetType][iJet];
          jetEtaArrayWTAOutput[iJetType][iJetOutput] = jetEtaArrayWTA[iJetType][iJet] ;
          jetRawPtArrayOutput[iJetType][iJetOutput] = jetRawPtArray[iJetType][iJet];
          jetMaxTrackPtArrayOutput[iJetType][iJetOutput] = jetMaxTrackPtArray[iJetType][iJet];
          
          if(isMC){
            jetRefPtArrayOutput[iJetType][iJetOutput] = jetRefPtArray[iJetType][iJet];
            jetRefFlavorArrayOutput[iJetType][iJetOutput] = jetRefFlavorArray[iJetType][iJet];
            jetRefFlavorForBArrayOutput[iJetType][iJetOutput] = jetRefFlavorForBArray[iJetType][iJet];
          }
          
          iJetOutput++;
        } else {
          nJetsOutput[iJetType]--;
        }
        
      } // Reconstructed jet loop
      
      if(isMC){
        
        iJetOutput = 0;
        nGenJetsOutput[iJetType] = nGenJets[iJetType];
        
        for(int iJet = 0; iJet < nGenJets[iJetType]; iJet++){
          
          passJetCuts = true;
          
          // Apply very basic jet cuts
          if(genJetPtArray[iJetType][iJet] < 25) passJetCuts = false;    // Minumum pT cut of 25 GeV
          if(genJetEtaArray[iJetType][iJet] > 2 && genJetEtaArrayWTA[iJetType][iJet] > 2) passJetCuts = false;    // Maximum eta cut of 2
          
          // Fill the jet arrays with generated jets
          if(passJetCuts){
            
            genJetPtArrayOutput[iJetType][iJetOutput] = genJetPtArray[iJetType][iJet];
            genJetPhiArrayOutput[iJetType][iJetOutput] = genJetPhiArray[iJetType][iJet];
            genJetPhiArrayWTAOutput[iJetType][iJetOutput] = genJetPhiArrayWTA[iJetType][iJet];
            genJetEtaArrayOutput[iJetType][iJetOutput] = genJetEtaArray[iJetType][iJet];
            genJetEtaArrayWTAOutput[iJetType][iJetOutput] = genJetEtaArrayWTA[iJetType][iJet];
            
            iJetOutput++;
            
          } else {
            nGenJetsOutput[iJetType]--;
          }
          
        } // Generator level jet loop
        
      } // If for filling generator jet loop
      
      jetTreeOutput[iJetType]->Fill();
      
    } // Loop over jet collections
    
    
    // Reco track loop
    for(int iTrack = 0; iTrack < nTracks; iTrack++){
      
      passTrackCuts = true;
      
      // Do basic track cuts for the reconstructed tracks
      if(trackHighPurityArray->at(iTrack) != 1) passTrackCuts = false;
      
      if(trackPtErrorArray->at(iTrack)/trackPtArray->at(iTrack) > 0.1) passTrackCuts = false;
      if(TMath::Abs(trackVertexDistanceZArray->at(iTrack)/trackVertexDistanceZErrorArray->at(iTrack)) > 3) passTrackCuts = false;
      if(TMath::Abs(trackVertexDistanceXYArray->at(iTrack)/trackVertexDistanceXYErrorArray->at(iTrack)) > 3) passTrackCuts = false;
      
      if(TMath::Abs(trackEtaArray->at(iTrack)) >= 2.4) passTrackCuts = false;  //acceptance of the tracker
      
      if(trackPtArray->at(iTrack) < 0.7) passTrackCuts = false;   // Minimum track pT
      
      if(passTrackCuts){
        trackPtOutput->push_back(trackPtArray->at(iTrack));
        trackPtErrorOutput->push_back(trackPtErrorArray->at(iTrack));
        trackPhiOutput->push_back(trackPhiArray->at(iTrack));
        trackEtaOutput->push_back(trackEtaArray->at(iTrack));
        trackHighPurityOutput->push_back(trackHighPurityArray->at(iTrack));
        trackVertexDistanceZOutput->push_back(trackVertexDistanceZArray->at(iTrack));
        trackVertexDistanceZErrorOutput->push_back(trackVertexDistanceZErrorArray->at(iTrack));
        trackVertexDistanceXYOutput->push_back(trackVertexDistanceXYArray->at(iTrack));
        trackVertexDistanceXYErrorOutput->push_back(trackVertexDistanceXYErrorArray->at(iTrack));
        trackNormalizedChi2Output->push_back(trackNormalizedChi2Array->at(iTrack));
        nHitsTrackerLayerOutput->push_back(nHitsTrackerLayerArray->at(iTrack));
        nHitsTrackOutput->push_back(nHitsTrackArray->at(iTrack));
        trackEnergyEcalOutput->push_back(trackEnergyEcalArray->at(iTrack));
        trackEnergyHcalOutput->push_back(trackEnergyHcalArray->at(iTrack));
        trackChargeOutput->push_back(trackChargeArray->at(iTrack));
      }
    }
    
    nTracksOutput = trackPtOutput->size();
    trackTreeOutput->Fill();
    
    // Clear the vectors before the next event! Otherwise all the tracks pile up cumulatively
    trackPtOutput->clear();
    trackPtErrorOutput->clear();
    trackPhiOutput->clear();
    trackEtaOutput->clear();
    trackHighPurityOutput->clear();
    trackVertexDistanceZOutput->clear();
    trackVertexDistanceZErrorOutput->clear();
    trackVertexDistanceXYOutput->clear();
    trackVertexDistanceXYErrorOutput->clear();
    trackNormalizedChi2Output->clear();
    nHitsTrackerLayerOutput->clear();
    nHitsTrackOutput->clear();
    trackEnergyEcalOutput->clear();
    trackEnergyHcalOutput->clear();
    trackChargeOutput->clear();
    
    // Gen track loop
    if(isMC){
      for(int iTrack = 0; iTrack < genTrackPtArray->size(); iTrack++){
        
        // Cut away low pT tracks and tracks with eta outside of tracker acceptance
        if(TMath::Abs(genTrackEtaArray->at(iTrack)) >= 2.4) continue; //acceptance of the tracker
        
        if(genTrackPtArray->at(iTrack) < 0.7) continue;   // Minimum track pT
        //if(genTrackPtArray->at(iTrack) > 300 ) continue;  // Maximum track pT
        
        // Fill the output vectors with gen particles surviving the cuts
        genTrackPtVector->push_back(genTrackPtArray->at(iTrack));
        genTrackPhiVector->push_back(genTrackPhiArray->at(iTrack));
        genTrackEtaVector->push_back(genTrackEtaArray->at(iTrack));
        genTrackPdgVector->push_back(genTrackPdgArray->at(iTrack));
        genTrackChargeVector->push_back(genTrackChargeArray->at(iTrack));
        genTrackSubeventVector->push_back(genTrackSubeventArray->at(iTrack));
      }
      
      genTrackTreeOutput->Fill();
      
    } // Filling gen tracks for MC
    
    // Particle flow candidate loop. Do only if requested
    if(includePfCandidates){
      for(int pfi = 0; pfi < particleFlowCandidateIdVector->size(); pfi++) {
        
        // pT cut for PF condidates
        if(particleFlowCandidatePtVector->at(pfi) < 1) continue;  // Minimum PF candidate pT cut of 1 GeV
        
        particleFlowCandidateIdOutputVector->push_back(particleFlowCandidateIdVector->at(pfi));
        particleFlowCandidatePtOutputVector->push_back(particleFlowCandidatePtVector->at(pfi));
        particleFlowCandidatePhiOutputVector->push_back(particleFlowCandidatePhiVector->at(pfi));
        particleFlowCandidateEtaOutputVector->push_back(particleFlowCandidateEtaVector->at(pfi));
        particleFlowCandidateMassOutputVector->push_back(particleFlowCandidateMassVector->at(pfi));
        
      } // particle flow candidate loop
      
      particleFlowCandidateTreeOutput->Fill();
    }
    
    // Clear the vectors before the next event! Otherwise all the tracks pile up cumulatively
    if(isMC){
      genTrackPtVector->clear();
      genTrackPhiVector->clear();
      genTrackEtaVector->clear();
      genTrackPdgVector->clear();
      genTrackChargeVector->clear();
      genTrackSubeventVector->clear();
    }
    
    if(includePfCandidates){
      particleFlowCandidateIdOutputVector->clear();
      particleFlowCandidatePtOutputVector->clear();
      particleFlowCandidatePhiOutputVector->clear();
      particleFlowCandidateEtaOutputVector->clear();
      particleFlowCandidateMassOutputVector->clear();
    }
    
  } // Event loop
  
  // Write the skimmed trees to the output file
  
  TFile *outputFile = new TFile(outputFileName, "RECREATE");
  
  gDirectory->mkdir("HiForestInfo");
  gDirectory->cd("HiForestInfo");
  
  hiForestInfoTreeOutput->Write();
  
  gDirectory->cd("../");
  gDirectory->mkdir("hiEvtAnalyzer");
  gDirectory->cd("hiEvtAnalyzer");
  
  heavyIonTreeOutput->Write();
  
  gDirectory->cd("../");
  gDirectory->mkdir("hltanalysis");
  gDirectory->cd("hltanalysis");
  
  hltTreeOutput->Write();
  
  gDirectory->cd("../");
  gDirectory->mkdir("skimanalysis");
  gDirectory->cd("skimanalysis");
  
  skimTreeOutput->Write();
  
  gDirectory->cd("../");
  
  const char *jetDirectories[] = {"akCs4PFJetAnalyzer","akFlowPuCs4PFJetAnalyzer"};
  
  for(int iJetType = 0; iJetType < nJetTrees; iJetType++){
    
    gDirectory->mkdir(jetDirectories[iJetType]);
    gDirectory->cd(jetDirectories[iJetType]);
    
    jetTreeOutput[iJetType]->Write();
    
    gDirectory->cd("../");
    
  } // Loop over jet types
  
  gDirectory->mkdir("PbPbTracks");
  gDirectory->cd("PbPbTracks");
  
  trackTreeOutput->Write();
  
  gDirectory->cd("../");
  
  // Generator particles only present in MC
  if(isMC){
    gDirectory->mkdir("HiGenParticleAna");
    gDirectory->cd("HiGenParticleAna");
    
    genTrackTreeOutput->Write();
    
    gDirectory->cd("../");
  }
  
  // Only write PF candidates if requested
  if(includePfCandidates){
    gDirectory->mkdir("particleFlowAnalyser");
    gDirectory->cd("particleFlowAnalyser");
    
    particleFlowCandidateTreeOutput->Write();
    
    gDirectory->cd("../");
  }
  
  outputFile->Close();
  
  return 0;
  
}
