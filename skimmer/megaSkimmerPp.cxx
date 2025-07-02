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

// Own includes for corrected multiplicity
#include "trackingEfficiency2017pp.h"
#include "TrackingEfficiencyInterface.h"

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
  const char* fileLocation[] = {"root://xrootd.rcac.purdue.edu/", "root://eoscms.cern.ch/", "root://xrootd-vanderbilt.sites.opensciencegrid.org/", "root://cmsxrootd.fnal.gov/"};
  
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
int main(int argc, char* argv[]){
  
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
  ReadFileList(fileNameVector,fileNameFile, 1,fileSearchIndex,runLocal);

  // Tracking efficiency correction needed for corrected multiplicity
  TrackingEfficiencyInterface* trackEfficiencyCorrector = new TrkEff2017pp(false, "trackCorrectionTables/pp2017/");

  // Open the output file
  TFile* outputFile = new TFile(outputFileName, "RECREATE");
  
  // Define trees to be read from the files
  TChain* heavyIonTree = new TChain("hiEvtAnalyzer/HiTree");
  TChain* skimTree = new TChain("skimanalysis/HltTree");
  TChain* trackTree = new TChain("ppTracks/trackTree");
  TChain* genTrackTree = new TChain("HiGenParticleAna/hi");
  
  // All the branches and leaves come in array of two, one for input and one for output

  // Branches for heavy ion tree
  TBranch* runBranch;             // Branch for run
  TBranch* eventBranch;           // Branch for event
  TBranch* lumiBranch;            // Branch for lumi
  TBranch* hiVzBranch;            // Branch for vertex z-position
  TBranch* hiHFplusBranch;        // Branch for energy in HF plus
  TBranch* hiHFminusBranch;       // Branch for energy in HF minus
  TBranch* ptHatBranch;           // Branch for pT hat
  TBranch* eventWeightBranch;     // Branch for jet weight for 2018 MC
  
  // Leaves for heavy ion tree
  UInt_t run;           // Run number
  ULong64_t event;      // Event number
  UInt_t lumi;          // Luminosity block
  Float_t vertexZ;      // Vertex z-position
  Float_t hiHFplus;     // Energy in HF plus
  Float_t hiHFminus;    // Energy in HF minus
  Float_t ptHat;        // pT hat
  Float_t eventWeight;  // jet weight in the 2018 MC tree
  
  // Branches for skim tree
  TBranch* primaryVertexBranch;             // Branch for primary vertex filter bit
  TBranch* pileupBranch;                    // Branch for pileup
  
  // Leaves for the skim tree
  Int_t primaryVertexFilterBit;             // Filter bit for primary vertex
  Int_t pileupFilterBit;                    // Filter bit for pileup
  
  // Branches for track tree
  TBranch* nTracksBranch;                    // Branch for number of tracks
  TBranch* trackPtBranch;                    // Branch for track pT
  TBranch* trackPtErrorBranch;               // Branch for track pT error
  TBranch* trackPhiBranch;                   // Branch for track phi
  TBranch* trackEtaBranch;                   // Branch for track eta
  TBranch* trackHighPurityBranch;            // Branch for high purity of the track
  TBranch* trackVertexDistanceZBranch;       // Branch for track distance from primary vertex in z-direction
  TBranch* trackVertexDistanceZErrorBranch;  // Branch for error for track distance from primary vertex in z-direction
  TBranch* trackVertexDistanceXYBranch;      // Branch for track distance from primary vertex in xy-direction
  TBranch* trackVertexDistanceXYErrorBranch; // Branch for error for track distance from primary vertex in xy-direction
  TBranch* trackChargeBranch;                // Branch for track charge
  
  // Leaves for the track tree
  Int_t nTracks;                                                 // Number of tracks
  std::vector<float>* trackPtArray = nullptr;                    // Array for track pT:s
  std::vector<float>* trackPtErrorArray = nullptr;               // Array for track pT errors
  std::vector<float>* trackPhiArray = nullptr;                   // Array for track phis
  std::vector<float>* trackEtaArray = nullptr;                   // Array for track etas
  std::vector<bool>* trackHighPurityArray = nullptr;             // Array for the high purity of tracks
  std::vector<float>* trackVertexDistanceZArray = nullptr;       // Array for track distance from primary vertex in z-direction
  std::vector<float>* trackVertexDistanceZErrorArray = nullptr;  // Array for error for track distance from primary vertex in z-direction
  std::vector<float>* trackVertexDistanceXYArray = nullptr;      // Array for track distance from primary vertex in xy-direction
  std::vector<float>* trackVertexDistanceXYErrorArray = nullptr; // Array for error for track distance from primary vertex in xy-direction
  std::vector<char>* trackChargeArray = nullptr;                 // Array for track charge
  
  // Branches for generator level track tree
  TBranch* genTrackPtBranch;         // Branch for generator level track pT:s
  TBranch* genTrackPhiBranch;        // Branch for generator level track phis
  TBranch* genTrackEtaBranch;        // Branch for generator level track etas
  TBranch* genTrackSubeBranch;       // Branch for generator level track subevent index
  TBranch* genTrackChargeBranch;     // Branch for generator level track charges
  
  // Leaves for generator level track tree
  std::vector<float>* genTrackPtArray = nullptr;       // Array for generator level track pT:s
  std::vector<float>* genTrackPhiArray = nullptr;      // Array for generator level track phis
  std::vector<float>* genTrackEtaArray = nullptr;      // Array for generator level track etas
  std::vector<int>* genTrackSubeArray = nullptr;       // Array for generator level track subevent index
  std::vector<int>* genTrackChargeArray = nullptr;     // Array for generator level track charges
  
  // Add all the files to the chain
  for(std::vector<TString>::iterator listIterator = fileNameVector.begin(); listIterator != fileNameVector.end(); listIterator++){
    
    cout << "Adding file " << *listIterator << " to the chains" << endl;
    
    heavyIonTree->Add(*listIterator);
    skimTree->Add(*listIterator);
    trackTree->Add(*listIterator);
    genTrackTree->Add(*listIterator);
    
  }

  
  // ========================================== //
  // Read all the branches from the input trees //
  // ========================================== //
  
  // Connect the branches of the heavy ion tree
  heavyIonTree->SetBranchStatus("*", 0);
  heavyIonTree->SetBranchStatus("run", 1);
  heavyIonTree->SetBranchAddress("run", &run, &runBranch);
  heavyIonTree->SetBranchStatus("evt", 1);
  heavyIonTree->SetBranchAddress("evt", &event, &eventBranch);
  heavyIonTree->SetBranchStatus("lumi", 1);
  heavyIonTree->SetBranchAddress("lumi", &lumi, &lumiBranch);
  heavyIonTree->SetBranchStatus("vz", 1);
  heavyIonTree->SetBranchAddress("vz", &vertexZ, &hiVzBranch);
  heavyIonTree->SetBranchStatus("hiHFPlus_pf", 1);
  heavyIonTree->SetBranchAddress("hiHFPlus_pf", &hiHFplus, &hiHFplusBranch);
  heavyIonTree->SetBranchStatus("hiHFMinus_pf", 1);
  heavyIonTree->SetBranchAddress("hiHFMinus_pf", &hiHFminus, &hiHFminusBranch);
  
  // ptHat and event weight only for MC
  if(isMC){
    heavyIonTree->SetBranchStatus("pthat", 1);
    heavyIonTree->SetBranchAddress("pthat", &ptHat, &ptHatBranch);
    heavyIonTree->SetBranchStatus("weight", 1);
    heavyIonTree->SetBranchAddress("weight", &eventWeight, &eventWeightBranch);
  }
  
  // Connect the branches to the skim tree
  skimTree->SetBranchStatus("*", 0);

  skimTree->SetBranchStatus("pprimaryVertexFilter", 1);
  skimTree->SetBranchAddress("pprimaryVertexFilter", &primaryVertexFilterBit, &primaryVertexBranch);

  skimTree->SetBranchStatus("pVertexFilterCutGplusUpsPP", 1);
  skimTree->SetBranchAddress("pVertexFilterCutGplusUpsPP", &pileupFilterBit, &pileupBranch);
  
  // Connect the branches to the track tree
  trackTree->SetBranchStatus("*", 0);
  
  trackTree->SetBranchStatus("trkPt", 1);
  trackTree->SetBranchAddress("trkPt", &trackPtArray, &trackPtBranch);
  trackTree->SetBranchStatus("trkPtError", 1);
  trackTree->SetBranchAddress("trkPtError", &trackPtErrorArray, &trackPtErrorBranch);
  trackTree->SetBranchStatus("trkPhi", 1);
  trackTree->SetBranchAddress("trkPhi", &trackPhiArray, &trackPhiBranch);
  trackTree->SetBranchStatus("trkEta", 1);
  trackTree->SetBranchAddress("trkEta", &trackEtaArray, &trackEtaBranch);
  trackTree->SetBranchStatus("nTrk", 1);
  trackTree->SetBranchAddress("nTrk", &nTracks, &nTracksBranch);
  trackTree->SetBranchStatus("highPurity", 1);
  trackTree->SetBranchAddress("highPurity", &trackHighPurityArray, &trackHighPurityBranch);
  trackTree->SetBranchStatus("trkDzFirstVtx", 1);
  trackTree->SetBranchAddress("trkDzFirstVtx", &trackVertexDistanceZArray, &trackVertexDistanceZBranch);
  trackTree->SetBranchStatus("trkDzErrFirstVtx", 1);
  trackTree->SetBranchAddress("trkDzErrFirstVtx", &trackVertexDistanceZErrorArray, &trackVertexDistanceZErrorBranch);
  trackTree->SetBranchStatus("trkDxyFirstVtx", 1);
  trackTree->SetBranchAddress("trkDxyFirstVtx", &trackVertexDistanceXYArray, &trackVertexDistanceXYBranch);
  trackTree->SetBranchStatus("trkDxyErrFirstVtx", 1);
  trackTree->SetBranchAddress("trkDxyErrFirstVtx", &trackVertexDistanceXYErrorArray, &trackVertexDistanceXYErrorBranch);
  
  // Additional saved branches
  trackTree->SetBranchStatus("trkCharge", 1);
  trackTree->SetBranchAddress("trkCharge", &trackChargeArray, &trackChargeBranch);
  
  // Generator level tracks only in Monte Carlo
  if(isMC){
    
    // Connect the branches to generator level track tree
    genTrackTree->SetBranchStatus("*",0);
    
    genTrackTree->SetBranchStatus("pt", 1);
    genTrackTree->SetBranchAddress("pt", &genTrackPtArray, &genTrackPtBranch);
    genTrackTree->SetBranchStatus("phi", 1);
    genTrackTree->SetBranchAddress("phi", &genTrackPhiArray, &genTrackPhiBranch);
    genTrackTree->SetBranchStatus("eta", 1);
    genTrackTree->SetBranchAddress("eta", &genTrackEtaArray, &genTrackEtaBranch);
    genTrackTree->SetBranchStatus("sube", 1);
    genTrackTree->SetBranchAddress("sube", &genTrackSubeArray, &genTrackSubeBranch);
    genTrackTree->SetBranchStatus("chg", 1);
    genTrackTree->SetBranchAddress("chg", &genTrackChargeArray, &genTrackChargeBranch);
    
  }
  
  
  // ========================================== //
  //           Define output trees
  // ========================================== //

  // Copy the heavy ion tree to the output
  TTree* heavyIonTreeOutput = new TTree("HiTree","");
  
  // Connect the branches of the heavy ion tree
  heavyIonTreeOutput->Branch("run" , &run, "run/i");
  heavyIonTreeOutput->Branch("evt", &event, "evt/l");
  heavyIonTreeOutput->Branch("lumi", &lumi, "lumi/i");
  heavyIonTreeOutput->Branch("vz", &vertexZ, "vz/F");
  heavyIonTreeOutput->Branch("hiHFPlus_pf", &hiHFplus, "hiHFPlus_pf/F");
  heavyIonTreeOutput->Branch("hiHFMinus_pf", &hiHFminus, "hiHFMinus_pf/F");
  
  // ptHat and event weight only for MC
  if(isMC){
    heavyIonTreeOutput->Branch("pthat", &ptHat, "pthat/F");
    heavyIonTreeOutput->Branch("weight", &eventWeight, "weight/F");
  }
  
  // Copy the track trees to the output
  TTree* trackTreeOutput = new TTree("trackTree","");
  
  Int_t nTracksOutput;              // Number of tracks
  std::vector<float>* trackPtOutput = new std::vector<float>(); trackPtOutput->clear();     // Track pT
  std::vector<float>* trackPhiOutput = new std::vector<float>(); trackPhiOutput->clear();   // Track phi
  std::vector<float>* trackEtaOutput = new std::vector<float>(); trackEtaOutput->clear();   // Track eta
  Float_t nTracksOutputCorrected;   // Number of efficiency corrected tracks
  
  trackTreeOutput->Branch("nTrk", &nTracksOutput, "nTrk/I");
  trackTreeOutput->Branch("nTrkCorrected", &nTracksOutputCorrected , "nTrkCorrected/F");
  trackTreeOutput->Branch("trkPt", "vector<float>", &trackPtOutput);
  trackTreeOutput->Branch("trkPhi", "vector<float>", &trackPhiOutput);
  trackTreeOutput->Branch("trkEta", "vector<float>", &trackEtaOutput);
  
  // Generator level tracks only in Monte Carlo
  TTree* genTrackTreeOutput = new TTree("hi","");
  
  std::vector<float> *genTrackPtVector = new std::vector<float>(); genTrackPtVector->clear();
  std::vector<float> *genTrackPhiVector = new std::vector<float>(); genTrackPhiVector->clear();
  std::vector<float> *genTrackEtaVector = new std::vector<float>(); genTrackEtaVector->clear();
  std::vector<int> *genTrackSubeVector = new std::vector<int>(); genTrackSubeVector->clear();
  
  // Connect the branches to generator level track tree
  if(isMC){
    genTrackTreeOutput->Branch("pt", "vector<float>", &genTrackPtVector);
    genTrackTreeOutput->Branch("phi", "vector<float>", &genTrackPhiVector);
    genTrackTreeOutput->Branch("eta", "vector<float>", &genTrackEtaVector);
    genTrackTreeOutput->Branch("sube", "vector<int>", &genTrackSubeVector);
  }
  
  // ========================================== //
  //          Loop over all events              //
  // ========================================== //
  
  int nEvents = heavyIonTree->GetEntries();
  cout << "There are " << nEvents << " events" << endl;
  
  int iTrackOutput;
  double iTrackOutputCorrected;

  for(int iEvent = 0; iEvent < nEvents; iEvent++) {
    
    if( iEvent % 1000 == 0 ) std::cout << "iEvent: " << iEvent <<  " of " << nEvents << std::endl;
    
    // ========================================== //
    //        Read the event to input trees
    // ========================================== //
    
    heavyIonTree->GetEntry(iEvent);
    skimTree->GetEntry(iEvent);
    trackTree->GetEntry(iEvent);
    if(isMC) genTrackTree->GetEntry(iEvent);

    // Do basic event selection based on skim tree variables
    if(primaryVertexFilterBit == 0) continue;
    if(pileupFilterBit == 0) continue;
    if(vertexZ < -15) continue;
    if(vertexZ > 15) continue;
    
    heavyIonTreeOutput->Fill();
    
    // Reco track loop
    iTrackOutput = 0;
    iTrackOutputCorrected = 0;
    for(int iTrack = 0; iTrack < nTracks; iTrack++){
            
      // Do basic track cuts for the reconstructed tracks
      if(trackHighPurityArray->at(iTrack) != 1) continue; // High purity cut

      if(trackPtArray->at(iTrack) < 0.7) continue;   // Minimum track pT
      if(trackPtArray->at(iTrack) > 500) continue; // Maximum track pT

      if(TMath::Abs(trackEtaArray->at(iTrack)) >= 2.4) continue;  //acceptance of the tracker
      
      // pT and vertex resolution cuts
      if(trackPtErrorArray->at(iTrack)/trackPtArray->at(iTrack) > 0.1) continue;
      if(TMath::Abs(trackVertexDistanceZArray->at(iTrack)/trackVertexDistanceZErrorArray->at(iTrack)) > 3) continue;
      if(TMath::Abs(trackVertexDistanceXYArray->at(iTrack)/trackVertexDistanceXYErrorArray->at(iTrack)) > 3) continue;
      
      trackPtOutput->push_back(trackPtArray->at(iTrack));
      trackPhiOutput->push_back(trackPhiArray->at(iTrack));
      trackEtaOutput->push_back(trackEtaArray->at(iTrack));

      iTrackOutput++;
      iTrackOutputCorrected += trackEfficiencyCorrector->getCorrection(trackPtArray->at(iTrack), trackEtaArray->at(iTrack), 0);

    }
    
    nTracksOutput = iTrackOutput;
    nTracksOutputCorrected = iTrackOutputCorrected;
    trackTreeOutput->Fill();

    // Clear the vectors before the next event! Otherwise all the tracks pile up cumulatively
    trackPtOutput->clear();
    trackPhiOutput->clear();
    trackEtaOutput->clear();
    
    // Gen track loop
    if(isMC){
      for(int iTrack = 0; iTrack < genTrackPtArray->size(); iTrack++){
        
        // Cut away low pT tracks and tracks with eta outside of tracker acceptance
        if(TMath::Abs(genTrackEtaArray->at(iTrack)) >= 2.4) continue; //acceptance of the tracker
        
        if(genTrackPtArray->at(iTrack) < 0.7) continue;   // Minimum track pT
        if(genTrackPtArray->at(iTrack) > 500 ) continue;  // Maximum track pT

        if(genTrackChargeArray->at(iTrack) == 0) continue; // Require the particle to be charged
        
        // Fill the output vectors with gen particles surviving the cuts
        genTrackPtVector->push_back(genTrackPtArray->at(iTrack));
        genTrackPhiVector->push_back(genTrackPhiArray->at(iTrack));
        genTrackEtaVector->push_back(genTrackEtaArray->at(iTrack));
        genTrackSubeVector->push_back(genTrackSubeArray->at(iTrack));
      }
      
      genTrackTreeOutput->Fill();

      // Clear the vectors before the next event! Otherwise all the tracks pile up cumulatively
      genTrackPtVector->clear();
      genTrackPhiVector->clear();
      genTrackEtaVector->clear();
      genTrackSubeVector->clear();
      
    } // Filling gen tracks for MC
    
  } // Event loop
  
  // Write the skimmed trees to the output file  
  gDirectory->mkdir("hiEvtAnalyzer");
  gDirectory->cd("hiEvtAnalyzer");
  
  heavyIonTreeOutput->Write();
  
  gDirectory->cd("../");
  
  gDirectory->mkdir("ppTracks");
  gDirectory->cd("ppTracks");
  
  trackTreeOutput->Write();
  
  gDirectory->cd("../");
  
  // Generator particles only present in MC
  if(isMC){
    gDirectory->mkdir("HiGenParticleAna");
    gDirectory->cd("HiGenParticleAna");
    
    genTrackTreeOutput->Write();
    
    gDirectory->cd("../");
  }
  
  outputFile->Close();
  
  return 0;
  
}
