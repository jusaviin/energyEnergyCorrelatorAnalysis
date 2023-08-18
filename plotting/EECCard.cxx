/*
 * Implementation of the EECCard class
 */

// Own includes
#include "EECCard.h"

/*
 * Contructor with input file
 *
 *  TFile *inFile = Input file
 */
EECCard::EECCard(TFile *inFile):
  fInputFile(inFile),
  fCardDirectory("JCard"),
  fDataType(-1),
  fDataTypeString(""),
  fAlternativeDataTypeString(""),
  fGitHash(0),
  fProjectionGitHash(0),
  fProcessGitHash(0),
  fUnfoldingGitHash(0)
{

  // Read the vectors from the input file
  fInputFile->cd(fCardDirectory.Data());
  ReadVectors();
  fDataType = (*fCardEntries[kDataType])[1];
  
  // Read the Monte Carlo type from the vector: 0 = RecoReco, 1 = RecoGen, 2 = GenReco, 3 = GenGen
  if(fCardEntries[kMcCorrelationType]) {
    fMonteCarloType = (*fCardEntries[kMcCorrelationType])[1];
  } else {
    fMonteCarloType = 0;
  }
  
  FindDataTypeString();
  
}

/*
 * Destructor
 */
EECCard::~EECCard(){

}

/*
 * Reader for all the vectors from the input file
 */
void EECCard::ReadVectors(){
  
  // Read the git hash
  fGitHash = (TObjString*) gDirectory->Get("GitHash");
  fProjectionGitHash = (TObjString*) gDirectory->Get("ProjectionGitHash");
  fProcessGitHash = (TObjString*) gDirectory->Get("ProcessGitHash");
  fUnfoldingGitHash = (TObjString*) gDirectory->Get("UnfoldingGitHash");

  // Read the TVectorT<float>:s
  for(int iEntry = 0; iEntry < knEntries; iEntry++){
    fCardEntries[iEntry] = (TVectorT<float>*) gDirectory->Get(fCardEntryNames[iEntry]);
  }
  
  // Read the file names
  for(int iFileName = 0; iFileName < knFileNames; iFileName++){
    fFileNames[iFileName] = (TObjString*) gDirectory->Get(fFileNameSaveName[iFileName]);
  }

}

/*
 * Construct a data type string based on information on the card
 */
void EECCard::FindDataTypeString(){
  
  // Define the different data types corresponding to certain indices
  TString dataTypes[5] = {"pp","PbPb","pp MC","PbPb MC","localTest"};
  TString alternativeDataTypes[5] = {"pp","PbPb","Pythia8","Pythia+Hydjet","localTest"};
  if(fDataType < 0 || fDataType > 4){
    fDataTypeString = "Unknown";
    fAlternativeDataTypeString = "Unknown";
    fDataTypeStringWithoutMCType = "Unknown";
    return;
  }
  
  // Remember the data type without the MC type
  fDataTypeStringWithoutMCType = alternativeDataTypes[fDataType];
  
  // Define Monte Carlo types and add them to MC productions, which are data types 2 and 3
  TString monteCarloString[4] = {" RecoReco"," RecoGen"," GenReco"," GenGen"};
  if(fMonteCarloType < 0 || fMonteCarloType > 3){
    fDataTypeString = "Unknown";
    fAlternativeDataTypeString = "Unknown";
    return;
  }
  
  for(int iDataType = 2; iDataType <=3; iDataType++){
    dataTypes[iDataType].Append(monteCarloString[fMonteCarloType]);
    alternativeDataTypes[iDataType].Append(monteCarloString[fMonteCarloType]);
  }
  
  // Remember the constructed data type string
  fDataTypeString = dataTypes[fDataType];
  fAlternativeDataTypeString = alternativeDataTypes[fDataType];
}

/*
 *  Getter for data type string
 */
TString EECCard::GetDataType() const{
  return fDataTypeString;
}

/*
 *  Getter for alternative data type string
 */
TString EECCard::GetAlternativeDataType(const bool includeMCtype) const{
  if(includeMCtype) return fAlternativeDataTypeString;
  return fDataTypeStringWithoutMCType;
}

// Getter for the first unfolded centrality bin index
int EECCard::GetFirstUnfoldedCentralityBin() const{
  return (*fCardEntries[kFirstUnfoldedCentralityBin])[1];
}

// Getter for the last unfolded centrality bin index
int EECCard::GetLastUnfoldedCentralityBin() const{
  return (*fCardEntries[kLastUnfoldedCentralityBin])[1];
}  

// Getter for the first unfolded track pT bin index
int EECCard::GetFirstUnfoldedTrackPtBin() const{
  return (*fCardEntries[kFirstUnfoldedTrackPtBin])[1];
}

// Getter for the last unfolded track pT bin index
int EECCard::GetLastUnfoldedTrackPtBin() const{
  return (*fCardEntries[kLastUnfoldedTrackPtBin])[1];
} 

// Getter for the first unfolded jet pT bin index
int EECCard::GetFirstUnfoldedJetPtBin() const{
  return (*fCardEntries[kFirstUnfoldedJetPtBin])[1];
} 

// Getter for the last unfolded jet pT bin index
int EECCard::GetLastUnfoldedJetPtBin() const{
  return (*fCardEntries[kLastUnfoldedJetPtBin])[1];
} 

/*
 *  Getter for subevent cut
 */
int EECCard::GetSubeventCut() const{
  return (*fCardEntries[kSubeventCut])[1];
}

/*
 *  Getter for jet type
 */
int EECCard::GetJetType() const{
  return (*fCardEntries[kJetType])[1];
}

// Getter for the minimum jet pT cut
double EECCard::GetJetPtCut() const{
  return (*fCardEntries[kMinPtCut])[1];
}

// Getter for the information if reflected cone histograms are filled
bool EECCard::GetDoReflectedCone() const{
  return !((*fCardEntries[kDoReflectedCone])[1] == 0);
}

// Getter for the information if reflected cone QA histograms are filled
bool EECCard::GetDoReflectedConeQA() const{
  return ((*fCardEntries[kDoReflectedCone])[1] >= 2);
}


/*
 * Get the number of bins for internal index
 * If no vector is found in the index, return 0.
 */
int EECCard::GetNBins(const int index) const{
  if(fCardEntries[index]) return fCardEntries[index]->GetNoElements()-1;
  return 0;
}

// Get the number of centrality bins
int EECCard::GetNCentralityBins() const{
  return GetNBins(kCentralityBinEdges);
}

// Get the number of track pT bins
int EECCard::GetNTrackPtBins() const{
  return GetNBins(kTrackPtBinEdges);
}

// Get the number of jet pT bins in energy-energy correlator analysis
int EECCard::GetNJetPtBinsEEC() const{
  return GetNBins(kJetPtBinEdgesEEC);
}

// Get the number of track pT bins in energy-energy correlator analysis
int EECCard::GetNTrackPtBinsEEC() const{
  return GetNBins(kTrackPtBinEdgesEEC);
}

// Get the number of reconstructed jet pT bins in the unfolding response matrix
int EECCard::GetNJetPtBinsUnfoldingReco() const{
  return GetNBins(kJetPtBinEdgesUnfoldingReco);
}

// Get the number of generator level jet pT bins in the unfolding response matrix
int EECCard::GetNJetPtBinsUnfoldingTruth() const{
  return GetNBins(kJetPtBinEdgesUnfoldingTruth);
}

/*
 * Get a bin index based on a given value.
 * If value is out of bounds, return -1
 */
int EECCard::GetBinIndex(const int index, const double value) const{
  
  // Find the number of bins in the array
  int nBins = GetNBins(index);
  
  // If the number given is smaller than the first value, it is out of bounds. Return -1 in this case.
  if(value < (*fCardEntries[index])[1]) return -1;
  
  // Find the bin in which the given value is and return it
  for(int iBin = 2; iBin <= nBins+1; iBin++){
    if(value < (*fCardEntries[index])[iBin]) return iBin-2;
  }
  
  // If a value was not found, it must be larger than anything in the array. Return -1 in this case.
  return -1;
}

// Get the bin index for a given centrality value
int EECCard::GetBinIndexCentrality(const double value) const{
  return GetBinIndex(kCentralityBinEdges,value);
}

// Get the bin index for a given track pT value
int EECCard::GetBinIndexTrackPt(const double value) const{
  return GetBinIndex(kTrackPtBinEdges,value);
}

// Get the bin index for a given jet pT value in energy-energy correlator analysis
int EECCard::GetBinIndexJetPtEEC(const double value) const{
  return GetBinIndex(kJetPtBinEdgesEEC,value);
}

// Get the bin index for a given track pT value in energy-energy correlator analysis
int EECCard::GetBinIndexTrackPtEEC(const double value) const{
  return GetBinIndex(kTrackPtBinEdgesEEC,value);
}

// Get the bin index for a given reconstructed jet pT value in the unfolding response matrix
int EECCard::GetBinIndexJetPtUnfoldingReco(const double value) const{
  return GetBinIndex(kJetPtBinEdgesUnfoldingReco,value);
}

// Get the bin index for a given generator level jet pT value in the unfolding response matrix
int EECCard::GetBinIndexJetPtUnfoldingTruth(const double value) const{
  return GetBinIndex(kJetPtBinEdgesUnfoldingTruth,value);
}

/*
 * Find if a bin with given borders exist in an array defined by internal index
 *
 *  const int index = Internal index defining the array to search the bin from
 *  const double lowBorder = Low bin border of the searched bin
 *  const double highBorder = High bin border of the searched bin
 *
 *  return: Return the index of the corresponding bin. If the bin is not found, return -1
 */
int EECCard::FindBinIndex(const int index, const double lowBorder, const double highBorder) const{
  
  // Find the number of bins in the array
  int nBins = GetNBins(index);

  // Define a small number
  double epsilon = 0.001;
 
  // Go through all the bins and check if there exist a bin that matches with in given bins
  for(int iBin = 0; iBin < nBins; iBin++){
    if(TMath::Abs(GetLowBinBorder(index, iBin) - lowBorder) < epsilon){
      if(TMath::Abs(GetHighBinBorder(index, iBin) - highBorder) < epsilon){
        return iBin;
      }
    }
  }
  
  // If a match for bin borders was not found, the given bin is not in the array. Return -1 to show that.
  return -1;
}

// Find if a centrality bin with given borders exists and return its index
int EECCard::FindBinIndexCentrality(const double lowBorder, const double highBorder) const{
  return FindBinIndex(kCentralityBinEdges,lowBorder,highBorder);
}

// Find if a centrality bin with given borders exists and return its index
int EECCard::FindBinIndexCentrality(const std::pair<double,double> binBorders) const{
  return FindBinIndex(kCentralityBinEdges,binBorders.first,binBorders.second);
}

// Find if a track pT bin with given borders exists and return its index
int EECCard::FindBinIndexTrackPt(const double lowBorder, const double highBorder) const{
  return FindBinIndex(kTrackPtBinEdges,lowBorder,highBorder);
}

// Find if a track pT bin with given borders exists and return its index
int EECCard::FindBinIndexTrackPt(const std::pair<double,double> binBorders) const{
  return FindBinIndex(kTrackPtBinEdges,binBorders.first,binBorders.second);
}

// Find if a jet pT bin in energy-energy correlator analysis with given borders exists and return its index
int EECCard::FindBinIndexJetPtEEC(const double lowBorder, const double highBorder) const{
  return FindBinIndex(kJetPtBinEdgesEEC,lowBorder,highBorder);
}

// Find if a jet pT bin in energy-energy correlator analysis with given borders exists and return its index
int EECCard::FindBinIndexJetPtEEC(const std::pair<double,double> binBorders) const{
  return FindBinIndex(kJetPtBinEdgesEEC,binBorders.first,binBorders.second);
}

// Find if a track pT bin in energy-energy correlator analysis with given borders exists and return its index
int EECCard::FindBinIndexTrackPtEEC(const double lowBorder, const double highBorder) const{
  return FindBinIndex(kTrackPtBinEdgesEEC,lowBorder,highBorder);
}

// Find if a track pT bin in energy-energy correlator analysis with given borders exists and return its index
int EECCard::FindBinIndexTrackPtEEC(const std::pair<double,double> binBorders) const{
  return FindBinIndex(kTrackPtBinEdgesEEC,binBorders.first,binBorders.second);
}

// Find if a reconstructed jet pT bin in unfolding response matrix with given borders exists and return its index
int EECCard::FindBinIndexJetPtUnfoldingReco(const double lowBorder, const double highBorder) const{
  return FindBinIndex(kJetPtBinEdgesUnfoldingReco,lowBorder,highBorder);
}

// Find if a reconstructed jet pT bin in unfolding response matrix with given borders exists and return its index
int EECCard::FindBinIndexJetPtUnfoldingReco(const std::pair<double,double> binBorders) const{
  return FindBinIndex(kJetPtBinEdgesUnfoldingReco,binBorders.first,binBorders.second);
}

// Find if a generator level jet pT bin in unfolding response matrix with given borders exists and return its index
int EECCard::FindBinIndexJetPtUnfoldingTruth(const double lowBorder, const double highBorder) const{
  return FindBinIndex(kJetPtBinEdgesUnfoldingTruth,lowBorder,highBorder);
}

// Find if a generator level jet pT bin in unfolding response matrix with given borders exists and return its index
int EECCard::FindBinIndexJetPtUnfoldingTruth(const std::pair<double,double> binBorders) const{
  return FindBinIndex(kJetPtBinEdgesUnfoldingTruth,binBorders.first,binBorders.second);
}

// Get the low border of i:th bin from internal index
double EECCard::GetLowBinBorder(const int index, const int iBin) const{
  
  // Sanity check for the input
  if(iBin < 0) return -1;
  if(iBin > GetNBins(index)) return -1;
  
  // Return the asked bin index
  if(fCardEntries[index]) return (*fCardEntries[index])[iBin+1];
  return -1;
}

// Get the low border of i:th centrality bin
double EECCard::GetLowBinBorderCentrality(const int iBin) const{
  return GetLowBinBorder(kCentralityBinEdges,iBin);
}

// Get the low border of i:th track pT bin
double EECCard::GetLowBinBorderTrackPt(const int iBin) const{
  return GetLowBinBorder(kTrackPtBinEdges,iBin);
}

// Get the low border of i:th jet pT bin in energy-energy correlator analysis
double EECCard::GetLowBinBorderJetPtEEC(const int iBin) const{
  return GetLowBinBorder(kJetPtBinEdgesEEC,iBin);
}

// Get the low border of i:th track pT bin in energy-energy correlator analysis
double EECCard::GetLowBinBorderTrackPtEEC(const int iBin) const{
  return GetLowBinBorder(kTrackPtBinEdgesEEC,iBin);
}

// Get the low border of i:th reconsturcted jet pT bin in the unfolding response matrix
double EECCard::GetLowBinBorderJetPtUnfoldingReco(const int iBin) const{
  return GetLowBinBorder(kJetPtBinEdgesUnfoldingReco,iBin);
}

// Get the low border of i:th generator level jet pT bin in the unfolding response matrix
double EECCard::GetLowBinBorderJetPtUnfoldingTruth(const int iBin) const{
  return GetLowBinBorder(kJetPtBinEdgesUnfoldingTruth,iBin);
}

// Get the high border of i:th bin from internal index
double EECCard::GetHighBinBorder(const int index, const int iBin) const{
  
  // Sanity check for the input
  if(iBin < 0) return -1;
  if(iBin >= GetNBins(index)) return -1;
  
  // Return the asked bin index
  if(fCardEntries[index]) return (*fCardEntries[index])[iBin+2];
  return -1;
}

// Get the high border of i:th centrality bin
double EECCard::GetHighBinBorderCentrality(const int iBin) const{
  return GetHighBinBorder(kCentralityBinEdges,iBin);
}

// Get the high border of i:th track pT bin
double EECCard::GetHighBinBorderTrackPt(const int iBin) const{
  return GetHighBinBorder(kTrackPtBinEdges,iBin);
}

// Get the high border of i:th jet pT bin in energy-energy correlator analysis
double EECCard::GetHighBinBorderJetPtEEC(const int iBin) const{
  return GetHighBinBorder(kJetPtBinEdgesEEC,iBin);
}

// Get the high border of i:th track pT bin in energy-energy correlator analysis
double EECCard::GetHighBinBorderTrackPtEEC(const int iBin) const{
  return GetHighBinBorder(kTrackPtBinEdgesEEC,iBin);
}

// Get the high border of i:th reconsturcted jet pT bin in the unfolding response matrix
double EECCard::GetHighBinBorderJetPtUnfoldingReco(const int iBin) const{
  return GetHighBinBorder(kJetPtBinEdgesUnfoldingReco,iBin);
}

// Get the high border of i:th generator level jet pT bin in the unfolding response matrix
double EECCard::GetHighBinBorderJetPtUnfoldingTruth(const int iBin) const{
  return GetHighBinBorder(kJetPtBinEdgesUnfoldingTruth,iBin);
}

// Get the bin borders of the i:th centrality bin
std::pair<double,double> EECCard::GetBinBordersCentrality(const int iBin) const{
  return std::make_pair(GetLowBinBorderCentrality(iBin), GetHighBinBorderCentrality(iBin)); 
}

// Get the bin borders of the i:th track pT bin
std::pair<double,double> EECCard::GetBinBordersTrackPt(const int iBin) const{
  return std::make_pair(GetLowBinBorderTrackPt(iBin), GetHighBinBorderTrackPt(iBin)); 
}

// Get the bin borders of the i:th jet pT bin in energy-energy correlator analysis
std::pair<double,double> EECCard::GetBinBordersJetPtEEC(const int iBin) const{
  return std::make_pair(GetLowBinBorderJetPtEEC(iBin), GetHighBinBorderJetPtEEC(iBin)); 
}

// Get the bin borders of the i:th track pT bin in energy-energy correlator analysis
std::pair<double,double> EECCard::GetBinBordersTrackPtEEC(const int iBin) const{
  return std::make_pair(GetLowBinBorderTrackPtEEC(iBin), GetHighBinBorderTrackPtEEC(iBin)); 
}

// Get the bin borders of the i:th reconsturcted jet pT bin in the unfolding response matrix
std::pair<double,double> EECCard::GetBinBordersJetPtUnfoldingReco(const int iBin) const{
  return std::make_pair(GetLowBinBorderJetPtUnfoldingReco(iBin), GetHighBinBorderJetPtUnfoldingReco(iBin)); 
}

// Get the bin borders of the i:th generator level jet pT bin in the unfolding response matrix
std::pair<double,double> EECCard::GetBinBordersJetPtUnfoldingTruth(const int iBin) const{
  return std::make_pair(GetLowBinBorderJetPtUnfoldingTruth(iBin), GetHighBinBorderJetPtUnfoldingTruth(iBin)); 
} 


/*
 * Add one-dimensional vector to the card
 *
 * Arguments:
 *  int entryIndex = Internal index for the added vector in entry array
 *  float entryContent = Content to be put into the vector
 */
void EECCard::AddOneDimensionalVector(int entryIndex, float entryContent){
  
  // Only allow addition to postprocessing vectors
  if(entryIndex <= kPtHatBinEdges) return;
  
  // Make a new one dimensional vector to the desired index with given content
  float contents[1] = {entryContent};
  fCardEntries[entryIndex] = new TVectorT<float>(1,1,contents);

}

/*
 * Add one-dimensional vector to the card
 *
 * Arguments:
 *  int entryIndex = Internal index for the added vector in entry array
 *  int dimension = Number of entries in the given array
 *  float *contents = Content to be put into the vector
 */
void EECCard::AddVector(int entryIndex, int dimension, double *contents){
  
  // Convert double pointer to float pointer
  float* convertedContents = new float[dimension];
  for(int i = 0; i < dimension; i++){
    convertedContents[i] = contents[i];
  }
  
  // Make a new vector to the desired index with given content
  fCardEntries[entryIndex] = new TVectorT<float>(1,dimension,convertedContents);
  
  // Delete the converted contents array
  delete[] convertedContents;
  
}

/*
 * Add file name to the card
 *
 * Arguments:
 *  int entryIndex = Internal index for the added file name in file name array
 *  TString fileName = Added file name
 */
void EECCard::AddFileName(int entryIndex, TString fileName){
  
  // Make convert the string to TObjString and add it to the file name array
  fFileNames[entryIndex] = new TObjString(fileName.Data());
  
}

/*
 * Add git hash for the projecting code used to project the histograms to the card
 *
 * Arguments:
 *  const char* gitHash = Git hash to be added for projection
 */
void EECCard::AddProjectionGitHash(const char* gitHash){
  
  // Convert the const char* to TObjString and assign it to projection git hash
  fProjectionGitHash = new TObjString(gitHash);
  
}

/*
 * Add a git hash used to process the histograms in the file to the card
 *
 * Arguments:
 *  const char* gitHash = Git hash to be added for processing
 */
void EECCard::AddProcessGitHash(const char* gitHash){
  
  // Convert the const char* to TObjString and assign it to processing git hash
  fProcessGitHash = new TObjString(gitHash);
  
}

/*
 * Add a git hash used to unfold the energy-energy correlators to the card
 *
 * Arguments:
 *  const char* gitHash = Git hash to be added for unfolding
 */
void EECCard::AddUnfoldingGitHash(const char* gitHash){
  
  // Convert the const char* to TObjString and assign it to unfolding git hash
  fUnfoldingGitHash = new TObjString(gitHash);
  
}

/*
 * Print the contents of the card to the console
 */
void EECCard::Print() const{

  const char* gitHash = "Unknown";
  if(fGitHash != NULL) gitHash = fGitHash->String().Data();
  
  std::cout<<std::endl<<"========================= EECCard =========================="<<std::endl;
  std::cout << Form("%25s","GitHash"); //print keyword
  std::cout << ": " << gitHash << std::endl;
  for(int iEntry = 0; iEntry < knEntries; iEntry++){
    if(fCardEntries[iEntry]){
      std::cout << Form("%25s",fCardEntryNames[iEntry]); //print keyword
      std::cout << " (dim = "<<fCardEntries[iEntry]->GetNoElements() << ") ";//print size of TVector
      for(int iElement = 1; iElement <= fCardEntries[iEntry]->GetNoElements(); iElement++){
        std::cout << (*fCardEntries[iEntry])[iElement] << " ";//TVector components
      }
      std::cout << std::endl;
    }
  }
  std::cout << std::endl;
  
  if(fProjectionGitHash != NULL) std::cout << "Git hash for projections: " << fProjectionGitHash->String().Data() << std::endl;
  if(fUnfoldingGitHash != NULL) std::cout << "Git hash for unfolding: " << fUnfoldingGitHash->String().Data() << std::endl;
  
  for(int iFileName = 0; iFileName < knFileNames; iFileName++){
    if(fFileNames[iFileName]){
      std::cout << "Used " << fFileNameType[iFileName] << " file: " << fFileNames[iFileName]->String().Data() << std::endl;
    }
  }
  
  if(fProcessGitHash != NULL) std::cout << "Git hash for histogram processing: " << fProcessGitHash->String().Data() << std::endl;
  
}

/*
 * Write the contents of the EECCard to a file
 */
void EECCard::Write(TDirectory* file){
  
  // Create a directory to store the card parameters
  if(!file->GetDirectory("JCard")) file->mkdir("JCard");
  file->cd("JCard");
  
  // Write the git hash to the file
  if(fGitHash) fGitHash->Write("GitHash");
  
  // Write all the vectors to the file. Not all of these exist in older versions of cards, thus check if exists before writing.
  for(int iEntry = 0; iEntry < knEntries; iEntry++){
    if(fCardEntries[iEntry])  fCardEntries[iEntry]->Write(fCardEntryNames[iEntry]);
  }
   
  // Write the git hash used for projections to the file
  if(fProjectionGitHash) fProjectionGitHash->Write("ProjectionGitHash");
  
  // Write all the data names to the file.
  for(int iFileName = 0; iFileName < knFileNames; iFileName++){
    if(fFileNames[iFileName]) fFileNames[iFileName]->Write(fFileNameSaveName[iFileName]);
  }
  
  // Write the git hash used to unfold the energy-energy correlators
  if(fUnfoldingGitHash) fUnfoldingGitHash->Write("UnfoldingGitHash");

  // Write the git hash used to process the histograms
  if(fProcessGitHash) fProcessGitHash->Write("ProcessGitHash");
  
  // Return back to the main directory
  file->cd("../");
}

/*
 * Write the git hash used for unfolding energy-energy correlators to the file
 */
void EECCard::WriteUnfoldInfo(TDirectory* file){
  
  // Go to a directory to store the card parameters
  if(!file->GetDirectory("JCard")) file->mkdir("JCard");
  file->cd("JCard");
  
  // Write the git hash used to unfold the histograms
  if(fUnfoldingGitHash) fUnfoldingGitHash->Write("UnfoldingGitHash");

  // Write the file name for the response matrix
  if(fFileNames[kResponseMatrixFile]) fFileNames[kResponseMatrixFile]->Write(fFileNameSaveName[kResponseMatrixFile]);

  // Write the unfold entries
  for(int iEntry = kFirstUnfoldedCentralityBin; iEntry <= kLastUnfoldedJetPtBin; iEntry++){
    if(fCardEntries[iEntry])  fCardEntries[iEntry]->Write(fCardEntryNames[iEntry]);
  }

  // Return back to the main directory
  file->cd("../");
  
}

/*
 * Write the git hash used for processing histograms to the file
 */
void EECCard::WriteProcessHash(TDirectory* file){
  
  // Go to a directory to store the card parameters
  if(!file->GetDirectory("JCard")) file->mkdir("JCard");
  file->cd("JCard");
  
  // Write the git hash used to process the histograms
  if(fProcessGitHash) fProcessGitHash->Write("ProcessGitHash");
  
  // Return back to the main directory
  file->cd("../");
  
}
