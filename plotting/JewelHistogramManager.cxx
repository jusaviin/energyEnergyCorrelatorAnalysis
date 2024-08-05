/*
 * Implementation of JewelHistogramManager
 */

// Own includes
#include "JewelHistogramManager.h"

/*
 * Default constructor
 */
JewelHistogramManager::JewelHistogramManager():
  fLowNormalizationDeltaR(0.008),
  fHighNormalizationDeltaR(0.39)
{  
  // Energy-energy correlator histograms
  for(int iJetPt = 0; iJetPt < kJetPtBins; iJetPt++){
    for(int iTrackPt = 0; iTrackPt < kTrackPtBins; iTrackPt++){
      for(int iEnergyWeight = 0; iEnergyWeight < kEnergyWeights; iEnergyWeight++){
        fEnergyEnergyCorrelatorPp[iJetPt][iTrackPt][iEnergyWeight] = NULL;
        for(int iCentrality = 0; iCentrality < kCentralityBins; iCentrality++){ 
          fEnergyEnergyCorrelatorPbPb[iCentrality][iJetPt][iTrackPt][iEnergyWeight] = NULL;
          fEnergyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt][iEnergyWeight] = NULL;
        } // Centrality loop
      } // Energy weight loop
    } // Track pT loop
  } // Jet pT loop
}

/*
 * Constructor with input file
 */
JewelHistogramManager::JewelHistogramManager(TString inputDirectory) :
  JewelHistogramManager()
{  
  // Load all the graphs from the given directory
  LoadHistograms(inputDirectory); 
}

/*
 * Copy constructor
 */
JewelHistogramManager::JewelHistogramManager(const JewelHistogramManager& in):
  fLowNormalizationDeltaR(in.fLowNormalizationDeltaR),
  fHighNormalizationDeltaR(in.fHighNormalizationDeltaR)
{
  // Copy constructor
    
  // Energy-energy correlator histograms
  for(int iJetPt = 0; iJetPt < kJetPtBins; iJetPt++){
    for(int iTrackPt = 0; iTrackPt < kTrackPtBins; iTrackPt++){
      for(int iEnergyWeight = 0; iEnergyWeight < kEnergyWeights; iEnergyWeight++){
        fEnergyEnergyCorrelatorPp[iJetPt][iTrackPt][iEnergyWeight] = in.fEnergyEnergyCorrelatorPp[iJetPt][iTrackPt][iEnergyWeight];
        for(int iCentrality = 0; iCentrality < kCentralityBins; iCentrality++){           
          fEnergyEnergyCorrelatorPbPb[iCentrality][iJetPt][iTrackPt][iEnergyWeight] = in.fEnergyEnergyCorrelatorPbPb[iCentrality][iJetPt][iTrackPt][iEnergyWeight];
          fEnergyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt][iEnergyWeight] = in.fEnergyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt][iEnergyWeight];
        } // Centrality loop
      } // Energy weight loop
    } // Track pT loop
  } // Jet pT loop
}

/*
 *  Load all the histograms from input directory. We assume that files for each configuration are present in the given directory.
 *
 *   Arguments: TString inputDirectory = Path to directory where .dat files are located
 */
void JewelHistogramManager::LoadHistograms(TString inputDirectory){

  // If there is a "/" sign at the end of the inputDirectory string, remove it
  if(inputDirectory.EndsWith("/")){
    inputDirectory.Remove(inputDirectory.Capacity()-2);
  }
   
  // Loop over all the files that are present in the input directory and create graphs from the information inside the files
  TFile* currentFile;
  for(int iJetPt = 0; iJetPt < kJetPtBins; iJetPt++){

    // Load the pp histograms
    currentFile = TFile::Open(Form("%s/pp_%s_100k.root", inputDirectory.Data(), fJetPtName[iJetPt]));

    for(int iTrackPt = 0; iTrackPt < kTrackPtBins; iTrackPt++){
      for(int iEnergyWeight = 0; iEnergyWeight < kEnergyWeights; iEnergyWeight++){

        // Read the histogram from the file
        fEnergyEnergyCorrelatorPp[iJetPt][iTrackPt][iEnergyWeight] = (TH1D*) currentFile->Get(Form("_%s_%s_1", fEnergyWeightName[iEnergyWeight], fTrackPtName[iTrackPt]));
        fEnergyEnergyCorrelatorPp[iJetPt][iTrackPt][iEnergyWeight]->SetName(Form("energyEnergyCorrelatorJewelPp%d%d%d", iJetPt, iTrackPt, iEnergyWeight));

        // Normalize the histogram to one
      } // Energy weight loop
    } // Track pT loop

    for(int iCentrality = 0; iCentrality < kCentralityBins; iCentrality++){

      // Load the PbPb histograms
      currentFile = TFile::Open(Form("%s/PbPb_Subtracted_%s_%s_100k.root", inputDirectory.Data(), fJetPtName[iJetPt], fCentralityName[iCentrality]));

      for(int iTrackPt = 0; iTrackPt < kTrackPtBins; iTrackPt++){
        for(int iEnergyWeight = 0; iEnergyWeight < kEnergyWeights; iEnergyWeight++){

          fEnergyEnergyCorrelatorPbPb[iCentrality][iJetPt][iTrackPt][iEnergyWeight] = (TH1D*) currentFile->Get(Form("_%s_%s_1", fEnergyWeightName[iEnergyWeight], fTrackPtName[iTrackPt]));
          fEnergyEnergyCorrelatorPbPb[iCentrality][iJetPt][iTrackPt][iEnergyWeight]->SetName(Form("energyEnergyCorrelatorJewelPbPb%d%d%d%d", iCentrality, iJetPt, iTrackPt, iEnergyWeight));

        } // Energy weight loop
      } // Track pT loop
    } // Centrality loop
  } // Jet pT loop

  NormalizeHistograms();
  
}

// Normalize all histograms and recalculate the ratios
void JewelHistogramManager::NormalizeHistograms(){

  int lowIntegralBin, highIntegralBin;

  for(int iJetPt = 0; iJetPt < kJetPtBins; iJetPt++){
    for(int iTrackPt = 0; iTrackPt < kTrackPtBins; iTrackPt++){
      for(int iEnergyWeight = 0; iEnergyWeight < kEnergyWeights; iEnergyWeight++){

        // Normalize pp distributions
        lowIntegralBin = fEnergyEnergyCorrelatorPp[iJetPt][iTrackPt][iEnergyWeight]->FindBin(fLowNormalizationDeltaR);
        highIntegralBin = fEnergyEnergyCorrelatorPp[iJetPt][iTrackPt][iEnergyWeight]->FindBin(fHighNormalizationDeltaR);
        fEnergyEnergyCorrelatorPp[iJetPt][iTrackPt][iEnergyWeight]->Scale(1.0 / fEnergyEnergyCorrelatorPp[iJetPt][iTrackPt][iEnergyWeight]->Integral(lowIntegralBin, highIntegralBin, "width"));

        for(int iCentrality = 0; iCentrality < kCentralityBins; iCentrality++){

          // Normalize PbPb distributions
          lowIntegralBin = fEnergyEnergyCorrelatorPbPb[iCentrality][iJetPt][iTrackPt][iEnergyWeight]->FindBin(fLowNormalizationDeltaR);
          highIntegralBin = fEnergyEnergyCorrelatorPbPb[iCentrality][iJetPt][iTrackPt][iEnergyWeight]->FindBin(fHighNormalizationDeltaR);
          fEnergyEnergyCorrelatorPbPb[iCentrality][iJetPt][iTrackPt][iEnergyWeight]->Scale(1.0 / fEnergyEnergyCorrelatorPbPb[iCentrality][iJetPt][iTrackPt][iEnergyWeight]->Integral(lowIntegralBin, highIntegralBin, "width"));

          // Recalculate the ratios with newly normalized histograms
          fEnergyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt][iEnergyWeight] = (TH1D*) fEnergyEnergyCorrelatorPbPb[iCentrality][iJetPt][iTrackPt][iEnergyWeight]->Clone(Form("jewelRatio%d%d%d%d", iCentrality, iJetPt, iTrackPt, iEnergyWeight));
          fEnergyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt][iEnergyWeight]->Divide(fEnergyEnergyCorrelatorPp[iJetPt][iTrackPt][iEnergyWeight]);
        }

      } // Energy weight loop
    } // Track pT loop
  } // Jet pT loop

}

// Set the normalization region for the analysis and immidiate apply the normalization
void JewelHistogramManager::SetNormalizationRegion(double lowDeltaR, double highDeltaR){
  fLowNormalizationDeltaR = lowDeltaR;
  fHighNormalizationDeltaR = highDeltaR;
  NormalizeHistograms();
}

// Getters for energy-energy correlator graphs

// Getter for energy-energy correlators from pp collisions using bin indices
TH1D* JewelHistogramManager::GetEnergyEnergyCorrelatorPp(const int iJetPt, const int iTrackPt, const int iEnergyWeight) const{

  // If any of the bins is negative, indicating the histogram does not exist, return null
  if(iJetPt < 0 || iTrackPt < 0 || iEnergyWeight < 0) return NULL;

  return fEnergyEnergyCorrelatorPp[iJetPt][iTrackPt][iEnergyWeight];
}

// Getter for energy-energy correlators from pp collisions using bin borders
TH1D* JewelHistogramManager::GetEnergyEnergyCorrelatorPp(std::pair<int,int> jetPtBin, double trackPtCut, double energyWeight) const{
  int iJetPt = FindJetPtBinIndex(jetPtBin);
  int iTrackPt = FindTrackPtBinIndex(trackPtCut);
  int iEnergyWeight = FindEnergyWeightIndex(energyWeight);
  return GetEnergyEnergyCorrelatorPp(iJetPt, iTrackPt, iEnergyWeight);
}

// Getter for energy-energy correlators from PbPb collisions using bin indices
TH1D* JewelHistogramManager::GetEnergyEnergyCorrelatorPbPb(const int iCentrality, const int iJetPt, const int iTrackPt, const int iEnergyWeight) const{

  // If any of the bins is negative, indicating the histogram does not exist, return null
  if(iCentrality < 0 || iJetPt < 0 || iTrackPt < 0 || iEnergyWeight < 0) return NULL;

  return fEnergyEnergyCorrelatorPbPb[iCentrality][iJetPt][iTrackPt][iEnergyWeight];
}

// Getter for energy-energy correlators from PbPb collisions using bin borders
TH1D* JewelHistogramManager::GetEnergyEnergyCorrelatorPbPb(std::pair<int,int> centralityBin, std::pair<int,int> jetPtBin, double trackPtCut, double energyWeight) const{
  int iCentrality = FindCentralityBinIndex(centralityBin);
  int iJetPt = FindJetPtBinIndex(jetPtBin);
  int iTrackPt = FindTrackPtBinIndex(trackPtCut);
  int iEnergyWeight = FindEnergyWeightIndex(energyWeight);
  return GetEnergyEnergyCorrelatorPbPb(iCentrality, iJetPt, iTrackPt, iEnergyWeight);
}

// Getter for PbPb to pp energy-energy correlator ratios using bin indices
TH1D* JewelHistogramManager::GetEnergyEnergyCorrelatorPbPbToPpRatio(const int iCentrality, const int iJetPt, const int iTrackPt, const int iEnergyWeight) const{

  // If any of the bins is negative, indicating the histogram does not exist, return null
  if(iCentrality < 0 || iJetPt < 0 || iTrackPt < 0 || iEnergyWeight < 0) return NULL;

  return fEnergyEnergyCorrelatorPbPbToPpRatio[iCentrality][iJetPt][iTrackPt][iEnergyWeight];

}

// Getter for PbPb to pp energy-energy correlator ratios using bin borders
TH1D* JewelHistogramManager::GetEnergyEnergyCorrelatorPbPbToPpRatio(std::pair<int,int> centralityBin, std::pair<int,int> jetPtBin, double trackPtCut, double energyWeight) const{
  int iCentrality = FindCentralityBinIndex(centralityBin);
  int iJetPt = FindJetPtBinIndex(jetPtBin);
  int iTrackPt = FindTrackPtBinIndex(trackPtCut);
  int iEnergyWeight = FindEnergyWeightIndex(energyWeight);
  return GetEnergyEnergyCorrelatorPbPbToPpRatio(iCentrality, iJetPt, iTrackPt, iEnergyWeight);
}

// Find bin indices from bin borders

// Get a centrality bin index from a given centrality bin borders
int JewelHistogramManager::FindCentralityBinIndex(std::pair<int,int> centralityBin) const{
  for(int iCentrality = 0; iCentrality < kCentralityBins; iCentrality++){
    if(centralityBin.first == fCentralityBinBorders[iCentrality].first){
      if(centralityBin.second == fCentralityBinBorders[iCentrality].second){
        return iCentrality;
      }
    }
  }

  // If matching bin is not found, return a negative number
  return -1;
}

// Get a jet pT bin index from a given jet pT bin borders
int JewelHistogramManager::FindJetPtBinIndex(std::pair<int,int> jetPtBin) const{
  for(int iJetPt = 0; iJetPt < kJetPtBins; iJetPt++){
    if(jetPtBin.first == fJetPtBinBorders[iJetPt].first){
      if(jetPtBin.second == fJetPtBinBorders[iJetPt].second){
        return iJetPt;
      }
    }
  }

  // If matching bin is not found, return a negative number
  return -1;
}

// Get a track pT bin index from a given track pT cut
int JewelHistogramManager::FindTrackPtBinIndex(double trackPtBin) const{
  for(int iTrackPt = 0; iTrackPt < kTrackPtBins; iTrackPt++){
    if(TMath::Abs(trackPtBin-fTrackPtCuts[iTrackPt]) < 0.0001){
      return iTrackPt;
    }
  }

  // If matching bin is not found, return a negative number
  return -1;
}

// Get an energy weight bin index from a given energy weight value
int JewelHistogramManager::FindEnergyWeightIndex(double energyWeight) const{
  for(int iEnergyWeight = 0; iEnergyWeight < kEnergyWeights; iEnergyWeight++){
    if(TMath::Abs(energyWeight-fEnergyWeights[iEnergyWeight]) < 0.0001){
      return iEnergyWeight;
    }
  }

  // If matching bin is not found, return a negative number
  return -1;
}
