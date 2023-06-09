#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "EECUnfoldConfiguration.h"
#include "JDrawer.h"

/*
 * Reser underflow and overflow bins for a one-dimensional histogram
 */
void removeOutOfRange(TH1D* histogramInNeedOfTrimming){
   histogramInNeedOfTrimming->SetBinContent(0, 0);
   histogramInNeedOfTrimming->SetBinContent(histogramInNeedOfTrimming->GetNbinsX()+1, 0);
}

/*
 * Reser underflow and overflow bins for a two-dimensional histogram
 */
void removeOutOfRange(TH2D* histogramInNeedOfTrimming)
{
   int NX = histogramInNeedOfTrimming->GetNbinsX();
   int NY = histogramInNeedOfTrimming->GetNbinsY();

   for(int iX = 0; iX <= NX + 1; iX++)
   {
      histogramInNeedOfTrimming->SetBinContent(iX, 0, 0);
      histogramInNeedOfTrimming->SetBinContent(iX, NY+1, 0);
   }
   for(int iY = 0; iY <= NY + 1; iY++)
   {
      histogramInNeedOfTrimming->SetBinContent(0, iY, 0);
      histogramInNeedOfTrimming->SetBinContent(NX+1, iY, 0);
   }
}

/*
 * Macro for doing the jet pT unfolding for the 
 *
 *  Arguments:
 *   TString fileName = File from which the histograms are read and to which the processed histograms are written
 *   TString outputFileName = If given, the processed histograms are written to this file. Otherwise the fileName file is updated.
 *   const int iSplit = The Monte Carlo split used to derive the response matrix
 *                      0: Use the full Monte Carlo statistics to construct the response matrix
 *                      1: Use the first half of the Monte Carlo statistics to construct the response matrix
 *                      2: Use the second half of the Monte Carlo statistics to construct the response matrix
 *   const int iSystematic = Index for unfolding systematic uncertainty study
 *                           0: Default results, no systematic study
 *                           1: Evaluate systematic uncertainties for jet pT resolution
 *                           2: Evaluate systematic uncertainties for jet energy scale
 *   const int iEnergyEnergyCorrelator = Energy-energy correlator index for the unfolded correlator. Indices are explained in EECHistogramManager.h
 */
void unfoldEEChistograms(TString dataFileName, TString outputFileName, const int iSplit = 0, const int iSystematic = 0, const int iEnergyEnergyCorrelator = EECHistogramManager::kEnergyEnergyCorrelator){

  // **********************************
  //       Open the input files
  // **********************************

  // First, open the data iput file and obtain the EECCard from the file
  TFile* dataInputFile = TFile::Open(dataFileName);

  if(dataInputFile == NULL) {
    cout << "Error! The file " << dataFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }

  EECCard* dataCard = new EECCard(dataInputFile);

  // Determine if we are dealing with pp or PbPb data
  TString collisionSystem = dataCard->GetDataType();
  bool isPbPbData = collisionSystem.Contains("PbPb");

  // Use the information from the dataCard to create an unfolding configuration
  EECUnfoldConfiguration* unfoldConfigurationProvider = new EECUnfoldConfiguration(dataCard, iSplit, iSystematic);

  // Open the response matric file based on the information on the unfolding configuration
  TFile* responseInputFile = TFile::Open(unfoldConfigurationProvider->GetResponseFileName());
  
  if(responseInputFile == NULL) {
    cout << "Error! The file " << unfoldConfigurationProvider->GetResponseFileName().Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }

  // Get the EECCard also from the response matric file
  EECCard* responseCard = new EECCard(responseInputFile);
  
  
  // **********************************************************************************
  //       Check that the centrality, jet pT and track pT bins match with the files
  // **********************************************************************************

  if(isPbPbData != (dataCard->GetDataType().Contains("PbPb"))){
    cout << "You are trying to unfold pp with PbPb unfolding histograms or vice versa!" << endl;
    cout << "Please ensure that the data and unfolding histograms are comparible!" << endl;
  }

  // Only check compatible centrality bins for PbPb
  if(isPbPbData){
    if(responseCard->GetNCentralityBins() != dataCard->GetNCentralityBins()) {
      cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
    for(int iCentrality = 0; iCentrality < responseCard->GetNCentralityBins(); iCentrality++){
      if(TMath::Abs(responseCard->GetLowBinBorderCentrality(iCentrality) - dataCard->GetLowBinBorderCentrality(iCentrality)) > 5){
        cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
        return;
      }
      if(TMath::Abs(responseCard->GetHighBinBorderCentrality(iCentrality) - dataCard->GetHighBinBorderCentrality(iCentrality)) > 5){
        cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
        return;
      }
    }
  }  // Checking compatible centrality bins

  if(responseCard->GetNTrackPtBinsEEC() != dataCard->GetNTrackPtBinsEEC()){
    cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
    return;
  }
  for(int iTrackPt = 0; iTrackPt < responseCard->GetNTrackPtBinsEEC(); iTrackPt++){
    if(TMath::Abs(responseCard->GetLowBinBorderTrackPtEEC(iTrackPt) - dataCard->GetLowBinBorderTrackPtEEC(iTrackPt)) > 0.01){
      cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
    if(TMath::Abs(responseCard->GetHighBinBorderTrackPtEEC(iTrackPt) - dataCard->GetHighBinBorderTrackPtEEC(iTrackPt)) > 0.01){
      cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
  }  // Checking the compatibility of track pT bins

  if(responseCard->GetNJetPtBinsUnfoldingReco() != dataCard->GetNJetPtBinsEEC()){
    cout << "Error! Measured jet pT bins do not match between the two files! Cannot execute the code!" << endl;
    return;
  }
  for(int iJetPt = 0; iJetPt < responseCard->GetNJetPtBinsUnfoldingReco(); iJetPt++){
    if(TMath::Abs(responseCard->GetLowBinBorderJetPtUnfoldingReco(iJetPt) - dataCard->GetLowBinBorderJetPtEEC(iJetPt)) > 0.01){
      cout << "Error! Measured jet pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
    if(TMath::Abs(responseCard->GetHighBinBorderJetPtUnfoldingReco(iJetPt) - dataCard->GetHighBinBorderJetPtEEC(iJetPt)) > 0.01){
      cout << "Error! Measured jet pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
  }  // Checking the compatibility of reconstructed jet pT bins

  // ********************************************************
  //       Binning configuration for the unfolding study
  // ********************************************************
  
  // Find the number of bins from the card
  const int nCentralityBins = (isPbPbData) ? responseCard->GetNCentralityBins() : 1;
  const int nTrackPtBins = responseCard->GetNTrackPtBinsEEC();
  const int nJetPtBins = dataCard->GetNJetPtBinsEEC();

  // Bin range to be studied
  int firstStudiedCentralityBin = 0;
  int lastStudiedCentralityBin = nCentralityBins-1;
  
  int firstStudiedTrackPtBinEEC = 3;
  int lastStudiedTrackPtBinEEC = 5;

  // Select explicitly the jet pT bins that we are going to unfold
  std::vector<std::pair<double,double>> unfoldedJetPtBins;
  unfoldedJetPtBins.push_back(std::make_pair(120,140));
  unfoldedJetPtBins.push_back(std::make_pair(140,160));
  unfoldedJetPtBins.push_back(std::make_pair(160,180));
  unfoldedJetPtBins.push_back(std::make_pair(180,200));

  const int nUnfoldedJetPtBins = unfoldedJetPtBins.size();

  bool saveFigures = false;
  TString saveComment = "_matrixInversion";
  TString figureFormat = "pdf";
    
  // ***************************************************************
  //    Create histogram managers and load the needed histograms
  // ***************************************************************

  // Load the histograms needed for unfolding from the unfolding histogram manager
  EECHistogramManager* responseHistograms = new EECHistogramManager(responseInputFile, responseCard);
  responseHistograms->SetLoadJetPtUnfoldingHistograms(true);
  responseHistograms->SetCentralityBinRange(firstStudiedCentralityBin,lastStudiedCentralityBin);
  responseHistograms->SetTrackPtBinRangeEEC(firstStudiedTrackPtBinEEC,lastStudiedTrackPtBinEEC);
  responseHistograms->LoadProcessedHistograms();

  // Load the data histograms to be unfolded from the data histogram manager
  EECHistogramManager* dataHistograms = new EECHistogramManager(dataInputFile, dataCard);
  dataHistograms->SetLoadEnergyEnergyCorrelators(iEnergyEnergyCorrelator == EECHistogramManager::kEnergyEnergyCorrelator);
  dataHistograms->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationPlus(iEnergyEnergyCorrelator == EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationPlus);
  dataHistograms->SetLoadEnergyEnergyCorrelatorsEfficiencyVariationMinus(iEnergyEnergyCorrelator == EECHistogramManager::kEnergyEnergyCorrelatorEfficiencyVariationMinus);
  dataHistograms->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationPlus(iEnergyEnergyCorrelator == EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationPlus);
  dataHistograms->SetLoadEnergyEnergyCorrelatorsPairEfficiencyVariationMinus(iEnergyEnergyCorrelator == EECHistogramManager::kEnergyEnergyCorrelatorPairEfficiencyVariationMinus);
  dataHistograms->SetCentralityBinRange(firstStudiedCentralityBin,lastStudiedCentralityBin);
  dataHistograms->SetTrackPtBinRangeEEC(firstStudiedTrackPtBinEEC,lastStudiedTrackPtBinEEC);
  dataHistograms->SetJetPtBinRangeEEC(0, dataCard->GetNJetPtBinsEEC());
  dataHistograms->LoadProcessedHistograms();

  // Histograms that are needed to create the unfolding response
  TH1D* hUnfoldingMeasured[nCentralityBins][nTrackPtBins];
  TH1D* hUnfoldingTruth[nCentralityBins][nTrackPtBins];
  TH2D* hUnfoldingResponse[nCentralityBins][nTrackPtBins];
  TH1D* hUnfoldedDistribution[nCentralityBins][nTrackPtBins];
  TH1D* energyEnergyCorrelatorForUnfolding[nCentralityBins][nTrackPtBins];
  TH1D* energyEnergyCorrelatorsFromData[nCentralityBins][nJetPtBins][nTrackPtBins];

  // RooUnfold response object
  RooUnfoldResponse* rooResponse[nCentralityBins][nTrackPtBins];

  // Initialize all the histograms to null
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
      hUnfoldingMeasured[iCentrality][iTrackPt] = NULL;
      hUnfoldingTruth[iCentrality][iTrackPt] = NULL;
      hUnfoldingResponse[iCentrality][iTrackPt] = NULL;
      energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt] = NULL;
      rooResponse[iCentrality][iTrackPt] = NULL;
      hUnfoldedDistribution[iCentrality][iTrackPt] = NULL;
      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){
        energyEnergyCorrelatorsFromData[iCentrality][iJetPt][iTrackPt] = NULL;
      }
    } // Track pT loop
  }  // Centrality loop

  // Read the histograms needed for the unfolding response and create the RooUnfold response objects
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      hUnfoldingMeasured[iCentrality][iTrackPt] = responseHistograms->GetHistogramJetPtUnfoldingMeasured(iCentrality, iTrackPt);
      hUnfoldingTruth[iCentrality][iTrackPt] = responseHistograms->GetHistogramJetPtUnfoldingTruth(iCentrality, iTrackPt);
      hUnfoldingResponse[iCentrality][iTrackPt] = responseHistograms->GetHistogramJetPtUnfoldingResponse(iCentrality, iTrackPt);

      // Revome the overflow and underflow bins from the histograms
      removeOutOfRange(hUnfoldingMeasured[iCentrality][iTrackPt]);
      removeOutOfRange(hUnfoldingTruth[iCentrality][iTrackPt]);
      removeOutOfRange(hUnfoldingResponse[iCentrality][iTrackPt]);

      // Create the response object with trimmed histograms
      rooResponse[iCentrality][iTrackPt] = new RooUnfoldResponse(hUnfoldingMeasured[iCentrality][iTrackPt], hUnfoldingTruth[iCentrality][iTrackPt], hUnfoldingResponse[iCentrality][iTrackPt]);

      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){
        energyEnergyCorrelatorsFromData[iCentrality][iJetPt][iTrackPt] = dataHistograms->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt);
      } // Measured jet pT loop
    } // Track pT loop
  } // Centrality loop

  // Next, we need to transform the data histograms into a format that can be read by RooUnfold
  // For this, we will need to combine the jet pT and deltaR axes
  int nDeltaRBinsData = energyEnergyCorrelatorsFromData[firstStudiedCentralityBin][0][firstStudiedTrackPtBinEEC]->GetNbinsX();
  double jetPtLowerBound;
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt] = (TH1D*) hUnfoldingMeasured[iCentrality][iTrackPt]->Clone(Form("dataCorrelatorForUnfolding%d%d", iCentrality, iTrackPt));
      energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt]->Reset();
      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){
        jetPtLowerBound = dataCard->GetLowBinBorderJetPtEEC(iJetPt);
        for(int iBin = 1; iBin <= nDeltaRBinsData; iBin++){
          // Check what happens if the measured jet pT is capped at 80
          if(jetPtLowerBound < 80){
            energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt]->SetBinContent(iBin + nDeltaRBinsData*iJetPt, 0);
            energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt]->SetBinError(iBin + nDeltaRBinsData*iJetPt, 0);
          } else {
            energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt]->SetBinContent(iBin + nDeltaRBinsData*iJetPt, energyEnergyCorrelatorsFromData[iCentrality][iJetPt][iTrackPt]->GetBinContent(iBin) * energyEnergyCorrelatorsFromData[iCentrality][iJetPt][iTrackPt]->GetBinWidth(iBin));
            energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt]->SetBinError(iBin + nDeltaRBinsData*iJetPt, energyEnergyCorrelatorsFromData[iCentrality][iJetPt][iTrackPt]->GetBinError(iBin) * energyEnergyCorrelatorsFromData[iCentrality][iJetPt][iTrackPt]->GetBinWidth(iBin));
          }
        } // DeltaR bin loop
      } // Jet pT loop
    } // Track pT loop
  }  // Centrality loop

  // ********************************************************
  //                     The actual unfolding
  // ********************************************************

  // Unfolding using Bayesian unfolding
  RooUnfoldBayes* bayesUnfold[nCentralityBins][nTrackPtBins];
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      bayesUnfold[iCentrality][iTrackPt] = new RooUnfoldBayes(rooResponse[iCentrality][iTrackPt], energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt], unfoldConfigurationProvider->GetNumberOfIterations(dataCard->GetBinBordersCentrality(iCentrality), dataCard->GetLowBinBorderTrackPtEEC(iTrackPt)));
      hUnfoldedDistribution[iCentrality][iTrackPt] = (TH1D*)bayesUnfold[iCentrality][iTrackPt]->Hunfold();
    } // Track pT loop
  }  // Centrality loop


  // *********************************************************************
  //     Transforming unfolded histograms back to one dimensional ones
  // *********************************************************************

  // Dissect the unfolded histograms in jet pT bins from the big histograms
  TH1D* hUnfolded[nCentralityBins][nUnfoldedJetPtBins][nTrackPtBins];

  // Initialize the dissected histograms to null
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < nUnfoldedJetPtBins; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
        hUnfolded[iCentrality][iJetPt][iTrackPt] = NULL;
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

  // Find the DeltaR binning within a single jet pT bin from the big histogram
  double previousBinWidth = 0;
  double currentBinWidth = 0;
  int nDeltaRBins = 0;
  for(int iBin = 1; iBin < hUnfoldingTruth[firstStudiedCentralityBin][firstStudiedTrackPtBinEEC]->GetNbinsX(); iBin++){
    currentBinWidth = hUnfoldingTruth[firstStudiedCentralityBin][firstStudiedTrackPtBinEEC]->GetBinWidth(iBin);
    if(currentBinWidth < previousBinWidth) break;
    previousBinWidth = currentBinWidth;
    nDeltaRBins++;
  }

  double deltaRBinning[nDeltaRBins+1];
  for(int iBin = 1; iBin <= nDeltaRBins+1; iBin++){
    deltaRBinning[iBin-1] = hUnfoldingTruth[firstStudiedCentralityBin][firstStudiedTrackPtBinEEC]->GetBinLowEdge(iBin);
  }

  // Create new histograms for the chosen bin range
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      for(int iJetPt = 0; iJetPt < nUnfoldedJetPtBins; iJetPt++){
        hUnfolded[iCentrality][iJetPt][iTrackPt] = new TH1D(Form("hUnfolded%d%d%d", iCentrality, iTrackPt, iJetPt), Form("hUnfolded%d%d%d", iCentrality, iTrackPt, iJetPt), nDeltaRBins, deltaRBinning);
      } // Jet pT loop
    } // Track pT loop
  } // Centrality loop

  // Dissect the big histograms and fill small histograms based on the information
  // Also normalize the bins to bin width
  int jetPtUnfoldIndex = 0;
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      for(int iJetPt = 0; iJetPt < nUnfoldedJetPtBins; iJetPt++){
        jetPtUnfoldIndex = dataCard->FindBinIndexJetPtEEC(unfoldedJetPtBins.at(iJetPt));

        for(int iBin = 1; iBin <= nDeltaRBins; iBin++){
          hUnfolded[iCentrality][iJetPt][iTrackPt]->SetBinContent(iBin, hUnfoldedDistribution[iCentrality][iTrackPt]->GetBinContent(iBin + nDeltaRBins*jetPtUnfoldIndex) / hUnfoldedDistribution[iCentrality][iTrackPt]->GetBinWidth(iBin));
          hUnfolded[iCentrality][iJetPt][iTrackPt]->SetBinError(iBin, hUnfoldedDistribution[iCentrality][iTrackPt]->GetBinError(iBin + nDeltaRBins*jetPtUnfoldIndex) / hUnfoldedDistribution[iCentrality][iTrackPt]->GetBinWidth(iBin));
        } // DeltaR bin loop

      } // Jet pT loop
    } // Track pT loop
  } // Centrality loop
  
  // *********************************************************************
  //                Save the unfolded histograms to a file
  // *********************************************************************
  
  // First, add the git hash used to unfold the histograms to the card
  int firstUnfoldedJetPtBin = dataCard->FindBinIndexJetPtEEC(unfoldedJetPtBins.at(0));
  int lastUnfoldedJetPtBin = dataCard->FindBinIndexJetPtEEC(unfoldedJetPtBins.at(nUnfoldedJetPtBins-1));
  const char* gitHash = "GITHASHHERE";
  dataCard->AddUnfoldingGitHash(gitHash);
  dataCard->AddFileName(EECCard::kResponseMatrixFile, unfoldConfigurationProvider->GetResponseFileName());
  dataCard->AddOneDimensionalVector(EECCard::kFirstUnfoldedCentralityBin, firstStudiedCentralityBin);
  dataCard->AddOneDimensionalVector(EECCard::kLastUnfoldedCentralityBin, lastStudiedCentralityBin);
  dataCard->AddOneDimensionalVector(EECCard::kFirstUnfoldedTrackPtBin, firstStudiedTrackPtBinEEC);
  dataCard->AddOneDimensionalVector(EECCard::kLastUnfoldedTrackPtBin, lastStudiedTrackPtBinEEC);
  dataCard->AddOneDimensionalVector(EECCard::kFirstUnfoldedJetPtBin, firstUnfoldedJetPtBin);
  dataCard->AddOneDimensionalVector(EECCard::kLastUnfoldedJetPtBin, lastUnfoldedJetPtBin);


  EECHistogramManager* histogramsSavedToFile = new EECHistogramManager(dataCard);
  histogramsSavedToFile->SetLoadEnergyEnergyCorrelators(true);
  histogramsSavedToFile->SetCentralityBinRange(firstStudiedCentralityBin,lastStudiedCentralityBin);
  histogramsSavedToFile->SetTrackPtBinRangeEEC(firstStudiedTrackPtBinEEC,lastStudiedTrackPtBinEEC);
  histogramsSavedToFile->SetJetPtBinRangeEEC(0, dataCard->GetNJetPtBinsEEC());

  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      for(int iJetPt = 0; iJetPt < nUnfoldedJetPtBins; iJetPt++){
        jetPtUnfoldIndex = dataCard->FindBinIndexJetPtEEC(unfoldedJetPtBins.at(iJetPt));
        histogramsSavedToFile->SetUnfoldedEnergyEnergyCorrelator(hUnfolded[iCentrality][iJetPt][iTrackPt], EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, jetPtUnfoldIndex, iTrackPt);
      } // Jet pT loop
    } // Track pT loop
  } // Centrality loop

  histogramsSavedToFile->WriteUnfoldedEnergyEnergyCorrelators(outputFileName, "UPDATE");
}
