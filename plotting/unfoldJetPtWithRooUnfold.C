#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Do the jet pT unfolding using RooUnfold
 */
void unfoldJetPtWithRooUnfold(){

  // **********************************
  //       Open the input files
  // **********************************

  // Define the name for the file containing histograms needed for unfolding
  TString unfoldingInputFileName = "data/ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_firstUnfoldingTest_noBinWidthNormalized_processed_2023-05-08.root";
  // ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_firstUnfoldingTest_processed_2023-05-08.root
  // ppMC2017_RecoReco_Pythia8_pfJets_wtaAxis_firstUnfoldingTest_noBinWidthNormalized_processed_2023-05-08.root

  // Name of the file containing the data that needs to be unfolded
  TString dataInputFileName = "data/ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_newBinning_wayMoreJetPtBins_processed_2023-04-13.root";

  // Open the input files
  TFile* unfoldingInputFile = TFile::Open(unfoldingInputFileName);
  TFile* dataInputFile = TFile::Open(dataInputFileName);

  if(unfoldingInputFile == NULL) {
    cout << "Error! The file " << unfoldingInputFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }

  if(dataInputFile == NULL) {
    cout << "Error! The file " << dataInputFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }

  // For binning information, read the card from the unfolding file
  EECCard* unfoldingCard = new EECCard(unfoldingInputFile);
  EECCard* dataCard = new EECCard(dataInputFile);

  // Determine if we are dealing with pp or PbPb data
  TString collisionSystem = unfoldingCard->GetDataType();
  bool isPbPbData = collisionSystem.Contains("PbPb");
  
  // *************************************************************************************
  //       Check that the centrality, jet and track pT bins match with the two files
  // *************************************************************************************

  bool dataSystemIsPbPb = dataCard->GetDataType().Contains("PbPb");
  if(isPbPbData != dataSystemIsPbPb){
    cout << "You are trying to unfold pp with PbPb unfolding histograms or vice versa!" << endl;
    cout << "Please ensure that the data and unfolding histograms are comparible!" << endl;
  }

  // Only check compatible centrality bins for PbPb
  if(isPbPbData){
    if(unfoldingCard->GetNCentralityBins() != dataCard->GetNCentralityBins()){
      cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
    for(int iCentrality = 0; iCentrality < unfoldingCard->GetNCentralityBins(); iCentrality++){
      if(TMath::Abs(unfoldingCard->GetLowBinBorderCentrality(iCentrality) - dataCard->GetLowBinBorderCentrality(iCentrality)) > 5){
        cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
        return;
      }
      if(TMath::Abs(unfoldingCard->GetHighBinBorderCentrality(iCentrality) - dataCard->GetHighBinBorderCentrality(iCentrality)) > 5){
        cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
        return;
      }
    }
  } // Checking compatible centrality bins

  if(unfoldingCard->GetNTrackPtBinsEEC() != dataCard->GetNTrackPtBinsEEC()){
    cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
    return;
  }
  for(int iTrackPt = 0; iTrackPt < unfoldingCard->GetNTrackPtBinsEEC(); iTrackPt++){
    if(TMath::Abs(unfoldingCard->GetLowBinBorderTrackPtEEC(iTrackPt) - dataCard->GetLowBinBorderTrackPtEEC(iTrackPt)) > 0.01){
      cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
    if(TMath::Abs(unfoldingCard->GetHighBinBorderTrackPtEEC(iTrackPt) - dataCard->GetHighBinBorderTrackPtEEC(iTrackPt)) > 0.01){
      cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
  } // Checking the compatibility of track pT bins

  if(unfoldingCard->GetNJetPtBinsEEC() != dataCard->GetNJetPtBinsEEC()){
    cout << "Error! Jet pT bins do not match between the two files! Cannot execute the code!" << endl;
    return;
  }
  for(int iJetPt = 0; iJetPt < unfoldingCard->GetNJetPtBinsEEC(); iJetPt++){
    if(TMath::Abs(unfoldingCard->GetLowBinBorderJetPtEEC(iJetPt) - dataCard->GetLowBinBorderJetPtEEC(iJetPt)) > 0.01){
      cout << "Error! Jet pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
    if(TMath::Abs(unfoldingCard->GetHighBinBorderJetPtEEC(iJetPt) - dataCard->GetHighBinBorderJetPtEEC(iJetPt)) > 0.01){
      cout << "Error! Jet pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
  } // Checking the compatibility of jet pT bins

  // ********************************************************
  //       Binning configuration for the unfolding study
  // ********************************************************
  
  // Find the number of bins from the card
  const int nCentralityBins = (isPbPbData) ? unfoldingCard->GetNCentralityBins() : 1;
  const int nJetPtBinsEEC = unfoldingCard->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = unfoldingCard->GetNTrackPtBinsEEC();

  // Bin range to be studied
  int firstStudiedCentralityBin = 0;
  int lastStudiedCentralityBin = 0;
  
  int firstStudiedTrackPtBinEEC = 1;
  int lastStudiedTrackPtBinEEC = 5;

  int firstStudiedJetPtBinEEC = 4;
  int lastStudiedJetPtBinEEC = 4;

  // Default binning ranges for reference
  // centrality = {0,10,30,50,90}

  bool saveFigures = false;
  TString saveComment = "_firstLook";
  TString figureFormat = "png";
    
  // ***************************************************************
  //    Create histogram managers and load the needed histograms
  // ***************************************************************

  // Load the histograms needed for unfolding from the unfolding histogram manager
  EECHistogramManager* unfoldingHistograms = new EECHistogramManager(unfoldingInputFile, unfoldingCard);
  unfoldingHistograms->SetLoadJetPtUnfoldingHistograms(true);
  unfoldingHistograms->SetCentralityBinRange(firstStudiedCentralityBin,lastStudiedCentralityBin);
  unfoldingHistograms->SetTrackPtBinRangeEEC(firstStudiedTrackPtBinEEC,lastStudiedTrackPtBinEEC);
  unfoldingHistograms->LoadProcessedHistograms();

  // Load the data histograms to be unfolded from the data histogram manager
  EECHistogramManager* dataHistograms = new EECHistogramManager(dataInputFile, dataCard);
  dataHistograms->SetLoadEnergyEnergyCorrelators(true);
  dataHistograms->SetCentralityBinRange(firstStudiedCentralityBin,lastStudiedCentralityBin);
  dataHistograms->SetTrackPtBinRangeEEC(firstStudiedTrackPtBinEEC,lastStudiedTrackPtBinEEC);
  dataHistograms->SetJetPtBinRangeEEC(0, dataCard->GetNJetPtBinsEEC());
  dataHistograms->LoadProcessedHistograms();

  // Histograms that are needed to create the unfolding response
  TH1D* hUnfoldingMeasured[nCentralityBins][nTrackPtBinsEEC];
  TH1D* hUnfoldingTruth[nCentralityBins][nTrackPtBinsEEC];
  TH2D* hUnfoldingResponse[nCentralityBins][nTrackPtBinsEEC];
  TH1D* hUnfoldedDistribution[nCentralityBins][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorForUnfolding[nCentralityBins][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorsFromData[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];

  // RooUnfold response object
  RooUnfoldResponse* rooResponse[nCentralityBins][nTrackPtBinsEEC];

  // Initialize all the histograms to null
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
      hUnfoldingMeasured[iCentrality][iTrackPt] = NULL;
      hUnfoldingTruth[iCentrality][iTrackPt] = NULL;
      hUnfoldingResponse[iCentrality][iTrackPt] = NULL;
      hUnfoldedDistribution[iCentrality][iTrackPt] = NULL;
      energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt] = NULL;
      rooResponse[iCentrality][iTrackPt] = NULL;
      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
        energyEnergyCorrelatorsFromData[iCentrality][iJetPt][iTrackPt] = NULL;
      }
    } // Track pT loop
  }  // Centrality loop

  // Read the histograms needed for the unfolding response and create the RooUnfold response objects
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      hUnfoldingMeasured[iCentrality][iTrackPt] = unfoldingHistograms->GetHistogramJetPtUnfoldingMeasured(iCentrality, iTrackPt);
      hUnfoldingTruth[iCentrality][iTrackPt] = unfoldingHistograms->GetHistogramJetPtUnfoldingTruth(iCentrality, iTrackPt);
      hUnfoldingResponse[iCentrality][iTrackPt] = unfoldingHistograms->GetHistogramJetPtUnfoldingResponse(iCentrality, iTrackPt);
      rooResponse[iCentrality][iTrackPt] = new RooUnfoldResponse(hUnfoldingMeasured[iCentrality][iTrackPt], hUnfoldingTruth[iCentrality][iTrackPt], hUnfoldingResponse[iCentrality][iTrackPt]);

      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
        energyEnergyCorrelatorsFromData[iCentrality][iJetPt][iTrackPt] = dataHistograms->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt);
      } // Jet pT loop
    } // Track pT loop
  }  // Centrality loop

  // Next, we need to transform the data histograms into a format that can be read by RooUnfold
  // For this, we will need to combine the jet pT and deltaR axes
  int nDeltaRBinsData = energyEnergyCorrelatorsFromData[firstStudiedCentralityBin][0][firstStudiedTrackPtBinEEC]->GetNbinsX();
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt] = (TH1D*) hUnfoldingMeasured[iCentrality][iTrackPt]->Clone(Form("dataCorrelatorForUnfolding%d%d", iCentrality, iTrackPt));
      energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt]->Reset();
      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
        for(int iBin = 1; iBin <= nDeltaRBinsData; iBin++){
          energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt]->SetBinContent(iBin + nDeltaRBinsData*iJetPt, energyEnergyCorrelatorsFromData[iCentrality][iJetPt][iTrackPt]->GetBinContent(iBin) * energyEnergyCorrelatorsFromData[iCentrality][iJetPt][iTrackPt]->GetBinWidth(iBin));
          energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt]->SetBinError(iBin + nDeltaRBinsData*iJetPt, energyEnergyCorrelatorsFromData[iCentrality][iJetPt][iTrackPt]->GetBinError(iBin) * energyEnergyCorrelatorsFromData[iCentrality][iJetPt][iTrackPt]->GetBinWidth(iBin));
        } // DeltaR bin loop
      } // Jet pT loop
    } // Track pT loop
  }  // Centrality loop


  // After the response it setup, as a first check, we can try to unfold the measured distribution
  RooUnfoldInvert* invertUnfold[nCentralityBins][nTrackPtBinsEEC];
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      invertUnfold[iCentrality][iTrackPt] = new RooUnfoldInvert(rooResponse[iCentrality][iTrackPt], energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt]);
      hUnfoldedDistribution[iCentrality][iTrackPt] = (TH1D*)invertUnfold[iCentrality][iTrackPt]->Hunfold();
    } // Track pT loop
  }  // Centrality loop

  // Dissect the big histograms back into small histograms that can be compared
  TH1D* hMeasured[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hTruth[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hUnfolded[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hMeasuredToTruthRatio[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hUnfoldedToTruthRatio[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];

  // Initialize the dissected histograms to null
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        hMeasured[iCentrality][iJetPt][iTrackPt] = NULL;
        hTruth[iCentrality][iJetPt][iTrackPt] = NULL;
        hUnfolded[iCentrality][iJetPt][iTrackPt] = NULL;
        hMeasuredToTruthRatio[iCentrality][iJetPt][iTrackPt] = NULL;
        hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt] = NULL;
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

  // Find the DeltaR binning within a single jet pT bin from the big histogram
  double previousBinWidth = 0;
  double currentBinWidth = 0;
  int nDeltaRBins = 0;
  for(int iBin = 1; iBin < hUnfoldingMeasured[firstStudiedCentralityBin][firstStudiedTrackPtBinEEC]->GetNbinsX(); iBin++){
    currentBinWidth = hUnfoldingMeasured[firstStudiedCentralityBin][firstStudiedTrackPtBinEEC]->GetBinWidth(iBin);
    if(currentBinWidth < previousBinWidth) break;
    previousBinWidth = currentBinWidth;
    nDeltaRBins++;
  }

  double deltaRBinning[nDeltaRBins+1];
  for(int iBin = 1; iBin <= nDeltaRBins+1; iBin++){
    deltaRBinning[iBin-1] = hUnfoldingMeasured[firstStudiedCentralityBin][firstStudiedTrackPtBinEEC]->GetBinLowEdge(iBin);
  }

  // Create new histograms for the chosen bin range
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
        hMeasured[iCentrality][iJetPt][iTrackPt] = new TH1D(Form("hMeasured%d%d%d", iCentrality, iTrackPt, iJetPt), Form("hMeasured%d%d%d", iCentrality, iTrackPt, iJetPt), nDeltaRBins, deltaRBinning);
        hTruth[iCentrality][iJetPt][iTrackPt] = new TH1D(Form("hTruth%d%d%d", iCentrality, iTrackPt, iJetPt), Form("hTruth%d%d%d", iCentrality, iTrackPt, iJetPt), nDeltaRBins, deltaRBinning);
        hUnfolded[iCentrality][iJetPt][iTrackPt] = new TH1D(Form("hUnfolded%d%d%d", iCentrality, iTrackPt, iJetPt), Form("hUnfolded%d%d%d", iCentrality, iTrackPt, iJetPt), nDeltaRBins, deltaRBinning);
      } // Jet pT loop
    } // Track pT loop
  } // Centrality loop

  // Dissect the big histograms and fill small histograms based on the information
  // Also normalize the bins to bin width and histograms to one in the defined range
  double normalizationRegionLow = 0.006;
  double normalizationRegionHigh = 0.4;
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
        for(int iBin = 1; iBin <= nDeltaRBins; iBin++){
          hMeasured[iCentrality][iJetPt][iTrackPt]->SetBinContent(iBin, hUnfoldingMeasured[iCentrality][iTrackPt]->GetBinContent(iBin + nDeltaRBins*iJetPt) / hUnfoldingMeasured[iCentrality][iTrackPt]->GetBinWidth(iBin));
          hMeasured[iCentrality][iJetPt][iTrackPt]->SetBinError(iBin, hUnfoldingMeasured[iCentrality][iTrackPt]->GetBinError(iBin + nDeltaRBins*iJetPt) / hUnfoldingMeasured[iCentrality][iTrackPt]->GetBinWidth(iBin));
          
          hTruth[iCentrality][iJetPt][iTrackPt]->SetBinContent(iBin, hUnfoldingTruth[iCentrality][iTrackPt]->GetBinContent(iBin + nDeltaRBins*iJetPt) / hUnfoldingTruth[iCentrality][iTrackPt]->GetBinWidth(iBin));
          hTruth[iCentrality][iJetPt][iTrackPt]->SetBinError(iBin, hUnfoldingTruth[iCentrality][iTrackPt]->GetBinError(iBin + nDeltaRBins*iJetPt) / hUnfoldingTruth[iCentrality][iTrackPt]->GetBinWidth(iBin));
          
          hUnfolded[iCentrality][iJetPt][iTrackPt]->SetBinContent(iBin, hUnfoldedDistribution[iCentrality][iTrackPt]->GetBinContent(iBin + nDeltaRBins*iJetPt) / hUnfoldedDistribution[iCentrality][iTrackPt]->GetBinWidth(iBin));
          hUnfolded[iCentrality][iJetPt][iTrackPt]->SetBinError(iBin, hUnfoldedDistribution[iCentrality][iTrackPt]->GetBinError(iBin + nDeltaRBins*iJetPt) / hUnfoldedDistribution[iCentrality][iTrackPt]->GetBinWidth(iBin));
        } // DeltaR bin loop

        // Do the normalization of total distribution to one
        hMeasured[iCentrality][iJetPt][iTrackPt]->Scale(1.0 / hMeasured[iCentrality][iJetPt][iTrackPt]->Integral(hMeasured[iCentrality][iJetPt][iTrackPt]->FindBin(normalizationRegionLow), hMeasured[iCentrality][iJetPt][iTrackPt]->FindBin(normalizationRegionHigh), "width"));
        hTruth[iCentrality][iJetPt][iTrackPt]->Scale(1.0 / hTruth[iCentrality][iJetPt][iTrackPt]->Integral(hTruth[iCentrality][iJetPt][iTrackPt]->FindBin(normalizationRegionLow), hTruth[iCentrality][iJetPt][iTrackPt]->FindBin(normalizationRegionHigh), "width"));
        hUnfolded[iCentrality][iJetPt][iTrackPt]->Scale(1.0 / hUnfolded[iCentrality][iJetPt][iTrackPt]->Integral(hUnfolded[iCentrality][iJetPt][iTrackPt]->FindBin(normalizationRegionLow), hUnfolded[iCentrality][iJetPt][iTrackPt]->FindBin(normalizationRegionHigh), "width"));
      } // Jet pT loop
    } // Track pT loop
  } // Centrality loop

  // Calculate the ratios of the other histograms to truth
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
        hMeasuredToTruthRatio[iCentrality][iJetPt][iTrackPt] = (TH1D*) hMeasured[iCentrality][iJetPt][iTrackPt]->Clone(Form("measuredToTruthRatio%d%d%d", iCentrality, iTrackPt, iJetPt));
        hMeasuredToTruthRatio[iCentrality][iJetPt][iTrackPt]->Divide(hTruth[iCentrality][iJetPt][iTrackPt]);

        hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt] = (TH1D*) hUnfolded[iCentrality][iJetPt][iTrackPt]->Clone(Form("unfoldedToTruthRatio%d%d%d", iCentrality, iTrackPt, iJetPt));
        hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt]->Divide(hTruth[iCentrality][iJetPt][iTrackPt]);
      } // Jet pT loop
    } // Track pT loop
  } // Centrality loop
  
  // ****************************************************************
  //         Check the unfolded distribution against truth
  // ****************************************************************

  // Prepare a JDrawer for drawing purposes
  JDrawer* drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);

  drawer->SetLogX(true);
  
  TLine* oneLine = new TLine(normalizationRegionLow, 1, normalizationRegionHigh, 1);
  oneLine->SetLineStyle(2);

  // Helper variables
  TLegend* legend;
  TString centralityString;
  TString compactCentralityString;
  TString jetPtString;
  TString compactJetPtString;
  TString trackPtString;
  TString compactTrackPtString;

  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){

    // Set the centrality information for legends and figure saving
    if(isPbPbData) {
      centralityString = Form("Pythia+Hydjet: %.0f-%.0f", unfoldingCard->GetLowBinBorderCentrality(iCentrality), unfoldingCard->GetHighBinBorderCentrality(iCentrality));
      compactCentralityString = Form("_C%.0f-%.0f", unfoldingCard->GetLowBinBorderCentrality(iCentrality), unfoldingCard->GetHighBinBorderCentrality(iCentrality));
    } else {
      centralityString = "Pythia8";
      compactCentralityString = "_pythia8";
    }

    for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
      // Set the jet pT information for legends and figure saving
      jetPtString = Form("%.0f < jet p_{T} < %.0f", unfoldingHistograms->GetJetPtBinBorderEEC(iJetPt), unfoldingHistograms->GetJetPtBinBorderEEC(iJetPt + 1));
      compactJetPtString = Form("_J=%.0f-%.0f", unfoldingHistograms->GetJetPtBinBorderEEC(iJetPt), unfoldingHistograms->GetJetPtBinBorderEEC(iJetPt + 1));

      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
        // Create a new canvas for the plot
        drawer->CreateSplitCanvas();

        // Set the track pT information for legends and figure saving
        trackPtString = Form("%.1f < track p_{T}", unfoldingHistograms->GetTrackPtBinBorderEEC(iTrackPt));
        compactTrackPtString = Form("_T%.0f", unfoldingHistograms->GetTrackPtBinBorderEEC(iTrackPt));

        // Draw first the generator level distribution
        hTruth[iCentrality][iJetPt][iTrackPt]->SetLineColor(kBlack);
        hTruth[iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(normalizationRegionLow, normalizationRegionHigh);
        drawer->SetLogY(true);
        drawer->DrawHistogramToUpperPad(hTruth[iCentrality][iJetPt][iTrackPt], "#Deltar", "EEC", " ", "");

        // Add the reconstructed and unfolded distributions to the same plot
        hMeasured[iCentrality][iJetPt][iTrackPt]->SetLineColor(kRed);
        hMeasured[iCentrality][iJetPt][iTrackPt]->Draw("same");
        hUnfolded[iCentrality][iJetPt][iTrackPt]->SetLineColor(kBlue);
        hUnfolded[iCentrality][iJetPt][iTrackPt]->Draw("same");

        // Add a legend to the figure
        legend = new TLegend(0.25, 0.15, 0.5, 0.5);
        legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);

        legend->AddEntry((TObject*)0, centralityString.Data(), "");
        legend->AddEntry((TObject*)0, jetPtString.Data(), "");
        legend->AddEntry((TObject*)0, trackPtString.Data(), "");

        legend->AddEntry(hTruth[iCentrality][iJetPt][iTrackPt], "Generator level reference", "l");
        legend->AddEntry(hMeasured[iCentrality][iJetPt][iTrackPt], "Reconstructed correlator", "l");
        legend->AddEntry(hUnfolded[iCentrality][iJetPt][iTrackPt], "Unfolded correlator", "l");

        legend->Draw();

        // Draw the ratios to lower pad
        drawer->SetLogY(false);
        hMeasuredToTruthRatio[iCentrality][iJetPt][iTrackPt]->SetLineColor(kRed);
        hMeasuredToTruthRatio[iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.8, 1.2);
        hMeasuredToTruthRatio[iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(normalizationRegionLow, normalizationRegionHigh);
        drawer->DrawHistogramToLowerPad(hMeasuredToTruthRatio[iCentrality][iJetPt][iTrackPt], "#Deltar", "Ratio to truth", " ", "");

        hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt]->SetLineColor(kBlue);
        hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt]->Draw("same");

        oneLine->Draw();

        // Save the figures to a file
        if(saveFigures) {
          gPad->GetCanvas()->SaveAs(Form("figures/jetPtUnfoldingWithRooUnfold%s%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
        }

      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
}
