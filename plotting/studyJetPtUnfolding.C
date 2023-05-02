#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Macro to study the performance of the STAR method for jet pT unfolding
 */
void studyJetPtUnfolding(){

  // Enumeration for the file types needed in the study
  enum enumDistributionTypes{kGeneratorLevel, kReconstructed, kResponseMatrix, knDistributionTypes};

  // **********************************
  //       Open the input files
  // **********************************

  // Define the files names to be opened
  TString inputFileName[knDistributionTypes];

  // File with the generator level distributions
  inputFileName[kGeneratorLevel] = "data/PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_newDeltaRBins_wayMoreJetPtBins_matchJets_processed_2023-04-19.root";
  // ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_newBinning_wayMoreJetPtBins_processed_2023-04-13.root
  // PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_newDeltaRBins_wayMoreJetPtBins_processed_2023-04-17.root
  // PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_newDeltaRBins_wayMoreJetPtBins_matchJets_processed_2023-04-18.root

  // File with reconstructed information
  inputFileName[kReconstructed] = "data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_newDeltaRBins_wayMoreJetPtBins_matchJets_processed_2023-04-19.root";
  // ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_newBinning_wayMoreJetPtBins_processed_2023-04-13.root
  // PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_newDeltaRBins_wayMoreJetPtBins_processed_2023-04-17.root
  // PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_newDeltaRBins_wayMoreJetPtBins_matchJets_processed_2023-04-17.root

  // Jet pT response matrix to be used for unfolding the reconstructed information
  inputFileName[kResponseMatrix] = "data/PbPbMC2018_RecoGen_akFlowJets_miniAOD_4pCentShift_noTrigger_jetPtClosure_newRatioBins_processed_2023-04-21.root";
  // ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_noCorrelations_jetPtResponseMatrixMoreBins_processed_2023-01-13.root
  // PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_jetPtResponseMatrix_finalMcWeight_processed_2023-03-06.root

  // Open the input files
  TFile* inputFile[knDistributionTypes];
  EECCard* card[knDistributionTypes]; 
  for(int iFile = 0; iFile < knDistributionTypes; iFile++){
    inputFile[iFile] = TFile::Open(inputFileName[iFile]);

    if(inputFile[iFile] == NULL) {
      cout << "Error! The file " << inputFileName[iFile].Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }

    // Load the cards from the file
    card[iFile] = new EECCard(inputFile[iFile]);

  } // File loop

  // DOUBLE FOLDING!!!
  // Start with constant reco pT bin
  // Look at which gen pT bins correcpond to this
  // For each gen pT bin, look at which reco pT bins correspond to these
  // Construct the reconstructed pT back from the reconstructed pT:s
  // This should correspond to the correct generator level distribution corresponding to the original jet pT selection!

  // Determine if we are dealing with pp or PbPb data
  TString collisionSystem = card[0]->GetDataType();
  bool isPbPbData = collisionSystem.Contains("PbPb");
  
  // ********************************************************
  //       Binning configuration for the unfolding study
  // ********************************************************
  
  // Find the number of bins from the card
  const int nCentralityBins = (isPbPbData) ? card[0]->GetNCentralityBins() : 1;
  const int nJetPtBinsEEC = card[kReconstructed]->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card[kReconstructed]->GetNTrackPtBinsEEC();

  // Bin range to be studied
  int firstStudiedCentralityBin = 0;
  int lastStudiedCentralityBin = 0;
  
  int firstStudiedJetPtBinEEC = card[kGeneratorLevel]->GetBinIndexJetPtEEC(130); // First studied bin is 120-140 in the true jet pT
  int lastStudiedJetPtBinEEC = card[kGeneratorLevel]->GetBinIndexJetPtEEC(190);  // First studied bin is 180-200 in the true jet pT
  
  int firstStudiedTrackPtBinEEC = 5;
  int lastStudiedTrackPtBinEEC = 5;

  const int nAnalysisJetPtBins = lastStudiedJetPtBinEEC - firstStudiedJetPtBinEEC + 1;

  // Default binning ranges for reference
  // centrality = {0,10,30,50,90}

  const bool drawUnfoldedToGeneratorLevelComparison = false;
  const bool drawFoldedToReconstructedComparison = true;

  bool saveFigures = false;
  TString saveComment = "_firstLook";
  TString figureFormat = "png";
    
  // ***************************************************************
  //    Create histogram managers and load the needed histograms
  // ***************************************************************

  EECHistogramManager* histograms[knDistributionTypes];

  // From the generator level file, load the energy-energy correlator histograms from the selected range
  histograms[kGeneratorLevel] = new EECHistogramManager(inputFile[kGeneratorLevel],card[kGeneratorLevel]);
  histograms[kGeneratorLevel]->SetLoadEnergyEnergyCorrelators(true);
  histograms[kGeneratorLevel]->SetCentralityBinRange(firstStudiedCentralityBin,lastStudiedCentralityBin);
  histograms[kGeneratorLevel]->SetJetPtBinRangeEEC(0,card[kGeneratorLevel]->GetNJetPtBinsEEC()-1);
  histograms[kGeneratorLevel]->SetTrackPtBinRangeEEC(firstStudiedTrackPtBinEEC,lastStudiedTrackPtBinEEC);
  histograms[kGeneratorLevel]->LoadProcessedHistograms();

  // From the reconstructed file, load the energy-energy correlators from all available jet pT bins 
  histograms[kReconstructed] = new EECHistogramManager(inputFile[kReconstructed],card[kReconstructed]);
  histograms[kReconstructed]->SetLoadEnergyEnergyCorrelators(true);
  histograms[kReconstructed]->SetCentralityBinRange(firstStudiedCentralityBin,lastStudiedCentralityBin);
  histograms[kReconstructed]->SetJetPtBinRangeEEC(0,card[kReconstructed]->GetNJetPtBinsEEC()-1);
  histograms[kReconstructed]->SetTrackPtBinRangeEEC(firstStudiedTrackPtBinEEC,lastStudiedTrackPtBinEEC);
  histograms[kReconstructed]->LoadProcessedHistograms();

  // Read all the response matrices from the response matrix file
  histograms[kResponseMatrix] = new EECHistogramManager(inputFile[kResponseMatrix],card[kResponseMatrix]);
  histograms[kResponseMatrix]->SetLoadJetPtResponseMatrix(true);
  histograms[kResponseMatrix]->SetCentralityBinRange(firstStudiedCentralityBin,lastStudiedCentralityBin);
  histograms[kResponseMatrix]->LoadProcessedHistograms();

  // Histograms for generator level, reconstructed, and unfolded energy-energy correlators
  TH1D* hGeneratorLevelEnergyEnergyCorrelator[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hReconstructedEnergyEnergyCorrelator[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hUnfoldedEnergyEnergyCorrelator[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hFoldedGeneratorLevelEnergyEnergyCorrelator[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];

  // Ratio histograms with respect to the generator level
  TH1D* hReconstructedToGeneratorLevelRatio[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hUnfoldedToGeneratorLevelRatio[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hFoldedToReconstructedRatio[nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];

  // Histograms for response matrix and projections for it
  TH2D* hJetPtResponseMatrix[nCentralityBins];
  TH1D* hJetPtResponseMatrixProjection[nCentralityBins][nJetPtBinsEEC];
  TH1D* hJetPtResponseMatrixProjectionForFolding[nCentralityBins][nJetPtBinsEEC];

  // Initialize all the histograms to null
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
        hGeneratorLevelEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt] = NULL; 
        hReconstructedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt] = NULL;
        hUnfoldedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt] = NULL;
        hFoldedGeneratorLevelEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt] = NULL;
        hReconstructedToGeneratorLevelRatio[iCentrality][iJetPt][iTrackPt] = NULL;
        hUnfoldedToGeneratorLevelRatio[iCentrality][iJetPt][iTrackPt] = NULL;
        hFoldedToReconstructedRatio[iCentrality][iJetPt][iTrackPt] = NULL;
      } // Full jet pT loop
    } // Track pT loop
    
    hJetPtResponseMatrix[iCentrality] = NULL;

    for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
      hJetPtResponseMatrixProjection[iCentrality][iJetPt] = NULL;
      hJetPtResponseMatrixProjectionForFolding[iCentrality][iJetPt] = NULL;
    } // Analysis jet pT loop

  }  // Centrality loop

  // Read the histograms from the files
  int projectedBin;
  double normalizationRegionLow = 0.006;   // Low range of the normalization region
  double normalizationRegionHigh = 0.4;    // High range of the normalization region
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
        hGeneratorLevelEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt] = histograms[kGeneratorLevel]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt);

        // Normalize the distribution to one within the normalization region
        hGeneratorLevelEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Scale(1.0 / hGeneratorLevelEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Integral(hGeneratorLevelEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->FindBin(normalizationRegionLow), hGeneratorLevelEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->FindBin(normalizationRegionHigh), "width"));


      } // Analysis jet pT loop
      for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
        hReconstructedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt] = histograms[kReconstructed]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt);

        // Normalize the distribution to one within the normalization region
        hReconstructedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Scale(1.0 / hReconstructedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Integral(hReconstructedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->FindBin(normalizationRegionLow), hReconstructedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->FindBin(normalizationRegionHigh), "width"));

      } // Full jet pT loop
    } // Track pT loop

    // Read the two dimensional response matrices
    hJetPtResponseMatrix[iCentrality] = histograms[kResponseMatrix]->GetHistogramJetPtResponseMatrix(iCentrality);

    // Do the projections from the two-dimensional response matrices
    for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){
      projectedBin = hJetPtResponseMatrix[iCentrality]->GetYaxis()->FindBin(card[kGeneratorLevel]->GetLowBinBorderJetPtEEC(iJetPt)+1);
      hJetPtResponseMatrixProjection[iCentrality][iJetPt] = hJetPtResponseMatrix[iCentrality]->ProjectionX(Form("constantGenProjection%d%d", iCentrality, iJetPt), projectedBin, projectedBin);

      // Normilize the projection to one to make it probability distribution
      hJetPtResponseMatrixProjection[iCentrality][iJetPt]->Scale(1.0 / hJetPtResponseMatrixProjection[iCentrality][iJetPt]->Integral());

      // Do the projection also with respect to the other axis
      projectedBin = hJetPtResponseMatrix[iCentrality]->GetXaxis()->FindBin(card[kReconstructed]->GetLowBinBorderJetPtEEC(iJetPt)+1);
      hJetPtResponseMatrixProjectionForFolding[iCentrality][iJetPt] = hJetPtResponseMatrix[iCentrality]->ProjectionY(Form("constantRecoProjection%d%d", iCentrality, iJetPt), projectedBin, projectedBin);

      // Normilize the projection to one to make it probability distribution
      hJetPtResponseMatrixProjectionForFolding[iCentrality][iJetPt]->Scale(1.0 / hJetPtResponseMatrixProjectionForFolding[iCentrality][iJetPt]->Integral());
    } // Jet pT loop 

  }  // Centrality loop

  // ***************************************************************************************************
  //      Do the jet pT unfolding using the STAR method and test folding of generator level results
  // ***************************************************************************************************
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){

        // First, get an empty histogram to start with
        hUnfoldedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt] = (TH1D*) hGeneratorLevelEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Clone(Form("unfoldedEEC%d%d%d", iCentrality, iTrackPt, iJetPt));
        hUnfoldedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Reset();

        // TODO: Some check that the energy-energy correlator and response matrix jet pT bins match

        // Unfold the distribution as a probability weighted sum of the reconstructed distributions
        for(int iUnfold = 1; iUnfold <= hJetPtResponseMatrixProjection[iCentrality][iJetPt]->GetNbinsX(); iUnfold++){
          if(hJetPtResponseMatrixProjection[iCentrality][iJetPt]->GetBinContent(iUnfold) < 0.001) continue;
          hUnfoldedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Add(hReconstructedEnergyEnergyCorrelator[iCentrality][iUnfold-1][iTrackPt], hJetPtResponseMatrixProjection[iCentrality][iJetPt]->GetBinContent(iUnfold));
        }

        // Calculate the ratios to generator level
        hReconstructedToGeneratorLevelRatio[iCentrality][iJetPt][iTrackPt] = (TH1D*) hReconstructedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Clone(Form("reconstructedToGeneratorLevelRatio%d%d%d", iCentrality, iJetPt, iTrackPt));
        hReconstructedToGeneratorLevelRatio[iCentrality][iJetPt][iTrackPt]->Divide(hGeneratorLevelEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]);

        hUnfoldedToGeneratorLevelRatio[iCentrality][iJetPt][iTrackPt] = (TH1D*) hUnfoldedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Clone(Form("unfoldedToGeneratorLevelRatio%d%d%d", iCentrality, iJetPt, iTrackPt));
        hUnfoldedToGeneratorLevelRatio[iCentrality][iJetPt][iTrackPt]->Divide(hGeneratorLevelEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]);

        // Get an empty histogram to get also folding started
        hFoldedGeneratorLevelEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt] = (TH1D*) hReconstructedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Clone(Form("foldedGeneratorLevelEEC%d%d%d", iCentrality, iTrackPt, iJetPt));
        hFoldedGeneratorLevelEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Reset();

        // Fold the generator level distributions with a probability weighted sum of the generator level distributions
        for(int iUnfold = 1; iUnfold <= hJetPtResponseMatrixProjectionForFolding[iCentrality][iJetPt]->GetNbinsX(); iUnfold++){
          if(hJetPtResponseMatrixProjectionForFolding[iCentrality][iJetPt]->GetBinContent(iUnfold) < 0.001) continue;
          hFoldedGeneratorLevelEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Add(hGeneratorLevelEnergyEnergyCorrelator[iCentrality][iUnfold-1][iTrackPt], hJetPtResponseMatrixProjectionForFolding[iCentrality][iJetPt]->GetBinContent(iUnfold));
        }

        // Calculate the ratios to reconstructed
        hFoldedToReconstructedRatio[iCentrality][iJetPt][iTrackPt] = (TH1D*) hFoldedGeneratorLevelEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Clone(Form("foldedToReconstructedRatio%d%d%d", iCentrality, iJetPt, iTrackPt));
        hFoldedToReconstructedRatio[iCentrality][iJetPt][iTrackPt]->Divide(hReconstructedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]);

      } // Jet pT loop
    } // Track pT loop
  } // Centrality loop
  
  // **********************************
  //         Draw the figures
  // **********************************

  // Prepare a JDrawer for drawing purposes
  JDrawer* drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);

  drawer->DrawHistogram(hJetPtResponseMatrixProjection[0][firstStudiedJetPtBinEEC], "pT", "counts", " ");

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

  // Compare the true distribution to reconstructed and unfolded distributions
  if(drawUnfoldedToGeneratorLevelComparison){

    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){

      // Set the centrality information for legends and figure saving
      if(isPbPbData){
        centralityString = Form("Pythia+Hydjet: %.0f-%.0f", card[kGeneratorLevel]->GetLowBinBorderCentrality(iCentrality), card[kGeneratorLevel]->GetHighBinBorderCentrality(iCentrality));
        compactCentralityString = Form("_C%.0f-%.0f", card[kGeneratorLevel]->GetLowBinBorderCentrality(iCentrality), card[kGeneratorLevel]->GetHighBinBorderCentrality(iCentrality));
      } else {
        centralityString = "Pythia8";
        compactCentralityString = "_pythia8";
      }

      for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){

        // Set the jet pT information for legends and figure saving
        jetPtString = Form("%.0f < jet p_{T} < %.0f", histograms[kGeneratorLevel]->GetJetPtBinBorderEEC(iJetPt), histograms[kGeneratorLevel]->GetJetPtBinBorderEEC(iJetPt + 1));
        compactJetPtString = Form("_J=%.0f-%.0f", histograms[kGeneratorLevel]->GetJetPtBinBorderEEC(iJetPt), histograms[kGeneratorLevel]->GetJetPtBinBorderEEC(iJetPt + 1));

        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++) {

          // Create a new canvas for the plot
          drawer->CreateSplitCanvas();

          // Set the track pT information for legends and figure saving
          trackPtString = Form("%.1f < track p_{T}", histograms[kGeneratorLevel]->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString = Form("_T%.0f", histograms[kGeneratorLevel]->GetTrackPtBinBorderEEC(iTrackPt));

          // Draw first the generator level distribution
          hGeneratorLevelEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->SetLineColor(kBlack);
          hGeneratorLevelEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(normalizationRegionLow, normalizationRegionHigh);
          drawer->SetLogY(true);
          drawer->DrawHistogramToUpperPad(hGeneratorLevelEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt], "#Deltar", "EEC", " ", "");

          // Add the reconstructed and unfolded distributions to the same plot
          hReconstructedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->SetLineColor(kRed);
          hReconstructedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Draw("same");
          hUnfoldedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->SetLineColor(kBlue);
          hUnfoldedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Draw("same");

          // Add a legend to the figure
          legend = new TLegend(0.25, 0.15, 0.5, 0.5);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);

          legend->AddEntry((TObject*)0, centralityString.Data(), "");
          legend->AddEntry((TObject*)0, jetPtString.Data(), "");
          legend->AddEntry((TObject*)0, trackPtString.Data(), "");

          legend->AddEntry(hGeneratorLevelEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt], "Generator level reference", "l");
          legend->AddEntry(hReconstructedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt], "Reconstructed correlator", "l");
          legend->AddEntry(hUnfoldedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt], "Unfolded correlator", "l");

          legend->Draw();

          // Draw the ratios to lower pad
          drawer->SetLogY(false);
          hReconstructedToGeneratorLevelRatio[iCentrality][iJetPt][iTrackPt]->SetLineColor(kRed);
          hReconstructedToGeneratorLevelRatio[iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.8, 1.2);
          hReconstructedToGeneratorLevelRatio[iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(normalizationRegionLow, normalizationRegionHigh);
          drawer->DrawHistogramToLowerPad(hReconstructedToGeneratorLevelRatio[iCentrality][iJetPt][iTrackPt], "#Deltar", "Ratio to gen", " ", "");

          hUnfoldedToGeneratorLevelRatio[iCentrality][iJetPt][iTrackPt]->SetLineColor(kBlue);
          hUnfoldedToGeneratorLevelRatio[iCentrality][iJetPt][iTrackPt]->Draw("same");

          oneLine->Draw();

          // Save the figures to a file
          if(saveFigures) {
            gPad->GetCanvas()->SaveAs(Form("figures/jetPtUnfolding%s%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
          }

        }  // Track pT loop
      }    // Jet pT loop
    }      // Centrality loop
  }

  // Compare the reconstructed distribution to folded generator level distribution
  if(drawFoldedToReconstructedComparison) {

    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){

      // Set the centrality information for legends and figure saving
      if(isPbPbData){
        centralityString = Form("Pythia+Hydjet: %.0f-%.0f", card[kReconstructed]->GetLowBinBorderCentrality(iCentrality), card[kReconstructed]->GetHighBinBorderCentrality(iCentrality));
        compactCentralityString = Form("_C%.0f-%.0f", card[kReconstructed]->GetLowBinBorderCentrality(iCentrality), card[kReconstructed]->GetHighBinBorderCentrality(iCentrality));
      } else {
        centralityString = "Pythia8";
        compactCentralityString = "_pythia8";
      }

      for(int iJetPt = firstStudiedJetPtBinEEC; iJetPt <= lastStudiedJetPtBinEEC; iJetPt++){

        // Set the jet pT information for legends and figure saving
        jetPtString = Form("%.0f < jet p_{T} < %.0f", histograms[kReconstructed]->GetJetPtBinBorderEEC(iJetPt), histograms[kReconstructed]->GetJetPtBinBorderEEC(iJetPt + 1));
        compactJetPtString = Form("_J=%.0f-%.0f", histograms[kReconstructed]->GetJetPtBinBorderEEC(iJetPt), histograms[kReconstructed]->GetJetPtBinBorderEEC(iJetPt + 1));

        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

          // Create a new canvas for the plot
          drawer->CreateSplitCanvas();

          // Set the track pT information for legends and figure saving
          trackPtString = Form("%.1f < track p_{T}", histograms[kReconstructed]->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString = Form("_T%.0f", histograms[kReconstructed]->GetTrackPtBinBorderEEC(iTrackPt));

          // Draw first the generator level distribution
          hReconstructedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->SetLineColor(kBlack);
          hReconstructedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(normalizationRegionLow, normalizationRegionHigh);
          drawer->SetLogY(true);
          drawer->DrawHistogramToUpperPad(hReconstructedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt], "#Deltar", "EEC", " ", "");

          // Add the reconstructed and unfolded distributions to the same plot
          hFoldedGeneratorLevelEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->SetLineColor(kRed);
          hFoldedGeneratorLevelEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt]->Draw("same");

          // Add a legend to the figure
          legend = new TLegend(0.25, 0.15, 0.5, 0.5);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);

          legend->AddEntry((TObject*)0, centralityString.Data(), "");
          legend->AddEntry((TObject*)0, jetPtString.Data(), "");
          legend->AddEntry((TObject*)0, trackPtString.Data(), "");

          legend->AddEntry(hReconstructedEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt], "Reconstructed correlator", "l");
          legend->AddEntry(hFoldedGeneratorLevelEnergyEnergyCorrelator[iCentrality][iJetPt][iTrackPt], "Folded correlator", "l");

          legend->Draw();

          // Draw the ratios to lower pad
          drawer->SetLogY(false);
          hFoldedToReconstructedRatio[iCentrality][iJetPt][iTrackPt]->SetLineColor(kRed);
          hFoldedToReconstructedRatio[iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.8, 1.2);
          hFoldedToReconstructedRatio[iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(normalizationRegionLow, normalizationRegionHigh);
          drawer->DrawHistogramToLowerPad(hFoldedToReconstructedRatio[iCentrality][iJetPt][iTrackPt], "#Deltar", "Ratio to gen", " ", "");

          oneLine->Draw();

          // Save the figures to a file
          if(saveFigures) {
            gPad->GetCanvas()->SaveAs(Form("figures/jetPtFolding%s%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
          }

        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Draw folded to reconstructed comparison
}
