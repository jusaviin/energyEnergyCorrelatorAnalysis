#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
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
 * Forward fold the generator level distribution using the response matrix to get the detector level distribution
 *
 * Arguments:
 *  TH1D* unfoldedHistogram = Histogram with unfolded results. This will be forward folded
 *  TH2D* responseMatrix = Response matrix from which the forward folding is calculated
 *
 *  return: Forward folded histogram
 */
TH1D* forwardFold(TH1D* unfoldedHistogram, TH2D* responseMatrix){

  // Cannot perform forward folding is the histograms do not exist
  if(unfoldedHistogram == NULL || responseMatrix == NULL) return NULL;

  // Find the number of bins in generator and reconstructed levels
  int nBinsGen = responseMatrix->GetNbinsY();
  int nBinsReco = responseMatrix->GetNbinsX();

  // Create a histogram to hold the forward folded distribution
  TH1D* foldedHistogram = (TH1D*)unfoldedHistogram->Clone(Form("%sFolded", unfoldedHistogram->GetName()));
  foldedHistogram->Reset();

  // Do the actual forward folding
  double recoSum;
  double recoFractionOfTotal, recoFractionOfTotalError;
  double truthContent, truthError;
  double valueUpdate, errorUpdate;

  // Loop over all generator level bins
  for(int iGenBin = 1; iGenBin <= nBinsGen; iGenBin++){

    // Sum all the reconstructed values corresponding to given generator level value
    recoSum = 0;
    for(int iRecoBin = 1; iRecoBin <= nBinsReco; iRecoBin++){
      recoSum = recoSum + responseMatrix->GetBinContent(iRecoBin, iGenBin);
    }

    // If there are no reconstructed entries, also the entry in the folded distribution will be 0
    if(recoSum == 0) continue;

    // Loop over all the reconstructed bins corresponding to the generator level bin again and update their values
    for(int iRecoBin = 1; iRecoBin <= nBinsReco; iRecoBin++){

      // Fraction of the total reconstructed value that is found from this generator level bin
      recoFractionOfTotal = responseMatrix->GetBinContent(iRecoBin, iGenBin) / recoSum;

      // Content of the generator level bin
      truthContent = unfoldedHistogram->GetBinContent(iGenBin);

      // Relative uncertianties of the two values above
      recoFractionOfTotalError = responseMatrix->GetBinError(iRecoBin, iGenBin) / responseMatrix->GetBinContent(iRecoBin, iGenBin);
      truthError = unfoldedHistogram->GetBinError(iGenBin) / unfoldedHistogram->GetBinContent(iGenBin);

      // If the bin content in matrix is 0, ensure that also error for this is 0
      if(responseMatrix->GetBinContent(iRecoBin, iGenBin) == 0) recoFractionOfTotalError = 0;

      // Take the fraction of the truth as an input for the folded histogram
      valueUpdate = recoFractionOfTotal * truthContent;
      errorUpdate = TMath::Sqrt(recoFractionOfTotalError * recoFractionOfTotalError + truthError * truthError) * valueUpdate;

      foldedHistogram->SetBinContent(iRecoBin, foldedHistogram->GetBinContent(iRecoBin) + valueUpdate);
      foldedHistogram->SetBinError(iRecoBin, TMath::Sqrt(foldedHistogram->GetBinError(iRecoBin)*foldedHistogram->GetBinError(iRecoBin) + errorUpdate*errorUpdate));
    } // Reconstructed level bin loop
  } // Generator level bin loop

  return foldedHistogram;
}

/*
 * Do the jet pT unfolding using RooUnfold for first 20 iterations and determine the chi2 of the truth-to-unfolded fit
 */
void determineNumberOfUnfoldingIterations(int iSplit = 1, int iSystematic = 0){

  // Enumeration variables
  enum enumEnergyEnergyCorrelatorFileType{kDataFile, kTruthReferenceFile, kNFileTypes};
  enum enumSystematicVariations{kNominal, kUncertaintySmearDown, kUncertaintySmearUp, kJECShiftDown, kJECShiftUp, kPriorShape, kNSystematicVariations};

  TString unfoldName = "Bayes, NIT iterations";

  // **********************************
  //       Open the input files
  // **********************************

  //int iSplit = 1; // Choose 1 or 2
  //int iSystematic = kNominal; // Systematic uncertainty index. 
  TString systematicName[kNSystematicVariations] = {"nominalSmear", "uncertaintySmearDown", "uncertaintySmearUp", "minusJECuncertainty", "plusJECuncertainty", "nominalSmear_jetPtWeight"};

  // Define the name for the file containing histograms needed for unfolding
  TString unfoldingInputFileName = Form("data/PbPbMC2018_GenGen_akFlowJets_4pCentShift_cutBadPhi_energyWeightSquared_nominalReflectedCone_%s_responseMatrix_part%d_processed_2023-12-02.root", systematicName[iSystematic].Data(), iSplit);
  // ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_32deltaRBins_responseMatrix_part%d_processed_2023-06-02.root
  // ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_32deltaRBins_smearJetPtResolution_responseMatrix_part%d_processed_2023-06-02.root
  // ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_32deltaRBins_smearJetPtUncertainty_responseMatrix_part%d_processed_2023-06-02.root
  // ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_32deltaRBins_jetPtWeight_responseMatrix_part%d_processed_2023-06-15.root
  // PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_moreLowPtBins_unfoldingHistograms_part%d_processed_2023-05-20.root
  // PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_noJetPtWeight_smearJetPtResolution_responseMatrix_part%d_processed_2023-06-02.root
  // PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_noJetPtWeight_smearJetPtUncertainty_responseMatrix_part%d_processed_2023-06-02.root
  // PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_jetPtWeight_responseMatrix_part%d_processed_2023-06-15.root

  // New files with JetMet jet energy resolution smearing factors
  // ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_32deltaRBins_%s_responseMatrix_part%d_processed_2023-06-21.root
  // ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_%s_responseMatrix_part%d_processed_2023-10-30.root
  // ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_%s_responseMatrix_part%d_processed_2023-11-08.root
  // PbPbMC2018_GenGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_%s_responseMatrix_part%d_processed_2023-06-23.root
  // PbPbMC2018_GenGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_%s_responseMatrix_part%d_processed_2023-07-06.root
  // PbPbMC2018_GenGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_energyWeightSquared_%s_responseMatrix_part%d_processed_2023-10-24.root

  // Name of the file containing the data that needs to be unfolded
  TString energyEnergyCorrelatorInputFileName[kNFileTypes];
  energyEnergyCorrelatorInputFileName[kDataFile] = Form("data/PbPbMC2018_RecoGen_akFlowJets_4pCentShift_cutBadPhi_energyWeightSquared_nominalReflectedCone_%s_reconstructedReference_part%d_processed_2023-12-02.root", systematicName[iSystematic].Data(), 3-iSplit);
  // ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_32deltaRBins_reconstructedReferenceForUnfolding_part%d_processed_2023-06-05.root
  // ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_32deltaRBins_jetPtWeight_reconstructedReference_part%d_processed_2023-06-15.root
  // PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_moreLowPtBins_reconstructedReferenceForUnfolding_part%d_processed_2023-05-20.root
  // PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_jetPtWeight_reconstructedReference_part%d_processed_2023-06-15.root

  // New files with JetMet jet energy resolution smearing factors
  // ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_32deltaRBins_%s_reconstructedReference_part%d_processed_2023-06-21.root
  // ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_%s_reconstructedReference_part%d_processed_2023-10-30.root
  // ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_%s_reconstructedReference_part%d_processed_2023-11-08.root
  // PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_%s_reconstructedReference_part%d_processed_2023-06-23.root
  // PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_%s_reconstructedReference_part%d_processed_2023-07-06.root
  // PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_energyWeightSquared_%s_reconstructedReference_part%d_processed_2023-10-24.root

  energyEnergyCorrelatorInputFileName[kTruthReferenceFile] = Form("data/PbPbMC2018_GenGen_akFlowJets_4pCentShift_cutBadPhi_energyWeightSquared_nominalReflectedCone_%s_truthReference_part%d_processed_2023-12-02.root", systematicName[iSystematic].Data(), 3-iSplit);
  // ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_32deltaRBins_truthReferenceForUnfolding_part%d_processed_2023-06-05.root
  // ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_32deltaRBins_jetPtWeight_truthReference_part%d_processed_2023-06-15.root
  // PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_moreLowPtBins_truthReferenceForUnfolding_part%d_processed_2023-05-20.root
  // PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_jetPtWeight_truthReference_part%d_processed_2023-06-15.root

  // New files with JetMet jet energy resolution smearing factors
  // ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_32deltaRBins_%s_truthReference_part%d_processed_2023-06-21.root
  // ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_%s_truthReference_part%d_processed_2023-10-30.root
  // ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_energyWeightSquared_%s_truthReference_part%d_processed_2023-11-08.root
  // PbPbMC2018_GenGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_%s_truthReference_part%d_processed_2023-06-23.root
  // PbPbMC2018_GenGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_%s_truthReference_part%d_processed_2023-07-06.root
  // PbPbMC2018_GenGen_eecAnalysis_akFlowJets_4pCentShift_cutBadPhi_energyWeightSquared_%s_truthReference_part%d_processed_2023-10-24.root

  // Option to ignore truth reference file. We might just want to do the regular unfolding without comparing results to truth.
  const bool ignoreTruthReferenceFile = false; 

  // Open the input files
  TFile* unfoldingInputFile = TFile::Open(unfoldingInputFileName);
  
  if(unfoldingInputFile == NULL) {
    cout << "Error! The file " << unfoldingInputFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }

  // For the binning information, read the card from the energy-energy correlator file
  EECCard* unfoldingCard = new EECCard(unfoldingInputFile);

  TFile* energyEnergyCorrelatorInputFile[kNFileTypes];
  EECCard* energyenergyCorrelatorCard[kNFileTypes]; 
  for(int iFileType = 0; iFileType < kNFileTypes-ignoreTruthReferenceFile; iFileType++){

    energyEnergyCorrelatorInputFile[iFileType] = TFile::Open(energyEnergyCorrelatorInputFileName[iFileType]);

    if(energyEnergyCorrelatorInputFile[iFileType] == NULL) {
      cout << "Error! The file " << energyEnergyCorrelatorInputFileName[iFileType].Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }

    // For binning information, read the card from the energy-energy correlator file
    energyenergyCorrelatorCard[iFileType] = new EECCard(energyEnergyCorrelatorInputFile[iFileType]);
  }

  // Determine if we are dealing with pp or PbPb data
  TString collisionSystem = unfoldingCard->GetDataType();
  bool isPbPbData = collisionSystem.Contains("PbPb");
  
  // **********************************************************************************
  //       Check that the centrality, jet pT and track pT bins match with the files
  // **********************************************************************************

  for(int iFileType = 0; iFileType < kNFileTypes - ignoreTruthReferenceFile; iFileType++) {

    if(isPbPbData != (energyenergyCorrelatorCard[iFileType]->GetDataType().Contains("PbPb"))) {
      cout << "You are trying to unfold pp with PbPb unfolding histograms or vice versa!" << endl;
      cout << "Please ensure that the data and unfolding histograms are comparible!" << endl;
    }

    // Only check compatible centrality bins for PbPb
    if(isPbPbData) {
      if(unfoldingCard->GetNCentralityBins() != energyenergyCorrelatorCard[iFileType]->GetNCentralityBins()) {
        cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
        return;
      }
      for(int iCentrality = 0; iCentrality < unfoldingCard->GetNCentralityBins(); iCentrality++) {
        if(TMath::Abs(unfoldingCard->GetLowBinBorderCentrality(iCentrality) - energyenergyCorrelatorCard[iFileType]->GetLowBinBorderCentrality(iCentrality)) > 5) {
          cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
          return;
        }
        if(TMath::Abs(unfoldingCard->GetHighBinBorderCentrality(iCentrality) - energyenergyCorrelatorCard[iFileType]->GetHighBinBorderCentrality(iCentrality)) > 5) {
          cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
          return;
        }
      }
    }  // Checking compatible centrality bins

    if(unfoldingCard->GetNTrackPtBinsEEC() != energyenergyCorrelatorCard[iFileType]->GetNTrackPtBinsEEC()) {
      cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
    for(int iTrackPt = 0; iTrackPt < unfoldingCard->GetNTrackPtBinsEEC(); iTrackPt++) {
      if(TMath::Abs(unfoldingCard->GetLowBinBorderTrackPtEEC(iTrackPt) - energyenergyCorrelatorCard[iFileType]->GetLowBinBorderTrackPtEEC(iTrackPt)) > 0.01) {
        cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
        return;
      }
      if(TMath::Abs(unfoldingCard->GetHighBinBorderTrackPtEEC(iTrackPt) - energyenergyCorrelatorCard[iFileType]->GetHighBinBorderTrackPtEEC(iTrackPt)) > 0.01) {
        cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
        return;
      }
    }  // Checking the compatibility of track pT bins
  }
  
  if(unfoldingCard->GetNJetPtBinsEEC() != energyenergyCorrelatorCard[kDataFile]->GetNJetPtBinsEEC()) {
    cout << "Error! Measured jet pT bins do not match between the two files! Cannot execute the code!" << endl;
    return;
  }
  for(int iJetPt = 0; iJetPt < unfoldingCard->GetNJetPtBinsEEC(); iJetPt++) {
    if(TMath::Abs(unfoldingCard->GetLowBinBorderJetPtEEC(iJetPt) - energyenergyCorrelatorCard[kDataFile]->GetLowBinBorderJetPtEEC(iJetPt)) > 0.01) {
      cout << "Error! Measured jet pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
    if(TMath::Abs(unfoldingCard->GetHighBinBorderJetPtEEC(iJetPt) - energyenergyCorrelatorCard[kDataFile]->GetHighBinBorderJetPtEEC(iJetPt)) > 0.01) {
      cout << "Error! Measured jet pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
  }  // Checking the compatibility of reconstructed jet pT bins

  if(unfoldingCard->GetNJetPtBinsEEC() != energyenergyCorrelatorCard[kTruthReferenceFile]->GetNJetPtBinsEEC()) {
    cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
    return;
  }
  for(int iJetPt = 0; iJetPt < unfoldingCard->GetNJetPtBinsEEC(); iJetPt++) {
    if(TMath::Abs(unfoldingCard->GetLowBinBorderJetPtEEC(iJetPt) - energyenergyCorrelatorCard[kTruthReferenceFile]->GetLowBinBorderJetPtEEC(iJetPt)) > 0.01) {
      cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
    if(TMath::Abs(unfoldingCard->GetHighBinBorderJetPtEEC(iJetPt) - energyenergyCorrelatorCard[kTruthReferenceFile]->GetHighBinBorderJetPtEEC(iJetPt)) > 0.01) {
      cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
  }  // Checking the compatibility of generator level jet pT bins

  // ********************************************************
  //       Binning configuration for the unfolding study
  // ********************************************************
  
  // Find the number of bins from the card
  const int nCentralityBins = (isPbPbData) ? unfoldingCard->GetNCentralityBins() : 1;
  const int nTrackPtBinsEEC = unfoldingCard->GetNTrackPtBinsEEC();
  const int nJetPtBinsMeasured = energyenergyCorrelatorCard[kDataFile]->GetNJetPtBinsEEC();
  const int nJetPtBinsTruth = energyenergyCorrelatorCard[kTruthReferenceFile]->GetNJetPtBinsEEC();

  // Bin range to be studied
  int firstStudiedCentralityBin = 0;
  int lastStudiedCentralityBin = nCentralityBins-1;
  
  int firstStudiedTrackPtBinEEC = 0;
  int lastStudiedTrackPtBinEEC = 5;

  // Select explicitly the jet pT bins that we are going to unfold
  std::vector<std::pair<double,double>> unfoldedJetPtBins;
  unfoldedJetPtBins.push_back(std::make_pair(120,140));
  unfoldedJetPtBins.push_back(std::make_pair(140,160));
  unfoldedJetPtBins.push_back(std::make_pair(160,180));
  unfoldedJetPtBins.push_back(std::make_pair(180,200));
  unfoldedJetPtBins.push_back(std::make_pair(200,220));

  const int nUnfoldedJetPtBins = unfoldedJetPtBins.size();
  double unfoldedJetPtBinBorders[nUnfoldedJetPtBins+1];
  for(int iBorder = 0; iBorder < nUnfoldedJetPtBins; iBorder++){
    unfoldedJetPtBinBorders[iBorder] = unfoldedJetPtBins.at(iBorder).first;
  }
  unfoldedJetPtBinBorders[nUnfoldedJetPtBins] = unfoldedJetPtBins.at(nUnfoldedJetPtBins-1).second;

  // Default binning ranges for reference
  // centrality = {0,10,30,50,90}

  const int nIterations = 20; // Number of iterations used with the regularized unfolding

  const bool drawChi2map = false;
  const bool drawUnfoldedToTruthComparison = false;    // Compare unfolded distribution to truth reference

  const bool writeChi2ToFile = true; // Write the chi2 histograms to file
  TString outputFileName = Form("chi2Files/chi2Histograms_PbPb_energyWeightSquared_nominalReflectedCone_split%d_%s_4pCentShift_2023-12-04.root", iSplit, systematicName[iSystematic].Data());

  bool saveFigures = false;
  TString saveComment = "_bayesSwapped";
  TString figureFormat = "pdf";
    
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
  EECHistogramManager* energyEnergyCorrelatorHistograms[kNFileTypes];
  for(int iFileType = 0; iFileType < kNFileTypes - ignoreTruthReferenceFile; iFileType++){
    energyEnergyCorrelatorHistograms[iFileType] = new EECHistogramManager(energyEnergyCorrelatorInputFile[iFileType], energyenergyCorrelatorCard[iFileType]);
    energyEnergyCorrelatorHistograms[iFileType]->SetLoadEnergyEnergyCorrelators(true);
    energyEnergyCorrelatorHistograms[iFileType]->SetCentralityBinRange(firstStudiedCentralityBin,lastStudiedCentralityBin);
    energyEnergyCorrelatorHistograms[iFileType]->SetTrackPtBinRangeEEC(firstStudiedTrackPtBinEEC,lastStudiedTrackPtBinEEC);
    energyEnergyCorrelatorHistograms[iFileType]->SetJetPtBinRangeEEC(0, energyenergyCorrelatorCard[iFileType]->GetNJetPtBinsEEC());
    energyEnergyCorrelatorHistograms[iFileType]->LoadProcessedHistograms();
  }

  // Histograms that are needed to create the unfolding response
  TH1D* hUnfoldingMeasured[nCentralityBins][nTrackPtBinsEEC];
  TH1D* hUnfoldingTruth[nCentralityBins][nTrackPtBinsEEC];
  TH2D* hUnfoldingResponse[nCentralityBins][nTrackPtBinsEEC];
  TH1D* hUnfoldedDistribution[nCentralityBins][nTrackPtBinsEEC][nIterations];
  TH1D* hForwardFoldedDistribution[nCentralityBins][nTrackPtBinsEEC][nIterations];
  TH1D* energyEnergyCorrelatorForUnfolding[nCentralityBins][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorsFromData[nCentralityBins][nJetPtBinsMeasured][nTrackPtBinsEEC];
  TH1D* energyEnergyCorrelatorsTruthReference[nCentralityBins][nJetPtBinsTruth][nTrackPtBinsEEC];

  // Histograms for chi2 and error maps
  TH1D* hChi2map[nCentralityBins][nUnfoldedJetPtBins][nTrackPtBinsEEC];
  TH1D* hChi2mapForwardFold[nCentralityBins][nUnfoldedJetPtBins][nTrackPtBinsEEC];
  TH1D* hErrorMap[nCentralityBins][nUnfoldedJetPtBins][nTrackPtBinsEEC];

  // RooUnfold response object
  RooUnfoldResponse* rooResponse[nCentralityBins][nTrackPtBinsEEC];

  // Initialize all the histograms to null
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
      hUnfoldingMeasured[iCentrality][iTrackPt] = NULL;
      hUnfoldingTruth[iCentrality][iTrackPt] = NULL;
      hUnfoldingResponse[iCentrality][iTrackPt] = NULL;
      energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt] = NULL;
      rooResponse[iCentrality][iTrackPt] = NULL;
      for(int iIteration = 0; iIteration < nIterations; iIteration++){
        hUnfoldedDistribution[iCentrality][iTrackPt][iIteration] = NULL;
        hForwardFoldedDistribution[iCentrality][iTrackPt][iIteration] = NULL;
      }
      for(int iJetPt = 0; iJetPt < nJetPtBinsMeasured; iJetPt++){
        energyEnergyCorrelatorsFromData[iCentrality][iJetPt][iTrackPt] = NULL;
      }
      for(int iJetPt = 0; iJetPt < nJetPtBinsTruth; iJetPt++){
        energyEnergyCorrelatorsTruthReference[iCentrality][iJetPt][iTrackPt] = NULL;
      }
      for(int iJetPt = 0; iJetPt < nUnfoldedJetPtBins; iJetPt++){
        hChi2map[iCentrality][iJetPt][iTrackPt] = NULL;
        hChi2mapForwardFold[iCentrality][iJetPt][iTrackPt] = NULL;
        hErrorMap[iCentrality][iJetPt][iTrackPt] = NULL;
      }
    } // Track pT loop
  }  // Centrality loop

  // Read the histograms needed for the unfolding response and create the RooUnfold response objects
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      hUnfoldingMeasured[iCentrality][iTrackPt] = unfoldingHistograms->GetHistogramJetPtUnfoldingMeasured(iCentrality, iTrackPt);
      hUnfoldingTruth[iCentrality][iTrackPt] = unfoldingHistograms->GetHistogramJetPtUnfoldingTruth(iCentrality, iTrackPt);
      hUnfoldingResponse[iCentrality][iTrackPt] = unfoldingHistograms->GetHistogramJetPtUnfoldingResponse(iCentrality, iTrackPt);

      // Revome the overflow and underflow bins from the histograms
      removeOutOfRange(hUnfoldingMeasured[iCentrality][iTrackPt]);
      removeOutOfRange(hUnfoldingTruth[iCentrality][iTrackPt]);
      removeOutOfRange(hUnfoldingResponse[iCentrality][iTrackPt]);

      // Create the response object with trimmed histograms
      rooResponse[iCentrality][iTrackPt] = new RooUnfoldResponse(hUnfoldingMeasured[iCentrality][iTrackPt], hUnfoldingTruth[iCentrality][iTrackPt], hUnfoldingResponse[iCentrality][iTrackPt]);

      for(int iJetPt = 0; iJetPt < nJetPtBinsMeasured; iJetPt++){
        energyEnergyCorrelatorsFromData[iCentrality][iJetPt][iTrackPt] = energyEnergyCorrelatorHistograms[kDataFile]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt);
      } // Measured jet pT loop
      if(!ignoreTruthReferenceFile){
        for(int iJetPt = 0; iJetPt < nJetPtBinsTruth; iJetPt++){
          energyEnergyCorrelatorsTruthReference[iCentrality][iJetPt][iTrackPt] = energyEnergyCorrelatorHistograms[kTruthReferenceFile]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt);
        } // Truth jet pT loop
      } // If for reading the truth file 
    } // Track pT loop
  }  // Centrality loop

  // Create new histograms for chi2 maps
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      for(int iJetPt = 0; iJetPt < nUnfoldedJetPtBins; iJetPt++){
        hChi2map[iCentrality][iJetPt][iTrackPt] = new TH1D(Form("chi2map%d%d%d", iCentrality, iTrackPt, iJetPt), Form("chi2map%d%d%d", iCentrality, iTrackPt, iJetPt), nIterations, 0.5, 0.5+nIterations);
        hChi2mapForwardFold[iCentrality][iJetPt][iTrackPt] = new TH1D(Form("chi2mapForwardFold%d%d%d", iCentrality, iTrackPt, iJetPt), Form("hci2mapForwardFold%d%d%d", iCentrality, iTrackPt, iJetPt), nIterations, 0.5, 0.5+nIterations);
        hErrorMap[iCentrality][iJetPt][iTrackPt] = new TH1D(Form("errorMap%d%d%d", iCentrality, iTrackPt, iJetPt), Form("errorMap%d%d%d", iCentrality, iTrackPt, iJetPt), nIterations, 0.5, 0.5+nIterations);
      } // Jet pT loop
    } // Track pT loop
  }  // Centrality loop

  // Next, we need to transform the data histograms into a format that can be read by RooUnfold
  // For this, we will need to combine the jet pT and deltaR axes
  int nDeltaRBinsData = energyEnergyCorrelatorsFromData[firstStudiedCentralityBin][0][firstStudiedTrackPtBinEEC]->GetNbinsX();
  double jetPtLowerBound;
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt] = (TH1D*) hUnfoldingMeasured[iCentrality][iTrackPt]->Clone(Form("dataCorrelatorForUnfolding%d%d", iCentrality, iTrackPt));
      energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt]->Reset();
      for(int iJetPt = 0; iJetPt < nJetPtBinsMeasured; iJetPt++){
        jetPtLowerBound = energyenergyCorrelatorCard[kDataFile]->GetLowBinBorderJetPtEEC(iJetPt);
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
  RooUnfoldBayes* bayesUnfold[nCentralityBins][nTrackPtBinsEEC][nIterations];
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      for(int iIteration = 0; iIteration < nIterations; iIteration++){
        bayesUnfold[iCentrality][iTrackPt][iIteration] = new RooUnfoldBayes(rooResponse[iCentrality][iTrackPt], energyEnergyCorrelatorForUnfolding[iCentrality][iTrackPt], iIteration+1);
        hUnfoldedDistribution[iCentrality][iTrackPt][iIteration] = (TH1D*)bayesUnfold[iCentrality][iTrackPt][iIteration]->Hunfold();

        // Forward fold the unfolded distribution
        hForwardFoldedDistribution[iCentrality][iTrackPt][iIteration] = forwardFold(hUnfoldedDistribution[iCentrality][iTrackPt][iIteration], hUnfoldingResponse[iCentrality][iTrackPt]);

      } // Iteration loop
    } // Track pT loop
  }  // Centrality loop


  // *********************************************************************
  //     Transforming unfolded histograms back to one dimensional ones
  // *********************************************************************

  // Dissect the big histograms back into small histograms that can be compared
  TH1D* hMeasured[nCentralityBins][nUnfoldedJetPtBins][nTrackPtBinsEEC];
  TH1D* hTruth[nCentralityBins][nUnfoldedJetPtBins][nTrackPtBinsEEC];
  TH1D* hUnfolded[nCentralityBins][nUnfoldedJetPtBins][nTrackPtBinsEEC][nIterations];
  TH1D* hForwardFolded[nCentralityBins][nUnfoldedJetPtBins][nTrackPtBinsEEC][nIterations];
  TH1D* hMeasuredToTruthRatio[nCentralityBins][nUnfoldedJetPtBins][nTrackPtBinsEEC];
  TH1D* hUnfoldedToTruthRatio[nCentralityBins][nUnfoldedJetPtBins][nTrackPtBinsEEC][nIterations];
  TH1D* hForwardFoldedToMeasuredRatio[nCentralityBins][nUnfoldedJetPtBins][nTrackPtBinsEEC][nIterations];

  // Initialize the dissected histograms to null
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < nUnfoldedJetPtBins; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        hMeasured[iCentrality][iJetPt][iTrackPt] = NULL;
        hTruth[iCentrality][iJetPt][iTrackPt] = NULL;
        hMeasuredToTruthRatio[iCentrality][iJetPt][iTrackPt] = NULL;
        for(int iIteration = 0; iIteration < nIterations; iIteration++){
          hUnfolded[iCentrality][iJetPt][iTrackPt][iIteration] = NULL;
          hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][iIteration] = NULL;
          hForwardFolded[iCentrality][iJetPt][iTrackPt][iIteration] = NULL;
          hForwardFoldedToMeasuredRatio[iCentrality][iJetPt][iTrackPt][iIteration] = NULL;
        }
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
        hMeasured[iCentrality][iJetPt][iTrackPt] = new TH1D(Form("hMeasured%d%d%d", iCentrality, iTrackPt, iJetPt), Form("hMeasured%d%d%d", iCentrality, iTrackPt, iJetPt), nDeltaRBins, deltaRBinning);
        hTruth[iCentrality][iJetPt][iTrackPt] = new TH1D(Form("hTruth%d%d%d", iCentrality, iTrackPt, iJetPt), Form("hTruth%d%d%d", iCentrality, iTrackPt, iJetPt), nDeltaRBins, deltaRBinning);
        for(int iIteration = 0; iIteration < nIterations; iIteration++){
          hUnfolded[iCentrality][iJetPt][iTrackPt][iIteration] = new TH1D(Form("hUnfolded%d%d%d%d", iCentrality, iTrackPt, iJetPt, iIteration), Form("hUnfolded%d%d%d%d", iCentrality, iTrackPt, iJetPt, iIteration), nDeltaRBins, deltaRBinning);
          hForwardFolded[iCentrality][iJetPt][iTrackPt][iIteration] = new TH1D(Form("hForwardFolded%d%d%d%d", iCentrality, iTrackPt, iJetPt, iIteration), Form("hForwardFolded%d%d%d%d", iCentrality, iTrackPt, iJetPt, iIteration), nDeltaRBins, deltaRBinning);
        }
      } // Jet pT loop
    } // Track pT loop
  } // Centrality loop

  // Dissect the big histograms and fill small histograms based on the information
  // Also normalize the bins to bin width and histograms to one in the defined range
  double normalizationRegionLow = 0.006;
  double normalizationRegionHigh = 0.39;
  int jetPtUnfoldIndexMeasured = 0;
  int jetPtUnfoldIndexTruth = 0;
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      for(int iJetPt = 0; iJetPt < nUnfoldedJetPtBins; iJetPt++){
        jetPtUnfoldIndexMeasured = energyenergyCorrelatorCard[kDataFile]->FindBinIndexJetPtEEC(unfoldedJetPtBins.at(iJetPt));
        jetPtUnfoldIndexTruth = energyenergyCorrelatorCard[kTruthReferenceFile]->FindBinIndexJetPtEEC(unfoldedJetPtBins.at(iJetPt));
        for(int iBin = 1; iBin <= nDeltaRBins; iBin++){

          hMeasured[iCentrality][iJetPt][iTrackPt]->SetBinContent(iBin, hUnfoldingMeasured[iCentrality][iTrackPt]->GetBinContent(iBin + nDeltaRBins*jetPtUnfoldIndexMeasured) / hUnfoldingMeasured[iCentrality][iTrackPt]->GetBinWidth(iBin));
          hMeasured[iCentrality][iJetPt][iTrackPt]->SetBinError(iBin, hUnfoldingMeasured[iCentrality][iTrackPt]->GetBinError(iBin + nDeltaRBins*jetPtUnfoldIndexMeasured) / hUnfoldingMeasured[iCentrality][iTrackPt]->GetBinWidth(iBin));
          
          hTruth[iCentrality][iJetPt][iTrackPt]->SetBinContent(iBin, hUnfoldingTruth[iCentrality][iTrackPt]->GetBinContent(iBin + nDeltaRBins*jetPtUnfoldIndexTruth) / hUnfoldingTruth[iCentrality][iTrackPt]->GetBinWidth(iBin));
          hTruth[iCentrality][iJetPt][iTrackPt]->SetBinError(iBin, hUnfoldingTruth[iCentrality][iTrackPt]->GetBinError(iBin + nDeltaRBins*jetPtUnfoldIndexTruth) / hUnfoldingTruth[iCentrality][iTrackPt]->GetBinWidth(iBin));
        } // DeltaR bin loop

        for(int iIteration = 0; iIteration < nIterations; iIteration++){
          for(int iBin = 1; iBin <= nDeltaRBins; iBin++){
            hUnfolded[iCentrality][iJetPt][iTrackPt][iIteration]->SetBinContent(iBin, hUnfoldedDistribution[iCentrality][iTrackPt][iIteration]->GetBinContent(iBin + nDeltaRBins*jetPtUnfoldIndexMeasured) / hUnfoldedDistribution[iCentrality][iTrackPt][iIteration]->GetBinWidth(iBin));
            hUnfolded[iCentrality][iJetPt][iTrackPt][iIteration]->SetBinError(iBin, hUnfoldedDistribution[iCentrality][iTrackPt][iIteration]->GetBinError(iBin + nDeltaRBins*jetPtUnfoldIndexMeasured) / hUnfoldedDistribution[iCentrality][iTrackPt][iIteration]->GetBinWidth(iBin));

            hForwardFolded[iCentrality][iJetPt][iTrackPt][iIteration]->SetBinContent(iBin, hForwardFoldedDistribution[iCentrality][iTrackPt][iIteration]->GetBinContent(iBin + nDeltaRBins*jetPtUnfoldIndexMeasured) / hForwardFoldedDistribution[iCentrality][iTrackPt][iIteration]->GetBinWidth(iBin));
            hForwardFolded[iCentrality][iJetPt][iTrackPt][iIteration]->SetBinError(iBin, hForwardFoldedDistribution[iCentrality][iTrackPt][iIteration]->GetBinError(iBin + nDeltaRBins*jetPtUnfoldIndexMeasured) / hForwardFoldedDistribution[iCentrality][iTrackPt][iIteration]->GetBinWidth(iBin));
          } // DeltaR bin loop
        } // Iteration loop

        // Do the normalization of total distribution to one
        hMeasured[iCentrality][iJetPt][iTrackPt]->Scale(1.0 / hMeasured[iCentrality][iJetPt][iTrackPt]->Integral(hMeasured[iCentrality][iJetPt][iTrackPt]->FindBin(normalizationRegionLow), hMeasured[iCentrality][iJetPt][iTrackPt]->FindBin(normalizationRegionHigh), "width"));
        hTruth[iCentrality][iJetPt][iTrackPt]->Scale(1.0 / hTruth[iCentrality][iJetPt][iTrackPt]->Integral(hTruth[iCentrality][iJetPt][iTrackPt]->FindBin(normalizationRegionLow), hTruth[iCentrality][iJetPt][iTrackPt]->FindBin(normalizationRegionHigh), "width"));
        energyEnergyCorrelatorsFromData[iCentrality][jetPtUnfoldIndexMeasured][iTrackPt]->Scale(1.0 / energyEnergyCorrelatorsFromData[iCentrality][jetPtUnfoldIndexMeasured][iTrackPt]->Integral(energyEnergyCorrelatorsFromData[iCentrality][jetPtUnfoldIndexMeasured][iTrackPt]->FindBin(normalizationRegionLow), energyEnergyCorrelatorsFromData[iCentrality][jetPtUnfoldIndexMeasured][iTrackPt]->FindBin(normalizationRegionHigh), "width"));
        energyEnergyCorrelatorsTruthReference[iCentrality][jetPtUnfoldIndexTruth][iTrackPt]->Scale(1.0 / energyEnergyCorrelatorsTruthReference[iCentrality][jetPtUnfoldIndexTruth][iTrackPt]->Integral(energyEnergyCorrelatorsTruthReference[iCentrality][jetPtUnfoldIndexTruth][iTrackPt]->FindBin(normalizationRegionLow), energyEnergyCorrelatorsTruthReference[iCentrality][jetPtUnfoldIndexTruth][iTrackPt]->FindBin(normalizationRegionHigh), "width"));

        for(int iIteration = 0; iIteration < nIterations; iIteration++){
          hUnfolded[iCentrality][iJetPt][iTrackPt][iIteration]->Scale(1.0 / hUnfolded[iCentrality][iJetPt][iTrackPt][iIteration]->Integral(hUnfolded[iCentrality][iJetPt][iTrackPt][iIteration]->FindBin(normalizationRegionLow), hUnfolded[iCentrality][iJetPt][iTrackPt][iIteration]->FindBin(normalizationRegionHigh), "width"));

          hForwardFolded[iCentrality][iJetPt][iTrackPt][iIteration]->Scale(1.0 / hForwardFolded[iCentrality][iJetPt][iTrackPt][iIteration]->Integral(hForwardFolded[iCentrality][iJetPt][iTrackPt][iIteration]->FindBin(normalizationRegionLow), hForwardFolded[iCentrality][iJetPt][iTrackPt][iIteration]->FindBin(normalizationRegionHigh), "width"));
        } // Iteration loop
      } // Jet pT loop
    } // Track pT loop
  } // Centrality loop

  // Calculate the ratios of the other histograms to truth
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      for(int iJetPt = 0; iJetPt < nUnfoldedJetPtBins; iJetPt++){
        jetPtUnfoldIndexMeasured = energyenergyCorrelatorCard[kDataFile]->FindBinIndexJetPtEEC(unfoldedJetPtBins.at(iJetPt));
        jetPtUnfoldIndexTruth = energyenergyCorrelatorCard[kTruthReferenceFile]->FindBinIndexJetPtEEC(unfoldedJetPtBins.at(iJetPt));

        hMeasuredToTruthRatio[iCentrality][iJetPt][iTrackPt] = (TH1D*) energyEnergyCorrelatorsFromData[iCentrality][jetPtUnfoldIndexMeasured][iTrackPt]->Clone(Form("measuredToTruthRatio%d%d%d", iCentrality, iTrackPt, iJetPt));
        hMeasuredToTruthRatio[iCentrality][iJetPt][iTrackPt]->Divide(energyEnergyCorrelatorsTruthReference[iCentrality][jetPtUnfoldIndexTruth][iTrackPt]);

        for(int iIteration = 0; iIteration < nIterations; iIteration++){
          hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][iIteration] = (TH1D*) hUnfolded[iCentrality][iJetPt][iTrackPt][iIteration]->Clone(Form("unfoldedToTruthRatio%d%d%d%d", iCentrality, iTrackPt, iJetPt, iIteration));
          hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][iIteration]->Divide(energyEnergyCorrelatorsTruthReference[iCentrality][jetPtUnfoldIndexTruth][iTrackPt]);

          hForwardFoldedToMeasuredRatio[iCentrality][iJetPt][iTrackPt][iIteration] = (TH1D*) hForwardFolded[iCentrality][iJetPt][iTrackPt][iIteration]->Clone(Form("unfoldedToTruthRatio%d%d%d%d", iCentrality, iTrackPt, iJetPt, iIteration));
          hForwardFoldedToMeasuredRatio[iCentrality][iJetPt][iTrackPt][iIteration]->Divide(energyEnergyCorrelatorsFromData[iCentrality][jetPtUnfoldIndexMeasured][iTrackPt]);
        } // Iteration loop
      } // Jet pT loop
    } // Track pT loop
  } // Centrality loop

  // After all the ratios have been calculated, fit the unfolded to truth ratios with a constant and read the chi2 of the fit
  TF1* testFunction = new TF1("testFunction","1",normalizationRegionLow,normalizationRegionHigh);
  double chiSquare;
  double errorSquare;
  int firstBin, lastBin;
  for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
    for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
      for(int iJetPt = 0; iJetPt < nUnfoldedJetPtBins; iJetPt++){
        for(int iIteration = 0; iIteration < nIterations; iIteration++){
          chiSquare = hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][iIteration]->Chisquare(testFunction, "R");
          hChi2map[iCentrality][iJetPt][iTrackPt]->SetBinContent(iIteration+1,chiSquare);

          chiSquare = hForwardFoldedToMeasuredRatio[iCentrality][iJetPt][iTrackPt][iIteration]->Chisquare(testFunction, "R");
          hChi2mapForwardFold[iCentrality][iJetPt][iTrackPt]->SetBinContent(iIteration+1,chiSquare);

          // Calculate also the error map of the unfolded histogram
          firstBin = hUnfolded[iCentrality][iJetPt][iTrackPt][iIteration]->FindBin(normalizationRegionLow);
          lastBin = hUnfolded[iCentrality][iJetPt][iTrackPt][iIteration]->FindBin(normalizationRegionHigh);
          errorSquare = 0;
          for(int iBin = firstBin; iBin <= lastBin; iBin++){
            errorSquare += TMath::Power(hUnfolded[iCentrality][iJetPt][iTrackPt][iIteration]->GetBinError(iBin),2);
          }
          hErrorMap[iCentrality][iJetPt][iTrackPt]->SetBinContent(iIteration+1,errorSquare);
        } // Iteration loop
      } // Jet pT loop
    } // Track pT loop
  } // Centrality loop
  
  // Save the chi2, truth, measured, and undolfed histograms to a file
  if(writeChi2ToFile){
    TFile* outputFile = TFile::Open(outputFileName, "RECREATE");

    // Create directories for all histograms
    gDirectory->mkdir("chi2map");
    gDirectory->mkdir("chi2mapForwardFold");
    gDirectory->mkdir("errorMap");
    gDirectory->mkdir("truth");
    gDirectory->mkdir("measured");
    gDirectory->mkdir("unfolded");
    gDirectory->mkdir("measuredToTruthRatio");
    gDirectory->mkdir("unfoldedToTruthRatio");
    gDirectory->mkdir("forwardFoldedToMeasuredRatio");
    
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
        for(int iJetPt = 0; iJetPt < nUnfoldedJetPtBins; iJetPt++){

          jetPtUnfoldIndexMeasured = energyenergyCorrelatorCard[kDataFile]->FindBinIndexJetPtEEC(unfoldedJetPtBins.at(iJetPt));
          jetPtUnfoldIndexTruth = energyenergyCorrelatorCard[kTruthReferenceFile]->FindBinIndexJetPtEEC(unfoldedJetPtBins.at(iJetPt));

          gDirectory->cd("chi2map");
          hChi2map[iCentrality][iJetPt][iTrackPt]->Write();
          gDirectory->cd("../");

          gDirectory->cd("chi2mapForwardFold");
          hChi2mapForwardFold[iCentrality][iJetPt][iTrackPt]->Write();
          gDirectory->cd("../");

          gDirectory->cd("errorMap");
          hErrorMap[iCentrality][iJetPt][iTrackPt]->Write();
          gDirectory->cd("../");

          gDirectory->cd("truth");
          energyEnergyCorrelatorsTruthReference[iCentrality][jetPtUnfoldIndexTruth][iTrackPt]->Write(Form("hTruth%d%d%d", iCentrality, iTrackPt, iJetPt));
          gDirectory->cd("../");

          gDirectory->cd("measured");
          energyEnergyCorrelatorsFromData[iCentrality][jetPtUnfoldIndexMeasured][iTrackPt]->Write(Form("hMeasured%d%d%d", iCentrality, iTrackPt, iJetPt));
          gDirectory->cd("../");

          gDirectory->cd("unfolded");
          for(int iIteration = 0; iIteration < nIterations; iIteration++){
            hUnfolded[iCentrality][iJetPt][iTrackPt][iIteration]->Write();
          }
          gDirectory->cd("../");

          gDirectory->cd("measuredToTruthRatio");
          hMeasuredToTruthRatio[iCentrality][iJetPt][iTrackPt]->Write();
          gDirectory->cd("../");

          gDirectory->cd("unfoldedToTruthRatio");
          for(int iIteration = 0; iIteration < nIterations; iIteration++){
            hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][iIteration]->Write();
          }
          gDirectory->cd("../");

          gDirectory->cd("forwardFoldedToMeasuredRatio");
          for(int iIteration = 0; iIteration < nIterations; iIteration++){
            hForwardFoldedToMeasuredRatio[iCentrality][iJetPt][iTrackPt][iIteration]->Write();
          }
          gDirectory->cd("../");
        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
    
    // Write the card to the file
    unfoldingCard->AddVector(EECCard::kJetPtBinEdgesEEC,nUnfoldedJetPtBins+1,unfoldedJetPtBinBorders); // Save the unfolded jet pT bin borders to the file
    unfoldingCard->Write(outputFile);
    outputFile->Close();
  }

  // ******************************************************************
  //         Draw different comparisons to check the perfoemance
  // ******************************************************************

  // Prepare a JDrawer for drawing purposes
  JDrawer* drawer = new JDrawer();
  drawer->SetRelativeCanvasSize(1.1,1.2);
  drawer->SetLeftMargin(0.15);
  drawer->SetTopMargin(0.07);
  drawer->SetRightMargin(0.05);
  drawer->SetTitleOffsetY(1.2);

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

  int iterationColor[] = {kBlue, kGreen+3, kMagenta, kCyan, kOrange, kViolet+3, kPink-7, kSpring+3, kAzure-7};

  if(drawChi2map){

    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){

      // Set the centrality information for legends and figure saving
      if(isPbPbData) {
        centralityString = Form("Pythia+Hydjet: %.0f-%.0f", unfoldingCard->GetLowBinBorderCentrality(iCentrality), unfoldingCard->GetHighBinBorderCentrality(iCentrality));
        compactCentralityString = Form("_C%.0f-%.0f", unfoldingCard->GetLowBinBorderCentrality(iCentrality), unfoldingCard->GetHighBinBorderCentrality(iCentrality));
      } else {
        centralityString = "Pythia8";
        compactCentralityString = "_pythia8";
      }

      for(int iJetPt = 0; iJetPt < nUnfoldedJetPtBins; iJetPt++){

        // Set the jet pT information for legends and figure saving
        jetPtString = Form("%.0f < jet p_{T} < %.0f", unfoldedJetPtBins.at(iJetPt).first, unfoldedJetPtBins.at(iJetPt).second);
        compactJetPtString = Form("_J=%.0f-%.0f", unfoldedJetPtBins.at(iJetPt).first, unfoldedJetPtBins.at(iJetPt).second);

        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

          // Set the track pT information for legends and figure saving
          trackPtString = Form("%.1f < track p_{T}", unfoldingHistograms->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString = Form("_T%.0f", unfoldingHistograms->GetTrackPtBinBorderEEC(iTrackPt));

          // Draw the chi2 map to the canvas
          hChi2map[iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(kFullCircle);
          hChi2map[iCentrality][iJetPt][iTrackPt]->SetMarkerColor(kBlack);
          drawer->DrawHistogram(hChi2map[iCentrality][iJetPt][iTrackPt], "Number of iterations", "#chi^{2}", " ", "p");

          // Add a legend to the figure
          legend = new TLegend(0.15, 0.67, 0.4, 0.87);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);

          legend->AddEntry((TObject*)0, centralityString.Data(), "");
          legend->AddEntry((TObject*)0, jetPtString.Data(), "");
          legend->AddEntry((TObject*)0, trackPtString.Data(), "");

          legend->Draw();

          // Save the figures to a file
          if(saveFigures) {
            gPad->GetCanvas()->SaveAs(Form("figures/chi2map%s%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
          }

        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop

  } // Drawing chi2 map

  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  drawer->SetLogX(true);

  if(drawUnfoldedToTruthComparison){

    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){

      // Set the centrality information for legends and figure saving
      if(isPbPbData) {
        centralityString = Form("Pythia+Hydjet: %.0f-%.0f", unfoldingCard->GetLowBinBorderCentrality(iCentrality), unfoldingCard->GetHighBinBorderCentrality(iCentrality));
        compactCentralityString = Form("_C%.0f-%.0f", unfoldingCard->GetLowBinBorderCentrality(iCentrality), unfoldingCard->GetHighBinBorderCentrality(iCentrality));
      } else {
        centralityString = "Pythia8";
        compactCentralityString = "_pythia8";
      }

      for(int iJetPt = 0; iJetPt < nUnfoldedJetPtBins; iJetPt++){

        // Set the unfolding jet pT indices to match potentially different binnings
        jetPtUnfoldIndexMeasured = energyenergyCorrelatorCard[kDataFile]->FindBinIndexJetPtEEC(unfoldedJetPtBins.at(iJetPt));
        jetPtUnfoldIndexTruth = energyenergyCorrelatorCard[kTruthReferenceFile]->FindBinIndexJetPtEEC(unfoldedJetPtBins.at(iJetPt));

        // Set the jet pT information for legends and figure saving
        jetPtString = Form("%.0f < jet p_{T} < %.0f", unfoldedJetPtBins.at(iJetPt).first, unfoldedJetPtBins.at(iJetPt).second);
        compactJetPtString = Form("_J=%.0f-%.0f", unfoldedJetPtBins.at(iJetPt).first, unfoldedJetPtBins.at(iJetPt).second);

        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){

          // Create a new canvas for the plot
          drawer->CreateSplitCanvas();

          // Set the track pT information for legends and figure saving
          trackPtString = Form("%.1f < track p_{T}", unfoldingHistograms->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString = Form("_T%.0f", unfoldingHistograms->GetTrackPtBinBorderEEC(iTrackPt));

          // Draw first the generator level distribution
          energyEnergyCorrelatorsTruthReference[iCentrality][jetPtUnfoldIndexTruth][iTrackPt]->SetLineColor(kBlack);
          energyEnergyCorrelatorsTruthReference[iCentrality][jetPtUnfoldIndexTruth][iTrackPt]->GetXaxis()->SetRangeUser(normalizationRegionLow, normalizationRegionHigh);
          drawer->SetLogY(true);
          drawer->DrawHistogramToUpperPad(energyEnergyCorrelatorsTruthReference[iCentrality][jetPtUnfoldIndexTruth][iTrackPt], "#Deltar", "EEC", " ", "");

          // Add the reconstructed and unfolded distributions to the same plot
          energyEnergyCorrelatorsFromData[iCentrality][jetPtUnfoldIndexMeasured][iTrackPt]->SetLineColor(kRed);
          energyEnergyCorrelatorsFromData[iCentrality][jetPtUnfoldIndexMeasured][iTrackPt]->Draw("same");

          for(int iIteration = 0; iIteration < nIterations; iIteration++){
            hUnfolded[iCentrality][iJetPt][iTrackPt][iIteration]->SetLineColor(iterationColor[iIteration]);
            hUnfolded[iCentrality][iJetPt][iTrackPt][iIteration]->Draw("same");
          } 

          // Add a legend to the figure
          legend = new TLegend(0.25, 0.18- 0.03*nIterations, 0.5, 0.48 + 0.02*nIterations);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);

          legend->AddEntry((TObject*)0, centralityString.Data(), "");
          legend->AddEntry((TObject*)0, jetPtString.Data(), "");
          legend->AddEntry((TObject*)0, trackPtString.Data(), "");

          legend->AddEntry(energyEnergyCorrelatorsTruthReference[iCentrality][jetPtUnfoldIndexTruth][iTrackPt], "Generator level reference", "l");
          legend->AddEntry(energyEnergyCorrelatorsFromData[iCentrality][jetPtUnfoldIndexMeasured][iTrackPt], "Reconstructed jets + gen particles", "l");
          for(int iIteration = 0; iIteration < nIterations; iIteration++){
            legend->AddEntry(hUnfolded[iCentrality][iJetPt][iTrackPt][iIteration], Form("Unfolded result (%s)", unfoldName.ReplaceAll("NIT",Form("%d",iIteration+1)).Data()), "l");
            unfoldName.ReplaceAll(Form("%d",iIteration+1),"NIT");
          }

          legend->Draw();

          // Draw the ratios to lower pad
          drawer->SetLogY(false);
          hMeasuredToTruthRatio[iCentrality][iJetPt][iTrackPt]->SetLineColor(kRed);
          hMeasuredToTruthRatio[iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.8, 1.2);
          hMeasuredToTruthRatio[iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(normalizationRegionLow, normalizationRegionHigh);
          drawer->DrawHistogramToLowerPad(hMeasuredToTruthRatio[iCentrality][iJetPt][iTrackPt], "#Deltar", "Ratio to truth", " ", "");

          for(int iIteration = 0; iIteration < nIterations; iIteration++){
            hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][iIteration]->SetLineColor(iterationColor[iIteration]);
            hUnfoldedToTruthRatio[iCentrality][iJetPt][iTrackPt][iIteration]->Draw("same");
          }

          oneLine->Draw();

          // Save the figures to a file
          if(saveFigures) {
            gPad->GetCanvas()->SaveAs(Form("figures/jetPtUnfoldingWithRooUnfoldTruthCheck%s%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat.Data()));
          }
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Unfolded to truth comparison
}
