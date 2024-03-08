#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"
#include "AlgorithmLibrary.h"
#include "SystematicUncertaintyOrganizer.h"

/*
 * Macro for comparing final energy-energy correlators results
 */
void comparePpToPythia(){

  enum enumDataType{kData, kPythia, knDataTypes};
  
  // Files for comparison
  const int weightExponent = 1;
  const int nWeightExponents = 2;

  TString fileName[knDataTypes][nWeightExponents];
  fileName[kPythia][0] = "data/ppData_pfJets_wtaAxis_nominalEnergyWeight_optimizedUnfoldingBins_fixedCovarianceMatrix_jet60or80triggers_unfoldingWithCovariance_processed_2024-01-23.root";
  fileName[kPythia][1] = "data/ppData_pfJets_wtaAxis_energyWeightSquared_optimizedUnfoldingBins_fixedCovarianceMatrix_jet60or80triggers_unfoldingWithCovariance_processed_2024-01-23.root";
  fileName[kData][0] = "data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_optimizedUnfoldingBins_nominalSmear_truthReference_processed_2024-01-11.root";
  fileName[kData][1] = "data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_optimizedUnfoldingBins_energyWeightSquared_nominalSmear_truthReference_processed_2024-01-10.root";

  TString uncertaintyFileName[nWeightExponents];
  uncertaintyFileName[0] = "systematicUncertainties/systematicUncertainties_pp_nominalEnergyWeight_includeMCnonClosure_2024-03-07.root";
  uncertaintyFileName[1] = "systematicUncertainties/systematicUncertainties_pp_energyWeightSquared_includeMCnonClosure_2024-03-07.root";

  TString fileDescription[knDataTypes];
  fileDescription[1] = "pp data";
  fileDescription[0] = "Pythia8 truth";
  
  // Open the files and check that they exist
  TFile* inputFile[knDataTypes];
  EECCard* card[knDataTypes];
  for(int iFile = 0; iFile < knDataTypes; iFile++){
    
    inputFile[iFile] = TFile::Open(fileName[iFile][weightExponent-1]);
    
    if(inputFile[iFile] == NULL){
      cout << "Error! The file " << fileName[iFile][weightExponent-1].Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }

    card[iFile] = new EECCard(inputFile[iFile]);
  }

  TFile* uncertaintyFile = TFile::Open(uncertaintyFileName[weightExponent-1]);

  if(uncertaintyFile == NULL){
    cout << "Error! The file " << uncertaintyFileName[weightExponent-1].Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }

  EECCard* uncertaintyCard = new EECCard(uncertaintyFile);
  
  // ====================================================
  //               Binning configuration
  // ====================================================
  
  // Find the number of bins from the card
  const int nJetPtBinsEEC = card[0]->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card[0]->GetNTrackPtBinsEEC();
  
  std::vector<std::pair<double,double>> comparedJetPtBin;
  comparedJetPtBin.push_back(std::make_pair(120,140));
  comparedJetPtBin.push_back(std::make_pair(140,160));
  comparedJetPtBin.push_back(std::make_pair(160,180));
  comparedJetPtBin.push_back(std::make_pair(180,200));

  std::vector<double> comparedTrackPtBin;
  if(card[0]->GetWeightExponent() == 2) comparedTrackPtBin.push_back(1.0);
  comparedTrackPtBin.push_back(2.0);
  comparedTrackPtBin.push_back(3.0);

  // Choose the type of draw energy-energy correlator
  // EECHistogramManager::kEnergyEnergyCorrelatorNormalized = Normalized energy-energy correlator
  // EECHistogramManager::kEnergyEnergyCorrelatorBackground = Estimated background
  // EECHistogramManager::kEnergyEnergyCorrelatorSignal = Background subtracted, but not unfolded energy-energy correlator
  // EECHistogramManager::kEnergyEnergyCorrelatorUnfolded = Unfolded energy-energy correlator
  // EECHistogramManager::kEnergyEnergyCorrelatorBackgroundAfterUnfolding = Estimated background after unfolding
  // EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal = Unfolded energy-energy correlator signal
  // EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels = Raw energy-energy correlator
  int drawnEnergyEnergyCorrelator = EECHistogramManager::kEnergyEnergyCorrelatorUnfoldedSignal;

  // ====================================================
  //                Drawing configuration
  // ====================================================
  
  // Figure saving
  const bool saveFigures = false;  // Save figures
  TString saveComment = "";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures

  if(weightExponent == 2){
    saveComment.Prepend("_energyWeightSquared");
  } else {
    saveComment.Prepend("_nominalEnergyWeight");
  }

  // Drawing configuration
  std::pair<double, double> ratioZoom = std::make_pair(0.7, 1.3);
  std::pair<double, double> eecZoom = std::make_pair(0.2, 30);
  const bool automaticZoom = true;

  // Sanity checks for input. Ensure that all the selected bins actually exist in the input files.
  // This check is only needed for unfolded bins, so skip it if only raw distribution is drawn.
  if(drawnEnergyEnergyCorrelator != EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels){
    for(int iFile = 0; iFile < knDataTypes; iFile++){

      // Sanity check for jet pT bins
      for(auto jetPtBin : comparedJetPtBin){
        if(card[iFile]->FindBinIndexJetPtEEC(jetPtBin) < card[iFile]->GetFirstUnfoldedJetPtBin() || card[iFile]->FindBinIndexJetPtEEC(jetPtBin) > card[iFile]->GetLastUnfoldedJetPtBin()){
          cout << "ERROR! Jet pT bin " << jetPtBin.first << "-" << jetPtBin.second << " does not exist in file " << fileName[iFile][weightExponent-1].Data() << endl;
          cout << "Please only choose jet pT bins that are included in the input files." << endl;
          return;
        }
      }

      // Sanity check for track pT bins
      for(auto trackPtBin : comparedTrackPtBin){
        if(card[iFile]->GetBinIndexTrackPtEEC(trackPtBin) < card[iFile]->GetFirstUnfoldedTrackPtBin() || card[iFile]->GetBinIndexTrackPtEEC(trackPtBin) > card[iFile]->GetLastUnfoldedTrackPtBin()){
          cout << "ERROR! Track pT cut > " << trackPtBin << " GeV does not exist in file " << fileName[iFile][weightExponent-1].Data() << endl;
          cout << "Please only choose track pT bins that are included in the input files." << endl;
          return;
        }
      } 
    } // File loop for input sanity check
  } // Unfolded distributions
  
  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager* histograms[knDataTypes];
  for(int iFile = 0; iFile < knDataTypes; iFile++){
    histograms[iFile] = new EECHistogramManager(inputFile[iFile], card[iFile]);

    // Choose the energy-energy correlator types to load
    histograms[iFile]->SetLoadEnergyEnergyCorrelators(true);

    // Choose the bin ranges
    histograms[iFile]->SetCentralityBinRange(0, card[iFile]->GetNCentralityBins() - 1);
    histograms[iFile]->SetJetPtBinRangeEEC(0, card[iFile]->GetNJetPtBinsEEC() - 1);
    histograms[iFile]->SetTrackPtBinRangeEEC(0, card[iFile]->GetNTrackPtBinsEEC() - 1);

    // Load the histograms from the file
    histograms[iFile]->LoadProcessedHistograms();
  }

  // Create systematic uncertainty organizer to illustrate to draw the pp systematic uncertainties to plots
  SystematicUncertaintyOrganizer* uncertaintyOrganizer = new SystematicUncertaintyOrganizer(uncertaintyFile);

  // Energy-energy correlator histograms
  TH1D* hEnergyEnergyCorrelator[knDataTypes][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hEnergyEnergyCorrelatorRatio[knDataTypes][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hRelativeUncertaintyPp[nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hAbsoluteUncertaintyPp[nJetPtBinsEEC][nTrackPtBinsEEC];
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
      hRelativeUncertaintyPp[iJetPt][iTrackPt] = NULL;
      hAbsoluteUncertaintyPp[iJetPt][iTrackPt] = NULL;
      for(int iFile = 0; iFile < knDataTypes; iFile++){
        hEnergyEnergyCorrelator[iFile][iJetPt][iTrackPt] = NULL;
        hEnergyEnergyCorrelatorRatio[iFile][iJetPt][iTrackPt] = NULL;
      } // File loop
    } // Track pT loop
  } // Jet pT loop
  
  // Transformer to transform absolute uncertainties to relative ones
  AlgorithmLibrary* optimusPrimeTheTransformer = new AlgorithmLibrary();
  
  // Helper variables
  std::pair<double, double> drawingRange = std::make_pair(0.008, 0.39);
  double epsilon = 0.0001;
  int lowNormalizationBin, highNormalizationBin;
  int iTrackPt, iTrackPtReference, iTrackPtUncertainty;
  int iJetPt, iJetPtReference, iJetPtUncertainty;

  // Get the histograms from the histogram manager and normalize the signal histograms to one
  for(int iFile = 0; iFile < knDataTypes; iFile++){
    for(auto jetPtBin : comparedJetPtBin){
      for(auto trackPtBin : comparedTrackPtBin){

        // Find the proper binning and express it in term of the bins in the first file
        iJetPt = card[iFile]->FindBinIndexJetPtEEC(jetPtBin);
        iTrackPt = card[iFile]->GetBinIndexTrackPtEEC(trackPtBin);
        iJetPtReference = card[0]->FindBinIndexJetPtEEC(jetPtBin);
        iTrackPtReference = card[0]->GetBinIndexTrackPtEEC(trackPtBin);

        // Load the selected energy-energy correlator histogram
        if(drawnEnergyEnergyCorrelator == EECHistogramManager::knEnergyEnergyCorrelatorProcessingLevels){
          hEnergyEnergyCorrelator[iFile][iJetPtReference][iTrackPtReference] = histograms[iFile]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPt, iTrackPt);
        } else {
          hEnergyEnergyCorrelator[iFile][iJetPtReference][iTrackPtReference] = histograms[iFile]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, 0, iJetPt, iTrackPt, drawnEnergyEnergyCorrelator);
        }

        // Normalize the distributions to one in the drawingRange
        lowNormalizationBin = hEnergyEnergyCorrelator[iFile][iJetPtReference][iTrackPtReference]->GetXaxis()->FindBin(drawingRange.first + epsilon);
        highNormalizationBin = hEnergyEnergyCorrelator[iFile][iJetPtReference][iTrackPtReference]->GetXaxis()->FindBin(drawingRange.second - epsilon);

        hEnergyEnergyCorrelator[iFile][iJetPtReference][iTrackPtReference]->Scale(1 / hEnergyEnergyCorrelator[iFile][iJetPtReference][iTrackPtReference]->Integral(lowNormalizationBin, highNormalizationBin, "width"));

        // Read also the uncertainty histograms
        if(iFile == 0){
          iJetPtUncertainty = uncertaintyCard->FindBinIndexJetPtEEC(jetPtBin);
          iTrackPtUncertainty = uncertaintyCard->GetBinIndexTrackPtEEC(trackPtBin);

          hRelativeUncertaintyPp[iJetPt][iTrackPt] = (TH1D*) uncertaintyOrganizer->GetSystematicUncertainty(0, iJetPtUncertainty, iTrackPtUncertainty)->Clone(Form("relativeUncertainty%d%d", iJetPt, iTrackPt));
          optimusPrimeTheTransformer->TransformToRelativeUncertainty(hRelativeUncertaintyPp[iJetPt][iTrackPt], true);
          hAbsoluteUncertaintyPp[iJetPt][iTrackPt] = (TH1D*) uncertaintyOrganizer->GetSystematicUncertainty(0, iJetPtUncertainty, iTrackPtUncertainty)->Clone(Form("absoluteUncertainty%d%d", iJetPt, iTrackPt));
        }

      } // Track pT loop
    } // Jet pT loop
  } // File loop

  // After all the histograms have been read, calculate the ratios
  for(int iFile = 0; iFile < knDataTypes; iFile++){
    for(auto jetPtBin : comparedJetPtBin){
      iJetPt = card[0]->FindBinIndexJetPtEEC(jetPtBin);
      for(auto trackPtBin : comparedTrackPtBin){
        iTrackPt = card[0]->GetBinIndexTrackPtEEC(trackPtBin);

        hEnergyEnergyCorrelatorRatio[iFile][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelator[iFile][iJetPt][iTrackPt]->Clone(Form("eecRatio%d%d%d", iFile, iJetPt, iTrackPt));
        hEnergyEnergyCorrelatorRatio[iFile][iJetPt][iTrackPt]->Divide(hEnergyEnergyCorrelator[0][iJetPt][iTrackPt]);

      } // Track pT loop
    } // Jet pT loop
  } // File loop
  
  // ==========================================================================
  //                Draw all the distribution in the same figure
  // ==========================================================================
  
  JDrawer* drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  
  drawer->SetLogX(true);

  TString compactJetPtString = "";
  TString compactTrackPtString = "";
  TString energyWeightString = (card[0]->GetWeightExponent() == 1) ? "Nominal energy weight" : "Energy weight squared";
  TString legendString;
  int markerStyle[5] = {kOpenSquare, kFullCircle, kOpenCross, kFullStar, kFullCross};
  int color[] = {kRed, kBlack,kBlue,kGreen+3,kMagenta,kCyan,kOrange,kViolet+3,kPink-7,kSpring+3,kAzure-7};

  // Binning vectors
  std::vector<int> currentJetPtIndices;
  std::vector<int> currentTrackPtIndices;
  int legendJetPtIndex = 0;
  int legendTrackPtIndex = 0;
  double minimumCandidate, maximumCandidate;

  for(auto jetPtBin : comparedJetPtBin){
    iJetPt = card[0]->FindBinIndexJetPtEEC(jetPtBin);
    for(auto trackPtBin : comparedTrackPtBin){
      iTrackPt = card[0]->GetBinIndexTrackPtEEC(trackPtBin);
          
      // Create a new canvas for the plot
      drawer->CreateSplitCanvas();

      // Logarithmic EEC axis
      drawer->SetLogY(true);

      // Create the legend and add jet and track pT information to it
      TLegend* legend = new TLegend(0.18,0.04,0.45,0.58);
      legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);

      legend->AddEntry((TObject*) 0, Form("%.0f < jet p_{T} < %.0f", jetPtBin.first, jetPtBin.second), "");
      legend->AddEntry((TObject*) 0, Form("%.1f < track p_{T}", trackPtBin), "");
      legend->AddEntry((TObject*) 0, energyWeightString, "");

      // Setup jet and track pT strings
      compactJetPtString = Form("_J=%.0f-%.0f", jetPtBin.first, jetPtBin.second);
      compactTrackPtString = Form("_T>%.1f",trackPtBin);
      compactTrackPtString.ReplaceAll(".","v");


      // Set drawing style for all histograms
      for(int iFile = 0; iFile < knDataTypes; iFile++){
        hEnergyEnergyCorrelator[iFile][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iFile]);
        hEnergyEnergyCorrelator[iFile][iJetPt][iTrackPt]->SetMarkerColor(color[iFile]);
        hEnergyEnergyCorrelator[iFile][iJetPt][iTrackPt]->SetLineColor(color[iFile]); 
        hEnergyEnergyCorrelatorRatio[iFile][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iFile]);
        hEnergyEnergyCorrelatorRatio[iFile][iJetPt][iTrackPt]->SetMarkerColor(color[iFile]);
        hEnergyEnergyCorrelatorRatio[iFile][iJetPt][iTrackPt]->SetLineColor(color[iFile]);
      } // File loop

      // Set the drawing style for systematic uncertainties
      hAbsoluteUncertaintyPp[iJetPt][iTrackPt]->SetFillColorAlpha(kBlack, 0.2);
      hAbsoluteUncertaintyPp[iJetPt][iTrackPt]->SetMarkerStyle(9);
      hAbsoluteUncertaintyPp[iJetPt][iTrackPt]->SetMarkerSize(0);
      hRelativeUncertaintyPp[iJetPt][iTrackPt]->SetFillColorAlpha(kBlack, 0.2);
      hRelativeUncertaintyPp[iJetPt][iTrackPt]->SetMarkerStyle(9);
      hRelativeUncertaintyPp[iJetPt][iTrackPt]->SetMarkerSize(0);

      // Automatic zooming for the drawn histograms
      if(automaticZoom){
        eecZoom.first = 10000;
        eecZoom.second = 0;
        for(int iFile = 0; iFile < knDataTypes; iFile++){
          hEnergyEnergyCorrelator[iFile][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
          minimumCandidate = hEnergyEnergyCorrelator[iFile][iJetPt][iTrackPt]->GetMinimum();
          maximumCandidate = hEnergyEnergyCorrelator[iFile][iJetPt][iTrackPt]->GetMaximum();
          if(minimumCandidate < eecZoom.first) eecZoom.first = minimumCandidate;
          if(maximumCandidate > eecZoom.second) eecZoom.second = maximumCandidate;
        } // File loop
        eecZoom.first = eecZoom.first / 2.0;
        eecZoom.second = eecZoom.second * 2.0;
      }

      // Set the x- and y-axis drawing ranges
      hAbsoluteUncertaintyPp[iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
      hAbsoluteUncertaintyPp[iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(eecZoom.first, eecZoom.second);
          
      // Draw the histograms to the upper canvas
      drawer->DrawHistogramToUpperPad(hAbsoluteUncertaintyPp[iJetPt][iTrackPt], "#Deltar", Form("EEC %s", histograms[0]->GetEnergyEnergyCorrelatorProcessSaveName(drawnEnergyEnergyCorrelator)), " ", "e2");

      for(int iFile = 0; iFile < knDataTypes; iFile++){
        hEnergyEnergyCorrelator[iFile][iJetPt][iTrackPt]->Draw("same");
      } // File loop

      // Add legends for drawn histograms
      legendJetPtIndex = 0;
      legendTrackPtIndex = 0;
      for(int iFile = 0; iFile < knDataTypes; iFile++){
        legend->AddEntry(hEnergyEnergyCorrelator[iFile][iJetPt][iTrackPt], fileDescription[iFile].Data(), "p");
      } // File loop
  
      // Draw the legends to the upper pad
      legend->Draw();
          
      // Linear scale for the ratio
      drawer->SetLogY(false);
          
      // Set the drawing ranges
      hRelativeUncertaintyPp[iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);          
      hRelativeUncertaintyPp[iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

      // Draw the histograms
      drawer->SetGridY(true);

      drawer->DrawHistogramToLowerPad(hRelativeUncertaintyPp[iJetPt][iTrackPt], "#Deltar", "#frac{pp}{Pythia8}", " ", "e2");
      for(int iFile = 0; iFile < knDataTypes; iFile++){
        hEnergyEnergyCorrelatorRatio[iFile][iJetPt][iTrackPt]->Draw("same");
      } // File loop
      drawer->SetGridY(false);
          
      // Save the figures to a file
      if(saveFigures){
      gPad->GetCanvas()->SaveAs(Form("figures/eecPpToPythiaComparison%s%s%s.%s", saveComment.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
      }
    } // Track pT loop
  } // Jet pT loop

}
