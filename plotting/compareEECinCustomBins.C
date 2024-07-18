#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"
#include "AlgorithmLibrary.h"

/*
 * Macro for comparing energy-energy correlators in custom bins
 */
void compareEECinCustomBins(){
  
  // Files for comparison
  const int nComparisonFiles = 2;
  TString fileName[nComparisonFiles];
  fileName[0] = "data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_nominalEnergyWeight_includeJetShape_processed_2024-07-17.root";
  fileName[1] = "data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_nominalEnergyWeight_includeJetShape_shiftedJetPt_processed_2024-07-17.root";
  
  
  // Open the files and check that they exist
  TFile* inputFile[nComparisonFiles];
  EECCard* card[nComparisonFiles];
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    
    inputFile[iFile] = TFile::Open(fileName[iFile]);
    
    if(inputFile[iFile] == NULL){
      cout << "Error! The file " << fileName[iFile].Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }

    card[iFile] = new EECCard(inputFile[iFile]);
  }
  
  // Check if we are using PbPb or pp data
  TString collisionSystem = card[0]->GetDataType();
  bool isPbPbData = collisionSystem.Contains("PbPb");
  
  // ====================================================
  //               Binning configuration
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = card[0]->GetNCentralityBins();
  const int nJetPtBinsEEC = card[0]->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card[0]->GetNTrackPtBinsEEC();

  // For legends, collect common and not common selections
  std::vector<TString> commonLegend;
  std::vector<TString> individualLegend;
  std::vector<TString> legendComment;

  // Define how much the jet pT is shifted in each file
  std::vector<double> jetPtShift;
  jetPtShift.push_back(0); // The first file in the list defines the unshifted bin
  jetPtShift.push_back(15); // The other files can contain results shifted with respect to the first file

  legendComment.push_back("(pp, no E-loss)");
  legendComment.push_back("(PbPb, before E-loss)");

  // Check that each file has a shift definition
  if(jetPtShift.size() < nComparisonFiles){
    cout << "ERROR! You have not defined a jet pT shift for every file!" << endl;
    cout << "Cannot run the code. Please define jet pT shifts for each compared file." << endl;
    return;
  }
  
  // Select explicitly which bins from the files are compared:
  std::vector<std::pair<double,double>> comparedCentralityBin;
  comparedCentralityBin.push_back(std::make_pair(0,100));

  std::vector<std::pair<double,double>> comparedJetPtBin;
  comparedJetPtBin.push_back(std::make_pair(120,140));
  comparedJetPtBin.push_back(std::make_pair(140,160));
  comparedJetPtBin.push_back(std::make_pair(160,180));
  comparedJetPtBin.push_back(std::make_pair(180,200));

  std::vector<double> comparedTrackPtBin;
  comparedTrackPtBin.push_back(1);

  // Create an object to easilty manipulate histograms
  AlgorithmLibrary* optimusPrimeTheTransformer = new AlgorithmLibrary();

  // ====================================================
  //                    Configuration
  // ====================================================
  
  // Figure saving
  const bool saveFigures = false;  // Save figures
  const char* saveComment = "_cumulantScaleExampleExtendedRegion";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures

  // Drawing configuration
  std::pair<double, double> ratioZoom = std::make_pair(0.65, 1.35);
  
  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager* histograms[nComparisonFiles];
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
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

  // Energy-energy correlator, jet shape, and cumulant histograms
  TH1D* hEnergyEnergyCorrelator[nComparisonFiles][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hEnergyEnergyCorrelatorRatio[nComparisonFiles][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hJetShape[nComparisonFiles][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hJetShapeCumulant[nComparisonFiles][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hCorrelatorCumulant[nComparisonFiles][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hJetShapeCumulantScaledRatio[nComparisonFiles][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* hCorrelatorCumulantScaledRatio[nComparisonFiles][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        for (int iFile = 0; iFile < nComparisonFiles; iFile++){
          hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt] = NULL;
          hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt] = NULL;
          hJetShape[iFile][iCentrality][iJetPt][iTrackPt] = NULL;
          hJetShapeCumulant[iFile][iCentrality][iJetPt][iTrackPt] = NULL;
          hCorrelatorCumulant[iFile][iCentrality][iJetPt][iTrackPt] = NULL;
          hJetShapeCumulantScaledRatio[iFile][iCentrality][iJetPt][iTrackPt] = NULL;
          hCorrelatorCumulantScaledRatio[iFile][iCentrality][iJetPt][iTrackPt] = NULL;
        } // File loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
  
  // Helper histograms
  std::pair<double, double> drawingRange = std::make_pair(0.008, 0.39);
  double epsilon = 0.0001;
  int lowNormalizationBin, highNormalizationBin;
  int iCentrality, iCentralityMatched;
  int iTrackPt, iTrackPtMatched;
  int iJetPt, iJetPtMatched;
  TH1D* helperHistogram;

  // Get the histograms from the histogram manager and normalize the signal histograms to one
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    for(auto centralityBin : comparedCentralityBin){
      iCentrality = isPbPbData ? card[0]->FindBinIndexCentrality(centralityBin) : 0;
      iCentralityMatched = isPbPbData ? card[iFile]->FindBinIndexCentrality(centralityBin) : 0;

      for(auto jetPtBin : comparedJetPtBin){
        iJetPt = card[0]->FindBinIndexJetPtEEC(jetPtBin);
        iJetPtMatched = card[iFile]->FindBinIndexJetPtEEC(jetPtBin.first + jetPtShift.at(iFile), jetPtBin.second + jetPtShift.at(iFile));

        for(auto trackPtBin : comparedTrackPtBin){
          iTrackPt = card[0]->GetBinIndexTrackPtEEC(trackPtBin);
          iTrackPtMatched = card[iFile]->GetBinIndexTrackPtEEC(trackPtBin);
      
          // Read the energy-energy correlator histograms
          hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt] = histograms[iFile]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentralityMatched, iJetPtMatched, iTrackPtMatched);

          // Normalize the distributions to one in the drawingRange
          lowNormalizationBin = hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(drawingRange.first + epsilon);
          highNormalizationBin = hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->GetXaxis()->FindBin(drawingRange.second - epsilon);

          hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->Scale(1 / hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->Integral(lowNormalizationBin, highNormalizationBin, "width"));

          // Calculate the ratio to first index
          hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt] = (TH1D*)hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->Clone(Form("eecRatio%d", iFile));
          hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->Divide(hEnergyEnergyCorrelator[0][iCentrality][iJetPt][iTrackPt]);

          // Read the jet shape histograms
          hJetShape[iFile][iCentrality][iJetPt][iTrackPt] = histograms[iFile]->GetHistogramJetShape(iCentralityMatched, iJetPtMatched, iTrackPtMatched);

          // Normalize the distributions to one in the drawingRange
          hJetShape[iFile][iCentrality][iJetPt][iTrackPt]->Scale(1 / hJetShape[iFile][iCentrality][iJetPt][iTrackPt]->Integral(lowNormalizationBin, highNormalizationBin, "width"));

          // Square the jet shape histograms
          optimusPrimeTheTransformer->SquareHistogram(hJetShape[iFile][iCentrality][iJetPt][iTrackPt]);

          // For the squared histogram, calculate the cumulant in the full DeltaR region and in the analysis region
          hJetShapeCumulant[iFile][iCentrality][iJetPt][iTrackPt] = optimusPrimeTheTransformer->GetCumulant(hJetShape[iFile][iCentrality][iJetPt][iTrackPt], lowNormalizationBin);
          hJetShapeCumulant[iFile][iCentrality][iJetPt][iTrackPt]->SetName(Form("jetShapeCumulant%d%d%d%d", iFile, iCentrality, iJetPt, iTrackPt));

          // Take the squared cumulant also for energy-energy correlators
          helperHistogram = (TH1D*) hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->Clone(Form("eecHelper%d%d%d%d", iFile, iCentrality, iJetPt, iTrackPt));
          optimusPrimeTheTransformer->SquareHistogram(helperHistogram);

          // Calculate the cumulant from the squared energy-energy correlator distribution
          hCorrelatorCumulant[iFile][iCentrality][iJetPt][iTrackPt] = optimusPrimeTheTransformer->GetCumulant(helperHistogram, lowNormalizationBin);
          hCorrelatorCumulant[iFile][iCentrality][iJetPt][iTrackPt]->SetName(Form("correlatorCumulant%d%d%d%d", iFile, iCentrality, iJetPt, iTrackPt));

        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop

  } // File loop

  // Helper histogram
  TH1D* scaleHistogram;

  // Scale the ratio histogram with the ratio of cumulants of squared jet shapes
  for(int iFile = 0; iFile < nComparisonFiles; iFile++){
    for(auto centralityBin : comparedCentralityBin){
      iCentrality = isPbPbData ? card[0]->FindBinIndexCentrality(centralityBin) : 0;
      for(auto jetPtBin : comparedJetPtBin){
        iJetPt = card[0]->FindBinIndexJetPtEEC(jetPtBin);
        for(auto trackPtBin : comparedTrackPtBin){
          iTrackPt = card[0]->GetBinIndexTrackPtEEC(trackPtBin);

          hJetShapeCumulantScaledRatio[iFile][iCentrality][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->Clone(Form("jetShapeCumulantScaledRatio%d%d%d%d", iFile, iCentrality, iJetPt, iTrackPt));
          scaleHistogram = (TH1D*) hJetShapeCumulant[0][iCentrality][iJetPt][iTrackPt]->Clone(Form("scaleForJetShapeCumulant%d%d%d%d", iFile, iCentrality, iJetPt, iTrackPt));
          scaleHistogram->Divide(hJetShapeCumulant[iFile][iCentrality][iJetPt][iTrackPt]);
          hJetShapeCumulantScaledRatio[iFile][iCentrality][iJetPt][iTrackPt]->Multiply(scaleHistogram);

          hCorrelatorCumulantScaledRatio[iFile][iCentrality][iJetPt][iTrackPt] = (TH1D*) hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->Clone(Form("correlatorCumulantScaledRatio%d%d%d%d", iFile, iCentrality, iJetPt, iTrackPt));
          scaleHistogram = (TH1D*) hCorrelatorCumulant[0][iCentrality][iJetPt][iTrackPt]->Clone(Form("scaleForCorrelatorCumulant%d%d%d%d", iFile, iCentrality, iJetPt, iTrackPt));
          scaleHistogram->Divide(hCorrelatorCumulant[iFile][iCentrality][iJetPt][iTrackPt]);
          hCorrelatorCumulantScaledRatio[iFile][iCentrality][iJetPt][iTrackPt]->Multiply(scaleHistogram);

        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // File loop
  
  
  // ==========================================================================
  //                    All the ratios in the same figure
  // ==========================================================================
  
  JDrawer* drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  
  drawer->SetLogX(true);

  TString compactCentralityString, compactTrackPtString, compactJetPtString;
  TString legendString;
  int markerStyle[2] = {kOpenSquare, kOpenCircle};
  int color[] = {kBlack, kRed, kBlue};
  TLegend* legend;

  for(auto centralityBin : comparedCentralityBin){
    iCentrality = isPbPbData ? card[0]->FindBinIndexCentrality(centralityBin) : 0;

    if(isPbPbData){
      compactCentralityString = Form("_C=%.0f-%.0f", centralityBin.first, centralityBin.second);
    } else {
      compactCentralityString = "_pp";  
    }

    for(auto jetPtBin : comparedJetPtBin){
      iJetPt = card[0]->FindBinIndexJetPtEEC(jetPtBin);
      compactJetPtString = Form("_J=%.0f-%.0f", jetPtBin.first, jetPtBin.second);
      for(auto trackPtBin : comparedTrackPtBin){
        iTrackPt = card[0]->GetBinIndexTrackPtEEC(trackPtBin);
        compactTrackPtString = Form("_T>%.1f",trackPtBin);

        // Create a new canvas for the plot
        drawer->CreateSplitCanvas();
          
        // Logarithmic EEC axis
        drawer->SetLogY(true);

        legend = new TLegend(0.18,0.04,0.45,0.52);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
        legend->AddEntry((TObject*) 0, Form("%s 5.02 TeV",card[0]->GetAlternativeDataType(false).Data()), "");

        // Add centrality and charged particle pT information
        if(isPbPbData){
          legend->AddEntry((TObject*) 0, Form("Cent: %.0f-%.0f%%", centralityBin.first, centralityBin.second), "");
        }
        legend->AddEntry((TObject*) 0, Form("p_{T}^{ch} > %.1f GeV", trackPtBin), "");

        // Set drawing style for all histograms
        for(int iFile = 0; iFile < nComparisonFiles; iFile++){
          hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iFile]);
          hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iFile]);
          hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iFile]);
          hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(markerStyle[iFile]);
          hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iFile]);
          hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iFile]);

          hJetShapeCumulantScaledRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(kFullDiamond);
          hJetShapeCumulantScaledRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(kGreen+3);
          hJetShapeCumulantScaledRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetLineColor(kGreen+3);

          hCorrelatorCumulantScaledRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(kFullCross);
          hCorrelatorCumulantScaledRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(kMagenta);
          hCorrelatorCumulantScaledRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetLineColor(kMagenta);
        }


        // Set the x-axis drawing range
        hEnergyEnergyCorrelator[0][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
          
        // Draw the histograms to the upper canvas
        drawer->DrawHistogramToUpperPad(hEnergyEnergyCorrelator[0][iCentrality][iJetPt][iTrackPt], "#Deltar", "EEC", " ");

        for(int iFile = 1; iFile < nComparisonFiles; iFile++){
          hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt]->Draw("same");
        }

        for(int iFile = 0; iFile < nComparisonFiles; iFile++){
          legend->AddEntry(hEnergyEnergyCorrelator[iFile][iCentrality][iJetPt][iTrackPt], Form("%.0f < jet p_{T} < %.0f GeV %s", jetPtBin.first + jetPtShift.at(iFile), jetPtBin.second + jetPtShift.at(iFile), legendComment.at(iFile).Data()), "p");
        }
        legend->AddEntry((TObject*) 0, "Cumulants calculated from analysis region", "");
  
        // Draw the legends to the upper pad
        legend->Draw();
          
        // Linear scale for the ratio
        drawer->SetLogY(false);
          
        // Set the drawing ranges
        hEnergyEnergyCorrelatorRatio[1][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);          
        hEnergyEnergyCorrelatorRatio[1][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);

        // Draw the histograms
        drawer->SetGridY(true);
        drawer->DrawHistogramToLowerPad(hEnergyEnergyCorrelatorRatio[1][iCentrality][iJetPt][iTrackPt], "#Deltar", Form("#frac{%.0f < jet p_{T} < %.0f }{%.0f < jet p_{T} < %.0f }", jetPtBin.first + jetPtShift.at(1), jetPtBin.second + jetPtShift.at(1), jetPtBin.first, jetPtBin.second), " ");

        for(int iFile = 2; iFile < nComparisonFiles; iFile++){
          hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->Draw("same");
        }
        hJetShapeCumulantScaledRatio[1][iCentrality][iJetPt][iTrackPt]->Draw("same");
        hCorrelatorCumulantScaledRatio[1][iCentrality][iJetPt][iTrackPt]->Draw("same");
        drawer->SetGridY(false);

        legend = new TLegend(0.25, 0.85, 0.95, 0.95);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.08);legend->SetTextFont(62);
        legend->SetNColumns(2);
        legend->AddEntry(hJetShapeCumulantScaledRatio[1][iCentrality][iJetPt][iTrackPt], "E1C cumulant scale", "pl");
        legend->AddEntry(hCorrelatorCumulantScaledRatio[1][iCentrality][iJetPt][iTrackPt], "E2C cumulant scale", "pl");
        legend->Draw();
          
        // Save the figures to a file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/eecCustomBinComparison%s%s%s%s.%s", saveComment, compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
        }
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

}
