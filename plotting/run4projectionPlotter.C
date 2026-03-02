 #include "AlgorithmLibrary.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
 #include "JDrawer.h"

/*
 * Macro for plotting run4 projected uncertainties
 */
void run4projectionPlotter(){

  // Luminosity scaling numbers in nb-1
  double currentLuminosity = 1.7;
  double run4luminosity = 7;
  double luminosityScale = (TMath::Sqrt(run4luminosity) / run4luminosity) / (TMath::Sqrt(currentLuminosity) / currentLuminosity);

  // Enumeration for current results and projections
  enum enumDataType{kCurrentResult, kProjection, knDataTypes};

  // Open the files for current results and projections
  TString fileName[knDataTypes] = {"run4projection/hepdata_energyEnergyCorrelatorRatio_hin-23-004.root",
                                   "run4projection/hepdata_energyEnergyCorrelatorRatio_run4projectionsForJES_hin-23-004.root"};

  TFile* inputFile[knDataTypes];
  for(int iFile = 0; iFile < knDataTypes; iFile++){
    
    inputFile[iFile] = TFile::Open(fileName[iFile]);
    
    if(inputFile[iFile] == NULL){
      cout << "Error! The file " << fileName[iFile].Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the run4projection/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }
  }
  
  // ====================================================
  //               Binning configuration
  // ====================================================

  // Choose which bins from the files are compared
  std::vector<std::pair<int,int>> comparedCentralityBin;
  comparedCentralityBin.push_back(std::make_pair(0,10));
  const int nCentralityBins = comparedCentralityBin.size();
  
  std::vector<std::pair<int,int>> comparedJetPtBin;
  comparedJetPtBin.push_back(std::make_pair(120,140));
  const int nJetPtBins = comparedJetPtBin.size();

  std::vector<int> comparedTrackPtBin;
  comparedTrackPtBin.push_back(1);
  const int nTrackPtBins = comparedTrackPtBin.size();

  // ====================================================
  //                Drawing configuration
  // ====================================================

  // Axis draw ranges
  std::pair<double, double> drawingRange = std::make_pair(0.008, 0.39);
  std::pair<double, double> ratioZoom = std::make_pair(0.15, 1.85);
  
  // Figure saving
  const bool saveFigures = false;  // Save figures
  TString saveComment = "";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures

  // ====================================================
  //       Read selected histograms from the files
  // ====================================================

  // Energy-energy correlator histograms
  TH1D* hEnergyEnergyCorrelatorRatio[knDataTypes][nCentralityBins][nJetPtBins][nTrackPtBins];
  TH1D* hEnergyEnergyCorrelatorRatioSystematics[knDataTypes][nCentralityBins][nJetPtBins][nTrackPtBins];
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iFile = 0; iFile < knDataTypes; iFile++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){
        for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){     
          hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt] = NULL;
          hEnergyEnergyCorrelatorRatioSystematics[iFile][iCentrality][iJetPt][iTrackPt] = NULL;
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // File loop
  
  // Create an uncertainty dealer to deal with uncertainties
  AlgorithmLibrary* uncertaintyDealer = new AlgorithmLibrary();

  // Get the histograms from the histogram manager and normalize the signal histograms to one
  TH1D* correlatedUncertainty;
  TH1D* uncorrelatedUncertainty;
  int iCentrality = -1;
  int iJetPt = -1;
  int iTrackPt = -1;
  for(int iFile = 0; iFile < knDataTypes; iFile++){
    iCentrality = -1;
    for(auto centralityBin : comparedCentralityBin){
      iCentrality++;
      iJetPt = -1;
      for(auto jetPtBin : comparedJetPtBin){
        iJetPt++;
        iTrackPt = -1;
        for(auto trackPtBin : comparedTrackPtBin){
          iTrackPt++;

          cout << "Doing bin " << iFile << iCentrality << iJetPt << iTrackPt << endl;

          // Read the histogram with statistical uncertainties from the file
          hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt] = (TH1D*) inputFile[iFile]->Get(Form("energyEnergyCorrelatorRatio_nominalEnergyWeight_C%d-%d_J%d-%d_T%d", centralityBin.first, centralityBin.second, jetPtBin.first, jetPtBin.second, trackPtBin));

          // If we are doing the run4 projection, we need to scale the statistical uncertainty to match predicted luminosity
          if(iFile == kProjection){
            uncertaintyDealer->ScaleUncertainties(hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt], luminosityScale);
          }

          // The systematic uncertainties are in two different histograms, one for correlated and one for uncorrelated systematic uncertainties
          // Read them both, and then combine them into one systematic uncertainty band
          correlatedUncertainty = (TH1D*) inputFile[iFile]->Get(Form("energyEnergyCorrelatorRatio_correlatedUncertainty_nominalEnergyWeight_C%d-%d_J%d-%d_T%d", centralityBin.first, centralityBin.second, jetPtBin.first, jetPtBin.second, trackPtBin));
          uncorrelatedUncertainty = (TH1D*) inputFile[iFile]->Get(Form("energyEnergyCorrelatorRatio_uncorrelatedUncertainty_nominalEnergyWeight_C%d-%d_J%d-%d_T%d", centralityBin.first, centralityBin.second, jetPtBin.first, jetPtBin.second, trackPtBin));

          hEnergyEnergyCorrelatorRatioSystematics[iFile][iCentrality][iJetPt][iTrackPt] = uncertaintyDealer->CombineUncertainties(correlatedUncertainty, uncorrelatedUncertainty);

        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // File loop

  cout << "All done" << endl;
  
  // ==========================================================================
  //       Draw the projected uncertainties together with the current ones
  // ==========================================================================
  
  JDrawer* drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  
  drawer->SetLogX(true);

  TString compactJetPtString;
  TString compactTrackPtString;
  TString compactCentralityString;
  TString jetPtString;
  TString trackPtString;
  TString centralityString;
  TString legendString;
  int color[] = {kRed, kBlue};

  iCentrality = -1;
  iJetPt = -1;
  iTrackPt = -1;
  for(auto centralityBin : comparedCentralityBin){
    iCentrality++;
    iJetPt = -1;
    centralityString = Form("C=%d-%d%%", centralityBin.first, centralityBin.second);
    compactCentralityString = Form("_C%d-%d", centralityBin.first, centralityBin.second);
    for(auto jetPtBin : comparedJetPtBin){
      iJetPt++;
      iTrackPt = -1;
      jetPtString = Form("%d < p_{T,jet} < %d GeV", jetPtBin.first, jetPtBin.second);
      compactJetPtString = Form("_J%d-%d", jetPtBin.first, jetPtBin.second);
      for(auto trackPtBin : comparedTrackPtBin){
        iTrackPt++;
        trackPtString = Form("p_{T}^{ch} > %d GeV", trackPtBin);
        compactTrackPtString = Form("_T%d", trackPtBin);
          
        // Create the legend and add binning information to it
        TLegend* legend = new TLegend(0.18,0.04,0.45,0.58);
        legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);

        legend->AddEntry((TObject*) 0, centralityString, "");
        legend->AddEntry((TObject*) 0, jetPtString, "");
        legend->AddEntry((TObject*) 0, trackPtString, "");


        // Set drawing style for all histograms
        for(int iFile = 0; iFile < knDataTypes; iFile++){
          hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerColor(color[iFile]);
          hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iFile]);
          hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->SetLineWidth(3);

          hEnergyEnergyCorrelatorRatioSystematics[iFile][iCentrality][iJetPt][iTrackPt]->SetFillColorAlpha(color[iFile], 0.2);
          hEnergyEnergyCorrelatorRatioSystematics[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(9);
          hEnergyEnergyCorrelatorRatioSystematics[iFile][iCentrality][iJetPt][iTrackPt]->SetMarkerSize(0);
        } // File loop


        // Set the x- and y-axis drawing ranges
        hEnergyEnergyCorrelatorRatioSystematics[kCurrentResult][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(drawingRange.first, drawingRange.second);
        hEnergyEnergyCorrelatorRatioSystematics[kCurrentResult][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
          
        // Draw first the uncertainties, and then the statistical uncertainties
        drawer->DrawHistogram(hEnergyEnergyCorrelatorRatioSystematics[kCurrentResult][iCentrality][iJetPt][iTrackPt], "#Deltar", "PbPb EEC / pp EEC", " ", "e3");
        hEnergyEnergyCorrelatorRatioSystematics[kProjection][iCentrality][iJetPt][iTrackPt]->Draw("same,e3");

        for(int iFile = 0; iFile < knDataTypes; iFile++){
          hEnergyEnergyCorrelatorRatio[iFile][iCentrality][iJetPt][iTrackPt]->Draw("same");
        }

        // Add legends for drawn histograms
        legend->AddEntry(hEnergyEnergyCorrelatorRatio[kCurrentResult][iCentrality][iJetPt][iTrackPt], "CMS-PAS-HIN-23-004", "lpf");
        legend->AddEntry(hEnergyEnergyCorrelatorRatio[kProjection][iCentrality][iJetPt][iTrackPt], "Run4 projection", "lpf");
  
        // Draw the legend
        legend->Draw();
          
        // Save the figures to a file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("run4projection/run4projectionEECratio%s%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
        }
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

}
