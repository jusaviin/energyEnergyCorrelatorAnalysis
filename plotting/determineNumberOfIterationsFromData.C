#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Read data files with different number of iterations already performed
 * Take the chi2 between two consequtive iterations
 * When things stabilize, we have found a good iteration
 *
 *  int iSystem = Studied system (0 = PbPb, 1 = pp)
 *. int weightExponent = Weight exponent for energy-energy correlators
 */
void determineNumberOfIterationsFromData(int iSystem = 0, int weightExponent = 2){

  // Check that use input is within reasonable range
  enum enumSystemindex{kPbPb, kPp, knSystems};
  const int minWeightExponent = 1;
  const int maxWeightExponent = 2;
  const int nIterations = 8;

  TString systemString[] = {"PbPb", "pp"};
  TString energyWeightString[] = {"nominalEnergyWeight", "energyWeightSquared"}; 
  TString energyWeightLegend[] = {"n=1", "n=2"};
  TString backgroundString[] = {"_allBackgrounds",""};
  TString fileString[] = {"eecAnalysis_akFlowJet","ppData_pfJets_wtaAxis"};

  if(iSystem >= knSystems || iSystem < 0){
    cout << "ERROR! Index of system is out of range! Please give one of the followinf inputs:" << endl;
    cout << "0 for PbPb, or 1 for pp." << endl;
    return;
  }

  if(weightExponent < minWeightExponent || weightExponent > maxWeightExponent){
    cout << "ERROR! Weight exponent is out of range! Current implemented exponents are:" << endl;
    cout << minWeightExponent << " - " << maxWeightExponent << endl; 
    return;
  }

  // Enumeration variables
  TString unfoldName = "D'Agostini, NIT iterations";

  // **********************************
  //       Open the input files
  // **********************************

  // Define the name for the file containing histograms needed for unfolding
  TString inputFileName[knSystems][maxWeightExponent - minWeightExponent + 1][nIterations];
  for(int jSystem = 0; jSystem < knSystems; jSystem++){
    for(int jWeightExponent = 0; jWeightExponent < maxWeightExponent - minWeightExponent + 1; jWeightExponent++){
      for(int jIteration = 0; jIteration < nIterations; jIteration++){
        inputFileName[jSystem][jWeightExponent][jIteration] = Form("data/%s_%s%s_iteration%d_processed_2024-06-24.root", fileString[jSystem].Data(), energyWeightString[jWeightExponent].Data(), backgroundString[jSystem].Data(), jIteration+1);
      }
    }
  }


  TFile* inputFile[nIterations];
  for(int iIteration = 0; iIteration < nIterations; iIteration++){
    inputFile[iIteration] = TFile::Open(inputFileName[iSystem][weightExponent-1][iIteration]);

    if(inputFile[iIteration] == NULL) {
      cout << "Error! The file " << inputFileName[iSystem][weightExponent-1][iIteration].Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }

  }

  // Card for input files
  EECCard* dataCard = new EECCard(inputFile[0]);

  // ********************************************************
  //       Binning configuration for the unfolding study
  // ********************************************************
  
  // Find the number of bins from the card
  const int nCentralityBins = (iSystem == kPbPb) ? dataCard->GetNCentralityBins() : 1;
  const int nTrackPtBinsEEC = dataCard->GetNTrackPtBinsEEC();
  const int nJetPtBins = dataCard->GetNJetPtBinsEEC();

  // Select explicitly which bins from the files are compared:
  std::vector<std::pair<double,double>> comparedCentralityBin;
  comparedCentralityBin.push_back(std::make_pair(0,10));
  comparedCentralityBin.push_back(std::make_pair(10,30));
  comparedCentralityBin.push_back(std::make_pair(30,50));
  comparedCentralityBin.push_back(std::make_pair(50,90));

  std::vector<std::pair<double,double>> comparedJetPtBin;
  comparedJetPtBin.push_back(std::make_pair(120,140));
  comparedJetPtBin.push_back(std::make_pair(140,160));
  comparedJetPtBin.push_back(std::make_pair(160,180));
  comparedJetPtBin.push_back(std::make_pair(180,200));

  std::vector<double> comparedTrackPtBin;
  comparedTrackPtBin.push_back(1.0);
  //comparedTrackPtBin.push_back(1.5);
  comparedTrackPtBin.push_back(2.0);
  //comparedTrackPtBin.push_back(2.5);
  //comparedTrackPtBin.push_back(3.0);

  // If we are dealing with pp data, reset the centrality vector
  if(iSystem == kPp){
    comparedCentralityBin.clear();
    comparedCentralityBin.push_back(std::make_pair(-1,100));
  }

  // Select which histograms to draw
  bool drawChi2Histograms = false;
  bool drawPvalueHistograms = true;

  // Figure saving
  const bool saveFigures = true;  // Save figures
  TString saveComment = "";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures

  // Add energy weight string to save comment
  saveComment.Prepend(Form("_%s", energyWeightString[weightExponent-1].Data()));

  // Find todays date from the oracle
  AlgorithmLibrary* theOracle = new AlgorithmLibrary();
  TString today = theOracle->GetToday();

  // ***************************************************************
  //    Create histogram managers and load the needed histograms
  // ***************************************************************

  // Load the data histograms to be unfolded from the data histogram manager
  EECHistogramManager* histograms[nIterations];
  for(int iIteration = 0; iIteration < nIterations; iIteration++){
    histograms[iIteration] = new EECHistogramManager(inputFile[iIteration]);
    histograms[iIteration]->SetLoadEnergyEnergyCorrelators(true);
    histograms[iIteration]->SetCentralityBinRange(0,nCentralityBins-1);
    histograms[iIteration]->SetTrackPtBinRangeEEC(0,nTrackPtBinsEEC-1);
    histograms[iIteration]->SetJetPtBinRangeEEC(0, nJetPtBins-1);
    histograms[iIteration]->LoadProcessedHistograms();
  }

  // Histograms that are needed to create the unfolding response
  TH1D* energyEnergyCorrelatorUnfolded[nCentralityBins][nJetPtBins][nTrackPtBinsEEC][nIterations+1];
  TH1D* energyEnergyCorrelatorIterationRatio[nCentralityBins][nJetPtBins][nTrackPtBinsEEC][nIterations];

  // Histogram for chi2 change between two iterations
  TH1D* hChi2map[nCentralityBins][nJetPtBins][nTrackPtBinsEEC];
  TH1D* hPvalue[nCentralityBins][nJetPtBins][nTrackPtBinsEEC];

  // Initialize all the histograms to null
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        hChi2map[iCentrality][iJetPt][iTrackPt] = NULL;
        hPvalue[iCentrality][iJetPt][iTrackPt] = NULL;
        
        for(int iIteration = 0; iIteration < nIterations; iIteration++){
          energyEnergyCorrelatorUnfolded[iCentrality][iJetPt][iTrackPt][iIteration] = NULL;
          energyEnergyCorrelatorIterationRatio[iCentrality][iJetPt][iTrackPt][iIteration] = NULL;
        }
        energyEnergyCorrelatorUnfolded[iCentrality][iJetPt][iTrackPt][nIterations] = NULL;
      } // Track pT loop
    } // Jet pT loop
  }  // Centrality loop

  // Normalization region for histograms
  double normalizationRegionLow = 0.008;
  double normalizationRegionHigh = 0.39;
  int firstNormalizationBin, lastNormalizationBin;

  // Read the energy-energy correlator histograms for different numbers of iterations
  int iCentrality, iJetPt, iTrackPt;
  for(auto centralityBin : comparedCentralityBin){
    iCentrality = dataCard->FindBinIndexCentrality(centralityBin);
    for(auto jetPtBin : comparedJetPtBin){
      iJetPt = dataCard->FindBinIndexJetPtEEC(jetPtBin);
      for(auto trackPtBin : comparedTrackPtBin){
        iTrackPt = dataCard->GetBinIndexTrackPtEEC(trackPtBin);

        energyEnergyCorrelatorUnfolded[iCentrality][iJetPt][iTrackPt][0] = histograms[0]->GetHistogramEnergyEnergyCorrelator(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt);
        firstNormalizationBin = energyEnergyCorrelatorUnfolded[iCentrality][iJetPt][iTrackPt][0]->GetXaxis()->FindBin(normalizationRegionLow);
        lastNormalizationBin = energyEnergyCorrelatorUnfolded[iCentrality][iJetPt][iTrackPt][0]->GetXaxis()->FindBin(normalizationRegionHigh);
        energyEnergyCorrelatorUnfolded[iCentrality][iJetPt][iTrackPt][0]->Scale(1.0 / energyEnergyCorrelatorUnfolded[iCentrality][iJetPt][iTrackPt][0]->Integral(firstNormalizationBin, lastNormalizationBin, "width"));

        for(int iIteration = 0; iIteration < nIterations; iIteration++){
          energyEnergyCorrelatorUnfolded[iCentrality][iJetPt][iTrackPt][iIteration+1] = histograms[iIteration]->GetHistogramEnergyEnergyCorrelatorProcessed(EECHistogramManager::kEnergyEnergyCorrelator, iCentrality, iJetPt, iTrackPt, EECHistogramManager::kEnergyEnergyCorrelatorUnfolded);
          energyEnergyCorrelatorUnfolded[iCentrality][iJetPt][iTrackPt][iIteration+1]->Scale(1.0 / energyEnergyCorrelatorUnfolded[iCentrality][iJetPt][iTrackPt][iIteration+1]->Integral(firstNormalizationBin, lastNormalizationBin, "width"));

          // Calculate the ratio between this iteration and the previous iteration
          energyEnergyCorrelatorIterationRatio[iCentrality][iJetPt][iTrackPt][iIteration] = (TH1D*) energyEnergyCorrelatorUnfolded[iCentrality][iJetPt][iTrackPt][iIteration+1]->Clone(Form("iterationRatio%d%d%d%D", iCentrality, iJetPt, iTrackPt, iIteration));
          energyEnergyCorrelatorIterationRatio[iCentrality][iJetPt][iTrackPt][iIteration]->Divide(energyEnergyCorrelatorUnfolded[iCentrality][iJetPt][iTrackPt][iIteration]);

        } // Iteration loop
      } // Track pT loop
    } // Jet pT loop
  }  // Centrality loop

  // Create new histograms for chi2 maps
  for(auto centralityBin : comparedCentralityBin){
    iCentrality = dataCard->FindBinIndexCentrality(centralityBin);
    for(auto jetPtBin : comparedJetPtBin){
      iJetPt = dataCard->FindBinIndexJetPtEEC(jetPtBin);
      for(auto trackPtBin : comparedTrackPtBin){
        iTrackPt = dataCard->GetBinIndexTrackPtEEC(trackPtBin);
        hChi2map[iCentrality][iJetPt][iTrackPt] = new TH1D(Form("chi2map%d%d%d", iCentrality, iTrackPt, iJetPt), Form("chi2map%d%d%d", iCentrality, iTrackPt, iJetPt), nIterations, 0.5, 0.5+nIterations);
        hPvalue[iCentrality][iJetPt][iTrackPt] = new TH1D(Form("pValue%d%d%d", iCentrality, iTrackPt, iJetPt), Form("pValue%d%d%d", iCentrality, iTrackPt, iJetPt), nIterations, 0.5, 0.5+nIterations);
      } // Track pT loop
    } // Jet pT loop
  }  // Centrality loop

  // After all the ratios have been calculated, fit the iteration to iteration ratios with a constant and read the chi2 of the fit
  TF1* testFunction = new TF1("testFunction","1",normalizationRegionLow,normalizationRegionHigh);
  double chiSquare;
  int ndf = lastNormalizationBin - firstNormalizationBin + 1;
  double pValue;
  for(auto centralityBin : comparedCentralityBin){
    iCentrality = dataCard->FindBinIndexCentrality(centralityBin);
    for(auto jetPtBin : comparedJetPtBin){
      iJetPt = dataCard->FindBinIndexJetPtEEC(jetPtBin);
      for(auto trackPtBin : comparedTrackPtBin){
        iTrackPt = dataCard->GetBinIndexTrackPtEEC(trackPtBin);
        for(int iIteration = 0; iIteration < nIterations; iIteration++){

          chiSquare = energyEnergyCorrelatorIterationRatio[iCentrality][iJetPt][iTrackPt][iIteration]->Chisquare(testFunction, "R");
          hChi2map[iCentrality][iJetPt][iTrackPt]->SetBinContent(iIteration+1,chiSquare);

          pValue = TMath::Prob(chiSquare, ndf);
          hPvalue[iCentrality][iJetPt][iTrackPt]->SetBinContent(iIteration+1,pValue);

        } // Iteration loop
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop

  // Draw the chi2 maps to canvases
  JDrawer* drawer = new JDrawer();
  TLegend* legend;

  TString centralityString, trackPtString, jetPtString;
  TString compactCentralityString, compactTrackPtString, compactJetPtString;

  TLine* veryGoodLine = new TLine(0.5, 0.05, nIterations+0.5, 0.05);
  veryGoodLine->SetLineColor(kRed);
  veryGoodLine->SetLineStyle(2);

  TLine* oneLine = new TLine(0.5, 1, nIterations+0.5, 1);
  oneLine->SetLineColor(kBlack);
  oneLine->SetLineStyle(2);

  TLine* zeroLine = new TLine(0.5, 0, nIterations+0.5, 0);
  zeroLine->SetLineColor(kBlack);
  zeroLine->SetLineStyle(2);

  if(drawChi2Histograms){

    for(auto centralityBin : comparedCentralityBin){

      if(iSystem == kPp){
        iCentrality = 0;
        centralityString = "pp";
        compactCentralityString = "_pp";
      } else {
        iCentrality = dataCard->FindBinIndexCentrality(centralityBin);
        centralityString = Form("PbPb %.0f-%.0f%%", centralityBin.first, centralityBin.second);
        compactCentralityString = Form("_C=%.0f-%.0f", centralityBin.first, centralityBin.second);
      }

      for(auto jetPtBin : comparedJetPtBin){
        iJetPt = dataCard->FindBinIndexJetPtEEC(jetPtBin);
        jetPtString = Form("%.0f < jet p_{T} < %.0f", jetPtBin.first, jetPtBin.second);
        compactJetPtString = Form("_J=%.0f-%.0f", jetPtBin.first, jetPtBin.second);

        for(auto trackPtBin : comparedTrackPtBin){
          iTrackPt = dataCard->GetBinIndexTrackPtEEC(trackPtBin);
          trackPtString = Form("%.1f < track p_{T}", trackPtBin);
          compactTrackPtString = Form("_T>%.1f", trackPtBin);
          compactTrackPtString.ReplaceAll(".","v");

          // Create a legend for the plot
          legend = new TLegend(0.23, 0.05, 0.53, 0.6);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);
          legend->AddEntry((TObject*)0, centralityString.Data(), "");
          legend->AddEntry((TObject*)0, jetPtString.Data(), "");
          legend->AddEntry((TObject*)0, trackPtString.Data(), "");
          legend->AddEntry((TObject*)0, energyWeightLegend[weightExponent-1].Data(), "");
          
          // Draw the histograms to the upper canvas
          drawer->DrawHistogram(hChi2map[iCentrality][iJetPt][iTrackPt], "Iteration", "#chi^{2}", " ");

          // Draw the legend to the upper pad
          legend->Draw();

          // Save the figures to a file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/unfoldingIterationChi2%s%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
          }

        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Drawing chi2 histograms

  if(drawPvalueHistograms){

    for(auto centralityBin : comparedCentralityBin){
      
      if(iSystem == kPp){
        iCentrality = 0;
        centralityString = "pp";
        compactCentralityString = "_pp";
      } else {
        iCentrality = dataCard->FindBinIndexCentrality(centralityBin);
        centralityString = Form("PbPb %.0f-%.0f%%", centralityBin.first, centralityBin.second);
        compactCentralityString = Form("_C=%.0f-%.0f", centralityBin.first, centralityBin.second);
      }

      for(auto jetPtBin : comparedJetPtBin){
        iJetPt = dataCard->FindBinIndexJetPtEEC(jetPtBin);
        jetPtString = Form("%.0f < jet p_{T} < %.0f", jetPtBin.first, jetPtBin.second);
        compactJetPtString = Form("_J=%.0f-%.0f", jetPtBin.first, jetPtBin.second);

        for(auto trackPtBin : comparedTrackPtBin){
          iTrackPt = dataCard->GetBinIndexTrackPtEEC(trackPtBin);
          trackPtString = Form("%.1f < track p_{T}", trackPtBin);
          compactTrackPtString = Form("_T>%.1f", trackPtBin);
          compactTrackPtString.ReplaceAll(".","v");

          // Create a legend for the plot
          legend = new TLegend(0.53, 0.25, 0.83, 0.6);
          legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);
          legend->AddEntry((TObject*)0, centralityString.Data(), "");
          legend->AddEntry((TObject*)0, jetPtString.Data(), "");
          legend->AddEntry((TObject*)0, trackPtString.Data(), "");
          legend->AddEntry((TObject*)0, energyWeightLegend[weightExponent-1].Data(), "");
          
          // Set a nice style for the histogram
          hPvalue[iCentrality][iJetPt][iTrackPt]->SetMarkerStyle(kFullCircle);
          hPvalue[iCentrality][iJetPt][iTrackPt]->SetMarkerColor(kBlack);

          // Set a nice drawing range
          hPvalue[iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(-0.02, 1.02);

          // Draw the histograms to the upper canvas
          drawer->DrawHistogram(hPvalue[iCentrality][iJetPt][iTrackPt], "Iteration", "P-value", " ", "p");

          // Draw a line showing p=0.05
          veryGoodLine->Draw();
          oneLine->Draw();
          zeroLine->Draw();

          // Draw the legend to the upper pad
          legend->Draw();

          // Save the figures to a file
          if(saveFigures){
            gPad->GetCanvas()->SaveAs(Form("figures/unfoldingIterationPvalues%s%s%s%s.%s", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
          }

        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Drawing chi2 histograms



}
