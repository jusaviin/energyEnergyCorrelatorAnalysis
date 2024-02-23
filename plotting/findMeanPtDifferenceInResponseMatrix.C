#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"
#include "../src/EECHistograms.h"

/*
 * Macro for finding how much mean pT changes when unfolding with a response matrix
 */
void findMeanPtDifferenceInResponseMatrix(){
  
  // Input file
  TString inputFileName = "data/PbPbMC2018_GenGen_akFlowJets_4pCentShift_cutBadPhi_optimizedUnfoldingBins_nominalSmear_oneDimensionalResponseMatrix_processed_2024-01-16.root";
  // PbPbMC2018_GenGen_akFlowJets_4pCentShift_cutBadPhi_optimizedUnfoldingBins_nominalSmear_oneDimensionalResponseMatrix_processed_2024-01-16.root
  // PbPbMC2018_GenGen_akFlowJets_4pCentShift_cutBadPhi_optimizedUnfoldingBins_energyWeightSquared_nominalSmear_oneDimensionalResponseMatrix_processed_2024-01-18.root

  TFile* inputFile = TFile::Open(inputFileName);

  // Check that the input file exists
  if(inputFile == NULL){
    cout << "Error! The file " << inputFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }

  // Read EECCard from the file
  EECCard* card = new EECCard(inputFile);
  
  // ====================================================
  //               Binning configuration
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = card->GetNCentralityBins();
  const int nJetPtBinsEEC = card->GetNJetPtBinsEEC();
  
  // Select explicitly which bins from the files are compared:
  std::vector<std::pair<double,double>> comparedCentralityBin;
  comparedCentralityBin.push_back(std::make_pair(4,14));
  comparedCentralityBin.push_back(std::make_pair(14,34));
  comparedCentralityBin.push_back(std::make_pair(34,54));
  comparedCentralityBin.push_back(std::make_pair(54,94));

  std::vector<std::pair<double,double>> comparedJetPtBin;
  comparedJetPtBin.push_back(std::make_pair(120,140));
  comparedJetPtBin.push_back(std::make_pair(140,160));
  comparedJetPtBin.push_back(std::make_pair(160,180));
  comparedJetPtBin.push_back(std::make_pair(180,200));

  const int nStudiedJetPtBins = comparedJetPtBin.size();
  
  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager* histograms;
  histograms = new EECHistogramManager(inputFile, card);

  // Choose the energy-energy correlator types to load
  histograms->SetLoadJetHistograms(true);
  histograms->SetLoadJetPtOneDimensionalUnfoldingHistograms(true);

  // Choose the bin ranges
  histograms->SetCentralityBinRange(0, card->GetNCentralityBins() - 1);
  histograms->SetJetPtBinRangeEEC(0, card->GetNJetPtBinsEEC() - 1);

  // Load the histograms from the file
  histograms->LoadProcessedHistograms();

  // Energy-energy correlator histograms
  TH1D* hJetPt[nCentralityBins];
  TH2D* hSimpleResponseMatrix[nCentralityBins];
  TH1D* hProjectionForMeanFit[nCentralityBins][nStudiedJetPtBins];
  TProfile* hMeasuredMean[nCentralityBins];

  // Arrays for collections numbers
  double trueMean[nCentralityBins][nStudiedJetPtBins];
  double measuredMean[nCentralityBins][nStudiedJetPtBins];
  double fittedMean[nCentralityBins][nStudiedJetPtBins];
  
  // Initialize the energy-energy correlator histogram arrays to NULL
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    hJetPt[iCentrality] = NULL;
    hSimpleResponseMatrix[iCentrality] = NULL;
    hMeasuredMean[iCentrality] = NULL;
    for(int iJetPt = 0; iJetPt < nStudiedJetPtBins; iJetPt++){
      hProjectionForMeanFit[iCentrality][iJetPt] = NULL;
      trueMean[iCentrality][iJetPt] = 0;
      measuredMean[iCentrality][iJetPt] = 0;
      fittedMean[iCentrality][iJetPt] = 0;
    }
  } // Centrality loop
  
  // Helper variables
  double epsilon = 0.001;
  int iCentrality, iJetPt;
  int lowIntegralBin, highIntegralBin;

  // JDrawer for checking that fits are reasonable
  JDrawer* drawer = new JDrawer();
  TLegend* legend;
  TString centralityString, compactCentralityString;
  TString jetPtString, compactJetPtString;

  // Get the histograms from the histogram manager and see how much the mean jet pT shifts with unfolding
  for(auto centralityBin : comparedCentralityBin){
    iCentrality = card->FindBinIndexCentrality(centralityBin);
    centralityString = Form("Cent: %.0f-%.0f%%", centralityBin.first, centralityBin.second);
    compactCentralityString = Form("_C=%.0f-%.0f", centralityBin.first, centralityBin.second);

    // Load the jet pT histogram for this centrality range
    hJetPt[iCentrality] = histograms->GetHistogramJetPt(iCentrality);

    // Load the one dimensional jet pT response matrix
    hSimpleResponseMatrix[iCentrality] = histograms->GetHistogramJetPtOneDimensionalUnfoldingResponse(iCentrality);

    // Get the profile from the response matrix to easily determine measured mean pT values
    hMeasuredMean[iCentrality] = hSimpleResponseMatrix[iCentrality]->ProfileX();

    // Reset the jet pT counter
    iJetPt = 0;

    // Loop over the selected jet pT ranges to see how much the mean jet pT is expected to change
    for(auto jetPtBin : comparedJetPtBin){

      jetPtString = Form("%.0f < jet p_{T} < %.0f GeV", jetPtBin.first, jetPtBin.second);
      compactJetPtString = Form("_J=%.0f-%.0f", jetPtBin.first, jetPtBin.second);

      // Calculate the true mean jet pT in a selected bin. This is the mean after unfolding.
      lowIntegralBin = hJetPt[iCentrality]->FindBin(jetPtBin.first + epsilon);
      highIntegralBin = hJetPt[iCentrality]->FindBin(jetPtBin.second - epsilon);
      hJetPt[iCentrality]->GetXaxis()->SetRange(lowIntegralBin, highIntegralBin);
      trueMean[iCentrality][iJetPt] = hJetPt[iCentrality]->GetMean();

      // Project the selected reconstructed pT bin from the response matrix
      lowIntegralBin = hSimpleResponseMatrix[iCentrality]->GetXaxis()->FindBin(jetPtBin.first + epsilon);
      hProjectionForMeanFit[iCentrality][iJetPt] = hSimpleResponseMatrix[iCentrality]->ProjectionY(Form("projection%d%d", iCentrality, iJetPt), lowIntegralBin, lowIntegralBin);
      hProjectionForMeanFit[iCentrality][iJetPt]->GetXaxis()->SetRangeUser(25,300);
      hProjectionForMeanFit[iCentrality][iJetPt]->Fit("gaus","0");
      drawer->DrawHistogram(hProjectionForMeanFit[iCentrality][iJetPt], "Gen p_{T}", "Frequency", " ");
      hProjectionForMeanFit[iCentrality][iJetPt]->GetFunction("gaus")->Draw("same");

      legend = new TLegend(0.18, 0.7, 0.58, 0.85);
      legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);
      legend->AddEntry((TObject*)0, centralityString.Data(), "");
      legend->AddEntry((TObject*)0, jetPtString.Data(), "");
      legend->Draw();

      // Read the projected mean value
      fittedMean[iCentrality][iJetPt] = hProjectionForMeanFit[iCentrality][iJetPt]->GetFunction("gaus")->GetParameter(1);

      // Calculate the measured mean jet pT in a selected bin. This is the mean before unfolding.
      measuredMean[iCentrality][iJetPt] = hMeasuredMean[iCentrality]->GetBinContent(hMeasuredMean[iCentrality]->FindBin(jetPtBin.first + epsilon));

      // Increment the jet pT counter
      iJetPt++;

    } // Jet pT loop
  } // Centrality loop

  // Print the results to console
  for(auto centralityBin : comparedCentralityBin){
    iCentrality = card->FindBinIndexCentrality(centralityBin);

    // Reset the jet pT counter
    iJetPt = 0;

    // Loop over the selected jet pT ranges to see how much the mean jet pT is expected to change
    for(auto jetPtBin : comparedJetPtBin){

      // Printing the stuff:
      cout << "Centrality bin: " << centralityBin.first << "-" << centralityBin.second << "  Jet pT bin: " << jetPtBin.first << "-" << jetPtBin.second << endl;
      cout << "True mean: " << trueMean[iCentrality][iJetPt] << "  Measured mean: " << measuredMean[iCentrality][iJetPt] << "  Fitted mean: " << fittedMean[iCentrality][iJetPt] << "  Difference to measured: " << trueMean[iCentrality][iJetPt] - measuredMean[iCentrality][iJetPt] << endl;

      // Increment the jet pT counter
      iJetPt++;

    } // Jet pT loop
  } // Centrality loop


}
