#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Macro studying specifically the fake+fake background in.
 * It is assumed in the code that the centrality and track pT bins match between the files
 */
void fakeFakeBackgroundPlotter(){

  enum enumDataType{kPythiaHydjet, kMinBiasHydjet, knDataTypes};
  
  // File containing the Pythia+Hydjet simulation result (index 0), and the one containing minimum bias Hydjet result (index 1)
  TString inputFileName[knDataTypes] = {"data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJet_fakeFakeReflectedCone_noTrigger_preprocessed_2022-09-23.root", "data/MinBiasHydjet_RecoGen_eecAnalysis_akFlowJet_firstMinBiasScan_noTrigger_preprocessed_2022-10-10.root"};
  // data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJet_fakeFakeReflectedCone_noTrigger_preprocessed_2022-09-23.root
  // data/PbPbMC2018_GenGen_eecAnalysis_genJet_fakeFakeReflectedCone_noTrigger_preprocessed_2022-09-30.root
  // "data/MinBiasHydjet_RecoGen_eecAnalysis_akFlowJet_firstMinBiasScan_noTrigger_preprocessed_2022-10-10.root"
  
  // Open the input files and read the cards
  TFile* inputFile[knDataTypes];
  EECCard* card[knDataTypes];
  for(int iDataType = 0; iDataType < knDataTypes; iDataType++){
    inputFile[iDataType] = TFile::Open(inputFileName[iDataType]);
    
    // Check that the files exist
    if(inputFile[iDataType] == NULL){
      cout << "Error! The file " << inputFileName[iDataType].Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }
    
    // Load the card from the file and read the collision system
    card[iDataType] = new EECCard(inputFile[iDataType]);
    
  }
    
  // Require matching centrality and track pT bins
  const double epsilon = 0.00001;
  
  const int nCentralityBins = card[kPythiaHydjet]->GetNCentralityBins();
  if(nCentralityBins != card[kMinBiasHydjet]->GetNCentralityBins()){
    cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
    return;
  }
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    if(TMath::Abs(card[kPythiaHydjet]->GetLowBinBorderCentrality(iCentrality) - card[kMinBiasHydjet]->GetLowBinBorderCentrality(iCentrality)) > epsilon){
      cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
    if(TMath::Abs(card[kPythiaHydjet]->GetHighBinBorderCentrality(iCentrality) - card[kMinBiasHydjet]->GetHighBinBorderCentrality(iCentrality)) > epsilon){
      cout << "Error! Centrality bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
  }
  
  const int nTrackPtBinsEEC = card[kPythiaHydjet]->GetNTrackPtBinsEEC();
  if(nTrackPtBinsEEC != card[kMinBiasHydjet]->GetNTrackPtBinsEEC()){
    cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
    return;
  }
  for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
    if(TMath::Abs(card[kPythiaHydjet]->GetLowBinBorderTrackPtEEC(iTrackPt) - card[kMinBiasHydjet]->GetLowBinBorderTrackPtEEC(iTrackPt)) > epsilon){
      cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
    if(TMath::Abs(card[kPythiaHydjet]->GetHighBinBorderTrackPtEEC(iTrackPt) - card[kMinBiasHydjet]->GetHighBinBorderTrackPtEEC(iTrackPt)) > epsilon){
      cout << "Error! Track pT bins do not match between the two files! Cannot execute the code!" << endl;
      return;
    }
  }
  
  // ====================================================
  //    Binning configuration for the background study
  // ====================================================
  
  // Find the number of jet pT bins
  const int nJetPtBinsEEC[knDataTypes] = {card[kPythiaHydjet]->GetNJetPtBinsEEC(), card[kMinBiasHydjet]->GetNJetPtBinsEEC()};
  
  // Default binning ranges for reference
  // centrality = {0,10,30,50,90}
  // track pT =  {0.7,1,2,3,4,8,12,300}
  // jet pT in energy-energy correlators = {120,140,160,180,200,300,500,5030}
  // track pT in energy-energy correlators = {0.7,1,2,3,4,6,300}
  
  // Bin range to be integrated
  int firstStudiedCentralityBin = 0;
  int lastStudiedCentralityBin = 0;
  
  int firstStudiedJetPtBinEEC[knDataTypes] = {nJetPtBinsEEC[kPythiaHydjet],nJetPtBinsEEC[kMinBiasHydjet]};
  int lastStudiedJetPtBinEEC[knDataTypes] = {nJetPtBinsEEC[kPythiaHydjet],nJetPtBinsEEC[kMinBiasHydjet]}; // Note: Jets integrated over all pT ranges are in nJetPtBinsEEC bin
  
  int firstStudiedTrackPtBinEEC = 5;
  int lastStudiedTrackPtBinEEC = 5;
  
  // Select the types of energy-energy correlators are studied
  bool studyEnergyEnergyCorrelator[EECHistogramManager::knEnergyEnergyCorrelatorTypes];
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator] = true;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPt] = false;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorUncorrected] = false;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPtUncorrected] = false;
  
  // Select which plots to draw
  const bool drawFakeFakeForConstantJetPt = false;
  const bool drawFakeFakeForConstantTrackPt = false;
  const bool drawMinBiasToFakeFakeComparison = true;
  
  // Logarithmic axes
  const bool logDeltaR = true;
  const bool logEEC = true;
  
  // Axis zooming
  std::pair<double,double> ratioZoom = std::make_pair(0, 2);
  
  // Figure saving
  const bool saveFigures = false;  // Save figures
  const char* saveComment = "";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures
  
  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager *histograms[knDataTypes];
  for(int iDataType = 0; iDataType < knDataTypes; iDataType++){
    histograms[iDataType] = new EECHistogramManager(inputFile[iDataType],card[iDataType]);
    
    // Choose the energy-energy correlator types to load
    histograms[iDataType]->SetLoadEnergyEnergyCorrelators(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator]);
    histograms[iDataType]->SetLoadEnergyEnergyCorrelatorsJetPt(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPt]);
    histograms[iDataType]->SetLoadEnergyEnergyCorrelatorsUncorrected(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorUncorrected]);
    histograms[iDataType]->SetLoadEnergyEnergyCorrelatorsJetPtUncorrected(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPtUncorrected]);
    
    // Choose the bin ranges
    histograms[iDataType]->SetCentralityBinRange(firstStudiedCentralityBin,lastStudiedCentralityBin);
    histograms[iDataType]->SetJetPtBinRangeEEC(firstStudiedJetPtBinEEC[iDataType],lastStudiedJetPtBinEEC[iDataType]);
    histograms[iDataType]->SetTrackPtBinRangeEEC(firstStudiedTrackPtBinEEC,lastStudiedTrackPtBinEEC);
    
    // Load the histograms from the file
    histograms[iDataType]->LoadProcessedHistograms();
  }
  
  // Energy-energy correlator histograms separated by subevents from the Pythia+Hydjet simulation
  TH1D* hFakeFake[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC[kPythiaHydjet]+1][nTrackPtBinsEEC];
  
  // Reflected cone energy-energy correlators
  TH1D* hMinBias[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC[kMinBiasHydjet]+1][nTrackPtBinsEEC];
  
  // Histograms for all different ratios
  TH1D* hFakeFakeRatioConstantTrackPt[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC[kPythiaHydjet]+1][nTrackPtBinsEEC];
  TH1D* hFakeFakeRatioConstantJetPt[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC[kPythiaHydjet]+1][nTrackPtBinsEEC];
  TH1D* hMinBiasToFakeFakeRatio[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC[kPythiaHydjet]+1][nJetPtBinsEEC[kMinBiasHydjet]+1][nTrackPtBinsEEC];
  
  // Initialize the energy-energy correlator histogram array to NULL
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        for(int iJetPtFakeFake = 0; iJetPtFakeFake < nJetPtBinsEEC[kPythiaHydjet]+1; iJetPtFakeFake++){
          hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPtFakeFake][iTrackPt] = NULL;
          hFakeFakeRatioConstantTrackPt[iEnergyEnergyCorrelator][iCentrality][iJetPtFakeFake][iTrackPt] = NULL;
          hFakeFakeRatioConstantJetPt[iEnergyEnergyCorrelator][iCentrality][iJetPtFakeFake][iTrackPt] = NULL;
          for(int iJetPtMinBias = 0; iJetPtMinBias < nJetPtBinsEEC[kMinBiasHydjet]+1; iJetPtMinBias++){
            hMinBiasToFakeFakeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPtFakeFake][iJetPtMinBias][iTrackPt] = NULL;
          }
        } // Pythia+Hydjet jet pT loop
        for(int iJetPtMinBias = 0; iJetPtMinBias < nJetPtBinsEEC[kMinBiasHydjet]+1; iJetPtMinBias++){
          hMinBias[iEnergyEnergyCorrelator][iCentrality][iJetPtMinBias][iTrackPt] = NULL;
        } // Hydjet only jet pT loop
      } // Track pT loop
    } // Centrality loop
  } // Energy-energy correlator loop
  
  // Find at least one energy-energy correlator index that is studied
  int studiedEnergyEnergyCorrelatorIndex = -1;
  
  // Get the histograms from the histogram manager
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only read the selected energy-energy correlator types
    if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
    studiedEnergyEnergyCorrelatorIndex = iEnergyEnergyCorrelator;
    
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
        for(int iJetPtFakeFake = firstStudiedJetPtBinEEC[kPythiaHydjet]; iJetPtFakeFake <= lastStudiedJetPtBinEEC[kPythiaHydjet]; iJetPtFakeFake++){
          
          // Read the fake+fake histogram and normalize everything to one
          hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPtFakeFake][iTrackPt] = histograms[kPythiaHydjet]->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPtFakeFake, iTrackPt, EECHistograms::kSameJetPair, EECHistogramManager::kHydjetHydjet);
          hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPtFakeFake][iTrackPt]->Scale(1/hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPtFakeFake][iTrackPt]->Integral("width"));
          
          // For constant jet pT ratios, calculate the ratio to the lowest shown track pT cut
          hFakeFakeRatioConstantJetPt[iEnergyEnergyCorrelator][iCentrality][iJetPtFakeFake][iTrackPt] = (TH1D*) hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPtFakeFake][iTrackPt]->Clone(Form("constantJetPtRatio%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPtFakeFake, iTrackPt));
          hFakeFakeRatioConstantJetPt[iEnergyEnergyCorrelator][iCentrality][iJetPtFakeFake][iTrackPt]->Divide(hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPtFakeFake][firstStudiedTrackPtBinEEC]);
          
          // For constant track pT ratios, calculate the ratio to the lowest shown jet pT cut
          hFakeFakeRatioConstantTrackPt[iEnergyEnergyCorrelator][iCentrality][iJetPtFakeFake][iTrackPt] = (TH1D*) hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPtFakeFake][iTrackPt]->Clone(Form("constantTrackPtRatio%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPtFakeFake, iTrackPt));
          hFakeFakeRatioConstantTrackPt[iEnergyEnergyCorrelator][iCentrality][iJetPtFakeFake][iTrackPt]->Divide(hFakeFake[iEnergyEnergyCorrelator][iCentrality][firstStudiedJetPtBinEEC[kPythiaHydjet]][iTrackPt]);

        } // Pythia+Hydjet jet pT loop
        
        for(int iJetPtMinBias = firstStudiedJetPtBinEEC[kMinBiasHydjet]; iJetPtMinBias <= lastStudiedJetPtBinEEC[kMinBiasHydjet]; iJetPtMinBias++){
          
          // Read the minimum bias histogram and normalize everything to one
          hMinBias[iEnergyEnergyCorrelator][iCentrality][iJetPtMinBias][iTrackPt] = histograms[kMinBiasHydjet]->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPtMinBias, iTrackPt, EECHistograms::kSameJetPair);
          hMinBias[iEnergyEnergyCorrelator][iCentrality][iJetPtMinBias][iTrackPt]->Scale(1/hMinBias[iEnergyEnergyCorrelator][iCentrality][iJetPtMinBias][iTrackPt]->Integral("width"));
          
          // Calculate the min bias to fake+fake ratio
          for(int iJetPtFakeFake = firstStudiedJetPtBinEEC[kPythiaHydjet]; iJetPtFakeFake <= lastStudiedJetPtBinEEC[kPythiaHydjet]; iJetPtFakeFake++){
            hMinBiasToFakeFakeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPtFakeFake][iJetPtMinBias][iTrackPt] = (TH1D*) hMinBias[iEnergyEnergyCorrelator][iCentrality][iJetPtMinBias][iTrackPt]->Clone(Form("minBiasToFakeFakeRatio%d%d%d%d%d", iEnergyEnergyCorrelator, iCentrality, iJetPtFakeFake, iJetPtMinBias, iTrackPt));
            hMinBiasToFakeFakeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPtFakeFake][iJetPtMinBias][iTrackPt]->Divide(hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPtFakeFake][iTrackPt]);
          } // Pythia+Hydjet jet pT loop
          
        } // Min bias Hydjet jet pT loop
        
      } // Track pT loop
    } // Centrality loop
  } // Energy-energy correlator loop
  
  // Find the binning information for DeltaR bins from a histogram
  const int nDeltaRBins = hFakeFake[studiedEnergyEnergyCorrelatorIndex][firstStudiedCentralityBin][firstStudiedJetPtBinEEC[kPythiaHydjet]][firstStudiedTrackPtBinEEC]->GetNbinsX();
  double deltaRBinBorders[nDeltaRBins+1];
  deltaRBinBorders[0] = hFakeFake[studiedEnergyEnergyCorrelatorIndex][firstStudiedCentralityBin][firstStudiedJetPtBinEEC[kPythiaHydjet]][firstStudiedTrackPtBinEEC]->GetXaxis()->GetBinLowEdge(1);
  for(int iBin = 1; iBin <= nDeltaRBins; iBin++){
    deltaRBinBorders[iBin] = hFakeFake[studiedEnergyEnergyCorrelatorIndex][firstStudiedCentralityBin][firstStudiedJetPtBinEEC[kPythiaHydjet]][firstStudiedTrackPtBinEEC]->GetXaxis()->GetBinUpEdge(iBin);
  }
  
  
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);

  if(logDeltaR) drawer->SetLogX(true);

  TLegend *legend, *ptLegend;
  TString centralityString, trackPtString, jetPtString, jetPtStringMinBias;
  TString compactCentralityString, compactTrackPtString, compactJetPtString, compactJetPtStringMinBias;
  int color[10] = {kBlack, kBlue, kRed, kGreen+2, kCyan, kMagenta, kOrange-1, kAzure+5, kOrange-2, kGray};

  // Change the legend size based on the number of histograms in the legend
  double trackLegendY1 = 0.03;
  double trackLegendY2 = trackLegendY1 + 0.06*(lastStudiedTrackPtBinEEC+1-firstStudiedTrackPtBinEEC);
  double jetLegendY1 = 0.03;
  double jetLegendY2 = trackLegendY1 + 0.06*(lastStudiedJetPtBinEEC[kPythiaHydjet]+1-firstStudiedJetPtBinEEC[kPythiaHydjet]);
  double minBiasLegendY1 = 0.03;
  double minBiasLegendY2 = minBiasLegendY1 + 0.06*(lastStudiedJetPtBinEEC[kMinBiasHydjet]+2-firstStudiedJetPtBinEEC[kMinBiasHydjet]);

  // Loop over all selected histograms
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only read the selected energy-energy correlator types
    if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
    
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      
      // Set the centrality information for legends and figure saving
      centralityString = Form("Cent: %.0f-%.0f%%", histograms[kPythiaHydjet]->GetCentralityBinBorder(iCentrality), histograms[kPythiaHydjet]->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f", histograms[kPythiaHydjet]->GetCentralityBinBorder(iCentrality), histograms[kPythiaHydjet]->GetCentralityBinBorder(iCentrality+1));
      
      // ====================================================
      // ===   Draw fake+fake plots for constant jet pT   ===
      // ====================================================
      if(drawFakeFakeForConstantJetPt){
        
        // Loop over jet pT bins
        for(int iJetPt = firstStudiedJetPtBinEEC[kPythiaHydjet]; iJetPt <= lastStudiedJetPtBinEEC[kPythiaHydjet]; iJetPt++){
          
          // Set the jet pT information for legends and figure saving
          if(iJetPt == histograms[kPythiaHydjet]->GetNJetPtBinsEEC()){
            jetPtString = Form("Jet p_{T} > %.0f", histograms[kPythiaHydjet]->GetCard()->GetJetPtCut());
            compactJetPtString = "";
          } else {
            jetPtString = Form("%.0f < jet p_{T} < %.0f", histograms[kPythiaHydjet]->GetJetPtBinBorderEEC(iJetPt), histograms[kPythiaHydjet]->GetJetPtBinBorderEEC(iJetPt+1));
            compactJetPtString = Form("_J=%.0f-%.0f", histograms[kPythiaHydjet]->GetJetPtBinBorderEEC(iJetPt), histograms[kPythiaHydjet]->GetJetPtBinBorderEEC(iJetPt+1));
          }
          
          // Create a legend for the figure
          legend = new TLegend(0.25,trackLegendY1,0.45,trackLegendY1+0.24);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          legend->AddEntry((TObject*) 0, histograms[kPythiaHydjet]->GetCard()->GetAlternativeDataType().Data(), "");
          legend->AddEntry((TObject*) 0, "Fake+fake pairs","");
          legend->AddEntry((TObject*) 0, centralityString.Data(),"");
          legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
          
          ptLegend = new TLegend(0.6,trackLegendY1,0.8,trackLegendY2);
          ptLegend->SetFillStyle(0);ptLegend->SetBorderSize(0);ptLegend->SetTextSize(0.05);ptLegend->SetTextFont(62);
          
          // Create a new canvas for the plot
          drawer->CreateSplitCanvas();
          
          // Logarithmic EEC axis
          if(logEEC) drawer->SetLogY(true);
          
          // Loop over track pT bins and draw the distributions to the upper pad
          for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
            
            trackPtString = Form("%.1f < track p_{T}",histograms[kPythiaHydjet]->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString = Form("_T%.1f",histograms[kPythiaHydjet]->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString.ReplaceAll(".","v");
            
            // For logarithmic drawing, cannot go down to zero
            if(logDeltaR){
              hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
            }
            
            // Set good y-ranges for plotting
            hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.001, 40);
            
            // Draw the background histogram to the upper canvas
            hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iTrackPt]);
            if(iTrackPt == firstStudiedTrackPtBinEEC){
              drawer->DrawHistogramToUpperPad(hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt], "#Deltar", histograms[kPythiaHydjet]->GetEnergyEnergyCorrelatorAxisName(iEnergyEnergyCorrelator), " ");
            } else {
              hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Draw("same");
            }
            ptLegend->AddEntry(hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt], trackPtString.Data(), "l");
            
          } // Track pT loop
          
          // Draw the legends to the upper pad
          legend->Draw();
          ptLegend->Draw();
          
          // Linear scale for the ratio
          drawer->SetLogY(false);
          
          // Loop again over the track pT bins and draw the ratios to the lower pad
          for(int iTrackPt = firstStudiedTrackPtBinEEC+1; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
            
            // For logarithmic drawing, cannot go down to zero
            if(logDeltaR){
              hFakeFakeRatioConstantJetPt[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
            }
            
            hFakeFakeRatioConstantJetPt[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iTrackPt]);
            if(iTrackPt == firstStudiedTrackPtBinEEC+1){
              drawer->SetGridY(true);
              hFakeFakeRatioConstantJetPt[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
                drawer->DrawHistogramToLowerPad(hFakeFakeRatioConstantJetPt[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt], "#Deltar", Form("#frac{Color}{%.1f < track p_{T}}",histograms[kPythiaHydjet]->GetTrackPtBinBorderEEC(firstStudiedTrackPtBinEEC)), " ");
              drawer->SetGridY(false);
            } else {
              hFakeFakeRatioConstantJetPt[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Draw("same");
            }
            
          } // Track pT loop
    
        } // Jet pT loop
      } // Drawing fake+fake plots for constant jet pT
      
      // ======================================================
      // ===   Draw fake+fake plots for constant track pT   ===
      // ======================================================
      if(drawFakeFakeForConstantTrackPt){
        
        // Loop over track pT bins
        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
          
          trackPtString = Form("%.1f < track p_{T}",histograms[kPythiaHydjet]->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString = Form("_T%.1f",histograms[kPythiaHydjet]->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString.ReplaceAll(".","v");
          
          // Create a legend for the figure
          legend = new TLegend(0.25,jetLegendY1,0.45,jetLegendY1+0.24);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          legend->AddEntry((TObject*) 0, histograms[kPythiaHydjet]->GetCard()->GetAlternativeDataType().Data(), "");
          legend->AddEntry((TObject*) 0, "Fake+fake pairs","");
          legend->AddEntry((TObject*) 0, centralityString.Data(),"");
          legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
          
          ptLegend = new TLegend(0.6,jetLegendY1,0.8,jetLegendY2);
          ptLegend->SetFillStyle(0);ptLegend->SetBorderSize(0);ptLegend->SetTextSize(0.05);ptLegend->SetTextFont(62);
          
          // Create a new canvas for the plot
          drawer->CreateSplitCanvas();
          
          // Logarithmic EEC axis
          if(logEEC) drawer->SetLogY(true);
          
          // Loop over jet pT bins and draw the distribution to the upped pad
          for(int iJetPt = firstStudiedJetPtBinEEC[kPythiaHydjet]; iJetPt <= lastStudiedJetPtBinEEC[kPythiaHydjet]; iJetPt++){
            
            // Set the jet pT information for legends and figure saving
            if(iJetPt == histograms[kPythiaHydjet]->GetNJetPtBinsEEC()){
              jetPtString = Form("Jet p_{T} > %.0f", histograms[kPythiaHydjet]->GetCard()->GetJetPtCut());
              compactJetPtString = "";
            } else {
              jetPtString = Form("%.0f < jet p_{T} < %.0f", histograms[kPythiaHydjet]->GetJetPtBinBorderEEC(iJetPt), histograms[kPythiaHydjet]->GetJetPtBinBorderEEC(iJetPt+1));
              compactJetPtString = Form("_J=%.0f-%.0f", histograms[kPythiaHydjet]->GetJetPtBinBorderEEC(iJetPt), histograms[kPythiaHydjet]->GetJetPtBinBorderEEC(iJetPt+1));
            }
            
            // For logarithmic drawing, cannot go down to zero
            if(logDeltaR){
              hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
            }
            
            // Set good y-ranges for plotting
            hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.001, 40);
            
            // Draw the background histogram to the upper canvas
            hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iJetPt]);
            if(iJetPt == firstStudiedJetPtBinEEC[kPythiaHydjet]){
              drawer->DrawHistogramToUpperPad(hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt], "#Deltar", histograms[kPythiaHydjet]->GetEnergyEnergyCorrelatorAxisName(iEnergyEnergyCorrelator), " ");
            } else {
              hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Draw("same");
            }
            ptLegend->AddEntry(hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt], jetPtString.Data(), "l");
            
          } // Jet pT loop
          
          // Draw the legends to the upper pad
          legend->Draw();
          ptLegend->Draw();
          
          // Linear scale for the ratio
          drawer->SetLogY(false);
          
          // Loop again over the jet pT bins and draw the ratios to the lower pad
          for(int iJetPt = firstStudiedJetPtBinEEC[kPythiaHydjet]+1; iJetPt <= lastStudiedJetPtBinEEC[kPythiaHydjet]; iJetPt++){
            
            // For logarithmic drawing, cannot go down to zero
            if(logDeltaR){
              hFakeFakeRatioConstantTrackPt[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
            }
            
            hFakeFakeRatioConstantTrackPt[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iJetPt]);
            if(iJetPt == firstStudiedJetPtBinEEC[kPythiaHydjet]+1){
              drawer->SetGridY(true);
              hFakeFakeRatioConstantTrackPt[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
                drawer->DrawHistogramToLowerPad(hFakeFakeRatioConstantTrackPt[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt], "#Deltar", Form("#frac{Color}{%.0f < jet p_{T} < %.0f}",histograms[kPythiaHydjet]->GetJetPtBinBorderEEC(firstStudiedJetPtBinEEC[kPythiaHydjet]), histograms[kPythiaHydjet]->GetJetPtBinBorderEEC(firstStudiedJetPtBinEEC[kPythiaHydjet]+1)), " ");
              drawer->SetGridY(false);
            } else {
              hFakeFakeRatioConstantTrackPt[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Draw("same");
            }
            
          } // Track pT loop
    
        } // Jet pT loop
      } // Drawing fake+fake plots for constant track pT
      
      // ======================================================
      // ===   Draw minimum bias to fake+fake comparison    ===
      // ======================================================
      if(drawMinBiasToFakeFakeComparison){
        
        // Loop over Pythia+Hydjet jet pT bins
        for(int iJetPt = firstStudiedJetPtBinEEC[kPythiaHydjet]; iJetPt <= lastStudiedJetPtBinEEC[kPythiaHydjet]; iJetPt++){
          
          // Set the jet pT information for legends and figure saving
          if(iJetPt == histograms[kPythiaHydjet]->GetNJetPtBinsEEC()){
            jetPtString = Form("Jet p_{T} > %.0f", histograms[kPythiaHydjet]->GetCard()->GetJetPtCut());
            compactJetPtString = "";
          } else {
            jetPtString = Form("%.0f < jet p_{T} < %.0f", histograms[kPythiaHydjet]->GetJetPtBinBorderEEC(iJetPt), histograms[kPythiaHydjet]->GetJetPtBinBorderEEC(iJetPt+1));
            compactJetPtString = Form("_J=%.0f-%.0f", histograms[kPythiaHydjet]->GetJetPtBinBorderEEC(iJetPt), histograms[kPythiaHydjet]->GetJetPtBinBorderEEC(iJetPt+1));
          }
          
          // Loop over track pT bins
          for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
            
            trackPtString = Form("%.1f < track p_{T}",histograms[kPythiaHydjet]->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString = Form("_T%.1f",histograms[kPythiaHydjet]->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString.ReplaceAll(".","v");
            
            // Create a legend for the figure
            legend = new TLegend(0.25,jetLegendY1,0.45,jetLegendY1+0.3);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->AddEntry((TObject*) 0, histograms[kPythiaHydjet]->GetCard()->GetAlternativeDataType().Data(), "");
            legend->AddEntry((TObject*) 0, centralityString.Data(),"");
            legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
            
            ptLegend = new TLegend(0.6,minBiasLegendY1,0.8,minBiasLegendY2);
            ptLegend->SetFillStyle(0);ptLegend->SetBorderSize(0);ptLegend->SetTextSize(0.05);ptLegend->SetTextFont(62);
            
            // Create a new canvas for the plot
            drawer->CreateSplitCanvas();
            
            // Logarithmic EEC axis
            if(logEEC) drawer->SetLogY(true);
            
            // For logarithmic drawing, cannot go down to zero
            if(logDeltaR){
              hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
            }
            
            // Set good y-ranges for plotting
            hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0.001, 40);
            
            // First, draw the fake+fake distribution to the upper pad
            hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[0]);
            drawer->DrawHistogramToUpperPad(hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt], "#Deltar", histograms[kPythiaHydjet]->GetEnergyEnergyCorrelatorAxisName(iEnergyEnergyCorrelator), " ");
            legend->AddEntry(hFakeFake[iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt], Form("Pythia+Hydjet: %s",jetPtString.Data()),"l");
            
            ptLegend->AddEntry((TObject*) 0, "Minimum bias Hydjet","");
            
            // Loop over minimum bias Hydjet jet pT bins and draw the distribution to the upper pad
            for(int iJetPtMinBias = firstStudiedJetPtBinEEC[kMinBiasHydjet]; iJetPtMinBias <= lastStudiedJetPtBinEEC[kMinBiasHydjet]; iJetPtMinBias++){
              
              // Set the jet pT information for legends and figure saving
              if(iJetPtMinBias == histograms[kMinBiasHydjet]->GetNJetPtBinsEEC()){
                jetPtStringMinBias = Form("Jet p_{T} > %.0f", histograms[kMinBiasHydjet]->GetCard()->GetJetPtCut());
                compactJetPtStringMinBias = "";
              } else {
                jetPtStringMinBias = Form("%.0f < jet p_{T} < %.0f", histograms[kMinBiasHydjet]->GetJetPtBinBorderEEC(iJetPtMinBias), histograms[kMinBiasHydjet]->GetJetPtBinBorderEEC(iJetPtMinBias+1));
                compactJetPtStringMinBias = Form("_J=%.0f-%.0f", histograms[kMinBiasHydjet]->GetJetPtBinBorderEEC(iJetPtMinBias), histograms[kMinBiasHydjet]->GetJetPtBinBorderEEC(iJetPtMinBias+1));
              }
              
              hMinBias[iEnergyEnergyCorrelator][iCentrality][iJetPtMinBias][iTrackPt]->SetLineColor(color[iJetPtMinBias+1]);
              hMinBias[iEnergyEnergyCorrelator][iCentrality][iJetPtMinBias][iTrackPt]->Draw("same");
              ptLegend->AddEntry(hMinBias[iEnergyEnergyCorrelator][iCentrality][iJetPtMinBias][iTrackPt], jetPtStringMinBias.Data(), "l");
              
            } // Minimum bias Hydjet jet pT bins
            
            // Draw the legends to the upper pad
            legend->Draw();
            ptLegend->Draw();
            
            // Linear scale for the ratio
            drawer->SetLogY(false);
            
            // Loop over minimum bias Hydjet jet pT bins and draw the ratios to the lower pad
            for(int iJetPtMinBias = firstStudiedJetPtBinEEC[kMinBiasHydjet]; iJetPtMinBias <= lastStudiedJetPtBinEEC[kMinBiasHydjet]; iJetPtMinBias++){
              
              // For logarithmic drawing, cannot go down to zero
              if(logDeltaR){
                hMinBiasToFakeFakeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
              }
              
              hMinBiasToFakeFakeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt]->SetLineColor(color[iJetPtMinBias+1]);
              if(iJetPtMinBias == firstStudiedJetPtBinEEC[kMinBiasHydjet]){
                drawer->SetGridY(true);
                hMinBiasToFakeFakeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
                drawer->DrawHistogramToLowerPad(hMinBiasToFakeFakeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt], "#Deltar", "#frac{Color}{Fake+fake}", " ");
                drawer->SetGridY(false);
              } else {
                hMinBiasToFakeFakeRatio[iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt]->Draw("same");
              }
              
            } // Minimum bias Hydjet jet pT bins
            
          } // Track pT bin loop
        } // Pythia+Hydjet jet pT bins
      } // Drawing minimum bias to fake+fake comparison
    } // Centrality loop
  } // Energy-energy correlator loop


}
