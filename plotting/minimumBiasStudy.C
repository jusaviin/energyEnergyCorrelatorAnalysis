#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Macro studying specifically the fake+fake background in.
 * It is assumed in the code that the centrality and track pT bins match between the files
 */
void minimumBiasStudy(){

  
  // File containing the Pythia+Hydjet simulation result (index 0), and the one containing minimum bias Hydjet result (index 1)
  const int nRegularFiles = 2;
  TString inputFileName[] = {"data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJet_updatedMultiplicityAndDensity_wtaAxis_noTrigger_preprocessed_2022-10-17.root", "data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_onlyFakeJets_wtaAxis_noTrigger_preprocessed_2022-11-22.root"};
  // data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJet_fakeFakeReflectedCone_noTrigger_preprocessed_2022-09-23.root
  // data/PbPbMC2018_GenGen_eecAnalysis_genJet_fakeFakeReflectedCone_noTrigger_preprocessed_2022-09-30.root
  // "data/MinBiasHydjet_RecoGen_eecAnalysis_akFlowJet_firstMinBiasScan_noTrigger_preprocessed_2022-10-10.root"
  TString minBiasFileName = "data/MinBiasHydjet_RecoGen_eecAnalysis_akFlowJet_MnD_wtaAxis_noTrigger_preprocessed_2022-10-19.root";
  
  TString legendString[] = {"Pythia+Hydjet", "Only fake"};
  
  // Open the input files and read the cards
  TFile* inputFile[nRegularFiles+1];
  EECCard* card[nRegularFiles+1];
  for(int iFile = 0; iFile <= nRegularFiles; iFile++){
    if(iFile < nRegularFiles){
      inputFile[iFile] = TFile::Open(inputFileName[iFile]);
    } else {
      inputFile[iFile] = TFile::Open(minBiasFileName);
    }
    
    // Check that the files exist
    if(inputFile[iFile] == NULL){
      if(iFile < nRegularFiles){
        cout << "Error! The file " << inputFileName[iFile].Data() << " does not exist!" << endl;
      } else {
        cout << "Error! The file " << minBiasFileName.Data() << " does not exist!" << endl;
      }
      cout << "Maybe you forgot the data/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }
    
    // Load the card from the file and read the collision system
    card[iFile] = new EECCard(inputFile[iFile]);
    
  }
    
  // Require matching centrality and track pT bins
  const double epsilon = 0.00001;
  
  const int nCentralityBins = card[0]->GetNCentralityBins();
  if(nCentralityBins != card[nRegularFiles]->GetNCentralityBins()){
    cout << "Error! Centrality bins do not match between the files! Cannot execute the code!" << endl;
    return;
  }
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    if(TMath::Abs(card[0]->GetLowBinBorderCentrality(iCentrality) - card[nRegularFiles]->GetLowBinBorderCentrality(iCentrality)) > epsilon){
      cout << "Error! Centrality bins do not match between the files! Cannot execute the code!" << endl;
      return;
    }
    if(TMath::Abs(card[0]->GetHighBinBorderCentrality(iCentrality) - card[nRegularFiles]->GetHighBinBorderCentrality(iCentrality)) > epsilon){
      cout << "Error! Centrality bins do not match between the files! Cannot execute the code!" << endl;
      return;
    }
  }
  
  const int nTrackPtBinsEEC = card[0]->GetNTrackPtBinsEEC();
  if(nTrackPtBinsEEC != card[nRegularFiles]->GetNTrackPtBinsEEC()){
    cout << "Error! Track pT bins do not match between the files! Cannot execute the code!" << endl;
    return;
  }
  for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
    if(TMath::Abs(card[0]->GetLowBinBorderTrackPtEEC(iTrackPt) - card[nRegularFiles]->GetLowBinBorderTrackPtEEC(iTrackPt)) > epsilon){
      cout << "Error! Track pT bins do not match between the files! Cannot execute the code!" << endl;
      return;
    }
    if(TMath::Abs(card[0]->GetHighBinBorderTrackPtEEC(iTrackPt) - card[nRegularFiles]->GetHighBinBorderTrackPtEEC(iTrackPt)) > epsilon){
      cout << "Error! Track pT bins do not match between the files! Cannot execute the code!" << endl;
      return;
    }
  }
  
  // ====================================================
  //    Binning configuration for the background study
  // ====================================================
  
  // Find the number of jet pT bins
  const int nJetPtBinsEEC[2] = {card[0]->GetNJetPtBinsEEC(), card[nRegularFiles]->GetNJetPtBinsEEC()};
  
  // Default binning ranges for reference
  // centrality = {0,10,30,50,90}
  // track pT =  {0.7,1,2,3,4,8,12,300}
  // jet pT in energy-energy correlators = {120,140,160,180,200,300,500,5030}
  // track pT in energy-energy correlators = {0.7,1,2,3,4,6,300}
  
  // Bin range to be integrated
  int firstStudiedCentralityBin = 0;
  int lastStudiedCentralityBin = 0;
  
  int firstStudiedJetPtBinEEC[2] = {nJetPtBinsEEC[0],nJetPtBinsEEC[1]-1};
  int lastStudiedJetPtBinEEC[2] = {nJetPtBinsEEC[0],nJetPtBinsEEC[1]-1}; // Note: Jets integrated over all pT ranges are in nJetPtBinsEEC bin
  
  int firstStudiedTrackPtBinEEC = 0;
  int lastStudiedTrackPtBinEEC = 5;
  
  // Select the types of energy-energy correlators are studied
  bool studyEnergyEnergyCorrelator[EECHistogramManager::knEnergyEnergyCorrelatorTypes];
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator] = true;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPt] = false;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorUncorrected] = false;
  studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPtUncorrected] = false;
  
  // Logarithmic axes
  const bool logDeltaR = true;
  const bool logEEC = true;
  
  // Axis zooming
  std::pair<double,double> ratioZoom = std::make_pair(0, 2);
  
  // Figure saving
  const bool saveFigures = true;  // Save figures
  const char* saveComment = "";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures
  
  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager *histograms[nRegularFiles+1];
  for(int iFile = 0; iFile <= nRegularFiles; iFile++){
    histograms[iFile] = new EECHistogramManager(inputFile[iFile],card[iFile]);
    
    // Choose the energy-energy correlator types to load
    histograms[iFile]->SetLoadEnergyEnergyCorrelators(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelator]);
    histograms[iFile]->SetLoadEnergyEnergyCorrelatorsJetPt(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPt]);
    histograms[iFile]->SetLoadEnergyEnergyCorrelatorsUncorrected(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorUncorrected]);
    histograms[iFile]->SetLoadEnergyEnergyCorrelatorsJetPtUncorrected(studyEnergyEnergyCorrelator[EECHistogramManager::kEnergyEnergyCorrelatorJetPtUncorrected]);
    
    // Choose the bin ranges
    histograms[iFile]->SetCentralityBinRange(firstStudiedCentralityBin,lastStudiedCentralityBin);
    histograms[iFile]->SetJetPtBinRangeEEC(firstStudiedJetPtBinEEC[iFile/nRegularFiles],lastStudiedJetPtBinEEC[iFile/nRegularFiles]);
    histograms[iFile]->SetTrackPtBinRangeEEC(firstStudiedTrackPtBinEEC,lastStudiedTrackPtBinEEC);
    
    // Load the histograms from the file
    histograms[iFile]->LoadProcessedHistograms();
  }
  
  // Energy-energy correlator histograms separated by subevents from the Pythia+Hydjet simulation
  TH1D* hRegular[nRegularFiles][EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC[0]+1][nTrackPtBinsEEC];
  
  // Reflected cone energy-energy correlators
  TH1D* hMinBias[EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC[1]+1][nTrackPtBinsEEC];
  
  // Histograms for all different ratios
  TH1D* hMinBiasToRegularRatio[nRegularFiles][EECHistogramManager::knEnergyEnergyCorrelatorTypes][nCentralityBins][nJetPtBinsEEC[0]+1][nJetPtBinsEEC[1]+1][nTrackPtBinsEEC];
  
  // Initialize the energy-energy correlator histogram array to NULL
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        for(int iJetPtRegular = 0; iJetPtRegular < nJetPtBinsEEC[0]+1; iJetPtRegular++){
          for(int iFile = 0; iFile < nRegularFiles; iFile++){
            hRegular[iFile][iEnergyEnergyCorrelator][iCentrality][iJetPtRegular][iTrackPt] = NULL;
            for(int iJetPtMinBias = 0; iJetPtMinBias < nJetPtBinsEEC[1]+1; iJetPtMinBias++){
              hMinBiasToRegularRatio[iFile][iEnergyEnergyCorrelator][iCentrality][iJetPtRegular][iJetPtMinBias][iTrackPt] = NULL;
            }
          }
        } // Pythia+Hydjet jet pT loop
        for(int iJetPtMinBias = 0; iJetPtMinBias < nJetPtBinsEEC[1]+1; iJetPtMinBias++){
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
    
    // Read the regular histograms
    for(int iFile = 0; iFile < nRegularFiles; iFile++){
      for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
          for(int iJetPtRegular = firstStudiedJetPtBinEEC[0]; iJetPtRegular <= lastStudiedJetPtBinEEC[0]; iJetPtRegular++){
            
            // Read the fake+fake histogram and normalize everything to one
            hRegular[iFile][iEnergyEnergyCorrelator][iCentrality][iJetPtRegular][iTrackPt] = histograms[iFile]->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPtRegular, iTrackPt, EECHistograms::kSameJetPair);
            hRegular[iFile][iEnergyEnergyCorrelator][iCentrality][iJetPtRegular][iTrackPt]->Scale(1/hRegular[iFile][iEnergyEnergyCorrelator][iCentrality][iJetPtRegular][iTrackPt]->Integral("width"));
            
          } // Regular jet pT loop
        } // Track pT loop
      } // Centrality loop
    } // File loop
    
    // Read the minimum bias histograms
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
        for(int iJetPtMinBias = firstStudiedJetPtBinEEC[1]; iJetPtMinBias <= lastStudiedJetPtBinEEC[1]; iJetPtMinBias++){
          
          // Read the minimum bias histogram and normalize everything to one
          hMinBias[iEnergyEnergyCorrelator][iCentrality][iJetPtMinBias][iTrackPt] = histograms[nRegularFiles]->GetHistogramEnergyEnergyCorrelator(iEnergyEnergyCorrelator, iCentrality, iJetPtMinBias, iTrackPt, EECHistograms::kSameJetPair);
          hMinBias[iEnergyEnergyCorrelator][iCentrality][iJetPtMinBias][iTrackPt]->Scale(1/hMinBias[iEnergyEnergyCorrelator][iCentrality][iJetPtMinBias][iTrackPt]->Integral("width"));
          
          // Calculate the min bias to fake+fake ratio
          for(int iFile = 0; iFile < nRegularFiles; iFile++){
            for(int iJetPtRegular = firstStudiedJetPtBinEEC[0]; iJetPtRegular <= lastStudiedJetPtBinEEC[0]; iJetPtRegular++){
              hMinBiasToRegularRatio[iFile][iEnergyEnergyCorrelator][iCentrality][iJetPtRegular][iJetPtMinBias][iTrackPt] = (TH1D*) hRegular[iFile][iEnergyEnergyCorrelator][iCentrality][iJetPtRegular][iTrackPt]->Clone(Form("minBiasToFakeFakeRatio%d%d%d%d%d%d", iFile, iEnergyEnergyCorrelator, iCentrality, iJetPtRegular, iJetPtMinBias, iTrackPt));
              hMinBiasToRegularRatio[iFile][iEnergyEnergyCorrelator][iCentrality][iJetPtRegular][iJetPtMinBias][iTrackPt]->Divide(hMinBias[iEnergyEnergyCorrelator][iCentrality][iJetPtMinBias][iTrackPt]);
            } // Regular jet pT loop
          } // Regular file loop
          
        } // Min bias jet pT loop
      } // Track pT loop
    } // Centrality loop
    
  } // Energy-energy correlator loop
  
  // Find the binning information for DeltaR bins from a histogram
  const int nDeltaRBins = hRegular[0][studiedEnergyEnergyCorrelatorIndex][firstStudiedCentralityBin][firstStudiedJetPtBinEEC[0]][firstStudiedTrackPtBinEEC]->GetNbinsX();
  double deltaRBinBorders[nDeltaRBins+1];
  deltaRBinBorders[0] = hRegular[0][studiedEnergyEnergyCorrelatorIndex][firstStudiedCentralityBin][firstStudiedJetPtBinEEC[0]][firstStudiedTrackPtBinEEC]->GetXaxis()->GetBinLowEdge(1);
  for(int iBin = 1; iBin <= nDeltaRBins; iBin++){
    deltaRBinBorders[iBin] = hRegular[0][studiedEnergyEnergyCorrelatorIndex][firstStudiedCentralityBin][firstStudiedJetPtBinEEC[0]][firstStudiedTrackPtBinEEC]->GetXaxis()->GetBinUpEdge(iBin);
  }
  
  
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);

  if(logDeltaR) drawer->SetLogX(true);

  TLegend *legend;
  TString centralityString, trackPtString, jetPtString, jetPtStringMinBias;
  TString compactCentralityString, compactTrackPtString, compactJetPtString, compactJetPtStringMinBias;
  int color[10] = {kBlack, kBlue, kRed, kGreen+2, kCyan, kMagenta, kOrange-1, kAzure+5, kOrange-2, kGray};

  // Loop over all selected histograms
  for(int iEnergyEnergyCorrelator = 0; iEnergyEnergyCorrelator < EECHistogramManager::knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelator++){
    
    // Only read the selected energy-energy correlator types
    if(!studyEnergyEnergyCorrelator[iEnergyEnergyCorrelator]) continue;
    
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      
      // Set the centrality information for legends and figure saving
      centralityString = Form("Cent: %.0f-%.0f%%", histograms[0]->GetCentralityBinBorder(iCentrality), histograms[0]->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f", histograms[0]->GetCentralityBinBorder(iCentrality), histograms[0]->GetCentralityBinBorder(iCentrality+1));
      
      // Loop over Pythia+Hydjet jet pT bins
      for(int iJetPt = firstStudiedJetPtBinEEC[0]; iJetPt <= lastStudiedJetPtBinEEC[0]; iJetPt++){
        
        // Set the jet pT information for legends and figure saving
        if(iJetPt == histograms[0]->GetNJetPtBinsEEC()){
          jetPtString = Form("Jet p_{T} > %.0f", histograms[0]->GetCard()->GetJetPtCut());
          compactJetPtString = "";
        } else {
          jetPtString = Form("%.0f < jet p_{T} < %.0f", histograms[0]->GetJetPtBinBorderEEC(iJetPt), histograms[0]->GetJetPtBinBorderEEC(iJetPt+1));
          compactJetPtString = Form("_J=%.0f-%.0f", histograms[0]->GetJetPtBinBorderEEC(iJetPt), histograms[0]->GetJetPtBinBorderEEC(iJetPt+1));
        }
        
        // Loop over track pT bins
        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
          
          trackPtString = Form("%.1f < track p_{T}",histograms[0]->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString = Form("_T%.1f",histograms[0]->GetTrackPtBinBorderEEC(iTrackPt));
          compactTrackPtString.ReplaceAll(".","v");
          
          // Create a legend for the figure
          legend = new TLegend(0.45,0.08,0.65,0.4);
          legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
          legend->AddEntry((TObject*) 0, centralityString.Data(),"");
          legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
          
          // Loop over minimum bias Hydjet jet pT bins and draw the distribution to the upper pad
          for(int iJetPtMinBias = firstStudiedJetPtBinEEC[1]; iJetPtMinBias <= lastStudiedJetPtBinEEC[1]; iJetPtMinBias++){
            
            // Create a new canvas for the plot
            drawer->CreateSplitCanvas();
            
            // Set the jet pT information for legends and figure saving
            if(iJetPtMinBias == histograms[nRegularFiles]->GetNJetPtBinsEEC()){
              jetPtStringMinBias = Form("Jet p_{T} > %.0f", histograms[nRegularFiles]->GetCard()->GetJetPtCut());
              compactJetPtStringMinBias = "";
            } else {
              jetPtStringMinBias = Form("%.0f < jet p_{T} < %.0f", histograms[nRegularFiles]->GetJetPtBinBorderEEC(iJetPtMinBias), histograms[nRegularFiles]->GetJetPtBinBorderEEC(iJetPtMinBias+1));
              compactJetPtStringMinBias = Form("_J=%.0f-%.0f", histograms[nRegularFiles]->GetJetPtBinBorderEEC(iJetPtMinBias), histograms[nRegularFiles]->GetJetPtBinBorderEEC(iJetPtMinBias+1));
            }
            
            // Logarithmic EEC axis
            if(logEEC) drawer->SetLogY(true);
            
            // For logarithmic drawing, cannot go down to zero
            if(logDeltaR){
              hMinBias[iEnergyEnergyCorrelator][iCentrality][iJetPtMinBias][iTrackPt]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
            }
            
            // Set good y-ranges for plotting
            hMinBias[iEnergyEnergyCorrelator][iCentrality][iJetPtMinBias][iTrackPt]->GetYaxis()->SetRangeUser(0.001, 40);
            
            // Draw the minimum bias distribution
            hMinBias[iEnergyEnergyCorrelator][iCentrality][iJetPtMinBias][iTrackPt]->SetLineColor(color[0]);
            drawer->DrawHistogramToUpperPad(hMinBias[iEnergyEnergyCorrelator][iCentrality][iJetPtMinBias][iTrackPt], "#Deltar", histograms[nRegularFiles]->GetEnergyEnergyCorrelatorAxisName(iEnergyEnergyCorrelator), " ");
            legend->AddEntry(hMinBias[iEnergyEnergyCorrelator][iCentrality][iJetPtMinBias][iTrackPt], Form("MinBias Hydjet, %s", jetPtStringMinBias.Data()), "l");
            
            
            // Draw all other distributions
            for(int iFile = 0; iFile < nRegularFiles; iFile++){
              
              // First, draw the fake+fake distribution to the upper pad
              hRegular[iFile][iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iFile+1]);
              hRegular[iFile][iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt]->Draw("same");
              legend->AddEntry(hRegular[iFile][iEnergyEnergyCorrelator][iCentrality][iJetPt][iTrackPt], Form("%s: %s", legendString[iFile].Data(), jetPtString.Data()),"l");
            } // File loop
                        
            // Draw the legends to the upper pad
            legend->Draw();
            
            // Linear scale for the ratio
            drawer->SetLogY(false);
            
            // Draw all the ratios
            for(int iFile = 0; iFile < nRegularFiles; iFile++){
            
              // For logarithmic drawing, cannot go down to zero
              if(logDeltaR){
                hMinBiasToRegularRatio[iFile][iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt]->GetXaxis()->SetRangeUser(0.001, deltaRBinBorders[nDeltaRBins]);
              }
              
              hMinBiasToRegularRatio[iFile][iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt]->SetLineColor(color[iFile+1]);
              if(iFile == 0){
                drawer->SetGridY(true);
                hMinBiasToRegularRatio[iFile][iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
                drawer->DrawHistogramToLowerPad(hMinBiasToRegularRatio[iFile][iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt], "#Deltar", "#frac{Color}{MinBias}", " ");
                drawer->SetGridY(false);
              } else {
                hMinBiasToRegularRatio[iFile][iEnergyEnergyCorrelator][iCentrality][iJetPt][iJetPtMinBias][iTrackPt]->Draw("same");
              }
              
            } // Regular file loop
            
            // Save the figures to a file
            if(saveFigures){
              gPad->GetCanvas()->SaveAs(Form("figures/regularToMinBiasComparison%s_%s%s%s%s%s.%s", saveComment, histograms[0]->GetEnergyEnergyCorrelatorHistogramName(iEnergyEnergyCorrelator), compactCentralityString.Data(), compactJetPtString.Data(), compactJetPtStringMinBias.Data(), compactTrackPtString.Data(), figureFormat));
            }
            
          } // Minimum bias Hydjet jet pT bins          
        } // Track pT bin loop
      } // Pythia+Hydjet jet pT bins
    } // Centrality loop
  } // Energy-energy correlator loop


}
