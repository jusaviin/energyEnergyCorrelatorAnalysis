#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "../src/EECHistograms.h"
#include "EECCard.h"
#include "JDrawer.h"
#include <tuple>

/*
 * Fit a Gauss function to a histogram and extract parameters from that
 *
 *  TH1* histogram = Histogram from which drawing range in searched
 *
 *  return: Gauss mean, Gauss sigma, Error for Gauss mean, Error for Gauss sigma
 */
std::tuple<double,double,double,double> fitGauss(TH1* histogram, TString title = "", TString jetTypeString = "", TString centralityBin = "", TString ptBin = "",  TString saveName = ""){
  histogram->Fit("gaus","","",0.5,1.5);
  TF1 * gaussFit = histogram->GetFunction("gaus");
  
  double gaussMean = 0;
  double gaussSigma = 0;
  double gaussMeanError = 0;
  double gaussSigmaError = 0;
  
  if(gaussFit){
    gaussMean = gaussFit->GetParameter(1);
    gaussSigma = gaussFit->GetParameter(2);
    gaussMeanError = gaussFit->GetParError(1);
    gaussSigmaError = gaussFit->GetParError(2);
  }
  
  // If title is given, print the fit
  if(!title.EqualTo("")){
    JDrawer *temporaryDrawer = new JDrawer();
    temporaryDrawer->DrawHistogram(histogram,"Reco p_{T} / Gen p_{T}","Counts", " ");
    TLegend *legend = new TLegend(0.57,0.68,0.8,0.93);
    legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
    legend->SetHeader(title);
    legend->AddEntry((TObject*)0, jetTypeString, "");
    legend->AddEntry((TObject*)0, centralityBin, "");
    legend->AddEntry((TObject*)0, ptBin, "");
    legend->Draw();
    
    if(!saveName.EqualTo("")){
      gPad->GetCanvas()->SaveAs(saveName);
    }
    
  }
  
  return std::make_tuple(gaussMean,gaussSigma,gaussMeanError,gaussSigmaError);
}

/*
 * Draw a closure histogram with quark/gluon discrimination to the plot
 *
 *  TH1D *histogram[EECHistograms::knClosureTypes+1] = Array of histograms for each closure particle type
 *  const char* xTitle = Title given to the x-axis
 *  const char* yTitle = Title given to the y-axis
 *  bool ppData = Are we using pp od PbPb data in closure plotting
 *  int iCentrality = Index of the centrality bin
 *  int legendNzoom = Define the y-axis zoom and the legend position in the plot
 *  bool includeQuarkGluon = Flag for including quark and gluon jet only closures
 *  const char* saveComment = Comment given to the save name file
 *  bool saveFigures = Choose whether to save the figures or not
 */
void drawClosureHistogram(TH1D *histogram[EECHistograms::knClosureParticleTypes+1], const char* xTitle, const char* yTitle, bool ppData, int iCentrality, int legendNzoom, bool includeQuarkGluon, const char* saveComment, bool saveFigures){
  
  // Create a new drawer and define bin borders and drawing style
  JDrawer *drawer = new JDrawer();
  drawer->SetCanvasSize(700,700);
  double centralityBinBorders[5] = {0,10,30,50,90};
  const char* centralityString;
  const char* centralitySaveName;
  const char* namer;
  int lineColors[2] = {kBlue,kRed};
  const char* particleNames[2] = {"Quark jets","Gluon jets"};
  
  // Zooming and legend position options
  double yZoomLow = 0.9;
  double yZoomHigh = 1.1;
  double legendX1 = 0.51;
  double legendX2 = 0.93;
  double legendY1 = 0.65;
  double legendY2 = 0.92;
  
  if(legendNzoom == 1){
    yZoomLow = 0;
    yZoomHigh = 0.34; // 0.24
    legendX1 = 0.51;
    legendX2 = 0.93;
    legendY1 = 0.65;
    legendY2 = 0.92;
    //legendX1 = 0.46;
    //legendX2 = 0.88;
    //legendY1 = 0.65;
    //legendY2 = 0.92;
  } else if(legendNzoom == 2){
    legendX1 = 0.16;
    legendX2 = 0.58;
    legendY1 = 0.78;
    legendY2 = 0.99;
  }
  
  // Set a good style for the inclusive histogram
  histogram[EECHistograms::knClosureParticleTypes]->SetLineColor(kBlack);
  histogram[EECHistograms::knClosureParticleTypes]->GetYaxis()->SetRangeUser(yZoomLow,yZoomHigh);
  drawer->DrawHistogram(histogram[EECHistograms::knClosureParticleTypes],xTitle,yTitle," ");
  
  // Create a legend to the plot
  TLegend *legend = new TLegend(legendX1,legendY1,legendX2,legendY2);
  legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.04);legend->SetTextFont(62);
  
  if(ppData){
    centralityString = ", pp";
    centralitySaveName = "_pp";
  } else {
    centralityString = Form(", Cent:%.0f-%.0f", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]);
    centralitySaveName = Form("_C=%.0f-%.0f", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]);
  }
  
  namer = Form("Inclusive jet%s", centralityString);
  legend->SetHeader(namer);
  legend->AddEntry(histogram[EECHistograms::knClosureParticleTypes],"All jets","l");
  
  if(includeQuarkGluon){
    for(int iClosureParticle = 0; iClosureParticle < EECHistograms::knClosureParticleTypes; iClosureParticle++){
      histogram[iClosureParticle]->SetLineColor(lineColors[iClosureParticle]);
      histogram[iClosureParticle]->Draw("same");
      legend->AddEntry(histogram[iClosureParticle],particleNames[iClosureParticle],"l");
    } // Closure particle loop (quark/gluon)
  }
  
  legend->Draw();
  
  // Draw a dashed line at one
  double startPoint = histogram[EECHistograms::knClosureParticleTypes]->GetXaxis()->GetBinLowEdge(1);
  double endPoint = histogram[EECHistograms::knClosureParticleTypes]->GetXaxis()->GetBinUpEdge(histogram[EECHistograms::knClosureParticleTypes]->GetNbinsX());
  TLine *oneLine = new TLine(startPoint, 1, endPoint, 1);
  oneLine->SetLineStyle(2);
  oneLine->SetLineColor(kBlack);
  oneLine->Draw("same");
  
  // Save the figures if selected to do so
  if(saveFigures){
    namer = Form("figures/jet%sInclusiveJet%s.pdf",saveComment,centralitySaveName);
    gPad->GetCanvas()->SaveAs(namer);
  }
  
}

/*
 * Macro for producing the JFF correction from PYTHIA simulation results
 */
void constructJetPtClosures(){
  
  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  TString closureFileName = "data/PbPbMC2018_RecoGen_akFlowJets_miniAOD_4pCentShift_noTrigger_jetPtClosure_newRatioBins_processed_2023-04-21.root";
  // data/PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_jetPtClosure_finalMcWeight_fixCentrality_processed_2023-03-06.root
  // data/PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_jetPtClosure_finalMcWeight_processed_2023-03-06.root
  // data/PbPbMC2018_GenGen_eecAnalysis_akFlowJets_mAOD_4pC_wtaAxis_noTrig_matchJetPt_closures_processed_2023-02-07.root
  // data/PbPbMC2018_GenGen_eecAnalysis_akFlowJets_mAODnewR_4pC_wtaAxis_noTrigger_jetPtClosure_processed_2023-01-30.root
  // data/ppMC2017_GenGen_Pythia8_pfJets_wtaAxis_noCorrelations_jetPtClosures_processed_2023-01-13.root
  
  bool drawPtClosure = true;
  bool drawEtaClosure = false;
  
  bool includeQuarkGluon = false; // Include only quark and only gluon jet curves
  bool drawGaussFitsPt = false;
    
  bool fitResolution = true;  // Fit the jet pT resolution histograms
  
  bool saveFigures = false;  // Save the figures to file
  
  // ==================================================================
  // =================== Configuration ready ==========================
  // ==================================================================
  
  // Open the input files
  TFile *closureFile = TFile::Open(closureFileName);
  
  // Check if we are using pp or PbPb data
  EECCard *card = new EECCard(closureFile);
  TString collisionSystem = card->GetDataType();
  bool ppData = (collisionSystem.Contains("pp") || collisionSystem.Contains("localTest"));
  
  // Create histogram managers to provide the histograms for the correction
  EECHistogramManager *closureHistograms = new EECHistogramManager(closureFile);
  closureHistograms->SetLoadJetPtClosureHistograms(true);
  if(ppData) closureHistograms->SetCentralityBinRange(0,0);  // Disable centrality binning for pp data
  closureHistograms->LoadProcessedHistograms();
  
  // Find the correct number of centrality and track pT bins
  const int nCentralityBins = ppData ? 1 : closureHistograms->GetNCentralityBins();
  double centralityBinBorders[5] = {0,10,30,50,90};  // Bin borders for centrality
  
  // Initialize reco/gen ratio and closure histograms
  TH1D *hRecoGenRatio[EECHistogramManager::knGenJetPtBins][nCentralityBins][EECHistograms::knClosureParticleTypes+1];
  TH1D *hRecoGenRatioEta[EECHistogramManager::knJetEtaBins][nCentralityBins][EECHistograms::knClosureParticleTypes+1];
  TH1D *hJetPtClosure[nCentralityBins][EECHistograms::knClosureParticleTypes+1];
  TH1D *hJetPtClosureSigma[nCentralityBins][EECHistograms::knClosureParticleTypes+1];
  TH1D *hJetPtClosureEta[nCentralityBins][EECHistograms::knClosureParticleTypes+1];
  TH1D *hJetPtClosureSigmaEta[nCentralityBins][EECHistograms::knClosureParticleTypes+1];
  const char* histogramNamer;
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iClosureParticle = 0; iClosureParticle < EECHistograms::knClosureParticleTypes+1; iClosureParticle++){
      for(int iGenJetPt = 0; iGenJetPt < EECHistogramManager::knGenJetPtBins; iGenJetPt++){
        hRecoGenRatio[iGenJetPt][iCentrality][iClosureParticle] = NULL;
      } // Generator level jet pT loop
      for(int iJetEta = 0; iJetEta < EECHistogramManager::knJetEtaBins; iJetEta++){
        hRecoGenRatioEta[iJetEta][iCentrality][iClosureParticle] = NULL;
      }
      histogramNamer = Form("jetPtClosure_Cent%d_Part%d", iCentrality, iClosureParticle);
      hJetPtClosure[iCentrality][iClosureParticle] = new TH1D(histogramNamer,histogramNamer,45,50,500);
      histogramNamer = Form("jetPtClosureSigma_Cent%d_Part%d", iCentrality, iClosureParticle);
      hJetPtClosureSigma[iCentrality][iClosureParticle] = new TH1D(histogramNamer,histogramNamer,45,50,500);
      histogramNamer = Form("jetPtClosureEta_Cent%d_Part%d", iCentrality, iClosureParticle);
      hJetPtClosureEta[iCentrality][iClosureParticle] = new TH1D(histogramNamer,histogramNamer,50,-2.5,2.5);
      histogramNamer = Form("jetPtClosureSigmaEta_Cent%d_Part%d", iCentrality, iClosureParticle);
      hJetPtClosureSigmaEta[iCentrality][iClosureParticle] = new TH1D(histogramNamer,histogramNamer,50,-2.5,2.5);
    } // Closure particle loop (quark/gluon/no selection)
  } // Centrality loop
  
  JDrawer *drawer = new JDrawer();
  TF1* gaussFit;
  double gaussMean = 0;
  double gaussSigma = 0;
  double gaussMeanError = 0;
  double gaussSigmaError = 0;
  int minGenPt = 7;  // Set this to 7 to skip bins below 120 GeV
  TString genPtString;
  TString centralityString;
  TString jetTypeName[EECHistograms::knClosureParticleTypes+1] = {"Quark", "Gluon", "All"};
  TString jetTypeString;
  TString gaussFitSaveString;
  
  // Read the reco/gen histograms from the file and fit them to construct the closure plots
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    for(int iClosureParticle = 0; iClosureParticle < EECHistograms::knClosureParticleTypes+1; iClosureParticle++){
      if(iClosureParticle < EECHistograms::knClosureParticleTypes && !includeQuarkGluon) continue; // Only create the quark and gluon curves if selected
      if(drawPtClosure){
        for(int iGenJetPt = minGenPt; iGenJetPt < EECHistogramManager::knGenJetPtBins; iGenJetPt++){

          // Read the reco/gen histogram from the file
          hRecoGenRatio[iGenJetPt][iCentrality][iClosureParticle] = closureHistograms->GetHistogramJetPtClosure(iGenJetPt, EECHistogramManager::knJetEtaBins, iCentrality, iClosureParticle);
          
          // Fit a gauss to the histogram
          if(drawGaussFitsPt){
            genPtString = Form("%d < Gen p_{T} < %d", 50+10*iGenJetPt, 60+10*iGenJetPt);
            centralityString = Form("Cent: %.0f-%.0f", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]);
            jetTypeString = Form("%s jets", jetTypeName[iClosureParticle].Data());
            gaussFitSaveString = Form("figures/jetPtClosureGaussFit%sJets_T%dC%d.pdf", jetTypeName[iClosureParticle].Data(), iGenJetPt, iCentrality);
            std::tie(gaussMean,gaussSigma,gaussMeanError,gaussSigmaError) = fitGauss(hRecoGenRatio[iGenJetPt][iCentrality][iClosureParticle], " ", jetTypeString, centralityString, genPtString, gaussFitSaveString);
          } else {
            std::tie(gaussMean,gaussSigma,gaussMeanError,gaussSigmaError) = fitGauss(hRecoGenRatio[iGenJetPt][iCentrality][iClosureParticle]);
          }
          
          // Fill the histogram with the fit parameters
          hJetPtClosure[iCentrality][iClosureParticle]->SetBinContent(iGenJetPt+1,gaussMean);
          hJetPtClosure[iCentrality][iClosureParticle]->SetBinError(iGenJetPt+1,gaussMeanError);
          hJetPtClosureSigma[iCentrality][iClosureParticle]->SetBinContent(iGenJetPt+1,gaussSigma);
          hJetPtClosureSigma[iCentrality][iClosureParticle]->SetBinError(iGenJetPt+1,gaussSigmaError);
          
        } // Generator level jet pT loop
      } // pT closure if
      
      if(drawEtaClosure){
        // For eta, bins from 9 to nBins-9 cover the region -1.6 < eta < 1.6
        for(int iJetEta = 9; iJetEta < EECHistogramManager::knJetEtaBins-9; iJetEta++){
          
          // Read the reco/gen histogram from the file
          hRecoGenRatioEta[iJetEta][iCentrality][iClosureParticle] = closureHistograms->GetHistogramJetPtClosure(EECHistogramManager::knGenJetPtBins, iJetEta, iCentrality, iClosureParticle);
          
          // Fit a gauss to the histogram
          std::tie(gaussMean,gaussSigma,gaussMeanError,gaussSigmaError) = fitGauss(hRecoGenRatioEta[iJetEta][iCentrality][iClosureParticle]);
          
          // Fill the histogram with the fit parameters
          hJetPtClosureEta[iCentrality][iClosureParticle]->SetBinContent(iJetEta+1,gaussMean);
          hJetPtClosureEta[iCentrality][iClosureParticle]->SetBinError(iJetEta+1,gaussMeanError);
          hJetPtClosureSigmaEta[iCentrality][iClosureParticle]->SetBinContent(iJetEta+1,gaussSigma);
          hJetPtClosureSigmaEta[iCentrality][iClosureParticle]->SetBinError(iJetEta+1,gaussSigmaError);
          
        } // Jet eta loop
      } // eta closure if
      
    } // Closure particle loop (quark/gluon/no selection)
  } // Centrality loop
  
  double minFitPt = 50+10*minGenPt;
  double maxFitPt = 500;
  
  // Fit the resolution plots with a polynomial function
  if(fitResolution){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      cout << "Fitting iCentrality: " << iCentrality << endl;
      if(iCentrality == 2) hJetPtClosureSigma[iCentrality][EECHistograms::knClosureParticleTypes]->SetBinError(16,100);
      hJetPtClosureSigma[iCentrality][EECHistograms::knClosureParticleTypes]->Fit("pol4","","",minFitPt,maxFitPt);
    } // Centrality loop
  } // Fitting the resolution
  
  // Draw the closure plots
  int lineColors[2] = {kBlue,kRed};
  const char* particleNames[2] = {"Quark","Gluon"};
  TLegend *legend;
  
  // Create the output file
  TFile *outputFile = new TFile("closureThingy.root","UPDATE");
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    if(drawPtClosure){
      drawClosureHistogram(hJetPtClosure[iCentrality], "Gen p_{T} (GeV)", "#mu(reco p_{T} / gen p_{T})", ppData, iCentrality, 0, includeQuarkGluon, "PtClosure", saveFigures);
      
      drawClosureHistogram(hJetPtClosureSigma[iCentrality], "Gen p_{T} (GeV)", "#sigma(reco p_{T} / gen p_{T})", ppData, iCentrality, 1, includeQuarkGluon, "PtResolution", saveFigures);
      
      hJetPtClosureSigma[iCentrality][2]->Write();
    }
    
    if(drawEtaClosure){
      drawClosureHistogram(hJetPtClosureEta[iCentrality], "#eta", "#mu(reco p_{T} / gen p_{T})", ppData, iCentrality, 0, includeQuarkGluon, "EtaClosure", saveFigures);
      
      drawClosureHistogram(hJetPtClosureSigmaEta[iCentrality], "#eta", "#sigma(reco p_{T} / gen p_{T})", ppData, iCentrality, 1, includeQuarkGluon, "EtaResolution", saveFigures);
      
    }
    
  } // Centrality loop
  
  outputFile->Close();
  
}
