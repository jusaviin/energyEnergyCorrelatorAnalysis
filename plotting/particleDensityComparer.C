#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Find the minimum and maximum values from a histogram. They must be more extreme than the current values
 *
 *  Arguments: TH1D* histogram = Histogram from which the minimum and maximum values are searched
 *             std::pair<double,double> currentMinMax = The found values need to be more extreme than these to be accepted
 */
std::pair<double,double> findHistogramMinMax(TH1D *histogram, std::pair<double,double> currentMinMax){
  
  // As initial guess, take the given minimum and maximum values
  std::pair<double,double> newMinMax = std::make_pair(currentMinMax.first, currentMinMax.second);
  
  // Loop through all the bins in the histogram and update the minimum and maximum values
  double currentValue, currentError;
  for(int iBin = 1; iBin <= histogram->GetNbinsX(); iBin++){
    currentValue = histogram->GetBinContent(iBin);
    currentError = histogram->GetBinError(iBin);
    if((currentValue-currentError) < newMinMax.first) newMinMax.first = currentValue-currentError;
    if((currentValue+currentError) > newMinMax.second) newMinMax.second = currentValue+currentError;
  }
  
  // Return the minimum and maximum values from the histogram
  return newMinMax;
  
}

/*
 * Macro studying particle densities in different situations
 */
void particleDensityComparer(){

  enum enumDataType{kPythiaHydjet, kMinBiasHydjet, knDataTypes};
  
  // File containing the Pythia+Hydjet simulation result (index 0), and the one containing minimum bias Hydjet result (index 1)
  TString inputFileName[knDataTypes] = {"data/PbPbMC2018_GenGen_eecAnalysis_akFlowJet_MnD_wtaAxis_noTrigger_preprocessed_2022-10-21.root", "data/MinBiasHydjet_RecoGen_eecAnalysis_akFlowJet_MnD_eschemeAxis_noTrigger_preprocessed_2022-10-19.root"};
  // data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJet_updatedMultiplicityAndDensity_eschemeAxis_noTrigger_preprocessed_2022-10-17.root
  // data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJet_updatedMultiplicityAndDensity_wtaAxis_noTrigger_preprocessed_2022-10-17.root
  // data/eecAnalysis_akFlowJets_removeBadAcceptance_wtaAxis_processed_2022-10-25.root
  // data/eecAnalysis_akFlowJets_updatedMultiplicityAndDensity_eschemeAxis_preprocessed_2022-10-17.root
  // data/eecAnalysis_akFlowJets_updatedMultiplicityAndDensity_wtaAxis_preprocessed_2022-10-17.root
  
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
  
  int firstStudiedJetPtBinEEC[knDataTypes] = {0,0};
  int lastStudiedJetPtBinEEC[knDataTypes] = {0,nJetPtBinsEEC[kMinBiasHydjet]}; // Note: Jets integrated over all pT ranges are in nJetPtBinsEEC bin
  
  int firstStudiedTrackPtBinEEC = 0;
  int lastStudiedTrackPtBinEEC = nTrackPtBinsEEC-1;
  
  // Select the types of energy-energy correlators are studied
  bool studyParticleDensityType[EECHistogramManager::knParticleDensityAroundJetAxisTypes];
  studyParticleDensityType[EECHistogramManager::kParticleDensityAroundJetAxis] = false;
  studyParticleDensityType[EECHistogramManager::kParticlePtDensityAroundJetAxis] = false;
  studyParticleDensityType[EECHistogramManager::kParticleDensityAroundJetAxisPtBinned] = true;
  studyParticleDensityType[EECHistogramManager::kParticlePtDensityAroundJetAxisPtBinned] = false;
  
  // Select which plots to draw
  const bool drawReflectedConeToHydjetRatio = true;
  const bool drawReflectedConeToSignalRatio = false;
  const bool drawSubeventWithinSignalConeRatio = false;
  
  bool drawMinBiasToRegularRatio[EECHistograms::knSubeventTypes+1];
  drawMinBiasToRegularRatio[EECHistograms::kPythia] = false;
  drawMinBiasToRegularRatio[EECHistograms::kHydjet] = false;
  drawMinBiasToRegularRatio[EECHistograms::knSubeventTypes] = false;
  
  // Logarithmic axes
  const bool logY = false;
  const double maxX = 0.6;
  
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
    
    // Choose the particle density types to load
    histograms[iDataType]->SetLoadParticleDensityAroundJets(studyParticleDensityType[EECHistogramManager::kParticleDensityAroundJetAxis]);
    histograms[iDataType]->SetLoadParticlePtDensityAroundJets(studyParticleDensityType[EECHistogramManager::kParticlePtDensityAroundJetAxis]);
    histograms[iDataType]->SetLoadParticleDensityAroundJetsPtBinned(studyParticleDensityType[EECHistogramManager::kParticleDensityAroundJetAxisPtBinned]);
    histograms[iDataType]->SetLoadParticlePtDensityAroundJetsPtBinned(studyParticleDensityType[EECHistogramManager::kParticlePtDensityAroundJetAxisPtBinned]);
    
    // Choose the bin ranges
    histograms[iDataType]->SetCentralityBinRange(firstStudiedCentralityBin,lastStudiedCentralityBin);
    histograms[iDataType]->SetJetPtBinRangeEEC(firstStudiedJetPtBinEEC[iDataType],lastStudiedJetPtBinEEC[iDataType]);
    histograms[iDataType]->SetTrackPtBinRangeEEC(firstStudiedTrackPtBinEEC,lastStudiedTrackPtBinEEC);
    
    // Load the histograms from the file
    histograms[iDataType]->LoadProcessedHistograms();
  }
  
  // Particle density histograms from the Pythia+Hydjet simulation
  TH1D* hParticleDensity[EECHistogramManager::knParticleDensityAroundJetAxisTypes][nCentralityBins][nJetPtBinsEEC[kPythiaHydjet]+1][nTrackPtBinsEEC][EECHistograms::knJetConeTypes][EECHistograms::knSubeventTypes+1];
  
  // Particle density histograms from the minimum bias file
  TH1D* hMinBias[EECHistogramManager::knParticleDensityAroundJetAxisTypes][nCentralityBins][nJetPtBinsEEC[kMinBiasHydjet]+1][nTrackPtBinsEEC][EECHistograms::knJetConeTypes][EECHistograms::knSubeventTypes+1];
  
  // Histograms for all different ratios
  TH1D* hReflectedConeToHydjetRatio[EECHistogramManager::knParticleDensityAroundJetAxisTypes][nCentralityBins][nJetPtBinsEEC[kPythiaHydjet]+1][nTrackPtBinsEEC];
  TH1D* hReflectedConeToSignalRatio[EECHistogramManager::knParticleDensityAroundJetAxisTypes][nCentralityBins][nJetPtBinsEEC[kPythiaHydjet]+1][nTrackPtBinsEEC];
  TH1D* hSubeventWithinSignalConeRatio[EECHistogramManager::knParticleDensityAroundJetAxisTypes][nCentralityBins][nJetPtBinsEEC[kPythiaHydjet]+1][nTrackPtBinsEEC][EECHistograms::knSubeventTypes];
  TH1D* hMinBiasToRegularRatio[EECHistogramManager::knParticleDensityAroundJetAxisTypes][nCentralityBins][nJetPtBinsEEC[kPythiaHydjet]+1][nJetPtBinsEEC[kMinBiasHydjet]+1][nTrackPtBinsEEC][EECHistograms::knSubeventTypes+1];
  
  // Initialize the energy-energy correlator histogram array to NULL
  for(int iParticleDensity = 0; iParticleDensity < EECHistogramManager::knParticleDensityAroundJetAxisTypes; iParticleDensity++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        for(int iJetCone = 0; iJetCone < EECHistograms::knJetConeTypes; iJetCone++){
          for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes+1; iSubevent++){
            for(int iJetPt = 0; iJetPt < nJetPtBinsEEC[kPythiaHydjet]+1; iJetPt++){
              hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone][iSubevent] = NULL;
            } // Pythia+Hydjet jet pT loop
            for(int iJetPtMinBias = 0; iJetPtMinBias < nJetPtBinsEEC[kMinBiasHydjet]+1; iJetPtMinBias++){
              hMinBias[iParticleDensity][iCentrality][iJetPtMinBias][iTrackPt][iJetCone][iSubevent] = NULL;
            } // Hydjet only jet pT loop
          } // Subevent loop
        } // Jet cone type loop
        for(int iJetPt = 0; iJetPt < nJetPtBinsEEC[kPythiaHydjet]+1; iJetPt++){
          hReflectedConeToHydjetRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt] = NULL;
          hReflectedConeToSignalRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt] = NULL;
          for(int iJetPtMinBias; iJetPtMinBias < nJetPtBinsEEC[kMinBiasHydjet]+1; iJetPtMinBias++){
            for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes+1; iSubevent++){
              hMinBiasToRegularRatio[iParticleDensity][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iSubevent] = NULL;
            }
          }
          for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes; iSubevent++){
            hSubeventWithinSignalConeRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt][iSubevent] = NULL;
          } // Subevent type loop
        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
  } // Particle density type loop
  
  double normalizationFactor;
  
  // Get the histograms from the histogram manager
  for(int iParticleDensity = 0; iParticleDensity < EECHistogramManager::knParticleDensityAroundJetAxisTypes; iParticleDensity++){
    
    // Only read the selected particle density types
    if(!studyParticleDensityType[iParticleDensity]) continue;
    
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
        for(int iJetCone = 0; iJetCone < EECHistograms::knJetConeTypes; iJetCone++){
          for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes+1; iSubevent++){
            
            // Load different subevents only for Pythia+Hydjet simulation
            if(!card[kPythiaHydjet]->GetAlternativeDataType().Contains("Hydjet") && iSubevent != EECHistograms::knSubeventTypes) continue;
              
            for(int iJetPt = firstStudiedJetPtBinEEC[kPythiaHydjet]; iJetPt <= lastStudiedJetPtBinEEC[kPythiaHydjet]; iJetPt++){
              
              // Find the jet pT normalization factor
              if(iJetPt == histograms[kPythiaHydjet]->GetNJetPtBinsEEC()){
                normalizationFactor = histograms[kPythiaHydjet]->GetJetPtIntegral(iCentrality);
              } else {
                normalizationFactor = histograms[kPythiaHydjet]->GetJetPtIntegral(iCentrality, histograms[kPythiaHydjet]->GetJetPtBinBorderEEC(iJetPt), histograms[kPythiaHydjet]->GetJetPtBinBorderEEC(iJetPt+1));
              }
              
              // Read the particle density histograms and normalize everything to the number of jets
              hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone][iSubevent] = histograms[kPythiaHydjet]->GetHistogramParticleDensityAroundJetAxis(iCentrality, iJetPt, iTrackPt, iJetCone, iParticleDensity, iSubevent);
              hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone][iSubevent]->Scale(1/normalizationFactor);
              
            } // Pythia+Hydjet jet pT loop
            
            for(int iJetPtMinBias = firstStudiedJetPtBinEEC[kMinBiasHydjet]; iJetPtMinBias <= lastStudiedJetPtBinEEC[kMinBiasHydjet]; iJetPtMinBias++){
              
              // Find the jet pT normalization factor
              if(iJetPtMinBias == histograms[kMinBiasHydjet]->GetNJetPtBinsEEC()){
                normalizationFactor = histograms[kMinBiasHydjet]->GetJetPtIntegral(iCentrality);
              } else {
                normalizationFactor = histograms[kMinBiasHydjet]->GetJetPtIntegral(iCentrality, histograms[kMinBiasHydjet]->GetJetPtBinBorderEEC(iJetPtMinBias), histograms[kMinBiasHydjet]->GetJetPtBinBorderEEC(iJetPtMinBias+1));
              }
              
              // Read the particle density histograms and normalize everything to the number of jets
              hMinBias[iParticleDensity][iCentrality][iJetPtMinBias][iTrackPt][iJetCone][iSubevent] = histograms[kMinBiasHydjet]->GetHistogramParticleDensityAroundJetAxis(iCentrality, iJetPtMinBias, iTrackPt, iJetCone, iParticleDensity, iSubevent);
              hMinBias[iParticleDensity][iCentrality][iJetPtMinBias][iTrackPt][iJetCone][iSubevent]->Scale(1/normalizationFactor);
                            
            } // Min bias Hydjet jet pT loop
          } // Subevent loop
        } // Jet cones loop
      } // Track pT loop
    } // Centrality loop
  } // Energy-energy correlator loop
  
  // After everything is loaded, calculate the defined ratios
  for(int iParticleDensity = 0; iParticleDensity < EECHistogramManager::knParticleDensityAroundJetAxisTypes; iParticleDensity++){
    
    // Only read the selected particle density types
    if(!studyParticleDensityType[iParticleDensity]) continue;
    
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
        for(int iJetPt = firstStudiedJetPtBinEEC[kPythiaHydjet]; iJetPt <= lastStudiedJetPtBinEEC[kPythiaHydjet]; iJetPt++){
          
          // Ratio between particle density in reflected cone and for the hydjet particles in the signal cone
          if(drawReflectedConeToHydjetRatio){
            hReflectedConeToHydjetRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt] = (TH1D*) hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kReflectedCone][EECHistograms::kHydjet]->Clone(Form("reflectedConeToHydjetRatio%d%d%d%d", iParticleDensity, iCentrality, iJetPt, iTrackPt));
            hReflectedConeToHydjetRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt]->Divide(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][EECHistograms::kHydjet]);
          }
          
          // Ratio between particle density in reflected cone and in the signal cone
          if(drawReflectedConeToSignalRatio){
            hReflectedConeToSignalRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt] = (TH1D*) hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kReflectedCone][EECHistograms::knSubeventTypes]->Clone(Form("reflectedConeToSignalRatio%d%d%d%d", iParticleDensity, iCentrality, iJetPt, iTrackPt));
            hReflectedConeToSignalRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt]->Divide(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][EECHistograms::knSubeventTypes]);
          }
          
          // Ratio between different subevent selections within the signal cone
          if(drawSubeventWithinSignalConeRatio){
            for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes; iSubevent++){
              hSubeventWithinSignalConeRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt][iSubevent] = (TH1D*) hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][iSubevent]->Clone(Form("reflectedConeToSignalRatio%d%d%d%d%d", iParticleDensity, iCentrality, iJetPt, iTrackPt, iSubevent));
              hSubeventWithinSignalConeRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt][iSubevent]->Divide(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][EECHistograms::knSubeventTypes]);
            } // Subevent loop
          }
          
          // Ratio between minimum bias and particles within the signal cone
          for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes+1; iSubevent++){
            
            // Only construct the ratios for selected subevents
            if(!drawMinBiasToRegularRatio[iSubevent]) continue;
            
            for(int iJetPtMinBias = firstStudiedJetPtBinEEC[kMinBiasHydjet]; iJetPtMinBias <= lastStudiedJetPtBinEEC[kMinBiasHydjet]; iJetPtMinBias++){
              hMinBiasToRegularRatio[iParticleDensity][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iSubevent] = (TH1D*) hMinBias[iParticleDensity][iCentrality][iJetPtMinBias][iTrackPt][EECHistograms::kSignalCone][EECHistograms::knSubeventTypes]->Clone(Form("minBiasToHydjetRatio%d%d%d%d%d", iParticleDensity, iCentrality, iJetPt, iJetPtMinBias, iTrackPt));
              hMinBiasToRegularRatio[iParticleDensity][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iSubevent]->Divide(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][iSubevent]);
            } // Min bias jet pT loop
          } // Loop over subevent types
          
        } // Jet pT loop
      } // Track pT loop
    } // Centrality loop
  } // Energy-energy correlator loop
  
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);

  TLegend *legend, *ptLegend;
  TString centralityString, trackPtString, jetPtString, jetPtStringMinBias;
  TString compactCentralityString, compactTrackPtString, compactJetPtString, compactJetPtStringMinBias;
  int color[10] = {kBlack, kBlue, kRed, kGreen+2, kCyan, kMagenta, kOrange-1, kAzure+5, kOrange-2, kGray};
  std::pair<double,double> histogramYrange;
  TString legendString;

  // Change the legend size based on the number of histograms in the legend
  double trackLegendY1 = 0.03;
  double trackLegendY2 = trackLegendY1 + 0.06*(lastStudiedTrackPtBinEEC+1-firstStudiedTrackPtBinEEC);
  double jetLegendY1 = 0.03;
  double jetLegendY2 = trackLegendY1 + 0.06*(lastStudiedJetPtBinEEC[kPythiaHydjet]+1-firstStudiedJetPtBinEEC[kPythiaHydjet]);
  double minBiasLegendY1 = 0.03;
  double minBiasLegendY2 = minBiasLegendY1 + 0.06*(lastStudiedJetPtBinEEC[kMinBiasHydjet]+2-firstStudiedJetPtBinEEC[kMinBiasHydjet]);
  
  // Loop over all selected histograms
  for(int iParticleDensity = 0; iParticleDensity < EECHistogramManager::knParticleDensityAroundJetAxisTypes; iParticleDensity++){
    
    // Only read the selected energy-energy correlator types
    if(!studyParticleDensityType[iParticleDensity]) continue;
    
    for(int iCentrality = firstStudiedCentralityBin; iCentrality <= lastStudiedCentralityBin; iCentrality++){
      
      // Set the centrality information for legends and figure saving
      centralityString = Form("Cent: %.0f-%.0f%%", histograms[kPythiaHydjet]->GetCentralityBinBorder(iCentrality), histograms[kPythiaHydjet]->GetCentralityBinBorder(iCentrality+1));
      compactCentralityString = Form("_C=%.0f-%.0f", histograms[kPythiaHydjet]->GetCentralityBinBorder(iCentrality), histograms[kPythiaHydjet]->GetCentralityBinBorder(iCentrality+1));
      
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
        
        // Loop over track pT bins
        for(int iTrackPt = firstStudiedTrackPtBinEEC; iTrackPt <= lastStudiedTrackPtBinEEC; iTrackPt++){
          
          // Track pT binning is different depending on the particle density type
          if(iParticleDensity == EECHistogramManager::kParticleDensityAroundJetAxisPtBinned || iParticleDensity == EECHistogramManager::kParticlePtDensityAroundJetAxisPtBinned){
            trackPtString = Form("%.1f < track p_{T} < %.1f",histograms[kPythiaHydjet]->GetTrackPtBinBorderEEC(iTrackPt), histograms[kPythiaHydjet]->GetTrackPtBinBorderEEC(iTrackPt+1));
            compactTrackPtString = Form("_T%.1f-%.1f",histograms[kPythiaHydjet]->GetTrackPtBinBorderEEC(iTrackPt), histograms[kPythiaHydjet]->GetTrackPtBinBorderEEC(iTrackPt+1));
            compactTrackPtString.ReplaceAll(".","v");
          } else {
            trackPtString = Form("%.1f < track p_{T}",histograms[kPythiaHydjet]->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString = Form("_T>%.1f",histograms[kPythiaHydjet]->GetTrackPtBinBorderEEC(iTrackPt));
            compactTrackPtString.ReplaceAll(".","v");
          }
          
          // =============================================================================================================
          // ===   Draw ratio between particle density in reflected cone and the hydjet particles in the signal cone   ===
          // =============================================================================================================
          if(drawReflectedConeToHydjetRatio){
            
            // Create a legend for the figure
            legend = new TLegend(0.58,0.48,0.83,0.84);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->AddEntry((TObject*) 0, histograms[kPythiaHydjet]->GetCard()->GetAlternativeDataType().Data(), "");
            legend->AddEntry((TObject*) 0, centralityString.Data(),"");
            legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
            legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
            
            // Create a new canvas for the plot
            drawer->CreateSplitCanvas();
            
            // Logarithmic EEC axis
            if(logY) drawer->SetLogY(true);

            // Set the x-axis ranges
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kReflectedCone][EECHistograms::kHydjet]->GetXaxis()->SetRangeUser(0.0, maxX);
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][EECHistograms::kHydjet]->GetXaxis()->SetRangeUser(0.0, maxX);
            
            // Find good y-ranges for plotting
            histogramYrange = std::make_pair(10e10, 0);
            histogramYrange = findHistogramMinMax(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kReflectedCone][EECHistograms::kHydjet], histogramYrange);
            histogramYrange = findHistogramMinMax(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][EECHistograms::kHydjet], histogramYrange);
            
            // Add some empty space to the top and the bottom of the histogram
            histogramYrange.first = histogramYrange.first - histogramYrange.first*0.1;
            histogramYrange.second = histogramYrange.second + histogramYrange.second*0.1;
            
            // Set the y-ranges for the histogram
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kReflectedCone][EECHistograms::kHydjet]->GetYaxis()->SetRangeUser(histogramYrange.first, histogramYrange.second);
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][EECHistograms::kHydjet]->GetYaxis()->SetRangeUser(histogramYrange.first, histogramYrange.second);
            
            // Draw the histograms to the upper canvas
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kReflectedCone][EECHistograms::kHydjet]->SetLineColor(color[0]);
            drawer->DrawHistogramToUpperPad(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kReflectedCone][EECHistograms::kHydjet], "#Deltar", histograms[kPythiaHydjet]->GetParticleDensityAroundJetAxisAxisName(iParticleDensity), " ");
            legend->AddEntry(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kReflectedCone][EECHistograms::kHydjet], "Reflected cone", "l");
            
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][EECHistograms::kHydjet]->SetLineColor(color[1]);
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][EECHistograms::kHydjet]->Draw("same");
            legend->AddEntry(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][EECHistograms::kHydjet], "Signal cone, Hydjet", "l");
            
            // Draw the legends to the upper pad
            legend->Draw();
            
            // Linear scale for the ratio
            drawer->SetLogY(false);
            
            // Set the x-axis range for the ratio histogram
            hReflectedConeToHydjetRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.0, maxX);
            
            hReflectedConeToHydjetRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[1]);
            drawer->SetGridY(true);
            hReflectedConeToHydjetRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
            drawer->DrawHistogramToLowerPad(hReflectedConeToHydjetRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt], "#Deltar", "#frac{Reflected cone}{Hydjet}", " ");
            drawer->SetGridY(false);
   
          } // Drawing reflected cone to hydjet ratio
          
          // =====================================================================================
          // ===   Draw ratio between particle density in reflected cone and the signal cone   ===
          // =====================================================================================
          if(drawReflectedConeToSignalRatio){
            
            // Create a legend for the figure
            legend = new TLegend(0.62,0.48,0.87,0.84);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->AddEntry((TObject*) 0, histograms[kPythiaHydjet]->GetCard()->GetAlternativeDataType().Data(), "");
            legend->AddEntry((TObject*) 0, centralityString.Data(),"");
            legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
            legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
            
            // Create a new canvas for the plot
            drawer->CreateSplitCanvas();
            
            // Logarithmic EEC axis
            if(logY) drawer->SetLogY(true);

            // Set the x-axis ranges
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kReflectedCone][EECHistograms::knSubeventTypes]->GetXaxis()->SetRangeUser(0.0, maxX);
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][EECHistograms::knSubeventTypes]->GetXaxis()->SetRangeUser(0.0, maxX);
            
            // Find good y-ranges for plotting
            histogramYrange = std::make_pair(10e10, 0);
            histogramYrange = findHistogramMinMax(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kReflectedCone][EECHistograms::knSubeventTypes], histogramYrange);
            histogramYrange = findHistogramMinMax(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][EECHistograms::knSubeventTypes], histogramYrange);
            
            // Add some empty space to the top and the bottom of the histogram
            histogramYrange.first = histogramYrange.first - histogramYrange.first*0.1;
            histogramYrange.second = histogramYrange.second + histogramYrange.second*0.1;
            
            // Set the y-ranges for the histogram
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kReflectedCone][EECHistograms::knSubeventTypes]->GetYaxis()->SetRangeUser(histogramYrange.first, histogramYrange.second);
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][EECHistograms::knSubeventTypes]->GetYaxis()->SetRangeUser(histogramYrange.first, histogramYrange.second);
            
            // Draw the histograms to the upper canvas
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kReflectedCone][EECHistograms::knSubeventTypes]->SetLineColor(color[0]);
            drawer->DrawHistogramToUpperPad(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kReflectedCone][EECHistograms::knSubeventTypes], "#Deltar", histograms[kPythiaHydjet]->GetParticleDensityAroundJetAxisAxisName(iParticleDensity), " ");
            legend->AddEntry(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kReflectedCone][EECHistograms::knSubeventTypes], "Reflected cone", "l");
            
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][EECHistograms::knSubeventTypes]->SetLineColor(color[1]);
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][EECHistograms::knSubeventTypes]->Draw("same");
            legend->AddEntry(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][EECHistograms::knSubeventTypes], "Signal cone", "l");
            
            // Draw the legends to the upper pad
            legend->Draw();
            
            // Linear scale for the ratio
            drawer->SetLogY(false);
            
            // Set the x-axis range for the ratio histogram
            hReflectedConeToSignalRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(0.0, maxX);
            
            hReflectedConeToSignalRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt]->SetLineColor(color[1]);
            drawer->SetGridY(true);
            hReflectedConeToSignalRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
            drawer->DrawHistogramToLowerPad(hReflectedConeToSignalRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt], "#Deltar", "#frac{Reflected cone}{Signal cone}", " ");
            drawer->SetGridY(false);
            
            // Save the figures to a file
            if(saveFigures){
              gPad->GetCanvas()->SaveAs(Form("figures/%sConeComparison%s%s%s%s.%s", histograms[kPythiaHydjet]->GetParticleDensityAroundJetAxisSaveName(iParticleDensity), saveComment, compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
            }
   
          } // Drawing reflected cone to signal cone ratio
          
          // =========================================================================
          // ===   Draw ratio between different subevents within the signal cone   ===
          // =========================================================================
          if(drawSubeventWithinSignalConeRatio){
            
            // Create a legend for the figure
            legend = new TLegend(0.58,0.42,0.83,0.84);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->AddEntry((TObject*) 0, histograms[kPythiaHydjet]->GetCard()->GetAlternativeDataType().Data(), "");
            legend->AddEntry((TObject*) 0, centralityString.Data(),"");
            legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
            legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
            
            // Create a new canvas for the plot
            drawer->CreateSplitCanvas();
            
            // Logarithmic EEC axis
            if(logY) drawer->SetLogY(true);
            
            // Set the x-axis ranges and find a good y-axis range
            histogramYrange = std::make_pair(10e10, 0);
            for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes+1; iSubevent++){
              hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][iSubevent]->GetXaxis()->SetRangeUser(0.0, maxX);
              
              // Find good y-ranges for plotting
              histogramYrange = findHistogramMinMax(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][iSubevent], histogramYrange);
            }
            
            // Add some empty space to the top and the bottom of the histogram
            histogramYrange.first = histogramYrange.first - histogramYrange.first*0.1;
            histogramYrange.second = histogramYrange.second + histogramYrange.second*0.1;
            
            // Set the y-ranges for the histograms
            for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes+1; iSubevent++){
              hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][iSubevent]->GetYaxis()->SetRangeUser(histogramYrange.first, histogramYrange.second);
            }
            
            // Draw first the histogram without discriminating by subevent
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][EECHistograms::knSubeventTypes]->SetLineColor(color[0]);
            drawer->DrawHistogramToUpperPad(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][EECHistograms::knSubeventTypes], "#Deltar", histograms[kPythiaHydjet]->GetParticleDensityAroundJetAxisAxisName(iParticleDensity), " ");
            legend->AddEntry(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][EECHistograms::knSubeventTypes], "All particles", "l");
            
            // Add to same canvas the densities from different subevents
            for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes; iSubevent++){
              hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][iSubevent]->SetLineColor(color[iSubevent+1]);
              hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][iSubevent]->Draw("same");
              legend->AddEntry(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][iSubevent], Form("Only %s", histograms[kPythiaHydjet]->GetSubeventType(iSubevent)), "l");
            } // Subevent loop
            
            // Draw the legends to the upper pad
            legend->Draw();
            
            // Linear scale for the ratio
            drawer->SetLogY(false);
            
            // Loop over all the subevents
            for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes; iSubevent++){
              
              // Set the x-axis range for the ratio histogram
              hSubeventWithinSignalConeRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt][iSubevent]->GetXaxis()->SetRangeUser(0.0, maxX);
              
              hSubeventWithinSignalConeRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt][iSubevent]->SetLineColor(color[iSubevent+1]);
              drawer->SetGridY(true);
              hSubeventWithinSignalConeRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt][iSubevent]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
              if(iSubevent == 0){
                drawer->DrawHistogramToLowerPad(hSubeventWithinSignalConeRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt][iSubevent], "#Deltar", "#frac{Subevent}{All particles}", " ");
              } else {
                hSubeventWithinSignalConeRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt][iSubevent]->Draw("same");
              }
              drawer->SetGridY(false);
            } // Subevent loop
            
            // Save the figures to a file
            if(saveFigures){
              gPad->GetCanvas()->SaveAs(Form("figures/%sSubeventDecomposition%s%s%s%s.%s", histograms[kPythiaHydjet]->GetParticleDensityAroundJetAxisSaveName(iParticleDensity), saveComment, compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data(), figureFormat));
            }
            
          } // Draw different subevents within a signal cone
          
          // ===================================================================================================================
          // ===   Draw ratio between minimum bias particle densities and Hydjet particle densities within the signal cone   ===
          // ===================================================================================================================
          for(int iSubevent = 0; iSubevent < EECHistograms::knSubeventTypes+1; iSubevent++){
            
            // Only construct the ratios for selected subevents
            if(!drawMinBiasToRegularRatio[iSubevent]) continue;
            
            // Create a legend for the figure
            legend = new TLegend(0.58,0.42,0.83,0.84);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->AddEntry((TObject*) 0, histograms[kPythiaHydjet]->GetCard()->GetAlternativeDataType().Data(), "");
            legend->AddEntry((TObject*) 0, centralityString.Data(),"");
            legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
            legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
            
            // Determine a good string to put to the legend to describe the subevent selection
            legendString = "All particles";
            if(iSubevent < EECHistograms::knSubeventTypes){
              legendString = Form("%s particles",histograms[kMinBiasHydjet]->GetSubeventType(iSubevent));
            }
            
            // Create a new canvas for the plot
            drawer->CreateSplitCanvas();
            
            // Logarithmic EEC axis
            if(logY) drawer->SetLogY(true);
            
            // Set the x-axis ranges and find a good y-axis range
            histogramYrange = std::make_pair(10e10, 0);
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][iSubevent]->GetXaxis()->SetRangeUser(0.0, maxX);
            histogramYrange = findHistogramMinMax(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][iSubevent], histogramYrange);
            
            for(int iJetPtMinBias = firstStudiedJetPtBinEEC[kMinBiasHydjet]; iJetPtMinBias <= lastStudiedJetPtBinEEC[kMinBiasHydjet]; iJetPtMinBias++){
              hMinBias[iParticleDensity][iCentrality][iJetPtMinBias][iTrackPt][EECHistograms::kSignalCone][EECHistograms::knSubeventTypes]->GetXaxis()->SetRangeUser(0.0, maxX);
              
              // Find good y-ranges for plotting
              histogramYrange = findHistogramMinMax(hMinBias[iParticleDensity][iCentrality][iJetPtMinBias][iTrackPt][EECHistograms::kSignalCone][EECHistograms::knSubeventTypes], histogramYrange);
            } // Minimum bias jet pT loop
            
            // Add some empty space to the top and the bottom of the histogram
            histogramYrange.first = histogramYrange.first - histogramYrange.first*0.1;
            histogramYrange.second = histogramYrange.second + histogramYrange.second*0.1;
            
            // Set the y-ranges for the histograms
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][iSubevent]->GetYaxis()->SetRangeUser(histogramYrange.first, histogramYrange.second);
            for(int iJetPtMinBias = firstStudiedJetPtBinEEC[kMinBiasHydjet]; iJetPtMinBias <= lastStudiedJetPtBinEEC[kMinBiasHydjet]; iJetPtMinBias++){
              hMinBias[iParticleDensity][iCentrality][iJetPtMinBias][iTrackPt][EECHistograms::kSignalCone][EECHistograms::knSubeventTypes]->GetYaxis()->SetRangeUser(histogramYrange.first, histogramYrange.second);
            }
            
            // Draw first the histogram from only hydjet particles
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][iSubevent]->SetLineColor(color[0]);
            drawer->DrawHistogramToUpperPad(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][iSubevent], "#Deltar", histograms[kPythiaHydjet]->GetParticleDensityAroundJetAxisAxisName(iParticleDensity), " ");
            legend->AddEntry(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][EECHistograms::kSignalCone][iSubevent], legendString.Data(), "l");
            
            // Add to same canvas the densities from different jet pT bins in minimum bias
            for(int iJetPtMinBias = firstStudiedJetPtBinEEC[kMinBiasHydjet]; iJetPtMinBias <= lastStudiedJetPtBinEEC[kMinBiasHydjet]; iJetPtMinBias++){
              
              // Set the jet pT information for legends and figure saving
              if(iJetPtMinBias == histograms[kMinBiasHydjet]->GetNJetPtBinsEEC()){
                jetPtStringMinBias = Form("Jet p_{T} > %.0f", histograms[kMinBiasHydjet]->GetCard()->GetJetPtCut());
              } else {
                jetPtStringMinBias = Form("%.0f < jet p_{T} < %.0f", histograms[kMinBiasHydjet]->GetJetPtBinBorderEEC(iJetPtMinBias), histograms[kMinBiasHydjet]->GetJetPtBinBorderEEC(iJetPtMinBias+1));
              }
              
              hMinBias[iParticleDensity][iCentrality][iJetPtMinBias][iTrackPt][EECHistograms::kSignalCone][EECHistograms::knSubeventTypes]->SetLineColor(color[iJetPtMinBias+1]);
              hMinBias[iParticleDensity][iCentrality][iJetPtMinBias][iTrackPt][EECHistograms::kSignalCone][EECHistograms::knSubeventTypes]->Draw("same");
              legend->AddEntry(hMinBias[iParticleDensity][iCentrality][iJetPtMinBias][iTrackPt][EECHistograms::kSignalCone][EECHistograms::knSubeventTypes], Form("MB: %s", jetPtStringMinBias.Data()), "l");
            } // Subevent loop
            
            // Draw the legends to the upper pad
            legend->Draw();
            
            // Linear scale for the ratio
            drawer->SetLogY(false);

            // Loop over all the subevents
            for(int iJetPtMinBias = firstStudiedJetPtBinEEC[kMinBiasHydjet]; iJetPtMinBias <= lastStudiedJetPtBinEEC[kMinBiasHydjet]; iJetPtMinBias++){
              
              // Set the x-axis range for the ratio histogram
              hMinBiasToRegularRatio[iParticleDensity][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iSubevent]->GetXaxis()->SetRangeUser(0.0, maxX);
              
              hMinBiasToRegularRatio[iParticleDensity][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iSubevent]->SetLineColor(color[iJetPtMinBias+1]);
              drawer->SetGridY(true);
              hMinBiasToRegularRatio[iParticleDensity][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iSubevent]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
              if(iJetPtMinBias == 0){
                drawer->DrawHistogramToLowerPad(hMinBiasToRegularRatio[iParticleDensity][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iSubevent], "#Deltar", Form("#frac{Minimum bias}{%s}", legendString.Data()), " ");
              } else {
                hMinBiasToRegularRatio[iParticleDensity][iCentrality][iJetPt][iJetPtMinBias][iTrackPt][iSubevent]->Draw("same");
              }
              drawer->SetGridY(false);
            } // Min bias jet pT loop
          } // Subevent loop for drawing minimum bias to hydjet ratio
        } // Track pT loop
      } // Jet pT loop
    } // Centrality loop
  } // Particle density type loop


}
