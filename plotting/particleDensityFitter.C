#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Define a piecewise linear function for three pieces.
 *  parameters: par[0] = First piece divider. Needs to be fixed.
 *              par[1] = Second piece divider. Needs to be fixed.
 *              par[2] = Constant for the first piece
 *              par[3] = Slope for the first piece
 *              par[4] = Constant for the second piece
 *              par[5] = Slope for the second piece
 *              par[6] = Constant for the third piece
 *              par[7] = Slope for the third piece
 */
double threePieceLinear(double *x, double *par){
  if(x[0] < par[0]) return par[2] + par[3]*x[0];
  if(x[0] > par[1]) return par[6] + par[7]*x[0];
  return par[4] + par[5]*x[0];
}

/*
 * Define a piecewise linear function for two pieces.
 *  parameters: par[0] = Piece divider. Needs to be fixed.
 *              par[1] = Constant for the first piece
 *              par[2] = Slope for the first piece
 *              par[3] = Constant for the second piece
 *              par[4] = Slope for the second piece
 */
double twoPieceLinear(double *x, double *par){
  if(x[0] < par[0]) return par[1] + par[2]*x[0];
  return par[3] + par[4]*x[0];
}

/*
 * Define a  linear function
 *  parameters: par[0] = Constant
 *              par[1] = Slope
 */
double linear(double *x, double *par){
  return par[0] + par[1]*x[0];
}

/*
 * Define a function that is exponential at low DeltaR and linear in high
 *  parameters: par[0] = Divider between linear and exponential components
 *              par[1] = Constant for the exponent
 *              par[2] = Scale for the exponent
 *              par[3] = Exponent
 *              par[4] = Contant for the linear part
 *              par[5] = Slope for the linear part
 */
double expoLinear(double *x, double *par){
  if(x[0] < par[0]) return par[1] + par[2] * TMath::Exp(-par[3] * x[0]);
  return par[4] + par[5] * x[0];
}

/*
 * Macro studying particle densities in different situations
 */
void particleDensityFitter(){

  // Open the input file
  TString inputFileName = "data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJet_updatedMultiplicityAndDensity_wtaAxis_noTrigger_preprocessed_2022-10-17.root";
  // eecAnalysis_akFlowJets_removeBadAcceptance_wtaAxis_processed_2022-10-25.root
  // PbPbMC2018_RecoGen_eecAnalysis_akFlowJet_updatedMultiplicityAndDensity_wtaAxis_noTrigger_preprocessed_2022-10-17.root
  // PbPbMC2018_RecoGen_eecAnalysis_akFlowJet_moreTrackPtBins_noEEC_wtaAxis_noTrigger_processed_2022-11-08.root
  TFile* inputFile = TFile::Open(inputFileName);
  
  // Check that the files exist
  if(inputFile == NULL){
    cout << "Error! The file " << inputFileName.Data() << " does not exist!" << endl;
    cout << "Maybe you forgot the data/ folder path?" << endl;
    cout << "Will not execute the code" << endl;
    return;
  }

  EECCard* card = new EECCard(inputFile);
  
  // ====================================================
  //    Binning configuration for the background study
  // ====================================================
  
  // Find the number of bins from the card
  const int nCentralityBins = card->GetNCentralityBins();
  const int nJetPtBinsEEC = card->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = card->GetNTrackPtBinsEEC();
  
  // Default binning ranges for reference
  // centrality = {0,10,30,50,90}
  // track pT =  {0.7,1,2,3,4,8,12,300}
  // jet pT in energy-energy correlators = {120,140,160,180,200,300,500,5030}
  // track pT in energy-energy correlators = {0.7,1,1.5,2,2.5,3,3.5,4,5,6,7,8,10,12,16,20,300}
  
  // Select which histograms are fitted
  int firstDrawnCentralityBin = 0;
  int lastDrawnCentralityBin = 0;
  
  int firstDrawnJetPtBinEEC = 0;
  int lastDrawnJetPtBinEEC = 0; // Note: Jets integrated over all pT ranges are in nJetPtBinsEEC bin
  
  int firstDrawnTrackPtBinEEC = 0;
  int lastDrawnTrackPtBinEEC = nTrackPtBinsEEC-1;
  
  // Select the types of particle density types that are drawn
  bool drawParticleDensityType[EECHistogramManager::knParticleDensityAroundJetAxisTypes];
  drawParticleDensityType[EECHistogramManager::kParticleDensityAroundJetAxis] = false;
  drawParticleDensityType[EECHistogramManager::kParticlePtDensityAroundJetAxis] = false;
  drawParticleDensityType[EECHistogramManager::kParticleDensityAroundJetAxisPtBinned] = true;
  drawParticleDensityType[EECHistogramManager::kParticlePtDensityAroundJetAxisPtBinned] = false;
  
  // Select the distribution that is fitted. Idea is that we can fit only Hydjet in MC for testing purposes. For data this must be knSubeventTypes
  int subeventIndex = EECHistograms::knSubeventTypes; // kPythia, kHydjet, knSubeventTypes
  
  // Flag for which distribution is fitted
  const bool fitFakeInSignalCone = (subeventIndex == EECHistograms::kHydjet);
  
  // Logarithmic axes
  const bool logY = false;
  
  // Axis zooming
  std::pair<double,double> ratioZoom = std::make_pair(0.9, 1.1);
  
  // Figure saving
  const bool saveFigures = false;  // Save figures
  const char* saveComment = "";   // Comment given for this specific file
  const char* figureFormat = "pdf"; // Format given for the figures
  
  // Create and setup a new histogram managers to project and handle the histograms
  EECHistogramManager *histograms;
  histograms = new EECHistogramManager(inputFile,card);
    
  // Choose the particle density types to load
  histograms->SetLoadParticleDensityAroundJets(true);
  histograms->SetLoadParticlePtDensityAroundJets(true);
  histograms->SetLoadParticleDensityAroundJetsPtBinned(true);
  histograms->SetLoadParticlePtDensityAroundJetsPtBinned(true);
    
  // Choose the bin ranges
  histograms->SetCentralityBinRange(0,nCentralityBins-1);
  histograms->SetJetPtBinRangeEEC(0,nJetPtBinsEEC);
  histograms->SetTrackPtBinRangeEEC(0,nTrackPtBinsEEC-1);
    
  // Load the histograms from the file
  histograms->LoadProcessedHistograms();
  
  // Particle density histograms from the data
  TH1D* hParticleDensity[EECHistogramManager::knParticleDensityAroundJetAxisTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC][EECHistograms::knJetConeTypes];
  TH1D* hParticleDensityToFitRatio[EECHistogramManager::knParticleDensityAroundJetAxisTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC][EECHistograms::knJetConeTypes];
  TF1 *fitFunction[EECHistogramManager::knParticleDensityAroundJetAxisTypes][nCentralityBins][nJetPtBinsEEC+1][nTrackPtBinsEEC][EECHistograms::knJetConeTypes];
  
  // Initialize the energy-energy correlator histogram array to NULL
  for(int iParticleDensity = 0; iParticleDensity < EECHistogramManager::knParticleDensityAroundJetAxisTypes; iParticleDensity++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        for(int iJetCone = 0; iJetCone < EECHistograms::knJetConeTypes; iJetCone++){
          for(int iJetPt = 0; iJetPt < nJetPtBinsEEC+1; iJetPt++){
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone] = NULL;
            hParticleDensityToFitRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone] = NULL;
            
            // Initialize the fit functions in case we are fitting the fake distributions in signal cone in Monte Carlo
            if(fitFakeInSignalCone){
              if(iTrackPt < 5){
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone] = new TF1(Form("fitFunction%d%d%d%d%d", iParticleDensity, iCentrality, iTrackPt, iJetPt, iJetPt), threePieceLinear, 0, 0.6, 8);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->FixParameter(0,0.34);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->FixParameter(1,0.45);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(2,1);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(3,0);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(4,1);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(5,0);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(6,1);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(7,0);
              } else if(iTrackPt < 10){
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone] = new TF1(Form("fitFunction%d%d%d%d%d", iParticleDensity, iCentrality, iTrackPt, iJetPt, iJetPt), twoPieceLinear, 0, 0.6, 5);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->FixParameter(0,0.1);
                if(iTrackPt == 6) fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->FixParameter(0,0.14);
                if(iTrackPt == 7 || iTrackPt == 8) fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->FixParameter(0,0.12);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(1,1);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(2,0);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(3,1);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(4,0);
              } else {
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone] = new TF1(Form("fitFunction%d%d%d%d%d", iParticleDensity, iCentrality, iTrackPt, iJetPt, iJetPt), expoLinear, 0, 0.6, 6);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->FixParameter(0,0.1);
                if(iTrackPt == nTrackPtBinsEEC - 1) fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->FixParameter(0,0.05);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(1,0);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(2,1);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(3,1);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(4,1);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(5,0);
                if(iTrackPt == nTrackPtBinsEEC - 1) {
                  fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->FixParameter(5,0);
                }
              }
            } else {
              // Initialize the fit functions in case we are fitting the total distribution as is done for data
              if(iTrackPt < 4 && iJetCone == EECHistograms::kSignalCone){
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone] = new TF1(Form("fitFunction%d%d%d%d%d", iParticleDensity, iCentrality, iTrackPt, iJetPt, iJetPt), threePieceLinear, 0, 0.6, 8);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->FixParameter(0,0.34);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->FixParameter(1,0.45);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(2,1);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(3,0);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(4,1);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(5,0);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(6,1);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(7,0);
                
              } else {
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone] = new TF1(Form("fitFunction%d%d%d%d%d", iParticleDensity, iCentrality, iTrackPt, iJetPt, iJetPt), linear, 0, 0.6, 2);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(0,1);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(1,0);
              }
            }
            } // Jet pT loop
          } // Jet cone type loop
        } // Track pT loop
      } // Centrality loop
  } // Particle density type loop
  
  double normalizationFactor;
  double lowFitRegion[2] = {0,0.05};
  double highFitRegion[2] = {0.6,0.5};
  double lowFitValue, highFitValue;
  
  // Get the histograms from the histogram manager and fit them with a constant
  for(int iParticleDensity = 0; iParticleDensity < EECHistogramManager::knParticleDensityAroundJetAxisTypes; iParticleDensity++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        for(int iJetCone = 0; iJetCone < EECHistograms::knJetConeTypes; iJetCone++){
          for(int iJetPt = 0; iJetPt <= nJetPtBinsEEC; iJetPt++){
            
            // Find the jet pT normalization factor
            if(iJetPt == histograms->GetNJetPtBinsEEC()){
              normalizationFactor = histograms->GetJetPtIntegral(iCentrality);
            } else {
              normalizationFactor = histograms->GetJetPtIntegral(iCentrality, histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1));
            }
            
            // Read the particle density histograms and normalize integrals to one
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone] = histograms->GetHistogramParticleDensityAroundJetAxis(iCentrality, iJetPt, iTrackPt, iJetCone, iParticleDensity, subeventIndex);
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->Scale(1/normalizationFactor);
            hParticleDensityToFitRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone] = (TH1D*) hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->Clone(Form("distributionToFitRatio%d%d%d%d%d", iParticleDensity, iCentrality, iTrackPt, iJetCone, iJetPt));
            
            // Fit a constant to the tail of the distribution, use region 0.45 < DeltaR < 0.6 for now. This can be optimized later if the approach works
            lowFitValue = lowFitRegion[iJetCone];
            highFitValue = highFitRegion[iJetCone];
            if(iJetCone == EECHistograms::kSignalCone && !fitFakeInSignalCone && iTrackPt >= 4){
              lowFitValue = 0.45;
              highFitValue = 0.6;
            }
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->Fit(fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone],"Q","",lowFitValue,highFitValue);
            
            // Divide the distribution with a with a fit the get the ratio
            hParticleDensityToFitRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->Divide(fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]);
            
          } // Jet pT loop
        } // Jet cone type loop
      } // Track pT loop
    } // Centrality loop
  } // Particle density type loop
  
  // After the fit, adjust the pieces of the fit function such that:
  // 1) The function evaluates to 1 at 0.4 by applying a total weight for the function
  // 2) The function is continuous in the places where the pieces are glued together. This is obtained by adjusting the contants.
  double currentValue, targetValue, neededAdjustment, currentParameter;
  double limitValue;
  double epsilon = 0.0000001;
  for(int iParticleDensity = 0; iParticleDensity < EECHistogramManager::knParticleDensityAroundJetAxisTypes; iParticleDensity++){
    for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
      for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
        for(int iJetCone = 0; iJetCone < EECHistograms::knJetConeTypes; iJetCone++){
          for(int iJetPt = 0; iJetPt <= nJetPtBinsEEC; iJetPt++){
            
            // Scaling the function parameters to be weights when doing the study fitting fake distribution in signal cone
            if(fitFakeInSignalCone){
              
              // 1) Ensure that function evaluates to 1 at 0.4
              targetValue = 1;
              currentValue = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->Eval(0.4);
              neededAdjustment = currentValue / targetValue;
              
              // 2) For three piece linear linear function, change the constant values for different pieces such that the function is continuous
              if(iTrackPt < 5){
                
                // Scale each parameter in the function with the needed adjustment
                for(int iParameter = 2; iParameter < 8; iParameter++){
                  currentParameter = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetParameter(iParameter);
                  fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(iParameter, currentParameter / neededAdjustment);
                }
                
                // Adjust the parts that do not contain 0.4 to match the part that contains it in the place where fit region changes
                limitValue = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetParameter(0);
                targetValue = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->Eval(limitValue + epsilon);
                currentValue = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->Eval(limitValue - epsilon);
                neededAdjustment = currentValue - targetValue;
                currentParameter = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetParameter(2);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(2, currentParameter - neededAdjustment);
                
                limitValue = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetParameter(1);
                targetValue = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->Eval(limitValue - epsilon);
                currentValue = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->Eval(limitValue + epsilon);
                neededAdjustment = currentValue - targetValue;
                currentParameter = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetParameter(6);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(6, currentParameter - neededAdjustment);
                
              } else if(iTrackPt < 10){
                // Mid-high pT behavior
                
                // Scale each parameter in the function with the needed adjustment
                for(int iParameter = 1; iParameter < 5; iParameter++){
                  currentParameter = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetParameter(iParameter);
                  fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(iParameter, currentParameter / neededAdjustment);
                }
                
                // Adjust the parts that do not contain 0.4 to match the part that contains it in the place where fit region changes
                limitValue = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetParameter(0);
                targetValue = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->Eval(limitValue + epsilon);
                currentValue = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->Eval(limitValue - epsilon);
                neededAdjustment = currentValue - targetValue;
                currentParameter = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetParameter(1);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(1, currentParameter - neededAdjustment);
                
              } else {
                // Do the same for the exponential+linear fit
                
                // Scale each parameter in the function with the needed adjustment
                for(int iParameter = 1; iParameter < 6; iParameter++){
                  if(iParameter == 3) continue; // No scaling for parameter number 3, since this is the exponent
                  currentParameter = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetParameter(iParameter);
                  fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(iParameter, currentParameter / neededAdjustment);
                }
                
                // Adjust the parts that do not contain 0.4 to match the part that contains it in the place where fit region changes
                limitValue = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetParameter(0);
                targetValue = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->Eval(limitValue + epsilon);
                currentValue = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->Eval(limitValue - epsilon);
                neededAdjustment = currentValue - targetValue;
                currentParameter = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetParameter(1);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(1, currentParameter - neededAdjustment);
                
              }
            } else {
              // Scaling the function parameters to be weights when doing the study fitting the tails of the signal cone as in data
              
              // 1) Ensure that function evaluates to 1 at 0.4
              targetValue = 1;
              currentValue = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->Eval(0.4);
              neededAdjustment = currentValue / targetValue;
              
              // 2) For three piece linear linear function, change the constant values for different pieces such that the function is continuous
              if(iTrackPt < 4 && iJetCone == EECHistograms::kSignalCone){
                
                // Scale each parameter in the function with the needed adjustment
                for(int iParameter = 2; iParameter < 8; iParameter++){
                  currentParameter = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetParameter(iParameter);
                  fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(iParameter, currentParameter / neededAdjustment);
                }
                
                // Set the slope of the linear function at low deltaR to the same slope it is at high deltaR
                currentParameter = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetParameter(7);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(3, currentParameter);
                
                // Adjust the parts that do not contain 0.4 to match the part that contains it in the place where fit region changes
                limitValue = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetParameter(0);
                targetValue = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->Eval(limitValue + epsilon);
                currentValue = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->Eval(limitValue - epsilon);
                neededAdjustment = currentValue - targetValue;
                currentParameter = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetParameter(2);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(2, currentParameter - neededAdjustment);
                
                limitValue = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetParameter(1);
                targetValue = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->Eval(limitValue - epsilon);
                currentValue = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->Eval(limitValue + epsilon);
                neededAdjustment = currentValue - targetValue;
                currentParameter = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetParameter(6);
                fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(6, currentParameter - neededAdjustment);
                
              } else {
                // Otherwise we are dealing with a simple linear fit, so just the initial scaling is needed
                
                // Scale each parameter in the function with the needed adjustment
                for(int iParameter = 0; iParameter < 2; iParameter++){
                  currentParameter = fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetParameter(iParameter);
                  fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetParameter(iParameter, currentParameter / neededAdjustment);
                }
                
              }
              
            }
          } // Jet pT loop
        } // Jet cone type loop
      } // Track pT loop
    } // Centrality loop
  } // Particle density type loop
  
  // =======================================================
  //    Collect the slope parameters from the linear fits
  // =======================================================
  for(int iParticleDensity = 0; iParticleDensity < EECHistogramManager::knParticleDensityAroundJetAxisTypes; iParticleDensity++){
    if(iParticleDensity != EECHistogramManager::kParticleDensityAroundJetAxisPtBinned) continue;
    for(int iJetCone = 0; iJetCone < EECHistograms::knJetConeTypes; iJetCone++){
      if(iJetCone > 0) continue;
      cout << " double " << histograms->GetParticleDensityAroundJetAxisSaveName(iParticleDensity) << histograms->GetJetConeTypeSaveName(iJetCone) << "Parameters[8][" << nCentralityBins << "][" << nJetPtBinsEEC << "][" << nTrackPtBinsEEC << "] = {" << endl;
      for(int iParameter = 0; iParameter < 8; iParameter++){
        cout << " // ============= Paramater " << iParameter << " =============" << endl << "{" << endl;
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          cout << " // " << Form("Centrality %.0f-%.0f%%", histograms->GetCentralityBinBorder(iCentrality), histograms->GetCentralityBinBorder(iCentrality+1)) << endl;
          cout << "{";
          for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
            cout << "{";
            for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
              if(fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone] == NULL){
                cout << "-1";
              } else {
                cout << fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetParameter(iParameter);
              }
              if(iTrackPt < nTrackPtBinsEEC-1){
                cout << ",";
              }
            } // Track pT loop
            cout << "}";
            if(iJetPt < nJetPtBinsEEC-1){
              cout << ", // " << Form("%.0f < jet pT < %.0f GeV", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1)) << endl;
            } else {
              cout << "}";
              if(iCentrality < nCentralityBins-1){
                cout << ", // " << Form("%.0f < jet pT < %.0f GeV", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1)) << endl;
              } else {
                cout << " // " << Form("%.0f < jet pT < %.0f GeV", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1)) << endl;
              }
            }
          } // Jet pT loop
        } // Centrality loop
        if(iParameter == 7){
          cout << "}" << endl;
        } else {
          cout << "}," << endl;
        }
      } // Parameter index loop
      cout << "};" << endl;
    } // Jet cone type loop
  } // Particle density type loop
  
//  // =======================================================
//  //    Collect the slope parameters from the linear fits
//  // =======================================================
//  for(int iParticleDensity = 0; iParticleDensity < EECHistogramManager::knParticleDensityAroundJetAxisTypes; iParticleDensity++){
//    if(iParticleDensity != EECHistogramManager::kParticleDensityAroundJetAxisPtBinned) continue;
//    for(int iJetCone = 0; iJetCone < EECHistograms::knJetConeTypes; iJetCone++){
//      if(iJetCone > 0) continue;
//      for(int iParameter = 0; iParameter < 8; iParameter++){
//        cout << " double " << histograms->GetParticleDensityAroundJetAxisSaveName(iParticleDensity) << histograms->GetJetConeTypeSaveName(iJetCone) << "Parameter" << iParameter << "[" << nCentralityBins << "][" << nJetPtBinsEEC << "][" << nTrackPtBinsEEC << "] = {" << endl;
//        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
//          cout << " // " << Form("Centrality %.0f-%.0f%%", histograms->GetCentralityBinBorder(iCentrality), histograms->GetCentralityBinBorder(iCentrality+1)) << endl;
//          cout << "{";
//          for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
//            cout << "{";
//            for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
//              if(fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone] == NULL){
//                cout << "-1";
//              } else if(iTrackPt == nTrackPtBinsEEC - 1 && iParameter > 5){
//                cout << "0";
//              } else {
//                cout << fitFunction[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetParameter(iParameter);
//              }
//              if(iTrackPt < nTrackPtBinsEEC-1){
//                cout << ",";
//              }
//            } // Track pT loop
//            cout << "}";
//            if(iJetPt < nJetPtBinsEEC-1){
//              cout << ", // " << Form("%.0f < jet pT < %.0f GeV", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1)) << endl;
//            } else {
//              cout << "}";
//              if(iCentrality < nCentralityBins-1){
//                cout << ", // " << Form("%.0f < jet pT < %.0f GeV", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1)) << endl;
//              } else {
//                cout << " // " << Form("%.0f < jet pT < %.0f GeV", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1)) << endl;
//              }
//            }
//          } // Jet pT loop
//        } // Centrality loop
//        cout << "};" << endl;
//      } // Parameter index loop
//    } // Jet cone type loop
//  } // Particle density type loop
  
  // ==========================================================================
  //    Draw the particle density histograms to check the quality of the fit
  // ==========================================================================
  
  JDrawer *drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);

  TLegend *legend;
  TString centralityString, trackPtString, jetPtString;
  TString compactCentralityString, compactTrackPtString, compactJetPtString;
  TString legendString;
  
  // Loop over all selected histograms
  for(int iParticleDensity = 0; iParticleDensity < EECHistogramManager::knParticleDensityAroundJetAxisTypes; iParticleDensity++){
    
    // Only read the selected energy-energy correlator types
    if(!drawParticleDensityType[iParticleDensity]) continue;
    
    for(int iJetCone = 0; iJetCone < EECHistograms::knJetConeTypes; iJetCone++){
      
      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        
        // Set the centrality information for legends and figure saving
        centralityString = Form("Cent: %.0f-%.0f%%", histograms->GetCentralityBinBorder(iCentrality), histograms->GetCentralityBinBorder(iCentrality+1));
        compactCentralityString = Form("_C=%.0f-%.0f", histograms->GetCentralityBinBorder(iCentrality), histograms->GetCentralityBinBorder(iCentrality+1));
        
        // Loop over jet pT bins
        for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
          
          // Set the jet pT information for legends and figure saving
          if(iJetPt == histograms->GetNJetPtBinsEEC()){
            jetPtString = Form("Jet p_{T} > %.0f", histograms->GetCard()->GetJetPtCut());
            compactJetPtString = "";
          } else {
            jetPtString = Form("%.0f < jet p_{T} < %.0f", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1));
            compactJetPtString = Form("_J=%.0f-%.0f", histograms->GetJetPtBinBorderEEC(iJetPt), histograms->GetJetPtBinBorderEEC(iJetPt+1));
          }
          
          // Loop over track pT bins
          for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
            
            // Track pT binning is different depending on the particle density type
            if(iParticleDensity == EECHistogramManager::kParticleDensityAroundJetAxisPtBinned || iParticleDensity == EECHistogramManager::kParticlePtDensityAroundJetAxisPtBinned){
              trackPtString = Form("%.1f < track p_{T} < %.1f",histograms->GetTrackPtBinBorderEEC(iTrackPt), histograms->GetTrackPtBinBorderEEC(iTrackPt+1));
              compactTrackPtString = Form("_T%.1f-%.1f",histograms->GetTrackPtBinBorderEEC(iTrackPt), histograms->GetTrackPtBinBorderEEC(iTrackPt+1));
              compactTrackPtString.ReplaceAll(".","v");
            } else {
              trackPtString = Form("%.1f < track p_{T}",histograms->GetTrackPtBinBorderEEC(iTrackPt));
              compactTrackPtString = Form("_T>%.1f",histograms->GetTrackPtBinBorderEEC(iTrackPt));
              compactTrackPtString.ReplaceAll(".","v");
            }
            
            // Create a legend for the figure
            legend = new TLegend(0.58,0.48,0.83,0.84);
            legend->SetFillStyle(0);legend->SetBorderSize(0);legend->SetTextSize(0.05);legend->SetTextFont(62);
            legend->AddEntry((TObject*) 0, histograms->GetCard()->GetAlternativeDataType().Data(), "");
            legend->AddEntry((TObject*) 0, centralityString.Data(),"");
            legend->AddEntry((TObject*) 0, jetPtString.Data(),"");
            legend->AddEntry((TObject*) 0, trackPtString.Data(),"");
            if(iJetCone == EECHistograms::kReflectedCone) legend->AddEntry((TObject*) 0,"Reflected cone","");
            
            // Create a new canvas for the plot
            drawer->CreateSplitCanvas();
            
            // Logarithmic EEC axis
            if(logY) drawer->SetLogY(true);
            
            // Set the x-axis range
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetXaxis()->SetRangeUser(0.0, 0.6);
            
            // Draw the histograms to the upper canvas
            hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetLineColor(kBlack);
            drawer->DrawHistogramToUpperPad(hParticleDensity[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone], "#Deltar", histograms->GetParticleDensityAroundJetAxisAxisName(iParticleDensity), " ");
            
            // Draw the legends to the upper pad
            legend->Draw();
            
            // Linear scale for the ratio
            drawer->SetLogY(false);
            
            // Set the x-axis range for the ratio histogram
            hParticleDensityToFitRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetXaxis()->SetRangeUser(0.0, 0.6);
            
            hParticleDensityToFitRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->SetLineColor(kBlack);
            drawer->SetGridY(true);
            hParticleDensityToFitRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone]->GetYaxis()->SetRangeUser(ratioZoom.first, ratioZoom.second);
            drawer->DrawHistogramToLowerPad(hParticleDensityToFitRatio[iParticleDensity][iCentrality][iJetPt][iTrackPt][iJetCone], "#Deltar", "#frac{Distribution}{Fit}", " ");
            drawer->SetGridY(false);
            
            
          } // Track pT loop
        } // Jet pT loop
      } // Centrality loop
    } // Jet cone loop
  } // Particle density type loop


}
