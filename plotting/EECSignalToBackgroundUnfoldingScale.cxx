/*
 * Implementation of the EECSignalToBackgroundUnfoldingScale class
 */

// Own includes
#include "EECSignalToBackgroundUnfoldingScale.h"


/*
 * Contructor
 */
EECSignalToBackgroundUnfoldingScale::EECSignalToBackgroundUnfoldingScale(){
  InitializeFunctions(1);
}

/*
 * Custom constructor
 *
 *  EECCard* card = Card from which the weight exponent is determined
 */
EECSignalToBackgroundUnfoldingScale::EECSignalToBackgroundUnfoldingScale(EECCard* card){
  TString dataType = card->GetDataType();
  int iWeightExponent = card->GetWeightExponent();
  InitializeFunctions(iWeightExponent);
}

/*
 * Custom constructor
 *
 *  const int iWeightExponent = Exponent given to the energy weight in energy-energy correlators
 *                              Currently implemented values: 1 and 2
 */
EECSignalToBackgroundUnfoldingScale::EECSignalToBackgroundUnfoldingScale(const int iWeightExponent){
  InitializeFunctions(iWeightExponent);
}

/*
 * Initialization for signal to background scale function parameter arrays
 *
 *  const int iWeightExponent = Exponent given to the energy weight in energy-energy correlators
 *                              Currently implemented values: 1 and 2
 */
void EECSignalToBackgroundUnfoldingScale::InitializeFunctions(const int iWeightExponent){
  
  // First, check that the input values are reasonable
  if(iWeightExponent < 1 || iWeightExponent > 2){
    std::cout << "EECSignalToBackgroundUnfoldingScale::ERROR! " << std::endl;
    std::cout << "You have defined your energy weight exponent as " << iWeightExponent << " which is outside of the allowed values 1 and 2." << std::endl;
    std::cout << "The code will now crash. Better luck next time!" << std::endl;
    throw std::invalid_argument("iWeightExponent in EECSignalToBackgroundUnfoldingScale is out of allowed range!");
  }

  // ======================================================== //
  //                                                          //
  //  M     M   EEEEEE   TTTTTTT   H    H    OOOOO    DDDD    //
  //  MMM MMM   E           T      H    H   O     O   D   D   //
  //  M  M  M   EEEEEE      T      HHHHHH   O     O   D    D  //
  //  M     M   E           T      H    H   O     O   D   D   //
  //  M     M   EEEEEE      T      H    H    OOOOO    DDDD    //
  //                                                          //
  // ======================================================== //
  
  // When jet pT is unfolded, the true jet pT in the given bin changes. As signal to background ratio depends on jet pT,
  // we must take into account when normalizing the reflected cone background after unfolding. There are two ingredients
  // that go to this calculation. First, to know how much the mean pT is expected to change, we determine this by
  // calculating the mean pT in the 120-140 bin for generator level jets, and mean generator level pT for 120-140 measured
  // jet pT bin in the jet pT response matrix. The difference in these values tells us how much the mean pT changes due
  // to unfolding. Then we determine from data how much the signal to background ratio depends on the jet pT. For this
  // study, we look at nominal energy weights and squared energy weights separately. We calculate the signal to background
  // ratio in the analysis jet pT bins and in the bins adjacent to these, and then fit a second order polynomial to
  // the determined points to find the trend. In the 50-90% bin a first order polynomial is used instead. When we shift
  // the measured mean pT by the amount determined in the first step and read the function value at that point, we know
  // how much the unfolding changes the signal to background ratio. Then we can return a scaling factor for the reflected
  // cone background to take this properly into account.
  //
  // Used configuration for step 1, determining the shift in mean pT from simulation
  // Macro: findMeanPtDifferenceInResponseMatrix.C
  // Input: PbPbMC2018_GenGen_akFlowJets_4pCentShift_cutBadPhi_optimizedUnfoldingBins_nominalSmear_oneDimensionalResponseMatrix_processed_2024-01-16.root
  // Git hash: 580ee453881f757c251e4352e15edf1b964a4bf8
  //
  // Used configuration for step 2, determining the signal to background trends in data
  // Macro: signalToBackgroundRatio.C
  //
  // For results with reflected cone:
  // Input for nominal weight: eecAnalysis_akFlowJet_nominalEnergyWeight_optimizedUnfoldingBins_fixedCovarianceMatrix_unfoldingWithCovariance_processed_2024-01-23.root
  // Input for squared weight: eecAnalysis_akFlowJet_energyWeightSquared_optimizedUnfoldingBins_fixedCovarianceMatrix_unfoldingWithCovariance_processed_2024-01-23.root
  // Git hash: 580ee453881f757c251e4352e15edf1b964a4bf8
  
  double signalToBackgroundFunctionParameters[2][4][5][3] =
  {

  // ==================================================================== //
  //    First part of the array is for nominal energy weights, pT1*pT2    //
  // ==================================================================== //
  {

   // Numbers using mixed cone method from file eecAnalysis_akFlowJet_nominalEnergyWeight_combinedMixedConeBackground_unfoldingWithNominalSmear_processed_2024-05-28.root
   // Centrality 0-10%
  {{0.212828,0.0102745,-6.82336e-06}, // track pT > 1.0 GeV
  {0.0176997,0.0168176,-6.2929e-06}, // track pT > 1.5 GeV
  {-0.433117,0.0297341,-2.22067e-06}, // track pT > 2.0 GeV
  {-1.37505,0.0544879,1.16051e-05}, // track pT > 2.5 GeV
  {-2.90534,0.0950086,6.19311e-05}}, // track pT > 3.0 GeV
   // Centrality 10-30%
  {{0.226719,0.01518,-5.91267e-06}, // track pT > 1.0 GeV
  {0.0300126,0.0250262,-2.05179e-06}, // track pT > 1.5 GeV
  {-0.377559,0.0432069,1.26082e-05}, // track pT > 2.0 GeV
  {-1.18232,0.0765366,4.77449e-05}, // track pT > 2.5 GeV
  {-3.74227,0.14793,6.87683e-05}}, // track pT > 3.0 GeV
   // Centrality 30-50%
  {{-0.057911,0.0346047,-1.02305e-05}, // track pT > 1.0 GeV
  {-0.687454,0.0617439,-1.10006e-05}, // track pT > 1.5 GeV
  {-1.24756,0.102603,3.05721e-05}, // track pT > 2.0 GeV
  {-1.64531,0.164534,0.000142551}, // track pT > 2.5 GeV
  {0.044957,0.22323,0.000502884}}, // track pT > 3.0 GeV
   // Centrality 50-90%
  {{-1.0163,0.11513,0}, // track pT > 1.0 GeV
  {-3.05475,0.215313,0}, // track pT > 1.5 GeV
  {-9.89126,0.418568,0}, // track pT > 2.0 GeV
  {-20.9565,0.750403,0}, // track pT > 2.5 GeV
  {-45.7433,1.30478,0}} // track pT > 3.0 GeV
  },

  // Old numbers from reflected cone method
  /*
  // Centrality 0-10%
  {{0.227696,0.00970569,-4.76595e-06}, // track pT > 1.0 GeV
   {0.0299163,0.0157925,-3.4004e-06}, // track pT > 1.5 GeV
   {-0.423722,0.028125,2.84121e-07}, // track pT > 2.0 GeV
   {-1.34799,0.0514155,1.11466e-05}, // track pT > 2.5 GeV
   {-3.06551,0.0925417,3.74225e-05}}, // track pT > 3.0 GeV
   // Centrality 10-30%
  {{0.257605,0.0132685,-6.7881e-07}, // track pT > 1.0 GeV
   {0.131351,0.0210007,7.3974e-06}, // track pT > 1.5 GeV
   {-0.299582,0.0379164,1.86976e-05}, // track pT > 2.0 GeV
   {-0.951508,0.0661299,4.96968e-05}, // track pT > 2.5 GeV
   {-2.53978,0.120272,7.91541e-05}}, // track pT > 3.0 GeV
   // Centrality 30-50%
  {{0.0806256,0.0302039,-5.56079e-06}, // track pT > 1.0 GeV
   {-0.30946,0.0517791,-2.52296e-06}, // track pT > 1.5 GeV
   {-0.850391,0.088588,1.2843e-05}, // track pT > 2.0 GeV
   {0.721572,0.114685,0.000155004}, // track pT > 2.5 GeV
   {3.42675,0.152728,0.000367973}}, // track pT > 3.0 GeV
   // Centrality 50-90%
  {{0.318311,0.0875095,0}, // track pT > 1.0 GeV
   {0.761381,0.143391,0}, // track pT > 1.5 GeV
   {2.566,0.221488,0}, // track pT > 2.0 GeV
   {6.56557,0.317407,0}, // track pT > 2.5 GeV
   {13.2054,0.417227,0}} // track pT > 3.0 GeV
  },*/

  // RecoReco MC. To be used only with MC test configuration
   // Centrality 4-14%
  /*{{0.601365,0.00477948,9.92827e-06}, // track pT > 1.0 GeV
   {0.687348,0.0073499,2.20725e-05}, // track pT > 1.5 GeV
   {0.691395,0.0139023,4.35853e-05}, // track pT > 2.0 GeV
   {0.617921,0.0257925,8.32832e-05}, // track pT > 2.5 GeV
   {0.764123,0.0407747,0.000163073}}, // track pT > 3.0 GeV
   // Centrality 14-34%
  {{0.526041,0.0103528,6.7249e-06}, // track pT > 1.0 GeV
   {0.493835,0.0179969,1.55891e-05}, // track pT > 1.5 GeV
   {0.184671,0.0359957,2.35665e-05}, // track pT > 2.0 GeV
   {-0.512562,0.0690165,3.34997e-05}, // track pT > 2.5 GeV
   {-2.44576,0.133254,1.52467e-06}}, // track pT > 3.0 GeV
   // Centrality 34-54%
  {{0.66406,0.0239817,1.14017e-05}, // track pT > 1.0 GeV
   {0.79706,0.0426432,2.40417e-05}, // track pT > 1.5 GeV
   {1.18239,0.0772479,4.6873e-05}, // track pT > 2.0 GeV
   {0.020684,0.162107,-2.57542e-05}, // track pT > 2.5 GeV
   {-1.80113,0.29167,-0.000185037}}, // track pT > 3.0 GeV
   // Centrality 54-94%
  {{1.23812,0.079532,0}, // track pT > 1.0 GeV
   {2.72657,0.133189,0}, // track pT > 1.5 GeV
   {6.57561,0.211957,0}, // track pT > 2.0 GeV
   {16.0042,0.287173,0}, // track pT > 2.5 GeV
   {30.422,0.342131,0}} // track pT > 3.0 GeV
  },*/

  // ===================================================================== //
  //    Second part of the array is for squared energy weights, pT1*pT2    //
  // ===================================================================== //

  {

    // Results with mixed cone method from file eecAnalysis_akFlowJet_energyWeightSquared_combinedMixedConeBackground_unfoldingWithNominalSmear_processed_2024-05-28.root
    // Centrality 0-10%
  {{0.506507,-0.00870152,0.000440677}, // track pT > 1.0 GeV
  {0.346352,-0.0126094,0.000650555}, // track pT > 1.5 GeV
  {-0.0851027,-0.0187053,0.0010397}, // track pT > 2.0 GeV
  {-0.23249,-0.0395869,0.00175816}, // track pT > 2.5 GeV
  {1.45974,-0.103526,0.00301874}}, // track pT > 3.0 GeV
   // Centrality 10-30%
  {{1.27033,-0.0192395,0.000760692}, // track pT > 1.0 GeV
  {2.08461,-0.0368002,0.00113573}, // track pT > 1.5 GeV
  {4.53884,-0.0841682,0.00186867}, // track pT > 2.0 GeV
  {8.50917,-0.163061,0.00308987}, // track pT > 2.5 GeV
  {10.5392,-0.227489,0.00473014}}, // track pT > 3.0 GeV
   // Centrality 30-50%
  {{6.05941,-0.110548,0.00212843}, // track pT > 1.0 GeV
  {10.3371,-0.18996,0.00318678}, // track pT > 1.5 GeV
  {21.081,-0.379161,0.00520371}, // track pT > 2.0 GeV
  {40.7054,-0.707387,0.00839907}, // track pT > 2.5 GeV
  {85.1368,-1.44112,0.0139829}}, // track pT > 3.0 GeV
   // Centrality 50-90%
  {{-117.002,1.69224,0}, // track pT > 1.0 GeV
  {-167.144,2.41799,0}, // track pT > 1.5 GeV
  {-251.738,3.60534,0}, // track pT > 2.0 GeV
  {-350.32,5.13027,0}, // track pT > 2.5 GeV
  {-492.506,7.22929,0}} // track pT > 3.0 GeV
  }

   // Old numbers from reflected cone method
   // Centrality 0-10%
  /*{{0.499383,-0.00837629,0.000428108}, // track pT > 1.0 GeV
   {0.263944,-0.0106703,0.000616011}, // track pT > 1.5 GeV
   {-0.357668,-0.0125132,0.000952519}, // track pT > 2.0 GeV
   {-2.01522,-0.00733771,0.00148347}, // track pT > 2.5 GeV
   {-6.13691,0.0231625,0.00219929}}, // track pT > 3.0 GeV
   // Centrality 10-30%
  {{1.65429,-0.0234934,0.000724683}, // track pT > 1.0 GeV
   {2.41124,-0.0381215,0.00104085}, // track pT > 1.5 GeV
   {3.51714,-0.0592093,0.00157079}, // track pT > 2.0 GeV
   {6.42176,-0.103571,0.00240646}, // track pT > 2.5 GeV
   {12.1206,-0.179236,0.00356544}}, // track pT > 3.0 GeV
   // Centrality 30-50%
  {{8.83531,-0.129564,0.00190212}, // track pT > 1.0 GeV
   {13.3796,-0.195066,0.00266809}, // track pT > 1.5 GeV
   {23.7676,-0.338373,0.00399028}, // track pT > 2.0 GeV
   {48.3489,-0.662604,0.00616829}, // track pT > 2.5 GeV
   {87.4468,-1.17521,0.00918269}}, // track pT > 3.0 GeV
   // Centrality 50-90%
  {{-50.2734,0.922068,0}, // track pT > 1.0 GeV
   {-53.8571,1.09745,0}, // track pT > 1.5 GeV
   {-51.1867,1.25769,0}, // track pT > 2.0 GeV
   {-36.5242,1.33926,0}, // track pT > 2.5 GeV
   {-15.7122,1.35529,0}} // track pT > 3.0 GeV
  }*/

  // RecoReco MC. To be used only with MC test configuration
   // Centrality 4-14%
  /*{{1.53417,-0.0192604,0.000408935}, // track pT > 1.0 GeV
   {1.79577,-0.0259716,0.000581697}, // track pT > 1.5 GeV
   {1.45241,-0.0255574,0.000838887}, // track pT > 2.0 GeV
   {0.431411,-0.0159333,0.00118728}, // track pT > 2.5 GeV
   {-2.18465,0.0205891,0.00154459}}, // track pT > 3.0 GeV
   // Centrality 14-34%
  {{0.03052,0.00600556,0.000557027}, // track pT > 1.0 GeV
   {-1.07382,0.0228738,0.000738911}, // track pT > 1.5 GeV
   {-3.89774,0.0695887,0.000964981}, // track pT > 2.0 GeV
   {-9.7489,0.170855,0.00113737}, // track pT > 2.5 GeV
   {-20.0201,0.353631,0.00110694}}, // track pT > 3.0 GeV
   // Centrality 34-54%
  {{-0.0916579,0.0279501,0.00116282}, // track pT > 1.0 GeV
   {-2.26112,0.0821343,0.00141608}, // track pT > 1.5 GeV
   {-6.24482,0.19242,0.00164656}, // track pT > 2.0 GeV
   {-14.2308,0.404424,0.00156501}, // track pT > 2.5 GeV
   {-19.5052,0.615917,0.00142732}}, // track pT > 3.0 GeV
   // Centrality 54-94%
  {{-38.3245,0.783749,0}, // track pT > 1.0 GeV
   {-39.5182,0.933375,0}, // track pT > 1.5 GeV
   {-35.4382,1.07683,0}, // track pT > 2.0 GeV
   {-23.5364,1.16107,0}, // track pT > 2.5 GeV
   {-8.84147,1.20042,0}} // track pT > 3.0 GeV
  }*/

  };

  // Copy the input values to the class array
  for(int iCentrality = 0; iCentrality < kNCentralityBins; iCentrality++){
    for(int iTrackPt = 0; iTrackPt < kNTrackPtBins; iTrackPt++){
      for(int iParameter = 0; iParameter < 3; iParameter++){
        fFunctionParameters[iCentrality][iTrackPt][iParameter] = signalToBackgroundFunctionParameters[iWeightExponent-1][iCentrality][iTrackPt][iParameter];
      }
    } // Track pT loop
  } // Centrality loop

  // Also initialize the mean jet pT array
  double meanJetPt[4][6] = {
    {108.846,128.777,148.873,168.985,189.065,209.116}, // Centrality 0-10%
    {108.911,128.777,148.895,168.953,189.024,209.128}, // Centrality 10-30%
    {108.935,128.791,148.881,168.952,189.103,209.21},  // Centrality 30-50%
    {108.843,128.78,148.872,168.999,189.117,209.084}   // Centrality 50-90%
  };

  // Mean jet pT array for RecoReco MC
  /*double meanJetPt[4][6] = {
    {108.572,128.58,148.669,168.766,188.944,208.971}, // Centrality 4-14%
    {108.5,128.574,148.752,168.873,188.98,209.017}, // Centrality 14-34%
    {108.465,128.588,148.758,168.914,188.991,209.055}, // Centrality 34-54%
    {108.458,128.622,148.787,168.902,189.009,209.078} // Centrality 54-94%
  };*/

  for(int iCentrality = 0; iCentrality < kNCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < kNJetPtBins; iJetPt++){
      fMeanJetPt[iCentrality][iJetPt] = meanJetPt[iCentrality][iJetPt];
    } // Jet pT loop
  } // Centrality loop

  // Initialize the central and paripheral functions
  fSignalToBackgroundFunction = new TF1("SignalToBackgroundFunction", "pol2", 100, 220);
}

/*
 * Get the background scaling factor for the given bin
 *
 *  const std::pair<double,double> centralityBinBorders = Centrality bin borders for the scaling factor
 *  const std::pair<double,double> jetPtBinBorders = Jet pT bin borders for the scaling factor
 *  const double trackPtLowBorder = Lower track pT bin border for the scaling factor. All tracks above this are used
 *  const bool isPbPbData = True for PbPb, false for pp
 *  int iSystematic = Systematic uncertainty index.
 *                    3: Lower scale for signal-to-background ratio uncertainty
 *                    4: Higher scale for signal-to-background ratio uncertainty
 *                    Any other index: Default value for signal-to-background ratio scaling factor
 *
 *  return: Scaling factor corresponding to the bin with the given bin borders
 */
double EECSignalToBackgroundUnfoldingScale::GetEECSignalToBackgroundUnfoldingScale(const std::pair<double,double> centralityBinBorders, const std::pair<double,double> jetPtBinBorders, const double trackPtLowBorder, const bool isPbPbData, int iSystematic) const{

  // No scaling is needed for pp, signal to background ratio is so good that the effect of scaling is negligible
  if(!isPbPbData) return 1;

  // Decode the systematic uncertainty index. Indices 1 and 2 are reserved for MC based scaling factor systematics.
  // Thus indices 3 and 4 are used to vary the background scaling factor coming from signal-to-background ratio changes
  // in unfolding. These need to be internally propagated back to indices 1 and 2 such that correct array indices are read.
  // Any other systematic uncertainty index is interpreted as default case, which translates to index 0 internally.
  if(iSystematic == 3){
    iSystematic = 1;
  } else if (iSystematic == 4){
    iSystematic = 2;
  } else {
    iSystematic = 0;
  }

  // ******************************************************************** //
  // First, find the bin indices that correspond to the input bin borders //
  // ******************************************************************** //

  // Small number
  const double epsilon = 0.0001;

  // For PbPb data, search if the given centrality bin borders are included in the background scale table
  int centralityIndex = -1;
  for(int iCentrality = 0; iCentrality < kNCentralityBins; iCentrality++){
    if(TMath::Abs(fCentralityBinBorderLow[iCentrality] - centralityBinBorders.first) < epsilon){
      if(TMath::Abs(fCentralityBinBorderHigh[iCentrality] - centralityBinBorders.second) < epsilon){
        centralityIndex = iCentrality;
        break;
      }
    }
    if(TMath::Abs(fCentralityBinBorderLowShifted[iCentrality] - centralityBinBorders.first) < epsilon){
      if(TMath::Abs(fCentralityBinBorderHighShifted[iCentrality] - centralityBinBorders.second) < epsilon){
        centralityIndex = iCentrality;
        break;
      }
    }
  } // Loop over centrality bins included in the scaling table

  // If centrality bin is not found, print an error message and return -1 to show we did not find a proper scaling factor.
  if(centralityIndex == -1){
    std::cout << "EECSignalToBackgroundUnfoldingScale::Error! Centrality bin " << centralityBinBorders.first << "-" << centralityBinBorders.second << " % was not found from the table." << std::endl;
    std::cout << "Cannot provide a background scaling factor. Please check your code for errors." << std::endl;
    return -1;
  }

  // Search if the given track pT bin borders are included in the background scale table
  int trackPtIndex = -1;
  for(int iTrackPt = 0; iTrackPt < kNTrackPtBins; iTrackPt++){
    if(TMath::Abs(fTrackPtBinBorderLow[iTrackPt] - trackPtLowBorder) < epsilon){
      trackPtIndex = iTrackPt;
      break;
    }
  } // Loop over track pT bins included in the scaling table

  // If track pT bin is not found, print an error message and return -1 to show we did not find a proper scaling factor.
  if(trackPtIndex == -1){
    std::cout << "EECSignalToBackgroundUnfoldingScale::Error! Track pT > " << trackPtLowBorder << " GeV bin was not found from the table." << std::endl;
    std::cout << "Cannot provide a background scaling factor. Please check your code for errors." << std::endl;
    return -1;
  }

  // Search if the given jet pT bin borders are included in the background scale table
  int jetPtIndex = -1;
  for(int iJetPt = 0; iJetPt < kNJetPtBins; iJetPt++){
    if(TMath::Abs(fJetPtBinBorderLow[iJetPt] - jetPtBinBorders.first) < epsilon){
      if(TMath::Abs(fJetPtBinBorderHigh[iJetPt] - jetPtBinBorders.second) < epsilon){
        jetPtIndex = iJetPt;
        break;
      }
    }
  } // Loop over jet pT bins included in the scaling table

  // If jet pT bin is not found, print an error message and return -1 to show we did not find a proper scaling factor.
  if(jetPtIndex == -1){
    std::cout << "EECSignalToBackgroundUnfoldingScale::Error! Jet pT bin " << jetPtBinBorders.first << "-" << jetPtBinBorders.second << " GeV was not found from the table." << std::endl;
    std::cout << "Cannot provide a background scaling factor. Please check your code for errors." << std::endl;
    return -1;
  }

  // ************************************************************** //
  // After the correct bin has been found, calculate scaling factor //
  // ************************************************************** //

  // Set the parameters for the signal to background function
  for(int iParameter = 0; iParameter < 3; iParameter++){
    fSignalToBackgroundFunction->SetParameter(iParameter, fFunctionParameters[centralityIndex][trackPtIndex][iParameter]);
  }

  // Calculate the difference on signal to background ratios in original and unfolded mean pT locations
  double shiftedMeanPt = fMeanJetPt[centralityIndex][jetPtIndex] + fRelativeUpshift[iSystematic][centralityIndex]*fMeanJetPt[centralityIndex][jetPtIndex];
  double originalSignalToBackgroundRatio = fSignalToBackgroundFunction->Eval(fMeanJetPt[centralityIndex][jetPtIndex]);
  double unfoldedSignalToBackgroundRatio = fSignalToBackgroundFunction->Eval(shiftedMeanPt);

  return originalSignalToBackgroundRatio / unfoldedSignalToBackgroundRatio;

}
