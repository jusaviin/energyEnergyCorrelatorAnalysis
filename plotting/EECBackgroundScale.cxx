/*
 * Implementation of the EECBackgroundScale class
 */

// Own includes
#include "EECBackgroundScale.h"


/*
 * Contructor
 */
EECBackgroundScale::EECBackgroundScale(){
  InitializeArrays();
}

/*
 * Destructor
 */
EECBackgroundScale::~EECBackgroundScale(){

}

/*
 * Initialization for seagull method arrays. Set all bins that need nonzero value
 */
void EECBackgroundScale::InitializeArrays(){
  
  // ======================================================== //
  //                                                          //
  //  M     M   EEEEEE   TTTTTTT   H    H    OOOOO    DDDD    //
  //  MMM MMM   E           T      H    H   O     O   D   D   //
  //  M  M  M   EEEEEE      T      HHHHHH   O     O   D    D  //
  //  M     M   E           T      H    H   O     O   D   D   //
  //  M     M   EEEEEE      T      H    H    OOOOO    DDDD    //
  //                                                          //
  // ======================================================== //
  
  // In the current method, we start with Pythia+Hydjet simulation where the reconstructed jets are matched with the
  // generator level jets, and it is required that the matched jet has at least 50% of the reconstructed jet pT.
  // Since generator level jets are clustered only from Pythia particles, it can with large confidence say that
  // matched reconstructed jets are also Pythia jets, instead of jets coming from the Hydjet simulation.
  // When doing thi smatching, we can trust that the Hydjet particles are mostly background, and Pythia particles mostly signal.
  // We determine the true background from Pythia+Hydjet and Hydjet+Hydjet pairings from the signal cone.
  // The reflected cone background estimate is obtained by combining all particles from signal cone with all particles from the reflected cone.
  // We know that in this case fake+fake contribution is not properly modeled, but this contribution is negligible in the region
  // Track pT > 2 GeV and DeltaR < 0.4. We also know that there is a bias in reconstruction jets on upwards fluctuations of background.
  // To take the fluctuations into account and to mitigate the mismodeling of the fake+fake part, we look at the integrals of
  // true background and the reflected cone background estimate in the region DeltaR < 0.4 and calculate their ratio.
  // This ratio now reflects excess fluctuations below the jet peak in Monte Carlo. Assuming the ratio is similar in data,
  // the reflected cone estimate from data can be scaled with this number to get a good estimate of the absolute background level.
  //
  // The file that was used for the study is PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_finalMcWeight_matchJets_processed_2023-03-06.root
  // The macro for the study is findBackgroundNormalizationScale.C
  // Git hash which produced these numbers is b6d4713c0216f54c5ef5f0632137575923442134
  
  double energyEnergyCorrelatorBackgroundScaleFactors[4][8][8] = {
  // Centrality 0-10%
  {{1.24339,1.16791,1.05333,0.972113,0.910191,0.86347,0.829653,0.808503}, // 120 < jet pT < 140 GeV
   {1.20838,1.13668,1.03133,0.954307,0.900422,0.853821,0.817682,0.797345}, // 140 < jet pT < 160 GeV
   {1.1826,1.11652,1.0222,0.951499,0.902771,0.86042,0.827675,0.802374}, // 160 < jet pT < 180 GeV
   {1.16575,1.10627,1.02144,0.964405,0.929485,0.90218,0.894321,0.890391}, // 180 < jet pT < 200 GeV
   {1.13862,1.08798,1.01654,0.969254,0.938357,0.915893,0.908892,0.907782}, // 200 < jet pT < 300 GeV
   {1.09833,1.06328,1.01886,0.994539,0.98487,0.989774,1.0158,1.05147}, // 300 < jet pT < 500 GeV
   {1.0666,1.04634,1.02522,1.02028,1.03206,1.06274,1.13249,1.21058}, // 500 < jet pT < 5020 GeV
   {1.21194,1.14184,1.03839,0.965045,0.912333,0.870805,0.84204,0.824816}}, // 120 GeV < jet pT
  // Centrality 10-30%
  {{1.18905,1.12023,1.02044,0.944802,0.890309,0.847696,0.818396,0.796878}, // 120 < jet pT < 140 GeV
   {1.16395,1.10174,1.01354,0.94719,0.902885,0.868585,0.841536,0.832731}, // 140 < jet pT < 160 GeV
   {1.14477,1.08868,1.01065,0.954218,0.913196,0.878826,0.862777,0.848241}, // 160 < jet pT < 180 GeV
   {1.12794,1.0769,1.0055,0.952528,0.914929,0.885585,0.865272,0.852716}, // 180 < jet pT < 200 GeV
   {1.10881,1.06702,1.01089,0.974611,0.956623,0.947139,0.954082,0.972533}, // 200 < jet pT < 300 GeV
   {1.07747,1.04928,1.0149,0.998819,1.00158,1.02584,1.07424,1.13769}, // 300 < jet pT < 500 GeV
   {1.05396,1.03879,1.02748,1.03259,1.05841,1.10699,1.19796,1.32016}, // 500 < jet pT < 5020 GeV
   {1.16515,1.10337,1.01569,0.95081,0.906026,0.871724,0.850157,0.838499}}, // 120 GeV < jet pT
  // Centrality 30-50%
  {{1.10536,1.05461,0.981911,0.927105,0.886759,0.854171,0.828703,0.824074}, // 120 < jet pT < 140 GeV
   {1.0922,1.04649,0.98409,0.938928,0.904429,0.875274,0.864164,0.85683}, // 140 < jet pT < 160 GeV
   {1.08467,1.04531,0.992115,0.956685,0.933917,0.919655,0.91831,0.928734}, // 160 < jet pT < 180 GeV
   {1.07648,1.04041,0.993401,0.964823,0.951769,0.944796,0.955022,0.980337}, // 180 < jet pT < 200 GeV
   {1.06958,1.04145,1.00864,0.993383,0.993525,1.01506,1.0569,1.10702}, // 200 < jet pT < 300 GeV
   {1.05704,1.04099,1.02754,1.03506,1.06998,1.13001,1.22063,1.33626}, // 300 < jet pT < 500 GeV
   {1.04685,1.0403,1.04544,1.07418,1.13509,1.24024,1.40037,1.59662}, // 500 < jet pT < 5020 GeV
   {1.09349,1.04909,0.987998,0.944806,0.915401,0.894741,0.886761,0.893141}}, // 120 GeV < jet pT
  // Centrality 50-90%
  {{1.04376,1.01547,0.979617,0.957955,0.950139,0.955351,0.973888,1.00828}, // 120 < jet pT < 140 GeV
   {1.04716,1.02514,1.00282,0.992133,1.00307,1.02457,1.06812,1.15755}, // 140 < jet pT < 160 GeV
   {1.04946,1.03152,1.01368,1.01335,1.03004,1.06758,1.12056,1.1681}, // 160 < jet pT < 180 GeV
   {1.05071,1.03654,1.02793,1.03669,1.0771,1.13267,1.2234,1.32576}, // 180 < jet pT < 200 GeV
   {1.05942,1.05222,1.05893,1.09465,1.16399,1.27751,1.42002,1.59819}, // 200 < jet pT < 300 GeV
   {1.07076,1.07541,1.10752,1.18214,1.3185,1.52103,1.78397,2.09363}, // 300 < jet pT < 500 GeV
   {1.08255,1.09649,1.15331,1.26584,1.47309,1.77735,2.19324,2.65178}, // 500 < jet pT < 5020 GeV
   {1.04791,1.02621,1.00343,0.996898,1.01039,1.04142,1.09119,1.16301}} // 120 GeV < jet pT
  };
  
//  // Alternative numbers for generator level Pythia+Hydjet study
//  // From file: PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_finalMcWeight_processed_2023-03-08.root
//  double energyEnergyCorrelatorBackgroundScaleFactors[4][8][8] = {
//  // Centrality 0-10%
//  {{1.26478,1.20645,1.12693,1.08358,1.07283,1.0808,1.10332,1.1416}, // 120 < jet pT < 140 GeV
//   {1.23447,1.17987,1.10976,1.07261,1.06104,1.07312,1.09625,1.13483}, // 140 < jet pT < 160 GeV
//   {1.21557,1.16781,1.10339,1.07051,1.06236,1.07833,1.10754,1.17682}, // 160 < jet pT < 180 GeV
//   {1.19826,1.15253,1.09747,1.0726,1.07283,1.09409,1.1359,1.18454}, // 180 < jet pT < 200 GeV
//   {1.16849,1.13015,1.08473,1.06827,1.07368,1.09574,1.1443,1.1988}, // 200 < jet pT < 300 GeV
//   {1.12183,1.0954,1.06793,1.06459,1.0826,1.12535,1.19264,1.29374}, // 300 < jet pT < 500 GeV
//   {1.08424,1.06995,1.05983,1.07047,1.10211,1.16079,1.26708,1.39992}, // 500 < jet pT < 5020 GeV
//   {1.23453,1.18156,1.11215,1.07654,1.06916,1.08248,1.11155,1.1581}}, // 120 GeV < jet pT
//  // Centrality 10-30%
//  {{1.22218,1.17152,1.10863,1.07499,1.06832,1.08457,1.11378,1.16138}, // 120 < jet pT < 140 GeV
//   {1.1996,1.15449,1.09938,1.07573,1.07917,1.10123,1.14595,1.21377}, // 140 < jet pT < 160 GeV
//   {1.1808,1.13966,1.09193,1.06946,1.074,1.09954,1.15137,1.21757}, // 160 < jet pT < 180 GeV
//   {1.1649,1.12782,1.08453,1.06886,1.07615,1.10353,1.14912,1.21658}, // 180 < jet pT < 200 GeV
//   {1.14204,1.11183,1.07952,1.07312,1.09286,1.12971,1.1915,1.27956}, // 200 < jet pT < 300 GeV
//   {1.10357,1.08373,1.06701,1.07326,1.1078,1.1742,1.2767,1.40589}, // 300 < jet pT < 500 GeV
//   {1.07442,1.06462,1.06414,1.08577,1.13746,1.21629,1.35567,1.53027}, // 500 < jet pT < 5020 GeV
//   {1.19709,1.15238,1.0986,1.07378,1.07577,1.09889,1.14123,1.20407}}, // 120 GeV < jet pT
//  // Centrality 30-50%
//  {{1.1606,1.1267,1.09032,1.0808,1.09467,1.12997,1.19044,1.26547}, // 120 < jet pT < 140 GeV
//   {1.14461,1.11577,1.08579,1.08336,1.10649,1.15093,1.22422,1.30981}, // 140 < jet pT < 160 GeV
//   {1.13243,1.10735,1.08401,1.08903,1.11675,1.1668,1.24425,1.35124}, // 160 < jet pT < 180 GeV
//   {1.12356,1.10204,1.08285,1.08955,1.11925,1.18031,1.27876,1.40628}, // 180 < jet pT < 200 GeV
//   {1.10868,1.09138,1.08198,1.0981,1.14311,1.22177,1.33707,1.47936}, // 200 < jet pT < 300 GeV
//   {1.08632,1.07824,1.08219,1.11359,1.18168,1.2934,1.44745,1.65153}, // 300 < jet pT < 500 GeV
//   {1.0706,1.07044,1.08974,1.13677,1.224,1.36979,1.59526,1.87562}, // 500 < jet pT < 5020 GeV
//   {1.14364,1.11512,1.08685,1.08577,1.10946,1.15741,1.23448,1.33111}}, // 120 GeV < jet pT
//  // Centrality 50-90%
//  {{1.11818,1.10864,1.11545,1.15384,1.23216,1.33466,1.46324,1.6471}, // 120 < jet pT < 140 GeV
//   {1.11291,1.10781,1.12257,1.16842,1.25304,1.38804,1.57155,1.80188}, // 140 < jet pT < 160 GeV
//   {1.11082,1.10918,1.13172,1.18314,1.28503,1.43661,1.64172,1.88243}, // 160 < jet pT < 180 GeV
//   {1.10771,1.10869,1.13447,1.19721,1.30633,1.48051,1.71875,2.01821}, // 180 < jet pT < 200 GeV
//   {1.10473,1.10922,1.14383,1.22247,1.35651,1.56217,1.84996,2.21101}, // 200 < jet pT < 300 GeV
//   {1.10617,1.12027,1.17405,1.28554,1.47385,1.76142,2.16066,2.6389}, // 300 < jet pT < 500 GeV
//   {1.11091,1.13303,1.21078,1.35645,1.60924,2.00994,2.55657,3.21134}, // 500 < jet pT < 5020 GeV
//   {1.11356,1.10885,1.12503,1.17484,1.26895,1.40669,1.59,1.82841}} // 120 GeV < jet pT
//  };
  
  // Copy the input values to the class array
  for(int iCentrality = 0; iCentrality < kNCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < kNJetPtBins+1; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < kNTrackPtBins; iTrackPt++){
        fBackgroundScale[iCentrality][iJetPt][iTrackPt] = energyEnergyCorrelatorBackgroundScaleFactors[iCentrality][iJetPt][iTrackPt];
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
  
}

/*
 * Check that same bin borders are used in the input EECCard as in the study used the get the numbers for background scale
 *
 *  EECCard* card = Card from which we check the similarity of bin borders to the analysis used to determine the background scales
 *
 *  return: True: bin borders are the same. False: Bin borders are not the same
 */
bool EECBackgroundScale::CheckBinBorders(EECCard* card) const{
  
  // Small number
  const double epsilon = 0.0001;
  
  // Check that centrality bins match
  const int nCentralityBins = card->GetNCentralityBins();
  if(nCentralityBins != kNCentralityBins) return false;
  
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    if(TMath::Abs(kCentralityBinBorders[iCentrality] - card->GetLowBinBorderCentrality(iCentrality)) > epsilon) return false;
    if(TMath::Abs(kCentralityBinBorders[iCentrality+1] - card->GetHighBinBorderCentrality(iCentrality)) > epsilon) return false;
  }
  
  // Check that jet pT bins match
  const int nJetPtBins = card->GetNJetPtBinsEEC();
  if(nJetPtBins != kNJetPtBins) return false;
  
  for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){
    if(TMath::Abs(kJetPtBinBorders[iJetPt] - card->GetLowBinBorderJetPtEEC(iJetPt)) > epsilon) return false;
    if(TMath::Abs(kJetPtBinBorders[iJetPt+1] - card->GetHighBinBorderJetPtEEC(iJetPt)) > epsilon) return false;
  }
  
  // Check that track pT bins match
  const int nTrackPtBins = card->GetNTrackPtBinsEEC();
  if(nTrackPtBins != kNTrackPtBins) return false;
  
  for(int iTrackPt = 0; iTrackPt < nTrackPtBins; iTrackPt++){
    if(TMath::Abs(kTrackPtBinBorders[iTrackPt] - card->GetLowBinBorderTrackPtEEC(iTrackPt)) > epsilon) return false;
    if(TMath::Abs(kTrackPtBinBorders[iTrackPt+1] - card->GetHighBinBorderTrackPtEEC(iTrackPt)) > epsilon) return false;
  }
  
  // If all checks pass, return true
  return true;;
}

/*
 * Get the background scaling factor for the given bin
 *
 *  const int iCentrality = Centrality bin for scaling factor
 *  const int iJetPt = Jet pT bin for the scaling factor
 *  const int iTrackPt = Track pT bin for the scaling factor
 *
 *  return: Scaling factor corresponding to the defined bin
 */
double EECBackgroundScale::GetEECBackgroundScale(const int iCentrality, const int iJetPt, const int iTrackPt) const{
  
  // Find the scaling factor in the bin corresponding to the given bin
  return fBackgroundScale[iCentrality][iJetPt][iTrackPt];
  
}
