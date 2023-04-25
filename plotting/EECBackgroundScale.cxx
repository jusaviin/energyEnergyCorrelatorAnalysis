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
 * Custom constructor
 */
EECBackgroundScale::EECBackgroundScale(EECCard* card){
  TString dataType = card->GetDataType();
  bool useGenJets = (dataType.Contains("GenReco") || dataType.Contains("GenGen"));
  InitializeArrays(useGenJets);
}

/*
 * Destructor
 */
EECBackgroundScale::~EECBackgroundScale(){

}

/*
 * Initialization for background scale arrays
 */
void EECBackgroundScale::InitializeArrays(const bool useGenJets){
  
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
  // When doing this matching, we can trust that the Hydjet particles are mostly background, and Pythia particles mostly signal.
  // We determine the true background from Pythia+Hydjet and Hydjet+Hydjet pairings from the signal cone.
  // The reflected cone background estimate is obtained by combining all particles from signal cone with all particles from the reflected cone.
  // We know that in this case fake+fake contribution is not properly modeled, but this contribution is negligible in the region
  // Track pT > 2 GeV and DeltaR < 0.4. We also know that there is a bias in reconstructing jets on upwards fluctuations of background.
  // To take the fluctuations into account and to mitigate the mismodeling of the fake+fake part, we look at the integrals of
  // true background and the reflected cone background estimate in the region DeltaR < 0.4 and calculate their ratio.
  // This ratio now reflects excess fluctuations below the jet peak in Monte Carlo. Assuming the ratio is similar in data,
  // the reflected cone estimate from data can be scaled with this number to get a good estimate of the absolute background level.
  //
  // The file that was used for the study is PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_finalMcWeight_matchJets_fixCentrality_processed_2023-03-06.root
  // The macro for the study is findBackgroundNormalizationScale.C
  // Git hash which produced these numbers is b6d4713c0216f54c5ef5f0632137575923442134
  
  double energyEnergyCorrelatorBackgroundScaleFactors[2][4][8][8] = 
  {
  // The first part of the array gives the results for reconstructed jets.
  // The file that was used for the study is PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_finalMcWeight_matchJets_fixCentrality_processed_2023-03-06.root
  // The macro for the study is findBackgroundNormalizationScale.C
  // Git hash which produced these numbers is b6d4713c0216f54c5ef5f0632137575923442134
  {
  // Centrality 4-14%
  {{1.18905,1.12023,1.02044,0.944802,0.890309,0.847696,0.818396,0.796878}, // 120 < jet pT < 140 GeV
   {1.16395,1.10174,1.01354,0.94719,0.902885,0.868585,0.841536,0.832731}, // 140 < jet pT < 160 GeV
   {1.14477,1.08868,1.01065,0.954218,0.913196,0.878826,0.862777,0.848241}, // 160 < jet pT < 180 GeV
   {1.12794,1.0769,1.0055,0.952528,0.914929,0.885585,0.865272,0.852716}, // 180 < jet pT < 200 GeV
   {1.10881,1.06702,1.01089,0.974611,0.956623,0.947139,0.954082,0.972533}, // 200 < jet pT < 300 GeV
   {1.07747,1.04928,1.0149,0.998819,1.00158,1.02584,1.07424,1.13769}, // 300 < jet pT < 500 GeV
   {1.05396,1.03879,1.02748,1.03259,1.05841,1.10699,1.19796,1.32016}, // 500 < jet pT < 5020 GeV
   {1.16515,1.10337,1.01569,0.95081,0.906026,0.871724,0.850157,0.838499}}, // 120 GeV < jet pT
   // Centrality 14-34%
  {{1.10536,1.05461,0.981911,0.927105,0.886759,0.854171,0.828703,0.824074}, // 120 < jet pT < 140 GeV
   {1.0922,1.04649,0.98409,0.938928,0.904429,0.875274,0.864164,0.85683}, // 140 < jet pT < 160 GeV
   {1.08467,1.04531,0.992115,0.956685,0.933917,0.919655,0.91831,0.928734}, // 160 < jet pT < 180 GeV
   {1.07648,1.04041,0.993401,0.964823,0.951769,0.944796,0.955022,0.980337}, // 180 < jet pT < 200 GeV
   {1.06958,1.04145,1.00864,0.993383,0.993525,1.01506,1.0569,1.10702}, // 200 < jet pT < 300 GeV
   {1.05704,1.04099,1.02754,1.03506,1.06998,1.13001,1.22063,1.33626}, // 300 < jet pT < 500 GeV
   {1.04685,1.0403,1.04544,1.07418,1.13509,1.24024,1.40037,1.59662}, // 500 < jet pT < 5020 GeV
   {1.09349,1.04909,0.987998,0.944806,0.915401,0.894741,0.886761,0.893141}}, // 120 GeV < jet pT
   // Centrality 34-54%
  {{1.04376,1.01547,0.979617,0.957955,0.950139,0.955351,0.973888,1.00828}, // 120 < jet pT < 140 GeV
   {1.04716,1.02514,1.00282,0.992133,1.00307,1.02457,1.06812,1.15755}, // 140 < jet pT < 160 GeV
   {1.04946,1.03152,1.01368,1.01335,1.03004,1.06758,1.12056,1.1681}, // 160 < jet pT < 180 GeV
   {1.05071,1.03654,1.02793,1.03669,1.0771,1.13267,1.2234,1.32576}, // 180 < jet pT < 200 GeV
   {1.05942,1.05222,1.05893,1.09465,1.16399,1.27751,1.42002,1.59819}, // 200 < jet pT < 300 GeV
   {1.07076,1.07541,1.10752,1.18214,1.3185,1.52103,1.78397,2.09363}, // 300 < jet pT < 500 GeV
   {1.08255,1.09649,1.15331,1.26584,1.47309,1.77735,2.19324,2.65178}, // 500 < jet pT < 5020 GeV
   {1.04791,1.02621,1.00343,0.996898,1.01039,1.04142,1.09119,1.16301}}, // 120 GeV < jet pT
   // Centrality 54-94%
  {{1.13556,1.14333,1.18819,1.28948,1.43372,1.64952,1.87327,2.12782}, // 120 < jet pT < 140 GeV
   {1.15545,1.16939,1.2346,1.36201,1.56897,1.84668,2.13731,2.42631}, // 140 < jet pT < 160 GeV
   {1.17116,1.19434,1.27858,1.45103,1.73068,2.10461,2.55429,2.94246}, // 160 < jet pT < 180 GeV
   {1.18055,1.2072,1.30953,1.47249,1.7461,2.10928,2.47767,2.82787}, // 180 < jet pT < 200 GeV
   {1.21108,1.25163,1.38978,1.65001,2.05774,2.61467,3.31257,4.05708}, // 200 < jet pT < 300 GeV
   {1.2636,1.32739,1.52738,1.89995,2.51681,3.43076,4.52688,5.84696}, // 300 < jet pT < 500 GeV
   {1.3183,1.40498,1.69352,2.21939,3.08818,4.38779,5.97017,7.95629}, // 500 < jet pT < 5020 GeV
   {1.15858,1.1757,1.24724,1.38973,1.60978,1.91682,2.2532,2.60344}} // 120 GeV < jet pT
  },
  // The second part of the array gives the same numbers for a generator level jet study
  // From file: PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_finalMcWeight_fixCentrality_processed_2023-03-08.root
  {
  // Centrality 0-10%
  {{1.22218,1.17152,1.10863,1.07499,1.06832,1.08457,1.11378,1.16138}, // 120 < jet pT < 140 GeV
   {1.1996,1.15449,1.09938,1.07573,1.07917,1.10123,1.14595,1.21377}, // 140 < jet pT < 160 GeV
   {1.1808,1.13966,1.09193,1.06946,1.074,1.09954,1.15137,1.21757}, // 160 < jet pT < 180 GeV
   {1.1649,1.12782,1.08453,1.06886,1.07615,1.10353,1.14912,1.21658}, // 180 < jet pT < 200 GeV
   {1.14204,1.11183,1.07952,1.07312,1.09286,1.12971,1.1915,1.27956}, // 200 < jet pT < 300 GeV
   {1.10357,1.08373,1.06701,1.07326,1.1078,1.1742,1.2767,1.40589}, // 300 < jet pT < 500 GeV
   {1.07442,1.06462,1.06414,1.08577,1.13746,1.21629,1.35567,1.53027}, // 500 < jet pT < 5020 GeV
   {1.19709,1.15238,1.0986,1.07378,1.07577,1.09889,1.14123,1.20407}}, // 120 GeV < jet pT
  // Centrality 10-30%
  {{1.1606,1.1267,1.09032,1.0808,1.09467,1.12997,1.19044,1.26547}, // 120 < jet pT < 140 GeV
   {1.14461,1.11577,1.08579,1.08336,1.10649,1.15093,1.22422,1.30981}, // 140 < jet pT < 160 GeV
   {1.13243,1.10735,1.08401,1.08903,1.11675,1.1668,1.24425,1.35124}, // 160 < jet pT < 180 GeV
   {1.12356,1.10204,1.08285,1.08955,1.11925,1.18031,1.27876,1.40628}, // 180 < jet pT < 200 GeV
   {1.10868,1.09138,1.08198,1.0981,1.14311,1.22177,1.33707,1.47936}, // 200 < jet pT < 300 GeV
   {1.08632,1.07824,1.08219,1.11359,1.18168,1.2934,1.44745,1.65153}, // 300 < jet pT < 500 GeV
   {1.0706,1.07044,1.08974,1.13677,1.224,1.36979,1.59526,1.87562}, // 500 < jet pT < 5020 GeV
   {1.14364,1.11512,1.08685,1.08577,1.10946,1.15741,1.23448,1.33111}}, // 120 GeV < jet pT
   // Centrality 30-50%
  {{1.11818,1.10864,1.11545,1.15384,1.23216,1.33466,1.46324,1.6471}, // 120 < jet pT < 140 GeV
   {1.11291,1.10781,1.12257,1.16842,1.25304,1.38804,1.57155,1.80188}, // 140 < jet pT < 160 GeV
   {1.11082,1.10918,1.13172,1.18314,1.28503,1.43661,1.64172,1.88243}, // 160 < jet pT < 180 GeV
   {1.10771,1.10869,1.13447,1.19721,1.30633,1.48051,1.71875,2.01821}, // 180 < jet pT < 200 GeV
   {1.10473,1.10922,1.14383,1.22247,1.35651,1.56217,1.84996,2.21101}, // 200 < jet pT < 300 GeV
   {1.10617,1.12027,1.17405,1.28554,1.47385,1.76142,2.16066,2.6389}, // 300 < jet pT < 500 GeV
   {1.11091,1.13303,1.21078,1.35645,1.60924,2.00994,2.55657,3.21134}, // 500 < jet pT < 5020 GeV
   {1.11356,1.10885,1.12503,1.17484,1.26895,1.40669,1.59,1.82841}}, // 120 GeV < jet pT
   // Centrality 50-90%
  {{1.22544,1.25723,1.36646,1.5672,1.90046,2.39018,2.99232,3.68625}, // 120 < jet pT < 140 GeV
   {1.23418,1.27066,1.39823,1.6352,2.0295,2.5826,3.28602,4.11905}, // 140 < jet pT < 160 GeV
   {1.24375,1.2887,1.42981,1.68919,2.12921,2.80657,3.73974,4.73199}, // 160 < jet pT < 180 GeV
   {1.25251,1.29934,1.45626,1.74812,2.22684,2.87438,3.70034,4.5799}, // 180 < jet pT < 200 GeV
   {1.26656,1.32405,1.50711,1.85368,2.42182,3.26258,4.38555,5.71303}, // 200 < jet pT < 300 GeV
   {1.30952,1.3866,1.63157,2.09052,2.85941,4.07983,5.64174,7.57191}, // 300 < jet pT < 500 GeV
   {1.35841,1.46074,1.78526,2.38641,3.40879,4.93738,6.97275,9.48935}, // 500 < jet pT < 5020 GeV
   {1.23819,1.27774,1.41011,1.65563,2.06399,2.66084,3.4269,4.30642}} // 120 GeV < jet pT
  }};
  
  // Copy the input values to the class array
  for(int iCentrality = 0; iCentrality < kNCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < kNJetPtBins+1; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < kNTrackPtBins; iTrackPt++){
        fBackgroundScale[iCentrality][iJetPt][iTrackPt] = energyEnergyCorrelatorBackgroundScaleFactors[useGenJets][iCentrality][iJetPt][iTrackPt];
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
