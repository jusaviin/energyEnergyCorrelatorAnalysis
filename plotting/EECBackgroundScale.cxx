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
  
  // Warning: the method has been checked only for 120 < jet pT < 140 GeV and 0-10 % centrality bin!
  // If it works in that bin, all the other bins will need to be rechecked.
  //
  // In the current method, we start with Pythia+Hydjet simulation from which the Hydjet initiating jets are removed by track pT cut.
  // In this simulation, we can trust that the Hydjet particles are mostly background, and Pythia particles mostly signal.
  // We determine the true background from Pythia+Hydjet and Hydjet+Hydjet pairings from the signal cone.
  // The reflected cone background estimate is obtained by combining all particles from signal cone with all particles from the reflected cone.
  // We know that in this case fake+fake contribution is not properly modeled, but this contribution is negligible in the region
  // Track pT > 2 GeV and DeltaR < 0.4. We also know that there is a bias in reconstruction jets on upwards fluctuations of background.
  // To take the fluctuations into account and to mitigate the mismodeling of the fake+fake part, we look at the integrals of
  // true background and the reflected cone background estimate in the region DeltaR < 0.4 and calculate theur ration.
  // This ratio now reflects excess fluctuations below the jet peak in Monte Carlo. Assuming the ratio is similar in data,
  // the reflected cone estimate from data can be scaled with this number to get a good estimate of the absolute background level,
  //
  // The file that was used for the study is data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_4pC_newSRC_cutBadPhiAndComb_wtaAxis_jetTrigger_preprocessed_2022-11-29.root
  // The macro for the study is findBackgroundNormalizationScale.C
  // Git hash which produced these numbers is f4918db708d8f3ed0ee9b5d8a4cdd40faee27508
  
  double energyEnergyCorrelatorBackgroundScaleFactors[4][8][8] = {
    // Centrality 0-10%
    {{1.23704,1.16744,1.06653,0.993784,0.943877,0.903095,0.869518,0.850716}, // 120 < jet pT < 140 GeV
      {1.20002,1.1335,1.03594,0.968634,0.918617,0.872191,0.845861,0.821166}, // 140 < jet pT < 160 GeV
      {1.17488,1.11325,1.02579,0.959541,0.911736,0.873274,0.852182,0.843595}, // 160 < jet pT < 180 GeV
      {1.15749,1.10005,1.01791,0.96239,0.91762,0.883312,0.86252,0.857055}, // 180 < jet pT < 200 GeV
      {1.13615,1.08882,1.0242,0.982212,0.95755,0.944951,0.943092,0.952309}, // 200 < jet pT < 300 GeV
      {1.09845,1.06539,1.02444,1.00358,1.00068,1.01533,1.0463,1.09176}, // 300 < jet pT < 500 GeV
      {1.07444,1.05303,1.0325,1.03943,1.04623,1.11102,1.19921,1.31657}, // 500 < jet pT < 5020 GeV
      {1.20208,1.13756,1.04494,0.979664,0.933861,0.896488,0.871957,0.858326}}, // 120 GeV < jet pT
    // Centrality 10-30%
    {{1.18651,1.12492,1.03667,0.971424,0.922675,0.887954,0.864709,0.853848}, // 120 < jet pT < 140 GeV
      {1.15479,1.09755,1.01515,0.952968,0.9137,0.876982,0.850657,0.832422}, // 140 < jet pT < 160 GeV
      {1.13687,1.08445,1.01214,0.958919,0.922108,0.888424,0.865302,0.846724}, // 160 < jet pT < 180 GeV
      {1.12513,1.0789,1.0145,0.971409,0.940166,0.9222,0.915541,0.911173}, // 180 < jet pT < 200 GeV
      {1.1061,1.06667,1.01514,0.982473,0.969136,0.962757,0.972652,0.995328}, // 200 < jet pT < 300 GeV
      {1.07685,1.05047,1.01974,1.00537,1.01129,1.03812,1.09317,1.16177}, // 300 < jet pT < 500 GeV
      {1.05568,1.0405,1.02954,1.036,1.06478,1.11214,1.20631,1.32478}, // 500 < jet pT < 5020 GeV
      {1.15828,1.1026,1.02413,0.967128,0.928308,0.898593,0.880683,0.872278}}, // 120 GeV < jet pT
    // Centrality 30-50%
    {{1.11296,1.06835,1.00664,0.963227,0.935016,0.914957,0.911703,0.920274}, // 120 < jet pT < 140 GeV
      {1.09092,1.04965,0.993018,0.952635,0.923766,0.903556,0.888267,0.885671}, // 140 < jet pT < 160 GeV
      {1.08221,1.04535,0.998229,0.968978,0.95106,0.944834,0.957095,0.968511}, // 160 < jet pT < 180 GeV
      {1.07853,1.0462,1.00657,0.983988,0.978828,0.987987,1.01681,1.05886}, // 180 < jet pT < 200 GeV
      {1.06965,1.04393,1.01481,1.00481,1.01121,1.03958,1.08776,1.14942}, // 200 < jet pT < 300 GeV
      {1.05711,1.04225,1.03113,1.04062,1.07619,1.13826,1.22899,1.34887}, // 300 < jet pT < 500 GeV
      {1.04911,1.04336,1.04793,1.07687,1.14002,1.25146,1.42424,1.632}, // 500 < jet pT < 5020 GeV
      {1.09571,1.05624,1.00362,0.968921,0.948548,0.939112,0.944441,0.960664}}, // 120 GeV < jet pT
    // Centrality 50-90%
    {{1.06068,1.03895,1.01718,1.01362,1.02771,1.05936,1.10195,1.17502}, // 120 < jet pT < 140 GeV
      {1.05475,1.03553,1.02042,1.02025,1.04096,1.07915,1.14348,1.24066}, // 140 < jet pT < 160 GeV
      {1.05552,1.0406,1.03252,1.04166,1.0723,1.12766,1.19662,1.27664}, // 160 < jet pT < 180 GeV
      {1.05843,1.04716,1.04363,1.06268,1.11031,1.19574,1.30304,1.4012}, // 180 < jet pT < 200 GeV
      {1.06181,1.05643,1.0661,1.10465,1.17753,1.29089,1.43808,1.619}, // 200 < jet pT < 300 GeV
      {1.0738,1.07934,1.11558,1.19524,1.33441,1.55237,1.83499,2.16146}, // 300 < jet pT < 500 GeV
      {1.08484,1.09889,1.15445,1.26709,1.47662,1.77611,2.19168,2.67354}, // 500 < jet pT < 5020 GeV
      {1.05885,1.04166,1.02914,1.03592,1.06504,1.11763,1.1893,1.28798}} // 120 GeV < jet pT
  };
  
//  // Alternative numbers for generator level Pythia+Hydjet study
//  // From file: PbPbMC2018_GenGen_eecAnalysis_akFlowJet_MnD_wtaAxis_noTrigger_preprocessed_2022-10-21.root
//  double energyEnergyCorrelatorBackgroundScaleFactors[4][8][8] = {
//    // Centrality 0-10%
//    {{1.24833,1.19255,1.11839,1.07752,1.06599,1.07518,1.09916,1.1413}, // 120 < jet pT < 140 GeV
//      {1.22193,1.17115,1.10711,1.0747,1.06981,1.08629,1.12001,1.17128}, // 140 < jet pT < 160 GeV
//      {1.20211,1.15629,1.09901,1.07061,1.06594,1.08314,1.12375,1.18547}, // 160 < jet pT < 180 GeV
//      {1.18522,1.14325,1.0922,1.07037,1.07329,1.09663,1.13678,1.19436}, // 180 < jet pT < 200 GeV
//      {1.15816,1.12301,1.0829,1.07049,1.08183,1.11091,1.16574,1.23718}, // 200 < jet pT < 300 GeV
//      {1.11491,1.09137,1.06849,1.06879,1.09368,1.1458,1.22535,1.33188}, // 300 < jet pT < 500 GeV
//      {1.08126,1.06803,1.06151,1.07598,1.1132,1.17707,1.29048,1.43664}, // 500 < jet pT < 5020 GeV
//      {1.22094,1.17082,1.10688,1.07459,1.06968,1.08567,1.11982,1.17249}}, // 120 GeV < jet pT
//    // Centrality 10-30%
//    {{1.17683,1.13794,1.09411,1.07634,1.08276,1.11182,1.1615,1.22112}, // 120 < jet pT < 140 GeV
//      {1.15876,1.12516,1.08869,1.07953,1.09607,1.1305,1.19191,1.2789}, // 140 < jet pT < 160 GeV
//      {1.14402,1.11443,1.08405,1.07977,1.0985,1.14326,1.21211,1.30134}, // 160 < jet pT < 180 GeV
//      {1.13412,1.10817,1.08351,1.08477,1.1101,1.16207,1.24468,1.35558}, // 180 < jet pT < 200 GeV
//      {1.11676,1.09576,1.07912,1.08884,1.12515,1.18841,1.28615,1.40962}, // 200 < jet pT < 300 GeV
//      {1.08986,1.07807,1.07604,1.09961,1.15619,1.2531,1.39579,1.58212}, // 300 < jet pT < 500 GeV
//      {1.07063,1.06725,1.08007,1.11892,1.19561,1.32337,1.525,1.78623}, // 500 < jet pT < 5020 GeV
//      {1.15784,1.12445,1.08891,1.07995,1.09592,1.13511,1.19973,1.28313}}, // 120 GeV < jet pT
//    // Centrality 30-50%
//    {{1.11932,1.10515,1.1031,1.13082,1.19136,1.27719,1.38949,1.54421}, // 120 < jet pT < 140 GeV
//      {1.11268,1.10305,1.10715,1.14148,1.20602,1.31254,1.46398,1.63687}, // 140 < jet pT < 160 GeV
//      {1.10957,1.10369,1.11525,1.15474,1.23319,1.35798,1.52638,1.73476}, // 160 < jet pT < 180 GeV
//      {1.10536,1.1014,1.11537,1.15947,1.24545,1.38424,1.57074,1.82358}, // 180 < jet pT < 200 GeV
//      {1.09889,1.09847,1.12217,1.18322,1.29174,1.45765,1.68554,1.97383}, // 200 < jet pT < 300 GeV
//      {1.09691,1.10611,1.14776,1.23926,1.39722,1.64241,1.97301,2.38533}, // 300 < jet pT < 500 GeV
//      {1.09733,1.11402,1.17259,1.28612,1.49102,1.81126,2.26866,2.83169}, // 500 < jet pT < 5020 GeV
//      {1.11301,1.10353,1.10939,1.14627,1.21926,1.33092,1.48239,1.67653}}, // 120 GeV < jet pT
//    // Centrality 50-90%
//    {{1.18122,1.20179,1.28279,1.43,1.67356,2.02258,2.43578,2.9503}, // 120 < jet pT < 140 GeV
//      {1.18733,1.21203,1.30333,1.48007,1.76714,2.18903,2.69846,3.35842}, // 140 < jet pT < 160 GeV
//      {1.18934,1.21873,1.31781,1.50644,1.82099,2.25648,2.84968,3.49}, // 160 < jet pT < 180 GeV
//      {1.19505,1.22692,1.33757,1.54575,1.88866,2.37181,3.06116,3.79745}, // 180 < jet pT < 200 GeV
//      {1.20777,1.24789,1.38147,1.6289,2.04402,2.658,3.47599,4.4359}, // 200 < jet pT < 300 GeV
//      {1.23583,1.29057,1.46626,1.79607,2.358,3.20339,4.29347,5.68616}, // 300 < jet pT < 500 GeV
//      {1.27183,1.34855,1.58075,2.01957,2.75531,3.86946,5.39865,7.20834}, // 500 < jet pT < 5020 GeV
//      {1.18861,1.21488,1.3104,1.4889,1.78436,2.21105,2.74798,3.3972}} // 120 GeV < jet pT
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
