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
