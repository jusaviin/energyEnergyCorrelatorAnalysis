/*
 * Implementation of the ReflectedConeWeight class
 */

// Own includes
#include "ReflectedConeWeight.h"

/*
 * Contructor
 */
ReflectedConeWeight::ReflectedConeWeight()
{
  InitializeArrays();
}

/*
 * Destructor
 */
ReflectedConeWeight::~ReflectedConeWeight(){

}

/*
 * Initialization for seagull method arrays. Set all bins that need nonzero value
 */
void ReflectedConeWeight::InitializeArrays(){
  
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
  // What is currently done, is that the signal particle density is divided into two regions, that are fitted
  // with a line. The regions are the "kink region", 0.34 < DeltaR < 0.44, and the "background region", 0.44 < DeltaR < 0.6
  // We assume that the particle density for the background throughout will follow the shape of the background region
  // in everywhere else except in the kink region, where it follows the shape of the kink region. We have also
  // fitted the particle density shape in the reflected cone within 0.05 < DeltaR < 0.5.
  //
  // To apply the weight, the reflected cone basic shape is taken into account, and additional weight is added on
  // top of that to make the reflected cone shape match with the background shape in the signal cone.
  // BGD = Background density. KD = Kink density. RCD = Reflected cone density
  //
  // The file that was used for the study is data/eecAnalysis_akFlowJets_removeBadAcceptance_wtaAxis_processed_2022-10-25.root
  // The macro for the study is particleDensityFitter.C
  
  // =======================================================================
  // == For signal cone, this is the BGD(DeltaR = 0) / BGD(DeltaR = 0.34) ==
  // =======================================================================
  
  double particleDensityPtBinned[4][7][8] = {
    // Centrality 0-10%
    {{1.03715,1.06107,1.09448,1.11818,1.16904,1.24356,1.30797,1.47422}, // 120 < jet pT < 140 GeV
      {1.03485,1.06486,1.09475,1.11883,1.17582,1.25746,1.33473,1.48767}, // 140 < jet pT < 160 GeV
      {1.0548,1.07493,1.07856,1.15852,1.17005,1.21852,1.27118,1.53213}, // 160 < jet pT < 180 GeV
      {1.04032,1.05537,1.09292,1.13319,1.11144,1.27421,1.25746,1.53812}, // 180 < jet pT < 200 GeV
      {1.04553,1.06528,1.08965,1.10083,1.13599,1.25247,1.24137,1.57262}, // 200 < jet pT < 300 GeV
      {1.05817,1.07609,1.06999,1.15767,1.15764,1.07978,1.31419,1.63778}, // 300 < jet pT < 500 GeV
      {1.13373,1.11026,1.12378,0.637375,0.980963,1.03702,1.92669,14.4891}}, // 500 < jet pT < 5020 GeV
    // Centrality 10-30%
    {{1.05156,1.08317,1.11909,1.15237,1.22908,1.25442,1.31647,1.55014}, // 120 < jet pT < 140 GeV
      {1.05598,1.08563,1.12896,1.16203,1.21976,1.27312,1.39177,1.63026}, // 140 < jet pT < 160 GeV
      {1.04015,1.07901,1.13128,1.17067,1.19572,1.32684,1.3499,1.5936}, // 160 < jet pT < 180 GeV
      {1.04827,1.0784,1.1361,1.17013,1.2098,1.23272,1.38883,1.66619}, // 180 < jet pT < 200 GeV
      {1.0325,1.08305,1.12816,1.16583,1.22409,1.35243,1.45061,1.65732}, // 200 < jet pT < 300 GeV
      {1.0743,1.08202,1.10218,1.08297,1.01662,1.37998,1.53247,1.75862}, // 300 < jet pT < 500 GeV
      {1.00842,1.28426,1.1829,1.50521,-0.115873,0.557432,1.66047,-0.681321}}, // 500 < jet pT < 5020 GeV
    // Centrality 30-50%
    {{1.06992,1.1141,1.17996,1.22621,1.28255,1.37186,1.54424,1.65992}, // 120 < jet pT < 140 GeV
      {1.09181,1.13607,1.21016,1.25016,1.33593,1.37461,1.55997,1.74541}, // 140 < jet pT < 160 GeV
      {1.06294,1.13414,1.18817,1.26763,1.37175,1.48958,1.52663,1.7646}, // 160 < jet pT < 180 GeV
      {1.02665,1.12624,1.15702,1.13182,1.49083,1.54615,1.59895,1.73059}, // 180 < jet pT < 200 GeV
      {1.11836,1.11096,1.21399,1.20476,1.29088,1.36134,1.39283,1.70273}, // 200 < jet pT < 300 GeV
      {1.13478,1.13885,1.1637,1.05514,1.5652,1.70816,1.77322,1.46233}, // 300 < jet pT < 500 GeV
      {0.822128,0.959372,1.84294,1.2499,1.95462,1.27786,-0.0557169,1.72203}}, // 500 < jet pT < 5020 GeV
    // Centrality 50-90%
    {{1.17063,1.23714,1.37509,1.48029,1.62406,1.58877,1.66833,1.85168}, // 120 < jet pT < 140 GeV
      {1.21553,1.28124,1.31239,1.4633,1.71336,1.72282,1.76285,1.87719}, // 140 < jet pT < 160 GeV
      {1.19425,1.23396,1.44234,1.54604,1.58989,1.29237,1.63028,1.95419}, // 160 < jet pT < 180 GeV
      {1.11314,1.24665,1.48834,1.65943,1.60079,1.86546,1.82767,2.04839}, // 180 < jet pT < 200 GeV
      {1.00121,1.21642,1.47881,1.43536,1.63387,1.62927,1.65165,1.82287}, // 200 < jet pT < 300 GeV
      {1.31225,1.32221,1.08278,1.11992,1.19511,0.338351,8.13283,1.37021}, // 300 < jet pT < 500 GeV
      {1.81099,1.07239,1.65551,2.0829,1.32638,1.45039,1.26522,2.07134}} // 500 < jet pT < 5020 GeV
  };
  
  // ==========================================================================
  // == For reflected cone, this is the RCD(DeltaR = 0) / RCD(DeltaR = 0.34) ==
  // ==========================================================================
  
  double particleDensityPtBinnedReflectedCone[4][7][8] = {
    // Centrality 0-10%
    {{1.02082,1.00862,1.01642,1.01515,1.01432,1.01606,1.00175,0.977148}, // 120 < jet pT < 140 GeV
      {1.01865,1.01014,1.01084,1.02183,1.01132,1.01743,1.01494,0.954175}, // 140 < jet pT < 160 GeV
      {1.02264,1.00819,1.00817,1.02303,1.01504,0.993528,1.0015,0.93633}, // 160 < jet pT < 180 GeV
      {1.01334,1.00689,1.00169,1.00844,1.00375,1.00659,0.924959,0.926645}, // 180 < jet pT < 200 GeV
      {1.01107,1.00693,1.00315,0.997835,1.02093,0.974936,0.995928,0.947067}, // 200 < jet pT < 300 GeV
      {1.01763,0.994131,0.979124,1.03165,1.01044,0.920156,0.889442,0.975322}, // 300 < jet pT < 500 GeV
      {1.03154,0.988376,0.891902,0.931592,0.785267,1.08857,1.79411,1.05325}}, // 500 < jet pT < 5020 GeV
    // Centrality 10-30%
    {{1.0238,1.01338,1.02008,1.01754,1.01021,1.01454,0.996103,0.948179}, // 120 < jet pT < 140 GeV
      {1.01828,1.01047,1.01417,1.02218,1.01015,0.998293,1.01471,0.902489}, // 140 < jet pT < 160 GeV
      {1.02516,1.01023,1.01009,1.0156,1.01338,1.00171,0.980917,0.904686}, // 160 < jet pT < 180 GeV
      {1.0227,1.01042,1.01133,1.03029,0.984111,1.0075,0.929582,0.864887}, // 180 < jet pT < 200 GeV
      {1.01775,1.01437,0.996733,1.01602,1.00932,1.01201,0.978764,0.839422}, // 200 < jet pT < 300 GeV
      {1.03797,1.00841,1.01821,1.03519,0.97356,0.942085,0.951721,0.784917}, // 300 < jet pT < 500 GeV
      {0.977463,1.00542,1.02743,1.00869,1.23815,1.50099,1.43687,1.68672}}, // 500 < jet pT < 5020 GeV
    // Centrality 30-50%
    {{1.01459,1.00239,1.00455,0.981497,0.989248,0.967323,0.936872,0.868135}, // 120 < jet pT < 140 GeV
      {1.03247,1.01844,0.99048,0.997131,0.962352,0.959492,0.992209,0.837536}, // 140 < jet pT < 160 GeV
      {1.02043,0.998105,0.99563,0.968784,1.02864,0.953518,0.905975,0.718988}, // 160 < jet pT < 180 GeV
      {1.01345,0.998425,0.997992,0.990203,0.958247,0.79677,0.85976,0.688165}, // 180 < jet pT < 200 GeV
      {0.992513,0.998578,0.985573,0.977705,0.972259,0.892676,0.84765,0.822296}, // 200 < jet pT < 300 GeV
      {0.951206,1.00538,0.990178,0.902898,1.01015,0.750032,0.280739,0.563559}, // 300 < jet pT < 500 GeV
      {0.891639,0.758953,0.883179,0.928523,1.85112,1.83924,2.12026,1.72923}}, // 500 < jet pT < 5020 GeV
    // Centrality 50-90%
    {{0.985088,1.00918,0.985538,0.940368,0.908459,0.767333,0.816034,0.619181}, // 120 < jet pT < 140 GeV
      {0.993021,0.980018,0.954558,0.895924,0.880789,0.681584,0.979857,0.578747}, // 140 < jet pT < 160 GeV
      {0.979499,1.02443,0.917511,0.953081,0.668039,0.749618,0.568544,0.650549}, // 160 < jet pT < 180 GeV
      {1.08008,0.959395,0.952771,0.942338,0.555029,0.499021,1.18723,0.486642}, // 180 < jet pT < 200 GeV
      {0.962675,0.915826,0.896399,0.965635,0.575113,0.956269,0.37644,0.283039}, // 200 < jet pT < 300 GeV
      {0.84009,0.842098,0.921152,1.09828,0.998449,1.29171,1.17712,1.10162}, // 300 < jet pT < 500 GeV
      {2.46241,1.73441,2.00884,2.35006,1.83857,2.48267,-1,5.1498}} // 500 < jet pT < 5020 GeV
  };
  
  // =======================================================================
  // == For signal cone, this is the KD(DeltaR = 0.34) / KD(DeltaR = 0.4) ==
  // =======================================================================
  double particleDensityPtBinnedKink[4][7][8] = {
    // Centrality 0-10%
    {{1.01956,1.03241,1.04839,1.06856,1.08588,1.11997,1.17462,1.28525}, // 120 < jet pT < 140 GeV
      {1.02127,1.03596,1.05271,1.07004,1.11401,1.15368,1.16052,1.34573}, // 140 < jet pT < 160 GeV
      {1.02387,1.03685,1.05011,1.07976,1.10573,1.14527,1.14068,1.34169}, // 160 < jet pT < 180 GeV
      {1.02582,1.02923,1.04774,1.07072,1.10056,1.14352,1.19858,1.40839}, // 180 < jet pT < 200 GeV
      {1.01785,1.03732,1.04849,1.07056,1.11744,1.14266,1.1558,1.37441}, // 200 < jet pT < 300 GeV
      {1.00627,1.03256,1.05228,1.05611,1.08365,1.23619,1.14564,1.46196}, // 300 < jet pT < 500 GeV
      {0.959969,0.982598,1.04342,1.30286,1.00475,1.17881,0.876961,1.83153}}, // 500 < jet pT < 5020 GeV
    // Centrality 10-30%
    {{1.02395,1.03951,1.06131,1.09065,1.12975,1.16352,1.21913,1.35333}, // 120 < jet pT < 140 GeV
      {1.02884,1.0405,1.07048,1.10318,1.12187,1.17087,1.22957,1.38244}, // 140 < jet pT < 160 GeV
      {1.02708,1.04445,1.07056,1.0936,1.12721,1.17368,1.23968,1.41247}, // 160 < jet pT < 180 GeV
      {1.01572,1.04758,1.0712,1.09351,1.13567,1.17539,1.21135,1.42362}, // 180 < jet pT < 200 GeV
      {1.03024,1.03774,1.0665,1.09791,1.11737,1.18229,1.26962,1.46026}, // 200 < jet pT < 300 GeV
      {1.01542,1.03429,1.05216,1.0916,1.14318,1.17209,1.2266,1.50431}, // 300 < jet pT < 500 GeV
      {0.999292,1.06331,1.13128,1.39252,0.860853,0.804866,1.00156,1.55059}}, // 500 < jet pT < 5020 GeV
    // Centrality 30-50%
    {{1.03259,1.06135,1.09904,1.13689,1.18594,1.20406,1.30674,1.46526}, // 120 < jet pT < 140 GeV
      {1.04369,1.06251,1.09863,1.1436,1.20815,1.25224,1.33174,1.51462}, // 140 < jet pT < 160 GeV
      {1.02079,1.05696,1.10183,1.10483,1.15322,1.29897,1.34689,1.50899}, // 160 < jet pT < 180 GeV
      {1.04224,1.05543,1.09219,1.15402,1.16793,1.24301,1.2559,1.55678}, // 180 < jet pT < 200 GeV
      {1.03616,1.05156,1.07223,1.17912,1.16249,1.21022,1.3317,1.49531}, // 200 < jet pT < 300 GeV
      {1.01021,1.08254,1.06674,1.12047,1.1412,1.25265,0.98907,1.6675}, // 300 < jet pT < 500 GeV
      {0.751401,1.18767,1.42633,1.04446,1.08159,1.22675,0.829435,1.52272}}, // 500 < jet pT < 5020 GeV
    // Centrality 50-90%
    {{1.08866,1.0919,1.18341,1.2572,1.30058,1.35086,1.41097,1.51184}, // 120 < jet pT < 140 GeV
      {1.07445,1.10666,1.19488,1.21595,1.2809,1.33918,1.40624,1.54237}, // 140 < jet pT < 160 GeV
      {1.06253,1.10385,1.11587,1.10688,1.35489,1.27738,1.51023,1.49006}, // 160 < jet pT < 180 GeV
      {1.07904,1.13127,1.1582,1.34338,1.33199,1.39426,1.33704,1.49161}, // 180 < jet pT < 200 GeV
      {1.06651,1.10501,1.17877,1.22793,1.38095,1.20649,1.5489,1.60173}, // 200 < jet pT < 300 GeV
      {0.932784,1.084,1.02406,1.16905,1.57546,1.44441,2.06656,1.49769}, // 300 < jet pT < 500 GeV
      {0.655846,1.0258,0.983349,0.570278,1.04107,-1,2.06656,1.21683}} // 500 < jet pT < 5020 GeV
  };
  
  // ===============================================
  // == Copy the input arrays to the class arrays ==
  // ===============================================
  
  for(int iCentrality = 0; iCentrality < kNCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < kNJetPtBins; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < kNTrackPtBins; iTrackPt++){
        fSignalConeMainShape[iCentrality][iJetPt][iTrackPt] = particleDensityPtBinned[iCentrality][iJetPt][iTrackPt];
        fSignalConeKink[iCentrality][iJetPt][iTrackPt] = particleDensityPtBinnedKink[iCentrality][iJetPt][iTrackPt];
        fReflectedConeShape[iCentrality][iJetPt][iTrackPt] = particleDensityPtBinnedReflectedCone[iCentrality][iJetPt][iTrackPt];
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
  
}

/*
 * Find a bin index from an array
 */
int ReflectedConeWeight::FindBinIndex(const double* array, const int nBins, const double value) const{
  
  // If the value is below the lowest index, return -1 to show an error
  if(value < array[0]) return -1;
  
  // Find the value from the arrey and return the bin index
  for(int iBin = 0; iBin < nBins; iBin++){
    if(value < array[iBin+1]) return iBin;
  }
  
  // If the value is higher than the upper bin border, return -1 to show an error
  return -1;
}


/*
 * Find centrality bin index
 */
int ReflectedConeWeight::FindCentralityBin(const double centrality) const{
  return FindBinIndex(kCentralityBinBorders, kNCentralityBins, centrality);
}

/*
 * Find jet pT bin index
 */
int ReflectedConeWeight::FindJetPtBin(const double jetPt) const{
  return FindBinIndex(kJetPtBinBorders, kNJetPtBins, jetPt);
}

/*
 * Find track pT bin index
 */
int ReflectedConeWeight::FindTrackPtBin(const double trackPt) const{
  return FindBinIndex(kTrackPtBinBorders, kNTrackPtBins, trackPt);
}

/*
 * Get the seagull veto flag for the given bin for RacoGen simulation
 *
 *  const double deltaR = DeltaR between the track and jet axis
 *  const double centrality = Centrality of the event
 *  const double jetPt = Signal jet pT
 *  const double trackPt = pT of the track in the reflected cone
 */
double ReflectedConeWeight::GetReflectedConeWeight(const double deltaR, const double centrality, const double jetPt, const double trackPt) const{
  
  // Find the bin which is used to determine the weight
  const int iCentrality = FindCentralityBin(centrality);
  const int iJetPt = FindJetPtBin(jetPt);
  const int iTrackPt = FindTrackPtBin(trackPt);
  
  // If we are out of bounds in some of the arrays, do not apply weight
  if(iCentrality == -1 || iJetPt == -1 || iTrackPt == -1) return 1;
  
  double relativeChange34, relativeChange06;
  double fractionOf34, fractionOf06;
  
  // If track pT is above 2.5 GeV, only weight by general shape. The kink is not important in this case.
  // Also, if we are before the kink region, use the weight by the general shape
  if(iTrackPt > 2.5 || deltaR < 0.34){
    
    // Find what the relative change up until 0.34 is
    relativeChange34 = fSignalConeMainShape[iCentrality][iJetPt][iTrackPt] - fReflectedConeShape[iCentrality][iJetPt][iTrackPt];
    
    // Find what fraction of this change we need to apply
    fractionOf34 = deltaR / 0.34;
    
    // Return the weight besed the the relative change and the fraction of it we need to apply
    return 1 - (relativeChange34 * fractionOf34);
    
  } else {
    // This is the kink region for low pT tracks
    
    // Find what the relative change up until 0.34 is
    relativeChange34 = fSignalConeMainShape[iCentrality][iJetPt][iTrackPt] - fReflectedConeShape[iCentrality][iJetPt][iTrackPt];
    
    // Find what the relative change from 0.34 to deltaR is
    relativeChange06 = fSignalConeKink[iCentrality][iJetPt][iTrackPt] - (1 + (fReflectedConeShape[iCentrality][iJetPt][iTrackPt]-1) * (0.06/0.34));
    
    // Find the fraction of the kink change we need to apply
    fractionOf06 = (deltaR - 0.34) / 0.06;
    
    // Return the weight based on the relative weight up to the kink and the kink weight
    return 1 - (relativeChange34 + relativeChange06*fractionOf06);
    
  }
  
}
