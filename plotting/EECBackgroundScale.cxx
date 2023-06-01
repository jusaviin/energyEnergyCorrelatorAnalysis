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
  
  double energyEnergyCorrelatorBackgroundScaleFactors[2][4][20][8] = 
  {
  // The first part of the array gives the results for reconstructed jets.
  // The file that was produced all the other numbers but the last two rows in each centrality bin is PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_matchJets_noJetPtWeight_processed_2023-05-26.root
  // The file that was used to produce the last two rows in each centrality bin is PbPbMC2018_RecoGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_finalMcWeight_matchJets_fixCentrality_processed_2023-03-06.root
  // The macro for the study is findBackgroundNormalizationScale.C
  // Git hash for all the other numbers except those in the last two rows is 4af9eb73e874758d159ccda143a265ed4abf29f1
  // Git hash which produced the numbers in the last two rows is b6d4713c0216f54c5ef5f0632137575923442134
  {

  // Centrality 0-10%
  {{1.39785,1.32176,1.20496,1.12722,1.09149,1.09222,1.12055,1.17931}, // 60 < jet pT < 70 GeV
   {1.35642,1.27865,1.16101,1.07808,1.03014,1.01009,1.01418,1.03222}, // 70 < jet pT < 80 GeV
   {1.31772,1.23843,1.11785,1.02908,0.975607,0.945008,0.927762,0.92732}, // 80 < jet pT < 90 GeV
   {1.28357,1.20542,1.08786,1.00208,0.947962,0.909092,0.885459,0.873589}, // 90 < jet pT < 100 GeV
   {1.25439,1.1782,1.06449,0.980728,0.925323,0.878733,0.850852,0.832934}, // 100 < jet pT < 110 GeV
   {1.23287,1.16058,1.04762,0.964998,0.90549,0.857908,0.829481,0.807518}, // 110 < jet pT < 120 GeV
   {1.20602,1.13465,1.03053,0.951582,0.894932,0.85086,0.820355,0.798587}, // 120 < jet pT < 140 GeV
   {1.17994,1.11509,1.02249,0.953231,0.907155,0.872099,0.845826,0.836501}, // 140 < jet pT < 160 GeV
   {1.15993,1.101,1.01848,0.959038,0.916504,0.881203,0.864311,0.848766}, // 160 < jet pT < 180 GeV
   {1.14232,1.0885,1.01317,0.95741,0.918475,0.889429,0.868711,0.85658}, // 180 < jet pT < 200 GeV
   {1.13197,1.08393,1.01764,0.973473,0.947287,0.930409,0.927996,0.92468}, // 200 < jet pT < 220 GeV
   {1.12197,1.0776,1.01706,0.978097,0.959847,0.955303,0.961111,0.993752}, // 220 < jet pT < 240 GeV
   {1.11304,1.07143,1.01655,0.98573,0.975402,0.967289,0.987284,1.02075}, // 240 < jet pT < 260 GeV
   {1.10703,1.06788,1.01774,0.986007,0.974744,0.975559,0.989997,1.02782}, // 260 < jet pT < 280 GeV
   {1.09869,1.06089,1.01515,0.984108,0.973707,0.969588,0.996037,1.03932}, // 280 < jet pT < 300 GeV
   {1.08748,1.05697,1.01946,1.00176,1.0035,1.02786,1.07663,1.14076}, // 300 < jet pT < 500 GeV
   {1.06178,1.04461,1.03088,1.03486,1.05944,1.10846,1.19897,1.32231}, // 500 < jet pT < 5020 GeV
   {1.31863,1.24129,1.12451,1.04123,0.991244,0.962887,0.952942,0.956806}, // 60 GeV < jet pT
   {1.10881,1.06702,1.01089,0.974611,0.956623,0.947139,0.954082,0.972533}, // 200 < jet pT < 300 GeV
   {1.16515,1.10337,1.01569,0.95081,0.906026,0.871724,0.850157,0.838499}}, // 120 GeV < jet pT
   // Centrality 10-30%
  {{1.27835,1.2153,1.13135,1.082,1.0699,1.08619,1.12645,1.1818}, // 60 < jet pT < 70 GeV
   {1.23324,1.16835,1.07736,1.01772,0.98722,0.978573,0.987761,1.00909}, // 70 < jet pT < 80 GeV
   {1.19765,1.13344,1.03887,0.973323,0.933304,0.912293,0.905537,0.910225}, // 80 < jet pT < 90 GeV
   {1.16968,1.10548,1.01371,0.947881,0.899393,0.865494,0.849882,0.850625}, // 90 < jet pT < 100 GeV
   {1.14908,1.08784,1.00283,0.936972,0.897243,0.864217,0.837092,0.824712}, // 100 < jet pT < 110 GeV
   {1.133,1.07484,0.99181,0.927048,0.881303,0.850816,0.827682,0.807845}, // 110 < jet pT < 120 GeV
   {1.1181,1.06495,0.988863,0.932039,0.890785,0.857578,0.832805,0.828497}, // 120 < jet pT < 140 GeV
   {1.10393,1.0558,0.990223,0.943023,0.908237,0.879383,0.869357,0.863023}, // 140 < jet pT < 160 GeV
   {1.09561,1.05388,0.997654,0.960358,0.936797,0.922042,0.920323,0.930456}, // 160 < jet pT < 180 GeV
   {1.08653,1.04831,0.998404,0.968269,0.954748,0.948386,0.95859,0.984393}, // 180 < jet pT < 200 GeV
   {1.08412,1.05065,1.01069,0.990582,0.981933,0.996771,1.03604,1.07133}, // 200 < jet pT < 220 GeV
   {1.07799,1.04768,1.01109,0.991818,0.99197,1.01111,1.04875,1.10821}, // 220 < jet pT < 240 GeV
   {1.07391,1.04662,1.01592,1.00377,1.01053,1.04102,1.08775,1.14963}, // 240 < jet pT < 260 GeV
   {1.06991,1.04334,1.01352,1.00146,1.00698,1.03198,1.07047,1.1213}, // 260 < jet pT < 280 GeV
   {1.07028,1.04905,1.02703,1.0249,1.04742,1.09501,1.16822,1.2578}, // 280 < jet pT < 300 GeV
   {1.06366,1.04597,1.03058,1.03729,1.07201,1.13227,1.22344,1.34008}, // 300 < jet pT < 500 GeV
   {1.05174,1.04379,1.04752,1.07491,1.13574,1.2406,1.40061,1.59647}, // 500 < jet pT < 5020 GeV
   {1.20989,1.14708,1.05938,1.00015,0.96746,0.952816,0.953488,0.965812}, // 60 GeV < jet pT
   {1.06958,1.04145,1.00864,0.993383,0.993525,1.01506,1.0569,1.10702}, // 200 < jet pT < 300 GeV
   {1.09349,1.04909,0.987998,0.944806,0.915401,0.894741,0.886761,0.893141}}, // 120 GeV < jet pT
   // Centrality 30-50%
  {{1.15134,1.11558,1.07919,1.07728,1.11252,1.17766,1.27808,1.39886}, // 60 < jet pT < 70 GeV
   {1.10531,1.06626,1.01781,0.996646,0.996264,1.01803,1.06944,1.13493}, // 70 < jet pT < 80 GeV
   {1.08107,1.04323,0.991731,0.962536,0.947407,0.950357,0.973088,1.00171}, // 80 < jet pT < 90 GeV
   {1.06424,1.02563,0.972614,0.937403,0.921653,0.918889,0.918388,0.936001}, // 90 < jet pT < 100 GeV
   {1.0561,1.02018,0.973948,0.942598,0.925156,0.921355,0.940857,0.977592}, // 100 < jet pT < 110 GeV
   {1.05259,1.02064,0.976748,0.949038,0.938371,0.943403,0.963727,0.995001}, // 110 < jet pT < 120 GeV
   {1.05058,1.02098,0.98349,0.961493,0.953467,0.959778,0.979675,1.01552}, // 120 < jet pT < 140 GeV
   {1.05348,1.03027,1.00656,0.995338,1.00722,1.03021,1.07484,1.16617}, // 140 < jet pT < 160 GeV
   {1.05524,1.03621,1.01728,1.01681,1.03469,1.07341,1.1283,1.1772}, // 160 < jet pT < 180 GeV
   {1.05595,1.04075,1.03113,1.03967,1.08085,1.13687,1.22962,1.33283}, // 180 < jet pT < 200 GeV
   {1.06479,1.0545,1.05775,1.08504,1.14229,1.24457,1.37672,1.5324}, // 200 < jet pT < 220 GeV
   {1.06229,1.05274,1.0566,1.09028,1.15713,1.2589,1.38709,1.5468}, // 220 < jet pT < 240 GeV
   {1.06308,1.05731,1.0658,1.11164,1.19481,1.32407,1.47483,1.69594}, // 240 < jet pT < 260 GeV
   {1.06471,1.06049,1.06856,1.11276,1.19173,1.33034,1.51113,1.72596}, // 260 < jet pT < 280 GeV
   {1.0662,1.06531,1.0841,1.14137,1.25624,1.42633,1.65717,1.93388}, // 280 < jet pT < 300 GeV
   {1.0741,1.07803,1.10965,1.18457,1.32182,1.52567,1.79058,2.10336}, // 300 < jet pT < 500 GeV
   {1.0847,1.09802,1.15441,1.26664,1.47394,1.77887,2.19516,2.65188}, // 500 < jet pT < 5020 GeV
   {1.10088,1.06475,1.02097,1.00202,1.00558,1.02873,1.07358,1.13262}, // 60 GeV < jet pT
   {1.05942,1.05222,1.05893,1.09465,1.16399,1.27751,1.42002,1.59819}, // 200 < jet pT < 300 GeV
   {1.04791,1.02621,1.00343,0.996898,1.01039,1.04142,1.09119,1.16301}}, // 120 GeV < jet pT
   // Centrality 50-90%
  {{1.17713,1.18053,1.22995,1.33769,1.51567,1.76498,2.13387,2.59061}, // 60 < jet pT < 70 GeV
   {1.1408,1.13641,1.16106,1.24259,1.3658,1.54251,1.77676,2.06929}, // 70 < jet pT < 80 GeV
   {1.12262,1.11612,1.13466,1.19295,1.30129,1.45365,1.64877,1.85582}, // 80 < jet pT < 90 GeV
   {1.11736,1.11213,1.132,1.19304,1.2973,1.43844,1.63387,1.83964}, // 90 < jet pT < 100 GeV
   {1.12286,1.11872,1.14732,1.21523,1.33743,1.5134,1.67082,1.86788}, // 100 < jet pT < 110 GeV
   {1.13069,1.13357,1.17006,1.2608,1.40093,1.56932,1.76113,1.98052}, // 110 < jet pT < 120 GeV
   {1.13938,1.14703,1.19253,1.29492,1.44109,1.65854,1.88541,2.14144}, // 120 < jet pT < 140 GeV
   {1.15876,1.17276,1.23901,1.36839,1.57998,1.8632,2.15946,2.45395}, // 140 < jet pT < 160 GeV
   {1.17428,1.19737,1.28219,1.45559,1.73759,2.11363,2.56425,2.95137}, // 160 < jet pT < 180 GeV
   {1.18383,1.21066,1.31402,1.4805,1.75656,2.12419,2.49648,2.85306}, // 180 < jet pT < 200 GeV
   {1.20237,1.23806,1.36368,1.5998,1.97996,2.49935,3.12101,3.86156}, // 200 < jet pT < 220 GeV
   {1.21235,1.25408,1.393,1.65535,2.05224,2.58806,3.301,3.99866}, // 220 < jet pT < 240 GeV
   {1.21416,1.25431,1.393,1.64616,2.03427,2.51036,3.0995,3.73129}, // 240 < jet pT < 260 GeV
   {1.24483,1.29904,1.47083,1.81265,2.35156,3.16484,4.2591,5.29721}, // 260 < jet pT < 280 GeV
   {1.23772,1.28349,1.45424,1.77095,2.28124,3.07379,4.02943,4.98414}, // 280 < jet pT < 300 GeV
   {1.26566,1.32965,1.53022,1.90496,2.5248,3.44323,4.54679,5.87069}, // 300 < jet pT < 500 GeV
   {1.3202,1.40733,1.69838,2.22605,3.09695,4.4036,6.00362,8.01026}, // 500 < jet pT < 5020 GeV
   {1.14785,1.14789,1.18411,1.27345,1.41899,1.62162,1.8826,2.18167}, // 60 GeV < jet pT
   {1.21108,1.25163,1.38978,1.65001,2.05774,2.61467,3.31257,4.05708}, // 200 < jet pT < 300 GeV
   {1.15858,1.1757,1.24724,1.38973,1.60978,1.91682,2.2532,2.60344}} // 120 GeV < jet pT
  },
  // The second part of the array gives the same numbers for a generator level jet study
  // All but last two rows from file: PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_cutBadPhi_moreLowPtBins_truthReferenceForUnfolding_processed_2023-05-20.root
  // Last two rows from file: PbPbMC2018_GenGen_eecAnalysis_akFlowJets_miniAOD_4pCentShift_noTrigger_finalMcWeight_fixCentrality_processed_2023-03-08.root
  {
   // Centrality 0-10%
  {{1.38632,1.31176,1.19974,1.12733,1.08949,1.07698,1.08529,1.10874}, // 60 < jet pT < 70 GeV
   {1.34945,1.27827,1.17336,1.10679,1.07283,1.06014,1.07207,1.08778}, // 70 < jet pT < 80 GeV
   {1.32211,1.2548,1.15912,1.09759,1.07252,1.06606,1.07945,1.09609}, // 80 < jet pT < 90 GeV
   {1.29802,1.23408,1.14671,1.09444,1.07891,1.07641,1.09454,1.13376}, // 90 < jet pT < 100 GeV
   {1.27754,1.21699,1.13411,1.08709,1.07195,1.07664,1.09635,1.13436}, // 100 < jet pT < 110 GeV
   {1.25577,1.19795,1.12077,1.07623,1.06055,1.0649,1.08217,1.10877}, // 110 < jet pT < 120 GeV
   {1.23902,1.18572,1.11883,1.08248,1.07472,1.0927,1.12392,1.17018}, // 120 < jet pT < 140 GeV
   {1.21563,1.16776,1.10839,1.08069,1.08403,1.1047,1.15116,1.2221}, // 140 < jet pT < 160 GeV
   {1.1959,1.15201,1.10044,1.0766,1.08078,1.1065,1.16148,1.22497}, // 160 < jet pT < 180 GeV
   {1.17867,1.13894,1.09219,1.07495,1.08327,1.1142,1.16193,1.23011}, // 180 < jet pT < 200 GeV
   {1.16629,1.13086,1.09037,1.07657,1.09122,1.12229,1.17036,1.24263}, // 200 < jet pT < 220 GeV
   {1.1549,1.12156,1.08488,1.07812,1.09914,1.13906,1.19859,1.28865}, // 220 < jet pT < 240 GeV
   {1.1453,1.11349,1.08006,1.07483,1.09515,1.14174,1.21901,1.33506}, // 240 < jet pT < 260 GeV
   {1.13671,1.10782,1.07747,1.07168,1.09587,1.14033,1.22939,1.3424}, // 260 < jet pT < 280 GeV
   {1.13099,1.10385,1.07799,1.0773,1.10162,1.15371,1.23885,1.35287}, // 280 < jet pT < 300 GeV
   {1.11305,1.09084,1.07091,1.07554,1.10865,1.17429,1.27523,1.40506}, // 300 < jet pT < 500 GeV
   {1.08236,1.07061,1.06867,1.08913,1.138,1.21488,1.3498,1.53014}, // 500 < jet pT < 5020 GeV
   {1.33584,1.26682,1.16776,1.10642,1.07933,1.07409,1.08987,1.11761}, // 60 GeV < jet pT
   {1.14204,1.11183,1.07952,1.07312,1.09286,1.12971,1.1915,1.27956}, // 200 < jet pT < 300 GeV
   {1.19709,1.15238,1.0986,1.07378,1.07577,1.09889,1.14123,1.20407}}, // 120 GeV < jet pT
  // Centrality 10-30%
  {{1.28655,1.22645,1.14635,1.0998,1.0816,1.08798,1.11018,1.13755}, // 60 < jet pT < 70 GeV
   {1.25884,1.2043,1.13254,1.0952,1.08691,1.09753,1.12702,1.16256}, // 70 < jet pT < 80 GeV
   {1.23471,1.18452,1.11984,1.08723,1.08727,1.10668,1.13481,1.17444}, // 80 < jet pT < 90 GeV
   {1.21756,1.17054,1.11331,1.08593,1.08016,1.10028,1.12619,1.1749}, // 90 < jet pT < 100 GeV
   {1.19959,1.15649,1.10554,1.08434,1.08841,1.10794,1.14446,1.20081}, // 100 < jet pT < 110 GeV
   {1.18567,1.14559,1.09813,1.07843,1.08799,1.12114,1.16888,1.23529}, // 110 < jet pT < 120 GeV
   {1.17321,1.1368,1.0968,1.08599,1.10024,1.13636,1.197,1.27227}, // 120 < jet pT < 140 GeV
   {1.15641,1.12533,1.0928,1.0892,1.11245,1.15647,1.23254,1.32013}, // 140 < jet pT < 160 GeV
   {1.1425,1.11495,1.08808,1.09018,1.11576,1.16525,1.2439,1.34708}, // 160 < jet pT < 180 GeV
   {1.13441,1.1108,1.08891,1.09579,1.1277,1.18863,1.28801,1.41329}, // 180 < jet pT < 200 GeV
   {1.12285,1.10073,1.08345,1.09428,1.1321,1.20086,1.30668,1.43032}, // 200 < jet pT < 220 GeV
   {1.11819,1.09855,1.08639,1.1011,1.14556,1.22916,1.34348,1.50499}, // 220 < jet pT < 240 GeV
   {1.11289,1.095,1.08681,1.10559,1.15685,1.24329,1.36041,1.50304}, // 240 < jet pT < 260 GeV
   {1.10728,1.09239,1.08656,1.10199,1.1542,1.24499,1.37962,1.54298}, // 260 < jet pT < 280 GeV
   {1.10102,1.08775,1.08363,1.10759,1.16454,1.26182,1.40848,1.59629}, // 280 < jet pT < 300 GeV
   {1.0927,1.08312,1.08471,1.11618,1.18522,1.29604,1.45237,1.66051}, // 300 < jet pT < 500 GeV
   {1.07596,1.07423,1.09108,1.13606,1.22184,1.36942,1.59333,1.86902}, // 500 < jet pT < 5020 GeV
   {1.24657,1.19409,1.12707,1.09299,1.08744,1.10427,1.13735,1.18154}, // 60 GeV < jet pT
   {1.10868,1.09138,1.08198,1.0981,1.14311,1.22177,1.33707,1.47936}, // 200 < jet pT < 300 GeV
   {1.14364,1.11512,1.08685,1.08577,1.10946,1.15741,1.23448,1.33111}}, // 120 GeV < jet pT
   // Centrality 30-50%
  {{1.17765,1.14814,1.11966,1.12211,1.15436,1.21149,1.29438,1.37825}, // 60 < jet pT < 70 GeV
   {1.16244,1.13649,1.11563,1.12713,1.1641,1.22336,1.3139,1.41572}, // 70 < jet pT < 80 GeV
   {1.15072,1.12834,1.11428,1.12993,1.16685,1.23223,1.3159,1.42168}, // 80 < jet pT < 90 GeV
   {1.14244,1.12377,1.11329,1.13278,1.17945,1.24958,1.34569,1.44076}, // 90 < jet pT < 100 GeV
   {1.13568,1.1198,1.11464,1.14004,1.19908,1.28572,1.40552,1.5607}, // 100 < jet pT < 110 GeV
   {1.12901,1.11542,1.11681,1.14358,1.21141,1.31483,1.45401,1.64396}, // 110 < jet pT < 120 GeV
   {1.12506,1.11466,1.12082,1.15849,1.23363,1.33569,1.46617,1.65172}, // 120 < jet pT < 140 GeV
   {1.11807,1.11154,1.12404,1.16919,1.25378,1.39288,1.57655,1.81604}, // 140 < jet pT < 160 GeV
   {1.11717,1.11474,1.13592,1.18671,1.28907,1.44308,1.64878,1.89172}, // 160 < jet pT < 180 GeV
   {1.11125,1.11088,1.13504,1.19714,1.31016,1.48738,1.72639,2.01935}, // 180 < jet pT < 200 GeV
   {1.11118,1.11226,1.14211,1.21572,1.33922,1.53104,1.80117,2.12328}, // 200 < jet pT < 220 GeV
   {1.109,1.11295,1.14311,1.21726,1.35376,1.54286,1.82591,2.18296}, // 220 < jet pT < 240 GeV
   {1.1066,1.11183,1.14934,1.23227,1.37361,1.59756,1.88549,2.27511}, // 240 < jet pT < 260 GeV
   {1.10424,1.11061,1.14634,1.23192,1.37793,1.61325,1.89952,2.28649}, // 260 < jet pT < 280 GeV
   {1.10969,1.11975,1.16881,1.26102,1.42611,1.68935,2.03155,2.46165}, // 280 < jet pT < 300 GeV
   {1.10965,1.12316,1.17639,1.28909,1.48248,1.7776,2.17915,2.66338}, // 300 < jet pT < 500 GeV
   {1.11354,1.1356,1.21514,1.36047,1.62087,2.02412,2.58572,3.25207}, // 500 < jet pT < 5020 GeV
   {1.15738,1.13423,1.11789,1.13291,1.1777,1.25023,1.35199,1.46979}, // 60 GeV < jet pT
   {1.10473,1.10922,1.14383,1.22247,1.35651,1.56217,1.84996,2.21101}, // 200 < jet pT < 300 GeV
   {1.11356,1.10885,1.12503,1.17484,1.26895,1.40669,1.59,1.82841}}, // 120 GeV < jet pT
   // Centrality 50-90%
  {{1.21944,1.23164,1.29256,1.41725,1.62002,1.87873,2.18221,2.57361}, // 60 < jet pT < 70 GeV
   {1.21407,1.22648,1.28958,1.42374,1.61992,1.89855,2.22607,2.59089}, // 70 < jet pT < 80 GeV
   {1.21274,1.22976,1.30287,1.44571,1.6663,1.97596,2.37649,2.77648}, // 80 < jet pT < 90 GeV
   {1.21738,1.2382,1.32278,1.47824,1.75435,2.12838,2.51666,2.91635}, // 90 < jet pT < 100 GeV
   {1.21902,1.24171,1.33283,1.51023,1.76633,2.13836,2.59744,3.15097}, // 100 < jet pT < 110 GeV
   {1.22841,1.25873,1.36303,1.55293,1.87618,2.3358,2.92476,3.59651}, // 110 < jet pT < 120 GeV
   {1.22915,1.26038,1.37046,1.5733,1.91025,2.41459,3.04212,3.77907}, // 120 < jet pT < 140 GeV
   {1.23666,1.27262,1.40099,1.64212,2.0416,2.59387,3.26911,4.07981}, // 140 < jet pT < 160 GeV
   {1.24773,1.29346,1.43459,1.69488,2.13996,2.81735,3.75066,4.76215}, // 160 < jet pT < 180 GeV
   {1.25506,1.30211,1.46002,1.75381,2.23157,2.88952,3.66253,4.53526}, // 180 < jet pT < 200 GeV
   {1.25409,1.30632,1.47174,1.78836,2.27256,2.98562,3.98609,5.10848}, // 200 < jet pT < 220 GeV
   {1.27015,1.32771,1.50801,1.84755,2.39034,3.18186,4.23529,5.53934}, // 220 < jet pT < 240 GeV
   {1.28668,1.35177,1.55965,1.95601,2.61671,3.62153,4.97261,6.59301}, // 240 < jet pT < 260 GeV
   {1.28346,1.34384,1.54878,1.93656,2.56677,3.52397,4.72041,6.02032}, // 260 < jet pT < 280 GeV
   {1.29112,1.36128,1.57741,1.9781,2.65157,3.65699,4.9165,6.58861}, // 280 < jet pT < 300 GeV
   {1.31249,1.39034,1.63676,2.09984,2.87619,4.09824,5.64267,7.55251}, // 300 < jet pT < 500 GeV
   {1.36059,1.46433,1.79656,2.40811,3.43729,5.00804,7.03842,9.53993}, // 500 < jet pT < 5020 GeV
   {1.2196,1.23741,1.31372,1.46476,1.7057,2.0389,2.44091,2.90668}, // 60 GeV < jet pT
   {1.26656,1.32405,1.50711,1.85368,2.42182,3.26258,4.38555,5.71303}, // 200 < jet pT < 300 GeV
   {1.23819,1.27774,1.41011,1.65563,2.06399,2.66084,3.4269,4.30642}} // 120 GeV < jet pT
  }};
   
  
  // Copy the input values to the class array
  for(int iCentrality = 0; iCentrality < kNCentralityBins; iCentrality++){
    for(int iJetPt = 0; iJetPt < kNJetPtBins; iJetPt++){
      for(int iTrackPt = 0; iTrackPt < kNTrackPtBins; iTrackPt++){
        fBackgroundScale[iCentrality][iJetPt][iTrackPt] = energyEnergyCorrelatorBackgroundScaleFactors[useGenJets][iCentrality][iJetPt][iTrackPt];
      } // Track pT loop
    } // Jet pT loop
  } // Centrality loop
  
}

/*
 * Get the background scaling factor for the given bin
 *
 *  const std::pair<double,double> centralityBinBorders = Centrality bin borders for the scaling factor
 *  const std::pair<double,double> jetPtBinBorders = Jet pT bin borders for the scaling factor
 *  const std::pair<double,double> trackPtBinBorders = Track pT bin borders for the scaling factor
 *
 *  return: Scaling factor corresponding to the bin with the given bin borders
 */
double EECBackgroundScale::GetEECBackgroundScale(const std::pair<double,double> centralityBinBorders, const std::pair<double,double> jetPtBinBorders, const std::pair<double,double> trackPtBinBorders) const{

  // ******************************************************************** //
  // First, find the bin indices that correspond to the input bin borders //
  // ******************************************************************** //

  // Small number
  const double epsilon = 0.0001;

  // Search if the given centrality bin borders are included in the background scale table
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
    std::cout << "EECBackgroundScale::Error! Centrality bin " << centralityBinBorders.first << "-" << centralityBinBorders.second << " % was not found from the table." << std::endl;
    std::cout << "Cannot provide a background scaling factor. Please check your code for errors." << std::endl;
    return -1;
  }

  // Search if the given track pT bin borders are included in the background scale table
  int trackPtIndex = -1;
  for(int iTrackPt = 0; iTrackPt < kNTrackPtBins; iTrackPt++){
    if(TMath::Abs(fTrackPtBinBorderLow[iTrackPt] - trackPtBinBorders.first) < epsilon){
      if(TMath::Abs(fTrackPtBinBorderHigh[iTrackPt] - trackPtBinBorders.second) < epsilon){
        trackPtIndex = iTrackPt;
        break;
      }
    }
  } // Loop over track pT bins included in the scaling table

  // If track pT bin is not found, print an error message and return -1 to show we did not find a proper scaling factor.
  if(trackPtIndex == -1){
    std::cout << "EECBackgroundScale::Error! Track pT bin " << trackPtBinBorders.first << "-" << trackPtBinBorders.second << " GeV was not found from the table." << std::endl;
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
    std::cout << "EECBackgroundScale::Error! Jet pT bin " << jetPtBinBorders.first << "-" << jetPtBinBorders.second << " GeV was not found from the table." << std::endl;
    std::cout << "Cannot provide a background scaling factor. Please check your code for errors." << std::endl;
    return -1;
  }

  // *************************************************************************** //
  // After the correct bin has been found, return the corresponding scale factor //
  // *************************************************************************** //

  return fBackgroundScale[centralityIndex][jetPtIndex][trackPtIndex];

}
