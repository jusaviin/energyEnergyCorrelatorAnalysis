/*
 * Implementation of the JetMetScalingFactorManager class
 */

// Own includes
#include "JetMetScalingFactorManager.h"


/*
 * Contructor
 */
JetMetScalingFactorManager::JetMetScalingFactorManager():
  fIsPbPbData(true),
  fSystematicIndex(0)
{

  // Find the scaling factors based on data type
  InitializeArrays();
}

/*
 * Custom constructor
 */
JetMetScalingFactorManager::JetMetScalingFactorManager(const bool isPbPbData, const int scalingFactorType):
  fIsPbPbData(isPbPbData)
{

  if(scalingFactorType < 0){
    fSystematicIndex = 0;
  } else {
    fSystematicIndex = scalingFactorType;
  }

  // Find the scaling factors based on data type
  InitializeArrays();
}

/*
 * Initialization for the scaling factor arrays
 */
void JetMetScalingFactorManager::InitializeArrays(){
  
  // ======================================================== //
  //                                                          //
  //  M     M   EEEEEE   TTTTTTT   H    H    OOOOO    DDDD    //
  //  MMM MMM   E           T      H    H   O     O   D   D   //
  //  M  M  M   EEEEEE      T      HHHHHH   O     O   D    D  //
  //  M     M   E           T      H    H   O     O   D   D   //
  //  M     M   EEEEEE      T      H    H    OOOOO    DDDD    //
  //                                                          //
  // ======================================================== //
  
  // The jet energy resolution scaling factors are determined by the JetMet group
  // We are using the scaling factors from the closest pp run periods of the PbPb data

  // The numbers for PbPb are taken from file: Autumn18_RunD_V7b_MC_SF_AK4PF.txt
  if(fIsPbPbData){

    fScalingFactors[0][kNominal] = 1.1742; // Nominal scaling factor for 0 < |jet eta| < 0.522
    fScalingFactors[1][kNominal] = 1.1930; // Nominal scaling factor for 0.522 < |jet eta| < 0.783
    fScalingFactors[2][kNominal] = 1.1451; // Nominal scaling factor for 0.783 < |jet eta| < 1.131
    fScalingFactors[3][kNominal] = 1.1618; // Nominal scaling factor for 1.131 < |jet eta| < 1.305
    fScalingFactors[4][kNominal] = 1.1455; // Nominal scaling factor for 1.305 < |jet eta| < 1.740

    fScalingFactors[0][kUncertaintyDown] = 1.1415; // Down uncertainty scaling factor for 0 < |jet eta| < 0.522
    fScalingFactors[1][kUncertaintyDown] = 1.1559; // Down uncertainty scaling factor for 0.522 < |jet eta| < 0.783
    fScalingFactors[2][kUncertaintyDown] = 1.0812; // Down uncertainty scaling factor for 0.783 < |jet eta| < 1.131
    fScalingFactors[3][kUncertaintyDown] = 1.1086; // Down uncertainty scaling factor for 1.131 < |jet eta| < 1.305
    fScalingFactors[4][kUncertaintyDown] = 1.0838; // Down uncertainty scaling factor for 1.305 < |jet eta| < 1.740

    fScalingFactors[0][kUncertaintyUp] = 1.2069; // Up uncertainty scaling factor for 0 < |jet eta| < 0.522
    fScalingFactors[1][kUncertaintyUp] = 1.2302; // Up uncertainty scaling factor for 0.522 < |jet eta| < 0.783
    fScalingFactors[2][kUncertaintyUp] = 1.2089; // Up uncertainty scaling factor for 0.783 < |jet eta| < 1.131
    fScalingFactors[3][kUncertaintyUp] = 1.2150; // Up uncertainty scaling factor for 1.131 < |jet eta| < 1.305
    fScalingFactors[4][kUncertaintyUp] = 1.2072; // Up uncertainty scaling factor for 1.305 < |jet eta| < 1.740


  // The numbers for pp are taken from file: Fall17_V3b_MC_SF_AK4PF.txt
  } else {

    fScalingFactors[0][kNominal] = 1.1432; // Nominal scaling factor for 0 < |jet eta| < 0.522
    fScalingFactors[1][kNominal] = 1.1815; // Nominal scaling factor for 0.522 < |jet eta| < 0.783
    fScalingFactors[2][kNominal] = 1.0989; // Nominal scaling factor for 0.783 < |jet eta| < 1.131
    fScalingFactors[3][kNominal] = 1.1137; // Nominal scaling factor for 1.131 < |jet eta| < 1.305
    fScalingFactors[4][kNominal] = 1.1307; // Nominal scaling factor for 1.305 < |jet eta| < 1.740

    fScalingFactors[0][kUncertaintyDown] = 1.1210; // Down uncertainty scaling factor for 0 < |jet eta| < 0.522
    fScalingFactors[1][kUncertaintyDown] = 1.1332; // Down uncertainty scaling factor for 0.522 < |jet eta| < 0.783
    fScalingFactors[2][kUncertaintyDown] = 1.0533; // Down uncertainty scaling factor for 0.783 < |jet eta| < 1.131
    fScalingFactors[3][kUncertaintyDown] = 0.9740; // Down uncertainty scaling factor for 1.131 < |jet eta| < 1.305
    fScalingFactors[4][kUncertaintyDown] = 0.9837; // Down uncertainty scaling factor for 1.305 < |jet eta| < 1.740

    fScalingFactors[0][kUncertaintyUp] = 1.1654; // Up uncertainty scaling factor for 0 < |jet eta| < 0.522
    fScalingFactors[1][kUncertaintyUp] = 1.2299; // Up uncertainty scaling factor for 0.522 < |jet eta| < 0.783
    fScalingFactors[2][kUncertaintyUp] = 1.1444; // Up uncertainty scaling factor for 0.783 < |jet eta| < 1.131
    fScalingFactors[3][kUncertaintyUp] = 1.2533; // Up uncertainty scaling factor for 1.131 < |jet eta| < 1.305
    fScalingFactors[4][kUncertaintyUp] = 1.2778; // Up uncertainty scaling factor for 1.305 < |jet eta| < 1.740

  }

  // After we have the numbers from the file, we need to process them as factors to be used with the smearing procedure
  // To do this, a following procedure is used:
  // We assume the jet energy resolution is a Gaussian distribution with some certain sigma. If you would like to add a Gaussian noise to make it worse, the sigma gets larger. The process obeys the random variable rule: X=Y+Z, where Y~N(y, sigmay) and Z~N(z,sigmaz), then X~N(y+z, sqrt(sigmay^2+sigmaz^2))). In this case, we assume that noise and the resolution are independent.
  // So let assume the sigmay is the jet energy resolution, then you want the sigmax = 1.2sigmay
  // which means that the sigmaz = sigmay * sqrt(1.2*1.2 - 1)

  for(int iEta = 0; iEta < kNJetEtaBins; iEta++){
    for(int iSystematic = 0; iSystematic < kNScalingFactorTypes; iSystematic++){
      fScalingFactors[iEta][iSystematic] = TMath::Sqrt(TMath::Max(fScalingFactors[iEta][iSystematic]*fScalingFactors[iEta][iSystematic] - 1.0, 0.0));
    }
  }

  
}

// Getter for the JetMet scaling factor corresponding to jet eta
double JetMetScalingFactorManager::GetScalingFactor(double jetEta) const{

  // Transform jet eta into positive number in the defined range
  if(jetEta < 0) jetEta = -jetEta;
  if(jetEta >= fJetEtaBinBorders[kNJetEtaBins]) jetEta = fJetEtaBinBorders[kNJetEtaBins] - 0.0001;

  // Find the eta bin corresponding to the input jet eta
  int iJetEta = 0;
  for(int iEta = 1; iEta < kNJetEtaBins+1; iEta++){
    if(jetEta < fJetEtaBinBorders[iEta]) break;
    iJetEta++;
  }

  // Return the scaling factor
  return fScalingFactors[iJetEta][fSystematicIndex];

}