#include "trackingEfficiency2016pPb.h"

inline bool TrkEff2016pPb::checkBounds(float pt, float eta){
  if( TMath::Abs(eta) > 2.4 ){
    if( ! fIsQuiet) std::cout << "TrkEff2016pPb: track outside |eta|<2.4, please apply this cut!  I am returning a correction factor of 0 for this track for now." << std::endl;
    return false;
  }
  
  if( pt< 0 || pt > 500 ){
    if( ! fIsQuiet) std::cout << "TrkEff2016pPb: pT is outside the range [0,500].  I am returning a correction factor of 0 for this track for now." << std::endl;
    return false;
  }

  return true;
}

float TrkEff2016pPb::getCorrection(float pt, float eta, int hiBin){
  return getCorrection(pt, eta);
}

float TrkEff2016pPb::getCorrection(float pt, float eta){
  if( !checkBounds(pt, eta) ) return 0;

  // Find track efficiency
  float efficiency = getEfficiency(pt, eta, true);

  //protect against dividing by 0
  if(efficiency > 0.001){
    return 1.0 / efficiency;
  } else {
    if( ! fIsQuiet ) std::cout << "TrkEff2016pPb: Warning! Tracking efficiency is very low for this track (close to dividing by 0).  Returning correction factor of 0 for this track for now." << std::endl;
    return 0;
  }

}

float TrkEff2016pPb::getEfficiency(float pt, float eta, bool passesCheck){
  if( !passesCheck){
    if(  !checkBounds(pt, eta) ) return 0;
  }

  return fEfficiencyMatrix->GetBinContent( fEfficiencyMatrix->FindBin(eta, pt) );
}


TrkEff2016pPb::TrkEff2016pPb(bool isQuiet, std::string filePath){
  fIsQuiet = isQuiet;
  
  std::cout << "Searching tracking corrections from folder: " << filePath << std::endl;
  
  if(!fIsQuiet) std::cout << "TrkEff2016pPb class opening in general tracks mode!" << std::endl;
  
  fEfficiencyFile = TFile::Open( (filePath + "Hijing_8TeV_dataBS.root").c_str(),"open");
  
  if( !(fEfficiencyFile->IsOpen() ) ){
    std::cout << "WARNING, COULD NOT FIND TRACK EFFICIENCY FILE FOR GENERAL TRACKS!" << std::endl;
  } else {
    fEfficiencyMatrix = (TH2F*) fEfficiencyFile->Get("rTotalEff3D_0");
  }
}

TrkEff2016pPb::~TrkEff2016pPb(){
  fEfficiencyFile->Close();
}
