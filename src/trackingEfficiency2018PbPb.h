#ifndef TRKEFF2018PBPB
#define TRKEFF2018PBPB

#include "TFile.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TMath.h"
#include <iostream>
#include <string>
#include "TrackingEfficiencyInterface.h"

class TrkEff2018PbPb : public TrackingEfficiencyInterface{
public:

  TrkEff2018PbPb( std::string collectionName = "general", bool isQuiet_ = false ,std::string filePath = "");
  virtual ~TrkEff2018PbPb();

  float getCorrection(float pt, float eta, int hiBin);
  float getEfficiency( float pt, float eta, int hiBin, bool passesCheck = false);
  float getFake( float pt, float eta, int hiBin, bool passesCheck = false);

private:

  inline bool checkBounds(float pt, float eta, int hiBin);

  std::string mode;
  bool isQuiet;

  TFile * trkEff;
  TFile * trkFake;
  TH3F * eff;
  TH3F * fake;

  TH2D * effPix[5];

};


#endif
