#ifndef TRKEFF2017PP
#define TRKEFF2017PP

#include "TFile.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TMath.h"
#include <iostream>
#include <string>
#include "TrackingEfficiencyInterface.h"

class TrkEff2017pp : public TrackingEfficiencyInterface{
public:

  TrkEff2017pp( bool isQuiet_ = false ,std::string filePath = "");
  virtual ~TrkEff2017pp();

  float getCorrection(float pt, float eta);
  float getCorrection(float pt, float eta, int hiBin);
  float getEfficiency( float pt, float eta, bool passesCheck = false);
  float getFake( float pt, float eta, bool passesCheck = false);
  float getSecondary( float pt, float eta, bool passesCheck = false);


private:

  inline bool checkBounds(float pt, float eta);

  bool isQuiet;

  TFile * trkEff;
  TH2F * eff;
  TH2F * fake;
  TH2F * sec;

};

#endif
