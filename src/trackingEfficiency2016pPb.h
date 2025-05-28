#ifndef TRKEFF2016PPB
#define TRKEFF2016PPB

#include "TFile.h"
#include "TH2D.h"
#include "TMath.h"
#include <iostream>
#include <string>
#include "TrackingEfficiencyInterface.h"

class TrkEff2016pPb : public TrackingEfficiencyInterface{
public:

  TrkEff2016pPb( bool isQuiet = false ,std::string filePath = "");
  virtual ~TrkEff2016pPb();

  float getCorrection(float pt, float eta);
  float getCorrection(float pt, float eta, int hiBin);
  float getEfficiency(float pt, float eta, bool passesCheck = false);

private:

  inline bool checkBounds(float pt, float eta);

  bool fIsQuiet;

  TFile* fEfficiencyFile;
  TH2F* fEfficiencyMatrix;

};

#endif
