#ifndef TRACKINGEFFICIENCYINTERFACE_H
#define TRACKINGEFFICIENCYINTERFACE_H

class TrackingEfficiencyInterface{
  
public:
  
  virtual float getCorrection(float pt, float eta, int hiBin) = 0;
  virtual ~TrackingEfficiencyInterface();
  
};

#endif
