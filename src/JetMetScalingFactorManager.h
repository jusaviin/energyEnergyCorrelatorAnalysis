#ifndef JETMETSCALINGFACTORMANAGER_H
#define JETMETSCALINGFACTORMANAGER_H

// Root includes
#include <TMath.h>

/*
 * JetMetScalingFactorManager class
 *
 * Class whose purpose is to provide unfolding configuration for the data
 */
class JetMetScalingFactorManager {
  
public:
 
  // Enumeration for different 
  enum enumScalingFactorType{kNominal, kUncertaintyDown, kUncertaintyUp, kNScalingFactorTypes};

  static const int kNJetEtaBins = 5;      // Number of jet eta bins for included in the scaling factor tables
  
  JetMetScalingFactorManager();                                                    // Constructor
  JetMetScalingFactorManager(const bool isPbPbData, const int scalingFactorType);  // Custom constructor
  ~JetMetScalingFactorManager() = default;                                         // Destructor

  // Getter for JetMet scaling factor
  double GetScalingFactor(double jetEta) const; // Getter for the scaling factor corresponding to jet eta
  
private:
  
  // Private variables
  bool fIsPbPbData;                // Flag for PbPb data, read from card
  int fSystematicIndex;            // Index for systematic uncertainty study

  // Binning information for the scaling tables
  const double fJetEtaBinBorders[kNJetEtaBins+1] = {0, 0.522, 0.783, 1.131, 1.305, 1.740};

  // Array that holds the best number of iterations for each response matrix
  double fScalingFactors[kNJetEtaBins][kNScalingFactorTypes];
  
  // Initialize the number of iterations array based on the predefined configuration index
  void InitializeArrays();
  
};

#endif
