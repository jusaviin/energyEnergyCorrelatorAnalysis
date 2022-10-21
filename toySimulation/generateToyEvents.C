/*
 * Get deltaR between two objects
 *
 *  Arguments:
 *   const double eta1 = Eta of the first object
 *   const double phi1 = Phi of the first object
 *   const double eta2 = Eta of the second object
 *   const double phi2 = Phi of the second object
 *
 *  return: DeltaR between the two objects
 */
double getDeltaR(const double eta1, const double phi1, const double eta2, const double phi2){

  double deltaEta = eta1 - eta2;
  double deltaPhi = phi1 - phi2;
  
  // Transform deltaPhi to interval [-pi,pi]
  while(deltaPhi > TMath::Pi()){deltaPhi += -2*TMath::Pi();}
  while(deltaPhi < -TMath::Pi()){deltaPhi += 2*TMath::Pi();}
  
  // Return the distance between the objects
  return TMath::Sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);
  
}

/*
 * Macro to generate toy events to study energy-energy correlators. The macro generates certain amount of events.
 * Each event has a multiplicity that follows the measured multiplicity dustribution within the jet cone.
 *
 *  Arguments:
 *    int nEvent = Number of events that is generated
 *    double slope = Slope for particle density. This says how much more probable it is to generate particle to deltaR = 0 compared to deltaR = 0.4
 *    double jetR = Used jet radius
 */
void generateToyEvents(int nEvent = 100000, double slope = 0.1, double jetR = 0.4){
  
  // Read the track pT and multiplicity distributions from a file
  TFile *inputFile = TFile::Open("../data/eecAnalysis_akFlowJets_updatedMultiplicityAndDensity_wtaAxis_preprocessed_2022-10-17.root");
  TH1D *hTrackPt = (TH1D*) inputFile->Get("track/trackPt_C0");
  TH1D *hMultiplicity = (TH1D*) inputFile->Get("multiplicityInReflectedCone/multiplicityInReflectedCone_C0T0");
  
  // Define the histograms that are filled in the event loop
  TH1D *hParticleDensity[2];
  hParticleDensity[0] = new TH1D("particleDensity", "particleDensity", TMath::Nint(jetR*100), 0, jetR);
  hParticleDensity[1] = new TH1D("particlePtDensity", "particlePtDensity", TMath::Nint(jetR*100), 0, jetR);
  
  // Logarithmic deltaR binning for energy-energy correlator histograms
  const int nDeltaRBinsEEC = 60;
  const double minDeltaREEC = 0;
  const double maxDeltaREEC = 0.8;
  const double binnerShift = 0.01;
  const double deltaRlogBinWidth = (TMath::Log(maxDeltaREEC+binnerShift) - TMath::Log(minDeltaREEC+binnerShift)) / nDeltaRBinsEEC;
  double deltaRBinsEEC[nDeltaRBinsEEC+1];
  for(int iDeltaR = 0; iDeltaR <= nDeltaRBinsEEC; iDeltaR++){
    deltaRBinsEEC[iDeltaR] = (minDeltaREEC+binnerShift)*TMath::Exp(iDeltaR*deltaRlogBinWidth)-binnerShift;
  }
  
  TH1D *hEnergyEnergyCorrelator = new TH1D("energyEnergyCorrelator", "energyEnergyCorrelator", nDeltaRBinsEEC, deltaRBinsEEC);
  TH1D *hMultiplicityOut = new TH1D("multiplicity", "multiplicity", 150, -0.5, 149.5);
  
  // Create a random number generator
  TRandom3 *rng = new TRandom3();
  rng->SetSeed(0);
  
  // Helper variables used in the loop
  int currentMultiplicity;
  double currentPt, currentEta, currentPhi;
  std::vector<double> particlePt;
  std::vector<double> particleEta;
  std::vector<double> particlePhi;
  double correlatorWeight, currentDeltaR;
  
  // Function for inducing a slope to the particle density distribution
  TF1* slopeFunction = new TF1("slopeFunction", "pol1", 0, 0.4);
  slopeFunction->SetParameter(0,0);
  slopeFunction->SetParameter(1,slope/jetR);
    
  // Generate nEvent toy events
  cout << "Generating " << nEvent << " events" << endl;
  for(int iEvent = 0; iEvent < nEvent; iEvent++){
    
    if(iEvent % 10000 == 0) cout << "Generating event " << iEvent << endl;
    
    // Clear the particle vectors from the previous event
    particlePt.clear();
    particlePhi.clear();
    particleEta.clear();
    
    // First, determine the particle multiplicity within the jet cone
    currentMultiplicity = TMath::Nint(hMultiplicity->GetRandom(rng));
    hMultiplicityOut->Fill(currentMultiplicity);
    
    // Create particles uniformly until we have the required number of particles
    for(int iParticle = 0; iParticle < currentMultiplicity; iParticle++){
      
      // Get the particle pT from the measured track pT distribution
      currentPt = hTrackPt->GetRandom(rng);
      particlePt.push_back(currentPt);
      
      // Get eta and phi from uniform distribution. Require that they are within the jet cone of 0.4
      do {
        currentEta = rng->Uniform(-jetR, jetR);
        currentPhi = rng->Uniform(-jetR, jetR);
        currentDeltaR = TMath::Sqrt(currentEta*currentEta + currentPhi*currentPhi);
      } while ((currentDeltaR > jetR) || (rng->Rndm() < slopeFunction->Eval(currentDeltaR)));
      particlePhi.push_back(currentPhi);
      particleEta.push_back(currentEta);
      
      // Fill the particle density histogram
      hParticleDensity[0]->Fill(TMath::Sqrt(currentEta*currentEta + currentPhi*currentPhi));
      hParticleDensity[1]->Fill(TMath::Sqrt(currentEta*currentEta + currentPhi*currentPhi), currentPt);
      
    } // End of particle loop
    
    // After the particles are generated for the event, calculate the energy-energy correlator
    for(int firstParticle = 0; firstParticle < currentMultiplicity; firstParticle++){
      for(int secondParticle = firstParticle+1; secondParticle < currentMultiplicity; secondParticle++){
        correlatorWeight = particlePt.at(firstParticle)*particlePt.at(secondParticle);
        currentDeltaR = getDeltaR(particleEta.at(firstParticle), particlePhi.at(firstParticle), particleEta.at(secondParticle), particlePhi.at(secondParticle));
        hEnergyEnergyCorrelator->Fill(currentDeltaR, correlatorWeight);
      }
    }
    
  } // End of event loop
  
  // Normalize the particle density histogram to the bin area
  double binLowBorder, binHighBorder;
  double binLowArea, binHighArea, binArea;
  double newBinContent, newBinError;
  
  // Loop over the bins and normalize each bin to the bin area
  for(int iParticleDensity = 0; iParticleDensity < 2; iParticleDensity++){
    for(int iBin = 1; iBin <= hParticleDensity[iParticleDensity]->GetNbinsX(); iBin++){
      
      // Calculate the bin area
      binLowBorder = hParticleDensity[iParticleDensity]->GetXaxis()->GetBinLowEdge(iBin);
      binHighBorder = hParticleDensity[iParticleDensity]->GetXaxis()->GetBinUpEdge(iBin);
      binLowArea = TMath::Pi() * binLowBorder * binLowBorder;
      binHighArea = TMath::Pi() * binHighBorder * binHighBorder;
      binArea = binHighArea - binLowArea;
      
      // Calculate the new bin content and error for particle density
      newBinContent = hParticleDensity[iParticleDensity]->GetBinContent(iBin) / binArea;
      newBinError = hParticleDensity[iParticleDensity]->GetBinError(iBin) / binArea;
      hParticleDensity[iParticleDensity]->SetBinContent(iBin, newBinContent);
      hParticleDensity[iParticleDensity]->SetBinError(iBin, newBinError);
      
    } // Bin loop for normalizing particle densities
  }
  
  // Normalize the energy-energy correlator by the bin width
  hEnergyEnergyCorrelator->Scale(1.0,"width");
  
  // After the histograms are normalized, we can save them to a file
  TFile *outputFile = new TFile("toySimulation100kevents10pSlope.root", "RECREATE");
  hParticleDensity[0]->Write();
  hParticleDensity[1]->Write();
  hEnergyEnergyCorrelator->Write();
  hMultiplicityOut->Write();
  outputFile->Close();
  
}
