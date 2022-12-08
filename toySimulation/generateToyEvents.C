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
 *    double ptCut = Apply a lower pT cut for the simulated particles
 *    double slope = Slope for particle density. This says how much more probable it is to generate particle to deltaR = 0 compared to deltaR = 0.4
 *    double jetR = Used jet radius
 */
void generateToyEvents(int nEvent = 100000, double ptCut = 0.0, double slope = 0.0, double jetR = 0.4){
  
  // Read the track pT and multiplicity distributions from a file
  TFile *inputFile = TFile::Open("../data/PbPbMC2018_RecoGen_eecAnalysis_akFlowJet_updatedMultiplicityAndDensity_wtaAxis_noTrigger_preprocessed_2022-10-17.root");
  TH1D *hTrackPt = (TH1D*) inputFile->Get("track/trackPt_C0");
  //TH1D *hMultiplicityJet = (TH1D*) inputFile->Get("multiplicityInJetCone/multiplicityInJetCone_C0T0S0");
  TH1D *hMultiplicityBackground = (TH1D*) inputFile->Get("multiplicityInJetCone/multiplicityInJetCone_C0T0S1");
  
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
  
  TH1D *hEnergyEnergyCorrelator = new TH1D("energyEnergyCorrelator", "energyEnergyCorrelator", nDeltaRBinsEEC, deltaRBinsEEC); hEnergyEnergyCorrelator->Sumw2();
  TH1D *hThreePointEnergyCorrelator = new TH1D("threePointEnergyCorrelator", "threePointEnergyCorrelator", nDeltaRBinsEEC, deltaRBinsEEC); hThreePointEnergyCorrelator->Sumw2();
  TH1D *hMultiplicityOut = new TH1D("multiplicity", "multiplicity", 150, -0.5, 149.5);
  TH1D *hTrackPtOut = (TH1D*) hTrackPt->Clone("trackPt");
  TH2D *hEtaPhi = new TH2D("etaPhi", "etaPhi", 100, -0.5, 0.5, 100, -0.5, 0.5);
  TH1D *hFlowPhase = new TH1D("flowPhase","flowPhase",64,-3.2,3.2);
  hTrackPtOut->Reset("M");
  
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
  double deltaR12, deltaR13, deltaR23;
  
  // Function for inducing a flow modulation for the particles
  TF1 *flowFunction = new TF1("flowFunction", "[1]*cos(2*x+[0])+1", -TMath::Pi(), TMath::Pi());
  flowFunction->SetParameter(0,0);
  flowFunction->SetParameter(1,0.1);
  
  // Function to induce a bias that the upwards flow fluctuation coindiced with the jet axis
  double biasPercentage = 0.1;
  TF1 *flowBiasFunction = new TF1("flowBias", "1-TMath::Abs(x)*[0]", -TMath::Pi(), TMath::Pi());
  flowBiasFunction->SetParameter(0, biasPercentage*2/TMath::Pi());
  
  // Function for inducing a slope to the particle density distribution
  TF1* slopeFunction = new TF1("slopeFunction", "pol1", 0, 0.4);
  slopeFunction->SetParameter(0,0);
  slopeFunction->SetParameter(1,slope/jetR);
    
  // Apply a lower pT cut to the generated particle distribution
  for(int iBin = 1; iBin <= hTrackPt->GetNbinsX(); iBin++){
    if(hTrackPt->GetXaxis()->GetBinUpEdge(iBin) < ptCut){
      hTrackPt->SetBinContent(iBin,0);
      hTrackPt->SetBinError(iBin,0);
    } else {
      break;
    }
  }
  
  // Do not allow 0 and 1 particles in the multiplicity distribution
  //hMultiplicityJet->SetBinContent(1,0);
  //hMultiplicityJet->SetBinError(1,0);
  hMultiplicityBackground->SetBinContent(1,0);
  hMultiplicityBackground->SetBinContent(2,0);
  hMultiplicityBackground->SetBinError(1,0);
  hMultiplicityBackground->SetBinError(2,0);
  
  // Generate nEvent toy events
  double flowPhase;
  for(int iEvent = 0; iEvent < nEvent; iEvent++){
    if(iEvent % 10000 == 0) cout << "Generating event " << iEvent << endl;
    
    // Clear the particle vectors from the previous event
    particlePt.clear();
    particlePhi.clear();
    particleEta.clear();
    
    // First, determine the particle multiplicity within the jet cone
    currentMultiplicity = TMath::Nint(hMultiplicityBackground->GetRandom(rng));  // Multiplicity sampled from data distribution
    //currentMultiplicity = TMath::Nint(rng->Uniform(2,140));          // Unoform multiplicity
    hMultiplicityOut->Fill(currentMultiplicity);
    
    // Find the phase for the flow for this event
    //flowPhase = rng->Uniform(-TMath::Pi()/2,TMath::Pi()/2);  // Completely random flow
    //flowPhase = flowBiasFunction->GetRandom(-TMath::Pi()/2,TMath::Pi()/2,rng);  // Make it more likely for an upwards fluctuation to be aligned with jet axis
    flowPhase = -TMath::Pi()/4;
    flowFunction->SetParameter(0, flowPhase);
    hFlowPhase->Fill(flowPhase);
    
    // Create particles uniformly until we have the required number of particles
    for(int iParticle = 0; iParticle < currentMultiplicity; iParticle++){
      
      // Get the particle pT from the measured track pT distribution
      currentPt = hTrackPt->GetRandom(rng);  // Track pT distribution sampled from data distribution
      //currentPt = rng->Uniform(0.7, 20);   // Uniform track pT
      particlePt.push_back(currentPt);
      hTrackPtOut->Fill(currentPt);
      
      // Get eta and phi from uniform distribution. Require that they are within the jet cone of 0.4
      do {
        currentEta = rng->Uniform(-jetR, jetR);  // Uniform eta
        currentPhi = rng->Uniform(-jetR, jetR);  // Uniform phi
        
        // Phi with random flow modulation
        // currentPhi = flowFunction->GetRandom(-jetR,jetR,rng);
        
        currentDeltaR = TMath::Sqrt(currentEta*currentEta + currentPhi*currentPhi);
      } while ((currentDeltaR > jetR)/* || (rng->Rndm() < slopeFunction->Eval(currentDeltaR))*/);
      particlePhi.push_back(currentPhi);
      particleEta.push_back(currentEta);
      hEtaPhi->Fill(currentPhi, currentEta);
      
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
    
    // Calculate also a three-point energy correlator
    for(int firstParticle = 0; firstParticle < currentMultiplicity; firstParticle++){
      for(int secondParticle = firstParticle+1; secondParticle < currentMultiplicity; secondParticle++){
        deltaR12 = getDeltaR(particleEta.at(firstParticle), particlePhi.at(firstParticle), particleEta.at(secondParticle), particlePhi.at(secondParticle));
        for(int thirdParticle = secondParticle+1; thirdParticle < currentMultiplicity; thirdParticle++){
          deltaR13 = getDeltaR(particleEta.at(firstParticle), particlePhi.at(firstParticle), particleEta.at(thirdParticle), particlePhi.at(thirdParticle));
          deltaR23 = getDeltaR(particleEta.at(secondParticle), particlePhi.at(secondParticle), particleEta.at(thirdParticle), particlePhi.at(thirdParticle));
          
          currentDeltaR = deltaR12;
          if(deltaR13 > currentDeltaR) currentDeltaR = deltaR13;
          if(deltaR23 > currentDeltaR) currentDeltaR = deltaR23;
          
          correlatorWeight = particlePt.at(firstParticle)*particlePt.at(secondParticle)*particlePt.at(thirdParticle);
          
          hThreePointEnergyCorrelator->Fill(currentDeltaR, correlatorWeight);
        }
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
  
  // Normalize the energy-energy correlators by the bin width
  hEnergyEnergyCorrelator->Scale(1.0,"width");
  hThreePointEnergyCorrelator->Scale(1.0,"width");
  
  // After the histograms are normalized, we can save them to a file
  TFile *outputFile = new TFile("toySimulation100keventsWith3Point.root", "RECREATE");
  hParticleDensity[0]->Write();
  hParticleDensity[1]->Write();
  hEnergyEnergyCorrelator->Write();
  hThreePointEnergyCorrelator->Write();
  hMultiplicityOut->Write();
  hTrackPtOut->Write();
  hEtaPhi->Write();
  hFlowPhase->Write();
  outputFile->Close();
  
}
