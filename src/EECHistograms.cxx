// Histograms needed in the dijet analysis

// C++ includes
#include <assert.h>

// Root includes
#include <TFile.h>
#include <TMath.h>

// Own includes
#include "EECHistograms.h"

/*
 * Default constructor
 */
EECHistograms::EECHistograms() :
  fhVertexZ(0),
  fhVertexZWeighted(0),
  fhEvents(0),
  fhTriggers(0),
  fhTriggersAfterSelection(0),
  fhTrackCuts(0),
  fhCentrality(0),
  fhCentralityWeighted(0),
  fhPtHat(0),
  fhPtHatWeighted(0),
  fhMultiplicity(0),
  fhInclusiveJet(0),
  fhTrack(0),
  fhTrackUncorrected(0),
  fhParticleDensityAroundJet(0),
  fhParticlePtDensityAroundJet(0),
  fhParticleMultiplicityInJet(0),
  fhParticleMultiplicityInReflectedCone(0),
  fhParticleMultiplicityInJetUncorrected(0),
  fhParticleMultiplicityInReflectedConeUncorrected(0),
  fhMaxPtParticleInJet(0),
  fhEnergyEnergyCorrelator(0),
  fhEnergyEnergyCorrelatorUncorrected(0),
  fhEnergyEnergyCorrelatorJetPt(0),
  fhEnergyEnergyCorrelatorJetPtUncorrected(0),
  fhJetPtClosure(0),
  fhUnfoldingMeasured(0),
  fhUnfoldingTruth(0),
  fhUnfoldingResponse(0),
  fCard(0)
{
  // Default constructor
  
}

/*
 * Custom constructor
 */
EECHistograms::EECHistograms(ConfigurationCard* newCard) :
  fhVertexZ(0),
  fhVertexZWeighted(0),
  fhEvents(0),
  fhTriggers(0),
  fhTriggersAfterSelection(0),
  fhTrackCuts(0),
  fhCentrality(0),
  fhCentralityWeighted(0),
  fhPtHat(0),
  fhPtHatWeighted(0),
  fhMultiplicity(0),
  fhInclusiveJet(0),
  fhTrack(0),
  fhTrackUncorrected(0),
  fhParticleDensityAroundJet(0),
  fhParticlePtDensityAroundJet(0),
  fhParticleMultiplicityInJet(0),
  fhParticleMultiplicityInReflectedCone(0),
  fhParticleMultiplicityInJetUncorrected(0),
  fhParticleMultiplicityInReflectedConeUncorrected(0),
  fhMaxPtParticleInJet(0),
  fhEnergyEnergyCorrelator(0),
  fhEnergyEnergyCorrelatorUncorrected(0),
  fhEnergyEnergyCorrelatorJetPt(0),
  fhEnergyEnergyCorrelatorJetPtUncorrected(0),
  fhJetPtClosure(0),
  fhUnfoldingMeasured(0),
  fhUnfoldingTruth(0),
  fhUnfoldingResponse(0),
  fCard(newCard)
{
  // Custom constructor

}

/*
 * Copy constructor
 */
EECHistograms::EECHistograms(const EECHistograms& in) :
  fhVertexZ(in.fhVertexZ),
  fhVertexZWeighted(in.fhVertexZWeighted),
  fhEvents(in.fhEvents),
  fhTriggers(in.fhTriggers),
  fhTriggersAfterSelection(in.fhTriggersAfterSelection),
  fhTrackCuts(in.fhTrackCuts),
  fhCentrality(in.fhCentrality),
  fhCentralityWeighted(in.fhCentralityWeighted),
  fhPtHat(in.fhPtHat),
  fhPtHatWeighted(in.fhPtHatWeighted),
  fhMultiplicity(in.fhMultiplicity),
  fhInclusiveJet(in.fhInclusiveJet),
  fhTrack(in.fhTrack),
  fhTrackUncorrected(in.fhTrackUncorrected),
  fhParticleDensityAroundJet(in.fhParticleDensityAroundJet),
  fhParticlePtDensityAroundJet(in.fhParticlePtDensityAroundJet),
  fhParticleMultiplicityInJet(in.fhParticleMultiplicityInJet),
  fhParticleMultiplicityInReflectedCone(in.fhParticleMultiplicityInReflectedCone),
  fhParticleMultiplicityInJetUncorrected(in.fhParticleMultiplicityInJetUncorrected),
  fhParticleMultiplicityInReflectedConeUncorrected(in.fhParticleMultiplicityInReflectedConeUncorrected),
  fhMaxPtParticleInJet(in.fhMaxPtParticleInJet),
  fhEnergyEnergyCorrelator(in.fhEnergyEnergyCorrelator),
  fhEnergyEnergyCorrelatorUncorrected(in.fhEnergyEnergyCorrelatorUncorrected),
  fhEnergyEnergyCorrelatorJetPt(in.fhEnergyEnergyCorrelatorJetPt),
  fhEnergyEnergyCorrelatorJetPtUncorrected(in.fhEnergyEnergyCorrelatorJetPtUncorrected),
  fhJetPtClosure(in.fhJetPtClosure),
  fhUnfoldingMeasured(in.fhUnfoldingMeasured),
  fhUnfoldingTruth(in.fhUnfoldingTruth),
  fhUnfoldingResponse(in.fhUnfoldingResponse),
  fCard(in.fCard)
{
  // Copy constructor

}

/*
 * Assingment operator
 */
EECHistograms& EECHistograms::operator=(const EECHistograms& in){
  // Assingment operator
  
  if (&in==this) return *this;
  
  fhVertexZ = in.fhVertexZ;
  fhVertexZWeighted = in.fhVertexZWeighted;
  fhEvents = in.fhEvents;
  fhTriggers = in.fhTriggers;
  fhTriggersAfterSelection = in.fhTriggersAfterSelection;
  fhTrackCuts = in.fhTrackCuts;
  fhCentrality = in.fhCentrality;
  fhCentralityWeighted = in.fhCentralityWeighted;
  fhPtHat = in.fhPtHat;
  fhPtHatWeighted = in.fhPtHatWeighted;
  fhMultiplicity = in.fhMultiplicity;
  fhInclusiveJet = in.fhInclusiveJet;
  fhTrack = in.fhTrack;
  fhTrackUncorrected = in.fhTrackUncorrected;
  fhParticleDensityAroundJet = in.fhParticleDensityAroundJet;
  fhParticlePtDensityAroundJet = in.fhParticlePtDensityAroundJet;
  fhParticleMultiplicityInJet = in.fhParticleMultiplicityInJet;
  fhParticleMultiplicityInReflectedCone = in.fhParticleMultiplicityInReflectedCone;
  fhParticleMultiplicityInJetUncorrected = in.fhParticleMultiplicityInJetUncorrected;
  fhParticleMultiplicityInReflectedConeUncorrected = in.fhParticleMultiplicityInReflectedConeUncorrected;
  fhMaxPtParticleInJet = in.fhMaxPtParticleInJet;
  fhEnergyEnergyCorrelator = in.fhEnergyEnergyCorrelator;
  fhEnergyEnergyCorrelatorUncorrected = in.fhEnergyEnergyCorrelatorUncorrected;
  fhEnergyEnergyCorrelatorJetPt = in.fhEnergyEnergyCorrelatorJetPt;
  fhEnergyEnergyCorrelatorJetPtUncorrected = in.fhEnergyEnergyCorrelatorJetPtUncorrected;
  fhJetPtClosure = in.fhJetPtClosure;
  fhUnfoldingMeasured = in.fhUnfoldingMeasured;
  fhUnfoldingTruth = in.fhUnfoldingTruth;
  fhUnfoldingResponse = in.fhUnfoldingResponse;
  fCard = in.fCard;
  
  return *this;
}

/*
 * Destructor
 */
EECHistograms::~EECHistograms(){
  // destructor
  delete fhVertexZ;
  delete fhVertexZWeighted;
  delete fhEvents;
  delete fhTriggers;
  delete fhTriggersAfterSelection;
  delete fhTrackCuts;
  delete fhCentrality;
  delete fhCentralityWeighted;
  delete fhPtHat;
  delete fhPtHatWeighted;
  delete fhMultiplicity;
  delete fhInclusiveJet;
  delete fhTrack;
  delete fhTrackUncorrected;
  delete fhParticleDensityAroundJet;
  delete fhParticlePtDensityAroundJet;
  delete fhParticleMultiplicityInJet;
  delete fhParticleMultiplicityInReflectedCone;
  delete fhParticleMultiplicityInJetUncorrected;
  delete fhParticleMultiplicityInReflectedConeUncorrected;
  delete fhMaxPtParticleInJet;
  delete fhEnergyEnergyCorrelator;
  delete fhEnergyEnergyCorrelatorUncorrected;
  delete fhEnergyEnergyCorrelatorJetPt;
  delete fhEnergyEnergyCorrelatorJetPtUncorrected;
  delete fhJetPtClosure;
  delete fhUnfoldingMeasured;
  delete fhUnfoldingTruth;
  delete fhUnfoldingResponse;
}

/*
 * Set the configuration card used for the histogram class
 */
void EECHistograms::SetCard(ConfigurationCard* newCard){
  fCard = newCard;
}

/*
 * Create the necessary histograms
 */
void EECHistograms::CreateHistograms(){
  
  // ======== Common binning information for histograms =========
  
  // Centrality
  const Double_t minCentrality = -1;   // Minimum centrality bin is negative since hiBin is -0.5 for pp
  const Double_t maxCentrality = 100;  // Maximum centrality bin
  const Int_t nCentralityBins = 202;      // Number of centrality bins
  
  // Jet pT
  const Double_t minPtJet = 0;     // Minimum jet pT
  const Double_t maxPtJet = 500;   // Maximum jet pT
  const Int_t nPtBinsJet = 100;    // Number of jet pT bins
  
  //Track pT
  const Double_t minPtTrack = 0;   // Minimum track pT for track histograms
  const Double_t maxPtTrack = 50;  // Maximum track pT for track histograms (Hallie's analysis = 20)
  const Int_t nPtBinsTrack = 500;  // Number of track pT bins for track histograms (Hallie's analysis = 500)
  
  // Track pT for maximum particle pT inside a jet
  const Double_t minParticlePt = 0;
  const Double_t maxParticlePt = 100;
  const Int_t nPtBinsParticle = 100;
  
  // Phi
  const Double_t minPhi = -TMath::Pi();  // Minimum phi
  const Double_t maxPhi = TMath::Pi();   // Maximum phi
  const Int_t nPhiBins = 64;             // Number of phi bins
  
  // Eta
  const Double_t minEta = -2.5;    // Minimum eta (current eta cut for tracks = 2.4)
  const Double_t maxEta = 2.5;     // Maximum eta (current eta cut for tracks = 2.4)
  const Int_t nEtaBins = 50;       // Number of eta bins
  
  // Vertex z-position
  const Double_t minVz = -20;   // Minimum vz
  const Double_t maxVz = 20;    // Maximum vz
  const Int_t nVzBins = 80;     // Number of vz bins
  
  // pT hat
  const Double_t minPtHat = 0;     // Minimum pT hat
  const Double_t maxPtHat = 460;   // Maximum pT hat
  const Int_t nFinePtHatBins = 230; // Number of fine pT hat bins
  
  // Generator level pT binning for closure histograms
  const Double_t minClosurePt = 50;                             // Minimum gen jet pT for closure plots
  const Double_t maxClosurePt = 500;                            // Maximum gen jet pT for closure plots
  const Int_t nClosurePtBins = (maxClosurePt-minClosurePt)/10;  // Bin width of 10 for the Gen pT in closure plots
  
  // Particle type for closure plots (0 = quark, 1 = gluon)
  const Double_t minClosureParticleType = -0.5;                        // Closure particle type indexing starts from zero
  const Double_t maxClosureParticleType = knClosureParticleTypes-0.5;  // Maximum closure particle type index
  const Int_t nClosureParticleTypeBins = knClosureParticleTypes;       // Bin width for particle type is 1
  
  // Track pairing type for energy-energy correlators (0 = same jet, 1 = signal-reflected cone 2 = only eta-reflected jet)
  const Double_t minTrackPairingType = -0.5;                // Track pairing type indexing starts from zero
  const Double_t maxTrackPairingType = knPairingTypes-0.5;  // Maximum track pairing type index
  const Int_t nTrackPairingTypeBins = knPairingTypes;       // Bin width for track pairing type is 1
  
  // Subevent index for particles (0 = Particle is from embedded Pythia, 1 = Particle is from Hydjet)
  const Double_t minSubeventType = -0.5;
  const Double_t maxSubeventType = knSubeventTypes-0.5;
  const Int_t nSubeventTypeBins = knSubeventTypes;
  
  // Subevent pairing index for particles (0 = both paritcles from Pythia, 1 = one particle from Pythia, one from Hydjet, 2 = Both particles from Hydjet)
  const Double_t minSubeventCombination = -0.5;
  const Double_t maxSubeventCombination = knSubeventCombinations-0.5;
  const Int_t nSubeventCombinationBins = knSubeventCombinations;
  
  // Jet cone type (0 = signal cone, 1 = reflected cone)
  const Double_t minJetConeType = -0.5;
  const Double_t maxJetConeType = knJetConeTypes-0.5;
  const Int_t nJetConeTypes = knJetConeTypes;
  
  // Binning for reco/gen ratio for closure histograms
  const Double_t minClosureRatio = 0.025;    // Minimum ratio for the closure plots
  const Double_t maxClosureRatio = 2.025;    // Maximum ratio for the closure plots
  const Int_t nClosureRatioBins = 40;    // Number of closure ratio bins
  
  // Binning for multiplicity
  const Double_t minMultiplicity = 0;
  const Double_t maxMultiplicity = 5000;
  const Double_t maxMultiplicityWeighted = 4000;
  const Int_t nMultiplicityBins = 500;
  const Int_t nMultiplicityBinsWeighted = 400;
  
  // Binning for multiplicity in within jet cone
  const Double_t minMultiplicityInJets = -0.5;
  const Double_t maxMultiplicityInJets = 149.5;
  const Int_t nMultiplicityInJetsBins = 150;
  
  // Centrality bins for THnSparses (We run into memory issues, if have all the bins)
  const Int_t nWideCentralityBins = fCard->GetNBin("CentralityBinEdges");
  Double_t wideCentralityBins[nWideCentralityBins+1];
  for(Int_t iCentrality = 0; iCentrality < nWideCentralityBins+1; iCentrality++){
    wideCentralityBins[iCentrality] = fCard->Get("CentralityBinEdges",iCentrality);
  }
  
  // Bins for the pT hat histogram
  const Int_t nPtHatBins = fCard->GetNBin("PtHatBinEdges");
  Double_t ptHatBins[nPtHatBins+1];
  for(Int_t iPtHat = 0; iPtHat < nPtHatBins+1; iPtHat++){
    ptHatBins[iPtHat] = fCard->Get("PtHatBinEdges",iPtHat);
  }
  
  // Jet pT binning for energy-energy correlator histograms
  const Int_t nJetPtBinsEEC = fCard->GetNBin("JetPtBinEdgesEEC");
  Double_t jetPtBinsEEC[nJetPtBinsEEC+1];
  for(Int_t iJetPt = 0; iJetPt < nJetPtBinsEEC+1; iJetPt++){
    jetPtBinsEEC[iJetPt] = fCard->Get("JetPtBinEdgesEEC",iJetPt);
  }
  const Double_t minJetPtEEC = jetPtBinsEEC[0];
  const Double_t maxJetPtEEC = jetPtBinsEEC[nJetPtBinsEEC];
  
  // Track pT binning for energy-energy correlator histograms
  const Int_t nTrackPtBinsEEC = fCard->GetNBin("TrackPtBinEdgesEEC");
  Double_t trackPtBinsEEC[nTrackPtBinsEEC+1];
  for(Int_t iTrackPt = 0; iTrackPt < nTrackPtBinsEEC+1; iTrackPt++){
    trackPtBinsEEC[iTrackPt] = fCard->Get("TrackPtBinEdgesEEC",iTrackPt);
  }
  const Double_t minTrackPtEEC = trackPtBinsEEC[0];
  const Double_t maxTrackPtEEC = trackPtBinsEEC[nTrackPtBinsEEC];
  
  // DeltaR binning for track density
  const Int_t nDeltaRBinsTrackDensity = 80;
  const Double_t minDeltaRTrackDensity = 0;
  const Double_t maxDeltaRTrackDensity = 0.8;
  
  // Logarithmic deltaR binning for energy-energy correlator histograms
  const Int_t nDeltaRBinsEEC = 40;
  const Double_t minDeltaREEC = 0;
  const Double_t maxDeltaREEC = 0.8;
  const Double_t binnerShift = 0.01;
  const Double_t deltaRlogBinWidth = (TMath::Log(maxDeltaREEC+binnerShift) - TMath::Log(minDeltaREEC+binnerShift)) / nDeltaRBinsEEC;
  Double_t deltaRBinsEEC[nDeltaRBinsEEC+1];
  for(int iDeltaR = 0; iDeltaR <= nDeltaRBinsEEC; iDeltaR++){
    deltaRBinsEEC[iDeltaR] = (minDeltaREEC+binnerShift)*TMath::Exp(iDeltaR*deltaRlogBinWidth)-binnerShift;
  }
  
  // Binning for the two-dimensional unfolding response
  const Int_t nUnfoldingBins = nDeltaRBinsEEC*nJetPtBinsEEC;
  const Double_t minUnfoldingBin = 0;
  const Double_t maxUnfoldingBin = nJetPtBinsEEC*maxDeltaREEC;
  Double_t fullUnfoldingBinning[nUnfoldingBins+1];
  for(int iDeltaR = 0; iDeltaR < nDeltaRBinsEEC; iDeltaR++){
    for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
      fullUnfoldingBinning[iDeltaR+iJetPt*nDeltaRBinsEEC] = deltaRBinsEEC[iDeltaR]+maxDeltaREEC*iJetPt;
    }
  }
  fullUnfoldingBinning[nUnfoldingBins] = maxUnfoldingBin;

  // Arrays for creating THnSparses
  const Int_t nAxesMultiplicity = 3;
  Int_t nBinsMultiplicity[nAxesMultiplicity];
  Double_t lowBinBorderMultiplicity[nAxesMultiplicity];
  Double_t highBinBorderMultiplicity[nAxesMultiplicity];
  
  const Int_t nAxesMultiplicityInJets = 5;
  Int_t nBinsMultiplicityInJets[nAxesMultiplicityInJets];
  Double_t lowBinBorderMultiplicityInJets[nAxesMultiplicityInJets];
  Double_t highBinBorderMultiplicityInJets[nAxesMultiplicityInJets];
  
  const Int_t nAxesJet = 5;
  Int_t nBinsJet[nAxesJet];
  Double_t lowBinBorderJet[nAxesJet];
  Double_t highBinBorderJet[nAxesJet];
  
  const Int_t nAxesTrack = 4;
  Int_t nBinsTrack[nAxesTrack];
  Double_t lowBinBorderTrack[nAxesTrack];
  Double_t highBinBorderTrack[nAxesTrack];
  
  const Int_t nAxesTrackDensity = 6;
  Int_t nBinsTrackDensity[nAxesTrackDensity];
  Double_t lowBinBorderTrackDensity[nAxesTrackDensity];
  Double_t highBinBorderTrackDensity[nAxesTrackDensity];
  
  const Int_t nAxesMaxParticlePtInJet = 4;
  Int_t nBinsMaxParticlePtInJet[nAxesMaxParticlePtInJet];
  Double_t lowBinBorderMaxParticlePtInJet[nAxesMaxParticlePtInJet];
  Double_t highBinBorderMaxParticlePtInJet[nAxesMaxParticlePtInJet];
  
  const Int_t nAxesEnergyEnergyCorrelator = 6;
  Int_t nBinsEnergyEnergyCorrelator[nAxesEnergyEnergyCorrelator];
  Double_t lowBinBorderEnergyEnergyCorrelator[nAxesEnergyEnergyCorrelator];
  Double_t highBinBorderEnergyEnergyCorrelator[nAxesEnergyEnergyCorrelator];
  
  const Int_t nAxesJetClosure = 6;
  Int_t nBinsJetClosure[nAxesJetClosure];
  Double_t lowBinBorderJetClosure[nAxesJetClosure];
  Double_t highBinBorderJetClosure[nAxesJetClosure];

  const Int_t nAxesJetPtUnfoldDistribution = 3;
  Int_t nBinsJetPtUnfoldDistribution[nAxesJetPtUnfoldDistribution];
  Double_t lowBinBorderJetPtUnfoldDistribution[nAxesJetPtUnfoldDistribution];
  Double_t highBinBorderJetPtUnfoldDistribution[nAxesJetPtUnfoldDistribution];

  const Int_t nAxesJetPtUnfoldResponse = 4;
  Int_t nBinsJetPtUnfoldResponse[nAxesJetPtUnfoldResponse];
  Double_t lowBinBorderJetPtUnfoldResponse[nAxesJetPtUnfoldResponse];
  Double_t highBinBorderJetPtUnfoldResponse[nAxesJetPtUnfoldResponse];
  
  
  // ======== Plain TH1 histograms ========
  
  fhVertexZ = new TH1F("vertexZ","vertexZ",nVzBins,minVz,maxVz); fhVertexZ->Sumw2();
  fhVertexZWeighted = new TH1F("vertexZweighted","vertexZweighted",nVzBins,minVz,maxVz); fhVertexZWeighted->Sumw2();
  fhEvents = new TH1F("nEvents","nEvents",knEventTypes,-0.5,knEventTypes-0.5); fhEvents->Sumw2();
  fhTriggers = new TH1F("triggers","triggers",knTriggerCombinations,-0.5,knTriggerCombinations-0.5); fhTriggers->Sumw2();
  fhTriggersAfterSelection = new TH1F("triggersAfterSelection","triggersAfterSelection",knTriggerCombinations,-0.5,knTriggerCombinations-0.5); fhTriggersAfterSelection->Sumw2();
  fhTrackCuts = new TH1F("trackCuts","trackCuts",knTrackCuts,-0.5,knTrackCuts-0.5); fhTrackCuts->Sumw2();
  fhCentrality = new TH1F("centrality","centrality",nCentralityBins,minCentrality,maxCentrality); fhCentrality->Sumw2();
  fhCentralityWeighted = new TH1F("centralityWeighted","centralityWeighted",nCentralityBins,minCentrality,maxCentrality); fhCentralityWeighted->Sumw2();
  fhPtHat = new TH1F("pthat","pthat",nPtHatBins,ptHatBins); fhPtHat->Sumw2();
  fhPtHatWeighted = new TH1F("pthatWeighted","pthatWeighted",nFinePtHatBins,minPtHat,maxPtHat); fhPtHatWeighted->Sumw2();
  
  // For the event histogram, label each bin corresponding to an event cut
  for(Int_t i = 0; i < knEventTypes; i++){
    fhEvents->GetXaxis()->SetBinLabel(i+1,kEventTypeStrings[i]);
  }
  
  // For the trigger histogram, label each bin corresponding to a trigger selection
  for(Int_t i = 0; i < knTriggerCombinations; i++){
    fhTriggers->GetXaxis()->SetBinLabel(i+1,kTriggerStrings[i]);
    fhTriggersAfterSelection->GetXaxis()->SetBinLabel(i+1,kTriggerStrings[i]);
  }
  
  // For the track cut histogram, label each bin corresponding to a track cut
  for(Int_t i = 0; i < knTrackCuts; i++){
    fhTrackCuts->GetXaxis()->SetBinLabel(i+1,kTrackCutStrings[i]);
  }
  
  // ======== THnSparses for multiplicity ========
  
  // Axis 0 for the multiplicity histogram: multiplicity
  nBinsMultiplicity[0] = nMultiplicityBins;       // nBins for multiplicity
  lowBinBorderMultiplicity[0] = minMultiplicity;  // low bin border for multiplicity
  highBinBorderMultiplicity[0] = maxMultiplicity; // high bin border for multiplicity
  
  // Axis 1 for the multiplicity histogram: weighted multiplicity
  nBinsMultiplicity[1] = nMultiplicityBinsWeighted;       // nBins for weighted multiplicity
  lowBinBorderMultiplicity[1] = minMultiplicity;          // low bin border for weighted multiplicity
  highBinBorderMultiplicity[1] = maxMultiplicityWeighted; // high bin border for weighted multiplicity
  
  // Axis 2 for the multiplicity histogram: centrality
  nBinsMultiplicity[2] = nCentralityBins;     // nBins for wide centrality bins
  lowBinBorderMultiplicity[2] = minCentrality;    // low bin border for centrality
  highBinBorderMultiplicity[2] = maxCentrality;   // high bin border for centrality
  
  // Create the histograms for leading and subleading jets using the above binning information
  fhMultiplicity = new THnSparseF("multiplicity", "multiplicity", nAxesMultiplicity, nBinsMultiplicity, lowBinBorderMultiplicity, highBinBorderMultiplicity); fhMultiplicity->Sumw2();
  
  // ======== THnSparse for multiplicity within the jet cone ========
  
  // Axis 0 for the multiplicity in jets histogram: multiplicity
  nBinsMultiplicityInJets[0] = nMultiplicityInJetsBins;       // nBins for non-corrected multiplicity
  lowBinBorderMultiplicityInJets[0] = minMultiplicityInJets;  // low bin border for non-corrected multiplicity
  highBinBorderMultiplicityInJets[0] = maxMultiplicityInJets; // high bin border for non-corrected multiplicity
  
  // Axis 1 for the multiplicity in jets histogram: jet pT
  nBinsMultiplicityInJets[1] = nJetPtBinsEEC;         // nBins for jet pT in energy-energy correlator histograms
  lowBinBorderMultiplicityInJets[1] = minJetPtEEC;    // low bin border for jet pT in energy-energy correlator histograms
  highBinBorderMultiplicityInJets[1] = maxJetPtEEC;   // high bin border for jet pT in energy-energy correlator histograms
  
  // Axis 2 for the multiplicity in jets histogram: track pT
  nBinsMultiplicityInJets[2] = nTrackPtBinsEEC;        // nBins for track pT in energy-energy correlator histograms
  lowBinBorderMultiplicityInJets[2] = minTrackPtEEC;   // low bin border for track pT in energy-energy correlator histograms
  highBinBorderMultiplicityInJets[2] = maxTrackPtEEC;  // high bin border for track pT in energy-energy correlator histograms
  
  // Axis 3 for the multiplicity in jets histogram: centrality
  nBinsMultiplicityInJets[3] = nWideCentralityBins;    // nBins for wide centrality bins
  lowBinBorderMultiplicityInJets[3] = minCentrality;   // low bin border for centrality
  highBinBorderMultiplicityInJets[3] = maxCentrality;  // high bin border for centrality
  
  // Axis 4 for the multiplicity in jets histogram: track subevent (only relevant for Monte Carlo)
  nBinsMultiplicityInJets[4] = nSubeventTypeBins+1;    // nBins for subevent types (Pythia/Hydjet/Combined)
  lowBinBorderMultiplicityInJets[4] = minSubeventType;       // low bin border for subevent types
  highBinBorderMultiplicityInJets[4] = maxSubeventType+1;    // high bin border for subevent types
  
  // Create the multiplicity in jets histrograms using the above binning information
  fhParticleMultiplicityInJet = new THnSparseF("multiplicityInJetCone", "multiplicityInJetCone", nAxesMultiplicityInJets, nBinsMultiplicityInJets, lowBinBorderMultiplicityInJets, highBinBorderMultiplicityInJets); fhParticleMultiplicityInJet->Sumw2();
  fhParticleMultiplicityInReflectedCone = new THnSparseF("multiplicityInReflectedCone", "multiplicityInReflectedCone", nAxesMultiplicityInJets, nBinsMultiplicityInJets, lowBinBorderMultiplicityInJets, highBinBorderMultiplicityInJets); fhParticleMultiplicityInReflectedCone->Sumw2();
  fhParticleMultiplicityInJetUncorrected = new THnSparseF("multiplicityInJetConeUncorrected", "multiplicityInJetConeUncorrected", nAxesMultiplicityInJets, nBinsMultiplicityInJets, lowBinBorderMultiplicityInJets, highBinBorderMultiplicityInJets); fhParticleMultiplicityInJetUncorrected->Sumw2();
  fhParticleMultiplicityInReflectedConeUncorrected = new THnSparseF("multiplicityInReflectedConeUncorrected", "multiplicityInReflectedConeUncorrected", nAxesMultiplicityInJets, nBinsMultiplicityInJets, lowBinBorderMultiplicityInJets, highBinBorderMultiplicityInJets); fhParticleMultiplicityInReflectedConeUncorrected->Sumw2();
  
  fhParticleMultiplicityInJet->SetBinEdges(1,jetPtBinsEEC);                              // Jet pT bins
  fhParticleMultiplicityInReflectedCone->SetBinEdges(1,jetPtBinsEEC);                    // Jet pT bins
  fhParticleMultiplicityInJetUncorrected->SetBinEdges(1,jetPtBinsEEC);                   // Jet pT bins
  fhParticleMultiplicityInReflectedConeUncorrected->SetBinEdges(1,jetPtBinsEEC);         // Jet pT bins
  
  fhParticleMultiplicityInJet->SetBinEdges(2,trackPtBinsEEC);                            // Track pT bins
  fhParticleMultiplicityInReflectedCone->SetBinEdges(2,trackPtBinsEEC);                  // Track pT bins
  fhParticleMultiplicityInJetUncorrected->SetBinEdges(2,trackPtBinsEEC);                 // Track pT bins
  fhParticleMultiplicityInReflectedConeUncorrected->SetBinEdges(2,trackPtBinsEEC);       // Track pT bins
  
  fhParticleMultiplicityInJet->SetBinEdges(3,wideCentralityBins);                        // Centrality bins
  fhParticleMultiplicityInReflectedCone->SetBinEdges(3,wideCentralityBins);              // Centrality bins
  fhParticleMultiplicityInJetUncorrected->SetBinEdges(3,wideCentralityBins);             // Centrality bins
  fhParticleMultiplicityInReflectedConeUncorrected->SetBinEdges(3,wideCentralityBins);   // Centrality bins
  
  // ======== THnSparse for all jets ========
  
  // Axis 0 for the any jet histogram: jet pT
  nBinsJet[0] = nPtBinsJet;         // nBins for any jet pT
  lowBinBorderJet[0] = minPtJet;    // low bin border for any jet pT
  highBinBorderJet[0] = maxPtJet;   // high bin border for any jet pT
  
  // Axis 1 for the any jet histogram: jet phi
  nBinsJet[1] = nPhiBins;        // nBins for any jet phi
  lowBinBorderJet[1] = minPhi;   // low bin border for any jet phi
  highBinBorderJet[1] = maxPhi;  // high bin border for any jet phi
  
  // Axis 2 for the any jet histogram: jet eta
  nBinsJet[2] = nEtaBins;        // nBins for any jet eta
  lowBinBorderJet[2] = minEta;   // low bin border for any jet eta
  highBinBorderJet[2] = maxEta;  // high bin border for any jet eta
  
  // Axis 3 for the any jet histogram: centrality
  nBinsJet[3] = nWideCentralityBins;   // nBins for wide centrality bins
  lowBinBorderJet[3] = minCentrality;  // low bin border for centrality
  highBinBorderJet[3] = maxCentrality; // high bin border for centrality
  
  // Axis 4 for the jet histogram: jet flavor (quark/gluon)
  nBinsJet[4] = nClosureParticleTypeBins;        // nBins for jet flavor
  lowBinBorderJet[4] = minClosureParticleType;   // low bin border for jet flavor
  highBinBorderJet[4] = maxClosureParticleType;  // high bin border for jet flavor
  
  // Create the histogram for all jets using the above binning information
  fhInclusiveJet = new THnSparseF("inclusiveJet","inclusiveJet",nAxesJet,nBinsJet,lowBinBorderJet,highBinBorderJet); fhInclusiveJet->Sumw2();

  // Set custom centrality bins for histograms
  fhInclusiveJet->SetBinEdges(3,wideCentralityBins);
  
  // ======== THnSparses for tracks and uncorrected tracks ========
  
  // Axis 0 for the track histogram: track pT
  nBinsTrack[0] = nPtBinsTrack;         // nBins for track pT
  lowBinBorderTrack[0] = minPtTrack;    // low bin border for track pT
  highBinBorderTrack[0] = maxPtTrack;   // high bin border for track pT
  
  // Axis 1 for the track histogram: track phi
  nBinsTrack[1] = nPhiBins;         // nBins for track phi
  lowBinBorderTrack[1] = minPhi;    // low bin border for track phi
  highBinBorderTrack[1] = maxPhi;   // high bin border for track phi
  
  // Axis 2 for the track histogram: track eta
  nBinsTrack[2] = nEtaBins;         // nBins for track eta
  lowBinBorderTrack[2] = minEta;    // low bin border for track eta
  highBinBorderTrack[2] = maxEta;   // high bin border for track eta
  
  // Axis 3 for the track histogram: centrality
  nBinsTrack[3] = nWideCentralityBins;   // nBins for wide centrality bins
  lowBinBorderTrack[3] = minCentrality;  // low bin border for centrality
  highBinBorderTrack[3] = maxCentrality; // high bin border for centrality
  
  // Create the histograms for tracks and uncorrected tracks using the above binning information
  fhTrack = new THnSparseF("track","track",nAxesTrack,nBinsTrack,lowBinBorderTrack,highBinBorderTrack); fhTrack->Sumw2();
  fhTrackUncorrected = new THnSparseF("trackUncorrected","trackUncorrected",nAxesTrack,nBinsTrack,lowBinBorderTrack,highBinBorderTrack); fhTrackUncorrected->Sumw2();

  // Set custom centrality bins for histograms
  fhTrack->SetBinEdges(3,wideCentralityBins);
  fhTrackUncorrected->SetBinEdges(3,wideCentralityBins);
  
  // ======== THnSparse for maximum particle pT within a jet ========
  
  // Axis 0 for the track density histogram: deltaR
  nBinsMaxParticlePtInJet[0] = nJetPtBinsEEC;         // nBins for jet pT
  lowBinBorderMaxParticlePtInJet[0] = minJetPtEEC;    // low bin border for jet pT
  highBinBorderMaxParticlePtInJet[0] = maxJetPtEEC;   // high bin border for jet pT
  
  // Axis 1 for the track density histogram: particle pT (signal)
  nBinsMaxParticlePtInJet[1] = nPtBinsParticle;        // nBins for max particle pT
  lowBinBorderMaxParticlePtInJet[1] = minParticlePt;   // low bin border for max particle pT
  highBinBorderMaxParticlePtInJet[1] = maxParticlePt;  // high bin border for max particle pT
  
  // Axis 2 for the track density histogram: particle pT (background)
  nBinsMaxParticlePtInJet[2] = nPtBinsParticle;        // nBins for max particle pT
  lowBinBorderMaxParticlePtInJet[2] = minParticlePt;   // low bin border for max particle pT
  highBinBorderMaxParticlePtInJet[2] = maxParticlePt;  // high bin border for max particle pT
  
  // Axis 3 for track density histogram: centrality
  nBinsMaxParticlePtInJet[3] = nWideCentralityBins;    // nBins for wide centrality bins
  lowBinBorderMaxParticlePtInJet[3] = minCentrality;   // low bin border for centrality
  highBinBorderMaxParticlePtInJet[3] = maxCentrality;  // high bin border for centrality
  
  // Create the histograms for track density and pT weighted track density using the above binning information
  fhMaxPtParticleInJet = new THnSparseF("maxParticlePtInJet", "maxParticlePtInJet", nAxesMaxParticlePtInJet, nBinsMaxParticlePtInJet, lowBinBorderMaxParticlePtInJet, highBinBorderMaxParticlePtInJet); fhMaxPtParticleInJet->Sumw2();
  
  // Set custom bin borders for histograms
  fhMaxPtParticleInJet->SetBinEdges(0,jetPtBinsEEC);              // Jet pT bins
  fhMaxPtParticleInJet->SetBinEdges(3,wideCentralityBins);        // Centrality bins
  
  // ======== THnSparse for track density around the jets ========
  
  // Axis 0 for the track density histogram: deltaR
  nBinsTrackDensity[0] = nDeltaRBinsTrackDensity;        // nBins for deltaR between the track and the jet axis
  lowBinBorderTrackDensity[0] = minDeltaRTrackDensity;   // low bin border for deltaR
  highBinBorderTrackDensity[0] = maxDeltaRTrackDensity;  // high bin border for deltaR
  
  // Axis 1 for the track density histogram: jet pT
  nBinsTrackDensity[1] = nJetPtBinsEEC;         // nBins for jet pT
  lowBinBorderTrackDensity[1] = minJetPtEEC;    // low bin border for jet pT
  highBinBorderTrackDensity[1] = maxJetPtEEC;   // high bin border for jet pT
  
  // Axis 2 for the track density histogram: track pT
  nBinsTrackDensity[2] = nTrackPtBinsEEC;        // nBins for track pT
  lowBinBorderTrackDensity[2] = minTrackPtEEC;   // low bin border for track pT
  highBinBorderTrackDensity[2] = maxTrackPtEEC;  // high bin border for track pT
  
  // Axis 3 for track density histogram: centrality
  nBinsTrackDensity[3] = nWideCentralityBins;    // nBins for wide centrality bins
  lowBinBorderTrackDensity[3] = minCentrality;   // low bin border for centrality
  highBinBorderTrackDensity[3] = maxCentrality;  // high bin border for centrality
  
  // Axis 4 for track density histogram: jet cone type (signal cone/reflected cone)
  nBinsTrackDensity[4] = nJetConeTypes;           // nBins for wide centrality bins
  lowBinBorderTrackDensity[4] = minJetConeType;   // low bin border for centrality
  highBinBorderTrackDensity[4] = maxJetConeType;  // high bin border for centrality
  
  // Axis 5 for track density histogram: track subevent (only relevant for Monte Carlo)
  nBinsTrackDensity[5] = nSubeventTypeBins;       // nBins for subevent types (Pythia/Hydjet)
  lowBinBorderTrackDensity[5] = minSubeventType;  // low bin border for subevent combinations
  highBinBorderTrackDensity[5] = maxSubeventType; // high bin border for subevent combinations
  
  // Create the histograms for track density and pT weighted track density using the above binning information
  fhParticleDensityAroundJet = new THnSparseF("particleDensity", "particleDensity", nAxesTrackDensity, nBinsTrackDensity ,lowBinBorderTrackDensity, highBinBorderTrackDensity); fhParticleDensityAroundJet->Sumw2();
  fhParticlePtDensityAroundJet = new THnSparseF("particlePtDensity", "particlePtDensity", nAxesTrackDensity, nBinsTrackDensity, lowBinBorderTrackDensity, highBinBorderTrackDensity); fhParticlePtDensityAroundJet->Sumw2();
  
  // Set custom bin borders for histograms
  fhParticleDensityAroundJet->SetBinEdges(1,jetPtBinsEEC);              // Jet pT bins
  fhParticlePtDensityAroundJet->SetBinEdges(1,jetPtBinsEEC);            // Jet pT bins
  
  fhParticleDensityAroundJet->SetBinEdges(2,trackPtBinsEEC);            // Track pT bins
  fhParticlePtDensityAroundJet->SetBinEdges(2,trackPtBinsEEC);          // Track pT bins
  
  fhParticleDensityAroundJet->SetBinEdges(3,wideCentralityBins);        // Centrality bins
  fhParticlePtDensityAroundJet->SetBinEdges(3,wideCentralityBins);      // Centrality bins
  
  // ======== THnSparses for energy-energy correlators ========
  
  // Axis 0 for the energy-energy correlator histogram: deltaR
  nBinsEnergyEnergyCorrelator[0] = nDeltaRBinsEEC;        // nBins for deltaR between the two tracks
  lowBinBorderEnergyEnergyCorrelator[0] = minDeltaREEC;   // low bin border for deltaR
  highBinBorderEnergyEnergyCorrelator[0] = maxDeltaREEC;  // high bin border for deltaR
  
  // Axis 1 for the energy-energy correlator histogram: jet pT
  nBinsEnergyEnergyCorrelator[1] = nJetPtBinsEEC;         // nBins for jet pT
  lowBinBorderEnergyEnergyCorrelator[1] = minJetPtEEC;    // low bin border for jet pT
  highBinBorderEnergyEnergyCorrelator[1] = maxJetPtEEC;   // high bin border for jet pT
  
  // Axis 2 for the energy-energy correlator histogram: track pT
  nBinsEnergyEnergyCorrelator[2] = nTrackPtBinsEEC;        // nBins for track pT
  lowBinBorderEnergyEnergyCorrelator[2] = minTrackPtEEC;   // low bin border for track pT
  highBinBorderEnergyEnergyCorrelator[2] = maxTrackPtEEC;  // high bin border for track pT
  
  // Axis 3 for the energy-energy correlator histogram: centrality
  nBinsEnergyEnergyCorrelator[3] = nWideCentralityBins;   // nBins for wide centrality bins
  lowBinBorderEnergyEnergyCorrelator[3] = minCentrality;  // low bin border for centrality
  highBinBorderEnergyEnergyCorrelator[3] = maxCentrality; // high bin border for centrality
  
  // Axis 4 for the energy-energy correlator histogram: track pairing type (same jet/reflected cone jet)
  nBinsEnergyEnergyCorrelator[4] = nTrackPairingTypeBins;       // nBins for track pairing types
  lowBinBorderEnergyEnergyCorrelator[4] = minTrackPairingType;  // low bin border for track pairing type
  highBinBorderEnergyEnergyCorrelator[4] = maxTrackPairingType; // high bin border for track pairing type
  
  // Axis 5 for the energy-energy correlator histogram: track subevent (only relevant for Monte Carlo)
  nBinsEnergyEnergyCorrelator[5] = nSubeventCombinationBins;       // nBins for subevent combinations (pythia-pythia, pythia-hydjet, hydjet-hydjet)
  lowBinBorderEnergyEnergyCorrelator[5] = minSubeventCombination;  // low bin border for subevent combinations
  highBinBorderEnergyEnergyCorrelator[5] = maxSubeventCombination; // high bin border for subevent combinations
  
  // Create the histograms for energy-energy correlators with and without track efficiency corrections
  fhEnergyEnergyCorrelator = new THnSparseF("energyEnergyCorrelator", "energyEnergyCorrelator", nAxesEnergyEnergyCorrelator, nBinsEnergyEnergyCorrelator, lowBinBorderEnergyEnergyCorrelator, highBinBorderEnergyEnergyCorrelator); fhEnergyEnergyCorrelator->Sumw2();
  fhEnergyEnergyCorrelatorUncorrected = new THnSparseF("energyEnergyCorrelatorUncorrected", "energyEnergyCorrelatorUncorrected", nAxesEnergyEnergyCorrelator, nBinsEnergyEnergyCorrelator, lowBinBorderEnergyEnergyCorrelator, highBinBorderEnergyEnergyCorrelator); fhEnergyEnergyCorrelatorUncorrected->Sumw2();
  fhEnergyEnergyCorrelatorJetPt = new THnSparseF("energyEnergyCorrelatorJetPt", "energyEnergyCorrelatorJetPt", nAxesEnergyEnergyCorrelator, nBinsEnergyEnergyCorrelator, lowBinBorderEnergyEnergyCorrelator, highBinBorderEnergyEnergyCorrelator); fhEnergyEnergyCorrelatorJetPt->Sumw2();
  fhEnergyEnergyCorrelatorJetPtUncorrected = new THnSparseF("energyEnergyCorrelatorJetPtUncorrected", "energyEnergyCorrelatorJetPtUncorrected", nAxesEnergyEnergyCorrelator, nBinsEnergyEnergyCorrelator, lowBinBorderEnergyEnergyCorrelator, highBinBorderEnergyEnergyCorrelator); fhEnergyEnergyCorrelatorJetPtUncorrected->Sumw2();
  
  // Set custom bin borders for histograms
  fhEnergyEnergyCorrelator->SetBinEdges(0,deltaRBinsEEC);                      // DeltaR bins
  fhEnergyEnergyCorrelatorUncorrected->SetBinEdges(0,deltaRBinsEEC);           // DeltaR bins
  fhEnergyEnergyCorrelatorJetPt->SetBinEdges(0,deltaRBinsEEC);                 // DeltaR bins
  fhEnergyEnergyCorrelatorJetPtUncorrected->SetBinEdges(0,deltaRBinsEEC);      // DeltaR bins
  
  fhEnergyEnergyCorrelator->SetBinEdges(1,jetPtBinsEEC);                       // Jet pT bins
  fhEnergyEnergyCorrelatorUncorrected->SetBinEdges(1,jetPtBinsEEC);            // Jet pT bins
  fhEnergyEnergyCorrelatorJetPt->SetBinEdges(1,jetPtBinsEEC);                  // Jet pT bins
  fhEnergyEnergyCorrelatorJetPtUncorrected->SetBinEdges(1,jetPtBinsEEC);       // Jet pT bins
  
  fhEnergyEnergyCorrelator->SetBinEdges(2,trackPtBinsEEC);                     // Track pT bins
  fhEnergyEnergyCorrelatorUncorrected->SetBinEdges(2,trackPtBinsEEC);          // Track pT bins
  fhEnergyEnergyCorrelatorJetPt->SetBinEdges(2,trackPtBinsEEC);                // Track pT bins
  fhEnergyEnergyCorrelatorJetPtUncorrected->SetBinEdges(2,trackPtBinsEEC);     // Track pT bins
  
  fhEnergyEnergyCorrelator->SetBinEdges(3,wideCentralityBins);                 // Centrality bins
  fhEnergyEnergyCorrelatorUncorrected->SetBinEdges(3,wideCentralityBins);      // Centrality bins
  fhEnergyEnergyCorrelatorJetPt->SetBinEdges(3,wideCentralityBins);            // Centrality bins
  fhEnergyEnergyCorrelatorJetPtUncorrected->SetBinEdges(3,wideCentralityBins); // Centrality bins
  
  // ======== THnSparses for jet pT closures ========
  
  // Axis 0 for the jet pT closure histogram: generator level jet pT
  nBinsJetClosure[0] = nClosurePtBins;       // nBins for generator level pT bins in closure plots
  lowBinBorderJetClosure[0] = minClosurePt;  // low bin border generator level pT in closure plots
  highBinBorderJetClosure[0] = maxClosurePt; // high bin border generator level pT in closure plots
  
  // Axis 1 for the jet pT closure histogram: reconstructed jet pT
  nBinsJetClosure[1] = nClosurePtBins;       // nBins for reconstructed jet pT bins in closure plots
  lowBinBorderJetClosure[1] = minClosurePt;  // low bin border for reconstructed jet pT in closure plots
  highBinBorderJetClosure[1] = maxClosurePt; // high bin border for reconstructed jet pT in closure plots
  
  // Axis 2 for the jet pT closure histogram: generator level jet eta
  nBinsJetClosure[2] = nEtaBins;             // nBins for jet eta
  lowBinBorderJetClosure[2] = minEta;        // low bin border for jet eta
  highBinBorderJetClosure[2] = maxEta;       // high bin border for jet eta
  
  // Axis 3 for the jet pT closure histogram: centrality
  nBinsJetClosure[3] = nWideCentralityBins;     // nBins for centrality
  lowBinBorderJetClosure[3] = minCentrality;    // low bin border for centrality
  highBinBorderJetClosure[3] = maxCentrality;   // high bin border for centrality
  
  // Axis 4 for the jet pT closure histogram: ref parton = quark/gluon
  nBinsJetClosure[4] = nClosureParticleTypeBins;         // nBins for reference parton
  lowBinBorderJetClosure[4] = minClosureParticleType;    // low bin border for reference parton
  highBinBorderJetClosure[4] = maxClosureParticleType;   // high bin border for reference parton
  
  // Axis 5 for the jet pT closure histogram: reco/gen ratio for closure
  nBinsJetClosure[5] = nClosureRatioBins;        // nBins for closure ratio
  lowBinBorderJetClosure[5] = minClosureRatio;   // low bin border for closure ratio
  highBinBorderJetClosure[5] = maxClosureRatio;  // high bin border for closure ratio
  
  // Create histograms for jet pT closure
  fhJetPtClosure = new THnSparseF("jetPtClosure", "jetPtClosure", nAxesJetClosure, nBinsJetClosure, lowBinBorderJetClosure, highBinBorderJetClosure); fhJetPtClosure->Sumw2();
  
  // Set custom centrality bins for histograms
  fhJetPtClosure->SetBinEdges(3,wideCentralityBins);
  
  // ======== RooUnfold compatible energy-energy correlator distributions. Assume only jet pT is unfolded ========

  // Axis 0 for the jet pT unfolding distribution histogram: deltaR in jet pT bins
  nBinsJetPtUnfoldDistribution[0] = nUnfoldingBins;          // nBins for the unfolded quantities in reconstructed level
  lowBinBorderJetPtUnfoldDistribution[0] = minUnfoldingBin;  // low bin border for the unfolded quantities in reconstructed level
  highBinBorderJetPtUnfoldDistribution[0] = maxUnfoldingBin; // high bin border for the unfolded quantities in reconstructed level

  // Axis 1 for the jet pT unfolding distribution histogram: track pT
  nBinsJetPtUnfoldDistribution[1] = nTrackPtBinsEEC;        // nBins for track pT in energy-energy correlator histograms
  lowBinBorderJetPtUnfoldDistribution[1] = minTrackPtEEC;   // low bin border for track pT in energy-energy correlator histograms
  highBinBorderJetPtUnfoldDistribution[1] = maxTrackPtEEC;  // high bin border for track pT in energy-energy correlator histograms

  // Axis 2 for the jet pT unfolding distribution histogram: centrality
  nBinsJetPtUnfoldDistribution[2] = nWideCentralityBins;     // nBins for centrality
  lowBinBorderJetPtUnfoldDistribution[2] = minCentrality;    // low bin border for centrality
  highBinBorderJetPtUnfoldDistribution[2] = maxCentrality;   // high bin border for centrality

  // Create the histogram for jet pT unfolding distributions
  fhUnfoldingMeasured = new THnSparseF("jetPtUnfoldingMeasured", "jetPtUnfoldingMeasured", nAxesJetPtUnfoldDistribution, nBinsJetPtUnfoldDistribution, lowBinBorderJetPtUnfoldDistribution, highBinBorderJetPtUnfoldDistribution); fhUnfoldingMeasured->Sumw2();
  fhUnfoldingTruth = new THnSparseF("jetPtUnfoldingTruth", "jetPtUnfoldingTruth", nAxesJetPtUnfoldDistribution, nBinsJetPtUnfoldDistribution, lowBinBorderJetPtUnfoldDistribution, highBinBorderJetPtUnfoldDistribution); fhUnfoldingTruth->Sumw2();

  // Set custom bin axes for the histograms
  fhUnfoldingMeasured->SetBinEdges(0, fullUnfoldingBinning);
  fhUnfoldingTruth->SetBinEdges(0, fullUnfoldingBinning);
  fhUnfoldingMeasured->SetBinEdges(1, trackPtBinsEEC);
  fhUnfoldingTruth->SetBinEdges(1, trackPtBinsEEC);
  fhUnfoldingMeasured->SetBinEdges(2, wideCentralityBins);
  fhUnfoldingTruth->SetBinEdges(2, wideCentralityBins);

  // ======== RooUnfold compatible response matrix. Assume only jet pT is unfolded ========

  // Axis 0 for the jet pT unfolding response histogram: deltaR in generator level jet pT bins
  nBinsJetPtUnfoldResponse[0] = nUnfoldingBins;          // nBins for the unfolded quantities in reconstructed level
  lowBinBorderJetPtUnfoldResponse[0] = minUnfoldingBin;  // low bin border for the unfolded quantities in reconstructed level
  highBinBorderJetPtUnfoldResponse[0] = maxUnfoldingBin; // high bin border for the unfolded quantities in reconstructed level
  
  // Axis 1 for the jet pT unfolding response histogram: deltaR in reconstructed jet pT bins
  nBinsJetPtUnfoldResponse[1] = nUnfoldingBins;          // nBins for the unfolded quantities in generator level
  lowBinBorderJetPtUnfoldResponse[1] = minUnfoldingBin;  // low bin border for the unfolded quantities in generator level
  highBinBorderJetPtUnfoldResponse[1] = maxUnfoldingBin; // high bin border for the unfolded quantities in generator level
  
  // Axis 2 for the jet pT unfolding response histogram: track pT
  nBinsJetPtUnfoldResponse[2] = nTrackPtBinsEEC;        // nBins for track pT in energy-energy correlator histograms
  lowBinBorderJetPtUnfoldResponse[2] = minTrackPtEEC;   // low bin border for track pT in energy-energy correlator histograms
  highBinBorderJetPtUnfoldResponse[2] = maxTrackPtEEC;  // high bin border for track pT in energy-energy correlator histograms

  // Axis 3 for the jet pT unfolding response histogram: centrality
  nBinsJetPtUnfoldResponse[3] = nWideCentralityBins;     // nBins for centrality
  lowBinBorderJetPtUnfoldResponse[3] = minCentrality;    // low bin border for centrality
  highBinBorderJetPtUnfoldResponse[3] = maxCentrality;   // high bin border for centrality

  // Create the histogram for jet pT unfolding response
  fhUnfoldingResponse = new THnSparseF("jetPtUnfoldingResponse", "jetPtUnfoldingResponse", nAxesJetPtUnfoldResponse, nBinsJetPtUnfoldResponse, lowBinBorderJetPtUnfoldResponse, highBinBorderJetPtUnfoldResponse); fhUnfoldingResponse->Sumw2();

  // Set custom bin axes for the histograms
  fhUnfoldingResponse->SetBinEdges(0, fullUnfoldingBinning);
  fhUnfoldingResponse->SetBinEdges(1, fullUnfoldingBinning);
  fhUnfoldingResponse->SetBinEdges(2, trackPtBinsEEC);
  fhUnfoldingResponse->SetBinEdges(3, wideCentralityBins);

}

/*
 * Write the histograms to file
 */
void EECHistograms::Write() const{
  
  // Write the histograms to file
  fhVertexZ->Write();
  fhVertexZWeighted->Write();
  fhEvents->Write();
  fhTriggers->Write();
  fhTriggersAfterSelection->Write();
  fhTrackCuts->Write();
  fhCentrality->Write();
  fhCentralityWeighted->Write();
  fhPtHat->Write();
  fhPtHatWeighted->Write();
  fhMultiplicity->Write();
  fhInclusiveJet->Write();
  fhTrack->Write();
  fhTrackUncorrected->Write();
  fhParticleDensityAroundJet->Write();
  fhParticlePtDensityAroundJet->Write();
  fhParticleMultiplicityInJet->Write();
  fhParticleMultiplicityInReflectedCone->Write();
  fhParticleMultiplicityInJetUncorrected->Write();
  fhParticleMultiplicityInReflectedConeUncorrected->Write();
  fhMaxPtParticleInJet->Write();
  fhEnergyEnergyCorrelator->Write();
  fhEnergyEnergyCorrelatorUncorrected->Write();
  fhEnergyEnergyCorrelatorJetPt->Write();
  fhEnergyEnergyCorrelatorJetPtUncorrected->Write();
  fhJetPtClosure->Write();
  fhUnfoldingMeasured->Write();
  fhUnfoldingTruth->Write();
  fhUnfoldingResponse->Write();
}

/*
 * Write the histograms to a given file
 */
void EECHistograms::Write(TString outputFileName) const{
  
  // Define the output file
  TFile *outputFile = new TFile(outputFileName, "RECREATE");
  
  // Write the histograms to file
  Write();
  
  // Close the file after everything is written
  outputFile->Close();
  
  // Delete the outputFile object
  delete outputFile;
}


