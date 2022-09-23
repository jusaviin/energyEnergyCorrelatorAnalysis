/*
 * Implementation of EECHistogramManager
 */

// Root includes
#include <THnSparse.h>
#include <TPad.h>

// Own includes
#include "EECHistogramManager.h"

/*
 * Default constructor
 */
EECHistogramManager::EECHistogramManager() :
  fInputFile(NULL),
  fCard(NULL),
  fSystemAndEnergy(""),
  fCompactSystemAndEnergy(""),
  fLoadEventInformation(false),
  fLoadJets(false),
  fLoad2DHistograms(false),
  fLoadJetPtClosureHistograms(false),
  fJetFlavor(0),
  fFirstLoadedCentralityBin(0),
  fLastLoadedCentralityBin(1),
  fFirstLoadedTrackPtBin(0),
  fLastLoadedTrackPtBin(1),
  fFirstLoadedJetPtBinEEC(0),
  fLastLoadedJetPtBinEEC(1),
  fFirstLoadedTrackPtBinEEC(0),
  fLastLoadedTrackPtBinEEC(1),
  fnCentralityBins(kMaxCentralityBins),
  fnTrackPtBins(kMaxTrackPtBins),
  fnJetPtBinsEEC(kMaxJetPtBinsEEC),
  fnTrackPtBinsEEC(kMaxTrackPtBinsEEC)
{
  
  // Do not draw anything by default
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    fLoadTracks[iTrackType] = false;
  }
  
  // Do not draw energy-energy correlator histograms for default
  for(int iEnergyEnergyCorrelatorType = 0; iEnergyEnergyCorrelatorType < knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelatorType++){
    fLoadEnergyEnergyCorrelatorHistograms[iEnergyEnergyCorrelatorType] = false;
  }
  
  // Default binning for centrality
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins + 1; iCentrality++){
    fCentralityBinIndices[iCentrality] = iCentrality+1;
    fCentralityBinBorders[iCentrality] = 0;
  }
  
  // Default binning for track pT
  for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins + 1; iTrackPt++){
    fTrackPtBinIndices[iTrackPt] = iTrackPt + 1;
    fTrackPtBinBorders[iTrackPt] = 0;
  }
  
  // Default binning for jet pT in energy-energy correlator histograms
  for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC + 1; iJetPt++){
    fJetPtIndicesEEC[iJetPt] = iJetPt+1;
    fJetPtBinBordersEEC[iJetPt] = 0;
  }
  
  // Default binning for track pT in energy-energy correlator histograms
  for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBinsEEC + 1; iTrackPt++){
    fTrackPtIndicesEEC[iTrackPt] = iTrackPt+1;
    fTrackPtBinBordersEEC[iTrackPt] = 0;
  }
  
  // Initialize all the other histograms to null
  fhVertexZ = NULL;            // Vertex z position
  fhVertexZWeighted = NULL;    // Weighted vertex z-position (only meaningfull for MC)
  fhEvents = NULL;             // Number of events surviving different event cuts
  fhTrackCuts = NULL;          // Number of tracks surviving different track cuts
  fhCentrality = NULL;         // Centrality of all events
  fhCentralityWeighted = NULL; // Weighted centrality distribution in all events (only meaningful for MC)
  fhPtHat = NULL;              // pT hat for MC events (only meaningful for MC)
  fhPtHatWeighted = NULL;      // Weighted pT hat distribution (only meaningful for MC)
  
  fhMultiplicityMap = NULL;              // Multiplicity vs. centrality map
  fhMultiplicityMapWeighted = NULL;      // Efficiency weighted multiplicity vs. centrality map
  
  // Centrality loop
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    
    fhMultiplicity[iCentrality] = NULL;               // Track multiplicity from all events
    fhMultiplicityWeighted[iCentrality] = NULL;       // Efficiency weighted track multiplicity from all events
    
    // Jet histograms
    fhJetPt[iCentrality] = NULL;      // Jet pT histograms
    fhJetPhi[iCentrality] = NULL;     // Jet phi histograms
    fhJetEta[iCentrality] = NULL;     // Jet eta histograms
    fhJetEtaPhi[iCentrality] = NULL;  // 2D eta-phi histogram for jets
    
    // Track histograms
    for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
      fhTrackPt[iTrackType][iCentrality] = NULL;   // Track pT histograms
      
      // Loop over track pT bins
      for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins + 1; iTrackPt++){
        fhTrackPhi[iTrackType][iCentrality][iTrackPt] = NULL;    // Track phi histograms
        fhTrackEta[iTrackType][iCentrality][iTrackPt] = NULL;    // Track eta histograms
        fhTrackEtaPhi[iTrackType][iCentrality][iTrackPt] = NULL; // 2D eta-phi histogram for track
      } // Track pT loop
      
    } // Track category loop
    
    // Energy-energy correlator histograms
    for(int iEnergyEnergyCorrelatorType = 0; iEnergyEnergyCorrelatorType < knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelatorType++){
      for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC; iJetPt++){
        for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBinsEEC; iTrackPt++){
          for(int iPairingType = 0; iPairingType < EECHistograms::knPairingTypes; iPairingType++){
            for(int iSubevent = 0; iSubevent < knSubeventTypes+1; iSubevent++){
              fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent] = NULL;
            } // Subevent loop
          } // Pairing type loop (same jet/reflected cone)
        } // Track pT bins for energy-energy correlators
      } // Jet pT bins for energy-energy correlators
    } // Energy-energy correlator type loop
    
    // Jet pT closure histograms
    for(int iGenJetPt = 0; iGenJetPt <= knGenJetPtBins; iGenJetPt++){
      for(int iJetEta = 0; iJetEta <= knJetEtaBins; iJetEta++){
        for(int iClosureParticle = 0; iClosureParticle < EECHistograms::knClosureParticleTypes+1; iClosureParticle++){
          fhJetPtClosure[iGenJetPt][iJetEta][iCentrality][iClosureParticle] = NULL;
        } // Closure particle loop
      } // Jet eta bin loop
    } // Gen jet pT loop
  } // Centrality loop
}

/*
 * Constructor
 */
EECHistogramManager::EECHistogramManager(TFile *inputFile) :
  EECHistogramManager()
{
  fInputFile = inputFile;
  
  // Read card from inputfile
  fCard = new EECCard(inputFile);
  
  // Initialize values using the information in card
  InitializeFromCard();
  
}

/*
 * Constructor
 */
EECHistogramManager::EECHistogramManager(TFile *inputFile, EECCard *card) :
  EECHistogramManager()
{
  fInputFile = inputFile;
  
  // Initialize values using the information in card
  fCard = card;
  InitializeFromCard();
  
}

/*
 * Initialize several member variables from EECCard
 */
void EECHistogramManager::InitializeFromCard(){
  
  // Read the collision system from the card
  TString collisionSystem = fCard->GetDataType();
  
  // Make a string for collision system based on information on the card
  fSystemAndEnergy = Form("%s 5.02 TeV",collisionSystem.Data());
  fCompactSystemAndEnergy = fSystemAndEnergy;
  fCompactSystemAndEnergy.ReplaceAll(" ","");
  fCompactSystemAndEnergy.ReplaceAll(".","v");
  
  // Read bins for centrality, track pT, jet pT for energy-energy correlators, and track pT for energy-energy correlators from the card
  fnCentralityBins = fCard->GetNCentralityBins();
  fnTrackPtBins = fCard->GetNTrackPtBins();
  fnJetPtBinsEEC = fCard->GetNJetPtBinsEEC();
  fnTrackPtBinsEEC = fCard->GetNTrackPtBinsEEC();
  
  // Centrality binning
  for(int iCentrality = 0; iCentrality <= fnCentralityBins; iCentrality++){
    fCentralityBinBorders[iCentrality] = fCard->GetLowBinBorderCentrality(iCentrality);
  }
  fLastLoadedCentralityBin = fnCentralityBins-1;
  
  // Track pT binning
  for(int iTrackPt = 0; iTrackPt <= fnTrackPtBins; iTrackPt++){
    fTrackPtBinBorders[iTrackPt] = fCard->GetLowBinBorderTrackPt(iTrackPt);
  }
  fLastLoadedTrackPtBin = fnTrackPtBins-1;
  
  // Jet pT binning for energy-energy correlators
  for(int iJetPt = 0; iJetPt <= fnJetPtBinsEEC; iJetPt++){
    fJetPtBinBordersEEC[iJetPt] = fCard->GetLowBinBorderJetPtEEC(iJetPt);
  }
  fLastLoadedJetPtBinEEC = fnJetPtBinsEEC-1;
  
  // Track pT binnins for energy-energy correlators
  for(int iTrackPt = 0; iTrackPt <= fnTrackPtBinsEEC; iTrackPt++){
    fTrackPtBinBordersEEC[iTrackPt] = fCard->GetLowBinBorderTrackPtEEC(iTrackPt);
  }
  fLastLoadedTrackPtBinEEC = fnTrackPtBinsEEC-1;
  
  // Remove centrality selection from pp data and local testing
  if(collisionSystem.Contains("pp") || collisionSystem.Contains("localTest")){
    fLastLoadedCentralityBin = 0;
    fCentralityBinBorders[0] = -0.5;
  }
  
}

/*
 * Copy constructor
 */
EECHistogramManager::EECHistogramManager(const EECHistogramManager& in) :
  fInputFile(in.fInputFile),
  fCard(in.fCard),
  fSystemAndEnergy(in.fSystemAndEnergy),
  fCompactSystemAndEnergy(in.fCompactSystemAndEnergy),
  fLoadEventInformation(in.fLoadEventInformation),
  fLoadJets(in.fLoadJets),
  fLoad2DHistograms(in.fLoad2DHistograms),
  fLoadJetPtClosureHistograms(in.fLoadJetPtClosureHistograms),
  fJetFlavor(in.fJetFlavor),
  fFirstLoadedCentralityBin(in.fFirstLoadedCentralityBin),
  fLastLoadedCentralityBin(in.fLastLoadedCentralityBin),
  fFirstLoadedTrackPtBin(in.fFirstLoadedTrackPtBin),
  fLastLoadedTrackPtBin(in.fLastLoadedTrackPtBin),
  fFirstLoadedJetPtBinEEC(in.fFirstLoadedJetPtBinEEC),
  fLastLoadedJetPtBinEEC(in.fLastLoadedJetPtBinEEC),
  fFirstLoadedTrackPtBinEEC(in.fFirstLoadedTrackPtBinEEC),
  fLastLoadedTrackPtBinEEC(in.fLastLoadedTrackPtBinEEC),
  fhVertexZ(in.fhVertexZ),
  fhVertexZWeighted(in.fhVertexZWeighted),
  fhEvents(in.fhEvents),
  fhTrackCuts(in.fhTrackCuts),
  fhCentrality(in.fhCentrality),
  fhCentralityWeighted(in.fhCentralityWeighted)
{
  // Copy constructor
  
  // Copy all values
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    fLoadTracks[iTrackType] = in.fLoadTracks[iTrackType];
  }
  
  // Copy energy-energy correlator histograms
  for(int iEnergyEnergyCorrelatorType = 0; iEnergyEnergyCorrelatorType < knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelatorType++){
    fLoadEnergyEnergyCorrelatorHistograms[iEnergyEnergyCorrelatorType] = in.fLoadEnergyEnergyCorrelatorHistograms[iEnergyEnergyCorrelatorType];
  }
  
  // Copy binning for centrality
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins + 1; iCentrality++){
    fCentralityBinIndices[iCentrality] = in.fCentralityBinIndices[iCentrality];
    fCentralityBinBorders[iCentrality] = in.fCentralityBinBorders[iCentrality];
  }
  
  // Copy binning for track pT
  for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins + 1; iTrackPt++){
    fTrackPtBinIndices[iTrackPt] = in.fTrackPtBinIndices[iTrackPt];
    fTrackPtBinBorders[iTrackPt] = in.fTrackPtBinBorders[iTrackPt];
  }
  
  // Copy binning for jet pT in energy-energy correlator histograms
  for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC+1; iJetPt++){
    fJetPtIndicesEEC[iJetPt] = in.fJetPtIndicesEEC[iJetPt];
    fJetPtBinBordersEEC[iJetPt] = in.fJetPtBinBordersEEC[iJetPt];
  }
  
  // Copy binning for track pT in energy-energy correlator histograms
  for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBinsEEC+1; iTrackPt++){
    fTrackPtIndicesEEC[iTrackPt] = in.fTrackPtIndicesEEC[iTrackPt];
    fTrackPtBinBordersEEC[iTrackPt] = in.fTrackPtBinBordersEEC[iTrackPt];
  }
  
  // Multiplicity maps
  fhMultiplicityMap = in.fhMultiplicityMap;
  fhMultiplicityMapWeighted = in.fhMultiplicityMapWeighted;
  
  // Centrality loop
  for(int iCentrality = 0; iCentrality < kMaxCentralityBins; iCentrality++){
    
    fhMultiplicity[iCentrality] = in.fhMultiplicity[iCentrality];                           // Track multiplicity from all events
    fhMultiplicityWeighted[iCentrality] = in.fhMultiplicityWeighted[iCentrality];           // Efficiency weighted track multiplicity from all events
    
    // Jet histograms
    fhJetPt[iCentrality] = in.fhJetPt[iCentrality];         // Jet pT histograms
    fhJetPhi[iCentrality] = in.fhJetPhi[iCentrality];       // Jet phi histograms
    fhJetEta[iCentrality] = in.fhJetEta[iCentrality];       // Jet eta histograms
    fhJetEtaPhi[iCentrality] = in.fhJetEtaPhi[iCentrality]; // 2D eta-phi histogram for jets

    
    // Track histograms
    for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
      fhTrackPt[iTrackType][iCentrality] = in.fhTrackPt[iTrackType][iCentrality];   // Track pT histograms
      
      // Loop over track pT bins
      for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBins + 1; iTrackPt++){
        fhTrackPhi[iTrackType][iCentrality][iTrackPt] = in.fhTrackPhi[iTrackType][iCentrality][iTrackPt];    // Track phi histograms
        fhTrackEta[iTrackType][iCentrality][iTrackPt] = in.fhTrackEta[iTrackType][iCentrality][iTrackPt];    // Track eta histograms
        fhTrackEtaPhi[iTrackType][iCentrality][iTrackPt] = in.fhTrackEtaPhi[iTrackType][iCentrality][iTrackPt]; // 2D eta-phi histogram for track
      } // Track pT loop
      
    } // Track category loop
    
    // Energy-energy correlator histograms
    for(int iEnergyEnergyCorrelatorType = 0; iEnergyEnergyCorrelatorType < knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelatorType++){
      for(int iJetPt = 0; iJetPt < kMaxJetPtBinsEEC; iJetPt++){
        for(int iTrackPt = 0; iTrackPt < kMaxTrackPtBinsEEC; iTrackPt++){
          for(int iPairingType = 0; iPairingType < EECHistograms::knPairingTypes; iPairingType++){
            for(int iSubevent = 0; iSubevent < knSubeventTypes+1; iSubevent++){
              fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent] = in.fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent];
            } // Subevent loop
          } // Pairing type loop
        } // Track pT bins for energy-energy correlators
      } // Jet pT bins for energy-energy correlators
    } // Energy-energy correlator type loop
    
    // Jet pT closure histograms
    for(int iGenJetPt = 0; iGenJetPt <= knGenJetPtBins; iGenJetPt++){
      for(int iJetEta = 0; iJetEta <= knJetEtaBins; iJetEta++){
        for(int iClosureParticle = 0; iClosureParticle < EECHistograms::knClosureParticleTypes+1; iClosureParticle++){
          fhJetPtClosure[iGenJetPt][iJetEta][iCentrality][iClosureParticle] = in.fhJetPtClosure[iGenJetPt][iJetEta][iCentrality][iClosureParticle];
        } // Closure particle loop
      } // Jet eta bin loop
    } // Gen jet pT loop
    
  } // Centrality loop
}

/*
 * Destructor
 */
EECHistogramManager::~EECHistogramManager(){
  delete fCard;
}

/*
 * Load all the selected histograms from the inputfile
 */
void EECHistogramManager::LoadHistograms(){
  
  // Always load the number of events histogram
  fhEvents = (TH1D*) fInputFile->Get("nEvents");                           // Number of events surviving different event cuts
  
  // Load the event information histograms
  if(fLoadEventInformation){
    fhVertexZ = (TH1D*) fInputFile->Get("vertexZ");                        // Vertex z position
    fhVertexZWeighted = (TH1D*) fInputFile->Get("vertexZweighted");        // MC weighted vertex z position
    fhTrackCuts = (TH1D*) fInputFile->Get("trackCuts");                    // Number of tracks surviving different track cuts
    fhCentrality = (TH1D*) fInputFile->Get("centrality");                  // Centrality in all events
    fhCentralityWeighted = (TH1D*) fInputFile->Get("centralityWeighted");  // MC weighted centrality in all events
    fhPtHat = (TH1D*) fInputFile->Get("pthat");                            // pT hat for MC events
    fhPtHatWeighted = (TH1D*) fInputFile->Get("pthatWeighted");            // Weighted pT hat for MC events
    
    LoadMultiplicityHistograms();
  }
  
  // Load jet histograms
  LoadJetHistograms();
  
  // Load track histograms
  LoadTrackHistograms();
  
  // Load energy-energy correlator histograms
  LoadEnergyEnergyCorrelatorHistograms();
  
  // Load jet pT closure histograms
  LoadJetPtClosureHistograms();
  
}

/*
 * Loader for multiplicity histograms
 *
 * THnSparse for multiplicity:
 *
 *     Histogram name        Axis index            Content of axis
 * --------------------------------------------------------------------------
 *      multiplicity           Axis 0             Track multiplicity
 *      multiplicity           Axis 1       Efficiency weighted multiplicity
 *      multiplicity           Axis 2                 Centrality
 */
void EECHistogramManager::LoadMultiplicityHistograms(){
    
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  
  // Multiplicity histograms have more denser centrality binning than other histograms, so need to determine them separately
  TH1D* hBinner;
  hBinner = FindHistogram(fInputFile,"multiplicity",2,0,0,0);
  double centralityBinIndicesMultiplicity[fnCentralityBins+1];
  for(int iBin = 0; iBin < fnCentralityBins+1; iBin++){
    centralityBinIndicesMultiplicity[iBin] = hBinner->GetXaxis()->FindBin(fCentralityBinBorders[iBin]);
  }
  
  for(int iCentralityBin = 0; iCentralityBin <= fnCentralityBins; iCentralityBin++){
    
    // Select the centrality bin indices
    if(iCentralityBin < fnCentralityBins){
      lowerCentralityBin = centralityBinIndicesMultiplicity[iCentralityBin];
      higherCentralityBin = centralityBinIndicesMultiplicity[iCentralityBin+1]+duplicateRemoverCentrality;
    } else {
      // Load the histograms without centrality selection to the index fnCentralityBins
      lowerCentralityBin = centralityBinIndicesMultiplicity[0];
      higherCentralityBin = centralityBinIndicesMultiplicity[fnCentralityBins]+duplicateRemoverCentrality;
    }
    
    fhMultiplicity[iCentralityBin] = FindHistogram(fInputFile, "multiplicity", 0, 2, lowerCentralityBin, higherCentralityBin);
    fhMultiplicityWeighted[iCentralityBin] = FindHistogram(fInputFile, "multiplicity", 1, 2, lowerCentralityBin, higherCentralityBin);
    
  } // Centrality loop
    
  fhMultiplicityMap = FindHistogram2D(fInputFile, "multiplicity", 0, 2, 1, 0, 0, 0);
  fhMultiplicityMapWeighted = FindHistogram2D(fInputFile, "multiplicity", 1, 2, 0, 0, 0, 0);
}


/*
 * Loader for jet histograms
 *
 * THnSparse for jets:
 *
 *   Histogram name: inclusiveJet
 *
 *     Axis index           Content of axis
 * --------------------------------------------------------
 *       Axis 0                 Jet pT
 *       Axis 1                 Jet phi
 *       Axis 2                 Jet eta
 *       Axis 3               Centrality
 *       Axis 4   Jet flavor (for MC) 1 = Quark, 2 = Gluon
 *
 */
void EECHistogramManager::LoadJetHistograms(){
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  
  // Define arrays to help find the histograms
  int axisIndices[2] = {0};
  int lowLimits[2] = {0};
  int highLimits[2] = {0};
  
  int nAxes = 1;           // Number of constraining axes for this iteration
  
  
  for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
    
    // Select the bin indices
    lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
    higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
    
    
    axisIndices[0] = 3; lowLimits[0] = lowerCentralityBin; highLimits[0] = higherCentralityBin;  // Centrality
    
    // If jet flavor is specified, only read jets of specific flavor
    if(fJetFlavor == 1 || fJetFlavor == 2){
      nAxes++;
      axisIndices[1] = 4; lowLimits[1] = fJetFlavor; highLimits[1] = fJetFlavor;  // Jet flavor
    }
    
    // Always load jet pT histograms
    fhJetPt[iCentralityBin] = FindHistogram(fInputFile,fJetHistogramName,0,nAxes,axisIndices,lowLimits,highLimits);
    
    if(!fLoadJets) continue;  // Only load the remaining jet histograms if selected
    
    fhJetPhi[iCentralityBin] = FindHistogram(fInputFile,fJetHistogramName,1,nAxes,axisIndices,lowLimits,highLimits);
    fhJetEta[iCentralityBin] = FindHistogram(fInputFile,fJetHistogramName,2,nAxes,axisIndices,lowLimits,highLimits);
    if(fLoad2DHistograms) fhJetEtaPhi[iCentralityBin] = FindHistogram2D(fInputFile,fJetHistogramName,1,2,nAxes,axisIndices,lowLimits,highLimits);
    
  } // Loop over centrality bins
}

/*
 * Loader for track histograms
 *
 * THnSparse for tracks:
 *
 *   Histogram name: track/trackUncorrected
 *
 *     Axis index       Content of axis
 * ----------------------------------------
 *       Axis 0            Track pT
 *       Axis 1            Track phi
 *       Axis 2            Track eta
 *       Axis 3            Centrality
 */
void EECHistogramManager::LoadTrackHistograms(){
  
  // Define arrays to help find the histograms
  int axisIndices[2] = {0};
  int lowLimits[2] = {0};
  int highLimits[2] = {0};
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  int duplicateRemoverTrackPt = -1;
  int lowerTrackPtBin = 0;
  int higherTrackPtBin = 0;
  
  // Loop over all track histograms
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    if(!fLoadTracks[iTrackType]) continue;  // Only load the selected track types
    
    for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
      
      // Select the bin indices
      lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
      higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
      
      // Setup axes with restrictions, (3 = centrality)
      axisIndices[0] = 3; lowLimits[0] = lowerCentralityBin; highLimits[0] = higherCentralityBin;
      
      fhTrackPt[iTrackType][iCentralityBin] = FindHistogram(fInputFile,fTrackHistogramNames[iTrackType],0,1,axisIndices,lowLimits,highLimits);
      fhTrackPhi[iTrackType][iCentralityBin][fnTrackPtBins] = FindHistogram(fInputFile,fTrackHistogramNames[iTrackType],1,1,axisIndices,lowLimits,highLimits);
      fhTrackEta[iTrackType][iCentralityBin][fnTrackPtBins] = FindHistogram(fInputFile,fTrackHistogramNames[iTrackType],2,1,axisIndices,lowLimits,highLimits);
      if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCentralityBin][fnTrackPtBins] = FindHistogram2D(fInputFile,fTrackHistogramNames[iTrackType],1,2,1,axisIndices,lowLimits,highLimits);
      
      for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
        
        // Select the bin indices for track pT
        lowerTrackPtBin = fTrackPtBinIndices[iTrackPtBin];
        higherTrackPtBin = fTrackPtBinIndices[iTrackPtBin+1]+duplicateRemoverTrackPt;
        
        // Add restriction for pT axis (0)
        axisIndices[1] = 0; lowLimits[1] = lowerTrackPtBin; highLimits[1] = higherTrackPtBin;
        
        // Read the angle histograms in track pT bins
        fhTrackPhi[iTrackType][iCentralityBin][iTrackPtBin] = FindHistogram(fInputFile,fTrackHistogramNames[iTrackType],1,2,axisIndices,lowLimits,highLimits);
        fhTrackEta[iTrackType][iCentralityBin][iTrackPtBin] = FindHistogram(fInputFile,fTrackHistogramNames[iTrackType],2,2,axisIndices,lowLimits,highLimits);
        if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCentralityBin][iTrackPtBin] = FindHistogram2D(fInputFile,fTrackHistogramNames[iTrackType],1,2,2,axisIndices,lowLimits,highLimits);
        
      } // Track pT loop
    } // Centrality loop
  } // Track category loop
}

/*
 * Loader for energy-energy correlator histograms
 *
 * THnSparse for energy-energy correlators:
 *
 *   Histogram name: energyEnergyCorrelator/energyEnergyCorrelatorJetPt/energyEnergyCorrelatorUncorrected/energyEnergyCorrelatorJetPtUncorrected
 *
 *     Axis index       Content of axis               Note
 * -----------------------------------------------------------------
 *       Axis 0              DeltaR
 *       Axis 1              Jet pT
 *       Axis 2           Track pT cut
 *       Axis 3            Centrality
 *       Axis 4           Pairing type       Same jet/reflected cone
 *       Axis 5             Subevent           Only relevant for MC
 */
void EECHistogramManager::LoadEnergyEnergyCorrelatorHistograms(){
  
  // Define arrays to help find the histograms
  int axisIndices[5] = {0};
  int lowLimits[5] = {0};
  int highLimits[5] = {0};
  
  // Define helper variables
  int duplicateRemover = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  int lowerJetPtBin = 0;
  int higherJetPtBin = 0;
  int lowerTrackPtBin = 0;
  int higherTrackPtBin = 0;
  THnSparseD *histogramArray;
  
  // Loop over all different energy-energy correlator histograms
  for(int iEnergyEnergyCorrelatorType = 0; iEnergyEnergyCorrelatorType < knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelatorType++){
    if(!fLoadEnergyEnergyCorrelatorHistograms[iEnergyEnergyCorrelatorType]) continue;  // Only load the selected energy-energy correlators
    
    // For track pT bins, we are looking at all the tracks above the lower threshold
    histogramArray = (THnSparseD*) fInputFile->Get(fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType]);
    higherTrackPtBin = histogramArray->GetAxis(2)->GetNbins()+1;
    
    // Loop over pairing types
    for(int iPairingType = 0; iPairingType < EECHistograms::knPairingTypes; iPairingType++){
      
      // If reflected cone histograms are not filled in the data file, do not try to load them
      if((iPairingType == EECHistograms::kSignalReflectedConePair || iPairingType == EECHistograms::kReflectedConePair) && !fCard->GetDoReflectedCone()) continue;
      
      // Setup axes with restrictions, (4 = pairing type)
      axisIndices[0] = 4; lowLimits[0] = iPairingType+1; highLimits[0] = iPairingType+1;
      
      // Loop over centrality bins
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
        
        // Select the centrality bin indices
        lowerCentralityBin = fCentralityBinIndices[iCentrality];
        higherCentralityBin = fCentralityBinIndices[iCentrality+1]+duplicateRemover;
        
        // Setup axes with restrictions, (3 = centrality)
        axisIndices[1] = 3; lowLimits[1] = lowerCentralityBin; highLimits[1] = higherCentralityBin;
        
        // Loop over track pT bins
        for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){
          
          // Select the track pT bin indices. Notice that we do not change the higher bin index
          lowerTrackPtBin = fTrackPtIndicesEEC[iTrackPt];
          
          // Add restriction for track pT axis (2 = track pT)
          axisIndices[2] = 2; lowLimits[2] = lowerTrackPtBin; highLimits[2] = higherTrackPtBin;
          
          // Read the energy-energy correlator histograms without jet pT restrictions
          fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][iPairingType][knSubeventTypes] = FindHistogram(fInputFile, fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], 0, 3, axisIndices, lowLimits, highLimits);
          
          // For PbPb MC, read the energy-energy correlator histograms without jet pT restrictions in subevent bins
          if(fSystemAndEnergy.Contains("PbPb MC")){
            for(int iSubevent = 0; iSubevent < knSubeventTypes; iSubevent++){
              
              // Add a restriction for the subevent axis (5 = subevent)
              axisIndices[3] = 5; lowLimits[3] = iSubevent+1; highLimits[3] = iSubevent+1;
              
              // Read the energy-energy correlator histograms
              fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][iPairingType][iSubevent] = FindHistogram(fInputFile, fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], 0, 4, axisIndices, lowLimits, highLimits);
              
            } // Subevent loop
          } // PbPb MC requirement
          
          // Loop over jet pT bins
          for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
            
            // Select the jet pT bin indices
            lowerJetPtBin = fJetPtIndicesEEC[iJetPt];
            higherJetPtBin = fJetPtIndicesEEC[iJetPt+1]+duplicateRemover;
            
            // Add restriction for jet pT axis (1 = jet pT)
            axisIndices[3] = 1; lowLimits[3] = lowerJetPtBin; highLimits[3] = higherJetPtBin;
            
            // Read the energy-energy correlator histograms
            fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][knSubeventTypes] = FindHistogram(fInputFile, fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], 0, 4, axisIndices, lowLimits, highLimits);
            
            // For PbPb MC, loop over subevent types
            if(fSystemAndEnergy.Contains("PbPb MC")){
              for(int iSubevent = 0; iSubevent < knSubeventTypes; iSubevent++){
                
                // Add a restriction for the subevent axis (5 = subevent)
                axisIndices[4] = 5; lowLimits[4] = iSubevent+1; highLimits[4] = iSubevent+1;
                
                // Read the energy-energy correlator histograms
                fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent] = FindHistogram(fInputFile, fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], 0, 5, axisIndices, lowLimits, highLimits);
                
              } // Subevent loop
            } // PbPb MC requirement
            
          } // Jet pT loop
        } // Track pT loop
      } // Centrality loop
    } // Pairing type loop (same jet/reflected cone)
  } // Energy-energy correlator type type loop
}

/*
 * Loader for jet pT closure histograms
 *
 * THnSparse for closure histograms:
 *
 *   Histogram name: jetPtClosure
 *
 *     Axis index                  Content of axis
 * -----------------------------------------------------------
 *       Axis 0              Matched generator level jet pT
 *       Axis 1               Matched reconstructed jet pT
 *       Axis 2                         Jet eta
 *       Axis 3                       Centrality
 *       Axis 4                      Quark / gluon
 *       Axis 5             Matched reco to gen jet pT ratio
 */
void EECHistogramManager::LoadJetPtClosureHistograms(){
  
  if(!fLoadJetPtClosureHistograms) return; // Do not load the histograms if they are not selected for loading
  
  // Define arrays to help find the histograms
  int axisIndices[4] = {0};
  int lowLimits[4] = {0};
  int highLimits[4] = {0};
  int nRestrictionAxes = 3;
  
  // Define helper variables
  int duplicateRemoverCentrality = -1;
  int lowerCentralityBin = 0;
  int higherCentralityBin = 0;
  
  // Load all the histograms from the file
  for(int iGenJetPt = 0; iGenJetPt < knGenJetPtBins; iGenJetPt++){
    for(int iClosureParticle = 0; iClosureParticle < EECHistograms::knClosureParticleTypes+1; iClosureParticle++){
      for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
        
        // Select the bin indices
        lowerCentralityBin = fCentralityBinIndices[iCentralityBin];
        higherCentralityBin = fCentralityBinIndices[iCentralityBin+1]+duplicateRemoverCentrality;
        
        // Setup the axes with restrictions
        axisIndices[0] = 1; lowLimits[0] = iGenJetPt+1;    highLimits[0] = iGenJetPt+1;             // Gen jet pT
        axisIndices[1] = 4; lowLimits[1] = lowerCentralityBin; highLimits[1] = higherCentralityBin; // Centrality
        axisIndices[2] = 5; lowLimits[2] = iClosureParticle+1; highLimits[2] = iClosureParticle+1;  // Quark/gluon
        
        // For the last closure particle bin no restrictions for quark/gluon jets
        if(iClosureParticle == EECHistograms::knClosureParticleTypes){
          
          // Remove the closure particle requirement from the restriction axes
          nRestrictionAxes--;
        }
        
        fhJetPtClosure[iGenJetPt][knJetEtaBins][iCentralityBin][iClosureParticle] = FindHistogram(fInputFile,"jetPtClosure",5,nRestrictionAxes,axisIndices,lowLimits,highLimits);
        
        
        // Eta binning for the closure histogram
        //            for(int iJetEta = 0; iJetEta < knJetEtaBins; iJetEta++){
        //
        //              // For the last closure particle bin no restrictions for quark/gluon jets
        //              /*if(iClosureParticle == EECHistograms::knClosureParticleTypes){
        //               nRestrictionAxes = 3;
        //               axisIndices[2] = 2; lowLimits[2] = iJetEta+1; highLimits[2] = iJetEta+1; // Jet eta
        //               } else {
        //               nRestrictionAxes = 4;
        //               axisIndices[3] = 2; lowLimits[3] = iJetEta+1; highLimits[3] = iJetEta+1; // Jet eta
        //               }
        //
        //               fhJetPtClosure[iClosureType][iGenJetPt][iJetEta][iCentralityBin][iClosureParticle] = FindHistogram(fInputFile,"jetPtClosure",5,nRestrictionAxes,axisIndices,lowLimits,highLimits);*/
        //
        //              // Fill the pT integrated eta slices only once
        //              if(iGenJetPt == 0){
        //
        //                // Setup the axes with restrictions
        //                nRestrictionAxes = 3;
        //                axisIndices[1] = 2; lowLimits[1] = iJetEta+1;    highLimits[1] = iJetEta+1;                 // Jet eta
        //                axisIndices[2] = 3; lowLimits[2] = lowerCentralityBin; highLimits[2] = higherCentralityBin; // Centrality
        //                axisIndices[3] = 4; lowLimits[3] = iClosureParticle+1; highLimits[3] = iClosureParticle+1;  // Quark/gluon
        //
        //                // For the last closure particle bin no restrictions for quark/gluon jets
        //                if(iClosureParticle == EECHistograms::knClosureParticleTypes){
        //
        //                  // Remove the last set array bin
        //                  nRestrictionAxes--;
        //                }
        //
        //                fhJetPtClosure[knGenJetPtBins][iJetEta][iCentralityBin][iClosureParticle] = FindHistogram(fInputFile,"jetPtClosure",5,nRestrictionAxes,axisIndices,lowLimits,highLimits);
        //              }
        //
        //            } // Jet eta bin loop
      } // Centrality loop
    } // Closure particle loop
  } // Gen jet pT loop
}

/*
 * Extract a 2D histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   TFile *inputFile = Inputfile containing the THnSparse to be read
 *   const char *name = Name of the THnSparse that is read
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH2D
 *   int yAxis = Index for the axis in THnSparse that is projected to y-axis for TH2D
 *   int nAxes = Number of axes that are restained for the projection
 *   int *axisNumber = Indices for the axes in THnSparse that are used as a restriction for the projection
 *   int *lowBinIndex = Indices of the lowest considered bins in the restriction axis
 *   int *highBinIndex = Indices of the highest considered bins in the restriction axis
 */
TH2D* EECHistogramManager::FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex){
  
  // Read the histogram with the given name from the file
  THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
  
  // If cannot find histogram, inform that it could not be found and return null
  if(histogramArray == nullptr){
    cout << "Could not find " << name << ". Skipping loading this histogram." << endl;
    return NULL;
  }
  
  // Apply the restrictions in the set of axes
  for(int i = 0; i < nAxes; i++) histogramArray->GetAxis(axisNumber[i])->SetRange(lowBinIndex[i],highBinIndex[i]);
  
  // Create a unique name for eeach histogram that is read from the file
  char newName[200];
  sprintf(newName,"%s",histogramArray->GetName());
  for(int iBinIndex = 0; iBinIndex < nAxes; iBinIndex++){
    sprintf(newName,"%s_%d=%d-%d",newName,axisNumber[iBinIndex],lowBinIndex[iBinIndex],highBinIndex[iBinIndex]);
  }
  
  // Project out the histogram and give it the created unique name
  TH2D *projectedHistogram = (TH2D*) histogramArray->Projection(yAxis,xAxis);
  projectedHistogram->SetName(newName);
  
  // Apply bin width normalization to the projected histogram
  projectedHistogram->Scale(1.0,"width");
  
  // Return the projected histogram
  return projectedHistogram;
}

/*
 * Extract a 2D histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   TFile *inputFile = Inputfile containing the THnSparse to be read
 *   const char *name = Name of the THnSparse that is read
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH2D
 *   int yAxis = Index for the axis in THnSparse that is projected to y-axis for TH2D
 *   int restrictionAxis = Index for the axis in THnSparse that is used as a restriction for the projection
 *   int lowBinIndex = Index of the lowest considered bin in the restriction axis
 *   int highBinIndex = Index of the highest considered bin in the restriction axis
 *   int restrictionAxis2 = Index for the axis in THnSparse that is used as a second restriction for the projection
 *   int lowBinIndex2 = Index of the lowest considered bin in the second restriction axis
 *   int highBinIndex2 = Index of the highest considered bin in the second restriction axis
 */
TH2D* EECHistogramManager::FindHistogram2D(TFile *inputFile, const char *name, int xAxis, int yAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2, int lowBinIndex2, int highBinIndex2){
  int restrictionAxes[2] = {restrictionAxis,restrictionAxis2};
  int lowBinIndices[2] = {lowBinIndex,lowBinIndex2};
  int highBinIndices[2] = {highBinIndex,highBinIndex2};
  int nAxes = 2;
  if(highBinIndex2 == 0 && lowBinIndex2 == 0) nAxes = 1;
  return FindHistogram2D(inputFile,name,xAxis,yAxis,nAxes,restrictionAxes,lowBinIndices,highBinIndices);
}

/*
 * Extract a histogram with given restrictions on other axes in THnSparse
 *
 *  Arguments:
 *   TFile *inputFile = Inputfile containing the THnSparse to be read
 *   const char *name = Name of the THnSparse that is read
 *   int xAxis = Index for the axis in THnSparse that is projected to x-axis for TH1D
 *   int nAxes = Number of axes that are restained for the projection
 *   int *axisNumber = Indices for the axes in THnSparse that are used as a restriction for the projection
 *   int *lowBinIndex = Indices of the lowest considered bins in the restriction axis
 *   int *highBinIndex = Indices of the highest considered bins in the restriction axis
 */
TH1D* EECHistogramManager::FindHistogram(TFile *inputFile, const char *name, int xAxis, int nAxes, int *axisNumber, int *lowBinIndex, int *highBinIndex){
  
  // Read the histogram with the given name from the file
  THnSparseD *histogramArray = (THnSparseD*) inputFile->Get(name);
  
  // If cannot find histogram, inform that it could not be found and return null
  if(histogramArray == nullptr){
    cout << "Could not find " << name << ". Skipping loading this histogram." << endl;
    return NULL;
  }
  
  // Apply the restrictions in the set of axes
  for(int i = 0; i < nAxes; i++) histogramArray->GetAxis(axisNumber[i])->SetRange(lowBinIndex[i],highBinIndex[i]);
  
  // Create a unique name for each histogram that is read from the file
  char newName[200];
  sprintf(newName,"%s",histogramArray->GetName());
  for(int iBinIndex = 0; iBinIndex < nAxes; iBinIndex++){
    sprintf(newName,"%s_%d=%d-%d",newName,axisNumber[iBinIndex],lowBinIndex[iBinIndex],highBinIndex[iBinIndex]);
  }

  // Project out the histogram and give it the created unique name
  TH1D *projectedHistogram = NULL;
  
  // Check that we are not trying to project a non-existing axis
  if(xAxis < histogramArray->GetNdimensions()){
    projectedHistogram = (TH1D*) histogramArray->Projection(xAxis);
    projectedHistogram->SetName(newName);
  
    // Apply bin width normalization to the projected histogram
    projectedHistogram->Scale(1.0,"width");
  }
  
  // Return the projected histogram
  return projectedHistogram;
}

/*
 * Extract a histogram from a given centrality bin from THnSparseD
 *
 *  Arguments:
 *   TFile *inputFile = Inputfile containing the THnSparse to be read
 *   const char *name = Name of the THnSparse that is read
 *   int xAxis = Index for the axis in THnSparse that is projected to TH1D
 *   int restrictionAxis = Index for the axis in THnSparse that is used as a restriction for the projection
 *   int lowBinIndex = Index of the lowest considered bin in the restriction axis
 *   int highBinIndex = Index of the highest considered bin in the restriction axis
 *   int restrictionAxis2 = Index for the axis in THnSparse that is used as a second restriction for the projection
 *   int lowBinIndex2 = Index of the lowest considered bin in the second restriction axis
 *   int highBinIndex2 = Index of the highest considered bin in the second restriction axis
 */
TH1D* EECHistogramManager::FindHistogram(TFile *inputFile, const char *name, int xAxis, int restrictionAxis, int lowBinIndex, int highBinIndex, int restrictionAxis2, int lowBinIndex2, int highBinIndex2){
  int restrictionAxes[2] = {restrictionAxis,restrictionAxis2};
  int lowBinIndices[2] = {lowBinIndex,lowBinIndex2};
  int highBinIndices[2] = {highBinIndex,highBinIndex2};
  int nAxes = 2;
  if(highBinIndex2 == 0 && lowBinIndex2 == 0) nAxes = 1;
  return FindHistogram(inputFile,name,xAxis,nAxes,restrictionAxes,lowBinIndices,highBinIndices);
}

/*
 * Write all the loaded histograms into a file
 *
 *  const char* fileName = Name of the file to which the histograms are written
 *  const char* fileOption = Option given to the file when it is loaded
 */
void EECHistogramManager::Write(const char* fileName, const char* fileOption){
  
  // Create the output file
  TFile *outputFile = new TFile(fileName,fileOption);
  
  // Helper variable for renaming the saved histograms
  char histogramNamer[200];
  
  // Write the event information histograms to the output file
  if(fLoadEventInformation){
    fhEvents->Write("",TObject::kOverwrite);             // Number of events surviving different event cuts
    fhVertexZ->Write("",TObject::kOverwrite);            // Vertex z position
    fhVertexZWeighted->Write("",TObject::kOverwrite);    // MC weighted vertex z position
    fhTrackCuts->Write("",TObject::kOverwrite);          // Number of tracks surviving different track cuts
    fhCentrality->Write("",TObject::kOverwrite);         // Centrality in all events
    fhCentralityWeighted->Write("",TObject::kOverwrite); // MC weighted centrality in all events
    fhPtHat->Write("",TObject::kOverwrite);              // pT hat for MC events (only meaningful for MC)
    fhPtHatWeighted->Write("",TObject::kOverwrite);      // Weighted pT hat distribution (only meaningful for MC)
    
    // Create a directory for the multiplicity histograms if it does not already exist
    if(!gDirectory->GetDirectory("multiplicity")) gDirectory->mkdir("multiplicity");
    gDirectory->cd("multiplicity");
    
    // Loop over centrality
    for(int iCentrality = 0; iCentrality < fnCentralityBins; iCentrality++){
      sprintf(histogramNamer, "multiplicity_C%d", iCentrality);
      if(fhMultiplicity[iCentrality]) fhMultiplicity[iCentrality]->Write(histogramNamer, TObject::kOverwrite);
      
      sprintf(histogramNamer, "multiplicityWeighted_C%d", iCentrality);
      if(fhMultiplicityWeighted[iCentrality]) fhMultiplicityWeighted[iCentrality]->Write(histogramNamer, TObject::kOverwrite);
    }
    
    // Histograms without cerntality binning
    sprintf(histogramNamer, "multiplicity");
    if(fhMultiplicity[fnCentralityBins]) fhMultiplicity[fnCentralityBins]->Write(histogramNamer, TObject::kOverwrite);
    
    sprintf(histogramNamer, "multiplicityWeighted");
    if(fhMultiplicityWeighted[fnCentralityBins]) fhMultiplicityWeighted[fnCentralityBins]->Write(histogramNamer, TObject::kOverwrite);
    
    // Multiplicity vs. centrality maps
    sprintf(histogramNamer, "multiplicityMap");
    if(fhMultiplicityMap) fhMultiplicityMap->Write(histogramNamer, TObject::kOverwrite);
    
    sprintf(histogramNamer, "multiplicityMapWeighted");
    if(fhMultiplicityMapWeighted) fhMultiplicityMapWeighted->Write(histogramNamer, TObject::kOverwrite);
    
    // Return back to main directory
    gDirectory->cd("../");
  }
 
  // Write the jet histograms to the output file
  WriteJetHistograms();
  
  // Write the track histograms to the output file
  WriteTrackHistograms();
  
  // Write the energy-energy correlator histograms to the output file
  WriteEnergyEnergyCorrelatorHistograms();
  
  // Write the jet pT closure histograms to a file
  WriteClosureHistograms();
  
  // Write the card to the output file if it is not already written
  if(!gDirectory->GetDirectory("JCard")) fCard->Write(outputFile);
  
  // Close the file after everything is written
  outputFile->Close();
  
  // Delete the outputFile object
  delete outputFile;
  
}

/*
 * Write the jet histograms to the file that is currently open
 */
void EECHistogramManager::WriteJetHistograms(){
  
  // Helper variable for histogram naming
  char histogramNamer[200];
  
  // Write the jet histograms to the output file
  if(!fLoadJets) return;  // Only write the jet histograms if they are loaded
  
  // Create a directory for the histograms if it does not already exist
  if(!gDirectory->GetDirectory(fJetHistogramName)) gDirectory->mkdir(fJetHistogramName);
  gDirectory->cd(fJetHistogramName);
  
  for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
    
    // Check that the histograms are actually there before trying to save them.
    if(fhJetPt[iCentralityBin] == NULL) {
      cout << "Could not find histograms of type " << fJetHistogramName << " to write. Will skip writing these." << endl;
      continue;
    }
    
    // Jet pT
    sprintf(histogramNamer,"%sPt_C%d",fJetHistogramName,iCentralityBin);
    if(fhJetPt[iCentralityBin]) fhJetPt[iCentralityBin]->Write(histogramNamer, TObject::kOverwrite);
    
    // Jet phi
    sprintf(histogramNamer,"%sPhi_C%d",fJetHistogramName,iCentralityBin);
    if(fhJetPhi[iCentralityBin]) fhJetPhi[iCentralityBin]->Write(histogramNamer, TObject::kOverwrite);
    
    // Jet eta
    sprintf(histogramNamer,"%sEta_C%d",fJetHistogramName,iCentralityBin);
    if(fhJetEta[iCentralityBin]) fhJetEta[iCentralityBin]->Write(histogramNamer, TObject::kOverwrite);
    
    // Jet eta-phi
    sprintf(histogramNamer,"%sEtaPhi_C%d",fJetHistogramName,iCentralityBin);
    if(fLoad2DHistograms && fhJetEtaPhi[iCentralityBin]) fhJetEtaPhi[iCentralityBin]->Write(histogramNamer, TObject::kOverwrite);
    
  } // Loop over centrality bins
  
  // Return back to main directory
  gDirectory->cd("../");
  
}

/*
 * Write the track histograms to the file that is currently open
 */
void EECHistogramManager::WriteTrackHistograms(){
  
  // Helper variable for histogram naming
  char histogramNamer[200];
  
  // Write the track histograms to the output file
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    if(!fLoadTracks[iTrackType]) continue;  // Only write the loaded track types
    
    // Create a directory for the histograms if it does not already exist
    if(!gDirectory->GetDirectory(fTrackHistogramNames[iTrackType])) gDirectory->mkdir(fTrackHistogramNames[iTrackType]);
    gDirectory->cd(fTrackHistogramNames[iTrackType]);
    
    for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
      
      // Track pT
      sprintf(histogramNamer, "%sPt_C%d", fTrackHistogramNames[iTrackType], iCentralityBin);
      fhTrackPt[iTrackType][iCentralityBin]->Write(histogramNamer, TObject::kOverwrite);
      
      // pT integrated track phi
      sprintf(histogramNamer, "%sPhi_C%dT%d", fTrackHistogramNames[iTrackType], iCentralityBin, fnTrackPtBins);
      fhTrackPhi[iTrackType][iCentralityBin][fnTrackPtBins]->Write(histogramNamer, TObject::kOverwrite);
      
      // pT integrated track eta
      sprintf(histogramNamer, "%sEta_C%dT%d", fTrackHistogramNames[iTrackType], iCentralityBin, fnTrackPtBins);
      fhTrackEta[iTrackType][iCentralityBin][fnTrackPtBins]->Write(histogramNamer, TObject::kOverwrite);
      
      // pT integrated track eta-phi
      sprintf(histogramNamer, "%sEtaPhi_C%dT%d", fTrackHistogramNames[iTrackType], iCentralityBin, fnTrackPtBins);
      if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCentralityBin][fnTrackPtBins]->Write(histogramNamer, TObject::kOverwrite);
      
      for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
        
        // Track phi in track pT bins
        sprintf(histogramNamer, "%sPhi_C%dT%d", fTrackHistogramNames[iTrackType], iCentralityBin, iTrackPtBin);
        fhTrackPhi[iTrackType][iCentralityBin][iTrackPtBin]->Write(histogramNamer, TObject::kOverwrite);
        
        // Track eta in track pT bins
        sprintf(histogramNamer, "%sEta_C%dT%d", fTrackHistogramNames[iTrackType], iCentralityBin, iTrackPtBin);
        fhTrackEta[iTrackType][iCentralityBin][iTrackPtBin]->Write(histogramNamer, TObject::kOverwrite);
        
        // Track eta-phi in track pT bins
        sprintf(histogramNamer, "%sEtaPhi_C%dT%d", fTrackHistogramNames[iTrackType], iCentralityBin, iTrackPtBin);
        if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCentralityBin][iTrackPtBin]->Write(histogramNamer, TObject::kOverwrite);
        
      } // Track pT loop
    } // Centrality loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Track category loop

}

/*
 * Write the energy-energy correlator histograms to the file that is currently open
 */
void EECHistogramManager::WriteEnergyEnergyCorrelatorHistograms(){
  
  // Helper variable for histogram naming
  char histogramNamer[200];
  
  // Loop over energy-energy correaltor histogram types
  for(int iEnergyEnergyCorrelatorType = 0; iEnergyEnergyCorrelatorType < knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelatorType++){
    
    // Only write the histograms that have been loaded
    if(!fLoadEnergyEnergyCorrelatorHistograms[iEnergyEnergyCorrelatorType]) continue;
    
    // Create a directory for the histograms if it does not already exist
    if(!gDirectory->GetDirectory(fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType])) gDirectory->mkdir(fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType]);
    gDirectory->cd(fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType]);
    
    // Loop over pairing types (same jet/reflected cone)
    for(int iPairingType = 0; iPairingType < EECHistograms::knPairingTypes; iPairingType++){
      
      // Loop over centrality
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
        
        // Loop over track pT
        for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){
          
          // Write histograms without jet pT binning
          sprintf(histogramNamer,"%s%s_C%dT%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fPairingTypeSaveName[iPairingType], iCentrality, iTrackPt);
          if(fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][iPairingType][knSubeventTypes]) fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][iPairingType][knSubeventTypes]->Write(histogramNamer, TObject::kOverwrite);
          
          // For PbPb MC, write histograms without jet pT and with subevent type binning
          if(fSystemAndEnergy.Contains("PbPb MC")){
            for(int iSubevent = 0; iSubevent < knSubeventTypes; iSubevent++){
              
              // Write the energy-energy correlator histograms with subevent binning
              sprintf(histogramNamer,"%s%s_C%dT%dS%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fPairingTypeSaveName[iPairingType], iCentrality, iTrackPt, iSubevent);
              if(fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][iPairingType][iSubevent]) fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][iPairingType][iSubevent]->Write(histogramNamer, TObject::kOverwrite);
              
            } // Subevent type loop
          } // Data is PbPb MC
          
          // Loop over jet pT
          for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
            
            // Write the energy-energy correlator histograms
            sprintf(histogramNamer,"%s%s_C%dT%dJ%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fPairingTypeSaveName[iPairingType], iCentrality, iTrackPt, iJetPt);
            if(fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][knSubeventTypes]) fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][knSubeventTypes]->Write(histogramNamer, TObject::kOverwrite);
            
            // For PbPb MC, loop over subevent types
            if(fSystemAndEnergy.Contains("PbPb MC")){
              for(int iSubevent = 0; iSubevent < knSubeventTypes; iSubevent++){
                
                // Write the energy-energy correlator histograms with subevent binning
                sprintf(histogramNamer,"%s%s_C%dT%dJ%dS%d",fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fPairingTypeSaveName[iPairingType], iCentrality, iTrackPt, iJetPt, iSubevent);
                if(fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent]) fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent]->Write(histogramNamer, TObject::kOverwrite);
                
              } // Subevent type loop
            } // Data is PbPb MC
            
          } // Loop over jet pT bins
        } // Loop over track pT bins
      } // Loop over centrality bins
    } // Pairing type loop (same jet/reflected cone)
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Loop over different energy-energy correlator types
  
}

/*
 * Write the closure histograms to the file that is currently open
 */
void EECHistogramManager::WriteClosureHistograms(){
  
  // Helper variable for histogram naming
  char histogramNamer[200];
  
  if(fLoadJetPtClosureHistograms){
    
    
    // Create a directory for the histograms if it does not already exist
    sprintf(histogramNamer,"jetPtClosure_%s",fJetHistogramName);
    if(!gDirectory->GetDirectory(histogramNamer)) gDirectory->mkdir(histogramNamer);
    gDirectory->cd(histogramNamer);
    
    // Centrality loop
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      
      // Loop over closure particles (quark/gluon/no selection)
      for(int iClosureParticle = 0; iClosureParticle < EECHistograms::knClosureParticleTypes+1; iClosureParticle++){
        
        // Loop over generator level jet pT bins
        for(int iGenJetPt = 0; iGenJetPt <= knGenJetPtBins; iGenJetPt++){
          
          // Loop ovet jet eta bins
          for(int iJetEta = 0; iJetEta <= knJetEtaBins; iJetEta++){
            
            // Only write histogram that are non-NULL
            if(fhJetPtClosure[iGenJetPt][iJetEta][iCentrality][iClosureParticle]){
              sprintf(histogramNamer, "jetPtClosure_%s%s_C%dT%dE%d", fJetHistogramName, fClosureParticleName[iClosureParticle], iCentrality, iGenJetPt, iJetEta);
              fhJetPtClosure[iGenJetPt][iJetEta][iCentrality][iClosureParticle]->Write(histogramNamer, TObject::kOverwrite);
            }
            
          } // Jet eta bin loop
        } // Generator level jet pT loop
      } // Closure particle type (quark/gluon) loop
    } // Centrality loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Writing jet pT closure histograms
  
}

/*
 * Load the selected histograms from a file containing readily processed histograms
 */
void EECHistogramManager::LoadProcessedHistograms(){
  
  // Helper variable for finding names of loaded histograms
  char histogramNamer[300];
  
  // Always load the number of events histogram
  fhEvents = (TH1D*) fInputFile->Get("nEvents");                           // Number of events surviving different event cuts
  
  // Load the event information histograms
  if(fLoadEventInformation){
    fhVertexZ = (TH1D*) fInputFile->Get("vertexZ");                        // Vertex z position
    fhVertexZWeighted = (TH1D*) fInputFile->Get("vertexZweighted");        // MC weighted vertex z position
    fhTrackCuts = (TH1D*) fInputFile->Get("trackCuts");                    // Number of tracks surviving different track cuts
    fhCentrality = (TH1D*) fInputFile->Get("centrality");                  // Centrality in all events
    fhCentralityWeighted = (TH1D*) fInputFile->Get("centralityWeighted");  // MC weighted centrality in all events
    fhPtHat = (TH1D*) fInputFile->Get("pthat");                            // pT hat for MC events
    fhPtHatWeighted = (TH1D*) fInputFile->Get("pthatWeighted");            // Weighted pT hat for MC events
    
    for(int iCentrality = 0; iCentrality < fnCentralityBins; iCentrality++){
      
      sprintf(histogramNamer, "multiplicity/multiplicity_C%d", iCentrality);
      fhMultiplicity[iCentrality] = (TH1D*) fInputFile->Get(histogramNamer);
      
      sprintf(histogramNamer, "multiplicity/multiplicityWeighted_C%d", iCentrality);
      fhMultiplicityWeighted[iCentrality] = (TH1D*) fInputFile->Get(histogramNamer);
      
    }
    
    // Multiplicity histograms without centrality bins
    sprintf(histogramNamer, "multiplicity/multiplicity");
    fhMultiplicity[fnCentralityBins] = (TH1D*) fInputFile->Get(histogramNamer);
    
    sprintf(histogramNamer, "multiplicity/multiplicityWeighted");
    fhMultiplicityWeighted[fnCentralityBins] = (TH1D*) fInputFile->Get(histogramNamer);
    
    // Multiplicity vs. centrality maps
    sprintf(histogramNamer, "multiplicity/multiplicityMap");
    fhMultiplicityMap = (TH2D*) fInputFile->Get(histogramNamer);
    
    sprintf(histogramNamer, "multiplicity/multiplicityMapWeighted");
    fhMultiplicityMapWeighted = (TH2D*) fInputFile->Get(histogramNamer);
    
  }
  
  // Load the jet histograms from the input file
  
  for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
    
    // Always load jet pT histograms
    sprintf(histogramNamer,"%s/%sPt_C%d", fJetHistogramName, fJetHistogramName, iCentralityBin);
    fhJetPt[iCentralityBin] = (TH1D*) fInputFile->Get(histogramNamer);
    
    if(!fLoadJets) continue;  // Only load the loaded the selected histograms
    
    // Jet phi
    sprintf(histogramNamer,"%s/%sPhi_C%d", fJetHistogramName, fJetHistogramName, iCentralityBin);
    fhJetPhi[iCentralityBin] = (TH1D*) fInputFile->Get(histogramNamer);
    
    // Jet eta
    sprintf(histogramNamer,"%s/%sEta_C%d", fJetHistogramName, fJetHistogramName, iCentralityBin);
    fhJetEta[iCentralityBin] = (TH1D*) fInputFile->Get(histogramNamer);
    
    // Jet eta-phi
    sprintf(histogramNamer,"%s/%sEtaPhi_C%d", fJetHistogramName, fJetHistogramName, iCentralityBin);
    if(fLoad2DHistograms) fhJetEtaPhi[iCentralityBin] = (TH2D*) fInputFile->Get(histogramNamer);
  } // Loop over centrality bins
  
  
  // Load the track histograms from the input file
  for(int iTrackType = 0; iTrackType < knTrackCategories; iTrackType++){
    if(!fLoadTracks[iTrackType]) continue;  // Only load the selected track types
      
      for(int iCentralityBin = fFirstLoadedCentralityBin; iCentralityBin <= fLastLoadedCentralityBin; iCentralityBin++){
        
        // Track pT
        sprintf(histogramNamer,"%s/%sPt_C%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],iCentralityBin);
        fhTrackPt[iTrackType][iCentralityBin] = (TH1D*) fInputFile->Get(histogramNamer);
        
        // pT integrated track phi
        sprintf(histogramNamer,"%s/%sPhi_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],iCentralityBin,fnTrackPtBins);
        fhTrackPhi[iTrackType][iCentralityBin][fnTrackPtBins] = (TH1D*) fInputFile->Get(histogramNamer);
        
        // pT integrated track eta
        sprintf(histogramNamer,"%s/%sEta_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],iCentralityBin,fnTrackPtBins);
        fhTrackEta[iTrackType][iCentralityBin][fnTrackPtBins] = (TH1D*) fInputFile->Get(histogramNamer);
        
        // pT integrated track eta-phi
        sprintf(histogramNamer,"%s/%sEtaPhi_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],iCentralityBin,fnTrackPtBins);
        if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCentralityBin][fnTrackPtBins] = (TH2D*) fInputFile->Get(histogramNamer);
        
        for(int iTrackPtBin = fFirstLoadedTrackPtBin; iTrackPtBin <= fLastLoadedTrackPtBin; iTrackPtBin++){
          
          // Track phi in track pT bins
          sprintf(histogramNamer,"%s/%sPhi_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],iCentralityBin,iTrackPtBin);
          fhTrackPhi[iTrackType][iCentralityBin][iTrackPtBin] = (TH1D*) fInputFile->Get(histogramNamer);
          
          // Track eta in track pT bins
          sprintf(histogramNamer,"%s/%sEta_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],iCentralityBin,iTrackPtBin);
          fhTrackEta[iTrackType][iCentralityBin][iTrackPtBin] = (TH1D*) fInputFile->Get(histogramNamer);
          
          // Track eta-phi in track pT bins
          sprintf(histogramNamer,"%s/%sEtaPhi_C%dT%d",fTrackHistogramNames[iTrackType],fTrackHistogramNames[iTrackType],iCentralityBin,iTrackPtBin);
          if(fLoad2DHistograms) fhTrackEtaPhi[iTrackType][iCentralityBin][iTrackPtBin] = (TH2D*) fInputFile->Get(histogramNamer);
          
        } // Track pT loop
      } // Centrality loop
  } // Track category loop
  
  // Load the energy-energy correlator histograms from the input file
  for(int iEnergyEnergyCorrelatorType = 0; iEnergyEnergyCorrelatorType < knEnergyEnergyCorrelatorTypes; iEnergyEnergyCorrelatorType++){
    
    // Only load the selected types of histograms
    if(!fLoadEnergyEnergyCorrelatorHistograms[iEnergyEnergyCorrelatorType]) continue;
    
    for(int iPairingType = 0; iPairingType < EECHistograms::knPairingTypes; iPairingType++){
      for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
        for(int iTrackPt = fFirstLoadedTrackPtBinEEC; iTrackPt <= fLastLoadedTrackPtBinEEC; iTrackPt++){
          
          // Load the histograms without jet pT binning
          sprintf(histogramNamer,"%s/%s%s_C%dT%d", fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fPairingTypeSaveName[iPairingType], iCentrality, iTrackPt);
          fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][iPairingType][knSubeventTypes] = (TH1D*) fInputFile->Get(histogramNamer);
          
          // For PbPb MC, load histograms without jet pT and with subevent type binning
          if(fSystemAndEnergy.Contains("PbPb MC")){
            for(int iSubevent = 0; iSubevent < knSubeventTypes; iSubevent++){
              
              // Write the energy-energy correlator histograms with subevent binning
              sprintf(histogramNamer,"%s/%s%s_C%dT%dS%d", fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fPairingTypeSaveName[iPairingType], iCentrality, iTrackPt, iSubevent);
              fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][fnJetPtBinsEEC][iTrackPt][iPairingType][iSubevent] = (TH1D*) fInputFile->Get(histogramNamer);
              
            } // Subevent type loop
          } // Data is PbPb MC
          
          // Loop over jet pT
          for(int iJetPt = fFirstLoadedJetPtBinEEC; iJetPt <= fLastLoadedJetPtBinEEC; iJetPt++){
            
            // Make sure that jet pT integrated histogram is not overwritten by null
            if(fLastLoadedJetPtBinEEC >= fnJetPtBinsEEC) continue;
            
            // Load the energy-energy correlator histograms
            sprintf(histogramNamer,"%s/%s%s_C%dT%dJ%d", fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fPairingTypeSaveName[iPairingType], iCentrality, iTrackPt, iJetPt);
            fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][knSubeventTypes] = (TH1D*) fInputFile->Get(histogramNamer);
            
            // For PbPb MC, loop over subevent types
            if(fSystemAndEnergy.Contains("PbPb MC")){
              for(int iSubevent = 0; iSubevent < knSubeventTypes; iSubevent++){
                
                // Write the energy-energy correlator histograms with subevent binning
                sprintf(histogramNamer,"%s/%s%s_C%dT%dJ%dS%d", fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType], fPairingTypeSaveName[iPairingType], iCentrality, iTrackPt, iJetPt, iSubevent);
                fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent] = (TH1D*) fInputFile->Get(histogramNamer);
                
              } // Subevent type loop
            } // Data is PbPb MC
            
          } // Jet pT loop
        } // Track pT loop
      } // Centrality loop
    } // Pairing type loop
  } // Load energy-energy correlator histograms
  
  // Load the jet pT closure histograms from a processed file
  if(fLoadJetPtClosureHistograms){
    
    // Centrality loop
    for(int iCentrality = fFirstLoadedCentralityBin; iCentrality <= fLastLoadedCentralityBin; iCentrality++){
      
      // Loop over closure particles (quark/gluon/no selection)
      for(int iClosureParticle = 0; iClosureParticle < EECHistograms::knClosureParticleTypes+1; iClosureParticle++){
        
        // Loop over generator level jet pT bins
        for(int iGenJetPt = 0; iGenJetPt <= knGenJetPtBins; iGenJetPt++){
          
          // Loop over jet eta bins
          for(int iJetEta = 0; iJetEta <= knJetEtaBins; iJetEta++){
            
            sprintf(histogramNamer, "jetPtClosure_%s/jetPtClosure_%s%s_C%dT%dE%d", fJetHistogramName, fJetHistogramName, fClosureParticleName[iClosureParticle], iCentrality, iGenJetPt, iJetEta);
            fhJetPtClosure[iGenJetPt][iJetEta][iCentrality][iClosureParticle] = (TH1D*) fInputFile->Get(histogramNamer);
            
          } // Jet eta bin loop
        } // Generator level jet pT loop
      } // Closure particle type (quark/gluon) loop
    } // Centrality loop
    
    // Return back to main directory
    gDirectory->cd("../");
    
  } // Opening jet pT closure histograms
  
}

/*
 * Read the bin indices for given bin borders
 *
 *  Arguments:
 *   const char* histogramName = Name of the histogram from which the bin indices are searched
 *   const int nBins = Number of bins for the indices
 *   int *binIndices = Array of integers to be filled with bin index information read from the file
 *   const double *binBorders = Array for bin borders that are searched from the file
 *   const int iAxis = Index of the axis used for reading bin indices
 */
void EECHistogramManager::SetBinIndices(const char* histogramName, const int nBins, int *binIndices, const double *binBorders, const int iAxis){
  TH1D* hBinner = FindHistogram(fInputFile,histogramName,iAxis,0,0,0);
  for(int iBin = 0; iBin < nBins+1; iBin++){
    binIndices[iBin] = hBinner->GetXaxis()->FindBin(binBorders[iBin]);
  }
}

/*
 * Read the bin indices for given bin borders
 *
 *  Arguments:
 *   const char* histogramName = Name of the histogram from which the bin indices are searched
 *   const int nBins = Number of bins for the indices
 *   double *copyBinBorders = Array to which a copy of bin borders is made
 *   int *binIndices = Array of integers to be filled with bin index information read from the file
 *   const double *binBorders = Array for bin borders that are searched from the file
 *   const int iAxis = Index of the axis used for reading bin indices
 *   const bool setIndices = true: Set both bin indices and bin borders, false: Set only bin borders
 */
void EECHistogramManager::SetBinBordersAndIndices(const char* histogramName, const int nBins, double *copyBinBorders, int *binIndices, const double *binBorders, const int iAxis, const bool setIndices){
  TH1D* hBinner;
  if(setIndices) hBinner = FindHistogram(fInputFile,histogramName,iAxis,0,0,0);
  for(int iBin = 0; iBin < nBins+1; iBin++){
    copyBinBorders[iBin] = binBorders[iBin];
    if(setIndices) binIndices[iBin] = hBinner->GetXaxis()->FindBin(binBorders[iBin]);
  }
}

/*
 * Set up generic bin borders and indices according to provided bin borders
 *
 *  const bool readBinsFromFile = True: Disregard given bins ans use the ones in fCard. False: Use given bins
 *  const char* histogramName = Name of the histogram from which the bin indices are searched
 *  const int iAxis = Axis from which the set indices can be found
 *  int nSetBins = Number of bins that is set
 *  double* setBinBorders = Bin borders that are set
 *  int* setBinIndices = Bin indices that are set
 *  const int nBins = New number of bins that is given
 *  const double *binBorders = New bin borders that are given
 *  const char* errorMessage = Type of the set bins to be printed in possible error message
 *  const int maxBins = Maximum number of allowed bins of this type
 *  const bool setIndices = Set the bin indices in THnSparse
 */
void EECHistogramManager::SetGenericBins(const bool readBinsFromFile, const char* histogramName, const int iAxis, int nSetBins, double* setBinBorders, int* setBinIndices, const int nBins, const double *binBorders, const char* errorMessage, const int maxBins, const bool setIndices){
  
  // If bins are read from file, do not use the given bin borders
  if(readBinsFromFile){
    if(setIndices) SetBinIndices(histogramName, nSetBins, setBinIndices, setBinBorders, iAxis);
  } else { // If the given bin borders are use, update the number of bins and bin borders according to input
    if(nBins <= maxBins){
      nSetBins = nBins;
      SetBinBordersAndIndices(histogramName, nSetBins, setBinBorders, setBinIndices, binBorders, iAxis, setIndices);
    } else {
      cout << "Error! Too many " << errorMessage << " bins given. Maximum number is " << maxBins << ". Will not set bins." << endl;
    }
  }
}

/*
 * Set up centrality bin borders and indices according to provided bin borders
 *
 *  const bool readBinsFromFile = True: Disregard given bins ans use the ones in fCard. False: Use given bins
 *  const int nBins = Number of given centrality bins
 *  const double *binBorders = New bin borders for centrality
 *  const bool setIndices = Set the bin indices in THnSparse
 */
void EECHistogramManager::SetCentralityBins(const bool readBinsFromFile, const int nBins, const double *binBorders, const bool setIndices){
  
  SetGenericBins(readBinsFromFile, fEnergyEnergyCorrelatorHistogramNames[kEnergyEnergyCorrelator], 3, fnCentralityBins, fCentralityBinBorders, fCentralityBinIndices, nBins, binBorders, "centrality", kMaxCentralityBins, setIndices);
  
}

/*
 * Set up track pT bin borders and indices according to provided bin borders
 *
 *  const bool readBinsFromFile = True: Disregard given bins ans use the ones in fCard. False: Use given bins
 *  const int nBins = Number of given track pT bins
 *  const double *binBorders = New bin borders for track pT
 *  const bool setIndices = Set the bin indices in THnSparse
 */
void EECHistogramManager::SetTrackPtBins(const bool readBinsFromFile, const int nBins, const double *binBorders, const bool setIndices){
  
  SetGenericBins(readBinsFromFile, "track", 0, fnTrackPtBins, fTrackPtBinBorders, fTrackPtBinIndices, nBins, binBorders, "track pT", kMaxTrackPtBins, setIndices);
  
}

/*
 * Set up jet pT bin indices for energy-energy correlator according to provided bin borders
 *
 *  const bool readBinsFromFile = True: Disregard given bins ans use the ones in fCard. False: Use given bins
 *  const int nBins = Number of given track pT bins
 *  const double *binBorders = New bin borders for track pT
 *  const bool setIndices = Set the bin indices in THnSparse
 */
void EECHistogramManager::SetJetPtBinsEEC(const bool readBinsFromFile, const int nBins, const double *binBorders, bool setIndices){
  
  SetGenericBins(readBinsFromFile, fEnergyEnergyCorrelatorHistogramNames[kEnergyEnergyCorrelator], 1, fnJetPtBinsEEC, fJetPtBinBordersEEC, fJetPtIndicesEEC, nBins, binBorders, "jet pT EEC", kMaxJetPtBinsEEC, setIndices);
  
}

/*
 * Set up track pT bin indices for energy-energy correlator according to provided bin borders
 *
 *  const bool readBinsFromFile = True: Disregard given bins ans use the ones in fCard. False: Use given bins
 *  const int nBins = Number of given track pT bins
 *  const double *binBorders = New bin borders for track pT
 *  const bool setIndices = Set the bin indices in THnSparse
 */
void EECHistogramManager::SetTrackPtBinsEEC(const bool readBinsFromFile, const int nBins, const double *binBorders, bool setIndices ){
 
  SetGenericBins(readBinsFromFile, fEnergyEnergyCorrelatorHistogramNames[kEnergyEnergyCorrelator], 2, fnTrackPtBinsEEC, fTrackPtBinBordersEEC, fTrackPtIndicesEEC, nBins, binBorders, "track pT EEC", kMaxTrackPtBinsEEC, setIndices);
  
}

// Setter for loading event information
void EECHistogramManager::SetLoadEventInformation(const bool loadOrNot){
  fLoadEventInformation = loadOrNot;
}

// Setter for loading jet histograms
void EECHistogramManager::SetLoadJetHistograms(const bool loadOrNot){
  fLoadJets = loadOrNot;
}

// Setter for loading tracks
void EECHistogramManager::SetLoadTracks(const bool loadOrNot){
  fLoadTracks[kTrack] = loadOrNot;
}

// Setter for loading uncorrected tracks
void EECHistogramManager::SetLoadTracksUncorrected(const bool loadOrNot){
  fLoadTracks[kUncorrectedTrack] = loadOrNot;
}

// Setter for loading all track histograms
void EECHistogramManager::SetLoadAllTracks(const bool loadTracks, const bool loadUncorrected){
  SetLoadTracks(loadTracks);
  SetLoadTracksUncorrected(loadUncorrected);
}

// Setter for loading energy-energy correlators
void EECHistogramManager::SetLoadEnergyEnergyCorrelators(const bool loadOrNot){
  fLoadEnergyEnergyCorrelatorHistograms[kEnergyEnergyCorrelator] = loadOrNot;
}

// Setter for loading jet pT weighted energy-energy correlators
void EECHistogramManager::SetLoadEnergyEnergyCorrelatorsJetPt(const bool loadOrNot){
  fLoadEnergyEnergyCorrelatorHistograms[kEnergyEnergyCorrelatorJetPt] = loadOrNot;
}

// Setter for loading energy-energy correlators without tracking efficiency
void EECHistogramManager::SetLoadEnergyEnergyCorrelatorsUncorrected(const bool loadOrNot){
  fLoadEnergyEnergyCorrelatorHistograms[kEnergyEnergyCorrelatorUncorrected] = loadOrNot;
}

// Setter for loading jet pT weighted energy-energy correlators without tracking efficiency
void EECHistogramManager::SetLoadEnergyEnergyCorrelatorsJetPtUncorrected(const bool loadOrNot){
  fLoadEnergyEnergyCorrelatorHistograms[kEnergyEnergyCorrelatorJetPtUncorrected] = loadOrNot;
}

// Setter for loading all energy-energy correlators
void EECHistogramManager::SetLoadAllEnergyEnergyCorrelators(const bool loadRegular, const bool loadJetPt, const bool loadUncorrected, const bool loadJetPtUncorrected){
  SetLoadEnergyEnergyCorrelators(loadRegular);
  SetLoadEnergyEnergyCorrelatorsJetPt(loadJetPt);
  SetLoadEnergyEnergyCorrelatorsUncorrected(loadUncorrected);
  SetLoadEnergyEnergyCorrelatorsJetPtUncorrected(loadJetPtUncorrected);
}

// Setter for flavor selection for jets in jet histograms
void EECHistogramManager::SetJetFlavor(const int iFlavor){
  fJetFlavor = iFlavor;
}

 // Setter for loading two-dimensional histograms
void EECHistogramManager::SetLoad2DHistograms(const bool loadOrNot){
  fLoad2DHistograms = loadOrNot;
}

// Setter for loading jet pT closure histograms
void EECHistogramManager::SetLoadJetPtClosureHistograms(const bool loadOrNot){
  fLoadJetPtClosureHistograms = loadOrNot;
}

// Setter for loaded centrality bins
void EECHistogramManager::SetCentralityBinRange(const int first, const int last){
  fFirstLoadedCentralityBin = first;
  fLastLoadedCentralityBin = last;
  
  // Sanity check for centrality bins
  BinSanityCheck(fnCentralityBins,fFirstLoadedCentralityBin,fLastLoadedCentralityBin);
}

// Setter for loaded track pT bins
void EECHistogramManager::SetTrackPtBinRange(const int first, const int last){
  fFirstLoadedTrackPtBin = first;
  fLastLoadedTrackPtBin = last;
  
  // Sanity check for track pT bins
  BinSanityCheck(fnTrackPtBins,fFirstLoadedTrackPtBin,fLastLoadedTrackPtBin);
}

// Setter for jet pT bin range in energy-energy correlator histograms
void EECHistogramManager::SetJetPtBinRangeEEC(const int first, const int last){
  fFirstLoadedJetPtBinEEC = first;
  fLastLoadedJetPtBinEEC = last;
  
  // Sanity check for jet pT bins in energy-energy correlator histograms
  BinSanityCheck(fnJetPtBinsEEC,fFirstLoadedJetPtBinEEC,fLastLoadedJetPtBinEEC);
}

// Setter for track pT bin range in energy-energy correlator histograms
void EECHistogramManager::SetTrackPtBinRangeEEC(const int first, const int last){
  fFirstLoadedTrackPtBinEEC = first;
  fLastLoadedTrackPtBinEEC = last;
  
  // Sanity check for track pT bins in energy-energy correlator histograms
  BinSanityCheck(fnTrackPtBinsEEC,fFirstLoadedTrackPtBinEEC,fLastLoadedTrackPtBinEEC);
}

// Sanity check for set bins
void EECHistogramManager::BinSanityCheck(const int nBins, int& first, int& last){
  if(first < 0) first = 0;
  if(last < first) last = first;
  if(last > nBins-1) last = nBins-1;
}

// Sanity check for input bin index
int EECHistogramManager::BinIndexCheck(const int nBins, const int binIndex) const{
  if(binIndex < 0) return 0;
  if(binIndex > nBins-1) return nBins-1;
  return binIndex;
}

// Getter for the number of centrality bins
int EECHistogramManager::GetNCentralityBins() const{
  return fnCentralityBins;
}

// Getter for the number of track pT bins
int EECHistogramManager::GetNTrackPtBins() const{
  return fnTrackPtBins;
}

// Getter for the number of jet pT bins in energy-energy correlator histograms
int EECHistogramManager::GetNJetPtBinsEEC() const{
  return fnJetPtBinsEEC;
}

// Getter for the number of track pT bins in energy-energy correlator histograms
int EECHistogramManager::GetNTrackPtBinsEEC() const{
  return fnTrackPtBinsEEC;
}

// Getter for track histogram name
const char* EECHistogramManager::GetTrackHistogramName(int iTrackType) const{
  iTrackType = BinIndexCheck(knTrackCategories,iTrackType);
  return fTrackHistogramNames[iTrackType];
}

// Getter for name suitable for x-axis in a given track histogram
const char* EECHistogramManager::GetTrackAxisName(int iTrackType) const{
  iTrackType = BinIndexCheck(knTrackCategories,iTrackType);
  return fTrackAxisNames[iTrackType];
}

// Getter for the jet histogram name
const char* EECHistogramManager::GetJetHistogramName() const{
  return fJetHistogramName;
}

// Getter for name suitable for x-axis in a given jet histogram
const char* EECHistogramManager::GetJetAxisName() const{
  return fJetAxisName;
}

// Getter for energy-energy correlator histogram name
const char* EECHistogramManager::GetEnergyEnergyCorrelatorHistogramName(int iEnergyEnergyCorrelatorType) const{
  iEnergyEnergyCorrelatorType = BinIndexCheck(knEnergyEnergyCorrelatorTypes, iEnergyEnergyCorrelatorType);
  return fEnergyEnergyCorrelatorHistogramNames[iEnergyEnergyCorrelatorType];
}

// Getter for energy-energy correlator axis name
const char* EECHistogramManager::GetEnergyEnergyCorrelatorAxisName(int iEnergyEnergyCorrelatorType) const{
  iEnergyEnergyCorrelatorType = BinIndexCheck(knEnergyEnergyCorrelatorTypes, iEnergyEnergyCorrelatorType);
  return fEnergyEnergyCorrelatorAxisNames[iEnergyEnergyCorrelatorType];
}

// Getter for subevent types
const char* EECHistogramManager::GetSubeventType(const int iSubeventType) const{
  if(iSubeventType < 0) return "All combinations";
  if(iSubeventType >= knSubeventTypes) return "All combinations";
  return fSubeventTypeName[iSubeventType];
}

// Getter for pairing type save names
const char* EECHistogramManager::GetPairingTypeSaveName(const int iPairingType) const{
  if(iPairingType < 0) return "NonsensicalIndex";
  if(iPairingType >= EECHistograms::knPairingTypes) return "NonsensicalIndex";
  return fPairingTypeSaveName[iPairingType];
}

// Getter for collision system
TString EECHistogramManager::GetSystem() const{
  return fCard->GetDataType();
}

// Getter for i:th centrality bin border
double EECHistogramManager::GetCentralityBinBorder(const int iCentrality) const{
  return fCentralityBinBorders[iCentrality];
}

// Getter for i:th track pT bin border
double EECHistogramManager::GetTrackPtBinBorder(const int iTrackPt) const{
  return fTrackPtBinBorders[iTrackPt];
}

// Getter for i:th jet pT bin border in energy-energy correlator histograms
double EECHistogramManager::GetJetPtBinBorderEEC(const int iJetPt) const{
  return fJetPtBinBordersEEC[iJetPt];
}

// Getter for i:th track pT bin border in energy-energy correlator histograms
double EECHistogramManager::GetTrackPtBinBorderEEC(const int iTrackPt) const{
  return fTrackPtBinBordersEEC[iTrackPt];
}

// Getters for event information histograms

// Getter for z-vertex histogram
TH1D* EECHistogramManager::GetHistogramVertexZ() const{
  return fhVertexZ;
}

// Getter for z-vertex histogram
TH1D* EECHistogramManager::GetHistogramVertexZWeighted() const{
  return fhVertexZWeighted;
}

// Getter for histogram for number of events surviving different event cuts
TH1D* EECHistogramManager::GetHistogramEvents() const{
  return fhEvents;
}

// Getter for histogram for number of tracks surviving different track cuts
TH1D* EECHistogramManager::GetHistogramTrackCuts() const{
  return fhTrackCuts;
}

// Getter for centrality histogram in all events
TH1D* EECHistogramManager::GetHistogramCentrality() const{
  return fhCentrality;
}

// Getter for weighted centrality histogram in all events
TH1D* EECHistogramManager::GetHistogramCentralityWeighted() const{
  return fhCentralityWeighted;
}

// Getter for multiplicity histogram in all events
TH1D* EECHistogramManager::GetHistogramMultiplicity(int iCentrality) const{
  if(fCard->GetDataType().Contains("pp",TString::kIgnoreCase)) iCentrality = 0;  // No centrality selection for pp
  return fhMultiplicity[iCentrality];
}

// Getter for track efficiency weighted multiplicity histogram in all events
TH1D* EECHistogramManager::GetHistogramMultiplicityWeighted(int iCentrality) const{
  if(fCard->GetDataType().Contains("pp",TString::kIgnoreCase)) iCentrality = 0;  // No centrality selection for pp
  return fhMultiplicityWeighted[iCentrality];
}

// Getter for multiplicity vs. centrality map
TH2D* EECHistogramManager::GetHistogramMultiplicityMap() const{
  return fhMultiplicityMap;
}

// Getter for multiplicity vs. centrality map
TH2D* EECHistogramManager::GetHistogramWeightedMultiplicityMap() const{
  return fhMultiplicityMapWeighted;
}

// Getters for jet histograms

// Getter for jet pT histograms
TH1D* EECHistogramManager::GetHistogramJetPt(int iCentrality) const{
  if(fCard->GetDataType().Contains("pp",TString::kIgnoreCase)) iCentrality = 0;  // No centrality selection for pp
  return fhJetPt[iCentrality];
}

// Getter for jet phi histograms
TH1D* EECHistogramManager::GetHistogramJetPhi(int iCentrality) const{
  if(fCard->GetDataType().Contains("pp",TString::kIgnoreCase)) iCentrality = 0;  // No centrality selection for pp
  return fhJetPhi[iCentrality];
}

// Getter for jet eta histograms
TH1D* EECHistogramManager::GetHistogramJetEta(int iCentrality) const{
  if(fCard->GetDataType().Contains("pp",TString::kIgnoreCase)) iCentrality = 0;  // No centrality selection for pp
  return fhJetEta[iCentrality];
}

// Getter for 2D eta-phi histogram for jets
TH2D* EECHistogramManager::GetHistogramJetEtaPhi(int iCentrality) const{
  if(fCard->GetDataType().Contains("pp",TString::kIgnoreCase)) iCentrality = 0;  // No centrality selection for pp
  return fhJetEtaPhi[iCentrality];
}

// Getters for histograms for tracks

// Getter for track pT histograms
TH1D* EECHistogramManager::GetHistogramTrackPt(const int iTrackType, const int iCentrality) const{
  return fhTrackPt[iTrackType][iCentrality];
}

// Getter for track phi histograms
TH1D* EECHistogramManager::GetHistogramTrackPhi(const int iTrackType, const int iCentrality, const int iTrackPt) const{
  return fhTrackPhi[iTrackType][iCentrality][iTrackPt];
}

// Getter for track eta histograms
TH1D* EECHistogramManager::GetHistogramTrackEta(const int iTrackType, const int iCentrality, const int iTrackPt) const{
  return fhTrackEta[iTrackType][iCentrality][iTrackPt];
}

// Getter for 2D eta-phi histogram for track
TH2D* EECHistogramManager::GetHistogramTrackEtaPhi(const int iTrackType, const int iCentrality, const int iTrackPt) const{
  return fhTrackEtaPhi[iTrackType][iCentrality][iTrackPt];
}

// Getter for energy-energy correlator histograms
TH1D* EECHistogramManager::GetHistogramEnergyEnergyCorrelator(const int iEnergyEnergyCorrelatorType, const int iCentrality, const int iJetPt, const int iTrackPt, const int iPairingType, const int iSubevent) const{
  return fhEnergyEnergyCorrelator[iEnergyEnergyCorrelatorType][iCentrality][iJetPt][iTrackPt][iPairingType][iSubevent];
}

// Getter for jet pT closure histograms
TH1D* EECHistogramManager::GetHistogramJetPtClosure(const int iGenPtBin, const int iEtaBin, const int iCentrality, const int iClosureParticle) const{
  return fhJetPtClosure[iGenPtBin][iEtaBin][iCentrality][iClosureParticle];
}

/*
 * Getter for any one-dimensional histogram based on input string
 *
 *  Arguments:
 *   TString name = Name corresponding to loaded histogram
 *   int bin1 = First bin index for the histogram
 *   int bin2 = Second bin index for the histogram
 *   int bin3 = Third bin index for the histogram
 *   int bin4 = Fourth bin index for the histogram
 *   int bin5 = Fifth bin index for the histogram
 *   int bin6 = Sixth bin index for the histogram
 *
 *   return: Histogram corresponding to given name and bins
 */
TH1D* EECHistogramManager::GetOneDimensionalHistogram(TString name, int bin1, int bin2, int bin3, int bin4, int bin5, int bin6) const{
  if(name.EqualTo("vertexz",TString::kIgnoreCase) || name.EqualTo("fhvertexz",TString::kIgnoreCase)) return GetHistogramVertexZ();
  if(name.EqualTo("events",TString::kIgnoreCase) || name.EqualTo("fhevents",TString::kIgnoreCase)) return GetHistogramEvents();
  if(name.EqualTo("trackcuts",TString::kIgnoreCase) || name.EqualTo("fhtrackcuts",TString::kIgnoreCase)) return GetHistogramTrackCuts();
  if(name.EqualTo("centrality",TString::kIgnoreCase) || name.EqualTo("fhcentrality",TString::kIgnoreCase)) return GetHistogramCentrality();
  if(name.EqualTo("centralityweighted",TString::kIgnoreCase) || name.EqualTo("fhcentralityweighted",TString::kIgnoreCase)) return GetHistogramCentralityWeighted();
  if(name.EqualTo("multiplicity",TString::kIgnoreCase) || name.EqualTo("fhmultiplicity",TString::kIgnoreCase)) return GetHistogramMultiplicity(bin1);
  if(name.EqualTo("multiplicityweighted",TString::kIgnoreCase) || name.EqualTo("fhmultiplicityweighted",TString::kIgnoreCase)) return GetHistogramMultiplicityWeighted(bin1);
  if(name.EqualTo("jetpt",TString::kIgnoreCase) || name.EqualTo("fhjetpt",TString::kIgnoreCase)) return GetHistogramJetPt(bin1);
  if(name.EqualTo("jetphi",TString::kIgnoreCase) || name.EqualTo("fhjetphi",TString::kIgnoreCase)) return GetHistogramJetPhi(bin1);
  if(name.EqualTo("jeteta",TString::kIgnoreCase) || name.EqualTo("fhjeteta",TString::kIgnoreCase)) return GetHistogramJetEta(bin1);
  if(name.EqualTo("trackpt",TString::kIgnoreCase) || name.EqualTo("fhtrackpt",TString::kIgnoreCase)) return GetHistogramTrackPt(bin1,bin2);
  if(name.EqualTo("trackphi",TString::kIgnoreCase) || name.EqualTo("fhtrackphi",TString::kIgnoreCase)) return GetHistogramTrackPhi(bin1,bin2,bin3);
  if(name.EqualTo("tracketa",TString::kIgnoreCase) || name.EqualTo("fhtracketa",TString::kIgnoreCase)) return GetHistogramTrackEta(bin1,bin2,bin3);
  if(name.EqualTo("energyenergycorrelator",TString::kIgnoreCase) || name.EqualTo("fhenergyenergycorrelator",TString::kIgnoreCase) || name.EqualTo("eec",TString::kIgnoreCase)) return GetHistogramEnergyEnergyCorrelator(bin1,bin2,bin3,bin4,bin5,bin6);
  return NULL;
}

/*
 * Getter for any two-dimensional histogram based on input string
 *
 *  Arguments:
 *   TString name = Name corresponding to loaded histogram
 *   int bin1 = First bin index for the histogram
 *   int bin2 = Second bin index for the histogram
 *   int bin3 = Third bin index for the histogram
 *   int bin4 = Fourth bin index for the histogram
 *   int bin5 = Fifth bin index for the histogram
 *
 *   return: Histogram corresponding to given name and bins
 */
TH2D* EECHistogramManager::GetTwoDimensionalHistogram(TString name, int bin1, int bin2, int bin3, int bin4, int bin5) const{
  if(name.EqualTo("jetetaphi",TString::kIgnoreCase) || name.EqualTo("fhjetetaphi",TString::kIgnoreCase)) return GetHistogramJetEtaPhi(bin1);
  if(name.EqualTo("tracketaphi",TString::kIgnoreCase) || name.EqualTo("fhtracketaphi",TString::kIgnoreCase)) return GetHistogramTrackEtaPhi(bin1,bin2,bin3);
  if(name.EqualTo("multiplicitymap",TString::kIgnoreCase) || name.EqualTo("fhmultiplicitymap",TString::kIgnoreCase)) return GetHistogramMultiplicityMap();
  if(name.EqualTo("multiplicitymapweighted",TString::kIgnoreCase) || name.EqualTo("fhmultiplicitymapweighted",TString::kIgnoreCase)) return GetHistogramMultiplicityMap();
  return NULL;
}

// Get the first loaded centrality bin
int EECHistogramManager::GetFirstCentralityBin() const{
  return fFirstLoadedCentralityBin;
}

// Get the last loaded centrality bin
int EECHistogramManager::GetLastCentralityBin() const{
  return fLastLoadedCentralityBin;
}

// Get the first loaded track pT bin
int EECHistogramManager::GetFirstTrackPtBin() const{
  return fFirstLoadedTrackPtBin;
}

// Get the last loaded track pT bin
int EECHistogramManager::GetLastTrackPtBin() const{
  return fLastLoadedTrackPtBin;
}

// Get the first loaded energy-energy correlator jet pT bin
int EECHistogramManager::GetFirstJetPtBinEEC() const{
  return fFirstLoadedJetPtBinEEC;
}

// Get the last loaded energy-energy correlator jet pT bin
int EECHistogramManager::GetLastJetPtBinEEC() const{
  return fLastLoadedJetPtBinEEC;
}

// Get the first loaded energy-energy correlator track pT bin
int EECHistogramManager::GetFirstTrackPtBinEEC() const{
  return fFirstLoadedTrackPtBinEEC;
}

// Get the last loaded energy-energy correlator track pT bin
int EECHistogramManager::GetLastTrackPtBinEEC() const{
  return fLastLoadedTrackPtBinEEC;
}

// Getter for the number of events passing the cuts
int EECHistogramManager::GetNEvents() const{
  return fhEvents->GetBinContent(fhEvents->FindBin(EECHistograms::kVzCut));
}

// Getter for the JCard
EECCard* EECHistogramManager::GetCard() const{
  return fCard;
}

// Getter for integral over inclusive jet pT. Include the overflow bin in the integral.
double EECHistogramManager::GetJetPtIntegral(const int iCentrality) const{
  return fhJetPt[iCentrality]->Integral(1,fhJetPt[iCentrality]->GetNbinsX()+1,"width");
}

/*
 * Getter for integral over inclusive jet pT over specified range
 *
 *  const int iCentrality = Centrality bin
 *  const double minPt = Lower pT range for integral calculation
 *  const double maxPt = Higher pT range for integral calculation
 */
double EECHistogramManager::GetJetPtIntegral(const int iCentrality, const double minPt, const double maxPt) const{
  return fhJetPt[iCentrality]->Integral(fhJetPt[iCentrality]->FindBin(minPt+0.001), fhJetPt[iCentrality]->FindBin(maxPt-0.001), "width");
}
