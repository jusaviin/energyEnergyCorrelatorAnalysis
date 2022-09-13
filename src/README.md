# Code to construct correlations

This folder has the code that constructs correlations from forests. The code in compiled with the Makefile in the parent directory. Below is given a brief explanation what each class does.

**EECAnalyzer.cxx**

The main correlation code can be find from here. The code loops over a list of files given into it. It goes through all the events in the files, and analyzes them. It fills track, inclusive jet and energy-energy correlator histograms. Most of the configuration for the analysis can be read from an configuration card.

**EECHistograms.cxx**

All the histograms that are filled by EECAnalyzer are managed by this class. If you ever need to change binning, add histograms or anything else regarding the histograms, it can be done in this class. THnSparses are used for most of the analysis histograms, since they can conveniently contain several binning options.

**ForestReader.cxx**

A class template defining what kind of functionality a class reading the trees from a forest must implement. There is a separate class reading the trees since the analysis should not depend on the exact implementation of the forest. If the forest changes, one should write a new class inheriting from the ForestReader class reading it and no changes should be made in the main analysis code.

**HighForestReader.cxx**

A forest reader that can read the standard CMS heavy ion forests.

**GeneratorLevelForestReader.cxx**

A forest reader that reads generator level jet and track information instead of reconstructed values.

**ConfigurationCard.cxx**

A class that extracts information from specifically formatted text file (.input) and stores it. This allows to change input parameters for the analysis without touching the source code. It also creates a reference on the parameter set that is used for each run that is stored in the root file with the correlations.

**JetCorrector.cxx**

Class for doing jet energy correction on the fly. The code is originally from https://twiki.cern.ch/twiki/pub/CMS/HiJetReco2019/JetCorrector.h. It has been separated to cxx and h files and some debugging options have been added.

**JetUncertainty.cxx**

Class for estimating an uncertainty for the jet energy. The code is originally from https://twiki.cern.ch/twiki/pub/CMS/HiJetReco2019/JetUncertainty.h. It has been separated to cxx and h files and some debugging options have been added.

**TrackingEfficiencyInterface.cxx**

Class interface that defines what functions tracking efficiency classes for new data must implement. Used to seemslessly switch between pp and PbPb tracking corrections in DijetAnalyzer.

**trackingEfficiency2018PbPb.cxx**

Class managing tracking efficiency correction for 2018 PbPb data. The code is originally from https://twiki.cern.ch/twiki/bin/viewauth/CMS/HITracking2018PbPb. It has been slightly modified for this analysis.

**trackingEfficiency2017pp.cxx**

Class managing tracking efficiency correction for 2017 pp data. The code is originally from https://twiki.cern.ch/twiki/bin/view/CMS/HiTracking2017pp. It has been slightly modified for this analysis.
