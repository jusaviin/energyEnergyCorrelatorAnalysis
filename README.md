# Energy-energy correlators in PbPb and pp collisions

## Producing final plots

If you just want to generate all the plots that are in the paper, you can use the script `makeAllPaperPlots.sh` provided you have access to all necessary data files. Contact Jussi Viinikainen to get these files.

## General structure of the reposity

The code is organized in folders, with different functionality in different folders. Below is a quick overview of this folder structure.

**Home folder**

This folder contains helpful processing scripts, together with a Makefile that compiles the code in the `src` folder that calculates the energy-energy correlators. Also the main program for the analysis, `eecAnalysis.cxx`, is located in this folder.

**src**

The code that calculates the energy-energy correlators is located in this folder.

**plotting**

Code that does post-processing for energy-energy correlator distributions and different plotting macros are located in this folder.

**crab**

The configuration to run the code on crab can be found from this folder.

**jetEnergyCorrections**

The jet energy corrections applied in this analysis are read from this folder.

**smearingFiles**

This folder contains particle level response matrices that can be used to smear the DeltaR values or energy weights to show the effect of detector resolution.

**trackCorrectionTables**

Corrections for tracks are read from the correction tables located in this folder.

**mixingFileList**

This folder contains file lists for the files from which the mixed events are read.

**skimmer**

This folder contains code to skim the minimum bias files from which the mixed events are read, leaving only bare minimum amount of information to perform the mixing. This greatly enhances the performance of the mixing.

**toySimulation**

Files to perform a toy simulation to study how energy-energy correlators look in completely uniform particle distribution.

## Overview of basic workflow

1. Test locally that everything works

Before sending any code to GRID processing, you should check that the code compiles and produces the histograms you want. First, compile the code

```
make
```

Then, you should check that the desired running configuration is in a card file. Wheck this with

```
vim cardEECPbPb.input
```

Newxt, run the code and see that the histograms you want are filled. Executing the main program without parameters tells you which parameters you need to give to it in order to run it.

```
./eecAnalysis
```

Once you have tested that everything works locally, it is time to run over larger datasets on CRAB.

2. Run the jobs on CRAB

For CRAB running, you will need a CMSSW area in your `lxplus` account. Any version will do, since no CMSSW functionality besised the CRAB job submitter will be used. Once you have a CMSSW area set up, copy the analysis code there and run it following the instructions in the `crab` folder in this repository. Once the jobs are finished, merge and download the produced root files.

3. Post-process the downloaded files

There are several post-processing steps before you can get the final energy-energy correlator distributions.
