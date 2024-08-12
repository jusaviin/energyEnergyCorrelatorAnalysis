# Energy-energy correlators in PbPb and pp collisions

## General structure of the reposity

The code is organized in folders, with different functionality in different folders. Below is a quick overview of this folder structure.

```Home folder```

This folder contains helpful processing scripts, together with a Makefile that compiles the code in the `src` folder that calculates the energy-energy correlators. Also the main program for the analysis, `eecAnalysis.cxx`, is located in this folder.

```src```

The code that calculates the energy-energy correlators is located in this folder.

```plotting```

Code that does post-processing for energy-energy correlator distributions and different plotting macros are located in this folder.

```crab```

The configuration to run the code on crab can be found from this folder.

```jetEnergyCorrections```

The jet energy corrections applied in this analysis are read from this folder.

```smearingFiles```

This folder contains particle level response matrices that can be used to smear the DeltaR values or energy weights to show the effect of detector resolution.

```trackCorrectionTables```

Corrections for tracks are read from the correction tables located in this folder.

```mixingFileList```

This folder contains file lists for the files from which the mixed events are read.

```skimmer```

This folder contains code to skim the minimum bias files from which the mixed events are read, leaving only bare minimum amount of information to perform the mixing. This greatly enhances the performance of the mixing.

```toySimulation```

Files to perform a toy simulation to study how energy-energy correlators look in completely uniform particle distribution.
