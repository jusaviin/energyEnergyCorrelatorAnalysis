#!/bin/bash

# Check that the necessary folders exist and create them if not
if [ ! -d "figures" ]; then
  mkdir figures
fi

if [ ! -d "paperFigures" ]; then
  mkdir paperFigures
fi

# First, run the finalResultPlotter.C macro to create big canvases with only data
root -l -b -q 'plotting/finalResultPlotter.C(true)'

# From the created figures, move the ones that are used in the paper to paper plots folder
mv figures/finalBigDualDistributionCanvas_T\>1v0.pdf paperFigures/Figure_002.pdf
mv figures/finalBigDualDistributionCanvas_T\>2v0.pdf paperFigures/Figure_003.pdf
mv figures/finalBigDualRatioCanvas_T\>1v0.pdf paperFigures/Figure_004.pdf
mv figures/finalBigDualRatioCanvas_T\>2v0.pdf paperFigures/Figure_005.pdf

# The rest of the figures are produced by the modelComparisons.C macro

# First, produce all distributions with Hybrid model comparisons
root -l -b -q 'plotting/modelComparison.C(1,0,0,true)'
root -l -b -q 'plotting/modelComparison.C(2,0,0,true)'

# Move the figures used in the paper to paperFigures folder
cp figures/energyEnergyCorrelator_distributionModelComparison_nominalEnergyWeight_pp_J\=120-140_T\>1v0.pdf paperFigures/Figure_006-a.pdf
cp figures/energyEnergyCorrelator_distributionModelComparison_energyWeightSquared_pp_J\=120-140_T\>1v0.pdf paperFigures/Figure_006-b.pdf

mv figures/energyEnergyCorrelator_distributionModelComparison_nominalEnergyWeight_pp_J\=120-140_T\>1v0.pdf paperFigures/Figure_A01-a.pdf
mv figures/energyEnergyCorrelator_distributionModelComparison_nominalEnergyWeight_pp_J\=120-140_T\>2v0.pdf paperFigures/Figure_A01-b.pdf
mv figures/energyEnergyCorrelator_distributionModelComparison_nominalEnergyWeight_pp_J\=140-160_T\>1v0.pdf paperFigures/Figure_A01-c.pdf
mv figures/energyEnergyCorrelator_distributionModelComparison_nominalEnergyWeight_pp_J\=140-160_T\>2v0.pdf paperFigures/Figure_A01-d.pdf
mv figures/energyEnergyCorrelator_distributionModelComparison_nominalEnergyWeight_pp_J\=160-180_T\>1v0.pdf paperFigures/Figure_A01-e.pdf
mv figures/energyEnergyCorrelator_distributionModelComparison_nominalEnergyWeight_pp_J\=160-180_T\>2v0.pdf paperFigures/Figure_A01-f.pdf
mv figures/energyEnergyCorrelator_distributionModelComparison_nominalEnergyWeight_pp_J\=180-200_T\>1v0.pdf paperFigures/Figure_A01-g.pdf
mv figures/energyEnergyCorrelator_distributionModelComparison_nominalEnergyWeight_pp_J\=180-200_T\>2v0.pdf paperFigures/Figure_A01-h.pdf

mv figures/energyEnergyCorrelator_distributionModelComparison_energyWeightSquared_pp_J\=120-140_T\>1v0.pdf paperFigures/Figure_A02-a.pdf
mv figures/energyEnergyCorrelator_distributionModelComparison_energyWeightSquared_pp_J\=120-140_T\>2v0.pdf paperFigures/Figure_A02-b.pdf
mv figures/energyEnergyCorrelator_distributionModelComparison_energyWeightSquared_pp_J\=140-160_T\>1v0.pdf paperFigures/Figure_A02-c.pdf
mv figures/energyEnergyCorrelator_distributionModelComparison_energyWeightSquared_pp_J\=140-160_T\>2v0.pdf paperFigures/Figure_A02-d.pdf
mv figures/energyEnergyCorrelator_distributionModelComparison_energyWeightSquared_pp_J\=160-180_T\>1v0.pdf paperFigures/Figure_A02-e.pdf
mv figures/energyEnergyCorrelator_distributionModelComparison_energyWeightSquared_pp_J\=160-180_T\>2v0.pdf paperFigures/Figure_A02-f.pdf
mv figures/energyEnergyCorrelator_distributionModelComparison_energyWeightSquared_pp_J\=180-200_T\>1v0.pdf paperFigures/Figure_A02-g.pdf
mv figures/energyEnergyCorrelator_distributionModelComparison_energyWeightSquared_pp_J\=180-200_T\>2v0.pdf paperFigures/Figure_A02-h.pdf

mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=50-90_J\=120-140_T\>1v0.pdf paperFigures/Figure_A05-a.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=50-90_J\=120-140_T\>2v0.pdf paperFigures/Figure_A05-b.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=50-90_J\=140-160_T\>1v0.pdf paperFigures/Figure_A05-c.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=50-90_J\=140-160_T\>2v0.pdf paperFigures/Figure_A05-d.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=50-90_J\=160-180_T\>1v0.pdf paperFigures/Figure_A05-e.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=50-90_J\=160-180_T\>2v0.pdf paperFigures/Figure_A05-f.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=50-90_J\=180-200_T\>1v0.pdf paperFigures/Figure_A05-g.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=50-90_J\=180-200_T\>2v0.pdf paperFigures/Figure_A05-h.pdf

mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=50-90_J\=120-140_T\>1v0.pdf paperFigures/Figure_A06-a.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=50-90_J\=120-140_T\>2v0.pdf paperFigures/Figure_A06-b.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=50-90_J\=140-160_T\>1v0.pdf paperFigures/Figure_A06-c.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=50-90_J\=140-160_T\>2v0.pdf paperFigures/Figure_A06-d.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=50-90_J\=160-180_T\>1v0.pdf paperFigures/Figure_A06-e.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=50-90_J\=160-180_T\>2v0.pdf paperFigures/Figure_A06-f.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=50-90_J\=180-200_T\>1v0.pdf paperFigures/Figure_A06-g.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=50-90_J\=180-200_T\>2v0.pdf paperFigures/Figure_A06-h.pdf

mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=30-50_J\=120-140_T\>1v0.pdf paperFigures/Figure_A07-a.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=30-50_J\=120-140_T\>2v0.pdf paperFigures/Figure_A07-b.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=30-50_J\=140-160_T\>1v0.pdf paperFigures/Figure_A07-c.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=30-50_J\=140-160_T\>2v0.pdf paperFigures/Figure_A07-d.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=30-50_J\=160-180_T\>1v0.pdf paperFigures/Figure_A07-e.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=30-50_J\=160-180_T\>2v0.pdf paperFigures/Figure_A07-f.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=30-50_J\=180-200_T\>1v0.pdf paperFigures/Figure_A07-g.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=30-50_J\=180-200_T\>2v0.pdf paperFigures/Figure_A07-h.pdf

mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=30-50_J\=120-140_T\>1v0.pdf paperFigures/Figure_A08-a.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=30-50_J\=120-140_T\>2v0.pdf paperFigures/Figure_A08-b.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=30-50_J\=140-160_T\>1v0.pdf paperFigures/Figure_A08-c.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=30-50_J\=140-160_T\>2v0.pdf paperFigures/Figure_A08-d.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=30-50_J\=160-180_T\>1v0.pdf paperFigures/Figure_A08-e.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=30-50_J\=160-180_T\>2v0.pdf paperFigures/Figure_A08-f.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=30-50_J\=180-200_T\>1v0.pdf paperFigures/Figure_A08-g.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=30-50_J\=180-200_T\>2v0.pdf paperFigures/Figure_A08-h.pdf

mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=10-30_J\=120-140_T\>1v0.pdf paperFigures/Figure_A09-a.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=10-30_J\=120-140_T\>2v0.pdf paperFigures/Figure_A09-b.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=10-30_J\=140-160_T\>1v0.pdf paperFigures/Figure_A09-c.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=10-30_J\=140-160_T\>2v0.pdf paperFigures/Figure_A09-d.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=10-30_J\=160-180_T\>1v0.pdf paperFigures/Figure_A09-e.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=10-30_J\=160-180_T\>2v0.pdf paperFigures/Figure_A09-f.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=10-30_J\=180-200_T\>1v0.pdf paperFigures/Figure_A09-g.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=10-30_J\=180-200_T\>2v0.pdf paperFigures/Figure_A09-h.pdf

mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=10-30_J\=120-140_T\>1v0.pdf paperFigures/Figure_A10-a.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=10-30_J\=120-140_T\>2v0.pdf paperFigures/Figure_A10-b.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=10-30_J\=140-160_T\>1v0.pdf paperFigures/Figure_A10-c.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=10-30_J\=140-160_T\>2v0.pdf paperFigures/Figure_A10-d.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=10-30_J\=160-180_T\>1v0.pdf paperFigures/Figure_A10-e.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=10-30_J\=160-180_T\>2v0.pdf paperFigures/Figure_A10-f.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=10-30_J\=180-200_T\>1v0.pdf paperFigures/Figure_A10-g.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=10-30_J\=180-200_T\>2v0.pdf paperFigures/Figure_A10-h.pdf

mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=0-10_J\=120-140_T\>1v0.pdf paperFigures/Figure_A11-a.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=0-10_J\=120-140_T\>2v0.pdf paperFigures/Figure_A11-b.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=0-10_J\=140-160_T\>1v0.pdf paperFigures/Figure_A11-c.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=0-10_J\=140-160_T\>2v0.pdf paperFigures/Figure_A11-d.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=0-10_J\=160-180_T\>1v0.pdf paperFigures/Figure_A11-e.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=0-10_J\=160-180_T\>2v0.pdf paperFigures/Figure_A11-f.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=0-10_J\=180-200_T\>1v0.pdf paperFigures/Figure_A11-g.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_nominalEnergyWeight_C\=0-10_J\=180-200_T\>2v0.pdf paperFigures/Figure_A11-h.pdf

mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=0-10_J\=120-140_T\>1v0.pdf paperFigures/Figure_A12-a.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=0-10_J\=120-140_T\>2v0.pdf paperFigures/Figure_A12-b.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=0-10_J\=140-160_T\>1v0.pdf paperFigures/Figure_A12-c.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=0-10_J\=140-160_T\>2v0.pdf paperFigures/Figure_A12-d.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=0-10_J\=160-180_T\>1v0.pdf paperFigures/Figure_A12-e.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=0-10_J\=160-180_T\>2v0.pdf paperFigures/Figure_A12-f.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=0-10_J\=180-200_T\>1v0.pdf paperFigures/Figure_A12-g.pdf
mv figures/energyEnergyCorrelator_hybridModelDistribution_energyWeightSquared_C\=0-10_J\=180-200_T\>2v0.pdf paperFigures/Figure_A12-h.pdf

# Next in the paper there are PbPb/pp ratios with Hybrid model comparisons
root -l -b -q 'plotting/modelComparison.C(1,0,1,true)'
root -l -b -q 'plotting/modelComparison.C(2,0,1,true)'

# Then, move all these plots to the paperFigures folder
cp figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=0-10_J\=120-140_T\>1v0.pdf paperFigures/Figure_007-a.pdf
cp figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=0-10_J\=180-200_T\>1v0.pdf paperFigures/Figure_007-b.pdf

mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=50-90_J\=120-140_T\>1v0.pdf paperFigures/Figure_A16-a.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=50-90_J\=120-140_T\>2v0.pdf paperFigures/Figure_A16-b.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=50-90_J\=140-160_T\>1v0.pdf paperFigures/Figure_A16-c.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=50-90_J\=140-160_T\>2v0.pdf paperFigures/Figure_A16-d.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=50-90_J\=160-180_T\>1v0.pdf paperFigures/Figure_A16-e.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=50-90_J\=160-180_T\>2v0.pdf paperFigures/Figure_A16-f.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=50-90_J\=180-200_T\>1v0.pdf paperFigures/Figure_A16-g.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=50-90_J\=180-200_T\>2v0.pdf paperFigures/Figure_A16-h.pdf

mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=50-90_J\=120-140_T\>1v0.pdf paperFigures/Figure_A17-a.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=50-90_J\=120-140_T\>2v0.pdf paperFigures/Figure_A17-b.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=50-90_J\=140-160_T\>1v0.pdf paperFigures/Figure_A17-c.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=50-90_J\=140-160_T\>2v0.pdf paperFigures/Figure_A17-d.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=50-90_J\=160-180_T\>1v0.pdf paperFigures/Figure_A17-e.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=50-90_J\=160-180_T\>2v0.pdf paperFigures/Figure_A17-f.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=50-90_J\=180-200_T\>1v0.pdf paperFigures/Figure_A17-g.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=50-90_J\=180-200_T\>2v0.pdf paperFigures/Figure_A17-h.pdf

mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=30-50_J\=120-140_T\>1v0.pdf paperFigures/Figure_A18-a.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=30-50_J\=120-140_T\>2v0.pdf paperFigures/Figure_A18-b.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=30-50_J\=140-160_T\>1v0.pdf paperFigures/Figure_A18-c.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=30-50_J\=140-160_T\>2v0.pdf paperFigures/Figure_A18-d.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=30-50_J\=160-180_T\>1v0.pdf paperFigures/Figure_A18-e.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=30-50_J\=160-180_T\>2v0.pdf paperFigures/Figure_A18-f.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=30-50_J\=180-200_T\>1v0.pdf paperFigures/Figure_A18-g.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=30-50_J\=180-200_T\>2v0.pdf paperFigures/Figure_A18-h.pdf

mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=30-50_J\=120-140_T\>1v0.pdf paperFigures/Figure_A19-a.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=30-50_J\=120-140_T\>2v0.pdf paperFigures/Figure_A19-b.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=30-50_J\=140-160_T\>1v0.pdf paperFigures/Figure_A19-c.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=30-50_J\=140-160_T\>2v0.pdf paperFigures/Figure_A19-d.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=30-50_J\=160-180_T\>1v0.pdf paperFigures/Figure_A19-e.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=30-50_J\=160-180_T\>2v0.pdf paperFigures/Figure_A19-f.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=30-50_J\=180-200_T\>1v0.pdf paperFigures/Figure_A19-g.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=30-50_J\=180-200_T\>2v0.pdf paperFigures/Figure_A19-h.pdf

mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=10-30_J\=120-140_T\>1v0.pdf paperFigures/Figure_A20-a.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=10-30_J\=120-140_T\>2v0.pdf paperFigures/Figure_A20-b.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=10-30_J\=140-160_T\>1v0.pdf paperFigures/Figure_A20-c.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=10-30_J\=140-160_T\>2v0.pdf paperFigures/Figure_A20-d.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=10-30_J\=160-180_T\>1v0.pdf paperFigures/Figure_A20-e.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=10-30_J\=160-180_T\>2v0.pdf paperFigures/Figure_A20-f.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=10-30_J\=180-200_T\>1v0.pdf paperFigures/Figure_A20-g.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=10-30_J\=180-200_T\>2v0.pdf paperFigures/Figure_A20-h.pdf

mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=10-30_J\=120-140_T\>1v0.pdf paperFigures/Figure_A21-a.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=10-30_J\=120-140_T\>2v0.pdf paperFigures/Figure_A21-b.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=10-30_J\=140-160_T\>1v0.pdf paperFigures/Figure_A21-c.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=10-30_J\=140-160_T\>2v0.pdf paperFigures/Figure_A21-d.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=10-30_J\=160-180_T\>1v0.pdf paperFigures/Figure_A21-e.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=10-30_J\=160-180_T\>2v0.pdf paperFigures/Figure_A21-f.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=10-30_J\=180-200_T\>1v0.pdf paperFigures/Figure_A21-g.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=10-30_J\=180-200_T\>2v0.pdf paperFigures/Figure_A21-h.pdf

mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=0-10_J\=120-140_T\>1v0.pdf paperFigures/Figure_A22-a.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=0-10_J\=120-140_T\>2v0.pdf paperFigures/Figure_A22-b.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=0-10_J\=140-160_T\>1v0.pdf paperFigures/Figure_A22-c.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=0-10_J\=140-160_T\>2v0.pdf paperFigures/Figure_A22-d.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=0-10_J\=160-180_T\>1v0.pdf paperFigures/Figure_A22-e.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=0-10_J\=160-180_T\>2v0.pdf paperFigures/Figure_A22-f.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=0-10_J\=180-200_T\>1v0.pdf paperFigures/Figure_A22-g.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_nominalEnergyWeight_C\=0-10_J\=180-200_T\>2v0.pdf paperFigures/Figure_A22-h.pdf

mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=0-10_J\=120-140_T\>1v0.pdf paperFigures/Figure_A23-a.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=0-10_J\=120-140_T\>2v0.pdf paperFigures/Figure_A23-b.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=0-10_J\=140-160_T\>1v0.pdf paperFigures/Figure_A23-c.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=0-10_J\=140-160_T\>2v0.pdf paperFigures/Figure_A23-d.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=0-10_J\=160-180_T\>1v0.pdf paperFigures/Figure_A23-e.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=0-10_J\=160-180_T\>2v0.pdf paperFigures/Figure_A23-f.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=0-10_J\=180-200_T\>1v0.pdf paperFigures/Figure_A23-g.pdf
mv figures/energyEnergyCorrelator_hybridModelRatio_energyWeightSquared_C\=0-10_J\=180-200_T\>2v0.pdf paperFigures/Figure_A23-h.pdf

# Next figures in the paper are the double ratios with Hybrid comparisons
root -l -b -q 'plotting/modelComparison.C(1,0,2,true)'
root -l -b -q 'plotting/modelComparison.C(2,0,2,true)'

# Then, move all these plots to the paperFigures folder
cp figures/energyEnergyCorrelator_hybridModelDoubleRatio_nominalEnergyWeight_C\=0-10_J\=120-140.pdf paperFigures/Figure_008-a.pdf
cp figures/energyEnergyCorrelator_hybridModelDoubleRatio_nominalEnergyWeight_C\=0-10_J\=180-200.pdf paperFigures/Figure_008-b.pdf

mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_nominalEnergyWeight_C\=50-90_J\=120-140.pdf paperFigures/Figure_A27-a.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_energyWeightSquared_C\=50-90_J\=120-140.pdf paperFigures/Figure_A27-b.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_nominalEnergyWeight_C\=50-90_J\=140-160.pdf paperFigures/Figure_A27-c.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_energyWeightSquared_C\=50-90_J\=140-160.pdf paperFigures/Figure_A27-d.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_nominalEnergyWeight_C\=50-90_J\=160-180.pdf paperFigures/Figure_A27-e.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_energyWeightSquared_C\=50-90_J\=160-180.pdf paperFigures/Figure_A27-f.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_nominalEnergyWeight_C\=50-90_J\=180-200.pdf paperFigures/Figure_A27-g.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_energyWeightSquared_C\=50-90_J\=180-200.pdf paperFigures/Figure_A27-h.pdf

mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_nominalEnergyWeight_C\=30-50_J\=120-140.pdf paperFigures/Figure_A28-a.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_energyWeightSquared_C\=30-50_J\=120-140.pdf paperFigures/Figure_A28-b.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_nominalEnergyWeight_C\=30-50_J\=140-160.pdf paperFigures/Figure_A28-c.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_energyWeightSquared_C\=30-50_J\=140-160.pdf paperFigures/Figure_A28-d.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_nominalEnergyWeight_C\=30-50_J\=160-180.pdf paperFigures/Figure_A28-e.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_energyWeightSquared_C\=30-50_J\=160-180.pdf paperFigures/Figure_A28-f.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_nominalEnergyWeight_C\=30-50_J\=180-200.pdf paperFigures/Figure_A28-g.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_energyWeightSquared_C\=30-50_J\=180-200.pdf paperFigures/Figure_A28-h.pdf

mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_nominalEnergyWeight_C\=10-30_J\=120-140.pdf paperFigures/Figure_A29-a.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_energyWeightSquared_C\=10-30_J\=120-140.pdf paperFigures/Figure_A29-b.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_nominalEnergyWeight_C\=10-30_J\=140-160.pdf paperFigures/Figure_A29-c.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_energyWeightSquared_C\=10-30_J\=140-160.pdf paperFigures/Figure_A29-d.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_nominalEnergyWeight_C\=10-30_J\=160-180.pdf paperFigures/Figure_A29-e.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_energyWeightSquared_C\=10-30_J\=160-180.pdf paperFigures/Figure_A29-f.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_nominalEnergyWeight_C\=10-30_J\=180-200.pdf paperFigures/Figure_A29-g.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_energyWeightSquared_C\=10-30_J\=180-200.pdf paperFigures/Figure_A29-h.pdf

mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_nominalEnergyWeight_C\=0-10_J\=120-140.pdf paperFigures/Figure_A30-a.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_energyWeightSquared_C\=0-10_J\=120-140.pdf paperFigures/Figure_A30-b.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_nominalEnergyWeight_C\=0-10_J\=140-160.pdf paperFigures/Figure_A30-c.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_energyWeightSquared_C\=0-10_J\=140-160.pdf paperFigures/Figure_A30-d.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_nominalEnergyWeight_C\=0-10_J\=160-180.pdf paperFigures/Figure_A30-e.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_energyWeightSquared_C\=0-10_J\=160-180.pdf paperFigures/Figure_A30-f.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_nominalEnergyWeight_C\=0-10_J\=180-200.pdf paperFigures/Figure_A30-g.pdf
mv figures/energyEnergyCorrelator_hybridModelDoubleRatio_energyWeightSquared_C\=0-10_J\=180-200.pdf paperFigures/Figure_A30-h.pdf

# The next figure in the paper has the perturbative calculation and JEWEL predictions for PbPb/pp ratio in the same plot
root -l -b -q 'plotting/modelComparison.C(1,6,1,true)'

# Move the relevant produced plots to the paperFigures folder
mv figures/energyEnergyCorrelator_holguinAndJewelRatio_nominalEnergyWeight_C\=0-10_J\=120-140_T\>1v0.pdf paperFigures/Figure_009-a.pdf
mv figures/energyEnergyCorrelator_holguinAndJewelRatio_nominalEnergyWeight_C\=0-10_J\=180-200_T\>1v0.pdf paperFigures/Figure_009-b.pdf

# Finally, in the main text of the paper we have the PbPb/pp ratio comparison with CoLBT model
root -l -b -q 'plotting/modelComparison.C(1,2,1,true)'

# Move the relevant plot to the paper figure folder
mv figures/energyEnergyCorrelator_colbtRatio_nominalEnergyWeight_C\=0-10_J\=120-140_T\>1v0.pdf paperFigures/Figure_010.pdf

# We still need to plot figures that are included in the appendix. First, make plots with pp and PbPb distributions compared with JEWEL predictions
root -l -b -q 'plotting/modelComparison.C(1,5,0,true)'
root -l -b -q 'plotting/modelComparison.C(2,5,0,true)'

# Move the produced figures to the paper figure folder
mv figures/energyEnergyCorrelator_distributionExtendedModelComparison_nominalEnergyWeight_pp_J\=120-140_T\>1v0.pdf paperFigures/Figure_A03-a.pdf
mv figures/energyEnergyCorrelator_distributionExtendedModelComparison_nominalEnergyWeight_pp_J\=120-140_T\>2v0.pdf paperFigures/Figure_A03-b.pdf
mv figures/energyEnergyCorrelator_distributionExtendedModelComparison_nominalEnergyWeight_pp_J\=140-160_T\>1v0.pdf paperFigures/Figure_A03-c.pdf
mv figures/energyEnergyCorrelator_distributionExtendedModelComparison_nominalEnergyWeight_pp_J\=140-160_T\>2v0.pdf paperFigures/Figure_A03-d.pdf
mv figures/energyEnergyCorrelator_distributionExtendedModelComparison_nominalEnergyWeight_pp_J\=160-180_T\>1v0.pdf paperFigures/Figure_A03-e.pdf
mv figures/energyEnergyCorrelator_distributionExtendedModelComparison_nominalEnergyWeight_pp_J\=160-180_T\>2v0.pdf paperFigures/Figure_A03-f.pdf
mv figures/energyEnergyCorrelator_distributionExtendedModelComparison_nominalEnergyWeight_pp_J\=180-200_T\>1v0.pdf paperFigures/Figure_A03-g.pdf
mv figures/energyEnergyCorrelator_distributionExtendedModelComparison_nominalEnergyWeight_pp_J\=180-200_T\>2v0.pdf paperFigures/Figure_A03-h.pdf

mv figures/energyEnergyCorrelator_distributionExtendedModelComparison_energyWeightSquared_pp_J\=120-140_T\>1v0.pdf paperFigures/Figure_A04-a.pdf
mv figures/energyEnergyCorrelator_distributionExtendedModelComparison_energyWeightSquared_pp_J\=120-140_T\>2v0.pdf paperFigures/Figure_A04-b.pdf
mv figures/energyEnergyCorrelator_distributionExtendedModelComparison_energyWeightSquared_pp_J\=140-160_T\>1v0.pdf paperFigures/Figure_A04-c.pdf
mv figures/energyEnergyCorrelator_distributionExtendedModelComparison_energyWeightSquared_pp_J\=140-160_T\>2v0.pdf paperFigures/Figure_A04-d.pdf
mv figures/energyEnergyCorrelator_distributionExtendedModelComparison_energyWeightSquared_pp_J\=160-180_T\>1v0.pdf paperFigures/Figure_A04-e.pdf
mv figures/energyEnergyCorrelator_distributionExtendedModelComparison_energyWeightSquared_pp_J\=160-180_T\>2v0.pdf paperFigures/Figure_A04-f.pdf
mv figures/energyEnergyCorrelator_distributionExtendedModelComparison_energyWeightSquared_pp_J\=180-200_T\>1v0.pdf paperFigures/Figure_A04-g.pdf
mv figures/energyEnergyCorrelator_distributionExtendedModelComparison_energyWeightSquared_pp_J\=180-200_T\>2v0.pdf paperFigures/Figure_A04-h.pdf

mv figures/energyEnergyCorrelator_jewelDistribution_nominalEnergyWeight_C\=0-10_J\=120-140_T\>1v0.pdf paperFigures/Figure_A14-a.pdf
mv figures/energyEnergyCorrelator_jewelDistribution_nominalEnergyWeight_C\=0-10_J\=120-140_T\>2v0.pdf paperFigures/Figure_A14-b.pdf
mv figures/energyEnergyCorrelator_jewelDistribution_nominalEnergyWeight_C\=0-10_J\=140-160_T\>1v0.pdf paperFigures/Figure_A14-c.pdf
mv figures/energyEnergyCorrelator_jewelDistribution_nominalEnergyWeight_C\=0-10_J\=140-160_T\>2v0.pdf paperFigures/Figure_A14-d.pdf
mv figures/energyEnergyCorrelator_jewelDistribution_nominalEnergyWeight_C\=0-10_J\=160-180_T\>1v0.pdf paperFigures/Figure_A14-e.pdf
mv figures/energyEnergyCorrelator_jewelDistribution_nominalEnergyWeight_C\=0-10_J\=160-180_T\>2v0.pdf paperFigures/Figure_A14-f.pdf
mv figures/energyEnergyCorrelator_jewelDistribution_nominalEnergyWeight_C\=0-10_J\=180-200_T\>1v0.pdf paperFigures/Figure_A14-g.pdf
mv figures/energyEnergyCorrelator_jewelDistribution_nominalEnergyWeight_C\=0-10_J\=180-200_T\>2v0.pdf paperFigures/Figure_A14-h.pdf

mv figures/energyEnergyCorrelator_jewelDistribution_energyWeightSquared_C\=0-10_J\=120-140_T\>1v0.pdf paperFigures/Figure_A15-a.pdf
mv figures/energyEnergyCorrelator_jewelDistribution_energyWeightSquared_C\=0-10_J\=120-140_T\>2v0.pdf paperFigures/Figure_A15-b.pdf
mv figures/energyEnergyCorrelator_jewelDistribution_energyWeightSquared_C\=0-10_J\=140-160_T\>1v0.pdf paperFigures/Figure_A15-c.pdf
mv figures/energyEnergyCorrelator_jewelDistribution_energyWeightSquared_C\=0-10_J\=140-160_T\>2v0.pdf paperFigures/Figure_A15-d.pdf
mv figures/energyEnergyCorrelator_jewelDistribution_energyWeightSquared_C\=0-10_J\=160-180_T\>1v0.pdf paperFigures/Figure_A15-e.pdf
mv figures/energyEnergyCorrelator_jewelDistribution_energyWeightSquared_C\=0-10_J\=160-180_T\>2v0.pdf paperFigures/Figure_A15-f.pdf
mv figures/energyEnergyCorrelator_jewelDistribution_energyWeightSquared_C\=0-10_J\=180-200_T\>1v0.pdf paperFigures/Figure_A15-g.pdf
mv figures/energyEnergyCorrelator_jewelDistribution_energyWeightSquared_C\=0-10_J\=180-200_T\>2v0.pdf paperFigures/Figure_A15-h.pdf

# Next, make the PbPb/pp ratios with only JEWEL results
root -l -b -q 'plotting/modelComparison.C(1,5,1,true)'
root -l -b -q 'plotting/modelComparison.C(2,5,1,true)'

# Move the relevant figures to the paper plots folder
mv figures/energyEnergyCorrelator_jewelRatio_nominalEnergyWeight_C\=0-10_J\=120-140_T\>1v0.pdf paperFigures/Figure_A25-a.pdf
mv figures/energyEnergyCorrelator_jewelRatio_nominalEnergyWeight_C\=0-10_J\=120-140_T\>2v0.pdf paperFigures/Figure_A25-b.pdf
mv figures/energyEnergyCorrelator_jewelRatio_nominalEnergyWeight_C\=0-10_J\=140-160_T\>1v0.pdf paperFigures/Figure_A25-c.pdf
mv figures/energyEnergyCorrelator_jewelRatio_nominalEnergyWeight_C\=0-10_J\=140-160_T\>2v0.pdf paperFigures/Figure_A25-d.pdf
mv figures/energyEnergyCorrelator_jewelRatio_nominalEnergyWeight_C\=0-10_J\=160-180_T\>1v0.pdf paperFigures/Figure_A25-e.pdf
mv figures/energyEnergyCorrelator_jewelRatio_nominalEnergyWeight_C\=0-10_J\=160-180_T\>2v0.pdf paperFigures/Figure_A25-f.pdf
mv figures/energyEnergyCorrelator_jewelRatio_nominalEnergyWeight_C\=0-10_J\=180-200_T\>1v0.pdf paperFigures/Figure_A25-g.pdf
mv figures/energyEnergyCorrelator_jewelRatio_nominalEnergyWeight_C\=0-10_J\=180-200_T\>2v0.pdf paperFigures/Figure_A25-h.pdf

mv figures/energyEnergyCorrelator_jewelRatio_energyWeightSquared_C\=0-10_J\=120-140_T\>1v0.pdf paperFigures/Figure_A26-a.pdf
mv figures/energyEnergyCorrelator_jewelRatio_energyWeightSquared_C\=0-10_J\=120-140_T\>2v0.pdf paperFigures/Figure_A26-b.pdf
mv figures/energyEnergyCorrelator_jewelRatio_energyWeightSquared_C\=0-10_J\=140-160_T\>1v0.pdf paperFigures/Figure_A26-c.pdf
mv figures/energyEnergyCorrelator_jewelRatio_energyWeightSquared_C\=0-10_J\=140-160_T\>2v0.pdf paperFigures/Figure_A26-d.pdf
mv figures/energyEnergyCorrelator_jewelRatio_energyWeightSquared_C\=0-10_J\=160-180_T\>1v0.pdf paperFigures/Figure_A26-e.pdf
mv figures/energyEnergyCorrelator_jewelRatio_energyWeightSquared_C\=0-10_J\=160-180_T\>2v0.pdf paperFigures/Figure_A26-f.pdf
mv figures/energyEnergyCorrelator_jewelRatio_energyWeightSquared_C\=0-10_J\=180-200_T\>1v0.pdf paperFigures/Figure_A26-g.pdf
mv figures/energyEnergyCorrelator_jewelRatio_energyWeightSquared_C\=0-10_J\=180-200_T\>2v0.pdf paperFigures/Figure_A26-h.pdf

# We can also compare double ratios with JEWEL predictions
root -l -b -q 'plotting/modelComparison.C(1,5,2,true)'
root -l -b -q 'plotting/modelComparison.C(2,5,2,true)'

# Move the relevant plots to the paper figure folder
mv figures/energyEnergyCorrelator_jewelDoubleRatio_nominalEnergyWeight_C\=0-10_J\=120-140.pdf paperFigures/Figure_A31-a.pdf
mv figures/energyEnergyCorrelator_jewelDoubleRatio_energyWeightSquared_C\=0-10_J\=120-140.pdf paperFigures/Figure_A31-b.pdf
mv figures/energyEnergyCorrelator_jewelDoubleRatio_nominalEnergyWeight_C\=0-10_J\=140-160.pdf paperFigures/Figure_A31-c.pdf
mv figures/energyEnergyCorrelator_jewelDoubleRatio_energyWeightSquared_C\=0-10_J\=140-160.pdf paperFigures/Figure_A31-d.pdf
mv figures/energyEnergyCorrelator_jewelDoubleRatio_nominalEnergyWeight_C\=0-10_J\=160-180.pdf paperFigures/Figure_A31-e.pdf
mv figures/energyEnergyCorrelator_jewelDoubleRatio_energyWeightSquared_C\=0-10_J\=160-180.pdf paperFigures/Figure_A31-f.pdf
mv figures/energyEnergyCorrelator_jewelDoubleRatio_nominalEnergyWeight_C\=0-10_J\=180-200.pdf paperFigures/Figure_A31-g.pdf
mv figures/energyEnergyCorrelator_jewelDoubleRatio_energyWeightSquared_C\=0-10_J\=180-200.pdf paperFigures/Figure_A31-h.pdf

# We still need the PbPb distribution and PbPb/pp ratio comparisons with different k-value in the perturbative calculation
root -l -b -q 'plotting/modelComparison.C(1,1,0,true)'
root -l -b -q 'plotting/modelComparison.C(1,1,1,true)'

# Move the relevent figures to the paper figures folder
mv figures/energyEnergyCorrelator_holguinDistribution_nominalEnergyWeight_C\=0-10_J\=120-140_T\>1v0.pdf paperFigures/Figure_A13-a.pdf
mv figures/energyEnergyCorrelator_holguinDistribution_nominalEnergyWeight_C\=0-10_J\=140-160_T\>1v0.pdf paperFigures/Figure_A13-b.pdf
mv figures/energyEnergyCorrelator_holguinDistribution_nominalEnergyWeight_C\=0-10_J\=160-180_T\>1v0.pdf paperFigures/Figure_A13-c.pdf
mv figures/energyEnergyCorrelator_holguinDistribution_nominalEnergyWeight_C\=0-10_J\=180-200_T\>1v0.pdf paperFigures/Figure_A13-d.pdf

mv figures/energyEnergyCorrelator_holguinRatio_nominalEnergyWeight_C\=0-10_J\=120-140_T\>1v0.pdf paperFigures/Figure_A24-a.pdf
mv figures/energyEnergyCorrelator_holguinRatio_nominalEnergyWeight_C\=0-10_J\=140-160_T\>1v0.pdf paperFigures/Figure_A24-b.pdf
mv figures/energyEnergyCorrelator_holguinRatio_nominalEnergyWeight_C\=0-10_J\=160-180_T\>1v0.pdf paperFigures/Figure_A24-c.pdf
mv figures/energyEnergyCorrelator_holguinRatio_nominalEnergyWeight_C\=0-10_J\=180-200_T\>1v0.pdf paperFigures/Figure_A24-d.pdf

# In the end, clean excess figures from the figure folder
rm figures/*.pdf
