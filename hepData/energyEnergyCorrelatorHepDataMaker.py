#
# Macro for creating HepData entry for energy-energy correlator analysis
#
# Code is writted based on examples in https://github.com/HEPData/hepdata_lib/blob/master/examples/
#

# Use hepdata_lib package for easier data handling
import hepdata_lib

# First, create a HepData submission object
# Later all the tables from different figures are added to this object
from hepdata_lib import Submission
submission = Submission()

# Add basic paper info to the submission object
submission.read_abstract("abstract.txt")
#submission.add_link("Webpage with preliminary figures and tables", "https://cms-results.web.cern.ch/cms-results/public-results/preliminary-results/HIN-23-004/index.html")
submission.add_link("Webpage with all figures and tables", "https://cms-results.web.cern.ch/cms-results/public-results/publications/HIN-23-004")
submission.add_link("DOI", "https://doi.org/10.1016/j.physletb.2025.139556")
submission.add_link("arXiv", "https://arxiv.org/abs/2503.19993")
submission.add_record_id(2904406, "inspire")

# All the figures are inserted as table objects to HepData
from hepdata_lib import Table

# Read the root files with the histograms
from hepdata_lib import RootFileReader

#########################################################################################
#  Read the root files containing final results from energy-energy correlator analysis  #
#########################################################################################

# Common naming convention in the final results files
centralityString = ["C0-10", "C10-30", "C30-50", "C50-90", "pp"]
trackPtString = ["T1","T2"]
jetPtString = ["J120-140","J140-160","J160-180","J180-200"]
energyWeightString = ["nominalEnergyWeight", "energyWeightSquared"]
centralityDescription = ["the 0-10 % centrality bin in PbPb", "the 10-30 % centrality bin in PbPb", "the 30-50 % centrality bin in PbPb", "the 50-90 % centrality bin in PbPb", "pp"]
trackPtDescription = ["charged particles with $p_{\\mathrm{T}} > 1$ GeV", "charged particles with $p_{\\mathrm{T}} > 2$ GeV"]
jetPtDescription = ["jet $p_{\\mathrm{T}}$ selection $120 < p_{\\mathrm{T,jet}} < 140$ GeV", "jet $p_{\\mathrm{T}}$ selection $140 < p_{\\mathrm{T,jet}} < 160$ GeV", "jet $p_{\\mathrm{T}}$ selection $160 < p_{\\mathrm{T,jet}} < 180$ GeV", "jet $p_{\\mathrm{T}}$ selection $180 < p_{\\mathrm{T,jet}} < 200$ GeV"]
energyWeightDescription = ["for energy weight $n=1$", "for energy weight $n=2$"]

# Numbers of bins
nCentrality = 5       # Number of centrality bins, including pp
nTrackPt = 2          # Number of track pT cuts for energy-energy correlators
nJetPt = 4            # Number of jet pT selections for energy-energy correlators
nEnergyWeight = 2     # Number of weight exponent for energy-energy correlators

# File for energy-energy correlator distributions. Paper figures 1 and 2.
distributionReader = RootFileReader("hepdata_energyEnergyCorrelatorDistribution_hin-23-004.root")

# Read the energy-energy correlator distributions from the root file
energyEnergyCorrelatorDistribution = []
energyEnergyCorrelatorDistributionCorrelatedUncertainty = []
energyEnergyCorrelatorDistributionUncorrelatedUncertainty = []

for iEnergyWeight in energyWeightString:
    for iTrackPt in trackPtString:
        for iJetPt in jetPtString:
            for iCentrality in centralityString:
                energyEnergyCorrelatorDistribution.append(distributionReader.read_hist_1d("energyEnergyCorrelator_" + iEnergyWeight + "_" + iCentrality + "_" + iJetPt + "_" + iTrackPt))
                energyEnergyCorrelatorDistributionCorrelatedUncertainty.append(distributionReader.read_hist_1d("energyEnergyCorrelator_correlatedUncertainty_" + iEnergyWeight + "_" + iCentrality + "_" + iJetPt + "_" + iTrackPt))
                energyEnergyCorrelatorDistributionUncorrelatedUncertainty.append(distributionReader.read_hist_1d("energyEnergyCorrelator_uncorrelatedUncertainty_" + iEnergyWeight + "_" + iCentrality + "_" + iJetPt + "_" + iTrackPt))

# File for PbPb/pp ratios of energy-energy correlators. For paper figures 3 and 4.
ratioReader = RootFileReader("hepdata_energyEnergyCorrelatorRatio_hin-23-004.root")

# Read the PbPb/pp ratios of energy-energy correlators from the root file
energyEnergyCorrelatorRatio = []
energyEnergyCorrelatorRatioCorrelatedUncertainty = []
energyEnergyCorrelatorRatioUncorrelatedUncertainty = []

for iEnergyWeight in energyWeightString:
    for iTrackPt in trackPtString:
        for iJetPt in jetPtString:
            for iCentrality in centralityString:
                # There is no ratio for pp/pp, so skip the pp centrality string for the ratio
                if iCentrality == "pp":
                    continue
                energyEnergyCorrelatorRatio.append(ratioReader.read_hist_1d("energyEnergyCorrelatorRatio_" + iEnergyWeight + "_" + iCentrality + "_" + iJetPt + "_" + iTrackPt))
                energyEnergyCorrelatorRatioCorrelatedUncertainty.append(ratioReader.read_hist_1d("energyEnergyCorrelatorRatio_correlatedUncertainty_" + iEnergyWeight + "_" + iCentrality + "_" + iJetPt + "_" + iTrackPt))
                energyEnergyCorrelatorRatioUncorrelatedUncertainty.append(ratioReader.read_hist_1d("energyEnergyCorrelatorRatio_uncorrelatedUncertainty_" + iEnergyWeight + "_" + iCentrality + "_" + iJetPt + "_" + iTrackPt))
        
# File for double ratios. Most of the figures are in supplemental material.
doubleRatioReader = RootFileReader("hepdata_energyEnergyCorrelatorDoubleRatio_hin-23-004.root")

# Read the double ratios from the root file
energyEnergyCorrelatorDoubleRatio = []
energyEnergyCorrelatorDoubleRatioUncertainty = []

for iEnergyWeight in energyWeightString:
    for iCentrality in centralityString:
        # There is no ratio for pp/pp, so skip the pp centrality string for the double ratio
        if iCentrality == "pp":
            continue
        for iJetPt in jetPtString:     
            energyEnergyCorrelatorDoubleRatio.append(doubleRatioReader.read_hist_1d("energyEnergyCorrelatorDoubleRatio_" + iEnergyWeight + "_" + iCentrality + "_" + iJetPt))
            energyEnergyCorrelatorDoubleRatioUncertainty.append(doubleRatioReader.read_hist_1d("energyEnergyCorrelatorDoubleRatio_systematicUncertainty_" + iEnergyWeight + "_" + iCentrality + "_" + iJetPt))



# Reaction labels to be added as keywords
pbpbReactionLabel = "PB PB --> JET CHARGED X"
ppReactionLabel = "P P --> JET CHARGED X"
reactionLabel = [pbpbReactionLabel, pbpbReactionLabel, pbpbReactionLabel, pbpbReactionLabel, ppReactionLabel]

# Function for finding the correct rounding information to a given number
#
# Arguments:
#   myNumber = Number which need determination of significant digits
#
# Return: The integer that should be given to the round() function in order to round the number to two significant digits
#
def decimalsToRound(myNumber):

    # Transform the number to a string and tokenize it from decimal point
    tokens = str(myNumber).split(".")

    # If there are more than one number before decimal place, there are enough significant digits before the decimal point that we can cut out everything below
    if len(tokens[0]) > 1:
        return 0

    # If there is exactly one non-zero number before decimal point, we need to include one digit after the decimal point
    if tokens[0] != "0":
        return 1

    # If the number before the decimal point is zero, count the number of leading zeros after the decimal point
    counter = 0
    for i in range(len(tokens[1])):
        if tokens[1][i] == "0":
            counter = counter+1
        else:
            break;

    # We are rounding to two significant digits, so the amount of numbers to keep after decimal point is two plus the first non-zero location
    return 2+counter;

# Read the variables from the histograms
from hepdata_lib import Variable, Uncertainty

# Function for finding x and y-axis variables from input histograms
#
# Arguments:
#  valueHistogram: Histogram containing the actual values and statistical errors
#  correlatedErrorHistogram: Histogram containing the correlated systematic errors
#  uncorrelatedErrorHistogram: Histogram containing the uncorrelated systematic errors
#  yAxisName: Label given to the histogram y-axis
#  includePp: True = Separate histogram exists for pp. False = Histograms are only available for PbPb
#  energyWeightBin: Given energy weight index
#  jetPtBin: Given track pT bin index
#  trackPtBin: Given track pT bin index
#
#  return: Returns HepData variables for x and y-axes from the input histograms
#
def findVariables(valueHistogram, correlatedErrorHistogram, uncorrelatedErrorHistogram, yAxisName, includePp, energyWeightBin, jetPtBin, trackPtBin):

    # Define reaction and centrality labels.
    centralityLabel = ["0-10%", "10-30%", "30-50%", "50-90%", "pp"]
    jetPtLabel = ["$120 < p_{\\mathrm{T,jet}} < 140$ GeV", "$140 < p_{\\mathrm{T,jet}} < 160$ GeV", "$160 < p_{\\mathrm{T,jet}} < 180$ GeV", "$180 < p_{\\mathrm{T,jet}} < 200$ GeV"]
    trackPtLabel = ["$> 1$ GeV", "$> 2$ GeV"]
    energyWeightLabel = ["$n=1$","n=2"]
        
    # Figure out the centrality bin range based on if separate histogram exists for pp
    firstCentralityBin = 0
    lastCentralityBin = 4
    
    if includePp:
        lastCentralityBin = 5
        
    # Find out the binning for DeltaR axis
    xAxis = Variable("$\\Delta r$", is_independent=True, is_binned=True, units="")
    xAxis.values = valueHistogram[0]["x_edges"]

    # Determine how many bins in the histogram is outside of the analysis range of [0.008, 0.4]

    # Find points above 0.4 that need to be removed
    overRange = 0      
    for (i,j) in xAxis.values:
        if j > 0.4:
            overRange = overRange + 1
                
    # Find points below 0.008 that need to be removed
    underRange = 0
    for (i,j) in xAxis.values:
        if j < 0.008:
            underRange = underRange + 1

    # Remove the points that are outside of the analysis range
    for i in range(overRange):
        xAxis.values.pop()
        
    for i in range(underRange):
        xAxis.values.pop(0)
        
    # Define a list for y-axis values
    yAxis = []

    # Read the values
    for iCentrality in range(firstCentralityBin,lastCentralityBin):
        myVariable = Variable(yAxisName, is_independent=False, is_binned=False, units="")
        myVariable.values = valueHistogram[(iCentrality-firstCentralityBin)]["y"]
        myVariable.add_qualifier("SQRT(S)/NUCLEON", "5.02 GeV")
        if includePp:
            myVariable.add_qualifier("RE", reactionLabel[iCentrality])
        else:
            myVariable.add_qualifier("RE", "{}  /  {}".format(pbpbReactionLabel, ppReactionLabel))
        myVariable.add_qualifier("Jet algorithm", "Anti-k$_{\\mathrm{T}}$ R = 0.4")
        myVariable.add_qualifier("Inclusive jet $p_{\\mathrm{T}}$", jetPtLabel[jetPtBin])
        myVariable.add_qualifier("$|\\eta^{\\mathrm{jet}}|$", "< 1.6")
        myVariable.add_qualifier("$p_{\\mathrm{T}}^{\\mathrm{ch}}$",trackPtLabel[trackPtBin])
        myVariable.add_qualifier("Centrality",centralityLabel[iCentrality])

        # If there are points outside of the analysis range, remove them from the values table:
        for j in range(overRange):
            myVariable.values.pop()
                
        for j in range(underRange):
            myVariable.values.pop(0)
        
        # Read the uncertainties
        statisticalUncertainty = Uncertainty("stat", is_symmetric=True)
        statisticalUncertainty.values = valueHistogram[(iCentrality-firstCentralityBin)]["dy"]

        correlatedSystematicUncertainty = Uncertainty("corr sys", is_symmetric=True)
        correlatedSystematicUncertainty.values = correlatedErrorHistogram[(iCentrality-firstCentralityBin)]["dy"]

        uncorrelatedSystematicUncertainty = Uncertainty("uncorr sys", is_symmetric=True)
        uncorrelatedSystematicUncertainty.values = uncorrelatedErrorHistogram[(iCentrality-firstCentralityBin)]["dy"]

        # If there are points outside of the analysis range, remove them from the values table:
        for j in range(overRange):
            statisticalUncertainty.values.pop()
            correlatedSystematicUncertainty.values.pop()
            uncorrelatedSystematicUncertainty.values.pop()
                
        for j in range(underRange):
            statisticalUncertainty.values.pop(0)
            correlatedSystematicUncertainty.values.pop(0)
            uncorrelatedSystematicUncertainty.values.pop(0)

        # Before defining the y-axis, round all numbers to specified precision
        for iBin in range(len(myVariable.values)):

            # Initialize the helper variables for this point
            roundInputValue = 0
            currentDecimals = 0

            # First we need to check which uncertainty is the smallest
            currentDecimals = decimalsToRound(statisticalUncertainty.values[iBin])
            if currentDecimals > roundInputValue:
                roundInputValue = currentDecimals

            currentDecimals = decimalsToRound(correlatedSystematicUncertainty.values[iBin])
            if currentDecimals > roundInputValue:
                roundInputValue = currentDecimals

            currentDecimals = decimalsToRound(uncorrelatedSystematicUncertainty.values[iBin])
            if currentDecimals > roundInputValue:
                roundInputValue = currentDecimals

            # Now that we have the desired precision determined, we can round the numbers accordingly
            myVariable.values[iBin] = round(myVariable.values[iBin], roundInputValue)
            statisticalUncertainty.values[iBin] = round(statisticalUncertainty.values[iBin], roundInputValue)
            correlatedSystematicUncertainty.values[iBin] = round(correlatedSystematicUncertainty.values[iBin], roundInputValue)
            uncorrelatedSystematicUncertainty.values[iBin] = round(uncorrelatedSystematicUncertainty.values[iBin], roundInputValue)
            

        # Add the rounded variables to the y-axis list
        yAxis.append(myVariable)
        yAxis[(iCentrality-firstCentralityBin)].add_uncertainty(statisticalUncertainty)
        yAxis[(iCentrality-firstCentralityBin)].add_uncertainty(correlatedSystematicUncertainty)
        yAxis[(iCentrality-firstCentralityBin)].add_uncertainty(uncorrelatedSystematicUncertainty)
        
    return xAxis, yAxis

# Function for finding x and y-axis variables from double ratio histograms
#
# Arguments:
#  valueHistogram: Histogram containing the actual values and statistical errors
#  errorHistogram: Histogram containing the systematic errors
#  energyWeightBin: Given energy weight index
#  centralityBin: Given centrality bin index
#
#  return: Returns HepData variables for x and y-axes from the input histograms
#
def findVariablesDoubleRatio(valueHistogram, errorHistogram, energyWeightBin, centralityBin):

    # Define labels for variables
    centralityLabel = ["0-10%", "10-30%", "30-50%", "50-90%"]
    jetPtLabel = ["$120 < p_{\\mathrm{T,jet}} < 140$ GeV", "$140 < p_{\\mathrm{T,jet}} < 160$ GeV", "$160 < p_{\\mathrm{T,jet}} < 180$ GeV", "$180 < p_{\\mathrm{T,jet}} < 200$ GeV"]
    energyWeightLabel = ["1","2"]
        
    # Find out the binning for DeltaR axis
    xAxis = Variable("$\\Delta r$", is_independent=True, is_binned=True, units="")
    xAxis.values = valueHistogram[0]["x_edges"]

    # Determine how many bins in the histogram is outside of the analysis range of [0.008, 0.4]

    # Find points above 0.4 that need to be removed
    overRange = 0      
    for (i,j) in xAxis.values:
        if j > 0.4:
            overRange = overRange + 1
                
    # Find points below 0.008 that need to be removed
    underRange = 0
    for (i,j) in xAxis.values:
        if j < 0.008:
            underRange = underRange + 1

    # Remove the points that are outside of the analysis range
    for i in range(overRange):
        xAxis.values.pop()
        
    for i in range(underRange):
        xAxis.values.pop(0)
        
    # Define a list for y-axis values
    yAxis = []

    # Loop is over all jet pT bins and read the variables
    for iJetPt in range(0, len(jetPtLabel)):
        myVariable = Variable("$\\frac{\\mathrm{PbPb/pp} (p_{\\mathrm{T}}^{\\mathrm{ch}} > 2 \\mathrm{GeV})}{\\mathrm{PbPb/pp} (p_{\\mathrm{T}}^{\\mathrm{ch}} > 1 \\mathrm{GeV})}$", is_independent=False, is_binned=False, units="")
        myVariable.values = valueHistogram[iJetPt]["y"]
        myVariable.add_qualifier("SQRT(S)/NUCLEON", "5.02 GeV")
        myVariable.add_qualifier("RE", "{}  /  {}".format(pbpbReactionLabel, ppReactionLabel))
        myVariable.add_qualifier("Jet algorithm", "Anti-k$_{\\mathrm{T}}$ R = 0.4")
        myVariable.add_qualifier("$|\\eta^{\\mathrm{jet}}|$", "< 1.6")
        myVariable.add_qualifier("Centrality",centralityLabel[centralityBin])
        myVariable.add_qualifier("$n$",energyWeightLabel[energyWeightBin])
        myVariable.add_qualifier("Inclusive jet $p_{\\mathrm{T}}$", jetPtLabel[iJetPt])

        # If there are points outside of the analysis range, remove them from the values table:
        for j in range(overRange):
            myVariable.values.pop()
                
        for j in range(underRange):
            myVariable.values.pop(0)
        
        statisticalUncertainty = Uncertainty("stat", is_symmetric=True)
        statisticalUncertainty.values = valueHistogram[iJetPt]["dy"]

        systematicUncertainty = Uncertainty("sys", is_symmetric=True)
        systematicUncertainty.values = errorHistogram[iJetPt]["dy"]

        # If there are points outside of the analysis range, remove them from the values table:
        for j in range(overRange):
            statisticalUncertainty.values.pop()
            systematicUncertainty.values.pop()
                
        for j in range(underRange):
            statisticalUncertainty.values.pop(0)
            systematicUncertainty.values.pop(0)

         # Before defining the y-axis, round all numbers to specified precision
        for iBin in range(len(myVariable.values)):

            # Initialize the helper variables for this point
            roundInputValue = 0
            currentDecimals = 0

            # First we need to check which uncertainty is the smallest
            currentDecimals = decimalsToRound(statisticalUncertainty.values[iBin])
            if currentDecimals > roundInputValue:
                roundInputValue = currentDecimals

            currentDecimals = decimalsToRound(systematicUncertainty.values[iBin])
            if currentDecimals > roundInputValue:
                roundInputValue = currentDecimals

            # Now that we have the desired precision determined, we can round the numbers accordingly
            myVariable.values[iBin] = round(myVariable.values[iBin], roundInputValue)
            statisticalUncertainty.values[iBin] = round(statisticalUncertainty.values[iBin], roundInputValue)
            systematicUncertainty.values[iBin] = round(systematicUncertainty.values[iBin], roundInputValue)

        # Add the rounded variables to the y-axis list
        yAxis.append(myVariable)
        yAxis[iJetPt].add_uncertainty(statisticalUncertainty)
        yAxis[iJetPt].add_uncertainty(systematicUncertainty)
        
    return xAxis, yAxis

# Different jet pT values are located in different columns in the figures. This needs to be added to the description of the tables.
jetPtLocation = ["leftmost column", "second column from the left", "second column from the right", "rightmost column"]
energyWeightLocation = ["left column", "right column"]

###############################################################################################
#     Table for energy-energy correlator distributions with pT cut 1 GeV. Paper figure 1.     #
###############################################################################################

# Make separate tables from each centrality and track pT bin
for iEnergyWeight in range(0, nEnergyWeight):
    for iJetPt in range(0, nJetPt):
        table1 = Table("Figure 2-{:d}".format(1+iJetPt+nJetPt*iEnergyWeight))
        table1.description = "The energy-energy correlator distributions constructed with {:s} {:s} and {:s}. The results are shown for different centrality bins in PbPb collisions and for pp collisions.".format(trackPtDescription[0], energyWeightDescription[iEnergyWeight], jetPtDescription[iJetPt])
        table1.location = "Data from figure 2, {:s}.".format(jetPtLocation[iJetPt])
        table1.keywords["observables"] = ["Energy-energy correlator"]
        table1.keywords["reactions"] = [pbpbReactionLabel, ppReactionLabel]
        table1.keywords["phrases"] = ["CMS", "jet", "energy", "correlator", "EEC", "E2C", "PbPb", "pp", "QGP"]
        #table1.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

        # Extract x- and y-axis information from the histograms
        xDeltaR, yDeltaR = findVariables(energyEnergyCorrelatorDistribution[nCentrality*iJetPt + nCentrality*nJetPt*nTrackPt*iEnergyWeight : nCentrality + nCentrality*iJetPt + nCentrality*nJetPt*nTrackPt*iEnergyWeight], energyEnergyCorrelatorDistributionCorrelatedUncertainty[nCentrality*iJetPt + nCentrality*nJetPt*nTrackPt*iEnergyWeight : nCentrality + nCentrality*iJetPt + nCentrality*nJetPt*nTrackPt*iEnergyWeight], energyEnergyCorrelatorDistributionUncorrelatedUncertainty[nCentrality*iJetPt + nCentrality*nJetPt*nTrackPt*iEnergyWeight : nCentrality + nCentrality*iJetPt + nCentrality*nJetPt*nTrackPt*iEnergyWeight], "EEC", True, iEnergyWeight, iJetPt, 0)

        # Add the variables to the table
        table1.add_variable(xDeltaR)
        for variable in yDeltaR:
            table1.add_variable(variable)

        # Add the table to submission object
        submission.add_table(table1)

###############################################################################################
#     Table for energy-energy correlator distributions with pT cut 2 GeV. Paper figure 2.     #
###############################################################################################

# Make separate tables from each centrality and track pT bin
for iEnergyWeight in range(0, nEnergyWeight):
    for iJetPt in range(0, nJetPt):
        table2 = Table("Figure 3-{:d}".format(1+iJetPt+nJetPt*iEnergyWeight))
        table2.description = "The energy-energy correlator distributions constructed with {:s} {:s} and {:s}. The results are shown for different centrality bins in PbPb collisions and for pp collisions.".format(trackPtDescription[1], energyWeightDescription[iEnergyWeight], jetPtDescription[iJetPt])
        table2.location = "Data from figure 3, {:s}.".format(jetPtLocation[iJetPt])
        table2.keywords["observables"] = ["Energy-energy correlator"]
        table2.keywords["reactions"] = [pbpbReactionLabel, ppReactionLabel]
        table2.keywords["phrases"] = ["CMS", "jet", "energy", "correlator", "EEC", "E2C", "PbPb", "pp", "QGP"]
        #table2.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

        # Extract x- and y-axis information from the histograms
        xDeltaR, yDeltaR = findVariables(energyEnergyCorrelatorDistribution[nCentrality*iJetPt + nCentrality*nJetPt + nCentrality*nJetPt*nTrackPt*iEnergyWeight : nCentrality + nCentrality*iJetPt + nCentrality*nJetPt  + nCentrality*nJetPt*nTrackPt*iEnergyWeight], energyEnergyCorrelatorDistributionCorrelatedUncertainty[nCentrality*iJetPt + nCentrality*nJetPt  + nCentrality*nJetPt*nTrackPt*iEnergyWeight : nCentrality + nCentrality*iJetPt + nCentrality*nJetPt  + nCentrality*nJetPt*nTrackPt*iEnergyWeight], energyEnergyCorrelatorDistributionUncorrelatedUncertainty[nCentrality*iJetPt + nCentrality*nJetPt  + nCentrality*nJetPt*nTrackPt*iEnergyWeight : nCentrality + nCentrality*iJetPt + nCentrality*nJetPt  + nCentrality*nJetPt*nTrackPt*iEnergyWeight], "EEC", True, iEnergyWeight, iJetPt, 1)

        # Add the variables to the table
        table2.add_variable(xDeltaR)
        for variable in yDeltaR:
            table2.add_variable(variable)

        # Add the table to submission object
        submission.add_table(table2)

# For the ratios, there is no separate pp distribution. This means that there will be only 4 centrality bins.
nCentrality = 4

# #############################################################################################
# #    Table for energy-energy correlator PbPb/pp ratios with pT cut 1 GeV. Paper figure 3.   #
# #############################################################################################

# Make separate tables from each centrality and track pT bin
for iEnergyWeight in range(0, nEnergyWeight):
    for iJetPt in range(0, nJetPt):
        table3 = Table("Figure 4-{:d}".format(1+iJetPt+nJetPt*iEnergyWeight))
        table3.description = "The energy-energy correlator PbPb/pp ratios constructed with {:s} {:s} and {:s}. The results are shown for different centrality bins in PbPb collisions.".format(trackPtDescription[0], energyWeightDescription[iEnergyWeight], jetPtDescription[iJetPt])
        table3.location = "Data from figure 4, {:s}.".format(jetPtLocation[iJetPt])
        table3.keywords["observables"] = ["Energy-energy correlator PbPb to pp ratio"]
        table3.keywords["reactions"] = [pbpbReactionLabel, ppReactionLabel]
        table3.keywords["phrases"] = ["CMS", "jet", "energy", "correlator", "EEC", "E2C", "PbPb", "pp", "QGP"]
        #table3.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

        # Extract x- and y-axis information from the histograms
        xDeltaR, yDeltaR = findVariables(energyEnergyCorrelatorRatio[nCentrality*iJetPt + nCentrality*nJetPt*nTrackPt*iEnergyWeight : nCentrality + nCentrality*iJetPt + nCentrality*nJetPt*nTrackPt*iEnergyWeight], energyEnergyCorrelatorRatioCorrelatedUncertainty[nCentrality*iJetPt + nCentrality*nJetPt*nTrackPt*iEnergyWeight : nCentrality + nCentrality*iJetPt + nCentrality*nJetPt*nTrackPt*iEnergyWeight], energyEnergyCorrelatorRatioUncorrelatedUncertainty[nCentrality*iJetPt + nCentrality*nJetPt*nTrackPt*iEnergyWeight : nCentrality + nCentrality*iJetPt + nCentrality*nJetPt*nTrackPt*iEnergyWeight], "PbPb/pp", False, iEnergyWeight, iJetPt, 0)

        # Add the variables to the table
        table3.add_variable(xDeltaR)
        for variable in yDeltaR:
            table3.add_variable(variable)

        # Add the table to submission object
        submission.add_table(table3)

# #############################################################################################
# #    Table for energy-energy correlator PbPb/pp ratios with pT cut 2 GeV. Paper figure 4.   #
# #############################################################################################

# Make separate tables from each centrality and track pT bin
for iEnergyWeight in range(0, nEnergyWeight):
    for iJetPt in range(0, nJetPt):
        table4 = Table("Figure 5-{:d}".format(1+iJetPt+nJetPt*iEnergyWeight))
        table4.description = "The energy-energy correlator PbPb/pp ratios constructed with {:s} {:s} and {:s}. The results are shown for different centrality bins in PbPb collisions.".format(trackPtDescription[1], energyWeightDescription[iEnergyWeight], jetPtDescription[iJetPt])
        table4.location = "Data from figure 5, {:s}.".format(jetPtLocation[iJetPt])
        table4.keywords["observables"] = ["Energy-energy correlator PbPb to pp ratio"]
        table4.keywords["reactions"] = [pbpbReactionLabel, ppReactionLabel]
        table4.keywords["phrases"] = ["CMS", "jet", "energy", "correlator", "EEC", "E2C", "PbPb", "pp", "QGP"]
        #table4.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

        # Extract x- and y-axis information from the histograms
        xDeltaR, yDeltaR = findVariables(energyEnergyCorrelatorRatio[nCentrality*iJetPt + nCentrality*nJetPt + nCentrality*nJetPt*nTrackPt*iEnergyWeight : nCentrality + nCentrality*iJetPt + nCentrality*nJetPt  + nCentrality*nJetPt*nTrackPt*iEnergyWeight], energyEnergyCorrelatorRatioCorrelatedUncertainty[nCentrality*iJetPt + nCentrality*nJetPt  + nCentrality*nJetPt*nTrackPt*iEnergyWeight : nCentrality + nCentrality*iJetPt + nCentrality*nJetPt  + nCentrality*nJetPt*nTrackPt*iEnergyWeight], energyEnergyCorrelatorRatioUncorrelatedUncertainty[nCentrality*iJetPt + nCentrality*nJetPt  + nCentrality*nJetPt*nTrackPt*iEnergyWeight : nCentrality + nCentrality*iJetPt + nCentrality*nJetPt  + nCentrality*nJetPt*nTrackPt*iEnergyWeight], "PbPb/pp", False, iEnergyWeight, iJetPt, 1)

        # Add the variables to the table
        table4.add_variable(xDeltaR)
        for variable in yDeltaR:
            table4.add_variable(variable)

        # Add the table to submission object
        submission.add_table(table4)

# ##############################################################################################
# #   Table for leading energy-energy correlator double ratios. Figures are in the appendix.   #
# ##############################################################################################

# Make separate tables for each jet pT bin
for iCentrality in range(nCentrality-1, -1, -1):
    for iEnergyWeight in range(0, nEnergyWeight):   
        table5 = Table("Figure A{:d}-{:d}".format(30 - iCentrality, iEnergyWeight))
        table5.description = "The double ratios of PbPb to pp single ratios with $p_{\\mathrm{T}}^{\\mathrm{ch}}$ > 2 GeV and $p_{\\mathrm{T}}^{\\mathrm{ch}}$ > 1 GeV."
        table5.location = "Data from figure A{:d}, {:s}.".format(30 - iCentrality, energyWeightLocation[iEnergyWeight])
        table5.keywords["observables"] = ["Energy-energy correlator double ratio"]
        table5.keywords["reactions"] = [pbpbReactionLabel, ppReactionLabel]
        table5.keywords["phrases"] = ["CMS", "jet", "energy", "correlator", "EEC", "E2C", "PbPb", "pp", "QGP"]
        # #table5.add_image("example_inputs/CMS-B2G-17-009_Figure_004-a.pdf") # Possibility to add image

        # Extract x- and y-axis information from the histograms
        xDeltaR, yDeltaR = findVariablesDoubleRatio(energyEnergyCorrelatorDoubleRatio[nJetPt*iCentrality + nJetPt*nCentrality*iEnergyWeight : nJetPt + nJetPt*iCentrality + nJetPt*nCentrality*iEnergyWeight], energyEnergyCorrelatorDoubleRatioUncertainty[nJetPt*iCentrality + nJetPt*nCentrality*iEnergyWeight : nJetPt + nJetPt*iCentrality + nJetPt*nCentrality*iEnergyWeight], iEnergyWeight, iCentrality)

        # Add the variables to the table
        table5.add_variable(xDeltaR)
        for variable in yDeltaR:
            table5.add_variable(variable)

        # Add the table to submission object
        submission.add_table(table5)


###########################################################################
#                      Finalize the submission                            #
###########################################################################


# Add common keywords for all tables
for table in submission.tables:
    table.keywords["cmenergies"] = [5020]   # Center-of-mass energy

# Create the submission file for upload
outputDirectory = "energyEnergyCorrelatorHepData"
submission.create_files(outputDirectory)
