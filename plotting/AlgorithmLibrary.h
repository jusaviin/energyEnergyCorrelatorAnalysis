#ifndef ALGORITHMLIBRARY_H
#define ALGORITHMLIBRARY_H

/*
 * This class is a collection of methods that are used to process the
 * results produced by the dijet analysis
 */

// C++ includes
#include <iostream>
#include <fstream>    // File stream for input/output to/from files
#include <assert.h>   // Standard c++ debugging tool. Terminates the program if expression given evaluates to 0.
#include <tuple>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>

// Root includes
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TF2.h>
#include <TMath.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TObjArray.h>
#include <TObjString.h>

class AlgorithmLibrary{
  
public:
    
  AlgorithmLibrary();   // Constructor
  AlgorithmLibrary(const AlgorithmLibrary& in); // Copy constructor
  ~AlgorithmLibrary();  // Destructor
  AlgorithmLibrary& operator=(const AlgorithmLibrary& in); // Equal sign operator
  
  std::pair<double,double> FindHistogramMinMax(TH1D* histogram, std::pair<double,double> currentMinMax); // Find the minimum and maximum values from a histogram. They must be more extreme than the current values, also given as input
  std::pair<double,double> FindHistogramMinMax(TH1D* histogram, std::pair<double,double> currentMinMax, std::pair<double,double> searchRange); // Find the minimum and maximum values from a histogram. They must be more extreme than the current values, also given as input
  std::pair<double,double> FindHistogramMinMax(TH1D* histogram, std::pair<double,double> currentMinMax, std::pair<int,int> searchRange); // Find the minimum and maximum values from a histogram. They must be more extreme than the current values, also given as input
  TH1D* RebinAsymmetric(TH1D* histogramInNeedOfRebinning, const int nBins, const double* binEdges); // Asymmetric rebinning for one-dimensional histograms
  TH2D* RebinHistogram(TH2D* histogramInNeedOfRebinning, const int nBinsX, const double* binBordersX, const int nBinsY, const double* binBordersY, const bool undoBinArea, const bool normalizeBinArea); // Rebin a two-dimensional histogram with given bin borders
  void NormalizeMatrix(TH2D* histogramInNeedOfNormalization, const double value = 1, const int direction = 1);  // Normalize rows or columns of a 2D histogram to a given value
  void NormalizeColumns(TH2D* histogramInNeedOfNormalization, const double value = 1);  // Normalize all the columns of a 2-D histogram to a given value
  void NormalizeRows(TH2D* histogramInNeedOfNormalization, const double value = 1);  // Normalize all the rows of a 2-D histogram to a given value
  TH2D* RotateHistogram(TH2D* originalHistogram); // Rotate two dimensional histogram 90 degrees
  void RemoveUncertainties(TH1* histogramInNeedForAccuracy, const double epsilon = 1e-10); // Remove uncertainties from a histogram
  void TransformToRelativeUncertainty(TH1D* transformedHistogram, const bool centerAtOne = false);
  void TransformToAbsoluteUncertainty(TH1D* transformedHistogram, TH1D* absoluteScaleHistogram, const bool centerAtOne = false);
  void SuppressSingleBinFluctuations(TH1D* fluctuatingHistogram, const double lowRange, const double highRange, const double threshold, const double suppressionLevel);
  void SquareHistogram(TH1* transformedHistogram); // Square the contents of a histogram
  void ExponentiateHistogram(TH1* transformedHistogram, const double exponent); // Exponentiate the contents of a histogram
  void ChangeBinWidthToBinAreaNormalization(TH1* normalizedHistogram); // Change the normalization of the histogram
  TH1D* GetCumulant(TH1* transformedHistogram, const int lowestBin = 1); // Transform histogram to cumulant
  TString GetToday(); // Getter for today's date
  std::vector<std::pair<TString, TGraph*>>  GetGraphsFromDatFile(TString fileName);
  TH1D* Histogrammify(TGraph* gSource, TH1D* hBinning);
  TH1D* HistogrammifyWithErrors(TGraphErrors* gSource, TH1D* hBinning);
  
private:
  
  bool CheckBinBoundaries(const int nCheckedBins, const double* checkedBins, const TAxis* originalAxis); // Checker that new bin boundaries are compatible with the old ones
  
};

#endif
