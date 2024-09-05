#ifndef SPLITCANVAS_H
#define SPLITCANVAS_H

// Root includes
#include "TCanvas.h"
#include "TPad.h"
#include "TString.h"

/*
 * SplitCanvas class
 *
 * Class whose purpose is to provide easy interface to create large split canvases
 */
class SplitCanvas: public TCanvas {

public :

  SplitCanvas();            // Default constructor
  SplitCanvas(TString canvasName, TString canvasTitle, double width, double height); // Custom constructor
  ~SplitCanvas() = default; // Default destructor
	void DivideNeatly(int nDivisionsHorizontal, int nDivisionsVertical);
  void DivideNeatly2x2();
	void CD(int iPad);
	double GetBottomRowScale();  // Getter for bottom row title size scaling factor
	double GetBottomPadMargin(); // Getter for the margin added to the bottom pad
	double GetLeftColumnScale(); // Getter for the left column horizontal scaling factor
	double GetLeftPadMargin();   // Getter for the margin added to the left pad

private :

  // Class variables	
	TString fCanvasName;     // Name given to the SplitCanvas
	TPad **fPad;             // Pad controlling to which pad within the split canvas we are drawing
	double fCanvasWidth;     // Width of the created canvas
	double fCanvasHeight;    // Height of the created canvas
	double fBottomRowScale;  // Scaling factor to be applied to the bottom row title sizes
	double fBottomPadMargin; // Get the extra margin added to the bottom pad
	double fLeftColumnScale; // Scaling factor to be applied to the leftmost column
	double fLeftPadMargin;   // Extra margin added to the left pad
};

#endif
