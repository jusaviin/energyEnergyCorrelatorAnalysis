#include "SplitCanvas.h"

/*
 * Implementation of the SplitCanvas class
 */

// Default constructor
SplitCanvas::SplitCanvas() :
  TCanvas("defaultName", "defaultTitle", 800, 600),
  fCanvasName("defaultName"),
  fCanvasWidth(800),
  fCanvasHeight(600),
  fBottomRowScale(1),
  fBottomPadMargin(0),
  fLeftPadMargin(0)
{
  fPad = NULL;
}

// Custom constructor
SplitCanvas::SplitCanvas(TString canvasName, TString canvasTitle, double width, double height):
  TCanvas(canvasName, canvasTitle, width, height),
  fCanvasName(canvasName),
  fCanvasWidth(width),
  fCanvasHeight(height),
  fBottomRowScale(1),
  fBottomPadMargin(0),
  fLeftPadMargin(0)
{
  fPad = NULL;
}

/* 
 * Divide the canvas to a defined number of divisions in horizontal (x) and vertical (y) directions
 *
 *  int nDivisionsHorizontal = Number of divisions in x-direction
 *  int nDivisionsVertical = Number of divisions in y-direction
 */
void SplitCanvas::DivideNeatly(int nDivisionsHorizontal, int nDivisionsVertical) {
  //if(nDivisionsHorizontal == 2 && nDivisionsVertical == 2){
  //  DivideNeatly2x2();
  //  return;
  //}

  // Remember the original margin, and then set the margin to zero
  double marginLeft = this->GetLeftMargin();
  double marginRight = this->GetRightMargin();
  double marginTop = this->GetTopMargin();
  double marginBottom = this->GetBottomMargin();
  this->SetMargin(0, 0, 0, 0);

  // Calculate the height and width of a single pad within the split canvas without margins
  double padHeightWithoutMargin = (1 - marginTop - marginBottom) / nDivisionsHorizontal;
  double padWidthWithoutMargin = (1 - marginLeft - marginRight) / nDivisionsVertical;

  // Calculate the absolute size of left and bottom margins
  double absoluteLeftMargin = fCanvasWidth * marginLeft;
  double absoluteBottomMargin = fCanvasHeight * marginBottom;

  // Create pads in correct positions with correct sizes
  double leftShift, bottomShift;
  fPad = new TPad*[nDivisionsHorizontal * nDivisionsVertical];
  for(int iRow = 0; iRow < nDivisionsHorizontal; iRow++) {
    for(int iColumn = 0; iColumn < nDivisionsVertical; iColumn++) {
      leftShift = 0;  bottomShift = 0;

      // If we are in the first column, add left margin to these pads in order to draw the y-axis
      if(iColumn == 0) leftShift = marginLeft;

      // If we are in the last row, add bottom margin to these pads in order to draw the x-axis
      if(iRow == nDivisionsHorizontal - 1) bottomShift = marginBottom;

      // Create the pad into correct position taking into account that left and bottom margins are only in the pads at the edge of the canvas
      fPad[iRow * nDivisionsVertical + iColumn] = new TPad(fCanvasName, "", padWidthWithoutMargin * iColumn + marginLeft - leftShift, padHeightWithoutMargin * (nDivisionsHorizontal - iRow - 1) + marginBottom - bottomShift, marginLeft + padWidthWithoutMargin * (iColumn + 1), padHeightWithoutMargin * (nDivisionsHorizontal - iRow) + marginBottom);

      // Set the size of the margin to be added for the first column and last row
      if(iColumn == 0) {
      	leftShift = absoluteLeftMargin / (absoluteLeftMargin + padWidthWithoutMargin * fCanvasWidth);
      	fLeftPadMargin = leftShift;
      }
      if(iRow == nDivisionsHorizontal - 1){
      	bottomShift = absoluteBottomMargin / (absoluteBottomMargin + padHeightWithoutMargin * fCanvasHeight);
      	fBottomRowScale = 1 - bottomShift;
        fBottomPadMargin = bottomShift;
      } 

      // Assign margins for the pads
      fPad[iRow * nDivisionsVertical + iColumn]->SetMargin(leftShift, 0, bottomShift, 0);
      fPad[iRow * nDivisionsVertical + iColumn]->Draw();
    }
  }
}

/*
 * Special settings if the division is done 2x2
 */
void SplitCanvas::DivideNeatly2x2(){

	// Do two division in both x- and y-directions
  int nDivisionsHorizontal = 2;
  int nDivisionsVertical = 2;

  // Remember the original margin, and then set the margin to zero
  double marginLeft = this->GetLeftMargin();
  double marginRight = this->GetRightMargin();
  double marginTop = this->GetTopMargin();
  double marginBottom = this->GetBottomMargin();
  this->SetMargin(0,0,0,0);

  // Calculate the height and width of a single pad within the split canvas without margins
  double padHeightWithoutMargin = (1-marginTop-marginBottom)/nDivisionsHorizontal;
  double padWidthWithoutMargin = (1-marginLeft-marginRight)/nDivisionsVertical;

  // Calculate the absolute size of left and bottom margins
  double absoluteLeftMargin = fCanvasWidth*marginLeft;
  double absoluteBottomMargin = fCanvasHeight*marginBottom;

  // Create pads in correct positions with correct sizes
  double leftShift, bottomShift;
  double extraMarginLeft, extraMarginRight, extraMarginBottom, extraMarginTop;
  fPad = new TPad*[nDivisionsHorizontal*nDivisionsVertical];
  for(int iRow =0 ; iRow < nDivisionsHorizontal; iRow++){
    for(int iColumn = 0; iColumn < nDivisionsVertical; iColumn++){
      leftShift=0; bottomShift =0;
      extraMarginLeft = 0; extraMarginRight = 0; extraMarginBottom = 0; extraMarginTop = 0;

      // If we are in the first column, add left margin to these pads in order to draw the y-axis
      if(iColumn == 0) leftShift = marginLeft;

      // If we are in the last row, add bottom margin to these pads in order to draw the x-axis
      if(iRow==nDivisionsHorizontal-1) bottomShift = marginBottom;

      // Add extra bottom margin for the top left pad
      if(iRow == 0 && iColumn == 0) extraMarginBottom = 0.1;

      // Add extra left margin for the bottom right pad
      if(iRow == 1 && iColumn == 1) extraMarginLeft = 0.1;

      // Add extra top and right margin for the bottom left pad
      if(iRow == 1 && iColumn == 0) {extraMarginTop = 0.1; extraMarginRight = 0.1;}

      // Create the pad into correct position taking into account extra margins and that left and bottom margins are only in the pads at the edge of the canvas
      fPad[iRow*nDivisionsVertical+iColumn] = new TPad(fCanvasName, "", padWidthWithoutMargin*iColumn+marginLeft-leftShift-extraMarginLeft, padHeightWithoutMargin*(nDivisionsHorizontal-iRow-1)+marginBottom-bottomShift-extraMarginBottom, marginLeft+padWidthWithoutMargin*(iColumn+1)-extraMarginRight, padHeightWithoutMargin*(nDivisionsHorizontal-iRow)+marginBottom-extraMarginTop);

      // Set the size of the margin to be added for the first column and last row
      if(iColumn ==0 ) leftShift = absoluteLeftMargin/(absoluteLeftMargin+padWidthWithoutMargin*fCanvasWidth);
      if(iRow == nDivisionsHorizontal-1) bottomShift = absoluteBottomMargin/(absoluteBottomMargin+padHeightWithoutMargin*fCanvasHeight);

      // Assign margins for the pads
      fPad[iRow*nDivisionsVertical+iColumn]->SetMargin(leftShift+extraMarginLeft*1.83,0,bottomShift+extraMarginBottom*1.96,0);
      fPad[iRow*nDivisionsVertical+iColumn]->Draw();
    }
  }
}

// Change to fPad index iPad within the canvas
void SplitCanvas::CD(int iPad){
		fPad[iPad]->cd();
}

// Getter for bottom row title size scaling factor
double SplitCanvas::GetBottomRowScale(){
	return fBottomRowScale;
}

// Getter for the margin added to the bottom pad
double SplitCanvas::GetBottomPadMargin(){
  return fBottomPadMargin;
}

// Getter for the margin added to the left pad
double SplitCanvas::GetLeftPadMargin(){
  return fLeftPadMargin;
} 