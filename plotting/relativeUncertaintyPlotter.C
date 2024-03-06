#include "SystematicUncertaintyOrganizer.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

/*
 * Transform an absolute uncertainty histogram to a relative uncertainty histogram
 */ 
void transformToRelativeUncertainty(TH1D* transformedHistogram){
  
  double binContent;
  double binError;
  double relativeError;

  for(int iBin = 1; iBin <= transformedHistogram->GetNbinsX(); iBin++){
    binContent = transformedHistogram->GetBinContent(iBin);
    binError = transformedHistogram->GetBinError(iBin);
    relativeError = binError / binContent;
    transformedHistogram->SetBinContent(iBin, relativeError);
    transformedHistogram->SetBinError(iBin, 0);
  }

}

/*
 * Macro for making final result plots comparing energy-energy correlators between pp and PbPb
 */
void relativeUncertaintyPlotter(){

  enum enumDataType{kPbPb, kPp, kNDataTypes};

  // ============= //
  // Configuration //
  // ============= //
  
  // Input files
  TString uncertaintyFileName[kNDataTypes];
  uncertaintyFileName[kPbPb] = "systematicUncertainties/systematicUncertainties_PbPb_nominalEnergyWeight_includeMCnonClosure_2024-03-06.root";
  // systematicUncertainties_PbPb_nominalEnergyWeight_includeMCnonClosure_2024-03-06.root
  // systematicUncertainties_PbPb_energyWeightSquared_includeMCnonClosure_2024-03-06.root
  // systematicUncertainties_PbPb_energyWeightSquared_includeMCnonClosure_2024-02-23.root
  // systematicUncertainties_PbPb_nominalEnergyWeight_includeMCnonClosure_2024-02-23.root
  uncertaintyFileName[kPp] = "systematicUncertainties/systematicUncertainties_pp_nominalEnergyWeight_includeMCnonClosure_2024-01-29.root";
  // systematicUncertainties_pp_nominalEnergyWeight_includeMCnonClosure_2024-01-29.root
  // systematicUncertainties_pp_energyWeightSquared_includeMCnonClosure_2024-01-29.root
  TFile* uncertaintyFile[kNDataTypes];
  EECCard* uncertaintyCard[kNDataTypes];
  for(int iFile = 0; iFile < kNDataTypes; iFile++){

    // Load the uncertainty file
    uncertaintyFile[iFile] = TFile::Open(uncertaintyFileName[iFile]);

    // Check that the uncertianty file exists
    if(uncertaintyFile[iFile] == NULL){
      cout << "Error! The file " << uncertaintyFileName[iFile].Data() << " does not exist!" << endl;
      cout << "Maybe you forgot the systematicUncertainties/ folder path?" << endl;
      cout << "Will not execute the code" << endl;
      return;
    }
    
    // Load the card from the file and read the collision system
    uncertaintyCard[iFile] = new EECCard(uncertaintyFile[iFile]);
  }

  // ==================================================================
  // ========================= Configuration ==========================
  // ==================================================================
  
  // Find the number of bins from the PbPb card
  const int nCentralityBins = uncertaintyCard[kPbPb]->GetNCentralityBins();
  const int nJetPtBinsEEC = uncertaintyCard[kPbPb]->GetNJetPtBinsEEC();
  const int nTrackPtBinsEEC = uncertaintyCard[kPbPb]->GetNTrackPtBinsEEC();
  
  // The final results are available for all the bins that are unfolded
  int firstDrawnCentralityBin = uncertaintyCard[kPbPb]->GetFirstUnfoldedCentralityBin();
  int lastDrawnCentralityBin = uncertaintyCard[kPbPb]->GetLastUnfoldedCentralityBin();
  
  int firstDrawnJetPtBinEEC = uncertaintyCard[kPbPb]->GetFirstUnfoldedJetPtBin();
  int lastDrawnJetPtBinEEC = uncertaintyCard[kPbPb]->GetLastUnfoldedJetPtBin();
  
  int firstDrawnTrackPtBinEEC = uncertaintyCard[kPbPb]->GetFirstUnfoldedTrackPtBin();
  int lastDrawnTrackPtBinEEC = uncertaintyCard[kPbPb]->GetLastUnfoldedTrackPtBin();

  //firstDrawnTrackPtBinEEC = 3;
  //lastDrawnTrackPtBinEEC = 5;
  
  // Save the plots
  const bool saveFigures = true;
  TString saveComment = "_optimizedUnfoldingBins";

  // Add a name describing the energy weight in the files
  if(uncertaintyCard[kPbPb]->GetWeightExponent() == 2){
    saveComment.Prepend("_energyWeightSquared");
  } else {
    saveComment.Prepend("_nominalEnergyWeight");
  }

  // Zoom settings
  std::pair<double, double> analysisDeltaR = std::make_pair(0.008, 0.39); // DeltaR span in which the analysis is done
  std::pair<double, double> relativeZoom[nCentralityBins+1][nTrackPtBinsEEC];

  for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
    relativeZoom[nCentralityBins][iTrackPt] = std::make_pair(0, 0.06); // Y-axis zoom for pp
    relativeZoom[0][iTrackPt] = std::make_pair(0, 0.26);   // Y-axis zoom for 0-10% PbPb with nominal energy weight
    relativeZoom[1][iTrackPt] = std::make_pair(0, 0.16);   // Y-axis zoom for 10-30% PbPb with nominal energy weight
    relativeZoom[2][iTrackPt] = std::make_pair(0, 0.1);  // Y-axis zoom for 30-50% PbPb with nominal energy weight
    relativeZoom[3][iTrackPt] = std::make_pair(0, 0.1);  // Y-axis zoom for 50-90% PbPb with nominal energy weight

    // Different zooming options for PbPb uncertainties with the energy weight squared
    if(uncertaintyCard[kPbPb]->GetWeightExponent() == 2){
      relativeZoom[0][iTrackPt] = std::make_pair(0, 0.5);   // Y-axis zoom for 0-10% PbPb with energy weight squared
      relativeZoom[1][iTrackPt] = std::make_pair(0, 0.3);   // Y-axis zoom for 10-30% PbPb with energy weight squared
      relativeZoom[2][iTrackPt] = std::make_pair(0, 0.18);  // Y-axis zoom for 30-50% PbPb with energy weight squared
      relativeZoom[3][iTrackPt] = std::make_pair(0, 0.15);  // Y-axis zoom for 50-90% PbPb with energy weight squared
    } 
  }

  relativeZoom[0][5] = std::make_pair(0, 0.22);   // Y-axis zoom for 0-10% and 3 GeV track cut PbPb with nominal energy weight
  relativeZoom[1][5] = std::make_pair(0, 0.12);   // Y-axis zoom for 10-30% and 3 GeV track cut PbPb with nomimal energy weight

  // Is some cases, we need to tune the zooming options for different track pT bins to better see the figures
  if(uncertaintyCard[kPbPb]->GetWeightExponent() == 2){
    relativeZoom[0][3] = std::make_pair(0, 0.5);   // Y-axis zoom for 0-10% and 2 GeV track cut PbPb with energy weight squared
    relativeZoom[0][5] = std::make_pair(0, 0.3);   // Y-axis zoom for 0-10% and 3 GeV track cut PbPb with energy weight squared
    relativeZoom[1][3] = std::make_pair(0, 0.25);   // Y-axis zoom for 10-30% and 2 GeV track cut PbPb with energy weight squared
    relativeZoom[1][5] = std::make_pair(0, 0.18);   // Y-axis zoom for 10-30% and 3 GeV track cut PbPb with energy weight squared
  } 

  // =================================================== //
  // Read the histograms from the uncertainty organizers //
  // =================================================== //

  // Create histogram managers for result files and systematic uncertainty organizers for systematic uncertainty files
  SystematicUncertaintyOrganizer* uncertainties[kNDataTypes];
  
  for(int iDataType = 0; iDataType < kNDataTypes; iDataType++){

    // Create a new systematic uncertainty organizer
    uncertainties[iDataType] = new SystematicUncertaintyOrganizer(uncertaintyFile[iDataType]);

  }
 
  // Define histograms for absolute and relative systematic uncertianties
  TH1D* systematicUncertaintyForPbPb[SystematicUncertaintyOrganizer::knUncertaintySources][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* systematicUncertaintyForPp[SystematicUncertaintyOrganizer::knUncertaintySources][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* relativeUncertaintyForPbPb[SystematicUncertaintyOrganizer::knUncertaintySources][nCentralityBins][nJetPtBinsEEC][nTrackPtBinsEEC];
  TH1D* relativeUncertaintyForPp[SystematicUncertaintyOrganizer::knUncertaintySources][nJetPtBinsEEC][nTrackPtBinsEEC];

  // Initialize histograms to NULL
  for(int iJetPt = 0; iJetPt < nJetPtBinsEEC; iJetPt++){
    for(int iTrackPt = 0; iTrackPt < nTrackPtBinsEEC; iTrackPt++){
      for(int iUncertainty = 0; iUncertainty < SystematicUncertaintyOrganizer::knUncertaintySources; iUncertainty++){
        systematicUncertaintyForPp[iUncertainty][iJetPt][iTrackPt] = NULL;
        relativeUncertaintyForPp[iUncertainty][iJetPt][iTrackPt] = NULL;
        for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
          systematicUncertaintyForPbPb[iUncertainty][iCentrality][iJetPt][iTrackPt] = NULL;
          relativeUncertaintyForPbPb[iUncertainty][iCentrality][iJetPt][iTrackPt] = NULL;
        } // Centrality loop
      } // Uncertainty type loop
    } // Track pT loop
  } // Jet pT loop

  // Read the histograms from managers
  int iTrackPtMatchedPp;
  int iJetPtMatchedPp;

  for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
    iJetPtMatchedPp = uncertaintyCard[kPp]->FindBinIndexJetPtEEC(uncertaintyCard[kPbPb]->GetBinBordersJetPtEEC(iJetPt));
    for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++){
      iTrackPtMatchedPp = uncertaintyCard[kPp]->FindBinIndexTrackPtEEC(uncertaintyCard[kPbPb]->GetBinBordersTrackPtEEC(iTrackPt));

      for(int iUncertainty = 0; iUncertainty < SystematicUncertaintyOrganizer::knUncertaintySources; iUncertainty++){

        // If the studied uncertainty source is relevent for pp, read the histograms and calculate relative uncertainty
        if(uncertainties[kPp]->GetSystematicUncertaintyRelevancyForPp(iUncertainty)) {
          systematicUncertaintyForPp[iUncertainty][iJetPt][iTrackPt] = uncertainties[kPp]->GetSystematicUncertainty(0, iJetPtMatchedPp, iTrackPtMatchedPp, iUncertainty);
          relativeUncertaintyForPp[iUncertainty][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPp[iUncertainty][iJetPt][iTrackPt]->Clone(Form("relativeUncertaintyPp%d%d%d", iUncertainty, iJetPt, iTrackPt));
          transformToRelativeUncertainty(relativeUncertaintyForPp[iUncertainty][iJetPt][iTrackPt]);
        }

        // Next, read the PbPb histograms and calculate relative uncertainties if the uncertainty source is relevant for PbPb
        if(!uncertainties[kPbPb]->GetSystematicUncertaintyRelevancyForPbPb(iUncertainty)) continue;

        for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++) {

          systematicUncertaintyForPbPb[iUncertainty][iCentrality][iJetPt][iTrackPt] = uncertainties[kPbPb]->GetSystematicUncertainty(iCentrality, iJetPt, iTrackPt, iUncertainty);
          relativeUncertaintyForPbPb[iUncertainty][iCentrality][iJetPt][iTrackPt] = (TH1D*) systematicUncertaintyForPbPb[iUncertainty][iCentrality][iJetPt][iTrackPt]->Clone(Form("relativeUncertaintyPbPb%d%d%d%d", iUncertainty, iCentrality, iJetPt, iTrackPt));
          transformToRelativeUncertainty(relativeUncertaintyForPbPb[iUncertainty][iCentrality][iJetPt][iTrackPt]);

        } // Centrality loop
      } // Uncertainty source loop
    } // Track pT loop
  } // Jet pT loop

  // ======================================================= //
  //   Draw all the relative uncertainties in the same plot  //
  // ======================================================= //

  JDrawer* drawer = new JDrawer();
  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);
  drawer->SetLogX(true); // Logarithmic deltaR axis

  TString centralityString;
  TString compactCentralityString;
  TString jetPtString;
  TString compactJetPtString;
  TString trackPtString;
  TString compactTrackPtString;

  double legendX1, legendX2, legendY1, legendY2;

  // Common variables for different plots
  TLegend* legend;

  // Do the drawing

  for(int iJetPt = firstDrawnJetPtBinEEC; iJetPt <= lastDrawnJetPtBinEEC; iJetPt++){
    jetPtString = Form("%.0f < jet p_{T} < %.0f", uncertaintyCard[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), uncertaintyCard[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));
    compactJetPtString = Form("_J=%.0f-%.0f", uncertaintyCard[kPbPb]->GetLowBinBorderJetPtEEC(iJetPt), uncertaintyCard[kPbPb]->GetHighBinBorderJetPtEEC(iJetPt));

    for(int iTrackPt = firstDrawnTrackPtBinEEC; iTrackPt <= lastDrawnTrackPtBinEEC; iTrackPt++) {
      trackPtString = Form("%.1f < track p_{T}", uncertaintyCard[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));
      compactTrackPtString = Form("_T>%.1f", uncertaintyCard[kPbPb]->GetLowBinBorderTrackPtEEC(iTrackPt));
      compactTrackPtString.ReplaceAll(".", "v");

      // Setup the legend for plots
      legendX1 = 0.41; legendX2 = 0.91; legendY1 = 0.45; legendY2 = 0.89;
      if(uncertaintyCard[kPp]->GetWeightExponent() == 2){
        legendX1 = 0.43; legendX2 = 0.93; legendY1 = 0.55; legendY2 = 0.98;
      }
      legend = new TLegend(legendX1, legendY1, legendX2, legendY2);
      legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.04); legend->SetTextFont(62);
      legend->AddEntry((TObject*)0, "pp data", "");
      legend->AddEntry((TObject*)0, jetPtString.Data(), "");
      legend->AddEntry((TObject*)0, trackPtString.Data(), "");

      // Set the drawing style for histogram containing all uncertainties
      relativeUncertaintyForPp[SystematicUncertaintyOrganizer::kAll][iJetPt][iTrackPt]->SetLineColor(uncertainties[kPp]->GetUncertaintyColor(SystematicUncertaintyOrganizer::kAll));

      // Set the x and y-axis drawing ranges
      relativeUncertaintyForPp[SystematicUncertaintyOrganizer::kAll][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
      relativeUncertaintyForPp[SystematicUncertaintyOrganizer::kAll][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(relativeZoom[nCentralityBins][iTrackPt].first, relativeZoom[nCentralityBins][iTrackPt].second);

      // Draw the relative uncertainties from all sources for pp
      drawer->DrawHistogram(relativeUncertaintyForPp[SystematicUncertaintyOrganizer::kAll][iJetPt][iTrackPt], "#Deltar", "Relative uncertainty", " ");
      legend->AddEntry(relativeUncertaintyForPp[SystematicUncertaintyOrganizer::kAll][iJetPt][iTrackPt], "Total uncertainty", "l");

      // Draw the different uncertainty sources to the same plot
      for(int iUncertainty = 0; iUncertainty < SystematicUncertaintyOrganizer::kAll; iUncertainty++){
        if(!uncertainties[kPp]->GetSystematicUncertaintyRelevancyForPp(iUncertainty)) continue;
        relativeUncertaintyForPp[iUncertainty][iJetPt][iTrackPt]->SetLineColor(uncertainties[kPp]->GetUncertaintyColor(iUncertainty));
        relativeUncertaintyForPp[iUncertainty][iJetPt][iTrackPt]->Draw("same");
        legend->AddEntry(relativeUncertaintyForPp[iUncertainty][iJetPt][iTrackPt], uncertainties[kPp]->GetSystematicUncertaintyLegendName(iUncertainty), "l");
      }

      legend->Draw();

      // If a plot name is given, save the plot in a file
      if(saveFigures){
        gPad->GetCanvas()->SaveAs(Form("figures/energyEnergyCorrelator_relativeUncertaintySummary%s_pp%s%s.pdf", saveComment.Data(), compactJetPtString.Data(), compactTrackPtString.Data()));
      }

      // Draw the same plots also for PbPb
      for(int iCentrality = firstDrawnCentralityBin; iCentrality <= lastDrawnCentralityBin; iCentrality++){
        centralityString = Form("PbPb data, %.0f-%.0f%%", uncertaintyCard[kPbPb]->GetLowBinBorderCentrality(iCentrality), uncertaintyCard[kPbPb]->GetHighBinBorderCentrality(iCentrality));
        compactCentralityString = Form("_C=%.0f-%.0f", uncertaintyCard[kPbPb]->GetLowBinBorderCentrality(iCentrality), uncertaintyCard[kPbPb]->GetHighBinBorderCentrality(iCentrality));

        // Setup the legend for plots
        legendX1 = 0.37; legendX2 = 0.87; legendY1 = 0.4; legendY2 = 0.9;
        if(uncertaintyCard[kPbPb]->GetWeightExponent() == 2){
          legendX1 = 0.34; legendX2 = 0.84; legendY1 = 0.4; legendY2 = 0.9;
        }
        legend = new TLegend(legendX1, legendY1, legendX2, legendY2);
        legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.04); legend->SetTextFont(62);
        legend->AddEntry((TObject*)0, centralityString.Data(), "");
        legend->AddEntry((TObject*)0, jetPtString.Data(), "");
        legend->AddEntry((TObject*)0, trackPtString.Data(), "");

        // Set the drawing style for histogram containing all uncertainties
        relativeUncertaintyForPbPb[SystematicUncertaintyOrganizer::kAll][iCentrality][iJetPt][iTrackPt]->SetLineColor(uncertainties[kPbPb]->GetUncertaintyColor(SystematicUncertaintyOrganizer::kAll));

        // Set the x and y-axis drawing ranges
        relativeUncertaintyForPbPb[SystematicUncertaintyOrganizer::kAll][iCentrality][iJetPt][iTrackPt]->GetXaxis()->SetRangeUser(analysisDeltaR.first, analysisDeltaR.second);
        relativeUncertaintyForPbPb[SystematicUncertaintyOrganizer::kAll][iCentrality][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(relativeZoom[iCentrality][iTrackPt].first, relativeZoom[iCentrality][iTrackPt].second);

        // Draw the relative uncertainties from all sources for PbPb
        drawer->DrawHistogram(relativeUncertaintyForPbPb[SystematicUncertaintyOrganizer::kAll][iCentrality][iJetPt][iTrackPt], "#Deltar", "Relative uncertainty", " ");
        legend->AddEntry(relativeUncertaintyForPbPb[SystematicUncertaintyOrganizer::kAll][iCentrality][iJetPt][iTrackPt], "Total uncertainty", "l");

        // Draw the different uncertainty sources to the same plot
        for(int iUncertainty = 0; iUncertainty < SystematicUncertaintyOrganizer::kAll; iUncertainty++){
          if(!uncertainties[kPbPb]->GetSystematicUncertaintyRelevancyForPbPb(iUncertainty)) continue;
          relativeUncertaintyForPbPb[iUncertainty][iCentrality][iJetPt][iTrackPt]->SetLineColor(uncertainties[kPbPb]->GetUncertaintyColor(iUncertainty));
          relativeUncertaintyForPbPb[iUncertainty][iCentrality][iJetPt][iTrackPt]->Draw("same");
          legend->AddEntry(relativeUncertaintyForPbPb[iUncertainty][iCentrality][iJetPt][iTrackPt], uncertainties[kPbPb]->GetSystematicUncertaintyLegendName(iUncertainty), "l");
        }

        legend->Draw();

        // If a plot name is given, save the plot in a file
        if(saveFigures){
          gPad->GetCanvas()->SaveAs(Form("figures/energyEnergyCorrelator_relativeUncertaintySummary%s%s%s%s.pdf", saveComment.Data(), compactCentralityString.Data(), compactJetPtString.Data(), compactTrackPtString.Data()));
        }

      }

    } // Track pT loop
  } // Jet pT loop
}
