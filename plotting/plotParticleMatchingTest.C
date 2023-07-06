#include "JDrawer.h"

/*
 * Macro for plotting histograms from particle matching test
 */
void plotParticleMatchingTest(){

  // ============= //
  // Configuration //
  // ============= //
  
  // Input files
  TFile* inputFile = TFile::Open("data/PbPbMC2018_RecoGen_akFlowJets_miniAOD_4pCentShift_noTrigger_nominalSmear_onlyTrackParticleMatching_2023-07-03.root");
  TFile* ppFile = TFile::Open("data/ppMC2017_RecoGen_Pythia8_pfJets_wtaAxis_32deltaRBins_nominalSmear_onlyTrackParticleMatching_2023-07-03.root");

  // Read the THnSparse from the file
  THnSparseF* histogramArray = (THnSparseF*) inputFile->Get("particlesCloseToTracks");
  THnSparseF* binaryHistogramArray = (THnSparseF*) inputFile->Get("tracksWithMatchedParticle");
  THnSparseF* ppHistogramArray = (THnSparseF*) ppFile->Get("particlesCloseToTracks");
  THnSparseF* ppBinaryHistogramArray = (THnSparseF*) ppFile->Get("tracksWithMatchedParticle");

  // Centrality and jet pT bins
  const int nCentralityBins = 4;
  const int nJetPtBins = 4;
  const int nTrackPtCuts = 2;
  int centralityBinBorders[nCentralityBins+1] = {0,10,30,50,90};
  int jetPtBinBorders[nJetPtBins+1] = {120, 140, 160, 180, 200};
  double trackPtCutValues[nTrackPtCuts] = {2, 3};
  int lowTrackBin = 0;
  int highTrackBin = 0;

  // Project the histograms from the array
  TH1D* particleCountNearTrack[nCentralityBins+1][nJetPtBins][nTrackPtCuts];
  TH1D* hasMatchingParticle[nCentralityBins+1][nJetPtBins][nTrackPtCuts];
  for(int iCentrality = 0; iCentrality < nCentralityBins; iCentrality++){
    histogramArray->GetAxis(3)->SetRange(iCentrality+1,iCentrality+1);
    binaryHistogramArray->GetAxis(3)->SetRange(iCentrality+1,iCentrality+1);
    for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){
      histogramArray->GetAxis(1)->SetRange(iJetPt+4,iJetPt+4);
      binaryHistogramArray->GetAxis(1)->SetRange(iJetPt+4,iJetPt+4);
      for(int iTrackPt = 0; iTrackPt < nTrackPtCuts; iTrackPt++){
        lowTrackBin = histogramArray->GetAxis(2)->FindBin(trackPtCutValues[iTrackPt]);
        highTrackBin = histogramArray->GetAxis(2)->GetNbins();
        histogramArray->GetAxis(2)->SetRange(lowTrackBin,highTrackBin);
        binaryHistogramArray->GetAxis(2)->SetRange(lowTrackBin,highTrackBin);
        particleCountNearTrack[iCentrality][iJetPt][iTrackPt] = (TH1D*) histogramArray->Projection(0);
        particleCountNearTrack[iCentrality][iJetPt][iTrackPt]->SetName(Form("myTracks%d%d%d", iCentrality, iJetPt, iTrackPt));
        hasMatchingParticle[iCentrality][iJetPt][iTrackPt] = (TH1D*) binaryHistogramArray->Projection(0);
        hasMatchingParticle[iCentrality][iJetPt][iTrackPt]->SetName(Form("myMatch%d%d%d", iCentrality, iJetPt, iTrackPt));
      }
    }
  }

  // Project the pp histograms
  for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++){
    ppHistogramArray->GetAxis(1)->SetRange(iJetPt+3,iJetPt+3);
    ppBinaryHistogramArray->GetAxis(1)->SetRange(iJetPt+3,iJetPt+3);
    for(int iTrackPt = 0; iTrackPt < nTrackPtCuts; iTrackPt++){
      lowTrackBin = ppHistogramArray->GetAxis(2)->FindBin(trackPtCutValues[iTrackPt]);
      highTrackBin = ppHistogramArray->GetAxis(2)->GetNbins();
      ppHistogramArray->GetAxis(2)->SetRange(lowTrackBin,highTrackBin);
      ppBinaryHistogramArray->GetAxis(2)->SetRange(lowTrackBin,highTrackBin);
      particleCountNearTrack[nCentralityBins][iJetPt][iTrackPt] = (TH1D*) ppHistogramArray->Projection(0);
      particleCountNearTrack[nCentralityBins][iJetPt][iTrackPt]->SetName(Form("myTracksPp%d%d", iJetPt, iTrackPt));
      hasMatchingParticle[nCentralityBins][iJetPt][iTrackPt] = (TH1D*) ppBinaryHistogramArray->Projection(0);
      hasMatchingParticle[nCentralityBins][iJetPt][iTrackPt]->SetName(Form("myMatchPp%d%d", iJetPt, iTrackPt));
    }
  }

  // ========================================= //
  //   Draw the histograms in separate plots   //
  // ========================================= //

  JDrawer* drawer = new JDrawer();
  TLegend* legend;
  int color[] = {kBlack, kRed, kBlue, kGreen+3, kMagenta};

  for(int iTrackPt = 0; iTrackPt < nTrackPtCuts; iTrackPt++){
    for(int iJetPt = 0; iJetPt < nJetPtBins; iJetPt++) {
      legend = new TLegend(0.53, 0.5, 0.83, 0.85);
      legend->SetFillStyle(0); legend->SetBorderSize(0); legend->SetTextSize(0.05); legend->SetTextFont(62);
      legend->AddEntry((TObject*)0, Form("%d < jet pT < %d", jetPtBinBorders[iJetPt], jetPtBinBorders[iJetPt+1]), "");
      legend->AddEntry((TObject*)0, Form("%.0f < Track pT", trackPtCutValues[iTrackPt]), "");

      for(int iCentrality = 0; iCentrality < nCentralityBins+1; iCentrality++){
        particleCountNearTrack[iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iCentrality]);
        particleCountNearTrack[iCentrality][iJetPt][iTrackPt]->Scale(1.0 / particleCountNearTrack[iCentrality][iJetPt][iTrackPt]->Integral());
      }

      particleCountNearTrack[0][iJetPt][iTrackPt]->GetYaxis()->SetRangeUser(0,1);
      drawer->DrawHistogram(particleCountNearTrack[0][iJetPt][iTrackPt], "Counts", "Probability", " ");
      legend->AddEntry(particleCountNearTrack[0][iJetPt][iTrackPt], Form("C: %d-%d", centralityBinBorders[0], centralityBinBorders[1]), "l");

      for(int iCentrality = 1; iCentrality < nCentralityBins; iCentrality++) {
        particleCountNearTrack[iCentrality][iJetPt][iTrackPt]->Draw("same");
        legend->AddEntry(particleCountNearTrack[iCentrality][iJetPt][iTrackPt], Form("C: %d-%d", centralityBinBorders[iCentrality], centralityBinBorders[iCentrality+1]), "l");
      }

      particleCountNearTrack[nCentralityBins][iJetPt][iTrackPt]->Draw("same");
      legend->AddEntry(particleCountNearTrack[nCentralityBins][iJetPt][iTrackPt], "pp", "l");

      // Also information if a matching particle is found to the same figure
      for(int iCentrality = 0; iCentrality < nCentralityBins+1; iCentrality++){
        hasMatchingParticle[iCentrality][iJetPt][iTrackPt]->SetLineColor(color[iCentrality]);
        hasMatchingParticle[iCentrality][iJetPt][iTrackPt]->SetLineStyle(2);
        hasMatchingParticle[iCentrality][iJetPt][iTrackPt]->Scale(1.0 / hasMatchingParticle[iCentrality][iJetPt][0]->Integral());
        hasMatchingParticle[iCentrality][iJetPt][iTrackPt]->Draw("same");
      }

      legend->Draw();
    }
  }
}
