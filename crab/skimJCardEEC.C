void skimJCardEEC(TString fileName, int nEntries){
  TFile *skimmedFile = new TFile(fileName,"UPDATE");
  skimmedFile->cd("JCard");
  char vectorName[100];
  const int nCardEntries = 46;
  const char *cardEntries[nCardEntries] = {"GitHash","DataType","McCorrelationType","MatchJets","TriggerSelection","FilledHistograms","JetType","JetAxis","JetEtaCut","MinJetPtCut","MaxJetPtCut","CutBadPhi","MinMaxTrackPtFraction","MaxMaxTrackPtFraction","JetUncertainty","TrackEtaCut","MinTrackPtCut","MaxTrackPtCut","MaxTrackPtRelativeError","VertexMaxDistance","CalorimeterSignalLimitPt","HighPtEtFraction","Chi2QualityCut","MinimumTrackHits","SubeventCut","TrackEfficiencyVariation","JetPtWeight","DisableTrackPairEfficiencyCorrection","ZVertexCut","LowPtHatCut","HighPtHatCut","MultiplicityMode","JetRadius","WeightExponent","SmearDeltaR","SmearEnergyWeight","JetPtBinEdgesEEC","TrackPtBinEdgesEEC","JetPtBinEdgesUnfoldingReco","JetPtBinEdgesUnfoldingTruth","CentralityBinEdges","TrackPtBinEdges","PtHatBinEdges","DoReflectedCone","ApplyReflectedConeWeight","DebugLevel"};
  for(int iVector = 0; iVector < nCardEntries; iVector++){
    for(int j = 2; j <= nEntries; j++){
      sprintf(vectorName,"%s;%d",cardEntries[iVector],j);
      gDirectory->Delete(vectorName);
    }
  }
  skimmedFile->Close();
}
