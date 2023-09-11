void checkNumberOfEventsMC(const char* fileName){
  TFile* dataFile = TFile::Open(fileName);
  TH1D* vzHistogram = (TH1D*) dataFile->Get("vertexZ");
  int nEntries = vzHistogram->GetEntries();

  int goodEntries[6] = {17764593, 17556674, 17764593+17556674, 21067938, 22279294, 21067938+22279294};
  
  cout << endl;
  cout << "Checking file " << fileName << endl;
  if(nEntries == goodEntries[0]){
    cout << "Yay! Entries consistent with part1 of Pythia+Hydjet simulation!" << endl;
  } else if (nEntries == goodEntries[1]){
    cout << "Yay! Entries consistent with part2 of Pythia+Hydjet simulation!" << endl;
  } else if (nEntries == goodEntries[2]){
    cout << "Yay! Entries consistent with whole Pythia+Hydjet simulation stats!" << endl;
  } else if (nEntries == goodEntries[3]){
    cout << "Yay! Entries consistent with part1 of Pythia8 simulation!" << endl;
  } else if (nEntries == goodEntries[4]){
    cout << "Yay! Entries consistent with part2 of Pythia8 simulation!" << endl;
  } else if (nEntries == goodEntries[5]){
    cout << "Yay! Entries consistent with whole Pythia8 simulation stats!" << endl;
  } else {
    cout << "ERROR! Entries do not match the expected!" << endl;
  }
  dataFile->Close();
}
