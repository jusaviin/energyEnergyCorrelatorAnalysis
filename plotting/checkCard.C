#include "EECCard.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)

/*
 * Macro for printing EECCard information from a given file
 */ 
 void checkCard(const char *fileName){
  TFile *file = TFile::Open(fileName);
  EECCard *card = new EECCard(file);
  card->Print();
  file->Close();
  delete card;
 }
