#ifndef RFSIGNAL_H
#define RFSIGNAL_H
#include "TGraph.h"
#include "FFTWComplex.h"

class RFFilter;

//!  This is a wrapper class for an RF Signal
/*!
  At the moment it doesn't do very much but tis might change in the future
*/
class RFSignal : public TGraph {
public:
  RFSignal(); ///<Default constructor
  RFSignal(TGraph *grWave, Int_t mvNs=0); ///<Assignnment constructor
  RFSignal(RFSignal *rfWave); ///<Assignnment constructor
  RFSignal(Int_t numPoints,Double_t *tVals,Double_t *vVals, Int_t mvNs=0);
  RFSignal(Int_t numFreqs, Double_t *freqVals,FFTWComplex *complexNums, Int_t mvNs=0);
  ~RFSignal(); ///<Destructor
  
  TGraph *getFreqMagGraph();
  Double_t *getFreqs();
  Double_t *getMags();
  Double_t *getPhases();
  FFTWComplex *getComplexNums();
  Int_t getNumFreqs();
  void addToSignal(RFSignal *grSignal);
  void applyFilter(RFFilter *theFilter);

 private:  
  void fillFreqStuff();
  void extractFromComplex();
  Int_t fGotFreqs;
  Int_t fNumFreqs;
  Int_t fMvNs;
  Double_t *fFreqs;  //[fNumFreqs]
  Double_t *fMags; //[fNumFreqs]
  Double_t *fPhases; //[fNumFreqs]
  FFTWComplex *fComplexNums; //[fNumFreqs]
  
ClassDef(RFSignal,2) ///< ROOT macro for persistence. 

};

#endif // RFSIGNAL_H
