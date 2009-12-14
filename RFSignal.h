#ifndef RFSIGNAL_H
#define RFSIGNAL_H
#include "TGraph.h"
#include "FFTWComplex.h"


//!  This is a wrapper class for a complex number
/*!
  And that's it.
*/
class RFSignal : public TGraph {
public:
  RFSignal(); ///<Default constructor
  RFSignal(TGraph *grWave); ///<Assignnment constructor
  RFSignal(Int_t numPoints,Double_t *tVals,Double_t *vVals);
  ~RFSignal(); ///<Destructor
  
  TGraph *getFreqMagGraph();
  Double_t *getFreqs();
  Double_t *getMags();
  Double_t *getPhases();
  FFTWComplex *getComplexNums();
  Int_t getNumFreqs();

 private:  
  void fillFreqStuff();
  Int_t fGotFreqs;
  Int_t fNumFreqs;
  Double_t *fFreqs;
  Double_t *fMags;
  Double_t *fPhases;
  FFTWComplex *fComplexNums;
  
ClassDef(RFSignal,2) ///< ROOT macro for persistence. 

};

#endif // RFSIGNAL_H
