
#include "RFSignal.h"
#include "FFTtools.h"

RFSignal::RFSignal() 
  : fGotFreqs(0)
{
//Default Constructor
}


RFSignal::RFSignal(TGraph *grWave)
 :TGraph(*grWave),fGotFreqs(0)
{
  ///<Assignnment constructor
}

RFSignal::RFSignal(Int_t numPoints,Double_t *tVals,Double_t *vVals)
  :TGraph(numPoints,tVals,vVals),fGotFreqs(0)
{
  ///<Assignnment constructor
}

RFSignal::~RFSignal() 
{
//Default destructor
  if(fGotFreqs) {
    delete [] fFreqs;
    delete [] fMags;
    delete [] fPhases;
    delete [] fComplexNums;
  }
}

TGraph *RFSignal::getFreqMagGraph()
{ 
  if(!fGotFreqs)
    fillFreqStuff();
  TGraph *grFreq = new TGraph(fNumFreqs,fFreqs,fMags);
  return grFreq;
}

Double_t *RFSignal::getFreqs()
{
  if(!fGotFreqs)
    fillFreqStuff();
  return fFreqs;
}

Double_t *RFSignal::getMags()
{
  if(!fGotFreqs)
    fillFreqStuff();
  return fMags;
}

Double_t *RFSignal::getPhases()
{
  if(!fGotFreqs)
    fillFreqStuff();
  return fPhases;
}


FFTWComplex *RFSignal::getComplexNums()
{
  if(!fGotFreqs)
    fillFreqStuff();
  return fComplexNums;
}

void RFSignal::fillFreqStuff()
{
  fGotFreqs=1;
  //Calcualte FFT and then...
  double deltaT=fX[1]-fX[0]; //From TGraph
  int length=fNpoints; //From TGraph
  fComplexNums=FFTtools::doFFT(length,fY);
  
  fNumFreqs=(length/2)+1;  
  fMags = new double [fNumFreqs];
  fFreqs = new double [fNumFreqs];

  //    double fMax = 1/(2*deltaT);  // In Hz
  double deltaF=1/(deltaT*length); //Hz
  deltaF*=1e-6; //MHz
  
  double tempF=0;
  for(int i=0;i<fNumFreqs;i++) {
    double power=fComplexNums[i].getAbsSq();
    if(i>0 && i<fNumFreqs-1) power*=2; //account for symmetry
    //	power*=deltaT/(length); //For time-integral squared amplitude
    //	power/=deltaF;//Just to normalise bin-widths
    //Ends up the same as dt^2, need to integrate the power (multiply by df)
    //to get a meaningful number out.    
    fFreqs[i]=tempF;
    fMags[i]=power;
    tempF+=deltaF;
  }
}

Int_t RFSignal::getNumFreqs()
{
  if(!fGotFreqs)
    fillFreqStuff();
  return fNumFreqs;
}
