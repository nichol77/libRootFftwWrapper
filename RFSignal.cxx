#include <iostream>
#include "RFSignal.h"
#include "FFTtools.h"

//Maybe need to add and RFSignal::RFSignal(RFSignal thingy)

RFSignal::RFSignal() 
  : fGotFreqs(0),fMvNs(0)
{
//Default Constructor
  fGotFreqs=0;
  fComplexNums=0;
  fFreqs=0;
  fPhases=0;
  fMags=0;
}


RFSignal::RFSignal(TGraph *grWave,Int_t mvNs)
  :TGraph(*grWave),fGotFreqs(0),fMvNs(mvNs)
{
  ///<Assignnment constructor
  fGotFreqs=0;
  fComplexNums=0;
  fFreqs=0;
  fPhases=0;
  fMags=0;
  fillFreqStuff();
}

RFSignal::RFSignal(Int_t numFreqs, Double_t *freqVals, FFTWComplex *complexNums,Int_t mvNs) 
  :TGraph(2*(numFreqs-1)),fMvNs(mvNs)
{
  ///  std::cerr << "Here\t" << numFreqs << "\t" << fNpoints << "\n";
  fNumFreqs=numFreqs;
  fComplexNums = new FFTWComplex [numFreqs];
  fFreqs = new Double_t[numFreqs];
  //  std::cerr << "New fFreqs RFSignal::RFSignal(numFreqs,...)\t" << fFreqs << "\n";

  fPhases = new Double_t[numFreqs];
  fMags = new Double_t[numFreqs];
  //  std::cerr << "Here still\t" << numFreqs << "\t" << fNpoints << "\n";
  for(int i=0;i<fNumFreqs;i++) {
    //    std::cerr << i << "\t" << fFreqs << "\t" << fComplexNums << ;
    fFreqs[i]=freqVals[i];
    fComplexNums[i]=complexNums[i];
    fPhases[i]=fComplexNums[i].getPhase();
    fMags[i]=fComplexNums[i].getAbsSq();
  }
  extractFromComplex();
}

RFSignal::RFSignal(Int_t numPoints,Double_t *tVals,Double_t *vVals,Int_t mvNs)
  :TGraph(numPoints,tVals,vVals),fGotFreqs(0),fMvNs(mvNs)
{
  ///<Assignnment constructor
  fGotFreqs=0;
  fComplexNums=0;
  fFreqs=0;
  fPhases=0;
  fMags=0;
  fillFreqStuff();
}

RFSignal::~RFSignal() 
{
//Default destructor
  if(fGotFreqs) {
    //    std::cerr << fFreqs << "\t" << fMags << "\t" << fPhases 
    //	      << "\t" << fComplexNums << "\n";
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
  if(fMvNs)
    deltaT/=1e9; //Convert to s
  int length=fNpoints; //From TGraph
  fComplexNums=FFTtools::doFFT(length,fY);
  
  fNumFreqs=(length/2)+1;  
  fMags = new double [fNumFreqs];
  fFreqs = new double [fNumFreqs];
  //  std::cerr << "New fFreqs RFSignal::fillFreqStuff\t" << fFreqs << "\n";
  fPhases = new Double_t[fNumFreqs];

  //    double fMax = 1/(2*deltaT);  // In Hz
  double deltaF=1/(deltaT*length); //Hz
  if(fMvNs)
    deltaF/=1e6; 

  //  deltaF*=1e-6; //MHz //Lets keep frequency in hz
  
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
    fPhases[i]=fComplexNums[i].getPhase();
    tempF+=deltaF;
  }
}

Int_t RFSignal::getNumFreqs()
{
  if(!fGotFreqs)
    fillFreqStuff();
  return fNumFreqs;
}


void RFSignal::extractFromComplex()
{
  //  std::cout << "Here also\n";
  fGotFreqs=1;
  double deltaF=fFreqs[1]-fFreqs[0];
  double deltaT=1./(deltaF*fNpoints);
  if(fMvNs)
    deltaT*=1e3;
  double temp=0;
  //  std::cout << "RFSignal: " << fNpoints << "\t" << fComplexNums[1].getAbs()
    //	    << "\n";
  double *fVoltVals = FFTtools::doInvFFT(fNpoints,fComplexNums);
  //  std::cout << "RFSignal: V " << fVoltVals[0] << "\t" << fVoltVals[1] 
  //	    << "\n";
  for(int i=0;i<fNpoints;i++) {
    fX[i]=temp;
    temp+=deltaT;
    fY[i]=fVoltVals[i];
//     if(i<fNpoints/2) {
//       fY[i]=fVoltVals[fNpoints/2+i];
//     }
//     else {
//       fY[i]=fVoltVals[i-fNpoints/2];
//     }
  }
}

void RFSignal::addToSignal(RFSignal *grSignal)
{
  if(grSignal->GetN()!=this->GetN()) {
    std::cout << "Different RFSignal sizes can't add\n";
    return;
  }
  FFTWComplex *otherNums = grSignal->getComplexNums();
  for(int i=0;i<fNumFreqs;i++) {
    fComplexNums[i]+=otherNums[i];
  }
  extractFromComplex();
}


