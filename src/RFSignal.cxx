#include <iostream>
#include "RFSignal.h"
#include "RFFilter.h"
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


RFSignal::RFSignal(RFSignal *rfWave)
  :TGraph(*((TGraph*)rfWave))
{
  ///<Assignnment constructor
  fGotFreqs=1;
  fNumFreqs=rfWave->getNumFreqs();
  fComplexNums = new FFTWComplex [fNumFreqs];
  fFreqs = new Double_t[fNumFreqs];
  fPhases = new Double_t[fNumFreqs];
  fMags = new Double_t[fNumFreqs];
  
  for(int i=0;i<fNumFreqs;i++) {
    fComplexNums[i]=rfWave->getComplexNums()[i];
    fFreqs[i]=rfWave->getFreqs()[i];
    fPhases[i]=rfWave->getPhases()[i];
    fMags[i]=rfWave->getMags()[i];
  }  
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
  //  fillFreqStuff();
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
    fMags[i]=fComplexNums[i].getAbs();
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
  //  fillFreqStuff();
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
    //Grrrh do we want mags or not, I say yes to mags
    fMags[i]=fComplexNums[i].getAbs();
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


void RFSignal::updateTimeDomain(){
  reExtractFromComplex();
}

void RFSignal::extractFromComplex()
{
//   std::cout << "EXTRACT:\n";
  //  std::cout << "Here also\n";
  fGotFreqs=1;
  double deltaF=fFreqs[1]-fFreqs[0];
  double deltaT=1./(deltaF*fNpoints);
  if(fMvNs)
    deltaT*=1e3;
  double temp=0;
//   std::cout << "RFSignal: " << fNpoints << "\t" << fComplexNums[1].getAbs()
//     	    << "\n";
//   std::cout << "RFSignal: dT " << deltaT << " dF " << deltaF
// 	    << "\n";
  double *fVoltVals = FFTtools::doInvFFT(fNpoints,fComplexNums);
//   std::cout << "RFSignal: V " << fVoltVals[0] << "\t" << fVoltVals[1] 
// 	    << "\n";
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
  delete [] fVoltVals;
}

void RFSignal::addToSignal(RFSignal *grSignal)
{
  if(grSignal->GetN()!=this->GetN()) {
    std::cout << "Different RFSignal sizes can't add\n";
    std::cout << "Adding " << grSignal->GetN() << " to " << this->GetN() << "\n";
    std::cout << "wave adding " << grSignal->getNumFreqs() << " to " << this->getNumFreqs() << "\n";
    return;
  }
  FFTWComplex *otherNums = grSignal->getComplexNums();
  for(int i=0;i<fNumFreqs;i++) {
    fComplexNums[i]+=otherNums[i];
    fPhases[i]=fComplexNums[i].getPhase();
    fMags[i]=fComplexNums[i].getAbs();
  }
  extractFromComplex();
}





void RFSignal::applyFilter(RFFilter *theFilter)
{
  Int_t filterNumFreqs=theFilter->getNumFreqs();
  if (filterNumFreqs==0) {
    std::cerr<<"[RFFilter::filter] Error!  Trying to use an uninitialized filter!\n";
    return;
  }
  Double_t *filterFreqs=theFilter->getFreqs();
  Double_t *filterMags=theFilter->getMags();
  Double_t *filterPhases=theFilter->getPhases();

  double temp_mag;
  double temp_phase;
  int curr_index=0;
//   std::cout << "numFreqs " << fNumFreqs << " numTime " << fNpoints << std::endl;
//   std::cout << "filterNumFreqs " << filterNumFreqs << std::endl;
  for (int i=0;i<fNumFreqs;++i)
    { 

//       if(filterFreqs[i]>200e6 && filterFreqs[i]<300e6)
// 	std::cout << "freq " << filterFreqs[i] << "\tsignal freq " << fFreqs[i] << " mag " << fMags[i] << " phase " << fPhases[i] << std::endl;

      //Find the signal frequency in the filter array
      while (curr_index+1 < filterNumFreqs && filterFreqs[curr_index+1] < fFreqs[i])
	curr_index++;

      if (curr_index+1 < filterNumFreqs)
	{
	  if (!(filterFreqs[curr_index] <= fFreqs[i] && filterFreqs[curr_index+1] >= fFreqs[i]))
	    { //Error checking.
	      std::cout<<"Error!  Found incorrect frequency inside RFFilter operator()!  Filter not applied.\n";
	      std::cout<<"curr_index = "<<curr_index<<", filterFreqs[curr_index] = "<<filterFreqs[curr_index]<<", filterNumFreqs = "<<filterNumFreqs<<", fFreqs["<<i<<"] = "<<fFreqs[i]<<std::endl;
	      return;
	    }

	  if (filterFreqs[curr_index+1] != fFreqs[i])
	    {
	      temp_mag = FFTtools::simpleInterploate(filterFreqs[curr_index],filterMags[curr_index],
						     filterFreqs[curr_index+1],filterMags[curr_index+1],
						     fFreqs[i]);
	      temp_phase = FFTtools::simpleInterploate(filterFreqs[curr_index],filterPhases[curr_index],
						       filterFreqs[curr_index+1],filterPhases[curr_index+1],
						       fFreqs[i]);
	      //This is probably wrong for filterPhases, will fix later
	  
	      fMags[i] *= temp_mag;
	      fPhases[i] += temp_phase;
	    }
	  else
	    {
	      fMags[i] *= filterMags[curr_index+1];
	      fPhases[i] += filterPhases[curr_index+1];
	    }
	} //if requested frequency is covered by the filter
      else
	{
	  //Set signal to 0 in all regions not covered by filter.
	  fMags[i] = 0;
	  fPhases[i] = 0;
	}
      fComplexNums[i].setMagPhase(fMags[i],fPhases[i]);

//       if(filterFreqs[i]>200e6 && filterFreqs[i]<300e6)
// 	std::cout << "\t\t\tsignal freq " << fFreqs[i] << " mag " << fMags[i] << " phase " << fPhases[i] << std::endl;
      //std::cout<<"sig_mag["<<i<<"] = "<<sig_mag[i]<<", sig_phase["<<i<<"] = "<<sig_phase[i]<<std::endl;
    } //for

//   std::cout << "numFreqs " << fNumFreqs << " numTime " << fNpoints << std::endl;
//   std::cout << "filterNumFreqs " << filterNumFreqs << std::endl;
  reExtractFromComplex();
}




void RFSignal::reExtractFromComplex()
{
//   std::cout << "REEXTRACT:\n";
  //  std::cout << "Here also\n";
  fGotFreqs=1;
//  double deltaF=fFreqs[1]-fFreqs[0];
//  double deltaT=1./(deltaF*fNpoints);
//   if(fMvNs)
//     deltaT*=1e3;
//  double temp=0;
//   std::cout << "RFSignal: " << fNpoints << "\t" << fComplexNums[1].getAbs()
//     	    << "\n";
//   std::cout << "RFSignal: dT " << deltaT << " dF " << deltaF
// 	    << "\n";
  double *fVoltVals = FFTtools::doInvFFT(fNpoints,fComplexNums);
//   std::cout << "RFSignal: V " << fVoltVals[0] << "\t" << fVoltVals[1] 
// 	    << "\n";
  for(int i=0;i<fNpoints;i++) {
//     fX[i]=temp;
//     temp+=deltaT;
    fY[i]=fVoltVals[i];
//     if(i<fNpoints/2) {
//       fY[i]=fVoltVals[fNpoints/2+i];
//     }
//     else {
//       fY[i]=fVoltVals[i-fNpoints/2];
//     }
  }
  delete [] fVoltVals;
}
