
//Stolen from
//FrequencyFilter.cxx  v1.3.2, 4 December 2009
//Stephen Hoover, hoover@physics.ucla.edu

#include <iostream>
#include <fstream>
#include <math.h>
#include "RFFilter.h"
#include "FFTtools.h"
#include "TGraph.h"


RFFilter::RFFilter()
  : fNumFreq(0),
    fNumTimes(0),
    fFrequencyInterval(0),
    fFrequency(0),
    fMagnitude(0),
    fPhase(0)
{}




RFFilter::RFFilter(int numFreqs,double *freqs, double *coefficients, int filter_type,int maxCoeff,int debugMode)
  : fNumFreq(numFreqs),
    fNumTimes(2*(numFreqs-1)),
    fFrequencyInterval(freqs[1]-freqs[0]),
    fDebugMode(debugMode)
{
  //Need to check here to make sure that the freqs array starts from zero
 
  this->constructFilter(numFreqs,freqs,coefficients,filter_type,maxCoeff);
}

//Filter is made causal by the following process:
//Start with all phases = 0, magnitudes given by input coefficients array, modified by Hanning window if desired.
//Transform to time domain, obtaining impulse response.
//Shift impulse response to later times, and truncate impulse response at t<0, and t>num_coeff.
//Transform shifted impulse response back to frequency domain.
void RFFilter::constructFilter(int numFreqs,double *freqs, double *coefficients, int filter_type, int maxCoeff)
{
  fNumFreq=numFreqs;
  fFrequency = new double [fNumFreq];
  fPhase = new double [fNumFreq];
  fMagnitude = new double [fNumFreq];
  fComplexNums = new FFTWComplex [fNumFreq];
  
  //Set the frequency and coefficients of each bin
  for(int i=0;i<fNumFreq;i++) {
    fFrequency[i]=freqs[i];
    fMagnitude[i]=coefficients[i];
    fComplexNums[i].re=coefficients[i];
    fPhase[i]=0;
  }

  /** Filter type 2: Acausal filter, square window. **/
  if (filter_type == 2)  {
    //Don't need to do anything else
    return;
  }

  double deltaT=1/(fFrequencyInterval*fNumTimes);  
  Double_t *fVoltVals = FFTtools::doInvFFT(fNumTimes,fComplexNums);
  Double_t *fTimeVals = new Double_t [fNumTimes];
  Double_t temp=0;
  for(int i=0;i<fNumTimes;i++) {
    fTimeVals[i]=temp;
    temp+=deltaT;
  }
 

  /* Apply windowing function to the filter's impulse response, and truncate series
     Window 0 : Square window, no change to nonzero coefficients
     Window 1 : Hanning window, where window function w(i) = 0.5 + 0.5cos(2pi*i / W) for -W/2<=i<=W/2, 0 elsewhere. */
  //First treat 0 as a special case.
  if(filter_type==1){} //Zeroth coefficient unchanged in Hanning window
  else if(filter_type==0){} //Unchanged in square window too
  //Now go through the rest of the coefficients, treating positive and negative cases simultaneously.
  //for the i=0 case, the coefficient is 1, no need to treat manually.
  for (int i=1;i<fNumTimes/2;i++) 
    {
      //      std::cout << i << "\t" << maxCoeff/2 << "\n";
    if (i <= maxCoeff/2) 
      {
      if(filter_type==1) //Hanning window
	{
	  //The t>0 coeffs
	  fVoltVals[i] *= (0.5 + 0.5 * cos(i * TMath::TwoPi() / maxCoeff));
	  //The t<0 coeffs
	  fVoltVals[fNumTimes-i] *= (0.5 + 0.5 * cos(i * TMath::TwoPi() / maxCoeff));
	} //if (Hanning window)
      else if (filter_type==0) //Do nothing for square window
	{
	} //else if (square window)
      } //if
    else //Restrict duration of impulse response by zeroing outside of specified number of parameters
      {
	fVoltVals[i]=0;
	fVoltVals[fNumTimes-i]=0;
      } //else
    } //for(apply window function)
  fVoltVals[fNumTimes/2]=0;

//   /* Now shift the impulse response function to make the filter causal. */
  for(int i=maxCoeff; i>=maxCoeff/2;i--) {
    fVoltVals[i]=fVoltVals[i-maxCoeff/2];
  }
  for(int i=0;i<maxCoeff/2;i++) {
    fVoltVals[i]=fVoltVals[(fNumTimes-maxCoeff/2)+i];
    fVoltVals[(fNumTimes-maxCoeff/2)+i]=0;
  }

  if(fDebugMode) {
    TGraph *grTimeDomain = new TGraph(fNumTimes,fTimeVals,fVoltVals);
    grTimeDomain->SetNameTitle("grTimeDomain","grTimeDomain");
    grTimeDomain->Write();
  }

  delete [] fComplexNums;
  fComplexNums = FFTtools::doFFT(fNumTimes,fVoltVals);
  
  for(int i=0;i<fNumFreq;i++) {
    fMagnitude[i]=fComplexNums[i].getAbs();
    fPhase[i]=fComplexNums[i].getPhase();
  }

  if(fDebugMode) {
    TGraph *grFreqMag = new TGraph(fNumFreq,fFrequency,fMagnitude);
    grFreqMag->SetNameTitle("grFreqMag","grFreqMag");
    grFreqMag->Write();
    

    TGraph *grFreqPhase = new TGraph(fNumFreq,fFrequency,fPhase);
    grFreqPhase->SetNameTitle("grFreqPhase","grFreqPhase");
    grFreqPhase->Write();
  }
  delete [] fVoltVals;
  delete [] fTimeVals;

} //RFFilter(double low, double high, int power_of_2, int maxCoeff, double df, int filter_type)


int RFFilter::updateFilter(int numFreqs,double *freqs,double *coefficients)
{
  if(numFreqs!=fNumFreq) {
    std::cout << "Have to make a new filter\n";
    return 0;
    //Could actually put something here that fixes this
  }
  for(int i=0;i<fNumFreq;i++) {
    fMagnitude[i]=coefficients[i];
    fComplexNums[i].setMagPhase(fMagnitude[i],fPhase[i]);
  }
  return 1;
  //Done
}

// RFFilter::RFFilter(const RFFilter & rhs) 
//   : N(rhs.N),
//     fNumFreq(rhs.fNumFreq),
//     fFrequencyInterval(rhs.fFrequencyInterval),
//     frequency(rhs.frequency),
//     magnitude(rhs.magnitude),
//     phase(rhs.phase)
// {} //RFFilter copy constructor


// RFFilter& RFFilter::operator=(const RFFilter &rhs)
// {
//   fFrequencyInterval = rhs.fFrequencyInterval;
//   N = rhs.N;
//   fNumFreq = rhs.fNumFreq;
//   frequency = Array(rhs.frequency);
//   magnitude = Array(rhs.magnitude);
//   phase = Array(rhs.phase);

//   return *this;
// } //RFFilter operator=


// RFFilter & RFFilter::operator*=(const RFFilter &rhs)
// {
//   if (this->N != rhs.N || this->fFrequencyInterval != rhs.fFrequencyInterval)
//     {
//       std::cerr<<"[RFFilter::operator*=]  Error!  Cannot multiply these two filters.  Filters must cover exactly the same frequencies to be multiplied!\n";
//       return *this;
//     } //end if

//   for (int index=0; index<frequency.size(); ++index)
//     {
//       magnitude[index] *= rhs.magnitude[index];
//       phase[index] += rhs.phase[index];
//     } //end for

//   return *this;
// } //RFFilter & RFFilter::operator*=(const RFFilter &rhs)


RFFilter::~RFFilter() 
{
  if(fNumFreq) {
    delete [] fFrequency;
    delete [] fMagnitude;
    delete [] fPhase;
    delete [] fComplexNums;
  }
 } //RFFilter destructor


// RFSignal RFFilter::operator()(RFSignal &signal)x
// {
//   Int_t numFreqs=signal.getNumFreqs();
//   Double_t *theFreqs= new Double_t [numFreqs];
//   FFTWComplex *theComplexVals=new FFTWComplex [numFreqs];
//   for(int i=0;i<numFreqs;i++) {
//     theFreqs[i]=signal.getFreqs()[i];
//     theComplexVals[i]=signal.getComplexNums()[i];
//   }
//   this->filter(numFreqs,theFreqs,theComplexVals);
//   return RFSignal(numFreqs,theFreqs,theComplexVals);
//   //Is this a memory leak?? Probably, will thinl about it
// } //RFSignal RFFilter::operator()(RFSignal &signal) const

RFSignal *RFFilter::filter(RFSignal *signal) 
{
  
  Int_t numFreqs=signal->getNumFreqs();
  Double_t *theFreqs= new Double_t [numFreqs];
  FFTWComplex *theComplexVals=new FFTWComplex [numFreqs];
  Double_t *inFreqs=signal->getFreqs();
  FFTWComplex *inComplex=signal->getComplexNums();
  for(int i=0;i<numFreqs;i++) {    
    theFreqs[i]=inFreqs[i];
    theComplexVals[i]=inComplex[i];
  }
  //HAck by RJN for testing
  this->filter(numFreqs,theFreqs,theComplexVals);
  RFSignal *retVal = new RFSignal(numFreqs,theFreqs,theComplexVals);
  //RFSignal *retVal = new RFSignal(signal->GetN(),signal->GetX(),signal->GetY());
  delete [] theFreqs;
  delete [] theComplexVals;
  return retVal;
} //RFSignal RFFilter::filter(RFSignal &signal) const



// FrequencyDomain RFFilter::filter(FrequencyDomain &signal) const
void RFFilter::filter(Int_t numFreqs,double *sig_freq,FFTWComplex *complexVals) ///This one here changes the input arrays
{
  if (fNumFreq==0) {
    std::cerr<<"[RFFilter::filter] Error!  Trying to use an uninitialized filter!\n";
    return;
  }
  double *sig_mag = new double [numFreqs];
  double *sig_phase = new double[numFreqs];

  double temp_mag;
  double temp_phase;
  int curr_index=0;
  for (int i=0;i<numFreqs;++i)
    { 
      sig_mag[i]=complexVals[i].getAbs();
      sig_phase[i]=complexVals[i].getPhase();

      //Find the signal frequency in the filter array
      while (curr_index+1 < fNumFreq && fFrequency[curr_index+1] < sig_freq[i])
	curr_index++;

      if (curr_index+1 < fNumFreq)
	{
	  if (!(fFrequency[curr_index] <= sig_freq[i] && fFrequency[curr_index+1] >= sig_freq[i]))
	    { //Error checking.
	      std::cout<<"Error!  Found incorrect frequency inside RFFilter operator()!  Filter not applied.\n";
	      std::cout<<"curr_index = "<<curr_index<<", fFrequency[curr_index] = "<<fFrequency[curr_index]<<", fNumFreq = "<<fNumFreq<<", sig_freq["<<i<<"] = "<<sig_freq[i]<<std::endl;
	      return;
	    }

	  if (fFrequency[curr_index+1] != sig_freq[i])
	    {
	      temp_mag = FFTtools::simpleInterploate(fFrequency[curr_index],fMagnitude[curr_index],
						     fFrequency[curr_index+1],fMagnitude[curr_index+1],
						     sig_freq[i]);
	      temp_phase = FFTtools::simpleInterploate(fFrequency[curr_index],fPhase[curr_index],
						       fFrequency[curr_index+1],fPhase[curr_index+1],
						       sig_freq[i]);
	      //This is probably wrong for fPhase, will fix later
	  
	      sig_mag[i] *= temp_mag;
	      sig_phase[i] += temp_phase;
	    }
	  else
	    {
	      sig_mag[i] *= fMagnitude[curr_index+1];
	      sig_phase[i] += fPhase[curr_index+1];
	    }
	} //if requested frequency is covered by the filter
      else
	{
	  //Set signal to 0 in all regions not covered by filter.
	  sig_mag[i] = 0;
	  sig_phase[i] = 0;
	}
      //std::cout<<"sig_mag["<<i<<"] = "<<sig_mag[i]<<", sig_phase["<<i<<"] = "<<sig_phase[i]<<std::endl;
      complexVals[i].setMagPhase(sig_mag[i],sig_phase[i]);
    } //for
  delete [] sig_phase;
  delete [] sig_mag;
} //FrequencyDomain RFFilter::filter(FrequencyDomain &signal) const
