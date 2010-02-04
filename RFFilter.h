#ifndef RFFILTER_H
#define RFFILTER_H
/////////////////////////////////////////////////
//  RFFilter  v1.3.2, 4 December 2009
//
//   This class builds a causal filter in the frequency domain and applies it
// to FrequencyDomain or RFWaveform objects. 
//
// Changes:
//   v1.1:   Added check that the lower frequency is positive to "constructFilter".
//           Added a constructor to mask a band between two passbands.
//   v1.2:   Fixed two array overflow bugs in "constructFilter".  Both bugs were in the for
//             loop that assigns coefficients to the "complex_value" Array.  Bug 1: The if
//             statement was <= freq_high instead of < freq_high.  Bug 2: I did not treat
//             the frequency=0 case separately, as I should.
//   v1.3:   Added filter type 2, which doesn't attempt to be causal.
//           Added operator *= to multiply the coefficients of two filters together.
//   v1.3.1: 1 July 2008 : Added some error checking to the second constructor - now checking if some of the
//             intervals don't make sense.
//   v1.3.2: 4 Dec 2009 : Fixed void constructFilter declaration (int coefficients=400 -> int num_coeff=400)
//
// Stephen Hoover, hoover@physics.ucla.edu
/////////////////////////////////////////////////

#include "RFSignal.h"
#include "FFTWComplex.h"

class RFFilter {
 public:
  RFFilter();
  virtual ~RFFilter();

  //This constructor takes an array of filter coefficients in the frequency domain, and turns it into a
  //causal filter. 
  //numFreqs is the number of frequencies
  //freqs is an evenly spaced array of frequencies
  //coefficients is the array of coefficients
    //filter_type: 0: Square window, 1: Hanning window (recommend using filter type 1).
  RFFilter(int numFreqs,double *freqs, double *coefficients, int filter_type=1, int maxCoeff=400, int debugMode=0);

  //This constructor builds a flat bandpass filter between low and high (in Hz).  Other arguements are as above.
  //  RFFilter(double low, double high, double df=1.0e6, int power_of_2=14, int num_coeff=400, int filter_type=1);

  //  RFFilter(double low1, double low2, double high1, double high2, double df=1.0e6, int power_of_2=14, int num_coeff=400, int filter_type=1);

/*   RFFilter(const RFFilter &rhs); //copy constructor */


/*   RFFilter& operator=(const RFFilter &rhs); */

/*   RFFilter & operator*=(const RFFilter &rhs); */

  //These two methods apply the filter to a signal.  
  RFSignal *filter(RFSignal *signal);
  void filter(Int_t numFreqs,double *freqs,FFTWComplex *complexVals); ///This one here changes the input arrays
  int getDebugMode() {return fDebugMode;}


 private:
  //Does the work of constructing the filter, so the two similar constructors don't have to repeat code.
  void constructFilter(int numFreqs,double *freqs, double *coefficients, int filter_type=1, int maxCoeff=400);

 private:
  int fNumFreq; //Size of arrays, only positive frequencies are stored
  int fNumTimes; // 2*(fNumFreqs-1)
  double fFrequencyInterval; //In Hz
  double *fFrequency;
  double *fMagnitude;
  double *fPhase;
  FFTWComplex *fComplexNums;
  Int_t fDebugMode;

};


#endif
