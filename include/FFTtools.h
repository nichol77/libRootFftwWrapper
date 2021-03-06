
#ifndef FFTTOOLS_H
#define FFTTOOLS_H



// If the following is uncommented, in and out arrays for FFTW are allocated contiguously. 
#define FFTTOOLS_ALLOCATE_CONTIGUOUS 
//#define FFTTOOLS_THREAD_SAFE 
//


// You may choose to define FFTTOOLS_COMPAT_MODE
// to get FFTtools to behave like it did several years ago (where it was a class with static methods). 
// This might help with legacy code. 


// c++ libraries thingies
#include <map>
#include <algorithm>
#include <iostream>
#include <fstream>

// ROOT
#include "TGraph.h"
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"

#ifdef FFTTOOLS_COMPAT_MODE
#define FFTTOOLS_METHOD static
#else
#define FFTTOOLS_METHOD  
#include "FFTWindow.h" 
#endif 


//My includes
#include "FFTWComplex.h"
#include "RFSignal.h"

class TRandom; 
class TH2; 

// FFTW
#include <complex>

/*! \mainpage ROOT FFTW Wrapper
 *
 * \section intro_sec Introduction
 *
 * This is the somewhat sketchy documentation for the ROOT FFTW Wrapper library (libRootFftwWrapper).
 *
 * \section prereq_sec Prerequisites
 *
 *  -# <A HREF="http://root.cern.ch">ROOT</A>
 *  -# <A HREF="http://www.fftw.org/">FFTW 3 -- Fastest Fourier Transform in the West</a>
 * 
 * \section install_sec Installation
 * -# Checkout the code from github https://github.com/nichol77/libRootFftwWrapper
 * -# Define the ANITA_UTIL_INSTALL_DIR (or ARA_UTIL_INSTALL_DIR) to point to the location you want the library installed (the library files will end up in (ANITA_UTIL_INSTALL_DIR)/lib and the header files in (ANITA_UTIL_INSTALL_DIR)/include).
 * -# Do <PRE>make</PRE><PRE>make install</PRE>
 * \section manual_sec Manual
 * If you are averse to reading web pages (and who wouldn't be) you can download a <a href="manual/libRootFftwWrapper.pdf">pdf copy of the reference material</a> but be warned it won't be a thrilling read.
 */

/*!
  This is the namespace that holds all the various useful FFT related functions. All of the functions canbe called using something like FFTtools::doYourFunkyThing(arg1,arg2).
*/


#ifdef FFTTOOLS_COMPAT_MODE
class 
#else
namespace 
#endif 
FFTtools
{
    
#ifdef FFTTOOLS_COMPAT_MODE
  public: 
#endif
  
  //! Interpolation Routines that use ROOT::Math::Interpolator 
  /*!
    \param grIn A pointer to the input TGraph.
    \param deltaT The desired period (1/rate) of the interpolated waveform.
    \return A pointer to ther interpolated TGraph, it is the users responsibility to delete this after use.
  */
   FFTTOOLS_METHOD TGraph *getInterpolatedGraph(const TGraph *grIn, Double_t deltaT);

  //! Interpolation Routines that zeropads the FFT
  /*!
    \param grIn A pointer to the input TGraph.
    \param deltaT The desired period (1/rate) of the interpolated waveform.
    \return A pointer to ther interpolated TGraph, it is the users responsibility to delete this after use.
  */
   FFTTOOLS_METHOD TGraph *getInterpolatedGraphFreqDom(const TGraph *grIn, Double_t deltaT);

  //! Convolution
  /*!
    \param grA A pointer to the input TGraph A.
    \param grB A pointer to the input TGraph B
    \return A pointer to the convolution of A and B stored in a TGraph, it is the users responsibility to delete this after use.
  */
   FFTTOOLS_METHOD TGraph *getConvolution(const TGraph *grA, const TGraph *grB);
  //! Convolution
  /*!
    \param grA A pointer to the input RFSignal A.
    \param grB A pointer to the input RFSignal B
    \return A pointer to the convolution of A and B stored in a TGraph, it is the users responsibility to delete this after use.
  */
   FFTTOOLS_METHOD RFSignal *getConvolution(RFSignal *grA, RFSignal *grB);

  //! Returns the magnitude of a complex number.
  /*!
    \param theNum The complex number.
    \return The magnitude of the complex number.
  */
   FFTTOOLS_METHOD double getAbs(FFTWComplex &theNum);

  //! Computes an inverse FFT
  /*!
    \param length The length of the output array
    \param theInput The input array of complex numbers of <i>(length/2 +1)</i>
    \return An array of <i>length</i> real numbers
  */
   FFTTOOLS_METHOD double *doInvFFT(int length, const FFTWComplex *theInput);
  //! Computes an FFT of an array of real numbers.
  /*!
    \param length The length of the input array.
    \param theInput The input array of <i>length</i> real numbers.
    \return An array of <i>length/2 + 1</i> complex numbers
  */
   FFTTOOLS_METHOD FFTWComplex *doFFT(int length,const double *theInput);    

   

   /** Version of doFFT don't require copying of memory. If these are not aligned properly (i.e. allocated with fftw_malloc, memalign or equivalent), bad things might happen. */ 
   FFTTOOLS_METHOD void doFFT(int length, const double * properly_aligned_input, FFTWComplex * properly_aligned_output); 

   /** Version of doInvFFt don't require copying of memory. If these are not aligned properly (i.e. allocated with fftw_malloc, memalign or equivalent), bad things might happen. Note that the input may be clobbered in this case.*/ 
   FFTTOOLS_METHOD void doInvFFTClobber(int length, FFTWComplex * properly_aligned_input_that_will_likely_be_clobbered, double * properly_aligned_output); 

   /** Version of doInvFFt that only requires copying of input. The input is
    * copied to a properly-aligned temporary on the stack.  If the output is
    * not  aligned properly (i.e. allocated with fftw_malloc, memalign or
    * equivalent), bad things might happen. Things will probably be a bit
    * faster if the input is aligned properly, as long as memcpy dispatches
    * correctly. 
    * */ 

   FFTTOOLS_METHOD void doInvFFTNoClobber(int length, const FFTWComplex * properly_aligned_input, double * properly_aligned_output); 

  //! Converts inputMag (linear units) to time domain by means of hilbert transform assuming R_signal = 1/sqrt(2) (R_mag - i R^_mag);
  /*!
    \param inputMag TGraph containg the inputMagnitude
    \return TGraph of the time domain signal
  */
   FFTTOOLS_METHOD TGraph *convertMagnitudeToTimeDomain(const TGraph *inputMag);    

  //! Computes the correlation of two subsets of TGraphs
  /*!
    \param gr1 The first TGraph in the correlation.
    \param gr2 The second TGraph in the correlation.
    \param firstIndex The index of the first sample in the sub traces.
    \param lastIndex The index of the last sample in the sub traces.
    \return The correlation as an array of <i>lastIndex-firstIndex</i> real numbers.
  */
   FFTTOOLS_METHOD double *getCorrelation(const TGraph *gr1, const TGraph *gr2,int firstIndex,int lastIndex);
  //! Computes the correlation of two arrays
  /*!
    \param length The length of the arrays
    \param oldY1 The first array in the correlation.
    \param oldY2 The second array in the correlation.
    \return The correlation as an array of <i>length</i> real numbers.
  */
   FFTTOOLS_METHOD double *getCorrelation(int length,const float *oldY1, const float *oldY2);
  //! Computes the correlation of two arrays
  /*!
    \param length The length of the arrays
    \param oldY1 The first array in the correlation.
    \param oldY2 The second array in the correlation.
    \return The correlation as an array of <i>length</i> real numbers.
  */
   FFTTOOLS_METHOD double *getCorrelation(int length,const double *oldY1,const  double *oldY2);


  //! This is designed for when you want to average a number of graphs of the same thing together. It uses a correlation to find the deltaT between graphs and then shifts the graphs and coherently sums them. The return is the average of the input graphs
  /*!
    \param numGraphs Number of graphs to average
    \param thePtrPtr A pointer to a two-dimensional array of TGraphs <i>*grPtrArray[numGraphs]</i>.
    \return A pointer to a graph containing the averaged waveform.
  */
   FFTTOOLS_METHOD TGraph *correlateAndAverage(Int_t numGraphs, TGraph **grPtrPtr);

  //! This is designed for when you want to average a number of graphs of the same thing together. It uses a correlation to find the deltaT between graphs and then shifts the graphs and coherently sums them. The return is the average of the input graphs
  //! This is intended to be used for a set of multiple waveforms where we have found some event-wide peak bin from previous correlation
  //! Using this method assures that the bin associated with the maximum correlation value is found within a limited range of the predetermined event wide peak bin
  /*!
    \param numGraphs Number of graphs to average
    \param thePtrPtr A pointer to a two-dimensional array of TGraphs <i>*grPtrArray[numGraphs]</i>.
    \param correlationBinPtr A pointer to the pre-determined average correlation on an event-by-event basis
    \param binWiggleRoom How many bins either side of the bin associated with the peak value in the reduced range used to find the max correlation,
    with the default set to 20 bins
    \return A pointer to a graph containing the averaged waveform.
  */

  FFTTOOLS_METHOD TGraph *correlateAndAveragePredeterminedZone(Int_t numGraphs, TGraph **grPtrPtr, Int_t *correlationBinPtr, Int_t binWiggleRoom);
  
  //! This is designed for when you want to average a number of graphs of the same thing together. It uses a correlation to find the deltaT between graphs and then shifts the graphs and coherently sums them. The return is the average of the input graphs
  /*!
    \param deltaTInt The time-step to interpolate to before doing the correlation and averaging
    \param numGraphs Number of graphs to average
    \param thePtrPtr A pointer to a two-dimensional array of TGraphs <i>*grPtrArray[numGraphs]</i>.
    \return A pointer to a graph containing the averaged waveform.
  */
  
   FFTTOOLS_METHOD TGraph *interpolateCorrelateAndAverage(Double_t deltaTInt,
						Int_t numGraphs,
						TGraph **grPtrPtr);
  
  //! Returns the time domain result of a frequency domain sum of a number of arrays. As of writing this documentation, I'm not sure why this would be interesting.
  /*!
    \param numArrays The number of arrays to sum.
    \param thePtrPtr A pointer to a two-dimensional array of doubles <i>[numArrays][eachLength]</i>.
    \param eachLength The length of each array.
    \return The time domain result of the frequency domain summation.
  */
  FFTTOOLS_METHOD  Double_t *combineValuesUsingFFTs(Int_t numArrays, Double_t **thePtrPtr, Int_t eachLength);

  //Higher level functions that take and return TGraphs
  //! Returns the power spectral density. Note the PSD is unormalised (or if you prefer is normalised to the sum squared amplitude of the time domain). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a>
  /*!
    \param grWave The input time domain waveform
    \return A pointer to a TGraph containing the power spectrum. 
  */
  FFTTOOLS_METHOD  TGraph *makePowerSpectrum(const TGraph *grWave);
  //! Returns the power spectral density. Note the PSD returned is the periodogram (or if you prefer is normalised to the mean squared amplitude of the time domain). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a>
  /*!
    \param grWave The input time domain waveform
    \return A pointer to a TGraph containing the power spectrum. 
  */
  FFTTOOLS_METHOD  TGraph *makePowerSpectrumPeriodogram(const TGraph *grWave);
  //! Returns the power spectral density. Note the PSD returned is normalised and divided by frequency bin width (or if you prefer it is normalised to the time-integral squared amplitude of the time domain and then divided by frequency bin width). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a> As the name suggests this function expects the input waveform to be a volts-seconds one.
  /*!
    \param grWave The input time domain waveform with units of volts-seconds.
    \return A pointer to a TGraph containing the power spectrum. 
  */  
  FFTTOOLS_METHOD  TGraph *makePowerSpectrumVoltsSeconds(const TGraph *grWave);
  //! Returns the power spectral density of the input waveform convolved with a Bartlett Window. Note the PSD returned is normalised and divided by frequency bin width (or if you prefer it is normalised to the time-integral squared amplitude of the time domain and then divided by frequency bin width). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a> As the name suggests this function expects the input waveform to be a volts-seconds one.
  /*!
    \param grWave The input time domain waveform with units of volts-seconds.
    \return A pointer to a TGraph containing the power spectrum. 
  */  
  FFTTOOLS_METHOD  TGraph *makePowerSpectrumVoltsSecondsBartlett(const TGraph *grWave);
  //! Returns the power spectral density of the input waveform. In this one we first zero pad the waveform and then split it up into overlapping segments and convolve each segment with the Bartlett window before summing the resulting PSD's. As the name suggests this function expects the input waveform to be a volts-seconds one. No idea if this one actually works, or where I read oabout this crazy method which is supposed to reduce the variance of the PSD estimator.
  /*!
    \param grWave The input time domain waveform with units of volts-seconds.
    \param padFactor The factor by which to zero pad the wave see padWave
    \param numFreqs The number of frequency bins required in the output, which is related to the length of overlapping segments.
    \return A pointer to a TGraph containing the power spectrum. 
  */    
  FFTTOOLS_METHOD  TGraph *makePSVSBartlettPaddedOverlap(const TGraph *grWave, Int_t padFactor=4, Int_t numFreqs=64);

  //! Returns the power spectral density in dB units. Note the PSD returned is normalised and divided by frequency bin width (or if you prefer it is normalised to the time-integral squared amplitude of the time domain and then divided by frequency bin width). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a> As the name suggests this function expects the input waveform to be a volts-seconds one.
  /*!
    \param grWave The input time domain waveform with units of volts-seconds.
    \return A pointer to a TGraph containing the power spectrum in dB units, the frequency units are MHz. 
  */  
  FFTTOOLS_METHOD  TGraph *makePowerSpectrumVoltsSecondsdB(const TGraph *grWave);
 //! Returns the power spectral density in dB units. Note the PSD returned is normalised and divided by frequency bin width (or if you prefer it is normalised to the time-integral squared amplitude of the time domain and then divided by frequency bin width). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a> As the name suggests this function expects the input waveform to be a millivolts-nanoseconds one.
  /*!
    \param grWave The input time domain waveform with units of millivolts-nanoseconds.
    \return A pointer to a TGraph containing the power spectrum in dB units, the frequency units are MHz. 
  */  
  FFTTOOLS_METHOD  TGraph *makePowerSpectrumMilliVoltsNanoSeconds(const TGraph *grWave);
  FFTTOOLS_METHOD  TGraph *makePowerSpectrumMilliVoltsNanoSecondsdB(const TGraph *grWave);
  //! Returns the power spectral density of the input waveform zero-padded by some factor. Note the PSD returned is normalised and divided by frequency bin width (or if you prefer it is normalised to the time-integral squared amplitude of the time domain and then divided by frequency bin width). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a> As the name suggests this function expects the input waveform to be a volts-seconds one.
  /*!
    \param grWave The input time domain waveform with units of volts-seconds.
    \param padFactor The factor by which to zero pad the wave.
    \return A pointer to a TGraph containing the power spectrum. 
  */  
  FFTTOOLS_METHOD  TGraph *makePowerSpectrumVoltsSecondsPadded(const TGraph *grWave, Int_t padFactor=4);
  //! Returns the power spectral density in dB units of the input waveform zero-padded by some factor. Note the PSD returned is normalised and divided by frequency bin width (or if you prefer it is normalised to the time-integral squared amplitude of the time domain and then divided by frequency bin width). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a> As the name suggests this function expects the input waveform to be a volts-seconds one.
  /*!
    \param grWave The input time domain waveform with units of volts-seconds.
    \param padFactor The factor by which to zero pad the wave.
    \return A pointer to a TGraph containing the power spectrum in dB units with MHz as the frequency unit. 
  */    
  FFTTOOLS_METHOD  TGraph *makePowerSpectrumVoltsSecondsPaddeddB(const TGraph *grWave, Int_t padFactor=4);
  
  //! Returns the power spectral density in completely unormalised unit (as in Parseval's theorem is not obeyed and there is an extra factor of N not removed form the PSD). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a>
  /*!
    \param grWave The input time domain waveform.
    \return A pointer to a TGraph containing the power spectrum. 
  */    
  FFTTOOLS_METHOD  TGraph *makeRawPowerSpectrum(const TGraph *grWave);

   //! Returns the correlation of two TGraphs
  /*!
    \param gr1 The first input TGraph
    \param gr2 The second input TGraph
    \return A pointer to a TGraph containing the correlation of <i>gr1</i> and <i>gr2</i>.
  */    
  FFTTOOLS_METHOD  TGraph *getCorrelationGraph(const TGraph *gr1, const TGraph *gr2, Int_t *zeroOffset=0);



   //! Returns the normalised correlation of two TGraphs
  /*!
    \param gr1 The first input TGraph (must be zero meaned)
    \param gr2 The second input TGraph (must be zero meaned)
    \return A pointer to a TGraph containing the correlation of <i>gr1</i> and <i>gr2</i> where each point is normalised by the number of valid samples in the correlation and by the product of the RMS of the input graphs.
  */    
  FFTTOOLS_METHOD  TGraph *getNormalisedCorrelationGraph(const TGraph *gr1, const TGraph *gr2, Int_t *zeroOffset=0);

   //! Returns the normalised correlation of two TGraphs
  /*!
    \param gr1 The first input TGraph (must be zero meaned)
    \param gr2 The second input TGraph (must be zero meaned)
    \param zeroOffset A pointer to an integer where the sample corresponding to zero offset will be stored
    \param useDtRange A flag to enable the setting of a limited range of deltat's to try
    \param dtMin The minimum delta-t to include in the correlation, the maximum delta-t to include in the correlation
    \return A pointer to a TGraph containing the correlation of <i>gr1</i> and <i>gr2</i> where each point is normalised by the number of valid samples in the correlation and by the product of the RMS of the input graphs.
  */  
  FFTTOOLS_METHOD  TGraph *getNormalisedCorrelationGraphTimeDomain(const TGraph *gr1, const TGraph *gr2, Int_t *zeroOffset=0, Int_t useDtRange=0, Double_t dtMin=-1000, Double_t dtMax=1000);

   //! Returns the correlation of two interpolated TGraphs
  /*!
    \param grIn1 The first input TGraph 
    \param grIn2 The second input TGraph 
    \param deltaT The desired time step for the interpolated graphs
    \return A pointer to a TGraph containing the correlation of <i>gr1</i> and <i>gr2</i>, after each is interpolated to have a timestep of <i>deltaT</i>.
  */      
  FFTTOOLS_METHOD  TGraph *getInterpolatedCorrelationGraph(const TGraph *grIn1, const TGraph *grIn2, Double_t deltaT);
  //! Returns the inverse FFT of the FFT of the input TGraph. Seems pointless.
  /*!
    \param grWave A pointer to the input TGraph
    \return A pointer to a TGraph containing the inverse FFT of the FFT of <i>grWave</i>.
  */        
  FFTTOOLS_METHOD  TGraph *makeInverseInverseSpectrum(const TGraph *grWave);
  
   //! Returns the time domain result of a frequency domain sum of a number of TGraphs. In the sum each TGraph is weighted by a value. As of writing this documentation, I'm not sure why this would be interesting.
  /*!
    \param numGraphs The number of TGraphs to sum.
    \param grPtr A pointer to an array of <i>numGraphs</i> TGraph pointers.
    \param theWeights An optional array of weights with which to sum the TGraphs.
    \return A pointer to a TGraph containing the inverse FFT of the weighted summed FFT of the input graphs.
  */
  FFTTOOLS_METHOD  TGraph *combineGraphsUsingFFTs(Int_t numGraphs, TGraph **grPtr,const double *theWeights=0);
  
   //! Smooth graph using box car smoothing
  /*!
    \param grWave A pointer to the input TGraph
    \param halfWidth The halfWidth in samples of the size of the box over which to smooth the input graph (eg. <i>halfWidth</i>=1 means that sample i is the average if i-1, i and i+1.
    \return A pointer to a TGraph containing the smoothed graph.
  */
   FFTTOOLS_METHOD TGraph *getBoxCar(const TGraph *grWave, Int_t halfWidth);
  //! The Hilbert transform of the input TGraph
  /*!
    \param grWave A pointer to the input TGraph
    \return A pointer to a TGraph containing the Hilbert transform
  */
   FFTTOOLS_METHOD TGraph *getHilbertTransform(const TGraph *grWave);


   /** Hilbert transform on an array y is input,N is length*/ 
   FFTTOOLS_METHOD double *  getHilbertTransform(int N, const double * y) ; 


  //! The Hilbert envelope of the input TGraph. This is defined as e_i=sqrt(v_i^2 + h_i^2), where e_i, v_i and h_i are the i-th sample of the envelope, input graph and hilbert transform of the input graph repsectively.
  /*!
    \param grWave A pointer to the input TGraph
    \return A pointer to a TGraph containing the Hilbert envelope function.
  */
   FFTTOOLS_METHOD TGraph *getHilbertEnvelope(const TGraph *grWave);

  //Utility functions (not necessarily FFT related but they'll live here for now

  //! The linear sum of the power in a TGraph (normally a PSD)
  /*!
    \param gr A pointer to the input TGraph (normally a PSD)
    \param firstBin The first bin to include in the sum.
    \param lastBin The last bin to include in the sum.
    \return The summed power
  */
   FFTTOOLS_METHOD Double_t sumPower(const TGraph *gr,Int_t firstBin=-1,Int_t lastBin=-1);
   //! The integral of the power in a TGraph (normally a PSD) (i.e the sum of bin content*bin width)
  /*!
    \param gr A pointer to the input TGraph (normally a PSD)
    \param firstBin The first bin to include in the integral.
    \param lastBin The last bin to include in the integral.
    \return The integral value.
  */
  FFTTOOLS_METHOD  Double_t integratePower(const TGraph *gr,Int_t firstBin=-1,Int_t lastBin=-1);
   //! The sum of the voltage squared in a waveform.
  /*!
    \param gr A pointer to the input TGraph 
    \param firstBin The first bin to include in the sum.
    \param lastBin The last bin to include in the sum.
    \return The value of the sum.
  */  
  FFTTOOLS_METHOD  Double_t sumVoltageSquared(const TGraph *gr,Int_t firstBin,Int_t lastBin); 
  //! The integral of the v^2*dt in a waveform. Now works for unevenly sampled waeforms
  /*!
    \param gr A pointer to the input TGraph 
    \param firstBin The first bin to include in the integral.
    \param lastBin The last bin to include in the integral.
    \return The value of the integral.
  */  
  FFTTOOLS_METHOD  Double_t integrateVoltageSquared(const TGraph *gr,Int_t firstBin=-1,Int_t lastBin=-1);
  //! Find the peak (maximum positive) bin in a TGraph
  /*!
    \param gr A pointer to the input TGraph 
    \return The index of the bin with peak value.
  */    
  FFTTOOLS_METHOD  Int_t getPeakBin(const TGraph *gr); 
  //! Find the peak (maximum positive) bin in a TGraph
  /*!
    \param gr A pointer to the input TGraph 
    \param firstBin The first bin to include.
    \param lastBin The last bin to include.
    \return The index of the bin with peak value.
  */    
  FFTTOOLS_METHOD  Int_t getPeakBin(const TGraph *gr, Int_t firstBin, Int_t lastBin);  
    //! Find the x value associated with the peak (maximum positive) in a TGraph
  /*!
    \param gr A pointer to the input TGraph 
    \return The x value of the bin with peak value.
  */    
  FFTTOOLS_METHOD Double_t getPeakXvalue(const TGraph *gr); 
  //! Find the peak (maximum positive) in a TGraph.
  /*!
    \param gr A pointer to the input TGraph.
    \param index An optional pointer in which the peak bin can be stored.
    \return The peak value.
  */    
  FFTTOOLS_METHOD  Double_t getPeakVal(const TGraph *gr, int *index=0);
  //! Find the peak (maximum positive) in a TGraph.
  /*!
    \param gr A pointer to the input TGraph.
    \param firstBin The first bin to include.
    \param lastBin The last bin to include.
    \param index An optional pointer in which the peak bin can be stored.
    \return The peak value.
  */    
  FFTTOOLS_METHOD Double_t getPeakVal(const TGraph *gr, Int_t firstBin, Int_t lastBin, int *index=0);
  //! Find the peak (v^2) in a TGraph.
  /*!
    \param gr A pointer to the input TGraph.
    \param index An optional pointer in which the peak bin can be stored.
    \return The peak (v^2)value.
  */
  
  FFTTOOLS_METHOD  Double_t getPeakSqVal(const TGraph *gr, int *index=0);
  //! Find the peak (v^2) and RMS (of v^2) in a TGraph.
  /*!
    \param gr A pointer to the input TGraph.
    \param peak A reference to a Double_t where the peak value will be stored.
    \param rms A reference to a Double_t where the rms value will be stored.
    \param index An optional pointer in which the peak bin can be stored.
  */      
  FFTTOOLS_METHOD  void getPeakRmsSqVal(const TGraph *gr, Double_t &peak, Double_t &rms, Int_t *index=0);
  
  //! Find the peak (v) and RMS (of v) of a rectified TGraph.
  /*!
    \param gr A pointer to the input TGraph.
    \param peak A reference to a Double_t where the peak value will be stored.
    \param rms A reference to a Double_t where the rms value will be stored.
    \param index An optional pointer in which the peak bin can be stored.
  */      
  FFTTOOLS_METHOD  void getPeakRmsRectified(const TGraph *gr, Double_t &peak, Double_t &rms, Int_t *index=0);

  //Graph returning utility funcs
  //! Returns the simple power envelope of a waveform. The waveform is first squared and then only local peaks are taken to form the power envelope.
  /*!
    \param gr A pointer to the input TGraph.
    \return A pointer to a TGraph containing the power envelope
  */     
  FFTTOOLS_METHOD  TGraph *getSimplePowerEnvelopeGraph(const TGraph *gr);
  //! Returns a smoothed FFT where N-bins are averaged to reduce variance.
  /*!
    \param gr A pointer to the input TGraph.
    \param factor The factor by which to smooth (eg. 2 returns a PSD with half as many frequency bins).
    \return A pointer to a TGraph containing the smoothed PSD.
  */     
  FFTTOOLS_METHOD  TGraph *smoothFFT(const TGraph *gr,Int_t factor) ;
  //! Returns the difference between two graphs (A-B).
  /*!
    \param grA A pointer to the first graph.
    \param grB A pointer to the second graph.
    \return A pointer to a TGraph containing the difference (A-B) between the input graphs.
  */     
  FFTTOOLS_METHOD  TGraph *subtractGraphs(const TGraph *grA, const TGraph *grB);
  //! Returns the ratio between two graphs (A/B).
  /*!
    \param grA A pointer to the first graph.
    \param grB A pointer to the second graph.
    \return A pointer to a TGraph containing the ratio (A/B) of the input graphs.
  */     
   FFTTOOLS_METHOD TGraph *divideGraphs(const TGraph *grA, const TGraph *grB);

  
  //! Returns a graph translated by deltaT. Such that t'=t+dt
  /*!
    \param grWave A pointer to the input graph.
    \param deltaT The amount to move time axis by
    \return A pointer to a TGraph containing the translated graph
  */     
  FFTTOOLS_METHOD  TGraph *translateGraph(const TGraph *grWave, const Double_t deltaT);

  //! Returns the one minus the ratio between two graphs (A/B).
  /*!
    \param grA A pointer to the first graph.
    \param grB A pointer to the second graph.
    \return A pointer to a TGraph containing the one minus the ratio (1 - A/B) of the input graphs.
  */     
  FFTTOOLS_METHOD  TGraph *ratioSubtractOneGraphs(const TGraph *grA, const TGraph *grB) ;
  //! Returns the ratio of two graphs as 10 *log(A/B)
  /*!
    \param grA A pointer to the first graph.
    \param grB A pointer to the second graph.
    \return A pointer to a TGraph containing 10 * log(A/B)
  */     
  FFTTOOLS_METHOD  TGraph *dbGraphs(const TGraph *grA, const TGraph *grB);
  //! Zero pad a wave making it a factor of N longer.
  /*!
    \param grA A pointer to the input graph.
    \param padFactor The factor by whcih to increas the length (eg. new length = <i>factor</i> * old length)
    \return A pointer to the zero padded TGraph.
  */     
  FFTTOOLS_METHOD  TGraph *padWave(const TGraph *grA, Int_t padFactor);
  //! Zero pad a wave making it up to newLength points.
  /*!
    \param grA A pointer to the input graph.
    \param newLength The new length the graph should be (it gets zero padded on either side)
    \return A pointer to the zero padded TGraph.
  */     
  FFTTOOLS_METHOD  TGraph *padWaveToLength(const TGraph *grA, Int_t newLength);
  //! Rectify a waveform, optionally returning an all negative waveform.
  /*!
    \param gr A pointer to the input graph.
    \param makeNeg An optional parameter which if true will return a negative waveform.
    \return A pointer to the rectified TGraph
  */     
  FFTTOOLS_METHOD  TGraph *rectifyWave(const TGraph *gr, Int_t makeNeg=0);
  
  //Window functions
  //! The Bartlett window function (it's basically a triangle that peaks in the middle at 1 and goes to zero at the end points)
  /*!
    \param j The index to evaluate the window at.
    \param n The length of sample.
    \return The value of the window function at index <i>j</i>.
  */       
  FFTTOOLS_METHOD  Double_t  bartlettWindow(Int_t j, Int_t n);
  //! The Welch window function similar to the Bartlett window except that it falls off from the middle as dx^2 (ie. less steeply).
  /*!
    \param j The index to evaluate the window at
    \param n The length of sample.
    \return The value of the window function at index <i>j</i>.
  */       
  FFTTOOLS_METHOD  Double_t welchWindow(Int_t j, Int_t n);



  //! This fucntions just calculates the simple bin by bin dv/dt derivative of the input data.
  /*!
    \param numPoints The number of points in the input data (the output will be numPoints-1 in length).
    \param inputX The input time values.
    \param inputY The input voltage values.
    \param outputX The output time values.
    \param outputY The output voltage values.
  */       
  FFTTOOLS_METHOD  void takeDerivative(Int_t numPoints, const Double_t *inputX, const Double_t *inputY, Double_t *outputX, Double_t *outputY);

 //! This returns a TGraph which is the derivative of the input graph
  /*!
    \param grIn The input graph.
    \return The derviative of grIn.
  */      
  FFTTOOLS_METHOD  TGraph *getDerviative(const TGraph *grIn);
 //! This returns a TGraph which is the derivative of the input graph
  /*!
    \param grIn The input graph.
    \return The derviative of grIn.
  */      
  FFTTOOLS_METHOD  TGraph *getDerivative(const TGraph *grIn);

 //! This returns a TGraph which has had a simple pass band filter applied
  /*!
    \param grWave The input graph.
    \param minFreq The lowest frequency to pass.
    \param maxFreq The highest frequency to pass.
    \return The derviative of grWave.
  */      
  FFTTOOLS_METHOD  TGraph *simplePassBandFilter(const TGraph *grWave, Double_t minFreq, Double_t maxFreq);

 //! This returns a TGraph which has had a simple notch band filter applied
  /*!
    \param grWave The input graph.
    \param minFreq The lower frequency of the notch.
    \param maxFreq The upper frequency of the notch.
    \return The derviative of grWave.
  */      
  FFTTOOLS_METHOD  TGraph *simpleNotchFilter(const TGraph *grWave, Double_t minFreq, Double_t maxFreq);


  //! This returns a TGraph which has had N simple notch band filters applied
  /*!
    \param grWave The input graph.
    \param numNotches The number of notch regiosn
    \param minFreq An array of lower frequency of the notches.
    \param maxFreq An array upper frequency of the notch.
    \return The derviative of grWave.
  */      
  FFTTOOLS_METHOD  TGraph *multipleSimpleNotchFilters(const TGraph *grWave, Int_t numNotches, const Double_t minFreq[], const Double_t maxFreq[]);

//! This returns a TGraph which has been cropped in.
  /*!
    \param grWave The input graph.
    \param minTime The lower edge of the time window.
    \param maxTime The upper edge of the time window.
    \return The derviative of grWave.
  */      
  FFTTOOLS_METHOD  TGraph *cropWave(const TGraph *grWave, Double_t minTime, Double_t maxTime);

  
  //! This returns the SNR ratio of the input waveform
  /*!
    \param gr The input graph.
    \return The SNR of the input waveform, where S is half of the peak-to-peak and N is the RMS of the first 25 ns.
  */      
  FFTTOOLS_METHOD  Double_t getWaveformSNR(const TGraph *gr); 
  //! This returns the SNR ratio of the input waveform
  /*!
    \param gr The input graph.
    \param peakToPeak A reference to a double which will be set to half of the peak-to-peak
    \param rms A reference to a double which will be set to the RMS of teh first 25 samples
    \return The SNR of the input waveform, where S is half of the peak-to-peak and N is the RMS of the first 25 samples.
  */     
  FFTTOOLS_METHOD  Double_t getWaveformSNR(const TGraph *gr,Double_t &peakToPeak,Double_t &rms);

  //! This returns the largest (i.e most positive, or least negative) value
  /*!
    \param gr The input graph.
    \return The peak
  */      
  FFTTOOLS_METHOD  Double_t getWaveformPeak(const TGraph *gr);
  FFTTOOLS_METHOD Double_t getEnvelopeSNR(const TGraph *gr);
  FFTTOOLS_METHOD Double_t getEnvelopeSNR(const TGraph *gr,Double_t &peakToPeak,Double_t &rms,Double_t &timeOfPeak);


  //! Linear interpolate to find the value at some point near two known points
  /*!
    \param x1 x value of point 1
    \param y1 y value of point 1
    \param x2 x value of point 2
    \param y2 y value of point 2
    \param x is the point of interest
    \return y
  */      
   FFTTOOLS_METHOD Double_t simpleInterploate(Double_t x1, Double_t y1, Double_t x2, Double_t y2,Double_t x);
  
   /*!
    * Correlation between two FFTs
    * bandpasses between min_i and max_i using butterworth filter of order order
    * work can be used for temporary to avoid allocation of new memory 
    */
   FFTTOOLS_METHOD double * FFTCorrelation(int waveformlength, const FFTWComplex * A, const FFTWComplex * B, FFTWComplex * work = 0, 
                           int min_i = 0, int max_i =0, int order=1);  


   /*! in place array rotation
    */
   FFTTOOLS_METHOD void inPlaceShift(int N, double *x); 


   /*! inplace graph rotation */ 
   FFTTOOLS_METHOD void rotate(TGraph * g, int rot); 

   /*! does fit to polynomial of order order then subtracts it */
   FFTTOOLS_METHOD void polySubtract(TGraph *g, int order=1); 


      
   /*! direct convolution for kernel can either assume zero's outside boundaries or repetition of first and last values */
   enum DirectConvolveEdgeBehavior
   {
        ZEROES_OUTSIDE, 
        REPEAT_OUTSIDE
   }; 


    /*!convolution without FFT of x with kernel h. y should have same size as x (it is cropped symmetrically) 
      if delay = 0, h will be centered around its middle sample, if - M/2 , will be purely causal, etc. 
    */
    FFTTOOLS_METHOD double * directConvolve(int N, const double *x, int M, const double * h, double *y = 0, int delay = 0,  DirectConvolveEdgeBehavior edge_behavior = ZEROES_OUTSIDE); 


    /*! wraps periodic array of doubles. center assumed to be period/2 */ 
   FFTTOOLS_METHOD void wrap(size_t N, double * vals, double period = 360); 

    /*! wraps periodic array of doubles. */ 
   FFTTOOLS_METHOD void wrap(size_t N, double * vals, double period, double center); 

    /*! wraps periodic array of floats */ 
   FFTTOOLS_METHOD void wrap(size_t N, float * vals, float period, float center); 

    /*! wraps periodic array of floats. center assumed to be period/2 */ 
   FFTTOOLS_METHOD void wrap(size_t N, float * vals, float period = 360); 


   /*! Faster than the libm floor, which cares about errno and the floating point environment and all that crap 
    **/ 
   FFTTOOLS_METHOD inline int fast_floor(double val) { return (int) val - (val < (int) val); }

   /*!Wraps a value to be within period centered at center.
    * This could be implemented in terms of the other wraps, but in practice it's better to inline since
    * the period and center are usually known at compile them. 
    *
    * */ 
   FFTTOOLS_METHOD inline double wrap(double val, double period, double center)
          { return val - period * fast_floor((val-center+period/2)/period); }

   /*! Wraps a value to be within period. Center assumed to be period/2 */ 
   FFTTOOLS_METHOD inline double wrap(double val, double period = 360)
          { return wrap(val,period,period/2); } 


    /*! unwraps periodic array of doubles. */ 
   FFTTOOLS_METHOD void unwrap(size_t N, double * vals, double period = 360); 


    /*! unwraps periodic array of floats. */ 
   FFTTOOLS_METHOD void unwrap(size_t N, float * vals, float period = 360); 


   /*! computes reasonable dt from unevenly sampled graph */ 
   FFTTOOLS_METHOD double getDt(const TGraph * g, int realN = 0);  

   /*! linearly interpolates value x in g . Similar to doing g->Eval() but much faster (since doesn't need to sort or binary search)*/ 
   FFTTOOLS_METHOD double evalEvenGraph(const TGraph * g, double x); 



#ifndef FFTTOOLS_COMPAT_MODE
   /*! applies window to graph */ 
   void applyWindow(TGraph *g, const FFTWindowType *w); 
#endif

   
   /*! random variable from rayleigh distribution */ 
   FFTTOOLS_METHOD double randomRayleigh(double sigma=1, TRandom * rng = 0); 

   /*! sinc function, if |x| < eps, sinc(x) = 1 */ 
   FFTTOOLS_METHOD double sinc(double x, double eps = 0); 
 
   FFTTOOLS_METHOD int loadWisdom(const char * file); 
   FFTTOOLS_METHOD int saveWisdom(const char * file); 

#ifndef FFTTOOLS_COMPAT_MODE
   /**
    *  Make a welch periodogram of evenly sampled input. A Bartlett periodogram can be emulated using overlap_fraction = 0 and window = &RECTANGULAR_WINDOW
    *
    *  @param gin evenly sampled graph to estimate power spectrum of
    *  @param segment_siz ethe number of samples in each segment
    *  @param overlap_fraction the fraction of overlap for each segment. 0 is no overlap. 1 is non-sensical full overlap which will probably cause some calamity. 
    *  @param window the window to use for each segment
    *  @param truncate_extra If true, any partial segment is discarded. Otherwise, it will be zero-padded 
    *  @param gout if non-zero, this TGraph will be used for output
    *  @return the Welch Periodogram 
    *
    */ 
   TGraph * welchPeriodogram(const TGraph * gin, int segment_size, double overlap_fraction = 0.5, const FFTWindowType * window = &GAUSSIAN_WINDOW , bool truncate_extra = true, TGraph * gout = 0); 



   
   /** "normal" lomb-scargle periodogram, not using the N log N algorithm in Press & Rybicki) */
   FFTTOOLS_METHOD double *  lombScarglePeriodogramSlow(int N, const double *x, const double * y, int nfreqs, const double * freqs, double * answer = 0); 

   /** fast periodogoram (as in Press & Rybicki) . Implementation in Periodogram.cxx */
   FFTTOOLS_METHOD TGraph * lombScarglePeriodogram(const TGraph * g, double dt = 0, double oversample_factor  = 4 , 
                       double high_factor = 2, TGraph * replaceme = 0, int extirpolation_factor =4)  ; 
   FFTTOOLS_METHOD TGraph * lombScarglePeriodogram(int N, double dt, const double * __restrict x, 
                                   const double * __restrict y, double oversample_factor  = 4 , 
                                   double high_factor = 2, TGraph * replaceme = 0, int extirpolation_factor = 4)  ; 


#endif

//   TH2 * getPowerVsTimeUsingLombScargle(const TGraph * g, int nbins, double sampling_dt = 0, double oversample_factor = 4, double high_factor = 2, TH2 * useme = 0); 

   /* Compute Stokes parameters from hpol / vpol and their hilbert transforms 
    *
    * This computes either the average, as well as optionally the instantaneous values.  
    *
    *
    * @param N number of samples
    * @param hpol hpol waveform (evenly sampled)
    * @param hpol_hat hilbert transform of hpol waveform (evenly sampled)
    * @param vpol vpol waveform (evenly sampled)
    * @param vpol_hat hilbert transform of vpol waveform (evenly sampled)
    * @param Iavg pointer to where Stokes I integral average will be stored (or null to not store). 
    * @param Qavg pointer to where Stokes Q integral average will be stored (or null to not store). 
    * @param Uavg pointer to where Stokes U integral average will be stored (or null to not store).
    * @param Vavg pointer to where Stokes V integral average will be stored (or null to not store). 
    * @param Iins pointer to where Stokes I instantaneous value will be stored (or null to not store). Must be of length N unless only_store_max_instantaneous is true 
    * @param Qins pointer to where Stokes Q instantaneous value will be stored (or null to not store).  Must be of length N unless only_store_max_instantaneous is true
    * @param Uins pointer to where Stokes U instantaneous value will be stored (or null to not store). Must be of length N unless only_store_max_instantaneous is true
    * @param Vins pointer to where Stokes V instantaneous value will be stored (or null to not store).  Must be of length N unless only_store_max_instantaneous is true
    * @param only_max_instantaneous  if true,  only returns the maximum of the instantaneous quantities (at the maximum I in the range) 
    */ 
   FFTTOOLS_METHOD void stokesParameters(int N, const double *  hpol, const double *  hpol_hat, const double *  vpol, const double *  vpol_hat, 
                         double *  Iavg = 0, double *  Qavg = 0, double *  Uavg = 0, double *  Vavg = 0, 
                         double *  Iins = 0, double *  Qins = 0, double *  Uins = 0, double *  Vins = 0 , 
                         bool only_max_instantaneous = true ); 



   /** Compute the DFT term  at a given frequency. Probably not want you want to do for a series of frequencies where one can 
    *  take advantage of trig relations, etc. 
    *
    *  @param g What to  look at
    *  @param freq the frequency
    *  @param phase pointer where phase will be stored or 0 to not store anything
    *  @param amp pointer where amplitude will be stored or 0 to not store anything
    *  @param real pointer where real will be stored or 0 to not store anything
    *  @param imag pointer where imag will be stored or 0 to not store anything
    */ 

   FFTTOOLS_METHOD void dftAtFreq(const TGraph * g, double freq, double * phase, double * amp = 0, double * real = 0, double * imag = 0); 


   /** Compute the DFT term at a given frequency and n multiples of it. Timebase does not have to be even. Phase and amp should be big enough to store nmultiples or 0 if you don't want to fill one. 
    * @param g waveform to look at 
    * @param freq base frequency
    * @param nmultiples number of multiples to consider
    * @param phase pointer to array where phase will be stored, or 0 to not store anything. 
    * @param amp pointer to array where amp will be stored, or 0 to not store anything. 
    * @param real pointer where real will be stored or 0 to not store anything
    * @param imag pointer where imag will be stored or 0 to not store anything
    *
    * */ 
   FFTTOOLS_METHOD void dftAtFreqAndMultiples(const TGraph * g, double freq, int nmultiples, double * phase, double * amp = 0, double * real = 0, double * imag = 0); 


   /** Given a desired amplitude  response G, makes a minimum phase version of the waveform. This takes advantage
    * of the fact that for a minimum phase signal, the phase response equals the hilbert transform of the log magnitude. 
    * @param N, the length of the amplitude response (typically the time domain length/ 2 + 1), 
    * @param F, the amplitude response G(w) with w in normalized units (first component is dc, last is nyquist). This is in linear units. 
    * @param mindb, an absolute zero in the frequency response cannot be represented with a minimum phase filter. Anything with a (power) magnitude smaller than mindb will be made equal to mindb 
    * @returns A minimum phase version of the waveform. 
    */ 
   FFTTOOLS_METHOD FFTWComplex * makeMinimumPhase(int N, const double * G,  double mindb = -100); 


   /** Compute a minimum phase version of a time-domain signal.
    * 
    * @param g the time-domain signal, probably needs to be zero-padded
    * @param mindb, an absolute zero in the frequency response cannot be represented with a minimum phase filter. Anything with a (power) magnitude smaller than mindb will be made equal to mindb 
    * @returns minimum phase version of g 
    */
   FFTTOOLS_METHOD TGraph * makeMinimumPhase(const TGraph *g, double mindb=-100); 


   /** Compute a fixed group delay version of a time-domain signal.
    * 
    * @param g the time-domain signal, probably needs to be zero-padded
    * @param phase the phase (in radians) of the output
    * @returns fixed phase version of g 
    */
   FFTTOOLS_METHOD TGraph * makeFixedDelay(const TGraph *g, double delay = 1); 


   /** Checks if a signal is causal or not by comparing the imaginary part to the hilbert transform of the real part and returning the square difference (scaled by N). 
    * A perfectly causal system will return 0, but numerical inaccuracy may make that not true. 
    * */ 
   FFTTOOLS_METHOD double checkCausality(int N, const FFTWComplex * signal) ; 


   /** Avg RMS envelope 
    *  Approximates envelope by computing rms in each window 
    *
    *  @param N  length of input (and output) arrays
    *  @param W  window length, in units of the x values (not all windows may contain the same length) 
    *  @param x  input x values
    *  @param y  input y values
    *  @param out output y values, if 0 output will be allocated
    *  @return output 
    */ 

   FFTTOOLS_METHOD double * rmsEnvelope(int N, double W, const double * x, const double * y, double * out = 0); 

   /** Peak envelope 
    * Approximates envelope by connecting peaks of |y|. 
    *
    * @param N the number of points in input / output arrays
    * @param min_distance  the minimum distance between peaks. This should be close to the half-period of the carrier wave you expect. 
    * @param x input x values (times) 
    * @param y input y values (voltages) 
    * @param out output y values (voltages) , or 0 to allocate new
    * @return output values 
    */
   FFTTOOLS_METHOD double * peakEnvelope(int N, double min_distance, const double * x, const double *y, double * out = 0); 
#ifdef FFTTOOLS_COMPAT_MODE
};
#else
}
#endif 
   
#endif //FFTTOOLS_H
