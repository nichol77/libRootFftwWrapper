
#ifndef FFTTOOLS_H
#define FFTTOOLS_H



// If the following is uncommented, in and out arrays for FFTW are allocated contiguously. 
#define FFTTOOLS_ALLOCATE_CONTIGUOUS 
//#define FFTTOOLS_THREAD_SAFE 

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
#include "FFTWindow.h" 


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
namespace FFTtools
{
    
  
  //! Interpolation Routines that use ROOT::Math::Interpolator 
  /*!
    \param grIn A pointer to the input TGraph.
    \param deltaT The desired period (1/rate) of the interpolated waveform.
    \return A pointer to ther interpolated TGraph, it is the users responsibility to delete this after use.
  */
   TGraph *getInterpolatedGraph(TGraph *grIn, Double_t deltaT);

  //! Interpolation Routines that zeropads the FFT
  /*!
    \param grIn A pointer to the input TGraph.
    \param deltaT The desired period (1/rate) of the interpolated waveform.
    \return A pointer to ther interpolated TGraph, it is the users responsibility to delete this after use.
  */
   TGraph *getInterpolatedGraphFreqDom(TGraph *grIn, Double_t deltaT);

  //! Convolution
  /*!
    \param grA A pointer to the input TGraph A.
    \param grB A pointer to the input TGraph B
    \return A pointer to the convolution of A and B stored in a TGraph, it is the users responsibility to delete this after use.
  */
   TGraph *getConvolution(TGraph *grA, TGraph *grB);
  //! Convolution
  /*!
    \param grA A pointer to the input RFSignal A.
    \param grB A pointer to the input RFSignal B
    \return A pointer to the convolution of A and B stored in a TGraph, it is the users responsibility to delete this after use.
  */
   RFSignal *getConvolution(RFSignal *grA, RFSignal *grB);

  //! Returns the magnitude of a complex number.
  /*!
    \param theNum The complex number.
    \return The magnitude of the complex number.
  */
   double getAbs(FFTWComplex &theNum);

  //! Computes an inverse FFT
  /*!
    \param length The length of the output array
    \param theInput The input array of complex numbers of <i>(length/2 +1)</i>
    \return An array of <i>length</i> real numbers
  */
   double *doInvFFT(int length, FFTWComplex *theInput);
  //! Computes an FFT of an array of real numbers.
  /*!
    \param length The length of the input array.
    \param theInput The input array of <i>length</i> real numbers.
    \return An array of <i>length/2 + 1</i> complex numbers
  */
   FFTWComplex *doFFT(int length,double *theInput);    

   
   /** Version of doFFT don't require copying of memory. If these are not aligned properly (i.e. allocated with fftw_malloc, memalign or equivalent), bad things might happen. */ 
   void doFFT(int length, const double * properly_aligned_input, FFTWComplex * properly_aligned_output); 

   /** Version of doInvFFt don't require copying of memory. If these are not aligned properly (i.e. allocated with fftw_malloc, memalign or equivalent), bad things might happen. Note that the input may be clobbered in this case.*/ 
   void doInvFFTClobber(int length, FFTWComplex * properly_aligned_input_that_will_likely_be_clobbered, double * properly_aligned_output); 

   /** Version of doInvFFt that only requires copying of input. The input is
    * copied to a properly-aligned temporary on the stack.  If the output is
    * not  aligned properly (i.e. allocated with fftw_malloc, memalign or
    * equivalent), bad things might happen. Things will probably be a bit
    * faster if the input is aligned properly, as long as memcpy dispatches
    * correctly. 
    * */ 

   void doInvFFTNoClobber(int length, const FFTWComplex * properly_aligned_input, double * properly_aligned_output); 

  //! Converts inputMag (linear units) to time domain by means of hilbert transform assuming R_signal = 1/sqrt(2) (R_mag - i R^_mag);
  /*!
    \param inputMag TGraph containg the inputMagnitude
    \return TGraph of the time domain signal
  */
   TGraph *convertMagnitudeToTimeDomain(TGraph *inputMag);    

  //! Computes the correlation of two subsets of TGraphs
  /*!
    \param gr1 The first TGraph in the correlation.
    \param gr2 The second TGraph in the correlation.
    \param firstIndex The index of the first sample in the sub traces.
    \param lastIndex The index of the last sample in the sub traces.
    \return The correlation as an array of <i>lastIndex-firstIndex</i> real numbers.
  */
   double *getCorrelation(TGraph *gr1, TGraph *gr2,int firstIndex,int lastIndex);
  //! Computes the correlation of two arrays
  /*!
    \param length The length of the arrays
    \param oldY1 The first array in the correlation.
    \param oldY2 The second array in the correlation.
    \return The correlation as an array of <i>length</i> real numbers.
  */
   double *getCorrelation(int length,float *oldY1, float *oldY2);
  //! Computes the correlation of two arrays
  /*!
    \param length The length of the arrays
    \param oldY1 The first array in the correlation.
    \param oldY2 The second array in the correlation.
    \return The correlation as an array of <i>length</i> real numbers.
  */
   double *getCorrelation(int length,double *oldY1, double *oldY2);


  //! This is designed for when you want to average a number of graphs of the same thing together. It uses a correlation to find the deltaT between graphs and then shifts the graphs and coherently sums them. The return is the average of the input graphs
  /*!
    \param numGraphs Number of graphs to avergae
    \param thePtrPtr A pointer to a two-dimensional array of TGraphs <i>*grPtrArray[numGraphs]</i>.
    \return A pointer to a graph containing the averaged waveform.
  */
   TGraph *correlateAndAverage(Int_t numGraphs, TGraph **grPtrPtr);

  //! This is designed for when you want to average a number of graphs of the same thing together. It uses a correlation to find the deltaT between graphs and then shifts the graphs and coherently sums them. The return is the average of the input graphs
  /*!
    \param deltaTInt The time-step to interpolate to before doing the correlation and averaging
    \param numGraphs Number of graphs to avergae
    \param thePtrPtr A pointer to a two-dimensional array of TGraphs <i>*grPtrArray[numGraphs]</i>.
    \return A pointer to a graph containing the averaged waveform.
  */
  
   TGraph *interpolateCorrelateAndAverage(Double_t deltaTInt,
						Int_t numGraphs,
						TGraph **grPtrPtr);
  
  //! Returns the time domain result of a frequency domain sum of a number of arrays. As of writing this documentation, I'm not sure why this would be interesting.
  /*!
    \param numArrays The number of arrays to sum.
    \param thePtrPtr A pointer to a two-dimensional array of doubles <i>[numArrays][eachLength]</i>.
    \param eachLength The length of each array.
    \return The time domain result of the frequency domain summation.
  */
   Double_t *combineValuesUsingFFTs(Int_t numArrays, Double_t **thePtrPtr, Int_t eachLength);

  //Higher level functions that take and return TGraphs
  //! Returns the power spectral density. Note the PSD is unormalised (or if you prefer is normalised to the sum squared amplitude of the time domain). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a>
  /*!
    \param grWave The input time domain waveform
    \return A pointer to a TGraph containing the power spectrum. 
  */
   TGraph *makePowerSpectrum(TGraph *grWave);
  //! Returns the power spectral density. Note the PSD returned is the periodogram (or if you prefer is normalised to the mean squared amplitude of the time domain). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a>
  /*!
    \param grWave The input time domain waveform
    \return A pointer to a TGraph containing the power spectrum. 
  */
   TGraph *makePowerSpectrumPeriodogram(TGraph *grWave);
  //! Returns the power spectral density. Note the PSD returned is normalised and divided by frequency bin width (or if you prefer it is normalised to the time-integral squared amplitude of the time domain and then divided by frequency bin width). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a> As the name suggests this function expects the input waveform to be a volts-seconds one.
  /*!
    \param grWave The input time domain waveform with units of volts-seconds.
    \return A pointer to a TGraph containing the power spectrum. 
  */  
   TGraph *makePowerSpectrumVoltsSeconds(TGraph *grWave);
  //! Returns the power spectral density of the input waveform convolved with a Bartlett Window. Note the PSD returned is normalised and divided by frequency bin width (or if you prefer it is normalised to the time-integral squared amplitude of the time domain and then divided by frequency bin width). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a> As the name suggests this function expects the input waveform to be a volts-seconds one.
  /*!
    \param grWave The input time domain waveform with units of volts-seconds.
    \return A pointer to a TGraph containing the power spectrum. 
  */  
   TGraph *makePowerSpectrumVoltsSecondsBartlett(TGraph *grWave);
  //! Returns the power spectral density of the input waveform. In this one we first zero pad the waveform and then split it up into overlapping segments and convolve each segment with the Bartlett window before summing the resulting PSD's. As the name suggests this function expects the input waveform to be a volts-seconds one. No idea if this one actually works, or where I read oabout this crazy method which is supposed to reduce the variance of the PSD estimator.
  /*!
    \param grWave The input time domain waveform with units of volts-seconds.
    \param padFactor The factor by which to zero pad the wave see padWave
    \param numFreqs The number of frequency bins required in the output, which is related to the length of overlapping segments.
    \return A pointer to a TGraph containing the power spectrum. 
  */    
   TGraph *makePSVSBartlettPaddedOverlap(TGraph *grWave, Int_t padFactor=4, Int_t numFreqs=64);

  //! Returns the power spectral density in dB units. Note the PSD returned is normalised and divided by frequency bin width (or if you prefer it is normalised to the time-integral squared amplitude of the time domain and then divided by frequency bin width). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a> As the name suggests this function expects the input waveform to be a volts-seconds one.
  /*!
    \param grWave The input time domain waveform with units of volts-seconds.
    \return A pointer to a TGraph containing the power spectrum in dB units, the frequency units are MHz. 
  */  
   TGraph *makePowerSpectrumVoltsSecondsdB(TGraph *grWave);
 //! Returns the power spectral density in dB units. Note the PSD returned is normalised and divided by frequency bin width (or if you prefer it is normalised to the time-integral squared amplitude of the time domain and then divided by frequency bin width). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a> As the name suggests this function expects the input waveform to be a millivolts-nanoseconds one.
  /*!
    \param grWave The input time domain waveform with units of millivolts-nanoseconds.
    \return A pointer to a TGraph containing the power spectrum in dB units, the frequency units are MHz. 
  */  
   TGraph *makePowerSpectrumMilliVoltsNanoSeconds(TGraph *grWave);
   TGraph *makePowerSpectrumMilliVoltsNanoSecondsdB(TGraph *grWave);
  //! Returns the power spectral density of the input waveform zero-padded by some factor. Note the PSD returned is normalised and divided by frequency bin width (or if you prefer it is normalised to the time-integral squared amplitude of the time domain and then divided by frequency bin width). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a> As the name suggests this function expects the input waveform to be a volts-seconds one.
  /*!
    \param grWave The input time domain waveform with units of volts-seconds.
    \param padFactor The factor by which to zero pad the wave.
    \return A pointer to a TGraph containing the power spectrum. 
  */  
   TGraph *makePowerSpectrumVoltsSecondsPadded(TGraph *grWave, Int_t padFactor=4);
  //! Returns the power spectral density in dB units of the input waveform zero-padded by some factor. Note the PSD returned is normalised and divided by frequency bin width (or if you prefer it is normalised to the time-integral squared amplitude of the time domain and then divided by frequency bin width). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a> As the name suggests this function expects the input waveform to be a volts-seconds one.
  /*!
    \param grWave The input time domain waveform with units of volts-seconds.
    \param padFactor The factor by which to zero pad the wave.
    \return A pointer to a TGraph containing the power spectrum in dB units with MHz as the frequency unit. 
  */    
   TGraph *makePowerSpectrumVoltsSecondsPaddeddB(TGraph *grWave, Int_t padFactor=4);
  
  //! Returns the power spectral density in completely unormalised unit (as in Parseval's theorem is not obeyed and there is an extra factor of N not removed form the PSD). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a>
  /*!
    \param grWave The input time domain waveform.
    \return A pointer to a TGraph containing the power spectrum. 
  */    
   TGraph *makeRawPowerSpectrum(TGraph *grWave);

   //! Returns the correlation of two TGraphs
  /*!
    \param gr1 The first input TGraph
    \param gr2 The second input TGraph
    \return A pointer to a TGraph containing the correlation of <i>gr1</i> and <i>gr2</i>.
  */    
   TGraph *getCorrelationGraph(TGraph *gr1, TGraph *gr2, Int_t *zeroOffset=0);



   //! Returns the normalised correlation of two TGraphs
  /*!
    \param gr1 The first input TGrap  (must be zero meaned)
    \param gr2 The second input TGraph (must be zero meaned)
    \return A pointer to a TGraph containing the correlation of <i>gr1</i> and <i>gr2</i> where each point is normalised by the number of valid samples in the correlation and by the product of the RMS of the input graphs.
  */    
   TGraph *getNormalisedCorrelationGraph(TGraph *gr1, TGraph *gr2, Int_t *zeroOffset=0);

   //! Returns the normalised correlation of two TGraphs
  /*!
    \param gr1 The first input TGrap  (must be zero meaned)
    \param gr2 The second input TGraph (must be zero meaned)
    \param zeroOffset A pointer to an integer where the sample corresponding to zero offset will be stored
    \param useDtRange A flag to enable the setting of a limited range of deltat's to try
    \param dtMin The minimum delta-t to include in the correlation, the maximum delta-t to include in the correlation
    \return A pointer to a TGraph containing the correlation of <i>gr1</i> and <i>gr2</i> where each point is normalised by the number of valid samples in the correlation and by the product of the RMS of the input graphs.
  */  
   TGraph *getNormalisedCorrelationGraphTimeDomain(TGraph *gr1, TGraph *gr2, Int_t *zeroOffset=0, Int_t useDtRange=0, Double_t dtMin=-1000, Double_t dtMax=1000);

   //! Returns the correlation of two interpolated TGraphs
  /*!
    \param grIn1 The first input TGraph 
    \param grIn2 The second input TGraph 
    \param deltaT The desired time step for the interpolated graphs
    \return A pointer to a TGraph containing the correlation of <i>gr1</i> and <i>gr2</i>, after each is interpolated to have a timestep of <i>deltaT</i>.
  */      
   TGraph *getInterpolatedCorrelationGraph(TGraph *grIn1, TGraph *grIn2, Double_t deltaT);
  //! Returns the inverse FFT of the FFT of the input TGraph. Seems pointless.
  /*!
    \param grWave A pointer to the input TGraph
    \return A pointer to a TGraph containing the inverse FFT of the FFT of <i>grWave</i>.
  */        
   TGraph *makeInverseInverseSpectrum(TGraph *grWave);
  
   //! Returns the time domain result of a frequency domain sum of a number of TGraphs. In the sum each TGraph is weighted by a value. As of writing this documentation, I'm not sure why this would be interesting.
  /*!
    \param numGraphs The number of TGraphs to sum.
    \param grPtr A pointer to an array of <i>numGraphs</i> TGraph pointers.
    \param theWeights An optional array of weights with which to sum the TGraphs.
    \return A pointer to a TGraph containing the inverse FFT of the weighted summed FFT of the input graphs.
  */
   TGraph *combineGraphsUsingFFTs(Int_t numGraphs, TGraph **grPtr,double *theWeights=0);
  
   //! Smooth graph using box car smoothing
  /*!
    \param grWave A pointer to the input TGraph
    \param halfWidth The halfWidth in samples of the size of the box over which to smooth the input graph (eg. <i>halfWidth</i>=1 means that sample i is the average if i-1, i and i+1.
    \return A pointer to a TGraph containing the smoothed graph.
  */
   TGraph *getBoxCar(TGraph *grWave, Int_t halfWidth);
  //! The Hilbert transform of the input TGraph
  /*!
    \param grWave A pointer to the input TGraph
    \return A pointer to a TGraph containing the Hilbert transform
  */
   TGraph *getHilbertTransform(TGraph *grWave);
  //! The Hilbert envelope of the input TGraph. This is defined as e_i=sqrt(v_i^2 + h_i^2), where e_i, v_i and h_i are the i-th sample of the envelope, input graph and hilbert transform of the input graph repsectively.
  /*!
    \param grWave A pointer to the input TGraph
    \return A pointer to a TGraph containing the Hilbert envelope function.
  */
   TGraph *getHilbertEnvelope(TGraph *grWave);

  //Utility functions (not necessarily FFT related but they'll live here for now

  //! The linear sum of the power in a TGraph (normally a PSD)
  /*!
    \param gr A pointer to the input TGraph (normally a PSD)
    \param firstBin The first bin to include in the sum.
    \param lastBin The last bin to include in the sum.
    \return The summed power
  */
   Double_t sumPower(TGraph *gr,Int_t firstBin=-1,Int_t lastBin=-1);
   //! The integral of the power in a TGraph (normally a PSD) (i.e the sum of bin content*bin width)
  /*!
    \param gr A pointer to the input TGraph (normally a PSD)
    \param firstBin The first bin to include in the integral.
    \param lastBin The last bin to include in the integral.
    \return The integral value.
  */
   Double_t integratePower(TGraph *gr,Int_t firstBin=-1,Int_t lastBin=-1);
   //! The sum of the voltage squared in a waveform.
  /*!
    \param gr A pointer to the input TGraph 
    \param firstBin The first bin to include in the sum.
    \param lastBin The last bin to include in the sum.
    \return The value of the sum.
  */  
   Double_t sumVoltageSquared(TGraph *gr,Int_t firstBin,Int_t lastBin); 
  //! The integral of the v^2*dt in a waveform. Now works for unevenly sampled waeforms
  /*!
    \param gr A pointer to the input TGraph 
    \param firstBin The first bin to include in the integral.
    \param lastBin The last bin to include in the integral.
    \return The value of the integral.
  */  
   Double_t integrateVoltageSquared(TGraph *gr,Int_t firstBin=-1,Int_t lastBin=-1);
  //! Find the peak (maximum positive) bin in a TGraph
  /*!
    \param gr A pointer to the input TGraph 
    \return The index of the bin with peak value.
  */    
   Int_t getPeakBin(TGraph *gr); 
 

  //! Find the peak (maximum positive) bin in a TGraph
  /*!
    \param gr A pointer to the input TGraph 
    \return The index of the bin with peak value.
  */    
   Int_t getPeakBin(TGraph *gr, Int_t firstBin, Int_t lastBin);  
  
  //! Find the peak (maximum positive) in a TGraph.
  /*!
    \param gr A pointer to the input TGraph.
    \param index An optional pointer in which the peak bin can be stored.
    \return The peak value.
  */    
   Double_t getPeakVal(TGraph *gr, int *index=0);
  //! Find the peak (v^2) in a TGraph.
  /*!
    \param gr A pointer to the input TGraph.
    \param index An optional pointer in which the peak bin can be stored.
    \return The peak (v^2)value.
  */    
   Double_t getPeakSqVal(TGraph *gr, int *index=0);
  //! Find the peak (v^2) and RMS (of v^2) in a TGraph.
  /*!
    \param gr A pointer to the input TGraph.
    \param peak A reference to a Double_t where the peak value will be stored.
    \param rms A reference to a Double_t where the rms value will be stored.
    \param index An optional pointer in which the peak bin can be stored.
  */      
   void getPeakRmsSqVal(TGraph *gr, Double_t &peak, Double_t &rms, Int_t *index=0);
  
  //! Find the peak (v) and RMS (of v) of a rectified TGraph.
  /*!
    \param gr A pointer to the input TGraph.
    \param peak A reference to a Double_t where the peak value will be stored.
    \param rms A reference to a Double_t where the rms value will be stored.
    \param index An optional pointer in which the peak bin can be stored.
  */      
   void getPeakRmsRectified(TGraph *gr, Double_t &peak, Double_t &rms, Int_t *index=0);

  //Graph returning utility funcs
  //! Returns the simple power envelope of a waveform. The waveform is first squared and then only local peaks are taken to form the power envelope.
  /*!
    \param gr A pointer to the input TGraph.
    \return A pointer to a TGraph containing the power envelope
  */     
   TGraph *getSimplePowerEnvelopeGraph(TGraph *gr);
  //! Returns a smoothed FFT where N-bins are averaged to reduce variance.
  /*!
    \param gr A pointer to the input TGraph.
    \param factor The factor by which to smooth (eg. 2 returns a PSD with half as many frequency bins).
    \return A pointer to a TGraph containing the smoothed PSD.
  */     
   TGraph *smoothFFT(TGraph *gr,Int_t factor) ;
  //! Returns the difference between two graphs (A-B).
  /*!
    \param grA A pointer to the first graph.
    \param grB A pointer to the second graph.
    \return A pointer to a TGraph containing the difference (A-B) between the input graphs.
  */     
   TGraph *subtractGraphs(TGraph *grA, TGraph *grB);
  //! Returns the ratio between two graphs (A/B).
  /*!
    \param grA A pointer to the first graph.
    \param grB A pointer to the second graph.
    \return A pointer to a TGraph containing the ratio (A/B) of the input graphs.
  */     
   TGraph *divideGraphs(TGraph *grA, TGraph *grB);

  
  //! Returns a graph translated by deltaT. Such that t'=t+dt
  /*!
    \param grWave A pointer to the input graph.
    \param deltaT The amount to move time axis by
    \return A pointer to a TGraph containing the translated graph
  */     
   TGraph *translateGraph(TGraph *grWave, Double_t deltaT);

  //! Returns the one minus the ratio between two graphs (A/B).
  /*!
    \param grA A pointer to the first graph.
    \param grB A pointer to the second graph.
    \return A pointer to a TGraph containing the one minus the ratio (1 - A/B) of the input graphs.
  */     
   TGraph *ratioSubtractOneGraphs(TGraph *grA, TGraph *grB) ;
  //! Returns the ratio of two graphs as 10 *log(A/B)
  /*!
    \param grA A pointer to the first graph.
    \param grB A pointer to the second graph.
    \return A pointer to a TGraph containing 10 * log(A/B)
  */     
   TGraph *dbGraphs(TGraph *grA, TGraph *grB);
  //! Zero pad a wave making it a factor of N longer.
  /*!
    \param grA A pointer to the input graph.
    \param padFactor The factor by whcih to increas the length (eg. new length = <i>factor</i> * old length)
    \return A pointer to the zero padded TGraph.
  */     
   TGraph *padWave(TGraph *grA, Int_t padFactor);
  //! Zero pad a wave making it up to newLength points.
  /*!
    \param grA A pointer to the input graph.
    \param newLength The new length the graph should be (it gets zero padded on either side)
    \return A pointer to the zero padded TGraph.
  */     
   TGraph *padWaveToLength(TGraph *grA, Int_t newLength);
  //! Rectify a waveform, optionally returning an all negative waveform.
  /*!
    \param gr A pointer to the input graph.
    \param makeNeg An optional parameter which if true will return a negative waveform.
    \return A pointer to the rectified TGraph
  */     
   TGraph *rectifyWave(TGraph *gr, Int_t makeNeg=0);
  
  //Window functions
  //! The Bartlett window function (it's basically a triangle that peaks in the middle at 1 and goes to zero at the end points)
  /*!
    \param j The index to evaluate the window at.
    \param n The length of sample.
    \return The value of the window function at index <i>j</i>.
  */       
   Double_t  bartlettWindow(Int_t j, Int_t n);
  //! The Welch window function similar to the Bartlett window except that it falls off from the middle as dx^2 (ie. less steeply).
  /*!
    \param j The index to evaluate the window at
    \param n The length of sample.
    \return The value of the window function at index <i>j</i>.
  */       
   Double_t welchWindow(Int_t j, Int_t n);



  //! This fucntions just calculates the simple bin by bin dv/dt derivative of the input data.
  /*!
    \param numPoints The number of points in the input data (the output will be numPoints-1 in length).
    \param inputX The input time values.
    \param inputY The input voltage values.
    \param outputX The output time values.
    \param outputY The output voltage values.
  */       
   void takeDerivative(Int_t numPoints, Double_t *inputX, Double_t *inputY, Double_t *outputX, Double_t *outputY);

 //! This returns a TGraph which is the derivative of the input graph
  /*!
    \param grIn The input graph.
    \return The derviative of grIn.
  */      
   TGraph *getDerviative(TGraph *grIn);
 //! This returns a TGraph which is the derivative of the input graph
  /*!
    \param grIn The input graph.
    \return The derviative of grIn.
  */      
   TGraph *getDerivative(TGraph *grIn);

 //! This returns a TGraph which has had a simple pass band filter applied
  /*!
    \param grWave The input graph.
    \param minFreq The lowest frequency to pass.
    \param maxFreq The highest frequency to pass.
    \return The derviative of grWave.
  */      
   TGraph *simplePassBandFilter(TGraph *grWave, Double_t minFreq, Double_t maxFreq);

 //! This returns a TGraph which has had a simple notch band filter applied
  /*!
    \param grWave The input graph.
    \param minFreq The lower frequency of the notch.
    \param maxFreq The upper frequency of the notch.
    \return The derviative of grWave.
  */      
   TGraph *simpleNotchFilter(TGraph *grWave, Double_t minFreq, Double_t maxFreq);


  //! This returns a TGraph which has had N simple notch band filters applied
  /*!
    \param grWave The input graph.
    \param numNotches The number of notch regiosn
    \param minFreq An array of lower frequency of the notches.
    \param maxFreq An array upper frequency of the notch.
    \return The derviative of grWave.
  */      
   TGraph *multipleSimpleNotchFilters(TGraph *grWave, Int_t numNotches, Double_t minFreq[], Double_t maxFreq[]);

//! This returns a TGraph which has been cropped in.
  /*!
    \param grWave The input graph.
    \param minTime The lower edge of the time window.
    \param maxTime The upper edge of the time window.
    \return The derviative of grWave.
  */      
   TGraph *cropWave(TGraph *grWave, Double_t minTime, Double_t maxTime);

  
  //! This returns the SNR ratio of the input waveform
  /*!
    \param gr The input graph.
    \return The SNR of the input waveform, where S is half of the peak-to-peak and N is the RMS of the first 25 ns.
  */      
   Double_t getWaveformSNR(TGraph *gr); 
  //! This returns the SNR ratio of the input waveform
  /*!
    \param gr The input graph.
    \param peakToPeak A reference to a double which will be set to half of the peak-to-peak
    \param rms A reference to a double which will be set to the RMS of teh first 25 samples
    \return The SNR of the input waveform, where S is half of the peak-to-peak and N is the RMS of the first 25 samples.
  */     
   Double_t getWaveformSNR(TGraph *gr,Double_t &peakToPeak,Double_t &rms);

  //! This returns the largest (i.e most positive, or least negative) value
  /*!
    \param gr The input graph.
    \return The peak
  */      
   Double_t getWaveformPeak(TGraph *gr);
   Double_t getEnvelopeSNR(TGraph *gr);
   Double_t getEnvelopeSNR(TGraph *gr,Double_t &peakToPeak,Double_t &rms,Double_t &timeOfPeak);


  //! Linear interpolate to find the value at some point near two known points
  /*!
    \param x1 x value of point 1
    \param y1 y value of point 1
    \param x2 x value of point 2
    \param y2 y value of point 2
    \param x is the point of interest
    \return y
  */      
   Double_t simpleInterploate(Double_t x1, Double_t y1, Double_t x2, Double_t y2,Double_t x);
  
   /*!
    * Correlation between two FFTs
    * bandpasses between min_i and max_i using butterworth filter of order order
    * work can be used for temporary to avoid allocation of new memory 
    */
   double * FFTCorrelation(int waveformlength, const FFTWComplex * A, const FFTWComplex * B, FFTWComplex * work = 0, 
                           int min_i = 0, int max_i =0, int order=1);  


   /*! in place array rotation
    */
   void inPlaceShift(int N, double *x); 


   /*! inplace graph rotation */ 
   void rotate(TGraph * g, int rot); 

   /*! does fit to polynomial of order order then subtracts it */
   void polySubtract(TGraph *g, int order=1); 


      
   /*! direct convolution for kernel can either assume zero's outside boundaries or repetition of first and last values */
   enum DirectConvolveEdgeBehavior
   {
        ZEROES_OUTSIDE, 
        REPEAT_OUTSIDE
   }; 


    /*!convolution without FFT of x with kernel h. y should have same size as x (it is cropped symmetrically) 
      if delay = 0, h will be centered around its middle sample, if - M/2 , will be purely causal, etc. 
    */
    double * directConvolve(int N, const double *x, int M, const double * h, double *y = 0, int delay = 0,  DirectConvolveEdgeBehavior edge_behavior = ZEROES_OUTSIDE); 


    /*! wraps periodic array of doubles. center assumed to be period/2 */ 
   void wrap(size_t N, double * vals, double period = 360); 

    /*! wraps periodic array of doubles. */ 
   void wrap(size_t N, double * vals, double period, double center); 

    /*! wraps periodic array of floats */ 
   void wrap(size_t N, float * vals, float period, float center); 

    /*! wraps periodic array of floats. center assumed to be period/2 */ 
   void wrap(size_t N, float * vals, float period = 360); 


   /*! Faster than the libm floor, which cares about errno and the floating point environment and all that crap 
    **/ 
   inline int fast_floor(double val) { return (int) val - (val < (int) val); }

   /*!Wraps a value to be within period centered at center.
    * This could be implemented in terms of the other wraps, but in practice it's better to inline since
    * the period and center are usually known at compile them. 
    *
    * */ 
   inline double wrap(double val, double period, double center)
          { return val - period * fast_floor((val-center+period/2)/period); }

   /*! Wraps a value to be within period. Center assumed to be period/2 */ 
   inline double wrap(double val, double period = 360)
          { return wrap(val,period,period/2); } 


    /*! unwraps periodic array of doubles. */ 
   void unwrap(size_t N, double * vals, double period = 360); 


    /*! unwraps periodic array of floats. */ 
   void unwrap(size_t N, float * vals, float period = 360); 


   /*! computes reasonable dt from unevenly sampled graph */ 
   double getDt(const TGraph * g, int realN = 0);  

   /*! linearly interpolates value x in g . Similar to doing g->Eval() but much faster (since doesn't need to sort or binary search)*/ 
   double evalEvenGraph(const TGraph * g, double x); 



   /*! applies window to graph */ 
   void applyWindow(TGraph *g, const FFTWindowType *w); 

   
   /*! random variable from rayleigh distribution */ 
   double randomRayleigh(double sigma=1, TRandom * rng = 0); 

   /*! sinc function, if |x| < eps, sinc(x) = 1 */ 
   double sinc(double x, double eps = 0); 
 
   int loadWisdom(const char * file); 
   int saveWisdom(const char * file); 


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
   double *  lombScarglePeriodogramSlow(int N, const double *x, const double * y, int nfreqs, const double * freqs, double * answer = 0); 

   /** fast periodogoram (as in Press & Rybicki) . Implementation in Periodogram.cxx */
   TGraph * lombScarglePeriodogram(const TGraph * g, double dt = 0, double oversample_factor  = 4 , 
                       double high_factor = 2, TGraph * replaceme = 0, int extirpolation_factor =4)  ; 
   TGraph * lombScarglePeriodogram(int N, double dt, const double * __restrict x, 
                                   const double * __restrict y, double oversample_factor  = 4 , 
                                   double high_factor = 2, TGraph * replaceme = 0, int extirpolation_factor = 4)  ; 


//   TH2 * getPowerVsTimeUsingLombScargle(const TGraph * g, int nbins, double sampling_dt = 0, double oversample_factor = 4, double high_factor = 2, TH2 * useme = 0); 

   /* Compute Stokes parameters from hpol / vpol and their hilbert transforms 
    * @param N number of samples
    * @param hpol hpol waveform (evenly sampled)
    * @param hpol_hat hilbert transform of hpol waveform (evenly sampled)
    * @param vpol vpol waveform (evenly sampled)
    * @param vpol_hat hilbert transform of vpol waveform (evenly sampled)
    * @param I pointer to where Stokes I will be stored (or null to not store)
    * @param Q pointer to where Stokes Q will be stored (or null to not store)
    * @param U pointer to where Stokes U will be stored (or null to not store)
    * @param V pointer to where Stokes V will be stored (or null to not store)
    */ 

   void stokesParameters(int N, const double * __restrict hpol, const double * __restrict hpol_hat, const double * __restrict vpol, const double * __restrict vpol_hat, 
                         double * I = 0, double * Q = 0, double * U = 0, double * V = 0); 



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

   void dftAtFreq(const TGraph * g, double freq, double * phase, double * amp = 0, double * real = 0, double * imag = 0); 


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
   void dftAtFreqAndMultiples(const TGraph * g, double freq, int nmultiples, double * phase, double * amp = 0, double * real = 0, double * imag = 0); 

}
   
#endif //FFTTOOLS_H
