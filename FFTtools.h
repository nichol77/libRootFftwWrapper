
#ifndef FFTTOOLS_H
#define FFTTOOLS_H

//#include <fftw3.h>
#include "TGraph.h"
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"

//My includes
#include "FFTWComplex.h"

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
 * -# Checkout the code from the SVN repository, eg.: <BR><PRE>svn co https://delos.mps.ohio-state.edu/anitaGround/libRootFftwWrapper/trunk myRootFftwWrapperDir</PRE>
 * -# Define the ANITA_UTIL_INSTALL_DIR to point to the location you want the library installed (the library files will end up in (ANITA_UTIL_INSTALL_DIR)/lib and the header files in (ANITA_UTIL_INSTALL_DIR)/include).
 * -# Do <PRE>make</PRE><PRE>make install</PRE>
 * \section manual_sec Manual
 * If you are averse to reading web pages (and who wouldn't be) you can download a <a href="manual/libRootFftwWrapper.pdf">pdf copy of the reference material</a> but be warned it won't be a thrilling read.
 */



//!  This is the class that holds all the various useful FFT related functions
/*!
  This is the class that holds all the various useful FFT related functions. All of the functions can (and should) be called statically using something like FFTtools::doYourFunkyThing(arg1,arg2).
*/
class FFTtools
{
public:
  FFTtools(); ///<Constructor, not used as all member functions are static
  ~FFTtools();  ///<Destructor, also not used. 

    
  
  //! Interpolation Routines that use ROOT::Math::Interpolator 
  /*!
    \param grIn A pointer to the input TGraph.
    \param deltaT The desired period (1/rate) of the interpolated waveform.
    \return A pointer to ther interpolated TGraph, it is the users responsibility to delete this after use.
  */
  static TGraph *getInterpolatedGraph(TGraph *grIn, Double_t deltaT);

  //! Returns the magnitude of a complex number.
  /*!
    \param theNum The complex number.
    \return The magnitude of the complex number.
  */
  static double getAbs(FFTWComplex &theNum);

  //! Computes an inverse FFT
  /*!
    \param length The length of the output array
    \param theInput The input array of complex numbers of <i>(length/2 +1)</i>
    \return An array of <i>length</i> real numbers
  */
  static double *doInvFFT(int length, FFTWComplex *theInput);


  //! Computes an FFT of an array of real numbers.
  /*!
    \param length The length of the input array.
    \param theInput The input array of <i>length</i> real numbers.
    \return An array of <i>length/2 + 1</i> complex numbers
  */
  static FFTWComplex *doFFT(int length,double *theInput);    

  //! Computes the correlation of two subsets of TGraphs
  /*!
    \param gr1 The first TGraph in the correlation.
    \param gr2 The second TGraph in the correlation.
    \param firstIndex The index of the first sample in the sub traces.
    \param lastIndex The index of the last sample in the sub traces.
    \return The correlation as an array of <i>lastIndex-firstIndex</i> real numbers.
  */
  static double *getCorrelation(TGraph *gr1, TGraph *gr2,int firstIndex,int lastIndex);
  //! Computes the correlation of two arrays
  /*!
    \param length The length of the arrays
    \param oldY1 The first array in the correlation.
    \param oldY2 The second array in the correlation.
    \return The correlation as an array of <i>length</i> real numbers.
  */
  static double *getCorrelation(int length,float *oldY1, float *oldY2);
  //! Computes the correlation of two arrays
  /*!
    \param length The length of the arrays
    \param oldY1 The first array in the correlation.
    \param oldY2 The second array in the correlation.
    \return The correlation as an array of <i>length</i> real numbers.
  */
  static double *getCorrelation(int length,double *oldY1, double *oldY2);

  //! Returns the time domain result of a frequency domain sum of a number of arrays. As of writing this documentation, I'm not sure why this would be interesting.
  /*!
    \param numArrays The number of arrays to sum.
    \param thePtrPtr A pointer to a two-dimensional array of doubles <i>[numArrays][eachLength]</i>.
    \param eachLength The length of each array.
    \return The time domain result of the frequency domain summation.
  */
  static Double_t *combineValuesUsingFFTs(Int_t numArrays, Double_t **thePtrPtr, Int_t eachLength);

  //Higher level functions that take and return TGraphs
  //! Returns the power spectral density. Note the PSD is unormalised (or if you prefer is normalised to the sum squared amplitude of the time domain). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a>
  /*!
    \param grWave The input time domain waveform
    \return A pointer to a TGraph containing the power spectrum. 
  */
  static TGraph *makePowerSpectrum(TGraph *grWave);
  //! Returns the power spectral density. Note the PSD returned is the periodogram (or if you prefer is normalised to the mean squared amplitude of the time domain). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a>
  /*!
    \param grWave The input time domain waveform
    \return A pointer to a TGraph containing the power spectrum. 
  */
  static TGraph *makePowerSpectrumPeriodogram(TGraph *grWave);
  //! Returns the power spectral density. Note the PSD returned is normalised and divided by frequency bin width (or if you prefer it is normalised to the time-integral squared amplitude of the time domain and then divided by frequency bin width). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a> As the name suggests this function expects the input waveform to be a volts-seconds one.
  /*!
    \param grWave The input time domain waveform with units of volts-seconds.
    \return A pointer to a TGraph containing the power spectrum. 
  */  
  static TGraph *makePowerSpectrumVoltsSeconds(TGraph *grWave);
  //! Returns the power spectral density of the input waveform convolved with a Bartlett Window. Note the PSD returned is normalised and divided by frequency bin width (or if you prefer it is normalised to the time-integral squared amplitude of the time domain and then divided by frequency bin width). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a> As the name suggests this function expects the input waveform to be a volts-seconds one.
  /*!
    \param grWave The input time domain waveform with units of volts-seconds.
    \return A pointer to a TGraph containing the power spectrum. 
  */  
  static TGraph *makePowerSpectrumVoltsSecondsBartlett(TGraph *grWave);
  //! Returns the power spectral density of the input waveform. In this one we first zero pad the waveform and then split it up into overlapping segments and convolve each segment with the Bartlett window before summing the resulting PSD's. As the name suggests this function expects the input waveform to be a volts-seconds one. No idea if this one actually works, or where I read oabout this crazy method which is supposed to reduce the variance of the PSD estimator.
  /*!
    \param grWave The input time domain waveform with units of volts-seconds.
    \param padFactor The factor by which to zero pad the wave see padWave
    \param numFreqs The number of frequency bins required in the output, which is related to the length of overlapping segments.
    \return A pointer to a TGraph containing the power spectrum. 
  */    
  static TGraph *makePSVSBartlettPaddedOverlap(TGraph *grWave, Int_t padFactor=4, Int_t numFreqs=64);

  //! Returns the power spectral density in dB units. Note the PSD returned is normalised and divided by frequency bin width (or if you prefer it is normalised to the time-integral squared amplitude of the time domain and then divided by frequency bin width). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a> As the name suggests this function expects the input waveform to be a volts-seconds one.
  /*!
    \param grWave The input time domain waveform with units of volts-seconds.
    \return A pointer to a TGraph containing the power spectrum in dB units, the frequency units are MHz. 
  */  
  static TGraph *makePowerSpectrumVoltsSecondsdB(TGraph *grWave);
 //! Returns the power spectral density in dB units. Note the PSD returned is normalised and divided by frequency bin width (or if you prefer it is normalised to the time-integral squared amplitude of the time domain and then divided by frequency bin width). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a> As the name suggests this function expects the input waveform to be a millivolts-nanoseconds one.
  /*!
    \param grWave The input time domain waveform with units of millivolts-nanoseconds.
    \return A pointer to a TGraph containing the power spectrum in dB units, the frequency units are MHz. 
  */  
  static TGraph *makePowerSpectrumMilliVoltsNanoSecondsdB(TGraph *grWave);
  //! Returns the power spectral density of the input waveform zero-padded by some factor. Note the PSD returned is normalised and divided by frequency bin width (or if you prefer it is normalised to the time-integral squared amplitude of the time domain and then divided by frequency bin width). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a> As the name suggests this function expects the input waveform to be a volts-seconds one.
  /*!
    \param grWave The input time domain waveform with units of volts-seconds.
    \param padFactor The factor by which to zero pad the wave.
    \return A pointer to a TGraph containing the power spectrum. 
  */  
  static TGraph *makePowerSpectrumVoltsSecondsPadded(TGraph *grWave, Int_t padFactor=4);
  //! Returns the power spectral density in dB units of the input waveform zero-padded by some factor. Note the PSD returned is normalised and divided by frequency bin width (or if you prefer it is normalised to the time-integral squared amplitude of the time domain and then divided by frequency bin width). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a> As the name suggests this function expects the input waveform to be a volts-seconds one.
  /*!
    \param grWave The input time domain waveform with units of volts-seconds.
    \param padFactor The factor by which to zero pad the wave.
    \return A pointer to a TGraph containing the power spectrum in dB units with MHz as the frequency unit. 
  */    
  static TGraph *makePowerSpectrumVoltsSecondsPaddeddB(TGraph *grWave, Int_t padFactor=4);
  
  //! Returns the power spectral density in completely unormalised unit (as in Parseval's theorem is not obeyed and there is an extra factor of N not removed form the PSD). <a href="http://www.hep.ucl.ac.uk/~rjn/saltStuff/fftNormalisation.pdf">See this short note for my terminology.</a>
  /*!
    \param grWave The input time domain waveform.
    \return A pointer to a TGraph containing the power spectrum. 
  */    
  static TGraph *makeRawPowerSpectrum(TGraph *grWave);
   //! Returns the correlation of two TGraphs
  /*!
    \param gr1 The first input TGraph
    \param gr2 The second input TGraph
    \return A pointer to a TGraph containing the correlation of <i>gr1</i> and <i>gr2</i>.
  */    
  static TGraph *getCorrelationGraph(TGraph *gr1, TGraph *gr2);
   //! Returns the correlation of two interpolated TGraphs
  /*!
    \param grIn1 The first input TGraph
    \param grIn2 The second input TGraph
    \param deltaT The desired time step for the interpolated graphs
    \return A pointer to a TGraph containing the correlation of <i>gr1</i> and <i>gr2</i>, after each is interpolated to have a timestep of <i>deltaT</i>.
  */      
  static TGraph *getInterpolatedCorrelationGraph(TGraph *grIn1, TGraph *grIn2, Double_t deltaT);
  //! Returns the inverse FFT of the FFT of the input TGraph. Seems pointless.
  /*!
    \param grWave A pointer to the input TGraph
    \return A pointer to a TGraph containing the inverse FFT of the FFT of <i>grWave</i>.
  */        
  static TGraph *makeInverseInverseSpectrum(TGraph *grWave);
  
   //! Returns the time domain result of a frequency domain sum of a number of TGraphs. In the sum each TGraph is weighted by a value. As of writing this documentation, I'm not sure why this would be interesting.
  /*!
    \param numGraphs The number of TGraphs to sum.
    \param grPtr A pointer to an array of <i>numGraphs</i> TGraph pointers.
    \param theWeights An optional array of weights with which to sum the TGraphs.
    \return A pointer to a TGraph containing the inverse FFT of the weighted summed FFT of the input graphs.
  */
  static TGraph *combineGraphsUsingFFTs(Int_t numGraphs, TGraph **grPtr,double *theWeights=0);
  
   //! Smooth graph using box car smoothing
  /*!
    \param grWave A pointer to the input TGraph
    \param halfWidth The halfWidth in samples of the size of the box over which to smooth the input graph (eg. <i>halfWidth</i>=1 means that sample i is the average if i-1, i and i+1.
    \return A pointer to a TGraph containing the smoothed graph.
  */
  static TGraph *getBoxCar(TGraph *grWave, Int_t halfWidth);
  //! The Hilbert transform of the input TGraph
  /*!
    \param grWave A pointer to the input TGraph
    \return A pointer to a TGraph containing the Hilbert transform
  */
  static TGraph *getHilbertTransform(TGraph *grWave);
  //! The Hilbert envelope of the input TGraph. This is defined as e_i=sqrt(v_i^2 + h_i^2), where e_i, v_i and h_i are the i-th sample of the envelope, input graph and hilbert transform of the input graph repsectively.
  /*!
    \param grWave A pointer to the input TGraph
    \return A pointer to a TGraph containing the Hilbert envelope function.
  */
  static TGraph *getHilbertEnvelope(TGraph *grWave);

  //Utility functions (not necessarily FFT related but they'll live here for now

  //! The linear sum of the power in a TGraph (normally a PSD)
  /*!
    \param gr A pointer to the input TGraph (normally a PSD)
    \param firstBin The first bin to include in the sum.
    \param lastBin The last bin to include in the sum.
    \return The summed power
  */
  static Double_t sumPower(TGraph *gr,Int_t firstBin=-1,Int_t lastBin=-1);
   //! The integral of the power in a TGraph (normally a PSD) (i.e the sum of bin content*bin width)
  /*!
    \param gr A pointer to the input TGraph (normally a PSD)
    \param firstBin The first bin to include in the integral.
    \param lastBin The last bin to include in the integral.
    \return The integral value.
  */
  static Double_t integratePower(TGraph *gr,Int_t firstBin=-1,Int_t lastBin=-1);
   //! The sum of the voltage squared in a waveform.
  /*!
    \param gr A pointer to the input TGraph 
    \param firstBin The first bin to include in the sum.
    \param lastBin The last bin to include in the sum.
    \return The value of the sum.
  */  
  static Double_t sumVoltageSquared(TGraph *gr,Int_t firstBin,Int_t lastBin); 
  //! The integral of the v^2*dt in a waveform.
  /*!
    \param gr A pointer to the input TGraph 
    \param firstBin The first bin to include in the integral.
    \param lastBin The last bin to include in the integral.
    \return The value of the integral.
  */  
  static Double_t integrateVoltageSquared(TGraph *gr,Int_t firstBin=-1,Int_t lastBin=-1);
  //! Find the peak (maximum positive) bin in a TGraph
  /*!
    \param gr A pointer to the input TGraph 
    \return The index of the bin with peak value.
  */    
  static Int_t getPeakBin(TGraph *gr); 
 

  //! Find the peak (maximum positive) bin in a TGraph
  /*!
    \param gr A pointer to the input TGraph 
    \return The index of the bin with peak value.
  */    
  static Int_t getPeakBin(TGraph *gr, Int_t firstBin, Int_t lastBin);  
  
  //! Find the peak (maximum positive) in a TGraph.
  /*!
    \param gr A pointer to the input TGraph.
    \param index An optional pointer in which the peak bin can be stored.
    \return The peak value.
  */    
  static Double_t getPeakVal(TGraph *gr, int *index=0);
  //! Find the peak (v^2) in a TGraph.
  /*!
    \param gr A pointer to the input TGraph.
    \param index An optional pointer in which the peak bin can be stored.
    \return The peak (v^2)value.
  */    
  static Double_t getPeakSqVal(TGraph *gr, int *index=0);
  //! Find the peak (v^2) and RMS (of v^2) in a TGraph.
  /*!
    \param gr A pointer to the input TGraph.
    \param peak A reference to a Double_t where the peak value will be stored.
    \param rms A reference to a Double_t where the rms value will be stored.
    \param index An optional pointer in which the peak bin can be stored.
  */      
  static void getPeakRmsSqVal(TGraph *gr, Double_t &peak, Double_t &rms, Int_t *index=0);
  
  //! Find the peak (v) and RMS (of v) of a rectified TGraph.
  /*!
    \param gr A pointer to the input TGraph.
    \param peak A reference to a Double_t where the peak value will be stored.
    \param rms A reference to a Double_t where the rms value will be stored.
    \param index An optional pointer in which the peak bin can be stored.
  */      
  static void getPeakRmsRectified(TGraph *gr, Double_t &peak, Double_t &rms, Int_t *index=0);

  //Graph returning utility funcs
  //! Returns the simple power envelope of a waveform. The waveform is first squared and then only local peaks are taken to form the power envelope.
  /*!
    \param gr A pointer to the input TGraph.
    \return A pointer to a TGraph containing the power envelope
  */     
  static TGraph *getSimplePowerEnvelopeGraph(TGraph *gr);
  //! Returns a smoothed FFT where N-bins are averaged to reduce variance.
  /*!
    \param gr A pointer to the input TGraph.
    \param factor The factor by which to smooth (eg. 2 returns a PSD with half as many frequency bins).
    \return A pointer to a TGraph containing the smoothed PSD.
  */     
  static TGraph *smoothFFT(TGraph *gr,Int_t factor) ;
  //! Returns the difference between two graphs (A-B).
  /*!
    \param grA A pointer to the first graph.
    \param grB A pointer to the second graph.
    \return A pointer to a TGraph containing the difference (A-B) between the input graphs.
  */     
  static TGraph *subtractGraphs(TGraph *grA, TGraph *grB);
  //! Returns the ratio between two graphs (A/B).
  /*!
    \param grA A pointer to the first graph.
    \param grB A pointer to the second graph.
    \return A pointer to a TGraph containing the ratio (A/B) of the input graphs.
  */     
  static TGraph *divideGraphs(TGraph *grA, TGraph *grB);
  //! Returns the one minus the ratio between two graphs (A/B).
  /*!
    \param grA A pointer to the first graph.
    \param grB A pointer to the second graph.
    \return A pointer to a TGraph containing the one minus the ratio (1 - A/B) of the input graphs.
  */     
  static TGraph *ratioSubtractOneGraphs(TGraph *grA, TGraph *grB) ;
  //! Returns the ratio of two graphs as 10 *log(A/B)
  /*!
    \param grA A pointer to the first graph.
    \param grB A pointer to the second graph.
    \return A pointer to a TGraph containing 10 * log(A/B)
  */     
  static TGraph *dbGraphs(TGraph *grA, TGraph *grB);
  //! Zero pad a wave making it a factor of N longer.
  /*!
    \param grA A pointer to the input graph.
    \param padFactor The factor by whcih to increas the length (eg. new length = <i>factor</i> * old length)
    \return A pointer to the zero padded TGraph.
  */     
  static TGraph *padWave(TGraph *grA, Int_t padFactor);
  //! Rectify a waveform, optionally returning an all negative waveform.
  /*!
    \param gr A pointer to the input graph.
    \param makeNeg An optional parameter which if true will return a negative waveform.
    \return A pointer to the rectified TGraph
  */     
  static TGraph *rectifyWave(TGraph *gr, Int_t makeNeg=0);
  
  //Window functions
  //! The Bartlett window function (it's basically a triangle that peaks in the middle at 1 and goes to zero at the end points)
  /*!
    \param j The index to evaluate the window at.
    \param n The length of sample.
    \return The value of the window function at index <i>j</i>.
  */       
  static Double_t  bartlettWindow(Int_t j, Int_t n);
  //! The Welch window function similar to the Bartlett window except that it falls off from the middle as dx^2 (ie. less steeply).
  /*!
    \param j The index to evaluate the window at
    \param n The length of sample.
    \return The value of the window function at index <i>j</i>.
  */       
  static Double_t welchWindow(Int_t j, Int_t n);



  //! This fucntions just calculates the simple bin by bin dv/dt derivative of the input data.
  /*!
    \param numPoints The number of points in the input data (the output will be numPoints-1 in length).
    \param inputX The input time values.
    \param inputY The input voltage values.
    \param outputX The output time values.
    \param outputY The output voltage values.
  */       
  static void takeDerivative(Int_t numPoints, Double_t *inputX, Double_t *inputY, Double_t *outputX, Double_t *outputY);

 //! This returns a TGraph which is the derivative of the input graph
  /*!
    \param grIn The input graph.
    \return The derviative of grIn.
  */      
  static TGraph *getDerviative(TGraph *grIn);

 //! This returns a TGraph which has had a simple pass band filter applied
  /*!
    \param grWave The input graph.
    \param minFreq The lowest frequency to pass.
    \param maxFreq The highest frequency to pass.
    \return The derviative of grWave.
  */      
  static TGraph *simplePassBandFilter(TGraph *grWave, Double_t minFreq, Double_t maxFreq);

 //! This returns a TGraph which has had a simple notch band filter applied
  /*!
    \param grWave The input graph.
    \param minFreq The lower frequency of the notch.
    \param maxFreq The upper frequency of the notch.
    \return The derviative of grWave.
  */      
  static TGraph *simpleNotchFilter(TGraph *grWave, Double_t minFreq, Double_t maxFreq);

};
   
#endif //FFTTOOLS_H
