#ifndef _SINE_SUBTRACT_H
#define _SINE_SUBTRACT_H

#ifdef DONT_HAVE_MINUIT2
#warning "Need Minuit2 for SineSubtract." 
#else


/**************************************************************************************
 * \class FFTtools::SineSubtract 
 *
 * \brief This class contains an independent implementation of the sine
 * subtraction CW removal algorithm, suggested by Andres/ Steph
 * (ANITA Elog 621).
 *
 * Cosmin Deaconu <cozzyd@kicp.uchicago.edu> is responsible for this mess. 
 *  
 * This, SineSubtract,  is the controlling class. Related classes are
 * SineFitter, which does the minimization, SineFitter::SineFitFn, which
 * provides the objective function and its derivative, SineFitterLimits, which
 * encodes parameter limits, and SineSubtractResult, which stores stuff. 
 *
 * Many options are available, some of which may even be documented. For
 * example, the part of the waveform used to train sine subtraction may be
 * limited (e.g. if you wanted to just use the non-trigger portion).  Multiple
 * waveforms may be fit for at once (a single frequency is fit for in each
 * iteration, with the amplitude and phase allowed to vary per waveform). It is
 * also possible to enable storage of intermediate graphs and power spectra at
 * each step.  This is useful for diagnostics and making pretty plots, but
 * incurs some overhead, so is not recommended for production use.
 * Gratuitously, there are methods to generate plots and even Beamer slides. 
 *
 * A sketch of the algorithm is: 
 *
 *    - Rectify waveform and estimate power spectrum (according to the power spectrum option). 
 *
 *    - Find peak, fmax, in power spectrum by finding max magnitude/(1+nfails)^nfails_exp) , 
 *      according to the peak finding option. See below for meaning of nfails. 
 *
 *    - Find sine ( A (sin wt + phi)) that when subtracted minimizes "power"
 *    (avg(v_i^2)); initialize with w = 2pi fmax, A = mag, phase from fft.
 *    Uninterpolated graph is used for this! 
 *
 *    - If power reduced by more than min_power_reduction, subtract 
 *
 *    - If power not reduced by more than min_power reduction, increment nfails
 *    for that bin 
 *
 *    - Repeat until max_iter_without_reduction iterations with not enough
 *    improvement   
 *
 * 
 * --- 
 *
 *  A lot of time has been spent profiling and trying to make this run as fast
 *  as possible. This results in code that is more difficult to understand than
 *  a naive implementation would allow. Any suggestions for performance
 *  enhancements are awelcome. 
 * 
 *  If ENABLE_VECTORIZE is defined, custom vectorized code (using Agner Fog's
 *  VCL, http://www.agner.org/optimize/#vectorclass) is used for the objective
 *  function and gradient.  At least until compilers get smart enough to
 *  autovectorize stuff better, this is substantially faster, but obviously
 *  only makes sense in a processor with SIMD instructions.  256-bit floating
 *  point registers are used (e.g. AVX) but the VCL will gracefully fallback to
 *  128-bit registers at compile time should your processor be too old. At some
 *  point AVX-512 will appear, so that might require a few changes. 
 *
 *  **WARNING**: You will not get the same "answers" in vectorized mode. This
 *  is expected, due to usage of different math functions and the general
 *  non-associativity of of floating point math. 
 *
 * ---
 *
 *  A number of additional compile-time options influence the behavior: 
 *
 *  - If SINE_SUBTRACT_USE_FLOATS is defined, the fitter will use floats
 *  instead of doubles. In principle, this should be twice the vectorization,
*  but in practice, the overhead of copying and casting seems to be
*  problematic. 
 *
 *  - If SINE_SUBTRACT_FORCE_ALIGNED is defined, the fitter will copy the
 *  doubles into an aligned array. This might make a performance difference on
 *  old processors, but seems to not help for the most part. 
 *  
 *  - If SINE_SUBTRACT_PROFILE is defined, the time spent calling subtractCW is
 *  printed. 
 *
 *  - if SINE_SUBTRACT_DONT_HORIZONTAL_ADD is defined, horizontal addition is
 *  disabled.  You probably don't want to do this, I added this when I was
 *  investigating some strange numerical behavior. 
 *
 *  See macros/testArtificialSubtract.C to see SineSubtract in action.  
 *
 * --- 
 *
 * TODO list: 
 *
 *  - Better  tuning of default step sizes / limits ? 
 *
 *  - Minuit2 adds quite a bit of overhead, and is one of the things preventing
 *  this from running even faster
 *
 ********************************************************************************************/ 

#include "Minuit2/Minuit2Minimizer.h" 
#include "FFTtools.h" 
#include "TObject.h"

#include <vector>

class TGraph; 
class TCanvas; 

namespace FFTtools
{

  /** The SineFitterLimits class holds a few constraints that are used to limit
   *  the parameter space of the fit. 
   *  
   *  The things that may be contrained are 
   *
   *    - How far each fit amplitude trace can deviate from the guess (two multiplicative factors) 
   *    - How far the fit frequency can deviate from the guess, in terms of fnyq / nsamples. 
   *
   * Any value set to zero means no limit on that parameter. 
   *
   *  I can't think of a good way to constrain the phase. 
   *
   *
   */
   class SineFitterLimits
   {
     public: 

       /** The maximum multiplicative factor the amplitudes can deviate from the guess. Negative to disable, 0 to fix amplitude. */ 
       double maxA_relative_to_guess; 

       /** The minimum multiplicative factor the amplitudes can deviate from the guess. Negative to disable, 0 to fix amplitude.  Note that if either min or max is 0, then it will be fixed.*/ 
       double minA_relative_to_guess; 

       /** The maximum multiplative factor of fnyq/nsamples the fit frequency can deviate from the guess. Negative to disable, 0 to fix frequency.  Note that if either min or max is 0, then it will be fixed.*/
       double max_n_df_relative_to_guess; 

       /** The starting error of the phase, in radians */
       double phase_start_error; 

       /** The fractional starting error of the amplitude */ 
       double amp_start_error; 

       /** The error, as a fraction of df, of the frequency estimate */ 
       double freq_start_error; 



       SineFitterLimits(double maxA_relative_to_guess_ = 4, double minA_relative_to_guess_ = 0.25, 
                        double max_n_df_relative_to_guess_ = 1, double phase_start_error_ = TMath::Pi() / 16, 
                        double amp_start_error_ = 0.1, double freq_start_error_  = 0.1) : 
                               maxA_relative_to_guess(maxA_relative_to_guess_), 
                               minA_relative_to_guess(minA_relative_to_guess_),
                               max_n_df_relative_to_guess(max_n_df_relative_to_guess_) , 
                               phase_start_error(phase_start_error_) , 
                               amp_start_error(amp_start_error_) , 
                               freq_start_error(freq_start_error_)  
                             {;} 

     
             
   }; 

  /** The SineFitter class handles the actual minimization of a set of points to 
   *  a sinusoid. It contains the Minuit2 object, the result, and a FitFn. It may 
   *  be used to fit for a single trace or multiple traces simultaneously (each with different amplitude and phase). 
   **/ 
    class SineFitter
    {

      public: 
        /** Initialize the SineFitter. There are no options. */ 
        SineFitter(); 

        /** Destructor */ 
        virtual ~SineFitter(); 

        /** Set the guess for the fitter. This must be done before doing the fit 
         * @param f The guessed frequency 
         * @param ntrace The number of traces we will next minimize
         * @param ph A pointer to a an array of guess phases, one for each trace. In practice, the phase is difficult to guess. 
         * @param amp the guess amplitude. 
         * */
        void setGuess(double f, int ntrace, const double *ph, double amp); 

        /** Perform the minimation to 1 or more traces. setGuess must be called before. 
         * @param ntrace the number of traces we will minimize
         * @param nsamples the number of samples per trace.
         * @param x x-values of the traces to use. x[0] should be first trace, etc. 
         * @param y y-values of the traces to use. x[0] should be first trace, etc. 
         * @param w weights for each of the traces. This is useful only when fitting multiple traces and you want one to be more important than another.  If 0, all are treated the same. 
         * @param env Envelopes for each of the traces. This is multipled against the sinusoid and and can be used to favor parts of the waveform over another. 
         *
         **/
        void doFit(int ntrace, const int * nsamples, const double ** x, const double **y, const double * w = 0, const double ** env = 0); 

        /** Get the best fit frequency. Only makes sense to call after doFit */ 
        double getFreq() const {return freq;} 

        /** Get pointer to array of best-fit phases. Only makes sense to call after doFit */
        const double* getPhase() const {return &phase[0];} 

        /** Get pointer to array of best-fit amplitudes. Only makes sense to call after doFit */
        const double* getAmp() const {return &amp[0];} 

        /** Get the parameter error on best-fit frequency. Since we don't take into account errors on x and y in the fit, this probably needs to be scaled to be meaninful*/
        double getFreqErr() const {return freq_err;} 

        /** Get the parameter error on best-fit phases. Since we don't take into account errors on x and y in the fit, this probably needs to be scaled to be meaninful*/
        const double *getPhaseErr() const {return &phase_err[0];} 

        /** Get the parameter error on best-fit amplitudes. Since we don't take into account errors on x and y in the fit, this probably needs to be scaled to be meaninful*/
        const double *getAmpErr() const {return &amp_err[0];} 

        /** Get the power (sum(v^2)/n after the fit. This is what is minimized  */
        double getPower() const {return min.MinValue(); } 

        /** Toggle verbose mode, in case looking at screen-fulls of Minuit output is something that excites you */ 
        void setVerbose(bool v) { verbose=v; } 

        /** Return minimzation status */ 
        int getStatus() const { return min.Status(); }

        /** Set the limit parameters */ 
        void setLimitOptions(const SineFitterLimits * lims) { limits = *lims; }

        ROOT::Minuit2::Minuit2Minimizer * minimizer() { return &min; }

      /** 
       * Access the i-th record of doEval in the minimization
       * 
       * @param i is the record to access, default is -1 (will try to find most recent for i < 0)
       * 
       * @return pointer to record graph if it exists, otherwise NULL. DO NOT DELETE
       */

      
        TGraph* getEvalRecordGraph(int i = -1){
          i = i < 0 ? grEvalRecords.size() - 1 : i;
          if(i >= 0 && i < grEvalRecords.size()){
            return grEvalRecords[i];
          }
          return NULL;
        }
        size_t nEvalRecords(){return grEvalRecords.size();}
        void SetDoEvalRecord(bool doEvalRecord){fDoEvalRecord = doEvalRecord;}
        Bool_t GetDoEvalRecord() const {return fDoEvalRecord;}
        void deleteEvalRecords();

      private:
        ROOT::Minuit2::Minuit2Minimizer min; 
        bool verbose; 
        double eval(const double *p); 
        double freq;
        std::vector<double> phase; 
        std::vector<double> amp; 
        double freq_err;
        std::vector<double> phase_err;
        std::vector<double> amp_err;
        bool fDoEvalRecord;
        std::vector<TGraph*> grEvalRecords;
      
        double freq_factor; 
        double min_amp_mean, max_amp_mean;  
        double min_amp_guess, max_amp_guess; 
        SineFitterLimits limits; 



#ifndef __CINT__ //CINT is stupid about subclasses 

        /** The SineFitFn implements a fast version of the objective function and
         * its derivative for minimzation by Minuit2. If ENABLE_VECTORIZE is defined
         * it uses Agner Fog's VCL for manual vectorization, although a lot of the 
         * performance improvement could just be from not using glibc's painfully slow sin(x)
         */
        class SineFitFn : public ROOT::Math::IGradientFunctionMultiDim 
        {
          public:

            /** Initialize. Optional pointer to parent, used to record eval function results. */ 
            SineFitFn(SineFitter* parent = NULL); 

            /** Destroy */
            virtual ~SineFitFn(); 

            /** Feed the FitFn the appropriate data. Note that if SINE_SUBTRACT_USE_FLOATS or SINE_SUBTRACT_FORCE_ALIGNED are defined, these are copied. 
             * Otherwise, we leave them where they are.
             * @param ntrace the number of traces we will minimize
             * @param nsamples the number of samples per trace.
             * @param x x-values of the traces to use. x[0] should be first trace, etc. 
             * @param y y-values of the traces to use. x[0] should be first trace, etc. 
             * @param w weights of each trace, or null to always use the same
             * @param env envelope for each trace, or null to have none (i.e. rectangular envelope)
             */
            void setXY(int ntraces, const int * nsamples, const double **xx, const double **yy, const double * w = 0, const double ** env = 0); 

            /** This evaluates the objective function for a parameter vector p */ 
            virtual double DoEval(const double *p) const; 

            /** This evaluates the coord-th derivative of the objective function for a parameter vector p */
            virtual double DoDerivative(const double *p, unsigned int coord) const; 

            /** Returns the number of dimensoins */ 
            unsigned NDim() const { return 1+2*nt; } 

            /** This needs to be implemented to be a proper subclass */ 
            virtual ROOT::Math::IBaseFunctionMultiDim * Clone() const;
          
          private:
            SineFitter* fContainer; // pointer to parent
          
            int *ns;
            int nt; 
#ifdef SINE_SUBTRACT_USE_FLOATS
            float **x; 
            float **y; 
            float **env; 
            const double **xp, **yp, ***envp; 
#elif defined(SINE_SUBTRACT_FORCE_ALIGNED)
            double **x; 
            double **y; 
            double **env; 
            const double **xp, **yp, **envp; 
#else
            const double ** x; 
            const double ** y; 
            const double ** env; 
#endif 
            const double * wgt; 

          
        } f;
#endif

    }; 

    /* This class stores the minimization result. It would be a struct if CINT weren't
     * stupid about structs */
    class SineSubtractResult
    {

      public: 

     /** Empty the bowels. */
      void clear(); 

      /** Add the result of another iteration */ 
      void append(const SineSubtractResult * r); 

      /** These destructors, they do nothing */ 
      virtual ~SineSubtractResult(){;} 


      /** Stores the power at each step. This vector will be bigger than the others since it includes the initial power.*/ 
      std::vector<double> powers; 

      /** Stores the fit phases at each step */ 
      std::vector<std::vector<double> >phases; 

      /** Stores the fit frequency at each step */ 
      std::vector<double> freqs; 

      /** Stores the fit amps at each step */ 
      std::vector<std::vector<double> > amps; 

      /** Stores the parameter errors on phases at each step */ 
      std::vector<std::vector<double > > phases_errs; 

      /** Stores the parameter error on frequency at each step */
      std::vector<double> freqs_errs; 

      /** Stores the parameter errors on amps at each step */ 
      std::vector<std::vector<double> >amps_errs; 


      ClassDef(SineSubtractResult,3); 
    }; 


    /** The SineSubtract class  handles the number of iterations
     * and initialization of the fitter. 
     */ 
    class SineSubtract 
    {

      public: 

        /** Create a new SineSubtract. A few options may be set here, although they can all be changed by setters later too 
         * (and there are other options that can only be set by setters. 
         * 
         *  @param max_iter_without_reduction  Set the maximum number of failed iterations before giving up. An iteration fails if the subtraction does not produce the required power reduction.
         *  @param min_power_reduction Set the minimum power reduction necessary to keep an iteration. Default is 5 percent. 
         *  @param store If true, intermediate steps will be saved. This incurs some overhead, so is mostly useful for making pretty plots and diagnostics. 
         *
         */
        SineSubtract(int max_iter_without_reduction = 3, double min_power_reduction = 0.05, bool store = false); 


        /** Create a new SineSubtract. A few options may be set here, although they can all be changed by setters later too 
         * (and there are other options that can only be set by setters. 
         * 
         *  @param freq_dependent_min_power_reductionmin_power_reduction Set the minimum power reduction necessary to keep an iteration, varying with frequency (freq in GHz). 
         *  @param max_iter_without_reduction  Set the maximum number of failed iterations before giving up. An iteration fails if the subtraction does not produce the required power reduction.
         *  @param store If true, intermediate steps will be saved. This incurs some overhead, so is mostly useful for making pretty plots and diagnostics. 
         *
         */

        SineSubtract(const TGraph * freq_dependent_min_power_reduction, int max_iter_without_reduction = 3, bool store = false); 


        /** Deallocate everything */ 
        virtual ~SineSubtract(); 

        /* Subtract CW from a single trace. Sine Subtraction can handle both evenly spaced and unevely spaced waveforms. 
         * For evenly spaced waveforms, you should pass dt <= 0, otherwise you should pass the nominal average. 
         * 
         * The input graph is not touched. A new graph is allocated for the subtraction. 
         *
         * @param g The input TGraph to subtract from 
         * @param dt The nominal sample rate for uneven waveforms. If <=0, then the graph is assumed to be even and dt is computed from the first time step. 
         * @param result if you already have the result, this just performs the transformations in the result 
         * @return The subtracted Graph.  
         */
        TGraph * subtractCW(const TGraph * g, double dt, const SineSubtractResult* result = NULL);

        /* Convenience method to subtract CW from an evenly-sampled array 
         *  (note that for now this just makes a TGraph copy internally, in the future this will be reversed for efficiency); 
         *  
         *  @param N number of samples
         *  @param y input y-values
         *  @param dt sampling period
         *  @param result if you already have the result, this just performs the transformations in the result 
         *  @returns waveform with subtracted CW. Either allocates new or uses y_out, depending on if y_out is NULL
         * */ 
        double * subtractCW(int N, const double * y, double dt, double * y_out = NULL, const SineSubtractResult *result = NULL ); 


        /** Subtract CW from one or more traces. Sine Subtraction can handle both evenly-spaced and unevenly-spaced waveforms.
         *  For evenly spaced waveforms, you should pass dt<-0, otherwise you should pass the nominal average. 
         *  Right now, all graphs are assumed to have the same sampling, although this limitation could be removed with some work. 
         *
         *  The input graphs are touched, you should make a copy if you want to preserve the originals. 
         *
         *  @param ng the number of graphs
         *  @param g The graphs. Will be modified. 
         *  @param dt The nominal sample rate for uneven waveforms. If <=0, then the graphs are assumed to be even and dt is computed from first time step. 
         *  @param w Scales for the input values. If the y-axis have different scales (different gains or units) you should pass an array here. If 0, everything is equally-weighted. 
         *
         */
        void subtractCW(int ng, TGraph ** g, double dt, const double * w = 0, const SineSubtractResult* result = NULL); 

        /** Set limits on the frequencies to try to subtract. If the units of the graph are in ns, the frequencies should be in GHz. 
         *
         * @param min The minimum frequency to try to remove
         * @param max The maximum frequency to try to remove
         */ 
        void setFreqLimits(double min, double max) { fmin.clear(); fmin.push_back(min); fmax.clear(); fmax.push_back(max); } 

        /** Set multiple allowed bands of frequencies to try to subtract. If the units of the graph are in ns, the frequencies should be in GHz. 
         *
         * @param nbands The number of frequency bands
         * @param min Array of minima of each band
         * @param max Array of maxima of each band
         */
        void setFreqLimits(int nbands, const double * min, const double * max) { fmin.clear(); fmax.clear(); for (int i = 0; i < nbands;i++) { fmin.push_back(min[i]); fmax.push_back(max[i]); }; } 

        /** Clear all frequency limits */ 
        void unsetFreqLimits() { fmin.clear(); fmax.clear(); high_factor = 1;  } 


        /** Toggle voluminous verbosity. It is a well-known fact that additional lines of text 
         * scrolling with great alacrity on your terminal have myriad benefits that this comment 
         * is far too brief to expound. */
        void setVerbose(bool v) {verbose = v; fitter.setVerbose(v); }


        /** In some cases, it may make sense to use only part of the trace for subtraction. For example, one may want to avoid
         * using the part near the trigger to avoid subtracting off a pulse. If trace limits are enabled, only the portion 
         * of the trace within limits is used for estimating the spectrum and finding sinusoids that reduce the power, however,
         * once those sinusoids are found, the entire trace has them subtracted. 
         *
         * While in principle, these could be different for each each waveform when subtracting multiple, right now
         * they are all the same. 
         *
         * @param min The minimum sample of the trace to use
         * @param max The maximum sample of the trace to use
         **/ 
        void setTraceLimits(int min, int max) { tmin = min; tmax = max; }



        /** 
         *  \enum PowerSpectrumEstimator
         *
         *  These are the algorithms for estimating the power spectrum, used to help estimate
         *  the starting parameters. 
         */
        enum PowerSpectrumEstimator
        {
          FFT, /// (default) Compute power spectrum using magnigude of FFT, see FFT_Param enum for parameter names and defaults
          LOMBSCARGLE /// Compute power spectrum using LombScargle Periodogram, see below for parameter names
        }; 


        /** The parameters for PowerSpectrumEstimator FFT */ 
        enum FFT_Params
        {
          FFT_NPAD, /// default =0 
          FFT_NPARAMS
        }; 

        /** The parameters for PowerSpectrumEstimator LOMBSCARGLE*/ 
        enum LOMBSCARGLE_Params
        {
          LS_OVERSAMPLE_FACTOR, ///default =2
          LS_NPARAMS
        }; 


        /** Each time an iteration fails, the power spectrum bin is weighted by
         * 1./(1+nfails)^exp to discourage further usage */ 
        void setNFailsExponent(double exp = 0.5) { nfail_exponent = exp ; }


        /** Set the power spectrum estimator used to initialize the SineSubtractFilter. 
         *  See  documentation for PowerSpectrumEstimator for what the types mean (or pass 0 to use the defaults). 
         **/ 
        void setPowerSpectrumEstimator(PowerSpectrumEstimator estimator, const double * params = 0); 


        /** 
         * \enum PeakFindingOption 
         *
         * These are the algorithms for finding the peak of the power spectrum.
         **/
        enum PeakFindingOption
        {
          GLOBALMAX, /// Just uses the global maximum (you probably don't want this) . No params. 
          NEIGHBORFACTOR, /// (default) requires being at least neighbor_factor bigger than one of the neighbors. See NEIGHBORFACTOR_Params for parametser
          TSPECTRUM, ///Use TSpectrum (probably super slow). See TSPECTRUM_Params for parameters. 
          SAVGOLSUB /// Subtract baseline computed with Savitzky Golay filter first. See SAVGOLSUB_Params for parameters. 
        }; 


        enum NEIGHBORFACTOR_Params
        {
          NF_NEIGHBOR_FACTOR, ///default = 0.15 
          NF_NPARAMS
        }; 

        /* These are basically based on the 
         * options to TSpectrum::searchHighRes
         */ 
        enum TSPECTRUM_Params
        {
          TS_SIGMA, /// default = 2
          TS_THRESHOLD, /// default = 0.05
          TS_NDECONV_ITERATIONS, //default = 3
          TS_AVERAGE_WINDOW, //default = 3
          TS_NPARAMS
        }; 

        enum SAVGOLSUB_Params
        {
          SGS_ORDER, ///default = 3 (treated as int) 
          SGS_WIDTH, ///default = 5 , really halfwidth, (treated as int) 
          SGS_NPARAMS
        }; 

        /** Set the peak finding option used to seed the sinusoid fitter. See PeakFindingOption for possibilities. 
         *  params depends on the particular peak finder you want (or pass 0 for defaults).
         **/ 

        void setPeakFindingOption(PeakFindingOption peak, const double * params = 0); 

        enum EnvelopeOption
        {
          ENV_NONE, /// No envelope is used
          ENV_HILBERT, //fit to hilbert envelope
          ENV_RMS, //fit to sliding RMS envelope 
          ENV_PEAK //fit to sliding RMS envelope 
        };  

        enum ENV_HILBERT_Params
        {
          ENV_HILBERT_FIT_ORDER, //default = 3, if < 0, no fit (just take values) 
          ENV_HILBERT_NPARAMS

        }; 

        enum ENV_RMS_Params
        {
          ENV_RMS_FIT_ORDER, //default = 3, if < 0, no fit (just take envelope values) 
          ENV_RMS_WINDOW, //default = 5
          ENV_RMS_NPARAMS
        }; 

       enum ENV_PEAK_Params
        {
          ENV_PEAK_FIT_ORDER, //default = 3, if < 0, no fit (just take envelope values) 
          ENV_PEAK_MINDISTANCE, //default = 5 (appropriate for ~200 MHz min freq) 
          ENV_PEAK_NPARAMS
        }; 

        /** Set the envelope option used to estimate an envelope to modify the sinusoids. See EnvelopeOption for possibilities. 
         *  params depends on the particular envelope option you want (or pass 0 for defaults).
         **/ 
        void setEnvelopeOption(EnvelopeOption env, const double * params = 0); 

        /** Returns a pointer to a vector of the sequence of powers vs. iteration */ 
        const std::vector<double> * getPowerSequence() const { return &r.powers; } 

        /** If storage of graphs is enabled, returns a pointer to a nested vector of stored graphs */ 
        const std::vector<std::vector<TGraph*> > * getStoredGraphs() const { return &gs; }

        /** If storage of graphs is enabled, returns a pointer to a nested vector of stored power spectra */ 
        const std::vector<TGraph*> * getStoredSpectra() const { return &spectra; }


        /** If storage of graphs is enabled, returns the number of stored graphs per channel (basically the number of iterations+1). Will be 0 if storage disabled. */ 
        size_t nStoredGraphsInChannel() const { return gs[0].size(); } 

        /** If storage of graphs is enabled, returns the number of channels stored (the number of traces that are being simultaneously fit for. Will be 0 if storage disabled. */ 
        size_t nChannels() const { return gs.size(); } 

        /** Returns the ith graph in the cth trace */ 
        TGraph* storedGraph(int i, int c) { return gs[c][i]; } 
        TGraph* storedEnvelope(int i, int c) { return env_gs[c][i]; } 

        /** If storage of graphs is enabled, returns the number of stored spectra. (Basically the number of iterations+1). Will be 0 if storage disabled. */ 
        size_t nStoredSpectra() const { return spectra.size(); } 

        /** Returns the power spectra of the ith iteration */
        TGraph* storedSpectra(int i) { return spectra[i]; } 

        /** Return the number of sines fit */ 
        int getNSines() const; 
        
        /* Return the number of traces fit for */
        int getNGraphs() const { return r.phases.size(); } 

        /** Get array of fit phases for the ith trace */
        const double * getPhases(int g) const { return &r.phases[g][0]; }

        /** Get array of fit frequencies*/
        const double * getFreqs() const { return &r.freqs[0]; }

        /** Get array of fit amplitudes for the ith trace */
        const double * getAmps(int g) const { return &r.amps[g][0]; }

        /** Get array of fit phase errors for the ith trace */
        const double * getPhaseErrs(int g) const { return &r.phases_errs[g][0]; }

        /** Get array of fit frequnecy errors */
        const double * getFreqErrs() const { return &r.freqs_errs[0]; }

        /** Get array of fit amplitude errors for the ith trace */
        const double * getAmpErrs(int g) const { return &r.amps_errs[g][0]; }

        /** Get a pointer to the SineSubtract result */
        const SineSubtractResult * getResult() const { return &r; } 

        /* Set the maximum number of failed iterations before giving up. */
        void setMaxIterationsWithoutReduction(int max) { maxiter = max;  } 

        /** Sets the maximum absolute number of iterations, even if power reduction is still happening. If max <= 0, there is no max */ 
        void setAbsoluteMaxIterations(int max) { abs_maxiter = max; } 

        /** Sets the maximum number of successful iterations, useful if you're trying to subtract a particular sinusoid. If max <=0, there is no max*/ 
        void setMaxSuccessfulIterations(int max) { max_successful_iter = max; } 

        /** Set the threshold for considering an iteration successful. This is a percentage of power reduction. */
        void setMinPowerReduction(double min) {min_power_reduction = min; g_min_power = 0;}; 

        /* Set the frequency-dependent threshold before giving up. */
        void setMinPowerReduction(const TGraph *g) { g_min_power= g; }  


        /** Togggle storage of intermediate traces */
        void setStore(bool s) { store = s; }

        /** Makes the fitter use this SineFitterLimits */
        void setFitterLimitOptions( const SineFitterLimits * lims) { fitter.setLimitOptions(lims); }

        //these are super auxilliary and could be moved out (they only need public methods anyway) 
        
        /**  Generate plots of traces and spectra.  
         * @param cpower Pointer to canvas to put traces + fits on. If NULL, new one is created. 
         * @param cspectra Pointer to canvas to put spectra on. If NULL, new one is created. 
         * @param ncols The number of columns per canvas  
         **/ 
        void makePlots(TCanvas * cpower = 0,TCanvas *cspectra = 0, int ncols = 4) const; 

        /**  Make Beamer slides showing all iterations, for those moments when you desperately need slides to show. 
         *
         * @param title Used in the frame title
         * @param prefix Used for naming files (e.g. sinsub will produce sinsub_0, sinsub_1, etc )
         * @param outdir Output directory
         * @param format The format of the plots (e.g. png, pdf) 
         * @param standalone If true, makes a self contained presentation, otherwise, can be copied and pasted in. 
         *
         * */ 
        void makeSlides(const char * title = "SineSubtract", const char * prefix = "sinsub", const char * outdir = ".", const char* format = "png", bool standalone = true) const; 

        SineFitter * sineFitter() { return &fitter; } 

      private: 
        int abs_maxiter; 
        int maxiter; 
        double nfail_exponent; 
        int max_successful_iter; 
        int tmin, tmax; 
        std::vector<double> fmin; 
        std::vector<double> fmax; 
        double min_power_reduction; 
        std::vector<std::vector<TGraph*> > gs; 
        std::vector<std::vector<TGraph*> > env_gs; 
        std::vector<TGraph*> spectra;
        SineSubtractResult r; 
        bool store; 
        double high_factor; 

        PowerSpectrumEstimator power_estimator;
        double * power_estimator_params; 
        PeakFindingOption peak_option; 
        double * peak_option_params; 


        EnvelopeOption envelope_option; 
        double * envelope_option_params; 

        const TGraph * g_min_power; 

        bool verbose; 
        void reset(); 
        unsigned id; 

        /** some utility stuff */ 
        int findMaxFreq(int Nfreq, const double * freq, const double * mag, const int * nfails) const; 
        bool allowedFreq(double f, double df) const; 

        SineFitter fitter; 
    }; 

}

#endif
#endif 
