#ifndef _SINE_SUBTRACT_H
#define _SINE_SUBTRACT_H

/**
 * \file SineSubtract.h 
 *
 *   Cosmin Deaconu <cozzyd@kicp.uchicago.edu> 
 *
 * Independent implementation of sine subtraction suggested by Andres/ Steph
 * with some differences. 
 *
 * SineSubtract is controlling class, SineFitter does the minimization,
 * SineFitFn provides the objective function.. 
 *
 *
 * The part of the waveform used by sine subtraction may be limited (e.g. if you wanted to just use the pretrigger portion). 
 * Also, multiple waveforms may be fit for at once (a single frequency is fit for in each iteration, with the amplitude and phase allowed to vary per waveform). 
 *
 * Algorithm is: 
 *
 *    -Rectify waveform and estimate power spectrum using lomb-scargle periodogram 
 *    -Find peak, fmax, in power spectrum by finding magnitude^2/(1+nfails) at
 *    least neighbor_factor^2 bigger than a neighboring magnitude^2/(1+nfails).
 *    (see below for meaning of nfails)
 *
 *    -Find sine ( A (sin wt + phi)) that when subtracted minimizes "power"
 *    (avg(v_i^2)); initialize with w = 2pi fmax, A = mag, phase from fft.
 *    Uninterpolated graph is used for this! 
 *
 *    -If power reduced by more than min_power_reduction, subtract 
 *
 *    -If power not reduced by more than min_power reduction, increment nfails
 *    for that bin 
 *
 *    -Repeat until max_iter_without_reduction iterations with not enough
 *    improvement   
 *
 * 
 *  If store is set to true, the graphs and power spectra at each step are
 *  stored so that they may be visually inspected (or using makePlots)
 *
 *
 *  SineFitter uses Minuit2 to minimize sum_i ( (g[i] - A (sin 2 pi f t +
 *  phase)) ^2)/n. The gradient is computed analytically. 
 *
 *  If ENABLE_VECTORIZE is defined, custom vectorized code (using Agner
 *  Fog's VCL, http://www.agner.org/optimize/#vectorclass) is used, which is
 *  significantly faster until compilers get smart enough to autovectorize the
 *  characteristic function and gradient. This only makes sense in a
 *  processor with SIMD. 256-bit floating point registers are used (e.g. AV2X)
 *  but the VCL will gracefully fallback to 128-bit registers at compile time
 *  should your processor suck. 
 *
 *  **WARNING**: You will not get the same "answers" in vectorized mode. This
 *  is expected.  Reasons for this are thought to include: 
 *
 *    - Horizontal addition (which breaks associativity, and can be disabled by
 *    defining SINE_SUBTRACT_DONT_HORIZONTAL_ADD) 
 *
 *    - Different implementations of e.g. sin (VCL appears to use a version
 *    based on libcephes, which is different from recent versions of 64-bit
 *    Linux) 
 *
 *    - FMA instructions may not be enabled, depending on compiler options, in
 *    non-vectorized mode 
 *
 *
 *  If SINE_SUBTRACT_USE_FLOATS is defined, the fitter will use floats instead
 *  of doubles. With vectorization, this should be roughly twice as fast (due
 *  to more floats fitting in registers). 
 *
 *  If SINE_SUBTRACT_PROFILE is defined, the time spent calling subtractCW is
 *  printed. 
 *
 *  See macros/testArtificialSubtract.C to see SineSubtract in action.  
 *
 *
 * TODO list: 
 *
 *  - Better peak finding algorithm (maybe smooth + non-maximum suppression or
 *    a wavelet-based peak-finder)  
 *
 *  - Better tuning of step sizes / limits ?  - Tuning of nfails scaling?  
 *
 */ 

#include "Minuit2/Minuit2Minimizer.h" 
#include "FFTtools.h" 
#include "TObject.h"

#include <vector>

class TGraph; 
class TCanvas; 

namespace FFTtools
{
    class SineFitter
    {

      public: 
        SineFitter(); 
        virtual ~SineFitter(); 
        void setGuess(double f, int ntrace, const double *ph, const double *amp); 
        void doFit(int ntrace, int nsamples, const double ** x, const double **y); 
        double getFreq() const {return freq;} 
        const double* getPhase() const {return &phase[0];} 
        const double* getAmp() const {return &amp[0];} 

        double getFreqErr() const {return freq_err;} 
        const double *getPhaseErr() const {return &phase_err[0];} 
        const double *getAmpErr() const {return &amp_err[0];} 

        double getPower() const {return min.MinValue(); } //"power" (sum (V^2)/n) 
        void setVerbose(bool v) { verbose=v; } 
        int getStatus() const { return min.Status(); }


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


#ifndef __CINT__ //CINT is stupid about subclasses 




        class SineFitFn : public ROOT::Math::IGradientFunctionMultiDim 
        {
          public:

            SineFitFn() {ns = 0;  nt = 0; x =0; y =0; }
            void setXY(int ntraces, int nsamples, const double **xx, const double **yy); 
            virtual double DoEval(const double *p) const; 
            virtual double DoDerivative(const double *p, unsigned int coord) const; 
            unsigned NDim() const { return 1+2*nt; } 
            virtual IBaseFunctionMultiDim * Clone() const; 


          private:
            int ns; 
            int nt; 
#ifdef SINE_SUBTRACT_USE_FLOATS
            float **x; 
            float **y; 
            const double **xp, **yp; 
#else
            const double ** x; 
            const double ** y; 
#endif 

        } f; 
#endif


    }; 

    class SineSubtractResult
    {

      public: 

      void clear(); 
      void append(const SineSubtractResult * r); 
      virtual ~SineSubtractResult(){;} 

      std::vector<double> powers; 
      std::vector<std::vector<double> >phases; 
      std::vector<double> freqs; 
      std::vector<std::vector<double> > amps; 
      std::vector<std::vector<double > > phases_errs; 
      std::vector<double> freqs_errs; 
      std::vector<std::vector<double> >amps_errs; 


      ClassDef(SineSubtractResult,3); 
    }; 


    class SineSubtract 
    {

      public: 
        SineSubtract(int max_iter_without_reduction = 3, double min_power_reduction = 0.01, bool store = false); 
        virtual ~SineSubtract(); 

        /* Subtract CW from a single trace. will be interpolated to dt if dt > 0 for power spectrum. */ 

        TGraph * subtractCW(const TGraph * g, double dt); 

        void subtractCW(int ng, TGraph ** g, double dt); 
        void unsetFreqLimits() { fmin.clear(); fmax.clear(); } 
        void setFreqLimits(double min, double max) { fmin.clear(); fmin.push_back(min); fmax.clear(); fmax.push_back(max); } 
        void setFreqLimits(int nbands, const double * min, const double * max) { fmin.clear(); fmax.clear(); for (int i = 0; i < nbands;i++) { fmin.push_back(min[i]); fmax.push_back(max[i]); }; } 
        void setOversampleFactor(double of) {oversample_factor = of;}
        void setHighFactor(double hf) {high_factor = hf;}

        void setVerbose(bool v) {verbose = v; fitter.setVerbose(v); }

        void setTraceLimits(int min, int max) { tmin = min; tmax = max; }


        const std::vector<double> * getPowerSequence() const { return &r.powers; } 

        const std::vector<std::vector<TGraph*> > * getStoredGraphs() const { return &gs; }

        const std::vector<TGraph*> * getStoredSpectra() const { return &spectra; }
        const std::vector<TGraph*> * getStoredFFTPhases() const { return &spectra; }

        size_t nStoredGraphsInChannel() const { return gs[0].size(); } 
        size_t nChannels() const { return gs.size(); } 
        TGraph* storedGraph(int i, int c) { return gs[c][i]; } 
        size_t nStoredSpectra() const { return spectra.size(); } 
        TGraph* storedSpectra(int i) { return spectra[i]; } 

        int getNSines() const; 
        int getNGraphs() const { return r.phases.size(); } 
        const double * getPhases(int g) const { return &r.phases[g][0]; }
        const double * getFreqs() const { return &r.freqs[0]; }
        const double * getAmps(int g) const { return &r.amps[g][0]; }
        const double * getPhaseErrs(int g) const { return &r.phases_errs[g][0]; }
        const double * getFreqErrs() const { return &r.freqs_errs[0]; }
        const double * getAmpErrs(int g) const { return &r.amps_errs[g][0]; }
        const SineSubtractResult * getResult() const { return &r; } 

        void setMaxIterationsWithoutReduction(int max) { maxiter = max; } 
        void setMinPowerReduction(double min) {min_power_reduction = min; }; 
        void setStore(bool s) { store = s; }
        void setNeighborFactor(double n) { neighbor_factor2 = n*n; } 

        //these are super auxilliary and could be moved out (they only need public methods anyway) 
        void makePlots(TCanvas * cpower = 0,TCanvas *cspectra = 0, int ncols = 4) const; 
        void makeSlides(const char * title = "SineSubtract", const char * prefix = "sinsub", const char * outdir = ".", const char* format = "png") const; 

      private: 
        int maxiter; 
        int tmin, tmax; 
        std::vector<double> fmin; 
        std::vector<double> fmax; 
        double min_power_reduction; 
        std::vector<std::vector<TGraph*> > gs; 
        std::vector<TGraph*> spectra; 
        SineSubtractResult r; 
        SineFitter fitter; 
        bool store; 
        double neighbor_factor2; 
        double oversample_factor ; 
        double high_factor; 

        bool verbose; 
        void reset(); 

    }; 

}

#endif 
