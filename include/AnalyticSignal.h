#ifndef FFTTOOLS_ANALYTIC_SIGNAL_H
#define FFTTOOLS_ANALYTIC_SIGNAL_H


#include <vector>
#include "TGraph.h" 

class TRandom; 

namespace FFTtools
{

    class FFTWindowType; 

    class AnalyticSignal 
    {
      public:
        AnalyticSignal(double extra_noise = 0) { setExtraNoise(extra_noise); isSetup = false; }
        virtual double eval(double t); 
        /* fill graph with signal */ 
        void fillWf(TGraph * g, bool zero_y = false); 


        /* Generate evenly sampled waveform with sample width dt */ 
        virtual TGraph * makeEvenWf(int N, double dt, double t0 = 0); 

        /* Generate waveform with Gaussian distributed dt */ 
        virtual TGraph * makeGaussianDtWf(int N, double mean_dt, double dt_sigma, double t0 = 0); 

        /* Generate waveform with alternating Gaussian distributed dt, separated by delta_t (like the DRS 4) */ 
        virtual TGraph * makeAlternatingGaussianDtWf(int N, double mean_dt, double delta_dt, double dt_sigma, double t0 = 0); 

        /* Generate waveform where each sample time is shifted from its evenly sampled location by a value generated uniformly by +/- max_jitter_frac */
        virtual TGraph * makeUniformJitterWf(int N, double dt, double max_jitter_frac = 0.25, double t0 = 0); 

        /* Generate waveform where each sample time is shifted from its evenly sampled location gaussianly  */
        TGraph * makeGaussianJitterWf(int N, double dt, double jitter_sigma, double t0 = 0); 

        void setExtraNoise(double val) { extra_noise = val; } 
        virtual void eval(int N, const double* t, double * v, bool zero_v = false); 
        void setRNG(TRandom  * rand) { rng = rand; } 

        virtual void setup() {;}
        virtual ~AnalyticSignal() {; }


      protected:
        TRandom * rng; 
        double extra_noise; 
        virtual double doEval(double t) {(void) t; return 0; }
        bool isSetup; 

    }; 

    class CompositeSignal : public AnalyticSignal 
    {

      public: 

        CompositeSignal(AnalyticSignal * sig1 = 0, AnalyticSignal * sig2 = 0, bool take_ownership = true) 
        {
          if (sig1) add(sig1,take_ownership); 
          if (sig2) add(sig2,take_ownership); 
          rng = 0; isSetup = false; extra_noise = 0; 
        }

        void add(AnalyticSignal * sig, bool take_ownership) { signals.push_back(sig); own.push_back(take_ownership);} 
        ~CompositeSignal(); 
        AnalyticSignal * getSignal(int i) { return signals[i]; }
        size_t nSignals() const { return signals.size(); }

        virtual void setup(); 
   
      protected:
        virtual double doEval(double t); 
      private:
        std::vector<AnalyticSignal *> signals; 
        std::vector<bool> own; 


    }; 


    class BandlimitedSampledSignal : public AnalyticSignal 
    {

      public: 
        BandlimitedSampledSignal(size_t n, const double *y, double dt, double t0 = 0, int nlobes = 0, const FFTWindowType * win = 0); 
        virtual void setup() {;}
      protected:
        TGraph g; 
        int max_lobes; 
        const FFTWindowType *win; 
        virtual double doEval(double t); 
    }; 



    class ThermalNoise : public AnalyticSignal 
    {

      public: 
        // generate ``analytic'' thermal noise
        // nosc is effective number of phasors (adjusted by bandwidth) 
        // fmin and fmax are frequency limits
        // rms is approximate RMS (might be different because of filtering) 
        // if order <= 0, a hard frequency cutoff is used, otherwise butterworth-like rolloff 
        ThermalNoise(int nosc, double fmin, double fmax, double rms, int order = 0);
        virtual ~ThermalNoise() { ; }
        void setup(); 
      protected:

        virtual double doEval(double t); 
      private:
        int nosc; 
        double fmin, fmax, rms; 
        std::vector<double> A;
        std::vector<double> phases;
        std::vector<double> freqs;
        int order; 


    }; 
}



#endif 
