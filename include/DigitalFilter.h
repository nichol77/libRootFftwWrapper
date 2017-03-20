#ifndef FFTTOOLS_DIGITAL_FILTER_HH
#define FFTTOOLS_DIGITAL_FILTER_HH

class TGraph; 
class TPad; 
#include "TString.h" 
#include <vector>
#include <complex>
#include <cstdlib>

/** 
 *
 * Digital Filter Implementation.
 * 
 *
 * Cosmin Deaconu <cozzyd@kicp.uchicago.edu> 
 *
 * */ 

namespace FFTtools
{
    class FFTWindowType; 

    /** Base class for all digital filters. Any implementing filters must implement filterOut and transfer.
    */
    class DigitalFilter
    {
     public: 

        /** 
         * Filter input, returning a newly allocated output. The caller is responsible for freeing the allocated memory.
         *
         * @param n size of waveform 
         * @param w waveform to filter
         * @returns newly allocated output 
         * */ 
        virtual double * filter(size_t n, const double * w) const; 

        /** Filter input with previously allocated output. 
         *
         * @param n size of waveform to filter
         * @param w input waveform
         * @param out output waveform (assumed to be allocated) 
         *
         *
         * */ 
        virtual void filterOut(size_t n, const double * w, double *out) const = 0; 

        /** Filter graph (and y errors, if available) ``in place'' (with temporary allocation)
         *
         * @param g graph to filter
         * @param filterErrors true if g->GetEY() should be filtered (if it exists) 
         *
         * */ 
        virtual void filterGraph(TGraph * g, bool filterErrors = false) const; 

        /**  Filter ``in place'' (with temporary allocation) 
         *
         * @param n size of waveform
         * @param w input waveform, output will also be written here 
         *
         * */
        virtual void filterReplace(size_t n, double * w) const; 


        /** Filter response to unit impulse. A delay may be given which is helpful for acausal filters. 
         * This works by applying filterOut on an input array with all zeroes and one one at position delay.  
         * 
         * @param n number of samples wanted in impulse response
         * @param out allocated memory where impulse response will be written
         * @param delay the time (in samples) of the unit impulse 
         * */ 
        virtual void impulse(size_t n, double * out, size_t delay = 0) const; 


        /** Filter response to unit impulse that allocates new memory. A delay
         * may be given which is helpful for acausal filters. The user is
         * responsible for deleting te allocated memory.   
         *
         * @param n number of samples wanted in impulse response 
         * @param delay the time (in samples) of the unit impulse 
         * @returns a newly allocated array of the impulse response
         *
         **/ 
        virtual double * impulse(size_t n, size_t delay = 0) const; 


        /** Get impulse response as a newly allocated TGraph. Caller is responsbile for deletion. 
         *
         * @param n number of samples wanted in impulse response
         * @param dt the sample rate of the graph
         * @param delay the ttime (in samples) of the unit impulse 
         * @returns a TGraph of the impulse response
         *
         **/ 
        TGraph * impulseGraph(size_t n  = 101, double dt = 1, size_t delay = 50) const; 


        /** Compute amplitude response, phase response and group delay directly from the transfer function, with all the gory ensuing numerical problems.
         * This method is provided becuase it is more efficient to calculate two (or three) of the responses at once than individually. 
         *
         * If any is 0, it is not computed. 
         * If any TGraph ** points to a TGraph * that points to 0, TGraph of the right size is allocated
         *
         */
        virtual void response(size_t n, TGraph ** amplitude_response, TGraph ** phase_response, TGraph ** group_delay = 0) const; 

        /* return amplitude response with n points */ 
        virtual TGraph* amplitudeResponse(size_t n = 101) const  {TGraph * amp = 0; response(n,&amp,0); return amp; }

        /* return phase response response with n points */ 
        virtual TGraph* phaseResponse(size_t n = 101) const { TGraph *ph = 0; response(n,0,&ph); return ph; }
        virtual TGraph * groupDelay(size_t n = 101) const { TGraph *gd = 0; response(n,0,0,&gd); return gd; }

        /** Draws the amplitude, phase, group delay and impulse response into the pad. Returns the pad.  */ 
        virtual TPad* drawResponse(TPad * c = 0, int n = 101, int delay = 50) const; 

        /* Computes transfer function */ 
        virtual std::complex<double> transfer(std::complex<double> z) const = 0;  

        virtual ~DigitalFilter() {;} 
    }; 


    /* A series of filters, computed by serially applying the filters (and multiplying the transfer functions) */ 
    class DigitalFilterSeries : public DigitalFilter
    {


      public: 

        /* Initialize filter series with a filter */ 
        DigitalFilterSeries(const DigitalFilter * f, bool takeOwnership = false) { add(f,takeOwnership); } 

        /* Empty filter series */ 
        DigitalFilterSeries() {; } 


        virtual void filterOut(size_t n, const double *w, double *out) const; 
        virtual std::complex<double> transfer(std::complex<double> z) const;  

        /* Add a filter to the series. Does NOT take ownership of it */ 
        virtual void add(const DigitalFilter *f, bool o) { series.push_back(f); own.push_back(o); }
        virtual void clear() 
        {
         for (size_t i =0; i< series.size(); i++) 
         { 
           if (own[i]) delete series[i]; 
         } 
         series.clear(); 
         own.clear(); 
        }
        virtual ~DigitalFilterSeries() {clear(); }

      protected:
        std::vector<const DigitalFilter*> series; 
        std::vector<bool> own; 

    }; 

    /* FIR filter
     *
     * This is effectively a convolution. This implementation has some non-standard options
     * that make an FIR filter NOT a subset of an IIR filter (Delay and extend). 
     *
     * */
    class FIRFilter : public DigitalFilter
    {
      public:
        //FIR filter of length N with values x. if extend is true, input to filter will be extended by copying first and last values N times before and after
        FIRFilter(size_t N, int delay = 0, bool extend =false) : coeffs(N),delay(delay),extend(extend) {; }
        FIRFilter(size_t N, const double * x, int delay = 0, bool extend = false) : coeffs(x,x+N), delay(delay), extend(extend) {; }
        virtual void filterOut(size_t n, const double * w, double * out) const; 
        virtual std::complex<double> transfer(std::complex<double> z) const ;  
        virtual void setDelay(int d) { delay = d; } 
        virtual void setExtend(bool ext) { extend = ext; } 
        size_t getOrder() const { return coeffs.size(); }
        int getDelay() const { return delay; } 
        const double* getCoeffs() const { return &(coeffs[0]); } 

      protected: 
        std::vector<double> coeffs; 
        int delay; 
        bool extend; 
    }; 

    /** A Savitzky-Golay Filter effectively fits an order n polynomial to 1+nleft+nright points, acting as a sort of
     * smoothing filter.  */ 

    class SavitzkyGolayFilter : public FIRFilter
    {
      public:
        //width_right < 0 means width_right = width_left 
        SavitzkyGolayFilter(int polynomial_order, int width_left, int width_right=-1, int derivative = 0);  
    }; 


    class GaussianFilter : public FIRFilter
    {
      public: 
        GaussianFilter(double sigma, double nsigma); 
    }; 

    class BoxFilter : public FIRFilter
    {
      public: 
        BoxFilter(int width); 
    }; 

    class DifferenceFilter : public FIRFilter
    {
      public: 
        DifferenceFilter(int order =1); 
    }; 


    class SincFilter : public FIRFilter 
    {
      public: 
        SincFilter(double w, int max_lobes, const FFTWindowType * win = 0, int delay = 0, bool extend = false); 
    }; 

    /** IIR filter implementation 
     *
     *  An IIR Filter is a rational function of z^-1 (which, of course, is e^(-jw), the unit delay) ) 
     *
     *           b0 + b1 z^-1 + b2 z^-2 + ... 
     *  H(z) =  -------------------------------
     *           a0 + a1 z^-1 + a2 z^-2 + ... 
     *
     *  a0 is typically 1 (otherwise b0 and a0 are degenerate), but this isn't enforced.  
     *
     *           
     * */ 
    class IIRFilter : public DigitalFilter 
    {
      public:
        /* Initialize filter from coeffs */ 
        IIRFilter(size_t order, const double * acoeffs, const double * bcoeffs) 
          : order(order), acoeffs(acoeffs, acoeffs+order+1), bcoeffs(bcoeffs, bcoeffs + order+1), gzero(0), gpole(0) 
        { computePolesAndZeroes(); } 

        /** initialize filter from digital zeroes/poles/gain */ 
        IIRFilter(std::complex<double> gain, size_t nzeroes, std::complex<double> * digiZeroes, size_t npoles, std::complex<double> * digiPoles) 
          : digi_gain(gain), digi_poles(digiPoles, digiPoles + npoles), digi_zeroes(digiZeroes, digiZeroes + nzeroes), gzero(0), gpole(0) 
        {
          computeCoeffsFromDigiPoles(); 
          order = npoles > nzeroes ? npoles : nzeroes; 

        }

        /** initialize filter from an FIR filter. 
         * If the FIR filter has a large delay, this might end up having quite a large order */ 
        IIRFilter(const FIRFilter& fir); 

        virtual void filterOut(size_t n, const double * w, double * out) const; 
        virtual std::complex<double> transfer(std::complex<double> z) const ;  

        /*analytic order, may not be number of coeffs if bandpass or notch. This is probably mostly useless. */
        size_t getOrder() const { return order; } 

        size_t nAcoeffs() const { return acoeffs.size(); } 
        size_t nBcoeffs() const { return bcoeffs.size(); } 
        virtual size_t nDigiPoles() const { return digi_poles.size(); }
        virtual size_t nDigiZeroes() const { return digi_zeroes.size(); }
        virtual const std::complex<double> * getDigiPoles() const { return &digi_poles[0]; } 
        virtual const std::complex<double> * getDigiZeroes()const  { return &digi_zeroes[0]; } 
        std::complex<double> getDigiGain() const { return digi_gain; } 
        virtual ~IIRFilter() {;}

        /* Make a string description of the filter. */
        TString asString() const; 
        /** Draw pole zero diagram 
         *
         *  Allowed options: 
         *    same  - draw on top 
         *    circ  - draw unit circle
         *    pol   - draw polar axis 
         * */ 
        virtual void Draw(const char * opt= "", int color=1) const; 


      protected:
        IIRFilter() { order = 0; gpole = 0; gzero = 0;} 
        size_t order; 
        std::vector<double> acoeffs; 
        std::vector<double> bcoeffs; 
        std::complex<double> digi_gain;
        std::vector<std::complex<double> > digi_poles; 
        void computePolesAndZeroes(); 
        void computeCoeffsFromDigiPoles(); 
        std::vector<std::complex<double> > digi_zeroes; 
        mutable TGraph * gzero; 
        mutable TGraph * gpole; 
    };


     
    /* Classic filter topologies */ 
    enum FilterTopology
    {
      LOWPASS, 
      HIGHPASS, 
      BANDPASS, 
      NOTCH 
    }; 


    /* Create a digital filter from bilinear transform of analog poles/zeroes/ gains */ 
    class TransformedZPKFilter : public IIRFilter
    {
      public:
        TransformedZPKFilter( int npoles, std::complex<double> *poles, int nzeroes, std::complex<double> * zeroes, double gain)
          : poles(poles, poles+npoles), zeroes(zeroes, zeroes+nzeroes), gain(gain) 
        {
          order = npoles > nzeroes ? npoles : nzeroes; 
          bilinearTransform(); 
        }

        TransformedZPKFilter(int npoles, std::complex<double> *poles, int nzeroes, std::complex<double> * zeroes, double gain, FilterTopology type, double w, double dw)
          : poles(poles, poles+npoles), zeroes(zeroes, zeroes+nzeroes), gain(gain) 
        {
          order = npoles > nzeroes ? npoles : nzeroes; 
          transform(type,w,dw); 
          bilinearTransform(); 
        }

        virtual ~TransformedZPKFilter() { ; }
        size_t nPoles() const { return poles.size(); }
        size_t nZeroes() const { return zeroes.size(); }
        const std::complex<double> * getPoles() { return &poles[0]; } 
        const std::complex<double> * getZeroes() { return &zeroes[0]; } 
        double getGain() const { return gain; } 



      protected:
        std::vector<std::complex<double> > poles; 
        std::vector<std::complex<double> > zeroes; 
        //transform from prototype filter to low/high/bandpass/notch with given frequency
        void transform(FilterTopology type, double w, double dw = 0); 
        void bilinearTransform(); 
        double gain; 
        TransformedZPKFilter() { ; } 

    }; 


    /* RC Filter  (aka any first-order filter)*/ 
    class RCFilter : public TransformedZPKFilter
    {
      public: 
        // w = 1./rc , but in normalized frequency
        //  For a Notch or Bandpass, w is the center frequency and dw is the ``width'' (i.e. the band is w-dw,w+dw)
        RCFilter(FilterTopology type = LOWPASS, double w = 0.5, double dw = 0); 
        virtual ~RCFilter() {;}

    }; 


    /*Butterworth Filter */ 
    class ButterworthFilter : public TransformedZPKFilter
    {

      public: 
        /* Create a butterworth filter of the given topology and order 
         *
         * For a Notch or Bandpass, w is the center frequency and dw is the ``width'' (i.e. the band is w-dw,w+dw)
         */
        ButterworthFilter(FilterTopology type, size_t order,  double w, double dw = 0); 
        virtual ~ButterworthFilter() {;}
    };


    /* Chebyshev Type 1 Filter */ 
    class ChebyshevIFilter : public TransformedZPKFilter
    {

      public: 
        // ripple_db is the size of the ripple 
        //  For a Notch or Bandpass, w is the center frequency and dw is the ``width'' (i.e. the band is w-dw,w+dw)
        ChebyshevIFilter(FilterTopology type, size_t order, double ripple_db, double w, double dw = 0); 
        virtual ~ChebyshevIFilter() {;}
    };


    /* Thiran AllPass filter. Can be used to implement a delay. */ 
    class ThiranFilter : public IIRFilter 
    {

      public:
        /** Construct a Thiran All pass filter with a given delay (in samples). An order can be requested or otherwise chosen automatically (to equal int(delay)) if 
         * a negative number is passed. delay must be positive. 
         */
        ThiranFilter(double delay, int order = -1) ; 


    }; 

    class MinimumPhaseFilter : public IIRFilter
    {
    
      public: 

        /* Construct a minimum phase filter from a general IIRFilter. This will reflect any zeroes and poles
         * around the unit circle. 
         *
         * If a pole or zero is exactly on the unit circle, it is moved inside by eps
         *
         * */ 
        MinimumPhaseFilter( const IIRFilter & other, double eps=1e-3); 



    }; 

    class InverseFilter : public IIRFilter
    {
      public: 

        /* Constructs the inverse of an IIRFilter (poles converted to zeroes and vice versa). It's
         * slightly more sophisticated than that because it will divide through by the gain to give a
         * 1 as the first a coeff. */
        InverseFilter(const IIRFilter &other); 
        InverseFilter(const FIRFilter &other); 
    }; 



}

#endif 
