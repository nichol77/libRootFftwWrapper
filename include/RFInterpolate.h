#ifndef __RF_INTERPOLATE_H 
#define __RF_INTERPOLATE_H 


/* Interpolation routines appropriate for band-limited waveforms 
 *
 * Cosmin Deaconu <cozzyd@kicp.uchicago.edu> 
 *
 */ 
class TGraph; 
class TGraphErrors; 
class TH2; 
class FFTWComplex; 

#include "FFTWindow.h" 

namespace FFTtools
{
    /* Transforms a graph with uneven spacing to one with even spacing
     * dt should be the effective sample rate (or 0 to estimate from the values). 
     *
     * Waveform will have g->GetN() * supersample values, i.e. supersample by passing supersample > 1, NOT by passing a smaller dt
     *
     * Inverts Shannon-Whitaker interpolation formula using matrix decomposition (and therefore takes O(n^3) time... sorry!)
     *
     */ 
    TGraph * getInterpolatedGraphInvert(const TGraph * g, double dt = 0, int nout = 0); 
    TGraph * getInterpolatedGraphWeightedInvert(const TGraph * g, double dt = 0, int nout = 0); 

    /* RECOMMENDED  
     * Sets all interpolator values farther than max_dist apart to 0 and use sparse operations. 
     * All values smaller than eps are set to 0; 
     * Values are weighted based on distance from even grid with exponent weight_exp (so weight_exp = 0 is no weighting) 
     * mu is regularization factor
     * lanczos_window windows weights by lanczos window (max width is max_dist) 
     * A will the sparse matrix as a histogram, if passed
     *
     * Errors on y will be set to error_scale/sqrt(sum(sincfactor^2)); 
     *
     *
     * */ 
    TGraphErrors * getInterpolatedGraphSparseInvert(const TGraph * g, double dt = 0, int nout = 0, 
                                                   double max_dist = 32, double eps = 0, double weight_exp = 0, 
                                                   double mu = 1e-3, int regularization_order = 0, double error_scale = 1, 
                                                  TH2 * A =0); 

    void getInterpolatedGraphSparseInvert(const TGraph * in, TGraph * out, 
                                                   double max_dist = 32, double eps = 0, double weight_exp = 0, 
                                                   double mu = 1e-3, int regularization_order = 0, double error_scale = 1, 
                                                   TH2 * A =0); 


    /* Same as getInterpolatedGraphInvert but split into overlapping sub signals */ 
    TGraph * getInterpolatedGraphInvertLapped(const TGraph * g, double dt = 0, int lapsize = 64, int nout = 0); 

    /* Same as getInterpolatedGraphInvert but split into sub signals with overlap of overlap */ 
    TGraph * getInterpolatedGraphInvertSplit(const TGraph * g, double dt = 0, int splitsize = 64, int overlap = 8, int out = 0); 

    /* Invert estimate of uneven DFT */ 
    TGraph * getInterpolatedGraphDFT(const TGraph *g, double dt = 0, int nout = 0, double maxF = 0); 
    /* Estimate DFT from uneven spacing... */
    FFTWComplex * getUnevenDFT(const TGraph *g, double df, int npts); 

    //supersample "exactly" using shannon whitaker interpolation (with FIR filter); 
    TGraph * supersample(const TGraph *g, int supersample_factor, int sinc_radius = 16); 

    TGraph * getInterpolatedGraphSincKernel(const TGraph * g, double dt = 0, int nout = 0); 
    


    /* Transforms a graph with uneven spacing to one with even spacing 
     * using lagrange interpolation
     *
     * dt should be the effective sample rate (or 0 to estimate from the values). 
     *
     * Waveform will have g->GetN() * supersample values, i.e. supersample by passing supersample > 1, NOT by passing a smaller dt
     *
     * THIS IS EXTREMELY INEFFICIENT AND GIVES SAME ANSWERS AS getInterpolatedGraphInvert 
     *
     * */ 
    TGraph * getInterpolatedGraphLagrange(const TGraph * g, double dt = 0, double supersample = 1); 



    /** Shannon-Whitaker interpolation of a regular graph */ 
    double shannonWhitakerInterpolateValue(double t, const TGraph * regular_graph, int max_lobe = 0, const FFTWindowType * win = 0); 
    double shannonWhitakerInterpolateValueAndError(double t, const TGraphErrors * regular_graph, double * err, int max_lobe = 0, const FFTWindowType* win = 0); 
    double linearInterpolateValueAndError(double t, const TGraph* regular_graph, double * err = 0); 

    /** Derivative of Shannon-Whitaker interpolation of a regular graph */ 
    double shannonWhitakerInterpolateDerivative(double t, const TGraph * regular_graph); 


}




#endif
