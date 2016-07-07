#include "FFTtools.h" 
#include <fftw3.h>
#include "TMath.h"
#include "TH2.h" 


#ifdef __APPLE__
#define SINCOS __sincos
#else
#define SINCOS sincos
#endif

#ifdef ENABLE_VECTORIZE
#include "vectormath_trig.h" 
#ifdef __AVX__
#define VEC Vec4d
#define VEC_T double 
#define VEC_N 4 
#define GATHER_REAL gather4d<0,2,4,6>
#define GATHER_IM gather4d<1,3,5,7>
static VEC INDEX (0,1,2,3);
#else
#define VEC Vec2d
#define VEC_T double 
#define VEC_N 2 
#define GATHER_REAL gather2d<0,2>
#define GATHER_IM gather2d<1,3>
static VEC INDEX (0,1);
#endif

#endif

/* Each point gets moved to 4 nearby evenly-spaced points */ 
static unsigned factorial_table[] = 
  { 0,
    1,
    1,
    2,
    6,
    24,
    120,
    720,
    5040,
    40320,
    362880,
    3628800,
    39916800,
    479001600
  };

/* lagrange interpolation.. I think. based on NR code*/ 
void extirpolate(double y, int n, double *yys, double x, int extirpolation_factor)
{
  // double * ys = (double*) __builtin_assume_aligned(yys, 32);
#if (__clang__ || (__GNUC__ >= 4 && __GNUC_MINOR__ >= 7))
  double * ys = (double*) __builtin_assume_aligned(yys, 32);  
#else
  double * ys = (double*) yys;
#endif  
  
  int ix = int(x); 

  extirpolation_factor = std::max(0,extirpolation_factor); 
  extirpolation_factor = std::min(sizeof(factorial_table)/sizeof(unsigned) - 1, size_t(extirpolation_factor));  

  // if it exactly equals a value, just set it 
  if (__builtin_expect(x == ix,0)) ys[ix-1] += y; 
  else
  {
    int ilo = std::min( std::max( int(x - 0.5 * extirpolation_factor),0),int(n-extirpolation_factor)); //make sure we dont' exceed bounds
    int ihi = ilo +extirpolation_factor ;

    double fac = x - ilo -1; 

    //lagrange interpolation is ugly 
    for (int j = ilo+1; j < ihi; j++) 
    {
      fac *= (x-j-1); 
    }

    int nden = factorial_table[extirpolation_factor]; 
    ys[ihi-1] += y * fac / (nden * (x- ihi)); 

    for (int j = ihi - 1; j > ilo; j--) 
    {
      nden = (nden / (j - ilo)) * (j -ihi); 
      ys[j-1] += y * fac / (nden * (x-j)); 
    }
  }
}


/* delegate so that can use optimization attributes */ 

static TGraph * _lombScarglePeriodogram(int n,  double dt, const double * __restrict  x, const double * __restrict y, double oversample_factor  = 4 , 
                       double high_factor = 2, TGraph * replaceme = 0, int extirpolation_factor = 4) 
#ifdef __clang__ 
   /* For OS X */
   // As best I can tell this would be the correct syntax for the fast math optimization with llvm
   // but it seems to not be supported :(
   // so for now we will just disable the fast math optimization for periodogram
   // [[gnu::optimize("fast-math")]];
    ; 
#else
  __attribute((__optimize__("fast-math","tree-vectorize"))); /* enable associativity and other things that help autovectorize */ 
#endif

TGraph * _lombScarglePeriodogram(int n, double dt,  const double *  __restrict x, const double * __restrict y , double oversample_factor, double high_factor,  TGraph * out, int extirpolation_factor) 
{


  int nout = n * oversample_factor * high_factor/2; 


  if (out) 
  {
    if (out->GetN() != nout) out->Set(nout); 
  }
  else
  {
    out = new TGraph(nout); 
  }

  int nfreq = n *oversample_factor * high_factor * extirpolation_factor; 
  int nwork = 2 * nfreq; 


  /* compute mean */
  double mean = 0;
  for (int i = 0; i < n; i++)  
  {
    mean += y[i]; 
  }


  mean /= n; 



  //create and zero workspace 
  double wk1[nwork] __attribute__((aligned(32)));
  double wk2[nwork] __attribute__((aligned(32)));
  memset(wk1,0, sizeof(wk1)); 
  memset(wk2,0, sizeof(wk2)); 
  
  dt = dt ? dt : (x[n-1] - x[0]) / n;
  double range = n * dt; 


  double scale_factor = 2 * nfreq / oversample_factor / range; 
  for (int i = 0; i < n; i++) 
  {
    double xx1 = fmod( (x[i] - x[0]) * scale_factor, 2 * nfreq); 
    double xx2 = fmod( 2 * (xx1++) , 2 * nfreq)+1; 

    extirpolate(y[i] - mean, nwork,wk1, xx1, extirpolation_factor); 
    extirpolate(1., nwork,wk2, xx2, extirpolation_factor); 
  }

  FFTWComplex * fft1 = (FFTWComplex* ) fftw_malloc(sizeof(FFTWComplex) * (nwork/2+1)); 
  FFTWComplex * fft2 = (FFTWComplex* ) fftw_malloc(sizeof(FFTWComplex) * (nwork/2+1)); 

  FFTtools::doFFT(nwork, wk1,fft1); 
  FFTtools::doFFT(nwork, wk2,fft2); 

  double df = 1./(range * oversample_factor); 
  double norm_factor = 4. / n; 

#ifdef ENABLE_VECTORIZE
  int leftover = nout % VEC_N; 
  int nit = nout / VEC_N + (leftover? 1: 0); 

  for (int i = 0; i < nit; i++)
  {
    VEC fft2_re=  GATHER_REAL((double*) &fft2[i*VEC_N+1]); 
    VEC fft2_im =  GATHER_IM((double*)  &fft2[i*VEC_N+1]); 
    VEC fft1_re=  GATHER_REAL((double*) &fft1[i*VEC_N+1]); 
    VEC fft1_im =  GATHER_IM((double*) &fft1[i*VEC_N+1]); 

    VEC invhypo = 1./sqrt((fft2_re* fft2_re+ fft2_im * fft2_im));

    VEC half_cos2wt = 0.5 * fft2_re  * invhypo; 
    VEC half_sin2wt = 0.5 * fft2_im  * invhypo; 

    VEC coswt = sqrt(0.5 + half_cos2wt); 
    VEC sinwt = sign_combine(sqrt(0.5 - half_cos2wt), half_sin2wt); 

    VEC density = 0.5 * n + half_cos2wt * fft2_re + half_sin2wt * fft2_im; 
    VEC costerm = square( coswt * fft1_re + sinwt * fft1_im)  / density; 
    VEC sinterm = square( coswt * fft1_re - sinwt * fft1_im)  / (n -density); 


    VEC index(i * VEC_N); 
    index += INDEX; 

    VEC ans_x = (index + 1) * df; 
    VEC ans_y = (costerm + sinterm)*norm_factor;
    if (i < nit-1 || leftover == 0)
    {
      ans_x.store(out->GetX() +i * VEC_N); 
      ans_y.store(out->GetY()+ i * VEC_N); 
    }
    else
    {
      ans_x.store_partial(leftover,out->GetX() +i * VEC_N); 
      ans_y.store_partial(leftover,out->GetY()+ i * VEC_N); 

    }
  }

#else
  for (int j = 0; j < nout; j++) 
  {
    double hypo = fft2[j+1].getAbs(); 
    double half_cos2wt = 0.5 * fft2[j+1].re / hypo; 
    double half_sin2wt = 0.5 * fft2[j+1].im / hypo; 
    double coswt = sqrt(0.5 + half_cos2wt); 
    double sinwt = TMath::Sign(sqrt(0.5 - half_cos2wt), half_sin2wt);

    double density = 0.5 * n + half_cos2wt * fft2[j+1].re + half_sin2wt * fft2[j+1].im; 

    double costerm = coswt * fft1[j+1].re + sinwt * fft1[j+1].im; 
    costerm *= costerm; 
    costerm /= density; 

    double sinterm = coswt * fft1[j+1].re - sinwt * fft1[j+1].im; 
    sinterm *= sinterm; 
    sinterm /= (n-density); 

    out->GetX()[j] = (j+1) * df; 
    out->GetY()[j] = (costerm + sinterm)  * norm_factor; 
  }
#endif
  free(fft1); 
  free(fft2); 

  return out; 

}
TGraph * FFTtools::lombScarglePeriodogram(int n, double dt, const double * __restrict x, const double * __restrict y,  double oversample_factor,
                     double high_factor, TGraph * out, int extirpolation_factor) 
{

  return _lombScarglePeriodogram(n,dt, x,y, oversample_factor, high_factor,out, extirpolation_factor); 
}



TGraph * FFTtools::lombScarglePeriodogram(const TGraph * g, double dt, double oversample_factor,
                     double high_factor, TGraph * out, int extirpolation_factor) 
{

  return _lombScarglePeriodogram(g->GetN(),dt, g->GetX(), g->GetY(), oversample_factor, high_factor,out, extirpolation_factor); 
}


TGraph * FFTtools::welchPeriodogram(const TGraph * gin, int segment_size, double overlap_fraction, const FFTWindowType * window, bool truncate_extra, TGraph * gout)
{

  double window_vals[segment_size]; 
  window->fill(segment_size, window_vals); 

  double window_weight = 0; 
  for (int i = 0; i < segment_size; i++) window_weight += window_vals[i] * window_vals[i] / segment_size; 

  double y[segment_size] __attribute__((aligned(32))); 



  TGraph * power = gout ? gout : new TGraph(segment_size/2+1); 

  if (gout)
  {
    power->Set(segment_size/2+1); 
    memset(power->GetY(),0, power->GetN() *sizeof(double)); 
  }

  double df = 1./(segment_size * (gin->GetX()[1] - gin->GetX()[0])); 
  for (int i = 0; i < power->GetN(); i++) 
  {
    power->GetX()[i] = df*i; 
  }


  int index = 0; 
  double nsegs = 0; 
  while (truncate_extra ? index + segment_size < gin->GetN() : index < gin->GetN())
  {

    for (int i = 0; i < segment_size; i++) 
    {
      y[i] = i + index < gin->GetN() ?  window_vals[i] * gin->GetY()[i + index] : 0; 
//      printf("%d %d %g %g\n", index, i, window_vals[i], y[i]); 
    }

    FFTWComplex * fft = doFFT(segment_size, y);

    power->GetY()[0] += fft[0].getAbsSq()/segment_size; 
    for (int j = 1; j < segment_size/2; j++)
    {
      power->GetY()[j] += fft[j].getAbsSq() *2 / segment_size; 
    }
    power->GetY()[segment_size/2] += fft[segment_size/2].getAbsSq() / segment_size; 

    delete [] fft; 

    index += (1.-overlap_fraction) * segment_size; 
    if (index + segment_size >= gin->GetN())
    {
      nsegs += double(gin->GetN() -index) / segment_size; 
    }
    else
    {
      nsegs +=1; 
    }
  }

  for (int i = 0; i <= segment_size/2; i++) power->GetY()[i] /= nsegs * window_weight; 

  return power; 
}




double * FFTtools::lombScarglePeriodogramSlow(int N, const double *x, const double *y, int nfreqs, const double * freqs, double * answer)
{
  if (!answer) answer = new double[nfreqs]; 

  for (int ifreq = 0; ifreq < nfreqs; ifreq++) 
  {

    double sin2wt = 0; 
    double cos2wt = 0; 
    double two_w = freqs[ifreq] * 2 * TMath::Pi(); 
    double mean = 0; 
    for (int i =0; i < N; i++)
    {
      double c,s; 
      SINCOS(x[i] * two_w, &s,&c); 
      sin2wt += s; 
      cos2wt += c; 
      mean += y[i]; 
    }

    mean /= N; 

    double tau = atan2(sin2wt, cos2wt) / two_w; 



    double ys = 0; 
    double yc = 0; 
    double s2 = 0; 
    double c2 = 0; 
    

    for (int i = 0; i < N; i++)
    {
      double yy = y[i] - mean; 

      double s,c; 

      SINCOS(two_w * (x[i] - tau), &s,&c); 

      ys += yy * s; 
      yc += yy * c; 

      s2 += s*s; 
      c2 += c*c; 
    }

    answer[ifreq] = ys*ys / s2 + yc*yc/c2; 
  }

  return answer; 
}



//TH2 *FFTtools::getPowerVsTimeUsingLombScargle(const TGraph * g, int nbins, double sampling_dt, double oversample, double 


