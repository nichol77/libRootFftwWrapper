#include "FFTtools.h" 


#include "TMath.h"

/* Each point gets moved to 4 nearby evenly-spaced points */ 
const int extirpolation_factor = 4; 

static int factorial_table[] = 
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
    479001600,
  };

/* lagrange interpolation.. I think. based on NR code*/ 
static void extirpolate(double y, int n, double *ys, double x)
{

  int ix = int(x); 

  // if it exactly equals a value, just set it 
  if (x == ix) ys[ix-1] += y; 
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

TGraph * FFTtools::lombScarglePeriodogram(const TGraph * g, double oversample_factor,
                     double high_factor, TGraph * out) 
{


  int n = g->GetN(); 
  int nout = n * oversample_factor * high_factor/2; 

  const double * x = g->GetX(); 
  const double * y = g->GetY(); 


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


  /* compute mean and variance... fast! */ 
  double mean = 0, var = 0; 

  for (int i = 0; i < n; i++)  
  {
    mean += y[i]; 
    var += y[i]*y[i]; 
  }


  mean /= n; 
  var = var/n - mean * mean; 



  //create and zero workspace 
  double wk1[nwork];
  double wk2[nwork];
  memset(wk1,0, sizeof(wk1)); 
  memset(wk2,0, sizeof(wk2)); 
  
  double range = x[n-1] - x[0]; 


  double scale_factor = 2 * nfreq / oversample_factor / range; 
  for (int i = 0; i < n; i++) 
  {
    double xx1 = fmod( (x[i] - x[0]) * scale_factor, 2 * nfreq); 
    double xx2 = fmod( 2 * (xx1++) , 2 * nfreq)+1; 

    extirpolate(y[i] - mean, nwork,wk1, xx1); 
    extirpolate(1., nwork,wk2, xx2); 

  }

  FFTWComplex * fft1 = FFTtools::doFFT(nwork, wk1); 
  FFTWComplex * fft2 = FFTtools::doFFT(nwork, wk2); 

  double df = 1./(range * oversample_factor); 

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
    out->GetY()[j] = (costerm + sinterm) / (2 * var); 
  }

  return out; 

}




TGraph * FFTtools::welchPeriodogram(const TGraph * gin, int segment_size, double overlap_fraction, const FFTWindowType * window, bool truncate_extra)
{

  double window_vals[segment_size]; 
  window->fill(segment_size, window_vals); 

  double window_weight = 0; 
  for (int i = 0; i < segment_size; i++) window_weight += window_vals[i] * window_vals[i] / segment_size; 

  double y[segment_size]; 

  TGraph * power = new TGraph(segment_size/2+1); 

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



