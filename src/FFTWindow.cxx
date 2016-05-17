#include "FFTWindow.h" 
#include "TGraph.h"
#include <cmath>
#include "TMath.h"


double* FFTtools::FFTWindowType::make(size_t N)  const
{
  double * ret = new double[N]; 
  fill(N,ret); 
  return ret; 
}

void FFTtools::FFTWindowType::fill(size_t N, double * x)  const
{

  size_t i; 
  int j; 
  for (i = 0, j = -(N-1)/2;  i < N; i++, j++) 
  {
    x[i] = value(j,N); 
  }
}

void FFTtools::FFTWindowType::apply(size_t N, double * x) const
{

  size_t i; 
  int j; 
  for (i = 0, j = -(N-1)/2 ; i < N; i++, j++) 
  {
    x[i] *= value(j,N); 
  }
}

double FFTtools::TriangularWindow::value(double i, size_t N) const
{
  double L = (N-1)/2;  
  return 1 - fabs(i / L); 
}

double FFTtools::HannWindow::value(double i, size_t N) const
{
  int j = i + (N-1)/2; 
  return 0.5 * (1. - cos(2.*M_PI * j / (N-1))) ; 

}

static const double hamming_alpha = 25./46; 
static const double hamming_beta = 21./46; 

double FFTtools::HammingWindow::value(double i, size_t N) const
{
  int j = i + (N-1)/2; 
  return hamming_alpha -hamming_beta* cos(2.*M_PI * j / (N-1)) ; 
}

static const double blackman_a0 = 7938./18308; 
static const double blackman_a1 = 9240./18308; 
static const double blackman_a2 = 1430./18308; 

double FFTtools::BlackmanWindow::value(double i, size_t N) const
{
  int j = i + (N-1)/2; 
   return blackman_a0 - blackman_a1 * cos(2*M_PI*j/(N-1)) + blackman_a2 * cos(4*M_PI*j/(N-1)); 
}

static const double blackmanharris_a0 = 0.35785;
static const double blackmanharris_a1 = 0.48829; 
static const double blackmanharris_a2 = 0.14128; 
static const double blackmanharris_a3 = 0.01168; 

double FFTtools::BlackmanHarrisWindow::value(double i, size_t N) const
{
  int j = i + (N-1)/2; 
   return blackmanharris_a0 - blackmanharris_a1 * cos(2*M_PI*j/(N-1)) + blackmanharris_a2 * cos(4*M_PI*j/(N-1)) - blackmanharris_a3 * cos(6*M_PI*j/(N-1)); ; 
}


double FFTtools::KaiserWindow::value(double i, size_t N) const
{
  int j = i + (N-1)/2; 
  return TMath::BesselI0(M_PI * alpha * TMath::Sqrt( 1 - TMath::Power(2.*j/(N-1)-1.,2))) / TMath::BesselI0(M_PI*alpha); 
}


double FFTtools::GaussianWindow::value(double i, size_t N) const
{
  int n = i + (N-1)/2; 
  return TMath::Exp(-2*TMath::Power(alpha*n/(N-1),2));
}



