#include "CWT.h" 
#include "TGraph.h" 
#include "FFTWComplex.h" 
#include "TH2.h" 
#include <assert.h>
#include "FFTtools.h"


const double PI_1_4 = sqrt(sqrt(M_PI)); 


void FFTtools::CWT::MotherWavelet::fill(int M, double w, double * wavelet) const
{
    for (int k = 0; k < M; k++) 
    {
        double x = M/2 -k; 
        wavelet[k] = eval(x,w); 
    }
}

double FFTtools::CWT::Ricker::eval(double x, double w) const 
{

    double winv2 = 1./(w*w); 
    double A = 2./sqrt(3*w) * PI_1_4; 
    return  A * (1. - x*x*winv2)  * exp(-x*x *winv2/2); 

}
void FFTtools::CWT::Ricker::fill(int M, double w, double * wavelet) const
{

    double winv2 = 1./(w*w); 
    double A = 2./sqrt(3*w) * PI_1_4; 
    for (int k = 0; k < M; k++) 
    {
        double x = M/2 -k; 
        wavelet[k] = A * (1. - x*x*winv2)  * exp(-x*x *winv2/2); 
    }
}

double FFTtools::CWT::Ridger::eval(double x, double w) const
{

    double winv2 = 1./(w*w); 
    double A = 1./sqrt(w); 
    return A * -x  * exp(-x*x *winv2/2); 
}


void FFTtools::CWT::Ridger::fill(int M, double w, double * wavelet) const
{

    double winv2 = 1./(w*w); 
    double A = 1./sqrt(w); 
    for (int k = 0; k < M; k++) 
    {
       double x = M/2 -k; 
       wavelet[k] = A * -x  * exp(-x*x *winv2/2); 
    }
}

double FFTtools::CWT::Morlet::eval(double x, double w) const
{

  double winv2 = 1./(w*w); 
  double A = 1./sqrt(w) * PI_1_4; 
  return A * exp(-x*x*winv2/2) * (cos(f*x/w) - exp(winv2/2))  ; 
}


void FFTtools::CWT::Morlet::fill(int M, double w, double * wavelet) const
{

  double winv2 = 1./(w*w); 
  double A = 1./sqrt(w) * PI_1_4; 
  for (int k = 0; k < M; k++) 
  {
     double x = M/2 -k; 
     wavelet[k] = A * exp(-x*x*winv2/2) * (cos(f*x/w) - exp(winv2/2))  ; 
  }
}



void FFTtools::CWT::setScales(int n, const double * new_scales) 
{
  scales.clear(); 
  scales.insert(scales.end(), new_scales, new_scales+n); 
}

void FFTtools::CWT::setScales(int n, double start_scale, double scale_step)
{
  scales.clear(); 
  scales.reserve(n); 
  for (int i = 0; i < n; i++) 
  {
    scales.push_back(start_scale); 
    start_scale += scale_step; 
  }
}



TH2 * FFTtools::CWT::getHist(const char * name) const
{

  std::vector<double> scale_centers(scales.size()+1); 

  for (unsigned i = 1; i < scales.size(); i++) 
  {
    scale_centers[i] = (scales[i]+scales[i-1])/2; 
  }

  scale_centers[0] = scales[0] - (scale_centers[1] - scales[0]); 
  scale_centers[scales.size()] = scales[scales.size()-1] - (scale_centers[scales.size()-1] - scales[scales.size()-1]); 

  
 TH2* h = new TH2D(name,name, size, xvec, scales.size(), &scale_centers[0]); 

 for (unsigned x = 0; x < size; x++) 
 {
   for (unsigned y = 0; y < scales.size(); y++) 
   {
     h->SetBinContent(x+1,y+1, yvecs[y][x]); 
   }
 }

 return h; 
}

double FFTtools::CWT::eval(double t, double min_abs_val, int min_scale, int max_scale) const
{
  double t0 = xvec[0]; 
  double dt = xvec[1] - t0; 

  if (t < t0 || t > t0 + size*dt) return 0; 

  int i = (0.5 +(t-t0)/dt); 
  double di = t - xvec[i]; 

  double val = 0; 

  int min = min_scale > 0 ? min_scale : 0;
  int max = max_scale > 0 ? max_scale : int(nScales())-1; 

  for (int s = min; s <= max; s++)
  {
    double scale = scales[s]; 
    int M = sw*scale; 
    for (int j = -M; j <= M; j++)
    {
      if (i+j <0 ||(i +j) > int(size)) continue; 
      double coeff = yvecs[s][i+j]; 
      if (fabs(coeff) > min_abs_val)
      {
          val += coeff* mw->eval(j+di, scale) / (scale*scale); 
      }
//      printf("%d %d %f\n",i,j,val); 
    }
  }

  return val; 
}


void FFTtools::CWT::eval(int N,  const double *t, double * out, double min_abs_val, int min_scale, int max_scale) const 
{
  for (int i = 0; i < N; i++) 
  {
    out[i] = eval(t[i], min_abs_val, min_scale, max_scale); 
  }
}

TGraph * FFTtools::CWT::getGraph(int scale) const
{
  return new TGraph (size, xvec, yvecs[scale]); 
}


FFTtools::CWT::CWT(const MotherWavelet  * w , double scaled_width, int nscales, double start_scale, double scale_step)
{
  setWidth(scaled_width); 
  setScales(nscales,start_scale,scale_step); 
  setMotherWavelet(w); 

  xvec = 0; 
}

void FFTtools::CWT::cleanup()
{

  if (xvec) delete xvec; 
  for (size_t i = 0; i < yvecs.size(); i++) { delete yvecs[i]; } 
  yvecs.clear(); 

}


void FFTtools::CWT::compute(const TGraph * win)
{
  cleanup(); 
  size = win->GetN(); 

  xvec = new double[size+1]; 
  memcpy(xvec, win->GetX(), sizeof(*xvec) * size); 
  xvec[size] = xvec[0] + size * (xvec[1]-xvec[0]); //add extra for hist 
  yvecs.reserve(scales.size()); 

  for (size_t j = 0; j < scales.size(); j++) 
  {

    double w = scales[j]; 
    int M = 2*w*sw+1; 
    double wavelet[M]; 
    mw->fill(M,w,wavelet); 
    yvecs.push_back(FFTtools::directConvolve(size, win->GetY(), M, wavelet)); 
    for (size_t i = 0; i < size; i++) 
    {
      yvecs[j][i] /= sqrt(w); 
    }
  }
}




