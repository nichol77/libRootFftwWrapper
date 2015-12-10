#include "AnalyticSignal.h" 
#include "TRandom.h" 
#include "TMath.h"
#include "RFInterpolate.h" 
#include "TGraph.h" 

#ifdef ENABLE_VECTORIZE
#include "vectormath_trig.h" 
#endif


FFTtools::ThermalNoise::ThermalNoise(int nosc, double fmin, double fmax, double rms, int filter_order) 
  :nosc(nosc/(fmax-fmin)), fmin(fmin), fmax(fmax), rms(rms) , order(filter_order) 
{

  A.reserve(nosc); 
  phases.reserve(nosc); 
  freqs.reserve(nosc); 
  rng = 0; 
  isSetup = false; 
  extra_noise = 0; 
}

void FFTtools::ThermalNoise::setup() 
{
  if (!rng) rng = gRandom;
  const double norm = 1./ sqrt(nosc) * rms;

  for (int i = 0; i <nosc; i++) 
  {
    double f = 0;
    if (order < 1) 
    {
      f = rng->Uniform(fmin,fmax); 
    }
    else
    {
      double sample_max =2 * fmax; 
      while(true) 
      {
        f = rng->Uniform(0,sample_max); 
        double prob = 1.; 
        prob /= (1 + TMath::Power(fmin/f,2*order)); 
        prob /= (1 + TMath::Power(f/fmax,2*order));     
        if (prob > rng->Uniform(0,1)) 
        {
            break;
        }
      }

    }

    freqs.push_back(2*TMath::Pi()*f); 
    A.push_back( norm * sqrt(-2*log(rng->Uniform(0,1))));
    phases.push_back(rng->Uniform(0,2*TMath::Pi())); 

  }
}


void FFTtools::AnalyticSignal::eval(int N, const double* t, double * v, bool zero_v) 
{

  for (int i = 0; i < N; i++) 
  {
   if (zero_v) v[i] = eval(t[i]); 
   else v[i] += eval(t[i]); 
  }

}
double FFTtools::AnalyticSignal::eval(double t) 
{

  if (!rng) rng = gRandom; 
  double val = extra_noise ? rng->Gaus(0,extra_noise) : 0; 
  if (!isSetup) 
  {
    setup(); 
    isSetup = true; 
  }
  val += doEval(t); 

  return val ; 
}

double FFTtools::ThermalNoise::doEval(double t) 
{
  double val = 0; 

#ifdef ENABLE_VECTORIZE // honestly I think this helps more because libmath has a super slow (but maybe more accurate) sin 
  Vec4d vf; 
  Vec4d vph; 
  Vec4d vA; 
  int leftover = A.size() % 4;
  int nit = A.size()/4 + (leftover ? 1 : 0); 

  for (int i = 0; i < nit; i++)
  {
    vf.load(&freqs[4*i]); 
    vph.load(&phases[4*i]); 
    vA.load(&A[4*i]); 

    if (i == nit-1 && leftover) //hopefully this gets unrolled? 
    {
        vA.cutoff(4-leftover); 
    }

    val += horizontal_add(vA * sin(mul_add(vf,t,vph))); 
  }
 
#else
  for (unsigned i = 0; i < A.size(); i++) 
  {
    val += A[i] * sin(freqs[i] * t + phases[i]); 
  }
#endif
  return val; 
}

TGraph * FFTtools::AnalyticSignal::makeAlternatingGaussianDtWf(int N, double mean_dt, double delta_dt, double sigma, double t0) 
{
  if(delta_dt == 0) return makeGaussianDtWf(N, mean_dt, sigma, t0); 
  if (!rng) rng = gRandom; 
  TGraph * g = new TGraph(N); 
  double t = t0; 
  double dt_odd, dt_even; 
  if (rng->Rndm() > 0.5) 
  {
    dt_odd = mean_dt - delta_dt; 
    dt_even = mean_dt + delta_dt; 
  }
  else
  {
    dt_odd = mean_dt + delta_dt; 
    dt_even = mean_dt - delta_dt; 
  }

  for (int i = 0; i < N; i++) 
  {
    g->SetPoint(i, t, 0); 
    double new_dt = 0; 
    while (new_dt <=0 ) 
    {
      new_dt= rng->Gaus((i % 2 == 0) ? dt_even : dt_odd,sigma); 
    }

    t += new_dt; 
  }

  fillWf(g); 
  return g; 
}


TGraph * FFTtools::AnalyticSignal::makeGaussianDtWf(int N, double dt, double sigma, double t0) 
{
  if (!rng) rng = gRandom; 
  TGraph * g = new TGraph(N); 
  double t = t0; 
  for (int i = 0; i < N; i++) 
  {
    g->SetPoint(i, t, 0); 
    double new_dt = 0; 
    while (new_dt <=0 ) 
    {
      new_dt= rng->Gaus(dt,sigma); 
    }

    t += new_dt; 
  }

  fillWf(g); 
  return g; 
}

TGraph * FFTtools::AnalyticSignal::makeUniformJitterWf(int N, double dt, double max_jitter, double t0) 
{
  TGraph * g = new TGraph(N); 
  if (!rng) rng = gRandom; 
  for (int i = 0; i < N; i++) 
  {
    double jitter = rng->Uniform(-max_jitter, max_jitter)*dt; 
    if (i == 0) t0 -= jitter;
    g->SetPoint(i, t0 + i* dt+jitter, 0); 
  }

  fillWf(g); 
  return g; 
}

TGraph * FFTtools::AnalyticSignal::makeGaussianJitterWf(int N, double dt, double jitter_sigma, double t0) 
{
  TGraph * g = new TGraph(N); 
  if (!rng) rng = gRandom; 
  for (int i = 0; i < N; i++) 
  {
    double jitter = rng->Gaus(0, jitter_sigma)*dt; 
    if (i == 0) t0 -= jitter; 
    g->SetPoint(i, t0 + i* dt+jitter, 0); 
  }

  fillWf(g); 
  return g; 
}

TGraph * FFTtools::AnalyticSignal::makeEvenWf(int N, double dt, double t0) 
{
  TGraph * g = new TGraph(N); 
  for (int i = 0; i < N; i++) 
  {
    g->SetPoint(i, t0 + i * dt, 0); 
  }

  fillWf(g); 
  return g; 
}

void FFTtools::AnalyticSignal::fillWf(TGraph * g, bool zero_y) 
{
  eval (g->GetN(), g->GetX(), g->GetY(), zero_y); 
}

FFTtools::CompositeSignal::~CompositeSignal() 
{

  for (size_t i = 0; i < nSignals(); i++) 
  { 
    if (own[i]) 
    {
      delete getSignal(i) ; 
    }
  }

}
void FFTtools::CompositeSignal::setup() 
{
   // for (size_t i = 0; i < nSignals(); i++) { getSignal(i)->setup(); } 
}

double FFTtools::CompositeSignal::doEval(double t) 
{
  double val = 0; 
  for (size_t i = 0; i < nSignals(); i++) { val += getSignal(i)->eval(t); } 
  return val; 
}


FFTtools::BandlimitedSampledSignal::BandlimitedSampledSignal(size_t n, const double * y, double dt, double t0)
  : g(n) 
{
  for (size_t i =0; i < n; i++)
  {
    g.SetPoint(i, t0+i * dt, y[i]); 
  }
  rng = 0; 
  extra_noise = 0; 
  isSetup = true; 
}

double FFTtools::BandlimitedSampledSignal::doEval(double t) 
{
  return shannonWhitakerInterpolateValue(t,&g); 
}





