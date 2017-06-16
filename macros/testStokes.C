#include "FFTtools.h" 
int nfreqs = 100; 
double dt = 1; 
double min_freq = 0.2; 
double max_freq = 0.5; 


void makePols(int N, TGraph * h, TGraph * v, double amplitude_offset , double phase_offset) 
{

  h->Set(N); 
  v->Set(N); 

  for (int j = 0; j < N; j++) 
  {
    h->GetX()[j] = j *dt; 
    v->GetX()[j] = j *dt; 
    h->GetY()[j] = 0;
    v->GetY()[j] = 0; 
  }

  for (int i= 0; i < nfreqs; i++) 
  {
    double f = gRandom->Uniform(min_freq, max_freq); 
    double phase = gRandom->Uniform(0,2*TMath::Pi()); 

    for (int j =0; j < N; j++)
    {
      double t = j *dt; 
      h->GetY()[j] += sin(amplitude_offset) *sin(TMath::Pi() * f*t + phase); 
      v->GetY()[j] += cos(amplitude_offset) *sin(TMath::Pi() * f*t + phase + phase_offset); 
    }

  }



}

void testStokes(double polarization_angle =0, double phase_offset = 0) 
{

  TGraph H;  
  TGraph V;  
  makePols(1024,&H,&V,polarization_angle, phase_offset); 

  TGraph *hH=FFTtools::getHilbertTransform(&H); 
  TGraph *hV=FFTtools::getHilbertTransform(&V); 

  double i,q,u,v; 

  printf("Polarization Angle: %g degrees\nPhase offset: %g degrees\n",polarization_angle *180/TMath::Pi(),phase_offset*180/TMath::Pi()); 
  FFTtools::stokesParameters(H.GetN(), H.GetY(), hH->GetY(), V.GetY(), hV->GetY(),&i,&q,&u,&v); 

  printf("I,Q,U,V=(%g,%g,%g,%g)\n",i,q,u,v); 


}

