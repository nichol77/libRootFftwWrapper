#include "FFTtools.h" 
int nfreqs = 100; 
double noise = 1; 
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
    h->GetY()[j] = gRandom->Gaus(0,noise);
    v->GetY()[j] = gRandom->Gaus(0,noise); 
  }

  for (int i= 0; i < nfreqs; i++) 
  {
    double f = gRandom->Uniform(min_freq, max_freq); 
    double phase = gRandom->Uniform(0,2*TMath::Pi()); 

    for (int j =0; j < N; j++)
    {
      double t = j *dt; 
      h->GetY()[j] += cos(amplitude_offset) *sin(TMath::Pi() * f*t + phase); 
      v->GetY()[j] += sin(amplitude_offset) *sin(TMath::Pi() * f*t + phase + phase_offset); 
    }

  }



}

void testStokes(double polarization_angle =0, double phase_offset = 0) 
{

  gRandom->SetSeed(0); 

  TGraph * H = new TGraph(1024);  
  TGraph * V = new TGraph(1024);  
  makePols(1024,H,V,polarization_angle, phase_offset); 

  TGraph *hH=FFTtools::getHilbertTransform(H); 
  TGraph *hV=FFTtools::getHilbertTransform(V); 

  TGraph *gI = new TGraph(H->GetN(),H->GetX(), H->GetY()); 
  TGraph *gQ = new TGraph(H->GetN(),H->GetX(), H->GetY()); 
  TGraph *gU = new TGraph(H->GetN(),H->GetX(), H->GetY()); 
  TGraph *gV = new TGraph(H->GetN(),H->GetX(), H->GetY()); 


  printf("Polarization Angle: %g degrees\nPhase offset: %g degrees\n",polarization_angle *180/TMath::Pi(),phase_offset*180/TMath::Pi()); 

  FFTtools::stokesParameters(H->GetN(), H->GetY(), hH->GetY(), V->GetY(), hV->GetY(),gI->GetY(),gQ->GetY(),gU->GetY(),gV->GetY(),true); 
  double i,q,u,v; 

  i = gI->GetY()[1023]; 
  q = gQ->GetY()[1023]; 
  u = gU->GetY()[1023]; 
  v = gV->GetY()[1023]; 

  printf("I,Q,U,V=(%g,%g,%g,%g)\n",i,q,u,v); 
  printf("measured polarization angle: %g\n", atan2(u,q)/2); 


}

