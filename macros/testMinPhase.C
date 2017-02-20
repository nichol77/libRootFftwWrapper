#include "FFTtools.h" 

void testMinPhase(int N = 256, double dt = 0.5) 
{

  int fftN = N/2+1; 


  double * G = new double[fftN]; 
  for (int i = 0; i < fftN; i++) 
  {
   G[i] = i < fftN/5  || i > fftN*4/5 ? 0 : sqrt(fftN); 
  }
  FFTWComplex * mps=FFTtools::makeMinimumPhase(fftN,G); 

  double * y = FFTtools::doInvFFT(N, mps); 


  TGraph *g = new TGraph(N); 

  for (int i = 0; i < N; i++) 
  {
    g->GetX()[i] = i * dt; 
    g->GetY()[i] = y[i]; 
  }

  TCanvas * c = new TCanvas("minphase","Minimum Phase"); 
  c->Divide(2,1); 
  c->cd(1); 

  g->Draw("alp"); 
  c->cd(2); 
  FFTtools::makePowerSpectrum(g)->Draw("alp"); 

}
