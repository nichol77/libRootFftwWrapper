#include "AnalyticSignal.h" 
#include "FFTtools.h" 
void testEnvelope()
{
  gRandom->SetSeed(0); 
 
  FFTtools::ThermalNoise noise(400,0.18/1.3,1.2/1.3,10); 

  TGraph * g = noise.makeGaussianDtWf(260, 1/2.6, 0.09,0); 

  double fcarrier = gRandom->Uniform(1.0,1.2); 
  double pcarrier = gRandom->Uniform(0, 2 * TMath::Pi()); 
  double f_mod = gRandom->Uniform(0.01,0.02); 
  double p_mod = gRandom->Uniform(0, 2 * TMath::Pi()); 
  double A = 50; 


  for (int i = 0; i < g->GetN(); i++)
  {
    printf("%g %g",g->GetX()[i],g->GetY()[i]); 
    g->GetY()[i] += A * sin(fcarrier *g->GetX()[i] + pcarrier) * sin(f_mod * g->GetX()[i] + p_mod); 

    printf("%g\n",g->GetY()[i]); 
  }


  TGraph * ig = FFTtools::getInterpolatedGraph(g,1/2.6); 

  TGraph * g_hilbert = FFTtools::getHilbertEnvelope(ig); 
  TGraph * g_rms = new TGraph(g->GetN(),g->GetX(),g->GetY()); 
  TGraph * g_peak = new TGraph(g->GetN(),g->GetX(),g->GetY()); 

  FFTtools::rmsEnvelope(g->GetN(), 10, g->GetX(), g->GetY(), g_rms->GetY()); 
  FFTtools::peakEnvelope(g->GetN(), 5, g->GetX(), g->GetY(), g_peak->GetY()); 

  g->Draw(); 

  g_hilbert->SetLineColor(2); 
  g_hilbert->Draw("lsame"); 

  g_rms->SetLineColor(3); 
  g_rms->Draw("lsame"); 
  g_peak->SetLineColor(4); 
  g_peak->Draw("lsame"); 


  





}
