#include "FFTtools.h"
#include <iostream>


#include <TMath.h>
#include <TApplication.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TLegend.h>

int main(int argc, char* argv[]){


  TApplication* theApp = new TApplication("App", &argc, argv);

  const int n = 256;
  double omega1 = 0.8;
  double omega2 = 0.6;


  double ts[n];
  double dt = 0.1;
  double phi = TMath::Pi()/2;
  double sineWave1[n];
  double sineWave2[n];


  for(int i=0; i<n; i++){
    double t = dt*i;
    ts[i] = t;
    sineWave1[i] = 3.3*TMath::Sin(TMath::TwoPi()*omega1*t + phi);
    sineWave2[i] = 1.6*TMath::Sin(TMath::TwoPi()*omega2*t);
  }

  TGraph* gr1 = new TGraph(n, ts, sineWave1);
  gr1->SetName("grSineWave");
  gr1->SetTitle("Sine waves; t (s); Amplitude (V)");
  TGraph* gr2 = new TGraph(n, ts, sineWave2);
  gr2->SetLineColor(kRed);

  TLegend* l1 = new TLegend(0.8, 0.8, 1, 1);
  l1->AddEntry(gr1, "y = 3.3 sin(2#pi 0.8t)", "l");
  l1->AddEntry(gr2, "y = 1.6 sin(2#pi 0.3t)", "l");

  TCanvas* c1 = new TCanvas("c1", "FFT Frivolities", 1200, 1600);
  c1->Divide(2, 2);
  c1->cd(1);
  gr1->Draw();
  gr2->Draw("l same");
  l1->Draw();


  TGraph* grCorr = FFTtools::getCorrelationGraph(gr1, gr2);
  grCorr->SetTitle("FFTtools::getCorrelationGraph;offset #deltat (s);Correlation");
  
  c1->cd(2);  
  grCorr->Draw();


  TGraph* grPs1 = FFTtools::makePowerSpectrumVoltsSecondsdB(gr1);
  for(int i=0; i<grPs1->GetN(); i++){
    grPs1->GetX()[i] *= 1e6; // Undo conversion of freqs to MHz
  }
  TGraph* grPs2 = FFTtools::makePowerSpectrumVoltsSecondsdB(gr2);
  for(int i=0; i<grPs2->GetN(); i++){
    grPs2->GetX()[i] *= 1e6; // Undo conversion of freqs to MHz
  }


  c1->cd(3);
  grPs1->Draw();
  grPs1->SetTitle("FFTtools::makePowerSpectrumVoltsSecondsdB; Frequency (Hz); Power (dB)");
  grPs2->SetLineColor(kRed);
  grPs2->Draw("lsame");
  l1->Draw();


  FFTWComplex* sineFFT1 = FFTtools::doFFT(n, sineWave1);
  FFTWComplex* sineFFT2 = FFTtools::doFFT(n, sineWave2);

  double* sineWave1Again = FFTtools::doInvFFT(n, sineFFT1);
  double* sineWave2Again = FFTtools::doInvFFT(n, sineFFT2);  
  
  TGraph* gr1Again = new TGraph(n, ts, sineWave1Again);
  TGraph* gr2Again = new TGraph(n, ts, sineWave2Again);  

  c1->cd(4);
  gr1Again->SetTitle("Doing a forward and then inverse FFT on the sine waves; t (s); Amplitude (V)");
  gr1Again->Draw();
  gr2Again->SetLineColor(kRed);
  gr2Again->Draw("lsame");
  
  
  /* Draw pretty things */
  
  c1->Update();
  
  c1->SaveAs(TString::Format("%sOutput.png", argv[0]));

  gSystem->ProcessEvents();
  std::cerr << "Select File->Quit to quit." << std::endl;
  theApp->Run();

  return 0;



}
