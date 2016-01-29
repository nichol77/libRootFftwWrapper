#include "FFTtools.h"
#include <iostream>


#include <TMath.h>
#include <TApplication.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TRandom3.h>

TCanvas* testFftFunctions();
TCanvas* testCorrelateInterpolateAndAverage();

int main(int argc, char* argv[]){

  TApplication* theApp = new TApplication("App", &argc, argv);

  // TCanvas* c1 = testFftFunctions();
  // TCanvas* c2 = testCorrelateInterpolateAndAverage();
  testFftFunctions();
  testCorrelateInterpolateAndAverage();
  
  gSystem->ProcessEvents();
  std::cerr << "Select File->Quit to quit." << std::endl;
  theApp->Run();

  return 0;  
}



TCanvas* testCorrelateInterpolateAndAverage(){

  const Int_t nSamp = 256;
  const Double_t deltaT = 1./2.6;

  const UInt_t mySeed = 1984;
  TRandom3 rand(mySeed);
  const Double_t gausMean = 0;
  const Double_t gausSigma = 1;  
  const Double_t fakePulseHeight = 10;
  const Double_t fakePulseWidth = 4;

  const Int_t numGraphs = 5;

  TCanvas* c1 = new TCanvas("c1", "Interpolate, Correlate and Average", 1200, 600); //1600);
  c1->Divide(2);
  c1->cd(1);
  TGraph* grs[numGraphs] = {NULL};
  
  for(int grInd=0; grInd < numGraphs; grInd++){
    const Double_t fakePulseTime = 30 + 10*grInd;

    Double_t volts[nSamp];
    Double_t times[nSamp];  
    for(int samp=0; samp<nSamp; samp++){
      times[samp] = deltaT*samp;
      volts[samp] = rand.Gaus(gausMean, gausSigma);

      Double_t expArg = -(times[samp]-fakePulseTime)*(times[samp]-fakePulseTime)/(fakePulseWidth*fakePulseWidth);
      volts[samp] += fakePulseHeight*exp(expArg);

    }
  
    grs[grInd] = new TGraph(nSamp, times, volts);
    grs[grInd]->SetTitle(TString::Format("Pulse time %2.0lf", fakePulseTime));    
    
    grs[grInd]->SetLineColor(grInd+1);
    grs[grInd]->SetFillColor(0); // For nice make legend function 

    TString drawOpt = grInd == 0 ? "alp" : "lpsame";
    grs[grInd]->Draw(drawOpt);

  }

  TLegend* l = gPad->BuildLegend();
  l->Draw();

 
  TGraph* grInterpCorrAve = FFTtools::interpolateCorrelateAndAverage(deltaT, numGraphs, grs);

  c1->cd(2);
  grInterpCorrAve->Draw("alp");
  grInterpCorrAve->SetTitle("from interpolateCorrelateAndAverage");
  grInterpCorrAve->SetLineStyle(2);
  grInterpCorrAve->SetFillColor(0);

  grs[0]->Draw("lpsame");
  
  TGraph* grDiff = (TGraph*) grInterpCorrAve->Clone("grDiff");
  grDiff->SetFillColor(0);
  grDiff->SetLineColor(kMagenta);  
  for(int samp=0; samp < grDiff->GetN(); samp++){
    grDiff->GetY()[samp] -= grs[0]->GetY()[samp];
  }
  grDiff->Draw("lpsame");
  grDiff->SetTitle("Difference");
  TLegend* l2 = gPad->BuildLegend();
  l2->Draw();


  grs[0]->SetTitle("Noisy gaussian pulses");  
  grInterpCorrAve->SetTitle("Interpolate, correlate and average");

  c1->Update();

  return c1;
}




TCanvas* testFftFunctions(){

  

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

  TCanvas* c2 = new TCanvas("c2", "FFT Frivolities", 1200, 1600);
  c2->Divide(2, 2);
  c2->cd(1);
  gr1->Draw("alp");
  gr2->Draw("l same");
  l1->Draw();


  TGraph* grCorr = FFTtools::getCorrelationGraph(gr1, gr2);
  grCorr->SetTitle("FFTtools::getCorrelationGraph;offset #deltat (s);Correlation");
  
  c2->cd(2);  
  grCorr->Draw("alp");


  TGraph* grPs1 = FFTtools::makePowerSpectrumVoltsSecondsdB(gr1);
  for(int i=0; i<grPs1->GetN(); i++){
    grPs1->GetX()[i] *= 1e6; // Undo conversion of freqs to MHz
  }
  TGraph* grPs2 = FFTtools::makePowerSpectrumVoltsSecondsdB(gr2);
  for(int i=0; i<grPs2->GetN(); i++){
    grPs2->GetX()[i] *= 1e6; // Undo conversion of freqs to MHz
  }


  c2->cd(3);
  grPs1->Draw("alp");
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

  c2->cd(4);
  gr1Again->SetTitle("Doing a forward and then inverse FFT on the sine waves; t (s); Amplitude (V)");
  gr1Again->Draw("alp");
  gr2Again->SetLineColor(kRed);
  gr2Again->Draw("lsame");
  
  
  /* Draw pretty things */
  
  c2->Update();
  
  return c2;



}
