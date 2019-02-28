#include "FFTtools.h" 
TMultiGraph * phasemod(TGraph *g, double mindb = -30, double delay = 20) 
{

  TGraph * m = FFTtools::makeMinimumPhase(g,mindb); 
  TGraph * f = FFTtools::makeFixedDelay(g,delay); 
  m->SetLineColor(2); 
  f->SetLineColor(4); 

  TGraph *ge =FFTtools::getHilbertEnvelope(g);
  TGraph *me =FFTtools::getHilbertEnvelope(m);
  TGraph *fe =FFTtools::getHilbertEnvelope(f);

  ge->SetLineStyle(2); 
  me->SetLineStyle(2); 
  fe->SetLineStyle(2); 

  me->SetLineColor(2); 
  fe->SetLineColor(4); 

  ge->SetTitle(TString("env ") + TString(g->GetTitle()));; 
  me->SetTitle(TString("env ") + TString(m->GetTitle()));; 
  fe->SetTitle(TString("env ") + TString(f->GetTitle()));; 


  TMultiGraph * mg = new TMultiGraph; 
  mg->Add(g);
  mg->Add(m);
  mg->Add(f);
  mg->Add(ge);
  mg->Add(me);
  mg->Add(fe);



  mg->Draw("al"); 
  gPad->BuildLegend(); 

  return mg; 
}
