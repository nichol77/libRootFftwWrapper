
#include "ResponseFitter.h" 

void testResponseFitter(int na=4, int nb=4, const char * rfile = "macros/resp.dat") 
{

  gRandom->SetSeed(0); 
  FFTtools::ResponseFitter f( FFTtools::ResponseFitter::Minuit2,  FFTtools::ResponseFitter::CrossCorrelation); 

  //test impulse response: 

  TGraph * resp = new TGraph(rfile); 
  resp->Draw(); 

  //put input impulse at first non-zero sample 


  TGraph imp(resp->GetN(), resp->GetX(), resp->GetY()); 
  memset(imp.GetY(),0, resp->GetN() * sizeof(double)); 

  int delay = -1; 
  for (int i = 0; i < resp->GetN(); i++) 
  { 
    if (resp->GetY()[i])
    {
      imp.GetY()[i] = 1; 
      delay = i; 
      break; 
    }
  }
  
  printf("%d\n", delay); 


  f.compute(&imp, resp, na, nb);    


  TGraph * filt = f.getBestFilter()->impulseGraph(resp->GetN(), resp->GetX()[1], delay); 
  //have to scale response by dt 
  

  filt->SetLineColor(2); 
  filt->SetMarkerColor(2); 

  filt->Draw("lsame"); 
//  filt->Draw(); 


}
