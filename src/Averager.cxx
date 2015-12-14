#include "Averager.h" 
#include "TGraph.h" 
#include "TTree.h"
#include "TEventList.h"
#include "FFTtools.h"


FFTtools::Averager::Averager() 
{
  n = 0; 
  avg = 0; 
  dirty = true; 
  t0 = 0; 
  t1 = 0; 
  max_shift = 0; 
  //h = new TH1F("h", "Delta t histo", 101, -50, 50);
}


void FFTtools::Averager::reset() 
{
  if (avg)
  {
    delete avg; 
    delete sum; 
    delete norm; 
  }
  n = 0; 
  dirty = true; 
}


void FFTtools::Averager::add(TTree*t, const char * branch_name, TEventList * list, double interpolation_dt) 
{

  TGraph * g = 0; 

  TBranch * b = t->GetBranch(branch_name); 
  char * oldaddr = b->GetAddress(); 
  b->SetAddress(&g); 

  int n = list ? list->GetN() : t->GetEntries(); 

  for (int i = 0; i <n; i++) 
  {
    //printf("%d/%d\n",i,n), 
    b->GetEntry( list ? list->GetEntry(i) : i); 
    //printf("%lld\n",list->GetEntry(i)), 
    add(g,interpolation_dt); 
  }


  delete g; 
  if (oldaddr) b->SetAddress(oldaddr); 

}

void FFTtools::Averager::add(const TGraph *g, double interpolate_dt) 
{

  if (!interpolate_dt)
  {
    add(g); 
    return; 
  }


  TGraph * ig = FFTtools::getInterpolatedGraph((TGraph*)g,interpolate_dt); 
  //zero mean 
  

  double mean = ig->GetMean(2); 
  if (fabs(mean) > 1e-6)
  {
    for (int i = 0; i < ig->GetN(); i++)
    {
      ig->GetY()[i]-=mean; 
    }
  }


  add(ig); 
  delete ig; 


}
void FFTtools::Averager::add(const TGraph *g) 
{
  dirty = true;

  TGraph * gg = (t0 || t1) ? FFTtools::cropWave((TGraph*)g, t0,t1) : (TGraph*)g; 


  if (!n) 
  {
    avg = new TGraph(gg->GetN(), gg->GetX(), gg->GetY()); 
    sum = new TGraph(avg->GetN(), avg->GetX(), avg->GetY()); 
    norm = new int[gg->GetN()]; 
    for (int i = 0; i < gg->GetN(); i++) { norm[i] = 1; } 
    n++;
    return; 
  }

  TGraph * cg = FFTtools::getCorrelationGraph(avg,gg); 

  int peak = -1; 
  double max = 0; 
  int shift = 0; 

  for (int i = 0; i <cg->GetN(); i++)
  {
    double shiftval = cg->GetX()[i]; 
    if (fabs(shiftval) > max_shift) continue; 
    if (peak < 0 || cg->GetY()[i] > max)
    {
      shift = peak - cg->GetN()/2; 
      max = cg->GetY()[i]; 
      peak = i; 
    }
  }



  delete cg; 
  
  for (int j = 0; j < sum->GetN(); j++)
  {
     int jj = j - shift; 
     if (jj < 0 || jj >= gg->GetN()) continue; 
     norm[j]++; 
     sum->GetY()[j] += gg->GetY()[jj]; 
  }

  n++; 
  
  if (t0 || t1)
  {
    delete gg; 
  }
}

void FFTtools::Averager::computeAverage()
{
  if (!dirty) return; 

  for (int j = 0; j < avg->GetN(); j++)
 {
    avg->GetY()[j] = norm[j] > 0 ? sum->GetY()[j] / norm[j] : 0; 
  }

  dirty = false; 
}
