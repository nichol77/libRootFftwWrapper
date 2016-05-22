#include "SineSubtract.h" 
#include "TGraph.h" 
#include "TCanvas.h" 
#include "TStyle.h" 
#include <malloc.h>

#ifdef SINE_SUBTRACT_PROFILE
#include "TStopwatch.h"
#endif

#include "TMath.h"
#include <set>
#include "FFTtools.h"
#include "TF1.h" 
#include "TH2.h"

#ifdef SINE_SUBTRACT_USE_FLOATS
#define VEC_T float
#else
#define VEC_T double
#endif

#ifdef ENABLE_VECTORIZE
#include "vectormath_trig.h" 

//assume AVX is present, otherwise it'll be simulated 

#ifdef SINE_SUBTRACT_USE_FLOATS
#define VEC Vec8f 
#define VEC_N 8
#else
#define VEC Vec4d 
#define VEC_N 4 
#endif


#endif 

extern int gErrorIgnoreLevel; // who ordered that? 



static double normalize_angle(double phi)
{
  return  phi - 2*TMath::Pi() * FFTtools::fast_floor((phi+ TMath::Pi()) / (2.*TMath::Pi())); 
}

FFTtools::SineFitter::SineFitter()
{
  min.SetFunction(f); 
  verbose = false; 

}

FFTtools::SineFitter::~SineFitter()
{


}

void FFTtools::SineFitter::setGuess(double fg, int ntrace, const double * phg, double ampg)
{
  freq = fg; 
  phase.clear(); 
  amp.clear(); 
  phase_err.clear(); 
  amp_err.clear(); 

  phase.insert(phase.end(), phg, phg + ntrace); 
  amp.insert(amp.end(), ntrace, ampg); 
  phase_err.insert(phase_err.end(), ntrace, 0); 
  amp_err.insert(amp_err.end(), ntrace, 0); 
}

ROOT::Math::IBaseFunctionMultiDim* FFTtools::SineFitter::SineFitFn::Clone() const
{
  SineFitFn * fn = new SineFitFn; 
#if defined(SINE_SUBTRACT_USE_FLOATS) || defined(SINE_SUBTRACT_FORCE_ALIGNED)
  fn->setXY(nt,ns,xp,yp); //incredibly inefficient, but quick kludge for now. We're probalby never going to clone anyway. 
#else
  fn->setXY(nt,ns,x,y); 
#endif
  return fn; 
}

FFTtools::SineFitter::SineFitFn::SineFitFn() 
{
  ns = 0; 
  nt = 0; 
  x= 0; 
  y = 0; 
#if defined(SINE_SUBTRACT_USE_FLOATS) || defined(SINE_SUBTRACT_FORCE_ALIGNED)
  xp = 0; 
  yp = 0; 
#endif

}

double FFTtools::SineFitter::SineFitFn::DoDerivative(const double * p, unsigned int coord) const
{

  double w = 2 *TMath::Pi() * p[0]; 
  double deriv = 0; 
  int type = coord == 0 ? 0 : (1 + ((coord-1) % 2)); 

  for (int ti = 0; ti < nt; ti++) 
  {

    if (type > 0 && (int(coord) - 1)/2 != ti) continue; 
    double ph = normalize_angle(p[1 + 2*ti]); 
    double A = p[2+2*ti]; 
    if (wgt) A*=wgt[ti]; 

#ifdef ENABLE_VECTORIZE 

    VEC vecx(0); 
    VEC vecy(0); 
    VEC vec_ph(ph); 
    VEC vec_w(w); 
    VEC vec_A(A); 
    VEC vec_cos(0); 
    VEC dYdp(0); 
    int leftover = ns % VEC_N;
    int nit = ns/VEC_N + (leftover ? 1 : 0); 

#ifdef SINE_SUBTRACT_FORCE_ALIGNED
    double * xx = (double*) __builtin_assume_aligned(x[ti], VEC_N * sizeof(VEC_T)); 
    double * yy = (double*) __builtin_assume_aligned(y[ti], VEC_N * sizeof(VEC_T)); 
#elif defined(SINE_SUBTRACT_USE_FLOATS)
    float * xx = (float*) __builtin_assume_aligned(x[ti], VEC_N * sizeof(VEC_T)); 
    float * yy = (float*) __builtin_assume_aligned(y[ti], VEC_N * sizeof(VEC_T)); 
#else
    double * xx = (double*) x[ti]; 
    double * yy = (double*) y[ti]; 
#endif
    

    for (int i = 0; i < nit; i++)
    {
#if defined(SINE_SUBTRACT_FORCE_ALIGNED) || defined(SINE_SUBTRACT_USE_FLOATS)
      vecx.load_a(xx+VEC_N*i); 
      vecy.load(yy+VEC_N*i); 
#else
      vecx.load(xx+VEC_N*i); 
      vecy.load(yy+VEC_N*i); 
#endif

      VEC vec_ang = mul_add(vecx, vec_w, vec_ph); 
      VEC vec_sin = sincos(&vec_cos, vec_ang); 
      VEC vec_Y = mul_sub(vec_sin, vec_A, vecy); 

      switch (type) 
      {
        case 0: 
          dYdp = vec_A *vecx*vec_cos * (2 * M_PI); 
          break; 
        case 1: 
          dYdp = vec_A *vec_cos; 
          break; 
        case 2: 
          dYdp = vec_sin; 
          break; 
      }
#ifndef SINE_SUBTRACT_DONT_HORIZONTAL_ADD
      if (i == nit-1 && leftover) //hopefully this gets unrolled? 
      {
        vec_Y.cutoff(leftover); 
      }

      deriv +=horizontal_add( vec_Y * dYdp)/ns; 
#else

      VEC ans = vec_Y * dYdp; 
      VEC_T ans_v[VEC_N] __attribute__((aligned(sizeof(VEC))); 
      ans.store_a(ans_v); 

      int vecn = (i == nit-1 && leftover) ? leftover : VEC_N;
      for (int j = 0; j < vecn; j++) 
      {
        deriv += ans_v[j]/ns; 
      }

#endif
    }

#else
    for (int i = 0; i < ns; i++)
    {
      double t = x[ti][i]; 
      double sinang = sin(w*t + ph); 
      double Y = A *sinang; 
      double dYdp = 0; 

      switch(type) 
      {
        case 0: 
          dYdp = A*t*cos(w*t + ph) * 2 * TMath::Pi(); 
          break; 
        case 1: 
          dYdp = A*cos(w*t + ph); 
          break; 
        case 2:
          dYdp = sinang; 
          break; 
      }

      deriv += ((Y-y[ti][i]) * dYdp)/ns; 
    }
#endif
  }

  deriv *=2; 

//  printf("dP/d%s (f=%f, ph=%f, A=%f) = %f\n", coord == 0 ? "f" : coord == 1? "ph" : "A", p[0], ph, A, deriv);  
  return deriv/nt; 
}



double FFTtools::SineFitter::SineFitFn::DoEval(const double * p) const
{

  VEC_T w = 2 *TMath::Pi() * p[0]; 
  VEC_T power = 0; 


  for (int ti = 0; ti < nt; ti++)
  {
    VEC_T ph = normalize_angle(p[1+2*ti]); 
    VEC_T A = p[2+2*ti]; 
    if (wgt) A*= wgt[ti]; 

#ifdef ENABLE_VECTORIZE 
    VEC vecx(0); 
    VEC vecy(0); 
    VEC vec_ph(ph); 
    VEC vec_w(w); 
    VEC vec_A(A); 

#ifdef SINE_SUBTRACT_FORCE_ALIGNED
    double * xx = (double*) __builtin_assume_aligned(x[ti], VEC_N * sizeof(VEC_T)); 
    double * yy = (double*) __builtin_assume_aligned(y[ti], VEC_N * sizeof(VEC_T)); 
#elif defined(SINE_SUBTRACT_USE_FLOATS)
    float * xx = (float*) __builtin_assume_aligned(x[ti], VEC_N * sizeof(VEC_T)); 
    float * yy = (float*) __builtin_assume_aligned(y[ti], VEC_N * sizeof(VEC_T)); 
#else
    double * xx = (double*) x[ti]; 
    double * yy = (double*) y[ti]; 
#endif


    int leftover = ns % VEC_N;
    int nit = ns/VEC_N + (leftover ? 1 : 0); 

    for (int i = 0; i < nit; i++)
    {
      vecx.load(xx+VEC_N*i); 
      vecy.load(yy+VEC_N*i); 
      VEC vec_ang = mul_add(vecx, vec_w, vec_ph); 
      VEC vec_sin = sin(vec_ang); 
      VEC vec_Y = mul_sub(vec_sin, vec_A, vecy); 

      VEC vec_Y2 = square(vec_Y); 
#ifndef SINE_SUBTRACT_DONT_HORIZONTAL_ADD
      if (i == nit-1 && leftover) //hopefully this gets unrolled? 
      {
        vec_Y2.cutoff(leftover); 
      }

      power += horizontal_add(vec_Y2)/ns; 
#else 
      int vecn = (i == nit-1 && leftover) ? leftover : VEC_N;
      for (int j = 0; j< vecn; j++)
      {
        power += vec_Y2[j]/ns; 
      }
#endif
    }
#else
    for (int i = 0; i < ns; i++)
    {
      VEC_T Y = A * sin(w*x[ti][i] + ph) -y[ti][i]; 
      power += Y*Y/ns; 
    }
#endif 

  }
//  printf("P(f=%f, ph=%f, A=%f) = %f\n", p[0], ph, A, power);  
  return power/nt;

}

FFTtools::SineFitter::SineFitFn::~SineFitFn()
{
#if defined(SINE_SUBTRACT_USE_FLOATS) || defined(SINE_SUBTRACT_FORCE_ALIGNED)
setXY(0,0,0,0); 
#endif
}

void FFTtools::SineFitter::SineFitFn::setXY(int ntraces, int nsamples, const double ** xx, const double ** yy, const double * wset) 
{
  wgt = wset; 
#ifdef SINE_SUBTRACT_USE_FLOATS
  // we have to convert everything to floats... 
  if (nt!=ntraces)
  {

    float ** newx = (float**) malloc(ntraces * sizeof(float*)); 
    float ** newy = (float**) malloc(ntraces * sizeof(float*)); 
    for (int i = 0; i < nt; i++) 
    {
      if (i < ntraces)
      {
        newx[i] = x[i]; 
        newy[i] = y[i]; 
      }
      else
      {
        free(x[i]); 
        free(y[i]); 
      }
    }

    for (int i = nt; i < ntraces; i++) 
    {
      newx[i] = 0; 
      newy[i] = 0; 
    }


    free(x); 
    free(y); 
    x = newx; 
    y = newy; 
  }

  if (ns!=nsamples)
  {
    for (int i = 0; i < ntraces; i++) 
    {
      if (x[i]) free(x[i]); 
      x[i] = (float*) memalign(VEC_N * sizeof(VEC_T), nsamples * sizeof(float)); 
      if (y[i]) free(y[i]); 
      y[i] = (float*) memalign(VEC_N * sizeof(VEC_T), nsamples * sizeof(float)); 
    }
  }


  for (int j = 0; j < ntraces; j++) 
  {
    for (int i = 0; i < nsamples; i++) 
    {
       x[j][i] = xx[j][i]; 
       y[j][i] = yy[j][i]; 
    }
  }

  xp = xx; 
  yp = yy; 
  
#elif SINE_SUBTRACT_FORCE_ALIGNED
  if (nt!=ntraces)
  {

    double ** newx = (double**) malloc(ntraces * sizeof(float*)); 
    double ** newy = (double**) malloc(ntraces * sizeof(float*)); 
    for (int i = 0; i < nt; i++) 
    {
      if (i < ntraces)
      {
        newx[i] = x[i]; 
        newy[i] = y[i]; 
      }
      else
      {
        free(x[i]); 
        free(y[i]); 
      }
    }

    for (int i = nt; i < ntraces; i++) 
    {
      newx[i] = 0; 
      newy[i] = 0; 
    }


    free(x); 
    free(y); 
    x = newx; 
    y = newy; 
  }

  if (ns!=nsamples)
  {
    for (int i = 0; i < ntraces; i++) 
    {
      free(x[i]); 
      x[i] = (double*) memalign(VEC_N * sizeof(VEC_T), nsamples * sizeof(double)); 
      y[i] = (double*) memalign(VEC_N * sizeof(VEC_T), nsamples * sizeof(double)); 
    }
  }


  for (int j = 0; j < ntraces; j++) 
  {
    memcpy(x[j], xx[j], nsamples * sizeof(double)); 
    memcpy(y[j], yy[j], nsamples * sizeof(double)); 
  }

  xp = xx; 
  yp = yy; 
 
#else
  x = xx; 
  y = yy;
#endif
  nt = ntraces; 
  ns = nsamples;
}


//stupid function that converts an array to a string 
static char * arrString(int n, const  double * arr, int sigfigs =8 ) 
{
  int size = 2 + (n-1) + (sigfigs+5) * n + 2;
  
  char * ret = new char[size]; 
  char format[8]; 
  sprintf(format,"%%.%dg",sigfigs); 

  int ctr = 0; 
  ctr += sprintf(ret + ctr, "("); 
  for (int i = 0; i < n; i++) 
  {
    ctr+=sprintf(ret + ctr, format, arr[i]); 
    if (i < n-1)  ctr += sprintf(ret + ctr, ","); 
  }
  ctr+= sprintf(ret+ctr, ")"); 

  assert(ctr < size); 

  return ret; 
}





void FFTtools::SineFitter::doFit(int ntraces, int nsamples, const double ** x, const double **y, const double * w)
{
  f.setXY(ntraces,nsamples,x,y,w); 

  if (verbose) 
  {
    char * Phstr =arrString(ntraces,&phase[0]); 
    printf("Guesses: f= %f, A=%f, ph=%s", freq, amp[0], Phstr); 
    delete [] Phstr; 
    double p[1 + 2*ntraces]; 
    p[0] = freq; 

    for (int i =0; i < ntraces; i++) 
    {
      p[1+2*i] = phase[i]; 
      p[2+2*i] = amp[0]; 
    }
    printf("Guess power is %f. Gradient is (", f.DoEval(p));
    for (int i = 0; i < 2 * ntraces+1; i++) 
    {
      printf("%f", f.DoDerivative(p,i)); 
      if (i < 2*ntraces) printf(","); 
    }
    printf(")\n"); 

  }

  double dt = (x[0][nsamples-1]  - x[0][0]) / (nsamples-1); 
  double fnyq = 1. / (2 * dt); 
  double df = fnyq / nsamples; 




  min.SetFunction(f); 


  if (limits.max_n_df_relative_to_guess)
  {
    min.SetLimitedVariable(0, "f",freq, df/10., freq-df * limits.max_n_df_relative_to_guess, freq+df * limits.max_n_df_relative_to_guess); 
  }
  else
  {
    min.SetVariable(0, "f",freq, df/10); 
  }


  double damp =0.1* amp[0]; 


  for (int i = 0; i < ntraces; i++) 
  {
    min.SetVariable(1+2*i, TString::Format("phi%d",i).Data(),phase[i], TMath::Pi()/64); 

    if (limits.maxA_relative_to_guess && limits.minA_relative_to_guess)
    {
       min.SetLimitedVariable(2+2*i, TString::Format("A%d",i).Data(), amp[0], damp, amp[0]*limits.minA_relative_to_guess, amp[0]*limits.maxA_relative_to_guess); 
    }
    else if (limits.maxA_relative_to_guess)
    {
       min.SetUpperLimitedVariable(2+2*i, TString::Format("A%d",i).Data(), amp[0], damp, amp[0]*limits.maxA_relative_to_guess); 
    }
    else if (limits.minA_relative_to_guess)
    {
       min.SetLowerLimitedVariable(2+2*i, TString::Format("A%d",i).Data(), amp[0], damp, amp[0]*limits.minA_relative_to_guess);
    }
    else
    {
       min.SetVariable(2+2*i, TString::Format("A%d",i).Data(), amp[0], damp); 
    }
  }


  min.SetPrintLevel(verbose ? 1 : 0); 

  int old_level = gErrorIgnoreLevel; 
  if (!verbose)
  {
    gErrorIgnoreLevel = 1001; 
  }
  min.Minimize(); 
  gErrorIgnoreLevel = old_level; 
  if(verbose)  min.PrintResults(); 

  freq = min.X()[0]; 
  freq_err = min.Errors()[0]; 
  for (int i = 0; i < ntraces; i++) 
  {
    phase[i] = normalize_angle(min.X()[1+2*i]); 
    amp[i] = min.X()[2+2*i]; 
    phase_err[i] = min.Errors()[1+2*i]; 
    amp_err[i] =  min.Errors()[2+2*i]; 
  }
}


FFTtools::SineSubtract::SineSubtract(int maxiter, double min_power_reduction, bool store)
  : maxiter(maxiter), min_power_reduction(min_power_reduction), store(store)
{

  power_estimator = LOMBSCARGLE; 
  oversample_factor = 2; 
  high_factor = 1; 
  neighbor_factor2 = 0.15; 
  verbose = false; 
  tmin = 0; 
  tmax = 0; 

#ifdef ENABLE_VECTORIZE
  no_subnormals();  // Unlikely to occur, but worth it if they do
#endif

  
}

void FFTtools::SineSubtractResult::append(const SineSubtractResult *r) 
{
  powers.insert(powers.end(), r->powers.begin(), r->powers.end()); 
  freqs.insert(freqs.end(), r->freqs.begin(), r->freqs.end()); 
  freqs_errs.insert(freqs_errs.end(), r->freqs_errs.begin(), r->freqs_errs.end()); 

  for (size_t i = 0; i < phases.size(); i++) 
  {
    phases[i].insert(phases[i].end(), r->phases[i].begin(), r->phases[i].end()); 
    amps[i].insert(amps[i].end(), r->amps[i].begin(), r->amps[i].end()); 
    phases_errs[i].insert(phases_errs[i].end(), r->phases_errs[i].begin(), r->phases_errs[i].end()); 
    amps_errs[i].insert(amps_errs[i].end(), r->amps_errs[i].begin(), r->amps_errs[i].end()); 
  }
}


void FFTtools::SineSubtractResult::clear() 
{
  powers.clear(); 
  freqs.clear(); 
  phases.clear(); 
  amps.clear(); 
  freqs_errs.clear(); 
  phases_errs.clear(); 
  amps_errs.clear(); 
}


void FFTtools::SineSubtract::reset() 
{

  r.clear(); 
  for (unsigned i = 0; i < gs.size(); i++)
  {
    for (unsigned j = 0; j < gs[i].size(); j++) 
    {
      delete gs[i][j]; 
    }
  }
  gs.clear(); 
  for (unsigned i = 0; i < spectra.size(); i++) delete spectra[i]; 
  spectra.clear(); 


}

TGraph * FFTtools::SineSubtract::subtractCW(const TGraph * g, double dt) 
{
  TGraph * gcopy = new TGraph(g->GetN(), g->GetX(), g->GetY()); 
  gcopy->SetTitle(TString::Format("%s (subtracted)", g->GetTitle())); 
  subtractCW(1,&gcopy,dt); 
  return gcopy; 
}

void FFTtools::SineSubtract::subtractCW(int ntraces, TGraph ** g, double dt, const double * w) 
{


#ifdef SINE_SUBTRACT_PROFILE
  TStopwatch sw; 
#endif
  reset(); 

  std::multiset<int> failed_bins; //nfails for each bin 


  int low = tmin < 0 || tmin >= g[0]->GetN() ? 0 : tmin; 
  int high = tmax <= 0 || tmax > g[0]->GetN() ? g[0]->GetN() : tmax; 


  int Nuse = high - low; 

  //zero mean and compute power at the same time. 
  double power = 0; 
  for (int ti = 0; ti < ntraces; ti++) 
  {
    double mean = g[ti]->GetMean(2); 
    for (int i = low; i <high ; i++)
    {
      g[ti]->GetY()[i] -=mean; 
      power += g[ti]->GetY()[i] * g[ti]->GetY()[i]; 
    }
  }

  r.powers.push_back(power/Nuse/ntraces); 

  if (store) 
  {
    gs.insert(gs.end(), ntraces, std::vector<TGraph*>()); 

    for (int i = 0; i < ntraces; i++) 
    {
      g[i]->SetTitle(TString::Format("Initial Waveform %d",i)); 
      gs[i].push_back(new TGraph(*g[i])); 
    }
  }

  int ntries = 0; 

  int spectrum_N = power_estimator == FFT ? Nuse/2 + 1 : Nuse * oversample_factor * high_factor / 2; 

  TGraph* power_spectra[ntraces]; 
  TGraph* fft_phases[ntraces]; // not all power estimation options use this
  for (int i = 0; i < ntraces; i++) 
  {
    power_spectra[i] = new TGraph(spectrum_N);  
    fft_phases[i] = 0; 
  }


  r.phases.insert(r.phases.end(),ntraces, std::vector<double>()); 
  r.phases_errs.insert(r.phases_errs.end(),ntraces, std::vector<double>()); 
  r.amps.insert(r.amps.end(),ntraces, std::vector<double>()); 
  r.amps_errs.insert(r.amps_errs.end(),ntraces, std::vector<double>()); 

  double fnyq = 1./(2*dt); 
  int nattempts = 0; 
  while(true) 
  {

    nattempts++; 
    for (int ti = 0; ti  < ntraces; ti++)
    {

      if (power_estimator == LOMBSCARGLE)
      {
        double hf = high_factor; 
        if (fmax.size()) hf = std::min(*std::max_element(fmax.begin(), fmax.end()) / fnyq, high_factor); 
        FFTtools::lombScarglePeriodogram(Nuse, dt, g[ti]->GetX() + low, g[ti]->GetY() + low, oversample_factor, hf, power_spectra[ti]);
      }
      else  //FFT 
      {
        double df = 1./(Nuse * dt); 
        if (fft_phases[ti] == 0)
        {
          fft_phases[ti] = new TGraph(spectrum_N); 
        }

        TGraph * ig = g[ti]; 
        if (dt > 0) // must interpolate
        {
          ig = FFTtools::getInterpolatedGraph(g[ti], dt); 
          ig->Set(g[0]->GetN());  // ensure same length
        }

        FFTWComplex * the_fft = FFTtools::doFFT(Nuse, ig->GetY() + low); 

        for (int i = 0; i < spectrum_N; i++)
        {
          power_spectra[ti]->GetY()[i] = the_fft[i].getAbsSq() / Nuse / 16; 
          if (i > 0 && i <spectrum_N-1) power_spectra[ti]->GetY()[i] *=2; 
          fft_phases[ti]->GetY()[i] = the_fft[i].getPhase(); 
          power_spectra[ti]->GetX()[i] = df *i; 
          fft_phases[ti]->GetX()[i] = df *i; 
        }


        delete [] the_fft; 

        if (dt > 0)
        {
          delete ig; 
        }
      }
    }


    int max_i = -1; 

    double mag_sum[spectrum_N] __attribute__((aligned(32))); 
    memset(mag_sum, 0, sizeof(mag_sum)); 
    double max_adj_mag2 = 0; 
    double max_f = 0; 

    double * spectra_x=0, *spectra_y=0; 

    if (store)
    {
      spectra_x  = new double[spectrum_N]; 
      spectra_y  = new double[spectrum_N]; 
    }

    //sum up all the spectra
    
    for (int ti = 0; ti < ntraces; ti++) 
    {
      for (int i = 0; i < spectrum_N; i++)
      {
         mag_sum[i] += power_spectra[ti]->GetY()[i];
      }
    }


    double df = power_spectra[0]->GetX()[1] - power_spectra[0]->GetX()[0]; 

    // find the maximum! 
    for (int i = 0; i < spectrum_N; i++)
    {
      double freq =  power_spectra[0]->GetX()[i]; 

      double adj_mag2 = mag_sum[i] / (1+failed_bins.count(i)); 
      double neigh_mag2_low = i == 0 ? DBL_MAX : mag_sum[i-1];; 
      double neigh_mag2_high = i == spectrum_N-1 ? DBL_MAX : mag_sum[i+1]; 

//      printf("%f %f\n", freq, adj_mag2); 

      if (store) 
      {
        spectra_x[i] = freq; 
        spectra_y[i] = mag_sum[i]/ntraces; 
      }

      //check if within allowed frequencies
      if (fmin.size())
      {
        bool ok = false; 
        for (unsigned i = 0; i < fmin.size(); i++) 
        {
          if (freq+df >= fmin[i] && freq-df <= fmax[i]) 
          {
            ok = true; 
            break; 
          }
        }

        if (!ok) 
        {
          continue; 
        }
      }
        

      //we need a smarter peak finder here... 
      if ( max_i < 0 ||
         ( adj_mag2 > max_adj_mag2 && 
           mag_sum[i] * neighbor_factor2 > TMath::Min(neigh_mag2_low, neigh_mag2_high)
           )
         )
      {
        
        max_i = i; 
        max_adj_mag2 =adj_mag2; 
        max_f = freq;
      }
    }


    double guess_ph[ntraces]; 
    double guess_A = 0; 

    for (int ti = 0; ti < ntraces; ti++)
    {
      guess_ph[ti] = fft_phases[ti] ? fft_phases[ti]->GetY()[max_i] : 0; 

      guess_A += sqrt(power_spectra[ti]->GetY()[max_i]) / ntraces; 
    }



    fitter.setGuess(max_f, ntraces, guess_ph, guess_A); 

    const double * x[ntraces];
    const double * y[ntraces];

    for (int i = 0; i < ntraces; i++) 
    {
      x[i] = g[i]->GetX()+low; 
      y[i] = g[i]->GetY()+low; 
    }

    fitter.doFit(ntraces, Nuse, x,y,w); 

//    printf("power before:%f\n", power); 
    power = fitter.getPower(); 
//    printf("power after:%f\n", power); 

    double ratio = 1. - power /r.powers[r.powers.size()-1]; 
    if (verbose) printf("Power Ratio: %f\n", ratio); 

    if(store && (ratio >= min_power_reduction || ntries == maxiter))
    {
      spectra.push_back(new TGraph(spectrum_N,spectra_x, spectra_y)); 
      if (r.powers.size() > 1)
      {
        spectra[spectra.size()-1]->SetTitle(TString::Format("Spectrum after %lu iterations",r.powers.size()-1)); 
      }
      else
      {
        spectra[spectra.size()-1]->SetTitle("Initial spectrum"); 
      }

      delete [] spectra_x; 
      delete [] spectra_y; 
    }



    if (ratio < min_power_reduction)
    {
      failed_bins.insert(max_i); 
      if (ntries++ > maxiter)   
      {
        break;
      }
      continue; 
    }
    ntries = 0; // restart clock

    r.powers.push_back(power); 
    r.freqs.push_back(fitter.getFreq()); 
    r.freqs_errs.push_back(fitter.getFreqErr()); 


    for (int i = 0; i < ntraces; i++) 
    {
      r.phases[i].push_back(fitter.getPhase()[i]); 
      r.phases_errs[i].push_back(fitter.getPhaseErr()[i]); 
      r.amps[i].push_back(fitter.getAmp()[i]); 
      r.amps_errs[i].push_back( fitter.getAmpErr()[i]); 
    }

    for (int ti = 0; ti < ntraces; ti++) 
    {
      double A = w ? w[ti] * fitter.getAmp()[ti] : fitter.getAmp()[ti]; 
      for (int i = 0; i < g[ti]->GetN(); i++) 
      {
        g[ti]->GetY()[i] -= A* sin(2*TMath::Pi() * fitter.getFreq() *g[ti]->GetX()[i] + fitter.getPhase()[ti]); 
      }

      if (store) 
      {
        //add sine we subtracted to previous graph 
        TF1 * fn = new TF1(TString::Format("fitsin%lu_%d",r.powers.size(),ti), "[0] * sin(2*pi*[1] * x + [2])",g[ti]->GetX()[0], g[ti]->GetX()[g[ti]->GetN()-1]); 
        fn->SetParameter(0,fitter.getAmp()[ti]); 
        fn->SetParameter(1,fitter.getFreq()); 
        fn->SetParameter(2,fitter.getPhase()[ti]); 
        fn->SetNpx(2*g[ti]->GetN()); 
        gs[ti][gs[ti].size()-1]->GetListOfFunctions()->Add(fn); 

        //add new graph 
        g[ti]->SetTitle(TString::Format("Wf %d after %lu iterations",ti, r.powers.size()-1)); 
        gs[ti].push_back(new TGraph(*g[ti])); 
      }

    }
  }

  for (int i = 0; i < ntraces; i++) 
  { 
    delete power_spectra[i];
    if (fft_phases[i]) delete fft_phases[i]; 
  }
#ifdef SINE_SUBTRACT_PROFILE
  printf("Time for SineSubtract::subtractCW(): "); 
  sw.Print("u"); 
  printf("nattempts: %d\n",nattempts); 

#endif
}


int FFTtools::SineSubtract::getNSines() const 
{
  return r.phases[0].size(); 
}

void FFTtools::SineSubtract::makeSlides(const char * title, const char * fileprefix, const char * outdir , const char * format) const 
{

  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);

  Int_t orig_msz = gStyle->GetMarkerSize();
  Int_t orig_mst = gStyle->GetMarkerStyle();
  Int_t orig_lt  = gStyle->GetLineWidth();

  gStyle->SetMarkerSize(2);
  gStyle->SetMarkerStyle(20);
  gStyle->SetLineWidth(orig_lt*4);

  int orig_stat = gStyle->GetOptStat(); 
  gStyle->SetOptStat(0); 

 
  TCanvas canvas("slidecanvas","SlideCanvas",4000,3000); 
  canvas.Divide(3,1); 

  TGraph g; 
  int niter = spectra.size(); 
  TH2I poweraxis("poweraxis","Iteration 0", 10, 0,niter, 10,0, r.powers[0]*1.1); 
  poweraxis.GetXaxis()->SetTitle("iteration"); 
  poweraxis.GetYaxis()->SetTitle("power"); 
  if (gs.size() > 1) 
  {
    canvas.cd(2)->Divide(1,gs.size()); 
  }

  FILE * texfile = fopen(TString::Format("%s/%s.tex", outdir, fileprefix),"w"); 

  for (int i = 0; i < niter; i++) 
  {

    canvas.cd(1); 
    poweraxis.SetTitle(i < niter -1 ? TString::Format("Iteration %d, Sub Freq = %g Ghz",i, r.freqs[i]) : TString::Format("Iteration %d (final)", i)); 

    poweraxis.Draw(); 
    g.SetPoint(i, i, r.powers[i]); 
    g.Draw("lp"); 
    canvas.cd(2); 

    int old_gs_width[gs.size()];

    for (size_t j = 0; j < gs.size(); j++) 
    {
      if (gs.size() > 1) 
      {
        canvas.cd(2)->cd(j+1); 
      }
      old_gs_width[j]= gs[j][i]->GetLineWidth(); 
      gs[j][i]->SetLineWidth(3); 
      gs[j][i]->Draw("al"); 
    }


    canvas.cd(3)->SetLogy(); 

    int old_spectrum_width = spectra[i]->GetLineWidth(); 
    spectra[i]->SetLineWidth(5); 
    spectra[i]->Draw("alp"); 


    TString canvfile = TString::Format("%s/%s_%d.%s", outdir, fileprefix, i, format); 
    canvas.SaveAs(canvfile); 

    spectra[i]->SetLineWidth(old_spectrum_width); 
    for (size_t j = 0; j < gs.size(); j++) 
    {
      gs[j][i]->SetLineWidth(old_gs_width[j]); 
    }

    fprintf(texfile, "\\begin{frame}\n"); 
    fprintf(texfile, "\t\\frametitle{%s (iteration %d)}\n", title, i); 
    fprintf(texfile, "\t\\begin{center}\n"); 
    fprintf(texfile, "\t\t\\includegraphics[width=4.2in]{%s_%d}\n",fileprefix,i); 
    fprintf(texfile, "\t\\end{center}\n"); 
    fprintf(texfile, "\\end{frame}\n\n"); 
    fflush(texfile); 
  }

  fclose(texfile); 

  gROOT->ForceStyle(kFALSE);
  gROOT->SetBatch(kFALSE);

  gStyle->SetMarkerSize(orig_msz);
  gStyle->SetMarkerStyle(orig_mst);
  gStyle->SetLineWidth(orig_lt);
  gStyle->SetOptStat(orig_stat); 
}

void FFTtools::SineSubtract::makePlots(TCanvas * cpower, TCanvas * cw, int ncols) const 
{

  int nplots = spectra.size() * (1 + gs.size()); 
  if (!cpower) cpower = new TCanvas("sinesubtractpowerplot","SineSubtract Power Evolution",1000,600); 
  cpower->Clear(); 

  if (!cw) cw = new TCanvas("sinesubtractplots","SineSubtract Waveforms/Spectra",1000,600); 
  cw->Clear(); 
  cw->Divide(ncols, (nplots-1)/ncols +1); 


  std::vector<double> iterations(r.powers.size()); 
  for (unsigned i = 0; i < r.powers.size(); i++) iterations[i] = i;
  TGraph *gpower = new TGraph(r.powers.size(), &iterations[0], &r.powers[0]); 
  cpower->cd(); 
  gpower->SetTitle("Power vs. iteration"); 
  gpower->Draw("alp"); 


  for (size_t i = 0; i < spectra.size(); i++) 
  {
    cw->cd((1+gs.size())*i+1)->SetLogy(); 
    ((TGraph*)spectra[i])->Draw("alp"); 
    for (size_t j = 0; j < gs.size(); j++)
    {
      cw->cd((1+gs.size())*i+2+j); 
      ((TGraph*)gs[j][i])->Draw("alp"); 
    }
  }

}


FFTtools::SineSubtract::~SineSubtract()
{
  reset(); 
}


ClassImp(FFTtools::SineSubtractResult); 


