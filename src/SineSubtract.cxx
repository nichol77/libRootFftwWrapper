#ifndef DONT_HAVE_MINUIT2
#include "SineSubtract.h" 


#include "TGraph.h" 
#include "TCanvas.h" 
#include "TStyle.h" 
#include "TMutex.h" 
#include "TLinearFitter.h" 
#include "DigitalFilter.h" 
#include "TSpectrum.h" 

#ifndef __APPLE__
#include <malloc.h>
#define SINCOS sincos 
#else
#include <malloc/malloc.h>
#define SINCOS __sincos
#endif

#ifdef SINE_SUBTRACT_PROFILE
#include "TStopwatch.h"
#endif

#include "TMath.h"
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



static TMutex fnLock; 
static int fnCount = 0; 
static int counter; 




#include "TError.h" 
//extern int gErrorIgnoreLevel; // who ordered that? 

static double guessPhase(const TGraph * g, double freq) 
{
  double phase; 
  FFTtools::dftAtFreq(g,freq,&phase,0); 
  return phase; 
}



static double normalize_angle(double phi)
{
  return  phi - 2*TMath::Pi() * FFTtools::fast_floor((phi+ TMath::Pi()) / (2.*TMath::Pi())); 
}


static __thread TLinearFitter * fitter = 0; 
static __thread int fitter_order = 0; 



static const char * order_strings[] = { "1", "1++x", "1++x++x*x","1++x++x*x++x*x*x","1++x++x*x++x*x*x++x*x*x*x" }; 


static void computeEnvelopeFit(int order, int N, const double * x, const double * y, double * out) 
{

  if (!fitter || fitter_order != order) 
  {
    fnLock.Lock(); 
    if (!fitter) 
    {
      fitter = new TLinearFitter(1,order_strings[order],""); 
    }
    else
    {
      fitter->SetFormula(order_strings[order]); 
    }
    fnLock.UnLock(); 
  }

  fitter->ClearPoints(); 
  fitter->AssignData(N,1,(double*)x,(double*)y); 
  fitter->Eval(); 

  double max = 0; 
  for (int i = 0; i < N; i++)
  {
    out[i] = fitter->GetParameter(0);
    double xx = 1; 
    for (int j = 1; j <= order; j++) 
    {
      xx *= x[i]; 
      out[i] += fitter->GetParameter(j) * xx; 
    }
    if (out[i] > max) max = out[i]; 
  }

  for (int i = 0; i < N; i++) 
  {
    out[i]/=max; 
  }
}

/* wrapper around hilbert envelope because need to interpolate to calculate it properly */
static void computeEnvelopeHilbert(const TGraph * g, double * out, double dt) 
{

  TGraph * ig = FFTtools::getInterpolatedGraph(g, dt); 
  TGraph * hilbert = FFTtools::getHilbertEnvelope(ig); 
  double xmax = ig->GetX()[ig->GetN()-1]; 

  for (int i = 0; i < g->GetN(); i++)
  {
    if (g->GetX()[i] < xmax)
    {
      out[i] = FFTtools::evalEvenGraph(hilbert, g->GetX()[i]); 
    }
    else
    {
      out[i] = ig->GetY()[ig->GetN()-1] ; 
    }
  }

  delete ig; 
  delete hilbert; 
}



FFTtools::SineFitter::SineFitter() : fDoEvalRecord(false), f(this)
{
  min.SetFunction(f);
  verbose = false; 

}


void FFTtools::SineFitter::deleteEvalRecords(){
  for(unsigned i=0; i < grEvalRecords.size(); i++){
    if(grEvalRecords[i]){
      delete grEvalRecords[i];
      grEvalRecords[i] = NULL;
    }
  }
  grEvalRecords.clear();
}

FFTtools::SineFitter::~SineFitter()
{

  deleteEvalRecords();
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
  fn->setXY(nt,ns,xp,yp,wgt,envp); //incredibly inefficient, but quick kludge for now. We're probalby never going to clone anyway. 
#else
  fn->setXY(nt,ns,x,y,wgt,env); 
#endif
  fn->fContainer = fContainer;
  return fn; 
}

FFTtools::SineFitter::SineFitFn::SineFitFn(SineFitter* parent)
    : ROOT::Math::IGradientFunctionMultiDim(), fContainer(parent)
{
  ns = NULL; 
  nt = 0; 
  x= 0; 
  y = 0; 
  env = 0; 
#if defined(SINE_SUBTRACT_USE_FLOATS) || defined(SINE_SUBTRACT_FORCE_ALIGNED)
  xp = 0; 
  yp = 0; 
  envp = 0; 
#endif

}


//here's hoping... 
__attribute__((optimize("unroll-loops")))
double FFTtools::SineFitter::SineFitFn::DoDerivative(const double * p, unsigned int coord) const
{

  double w = 2 *TMath::Pi() * p[0]; 
  double deriv = 0; 
  int type = coord == 0 ? 0 : (1 + ((coord-1) % 2)); 

  for (int ti = 0; ti < nt; ti++) 
  {

    if (type > 0 && (int(coord) - 1)/2 != ti) continue; 
    __builtin_prefetch(x[ti]); 
    __builtin_prefetch(y[ti]); 
    double ph = normalize_angle(p[1 + 2*ti]); 
    double A = p[2+2*ti]; 
    if (wgt) A*=wgt[ti]; 

#ifdef ENABLE_VECTORIZE 

    VEC vecx(0); 
    VEC vecy(0); 
    VEC vece(1); 
    VEC vec_ph(ph); 
    VEC vec_w(w); 
    VEC vec_A(A); 
    VEC vec_cos(0); 
    VEC dYdp(0); 
    int leftover = ns[ti] % VEC_N;
    int nit = ns[ti]/VEC_N + (leftover ? 1 : 0); 

#ifdef SINE_SUBTRACT_FORCE_ALIGNED
    double * xx = (double*) __builtin_assume_aligned(x[ti], VEC_N * sizeof(VEC_T)); 
    double * yy = (double*) __builtin_assume_aligned(y[ti], VEC_N * sizeof(VEC_T)); 
    double * envelope =   (env && env[ti]) ?  (double*) __builtin_assume_aligned(env[ti], VEC_N * sizeof(VEC_T) : 0; 
#elif defined(SINE_SUBTRACT_USE_FLOATS)
    float * xx = (float*) __builtin_assume_aligned(x[ti], VEC_N * sizeof(VEC_T)); 
    float * yy = (float*) __builtin_assume_aligned(y[ti], VEC_N * sizeof(VEC_T)); 
    float * envelope =   (env && env[ti]) ?  (float*) __builtin_assume_aligned(env[ti], VEC_N * sizeof(VEC_T) : 0; 
#else
    double * xx = (double*) x[ti]; 
    double * yy = (double*) y[ti]; 
    double * envelope =   (env && env[ti]) ? (double*) env[ti] : 0; 
#endif
    

    for (int i = 0; i < nit; i++)
    {
      if (i < nit-1 || !leftover)
      {
#if defined(SINE_SUBTRACT_FORCE_ALIGNED) || defined(SINE_SUBTRACT_USE_FLOATS)
        vecx.load_a(xx+VEC_N*i); 
        vecy.load(yy+VEC_N*i); 
#else
        vecx.load(xx+VEC_N*i); 
        vecy.load(yy+VEC_N*i); 
#endif
      }
      else
      {
        vecx.load_partial(leftover, xx+VEC_N*i); 
        vecy.load_partial(leftover, yy+VEC_N*i); 
      }

      VEC vec_ang = mul_add(vecx, vec_w, vec_ph); 
      VEC vec_sin = sincos(&vec_cos, vec_ang); 

      if (envelope)
      {

        if (i < nit -1 || !leftover)
        {
          vece.load(envelope + VEC_N*i); 
        }
        else
        {
          vece.load_partial(leftover, envelope + VEC_N*i); 
        }
        vec_sin *= vece; 
        vec_cos *= vece; 
      }


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

      deriv +=horizontal_add( vec_Y * dYdp)/ns[ti]; 
#else

      VEC ans = vec_Y * dYdp; 
      VEC_T ans_v[VEC_N] __attribute__((aligned(sizeof(VEC))); 
      ans.store_a(ans_v); 

      int vecn = (i == nit-1 && leftover) ? leftover : VEC_N;
      for (int j = 0; j < vecn; j++) 
      {
        deriv += ans_v[j]/ns[ti]; 
      }

#endif
    }

#else
    for (int i = 0; i < ns[ti]; i++)
    {
      double t = x[ti][i]; 
      double sinang = sin(w*t + ph); 
      double Y = A *sinang; 
      if (env && env[ti]) Y*= env[ti][i]; 

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
      if (env && env[ti]) sinang*= env[ti][i]; 

      deriv += ((Y-y[ti][i]) * dYdp)/ns[ti]; 
    }
#endif
  }

  deriv *=2; 

//  printf("dP/d%s (f=%f, ph=%f, A=%f) = %f\n", coord == 0 ? "f" : coord == 1? "ph" : "A", p[0], ph, A, deriv);  
  return deriv/nt; 
}



// yeah, like the compiler is really going to do what I want ... may have to unroll this manually 
__attribute__((optimize("unroll-loops")))
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
    VEC vece(1); 
    VEC vec_ph(ph); 
    VEC vec_w(w); 
    VEC vec_A(A); 

#ifdef SINE_SUBTRACT_FORCE_ALIGNED
    double * xx = (double*) __builtin_assume_aligned(x[ti], VEC_N * sizeof(VEC_T)); 
    double * yy = (double*) __builtin_assume_aligned(y[ti], VEC_N * sizeof(VEC_T)); 
    double * envelope =   (env && env[ti]) ?  (double*) __builtin_assume_aligned(env[ti], VEC_N * sizeof(VEC_T) : 0; 
#elif defined(SINE_SUBTRACT_USE_FLOATS)
    float * xx = (float*) __builtin_assume_aligned(x[ti], VEC_N * sizeof(VEC_T)); 
    float * yy = (float*) __builtin_assume_aligned(y[ti], VEC_N * sizeof(VEC_T)); 
    float * envelope =   (env && env[ti]) ?  (float*) __builtin_assume_aligned(env[ti], VEC_N * sizeof(VEC_T) : 0; 
#else
    double * xx = (double*) x[ti]; 
    double * yy = (double*) y[ti]; 
    double * envelope =   (env && env[ti]) ?  (double*) env[ti] : 0; 
#endif


    int leftover = ns[ti] % VEC_N;
    int nit = ns[ti]/VEC_N + (leftover ? 1 : 0); 

    for (int i = 0; i < nit; i++)
    {
      if (i < nit -1 || !leftover)
      {
        vecx.load(xx+VEC_N*i); 
        vecy.load(yy+VEC_N*i); 
      }
      else
      {
        vecx.load_partial(leftover, xx+VEC_N*i); 
        vecy.load_partial(leftover, yy+VEC_N*i); 
      }


      VEC vec_ang = mul_add(vecx, vec_w, vec_ph); 
      VEC vec_sin = sin(vec_ang); 

      //TODO if branch-prediction hurts, we should have two versions of this loop 
      if (envelope)
      {

        if (i < nit -1 || !leftover)
        {
          vece.load(envelope + VEC_N*i); 
        }
        else
        {
          vece.load_partial(leftover, envelope + VEC_N*i); 
        }
        vec_sin *= vece; 
      }

      VEC vec_Y = mul_sub(vec_sin, vec_A, vecy); 

      VEC vec_Y2 = square(vec_Y); 
#ifndef SINE_SUBTRACT_DONT_HORIZONTAL_ADD
      if (i == nit-1 && leftover) //hopefully this gets unrolled? 
      {
        vec_Y2.cutoff(leftover); 
      }

      power += horizontal_add(vec_Y2)/ns[ti]; 
#else 
      int vecn = (i == nit-1 && leftover) ? leftover : VEC_N;
      for (int j = 0; j< vecn; j++)
      {
        power += vec_Y2[j]/ns[ti]; 
      }
#endif
    }
#else
    for (int i = 0; i < ns[ti]; i++)
    {
      VEC_T Y = A * sin(w*x[ti][i] + ph); 
      if (env && env[ti]) Y*= env[ti][i]; 
      Y-=y[ti][i]; 
      power += Y*Y/ns[ti]; 
    }
#endif 

  }


  double evalResult = power/nt;
  if(fContainer && fContainer->GetDoEvalRecord()){
    TGraph* grEval = fContainer->getEvalRecordGraph();
    if(grEval){
      int numEvals = grEval->GetN();
      grEval->SetPoint(numEvals, numEvals, evalResult);
    }
  }

//  printf("P(f=%f, ph=%f, A=%f) = %f\n", p[0], ph, A, power);

  return evalResult;

}

FFTtools::SineFitter::SineFitFn::~SineFitFn()
{
#if defined(SINE_SUBTRACT_USE_FLOATS) || defined(SINE_SUBTRACT_FORCE_ALIGNED)
setXY(0,0,0,0); 
#endif

if (ns) delete [] ns; 

}

void FFTtools::SineFitter::SineFitFn::setXY(int ntraces, const int * nsamples, const double ** xx, const double ** yy, const double * wset, const double ** envset) 
{
  wgt = wset; 
#ifdef SINE_SUBTRACT_USE_FLOATS
  // we have to convert everything to floats... 
  if (nt!=ntraces)
  {

    float ** newx = (float**) malloc(ntraces * sizeof(float*)); 
    float ** newy = (float**) malloc(ntraces * sizeof(float*)); 
    float ** newenv = 0; 
    if (envset) 
      new_env = (float**)  emalloc(ntraces * sizeof(float*)); 
    for (int i = 0; i < nt; i++) 
    {
      if (i < ntraces)
      {
        newx[i] = x[i]; 
        newy[i] = y[i]; 
        if (env) 
          newenv[i] = env[i]; 
      }
      else
      {
        free(x[i]); 
        free(y[i]); 
        if (env) 
          free(env[i]); 
      }
    }

    for (int i = nt; i < ntraces; i++) 
    {
      newx[i] = 0; 
      newy[i] = 0; 
      if (newenv) 
        newenv[i] = 0; 
    }


    free(x); 
    free(y); 
    if (env) 
      free(env); 
    x = newx; 
    y = newy; 
    env = newenv; 
  }
   
  for (int i = 0; i < ntraces; i++) 
  {
    if (i >= nt || ns[i]!=nsamples[i])
    {
       if (x[i]) free(x[i]); 
       x[i] = (float*) memalign(VEC_N * sizeof(VEC_T), nsamples[i] * sizeof(float)); 
       if (y[i]) free(y[i]); 
       y[i] = (float*) memalign(VEC_N * sizeof(VEC_T), nsamples[i] * sizeof(float)); 
       if (env && env[i]) free(env[i]); 
       if (env) env[i] = (float*) memalign(VEC_N * sizeof(VEC_T), nsamples[i] * sizeof(float)); 
    }
  }


  for (int j = 0; j < ntraces; j++) 
  {
    for (int i = 0; i < nsamples[j]; i++) 
    {
       x[j][i] = xx[j][i]; 
       y[j][i] = yy[j][i]; 
       if (env && env[j])
         env[j][i] = envset[j][i];  
    }
  }

  xp = xx; 
  yp = yy; 
  envp = envset; 
  
#elif SINE_SUBTRACT_FORCE_ALIGNED
  if (nt!=ntraces)
  {

    double ** newx = (double**) malloc(ntraces * sizeof(float*)); 
    double ** newy = (double**) malloc(ntraces * sizeof(float*)); 

    double ** newenv = 0;

    if (envset) 
      newenv = (double**) malloc(ntraces * sizeof(float*)); 

    for (int i = 0; i < nt; i++) 
    {
      if (i < ntraces)
      {
        newx[i] = x[i]; 
        newy[i] = y[i]; 
        if (env) newenv[i] = env[i]; 
      }
      else
      {
        free(x[i]); 
        free(y[i]); 
        if (env) free(env[i]); 
      }
    }

    for (int i = nt; i < ntraces; i++) 
    {
      newx[i] = 0; 
      newy[i] = 0; 
      if (newenv) newenv[i] = 0; 
    }


    free(x); 
    free(y); 
    if (env) 
      free(env); 
    x = newx; 
    y = newy; 

    env = newenv; 
  }

  for (int i = 0; i < ntraces; i++) 
  {
    if (i >= nt || ns[i]!=nsamples[i])
    {
       if (x[i]) free(x[i]); 
       x[i] = (double*) memalign(VEC_N * sizeof(VEC_T), nsamples[i] * sizeof(double)); 
       if (y[i]) free(y[i]); 
       y[i] = (double*) memalign(VEC_N * sizeof(VEC_T), nsamples[i] * sizeof(double)); 
       if (env && env[i]) free(env[i]); 
       if (env) env[i] = (double*) memalign(VEC_N * sizeof(VEC_T), nsamples[i] * sizeof(double)); 
    }
  }



  for (int j = 0; j < ntraces; j++) 
  {
    memcpy(x[j], xx[j], nsamples[j] * sizeof(double)); 
    memcpy(y[j], yy[j], nsamples[j] * sizeof(double)); 
    if (envset && envset[j])
    {
      memcpy(env[j], envset[j], nsamples[j] * sizeof(double)); 
    }
  }

  xp = xx; 
  yp = yy; 
  envp = env; 
 
#else
  x = xx; 
  y = yy;
  env = envset; 

#endif

  if (nsamples)
  {
    if (!ns) 
    {
      ns = new int[ntraces]; 
    }
    else if (nt != ntraces)
    {
      delete [] ns; 
      ns = new int[ntraces];  
    }
    memcpy(ns,nsamples, ntraces*sizeof(int)); 
  }
  else 
  {
    if (ns) 
    {
      delete [] ns; 
    }

    ns=0; 
  }

  nt = ntraces; 
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





void FFTtools::SineFitter::doFit(int ntraces, const int * nsamples, const double ** x, const double **y, const double * w, const double ** env)
{
  f.setXY(ntraces,nsamples,x,y,w,env); 

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

  //Estimate Nyquist to set frequency cut 
  // Right now this just uses the first graph, although it would probably be smarter
  // to use the average or something 
  
  double dt = (x[0][nsamples[0]-1]  - x[0][0]) / (nsamples[0]-1); 
  double fnyq = 1. / (2 * dt); 
  double df = fnyq / nsamples[0]; 



  min.Clear(); 

  min.SetFunction(f); 


  if (limits.max_n_df_relative_to_guess > 0)
  {
    min.SetLimitedVariable(0, "f",freq, df  * limits.freq_start_error, freq-df * limits.max_n_df_relative_to_guess, freq+df * limits.max_n_df_relative_to_guess); 
  }
  else if (limits.max_n_df_relative_to_guess == 0)
  {
    min.SetFixedVariable(0,"f",freq); 
  }
  else 
  {
    min.SetVariable(0, "f",freq, df * limits.freq_start_error); 
  }


  double damp = limits.amp_start_error* amp[0]; 


  for (int i = 0; i < ntraces; i++) 
  {
    min.SetVariable(1+2*i, TString::Format("phi%d",i).Data(),phase[i], limits.phase_start_error); 

    if (limits.maxA_relative_to_guess > 0 && limits.minA_relative_to_guess > 0)
    {
       min.SetLimitedVariable(2+2*i, TString::Format("A%d",i).Data(), amp[0], damp, amp[0]*limits.minA_relative_to_guess, amp[0]*limits.maxA_relative_to_guess); 
    }
    else if (limits.maxA_relative_to_guess > 0)
    {
       min.SetUpperLimitedVariable(2+2*i, TString::Format("A%d",i).Data(), amp[0], damp, amp[0]*limits.maxA_relative_to_guess); 
    }
    else if (limits.minA_relative_to_guess > 0)
    {
       min.SetLowerLimitedVariable(2+2*i, TString::Format("A%d",i).Data(), amp[0], damp, amp[0]*limits.minA_relative_to_guess);
    }
    else if (limits.minA_relative_to_guess < 0 && limits.maxA_relative_to_guess > 0)
    {
       min.SetVariable(2+2*i, TString::Format("A%d",i).Data(), amp[0], damp); 
    }
    else
    {
      min.SetFixedVariable(2 + 2 * i, TString::Format("A%d",i).Data(), amp[0]); 
    }
  }


  min.SetPrintLevel(verbose ? 1 : 0); 

  int old_level = gErrorIgnoreLevel; 
  if (!verbose)
  {
    gErrorIgnoreLevel = 1001; 
  }
  if(fDoEvalRecord){
    grEvalRecords.push_back(new TGraph());
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

FFTtools::SineSubtract::SineSubtract(const TGraph * gmp, int maxiter, bool store)
  : abs_maxiter(0), maxiter(maxiter), max_successful_iter(0),
  min_power_reduction(gmp->GetMean(2)), store(store),
  power_estimator_params(0), peak_option_params(0), envelope_option_params(0) 

{
  g_min_power = gmp; 
  setPowerSpectrumEstimator(FFT); 
  setPeakFindingOption(NEIGHBORFACTOR); 
  setEnvelopeOption(ENV_NONE); 
  high_factor = 1; 
  verbose = false; 
  tmin = 0; 
  tmax = 0; 
  nfail_exponent =0.5; 
  id = counter++; 
#ifdef ENABLE_VECTORIZE
  no_subnormals();  // Unlikely to occur, but worth it if they do
#endif

  
}



FFTtools::SineSubtract::SineSubtract(int maxiter, double min_power_reduction, bool store)
  : abs_maxiter(0), maxiter(maxiter), max_successful_iter(0),
    min_power_reduction(min_power_reduction), store(store),
    power_estimator_params(0), peak_option_params(0), envelope_option_params(0) 

{

  setPowerSpectrumEstimator(FFT); 
  setPeakFindingOption(NEIGHBORFACTOR); 
  setEnvelopeOption(ENV_NONE); 
  high_factor = 1; 
  verbose = false; 
  tmin = 0; 
  tmax = 0; 
  g_min_power =0; 
  nfail_exponent =0.5; 
  id = counter++; 
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
    for (unsigned j = 0; j < env_gs[i].size(); j++) 
    {
      delete env_gs[i][j]; 
    }
  }
  gs.clear(); 
  env_gs.clear(); 
  for (unsigned i = 0; i < spectra.size(); i++) delete spectra[i]; 
  spectra.clear(); 


}

TGraph * FFTtools::SineSubtract::subtractCW(const TGraph * g, double dt, const SineSubtractResult* result) 
{
  TGraph * gcopy = new TGraph(g->GetN(), g->GetX(), g->GetY()); 
  gcopy->SetTitle(TString::Format("%s (subtracted)", g->GetTitle())); 
  subtractCW(1,&gcopy,dt, NULL, result);
  return gcopy; 
}


double * FFTtools::SineSubtract::subtractCW(int N, const double * y, double dt, double * yout,  const SineSubtractResult* result) 
{
  TGraph g(N); 
  for (int i = 0; i < N;i ++) 
  {
    g.GetX()[i] = i*dt; 
    g.GetY()[i] = y[i]; 
  }
  TGraph *gptr = &g;
  subtractCW(1,&gptr,-1,NULL,result); 
  if (!yout) 
  {
    yout = new double[N]; 
  }
  memcpy(yout, g.GetY(), N*sizeof(double)); 
  return yout; 
}


void FFTtools::SineSubtract::subtractCW(int ntraces, TGraph ** g, double dt, const double * w, const SineSubtractResult* result) 
{


#ifdef SINE_SUBTRACT_PROFILE
  TStopwatch sw;
  fitter.SetDoEvalRecord(true);
#endif
  fitter.deleteEvalRecords();
  reset(); 
  double orig_dt = dt;

  if (dt<=0) dt = g[0]->GetX()[1] - g[0]->GetX()[0]; 

  int low = tmin < 0 || tmin >= g[0]->GetN() ? 0 : tmin; 
  int high[ntraces];  

  int Nuse[ntraces];

  int NuseMax = 0; 

  //zero mean and compute power at the same time. 
  double power = 0; 
  for (int ti = 0; ti < ntraces; ti++) 
  {
    high[ti] = tmax <= 0 || tmax > g[ti]->GetN() ? g[ti]->GetN() : tmax; 
    Nuse[ti] = high[ti] - low; 

    if (Nuse[ti] > NuseMax) NuseMax = Nuse[ti]; 

    double mean = g[ti]->GetMean(2); 
    for (int i = low; i <high[ti] ; i++)
    {
      g[ti]->GetY()[i] -=mean; 
      power += g[ti]->GetY()[i] * g[ti]->GetY()[i] / Nuse[ti]; 
    }
  }

  r.powers.push_back(power/ntraces);

  if(result){
    double diff = TMath::Abs(power/ntraces - result->powers.at(0));
    if(diff > 1e-12){
      static int numError = 0;
      const int maxError = 100;
      if(numError < maxError){
        std::cerr << "Warning " <<  (numError+1) << " of " << maxError << " in " << __PRETTY_FUNCTION__
                  << ",  potential mismatch between input result and calculated power. Difference is "
                  << diff << std::endl;
        std::cerr << "Will do recalculation!" << std::endl;
        numError++;
      }
    }
    result = NULL;
  }

  if (store) 
  {
    gs.insert(gs.end(), ntraces, std::vector<TGraph*>()); 
    env_gs.insert(env_gs.end(), ntraces, std::vector<TGraph*>()); 
    for (int i = 0; i < ntraces; i++) 
    {
      g[i]->SetTitle(TString::Format("Initial Waveform %d",i)); 
      gs[i].push_back(new TGraph(g[i]->GetN(), g[i]->GetX(), g[i]->GetY())); 
    }
  }

  int ntries = 0; 

  double hf = high_factor; 
  double fnyq = 1./(2*dt); 
  if (fmax.size()) hf = std::min(*std::max_element(fmax.begin(), fmax.end()) / fnyq, high_factor); 

  double pad_scale_factor = 1; 
  if (power_estimator == FFT) pad_scale_factor = (1+power_estimator_params[FFT_NPAD]); 


  int spectrum_N = power_estimator == FFT ? (pad_scale_factor*NuseMax)/2 + 1
                                          : NuseMax * power_estimator_params[LS_OVERSAMPLE_FACTOR] * hf / 2; 


  std::vector<int> nfails(spectrum_N); 

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

  int nattempts = 0; 


  /** We need padded copies of traces if they don't all have the same length
   * so that the spectrum has the same bins (or if we want to pad for extra resolution)
   *
   **/ 
  TGraph * gPadded[ntraces]; 
  memset(gPadded,0,sizeof(gPadded)); 

  // double * envelopes[ntraces] = {0};
  std::vector<double*> envelopes(ntraces, NULL);  

  for (int ti = 0; ti < ntraces; ti++)
  {
    if (Nuse[ti] < NuseMax || pad_scale_factor > 1.)
    {
      gPadded[ti] = new TGraph(NuseMax * pad_scale_factor); 

      memcpy(gPadded[ti]->GetX(), g[ti]->GetX(), Nuse[ti] *sizeof(double));
      memcpy(gPadded[ti]->GetY(), g[ti]->GetY(), Nuse[ti] *sizeof(double));

      for (int i = Nuse[ti]; i < gPadded[ti]->GetN(); i++) 
      {
        gPadded[ti]->GetX()[i] = gPadded[ti]->GetX()[i-1] + dt; 
        gPadded[ti]->GetY()[i] = 0; //should be unnecessary 
      }
    }

    if (envelope_option != ENV_NONE) 
    {
      envelopes[ti] = new double[g[ti]->GetN()]; 
    }
    else
    {
      envelopes[ti] = 0; 

    }
  }


  int rIter = 0; // only gets incremented inside this while loop in the case of non-NULL input result
  const int nRIter = result ? result->freqs.size() : 1; // effectively making this a while(true) loop for NULL input result
  while(rIter < nRIter) 
  {
    nattempts++; 

    for (int ti = 0; ti  < ntraces; ti++)
    {


      if (envelope_option == ENV_HILBERT) 
      {
         computeEnvelopeHilbert(g[ti], envelopes[ti],dt);  
      }
      else if (envelope_option == ENV_RMS)
      {
        FFTtools::rmsEnvelope(g[ti]->GetN(), envelope_option_params[ENV_RMS_WINDOW], g[ti]->GetX(), g[ti]->GetY(), envelopes[ti]);
      }
      else if (envelope_option == ENV_PEAK)
      {

        FFTtools::peakEnvelope(g[ti]->GetN(), envelope_option_params[ENV_PEAK_MINDISTANCE], g[ti]->GetX(), g[ti]->GetY(), envelopes[ti]);
      }

      /*all of the params have the peak order as the first paramter (by construction!) */
      if (envelope_option != ENV_NONE && envelope_option_params[0] >=0) 
      {
        computeEnvelopeFit((int) envelope_option_params[0], g[ti]->GetN(), g[ti]->GetX(), envelopes[ti],  envelopes[ti]); 
      }


      TGraph * take_spectrum_of_this = gPadded[ti] ? gPadded[ti] : g[ti];  

      if (power_estimator == LOMBSCARGLE)
      {

         FFTtools::lombScarglePeriodogram(NuseMax, dt, take_spectrum_of_this->GetX() + low, take_spectrum_of_this->GetY() + low, power_estimator_params[LS_OVERSAMPLE_FACTOR], hf, power_spectra[ti]);
      }
      else //FFT 
      {
        double df = 1./(NuseMax * dt * pad_scale_factor); 
        if (fft_phases[ti] == 0)
        {
          fft_phases[ti] = new TGraph(spectrum_N); 
        }

        TGraph * ig = take_spectrum_of_this; 
        if (orig_dt > 0) // must interpolate
        {
          ig = FFTtools::getInterpolatedGraph(take_spectrum_of_this, dt); 
          ig->Set(NuseMax * pad_scale_factor);  // ensure same length
        }

        FFTWComplex * the_fft = FFTtools::doFFT(NuseMax * pad_scale_factor, ig->GetY() + low); 

        for (int i = 0; i < spectrum_N; i++)
        {
          power_spectra[ti]->GetY()[i] = the_fft[i].getAbsSq() / NuseMax / pad_scale_factor; 
          if (i > 0 && i <spectrum_N-1) power_spectra[ti]->GetY()[i] *=2; 
          fft_phases[ti]->GetY()[i] = the_fft[i].getPhase(); 
          power_spectra[ti]->GetX()[i] = df *i; 
          fft_phases[ti]->GetX()[i] = df *i; 
        }


        delete [] the_fft; 

        if (orig_dt > 0)
        {
          delete ig; 
        }
      }

      
    } // got a power estimate for it in ntraces, and calculated envelope parameters!



    std::vector<double> mag_sum(spectrum_N); 

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

    
    if (store)
    {
      for (int i = 0; i < spectrum_N; i++)
      {
         spectra_x[i] = power_spectra[0]->GetX()[i]; 
         spectra_y[i] = mag_sum[i]/ntraces; 
      }
    }

    // find the maximum 
    int max_i = findMaxFreq(spectrum_N, power_spectra[0]->GetX(), &mag_sum[0], &nfails[0]); 
    double max_f = power_spectra[0]->GetX()[max_i]; 

    if (verbose) 
      printf("max_freq (ntries: %d nattempts: %d): %d %g\n",ntries, nattempts, max_i, max_f); 



    if ( (abs_maxiter > 0 && nattempts > abs_maxiter)
        || (max_successful_iter > 0 && int(r.freqs.size()) >= max_successful_iter)
        || max_i < 0) 
    {
      if (store) 
      {
        spectra.push_back(new TGraph(spectrum_N,spectra_x, spectra_y)); 
        spectra[spectra.size()-1]->SetTitle(TString::Format("Spectrum after %lu iterations",r.powers.size()-1)); 
        delete [] spectra_x; 
        delete [] spectra_y; 
      }

 
      break; 
    }



    double guess_ph[ntraces]; 
    double guess_A = 0; 

    for (int ti = 0; ti < ntraces; ti++)
    {
      guess_ph[ti] = fft_phases[ti] ? fft_phases[ti]->GetY()[max_i] : guessPhase(g[ti], max_f); 

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

    // fitter.doFit(ntraces, Nuse, x,y,w, envelope_option == ENV_NONE? 0 : (const double **) envelopes);

    if(!result){
      fitter.doFit(ntraces, Nuse, x,y,w, envelope_option == ENV_NONE? 0 : (const double **) &envelopes[0]);

      //    printf("power before:%f\n", power); 
      power = fitter.getPower(); 
      //    printf("power after:%f\n", power);


      double ratio = 1. - power /r.powers[r.powers.size()-1]; 
      if (verbose) printf("Power Ratio: %f\n", ratio); 


      double mpr = g_min_power ? g_min_power->Eval(max_f) : min_power_reduction; 


      if(store && (ratio >= mpr || ntries == maxiter))
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


    


      if (ratio < mpr || std::isnan(mpr))
      {
        nfails[max_i]++; 
        if ((++ntries) > maxiter)   
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
    } // if(!result)
    else{ // we have an input result, so use it
      
      r.powers.push_back(result->powers[rIter+1]);
      r.freqs.push_back(result->freqs[rIter]);
      r.freqs_errs.push_back(result->freqs_errs[rIter]); 

      for (int i = 0; i < ntraces; i++) 
      {
        r.phases[i].push_back(result->phases[i][rIter]); 
        r.phases_errs[i].push_back(result->phases_errs[i][rIter]); 
        r.amps[i].push_back(result->amps[i][rIter]);
        r.amps_errs[i].push_back(result->amps_errs[i][rIter]);
      }
      rIter++; // increment this iter
    }
    

    for (int ti = 0; ti < ntraces; ti++) 
    {
      double A = w ? w[ti] * r.amps[ti].back() : r.amps[ti].back(); 
      for (int i = 0; i < g[ti]->GetN(); i++) 
      {
        double AA = envelopes[ti] ?  A * envelopes[ti][i] : A; 
        g[ti]->GetY()[i] -= AA* sin(2*TMath::Pi() * r.freqs.back() *g[ti]->GetX()[i] + r.phases[ti].back()); 
      }

      if (store)  
      {
        //add sine we subtracted to previous graph 
        
        
        fnLock.Lock(); 
        TF1 * fn = new TF1(TString::Format("fitsin_%d_%lu_%d",fnCount++, r.powers.size(),ti), "[0] * sin(2*pi*[1] * x + [2])",g[ti]->GetX()[0], g[ti]->GetX()[g[ti]->GetN()-1]); 
        fnLock.UnLock(); 
        fn->SetParameter(0,r.amps[ti].back()); 
        fn->SetParameter(1,r.freqs.back()); 
        fn->SetParameter(2,r.phases[ti].back()); 
        fn->SetNpx(2*g[ti]->GetN()); 
        gs[ti][gs[ti].size()-1]->GetListOfFunctions()->Add(fn); 

        //add new graph 
        g[ti]->SetTitle(TString::Format("Wf %d after %lu iterations",ti, r.powers.size()-1)); 
        gs[ti].push_back(new TGraph(*g[ti])); 

        if (envelopes[ti]) 
          env_gs[ti].push_back(new TGraph(g[ti]->GetN(), g[ti]->GetX(), envelopes[ti])); 
      }

    }


  }

  for (int i = 0; i < ntraces; i++) 
  { 
    if (gPadded[i]) delete gPadded[i]; 
    delete power_spectra[i];
    if (envelope_option != ENV_NONE) delete envelopes[i]; 
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
  return r.freqs.size(); 
}

void FFTtools::SineSubtract::makeSlides(const char * title, const char * fileprefix, const char * outdir , const char * format, bool standalone) const 
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

  if(standalone)
  {
    fprintf(texfile,"\\documentclass[hyperref={pdfpagelabels=false}]{beamer} \\mode<presentation> { \\usetheme{Boadilla} \\usecolortheme{beaver} }\n"); 
    fprintf(texfile, "\\setbeamertemplate{navigation symbols}{}\n"); 
    fprintf(texfile,"\\begin{document}\n"); 
  }

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

      if ((int) env_gs[j].size() > i) 
      {
        env_gs[j][i]->Draw("lsame"); 
      }
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

  if (standalone)
  {
    fprintf(texfile,"\\end{document}\n"); 
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
      if (env_gs[j].size() > i) 
        ((TGraph*)env_gs[j][i])->Draw("lsame"); 
    }
  }

}


FFTtools::SineSubtract::~SineSubtract()
{
  delete [] power_estimator_params; 
  delete [] peak_option_params; 
  reset(); 
}




bool FFTtools::SineSubtract::allowedFreq(double freq, double df) const
{
  if (!fmin.size()) return true; 
 
  for (unsigned j = 0; j < fmin.size(); j++) 
  {
     if (freq+df >= fmin[j] && freq-df <= fmax[j]) 
     {
       return true; 
     }
  }

  return false; 
}

static __thread FFTtools::SavitzkyGolayFilter* savgol = 0; 
static __thread TSpectrum* tspect = 0; 



int FFTtools::SineSubtract::findMaxFreq(int Nfreq, const double * freq, const double * mag, const int * nfails) const
{

  int max_i =-1; 
  double max = 0; 
  double df = freq[1]-freq[0]; 
       
  if (peak_option == GLOBALMAX || peak_option == NEIGHBORFACTOR)   // we can just iterate through and do this
  {
    for (int i = 0; i < Nfreq;i++)
    {

      //TODO this is a very inefficient way to do this 
      if (!allowedFreq(freq[i],df)) continue; 

      double val =mag[i]; 
      //check neighbors if we are doing that
      if (peak_option == NEIGHBORFACTOR) 
      {

        double neigh_low  = i == 0 ? DBL_MAX : mag[i-1]; 
        double neigh_high  = i == Nfreq-1 ? DBL_MAX : mag[i+1]; 
        if (val * peak_option_params[NF_NEIGHBOR_FACTOR] <= TMath::Min(neigh_low, neigh_high)) 
        {
          continue;
        }
      }

      double adjval = val;
      /* adjust the value */ 
      if (nfails[i]) adjval /= pow(1 + nfails[i],nfail_exponent); 


      if (adjval > max) 
      {
        max = adjval; 
        max_i = i; 
      }
    }
  }

  else if (peak_option == TSPECTRUM)
  {
#if ROOT_VERSION_CODE <= ROOT_VERSION(6,0,0)
    fprintf(stderr,"TSPECTRUM option not supported for ROOT < 6 due to API change\n"); 
#else
    if (!tspect) tspect = new TSpectrum(32); //32 peaks should be good enough for anyone! 

    double out[Nfreq]; 

    int npeaks = tspect->SearchHighRes((double*)mag,out,Nfreq,
        peak_option_params[TS_SIGMA], peak_option_params[TS_THRESHOLD], true,
        (int) peak_option_params[TS_NDECONV_ITERATIONS], true, (int)
        peak_option_params[TS_AVERAGE_WINDOW]) ; 


    for (int p = 0; p < npeaks; p++) 
    {
       int i= tspect->GetPositionX()[p] +0.5;  // I guess round to nearest bin? 
       if (!allowedFreq(freq[i],df)) continue; 

        double val =mag[i]; 
        if (nfails[i]) val /= pow(1 + nfails[i],nfail_exponent); 

        if (val > max) 
        {
          max = val; 
          max_i = i; 
        }
    }

    // it's possible that none of our peaks are valid, in which case we will just pick the maximum from the
    // deconvolved array :( 

    if (max_i < 0) 
    {
      for (int i = 0; i < Nfreq; i++) 
      {
        if (!allowedFreq(freq[i],df)) continue; 
        double val =mag[i]; 
        if (nfails[i]) val /= pow(1 + nfails[i],nfail_exponent); 
        if (val > max) 
        {
          max = val; 
          max_i = i; 
        }
      }
    }
#endif
  }

  else if (peak_option == SAVGOLSUB)
  {

    //check if we need a new savgol 
    if (!savgol || 
        savgol->getOrder() != (int) peak_option_params[SGS_ORDER] || 
        savgol->getWidthLeft() != (int) peak_option_params[SGS_WIDTH] )
    {
      if (savgol) delete savgol; 
      savgol = new FFTtools::SavitzkyGolayFilter(peak_option_params[SGS_ORDER], peak_option_params[SGS_WIDTH]); 
    }
    
     
    double bg[Nfreq]; 
    savgol->filterOut(Nfreq, mag, bg); 
    
    for (int i = 0; i < Nfreq; i++)
    {
      //TODO this is a very inefficient way of doing this 
      if (!allowedFreq(freq[i],df)) continue; 
      double val = mag[i] - bg[i]; 

      if (nfails[i]) val /= pow(1 + nfails[i], nfail_exponent); 
      if (val > max) 
      {
        max = val; 
        max_i = i; 
      }
    }
  }


  return max_i; 
}


static const double default_FFT_params[] = {0.} ; 
static const double default_LOMBSCARGLE_params[] = {2.} ; 
static double dummy = 0; 

void FFTtools::SineSubtract::setPowerSpectrumEstimator(PowerSpectrumEstimator estimator, const double * p) 
{
  power_estimator = estimator; 
  int N = estimator == FFT? FFT_NPARAMS : 
          estimator == LOMBSCARGLE? LS_NPARAMS :
          0; 
  
  if (power_estimator_params) delete [] power_estimator_params; 
  power_estimator_params = new double[N]; 

  if (N) 
  {
    memcpy(power_estimator_params,
           p? p : 
           estimator == FFT?  default_FFT_params  : 
           estimator == LOMBSCARGLE?  default_LOMBSCARGLE_params  : 
           &dummy,//to avoid compiler warning 
           sizeof(double) *N); 
  }

}

static const double default_NEIGHBORFACTOR_params[] = {0.15}; 
static const double default_TSPECTRUM_params[] = {2,0.05,3,3}; 
static const double default_SAVGOLSUB_params[] = {3,10}; 



void FFTtools::SineSubtract::setPeakFindingOption(PeakFindingOption option, const double * p)
{
  peak_option = option; 

  int N = option == GLOBALMAX? 0 : 
          option == NEIGHBORFACTOR? NF_NPARAMS :         
          option == TSPECTRUM? TS_NPARAMS : 
          option == SAVGOLSUB ? SGS_NPARAMS : 
          0; 

  if (peak_option_params) delete[] peak_option_params; 

  peak_option_params = new double[N]; 

  if (N) 
  {
    memcpy(peak_option_params, 
           p? p : 
           option == NEIGHBORFACTOR?  default_NEIGHBORFACTOR_params : 
           option == TSPECTRUM?  default_TSPECTRUM_params : 
           option == SAVGOLSUB?  default_SAVGOLSUB_params : 
           &dummy ,  //to avoid compiler warning
           sizeof(double) *N); 
  }

}

static const double default_ENV_HILBERT_params[] = {3}; 
static const double default_ENV_RMS_params[] = {3,5}; 
static const double default_ENV_PEAK_params[] = {3,5}; 

void FFTtools::SineSubtract::setEnvelopeOption(EnvelopeOption option, const double * p) 
{
  envelope_option = option; 

  int N = option == ENV_NONE? 0 : 
          option == ENV_HILBERT? ENV_HILBERT_NPARAMS :         
          option == ENV_RMS ? ENV_RMS_NPARAMS : 
          option == ENV_PEAK ? ENV_PEAK_NPARAMS : 
          0; 

  if (envelope_option_params) delete[] envelope_option_params; 

  envelope_option_params = new double[N]; 

  if (N) 
  {
    memcpy(envelope_option_params, 
           p? p : 
           option == ENV_HILBERT?  default_ENV_HILBERT_params : 
           option == ENV_RMS?  default_ENV_RMS_params : 
           option == ENV_PEAK?  default_ENV_PEAK_params : 
           &dummy ,  //to avoid compiler warning
           sizeof(double) *N); 
  }
}




ClassImp(FFTtools::SineSubtractResult); 

#endif

