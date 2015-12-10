{
  const int N = 256;  //nsamples
  const int ntraces = 1; 
  double dt = 0.5; // mean sample period
  double jitter = 0.1; //sample timing jitter in nanoseconds
  double rms = 3; //noise rms 
  const int nfreq = 30; //number of CW
  double min_freq = 0.24; 
  double max_freq = 0.25; 
  double max_amp = 3.5; //max CW amplitude 
  double min_amp = 6.5; //min CW amplitude 
  int max_failed_iterations = 15; //settings for sinesubtract 
  double min_power_ratio = 0.01; 


  bool use_freq_limits = true; 
  bool enable_trace_limits = false; 
  int trace_min = 0; 
  int trace_max = 450; 


  bool verbose = false; 
  gRandom->SetSeed(0); // Random seed (0 will pick a new one, do something else if you want something reproducible) 


  //end config



  double x[ntraces][N]; 
  double y[ntraces][N]; 


  double f[nfreq]; 
  double A[ntraces][nfreq]; 
  double ph[ntraces][nfreq]; 

  for (int i = 0; i < nfreq; i++) 
  {
    f[i] = gRandom->Uniform( min_freq / (2*dt), max_freq/(2*dt)); 
    for (int j = 0; j < ntraces; j++)
    {
      A[j][i] = gRandom->Uniform(min_amp,max_amp); 
      ph[j][i] = gRandom->Uniform(-TMath::Pi(),TMath::Pi()); 
  //    printf("True CW %d: f = %f, A = %f, ph = %f\n", i, f[i], A[i], ph[i]); 
    }
  }

  TGraph * g[ntraces]; 
  for (int ti = 0; ti < ntraces; ti++)
  {
    x[ti][0] = 0; 
    for (int i = 0; i < N; i++) 
    {
      if (i > 0) 
      {
        x[ti][i] = dt * i + gRandom->Uniform(-jitter,jitter); 
      }

      y[ti][i] = rms ? gRandom->Gaus(0,rms) : 0 ; 

      for (int j = 0; j < nfreq; j++) 
      {
        y[ti][i] += A[ti][j] * TMath::Sin( 2*TMath::Pi()*f[j]* (x[ti][i] + ph[ti][j])); 
      }

    }


    g[ti] = new TGraph(N,x[ti],y[ti]); 
  }

  FFTtools::SineSubtract * sub = new FFTtools::SineSubtract(max_failed_iterations,min_power_ratio,&FFTtools::HAMMING_WINDOW, true); 
  sub->setVerbose(verbose); 

  if (use_freq_limits)
  {
    sub->setFreqLimits(min_freq, max_freq); 
  }
  sub->subtractCW(ntraces, g,dt); 

  sub->makePlots(); 

  for (int i = 0; i < nfreq; i++) 
  {
    double trueA[ntraces]; 
    double truePh[ntraces]; 
    for (int j = 0; j < ntraces; j++) 
    {
      trueA[j] = A[j][i]; 
      truePh[j] = ph[j][i]; 
    }


    if (ntraces == 1) 
    {
      printf("True CW %d: f = %f, A = %f, ph = %f\n", i, f[i], trueA[0],truePh[0]); 
    }
    else
    {

      printf("True CW %d: f = %f\n", i, f[i]); 
      for (int j = 0; j < ntraces; j++) 
      {
        printf("\t(trace %d) A= %f , ph = %f\n", j+1, trueA[j], truePh[j]); 
      }
    }
  }

  for (int i = 0; i < sub->getNSines(); i++) 
  {
    double recoA[ntraces]; 
    double recoPh[ntraces]; 
    double recoAErr[ntraces]; 
    double recoPhErr[ntraces]; 

    for (int j = 0; j < ntraces; j++) 
    {
      recoA[j] = sub.getAmps(j)[i]; 
      recoAErr[j] = sub.getAmpErrs(j)[i]; 
      recoPh[j] = sub.getPhases(j)[i]; 
      recoPhErr[j] = sub.getPhaseErrs(j)[i]; 
    }

    
    if (ntraces == 1) 
    {
      printf("Reconstructed CW %d: f =%f+/-%f, A= %f+/-%f, ph = %f+/-%f\n", i, sub.getFreqs()[i], sub.getFreqErrs()[i], recoA[0], recoAErr[0], recoPh[0], recoPhErr[0]); 
    }
    else
    {
      printf("Reconstructed CW %d: f =%f+/-%f\n", i, sub.getFreqs()[i], sub.getFreqErrs()[i]); 
      for (int j = 0; j < ntraces; j++)
      {
        printf("\t(trace %d) A= %f+/-%f, ph = %f+/-%f\n", j+1, recoA[j], recoAerr[j], recoPh[j], recoPhErr[j]); 
      }
    }
  }
}
