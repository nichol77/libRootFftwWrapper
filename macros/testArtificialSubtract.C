{
  FFTtools::loadWisdom("wisdom.dat"); 

  bool reinterpolate = true;  
  const int N = 256;  //nsamples
  const int ntraces = 3; 
  double dt = 1./2.6; // mean sample period
  double jitter = 0.1; //sample timing jitter in nanoseconds
  double rms = 3; //noise rms 
  const int nfreq = 5; //number of CW
  double min_freq = 0.20; 
  double max_freq = 0.5; 
  double min_noise_freq = 0.2;
  double max_noise_freq = 1.2; 
  double min_amp = 3.5; //max CW amplitude 
  double max_amp = 6.5; //min CW amplitude 
  int max_failed_iterations = 5; //settings for sinesubtract 
  double min_power_ratio = 0.05; 


  bool use_freq_limits = false; 
  bool enable_trace_limits = false; 
  int trace_min = 0; 
  int trace_max = 450; 


  bool verbose = true; 
  gRandom->SetSeed(11); // Random seed (0 will pick a new one, do something else if you want something reproducible) 

  //end config



  double f[nfreq]; 
  double A[ntraces][nfreq]; 
  double ph[ntraces][nfreq]; 
  double fnyq = 1./(2*dt); 

  for (int i = 0; i < nfreq; i++) 
  {
    f[i] = gRandom->Uniform( min_freq, max_freq); 
    double  amp = gRandom->Uniform(min_amp, max_amp); 
    for (int j = 0; j < ntraces; j++)
    {
      A[j][i] = amp * gRandom->Uniform(0.9,1.1); 
      ph[j][i] = gRandom->Uniform(-TMath::Pi(),TMath::Pi()); 
  //    printf("True CW %d: f = %f, A = %f, ph = %f\n", i, f[i], A[i], ph[i]); 
    }
  }

  TGraph * g[ntraces]; 

  for (int ti = 0; ti < ntraces; ti++)
  {

    FFTtools::ThermalNoise noise(N*10, min_noise_freq/fnyq, max_noise_freq/fnyq, rms, 2); 
    g[ti] = reinterpolate ? noise.makeGaussianDtWf(N, dt, jitter) : noise.makeEvenWf(N,dt); 
    for (int i = 0; i < N; i++) 
    {
      for (int j = 0; j < nfreq; j++) 
      {
         g[ti]->GetY()[i] += A[ti][j] * TMath::Sin( 2*TMath::Pi()*f[j]* (g[ti]->GetX()[i] + ph[ti][j])); 
      }
    }

  }


  FFTtools::SineSubtract * sub = new FFTtools::SineSubtract(max_failed_iterations,min_power_ratio,true); 
  sub->setVerbose(verbose); 

  if (use_freq_limits)
  {
    sub->setFreqLimits(min_freq, max_freq); 
  }

  sub->subtractCW(ntraces, g, reinterpolate ?  dt : 0); 

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
        printf("\t(trace %d) A= %f+/-%f, ph = %f+/-%f\n", j+1, recoA[j], recoAErr[j], recoPh[j], recoPhErr[j]); 
      }
    }
  }
  FFTtools::saveWisdom("wisdom.dat"); 
}
