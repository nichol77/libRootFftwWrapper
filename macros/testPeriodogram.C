TGraph * testPeriodogram(int oversample_factor = 4, double high_factor =1)
{
  FFTtools::loadWisdom("wisdom.dat"); 
  const int N = 256;  //nsamples
  double dt = 1./2.6; // mean sample period
  double jitter = 0.1; //sample timing jitter in nanoseconds
  double rms = 3; //noise rms 
  const int nfreq = 20; //number of CW
  double min_freq = 0.20; 
  double max_freq = 0.5; 
  double min_noise_freq = 0.2;
  double max_noise_freq = 1.2; 
  double max_amp = 3.5; //max CW amplitude 
  double min_amp = 6.5; //min CW amplitude 
  gRandom->SetSeed(11); // Random seed (0 will pick a new one, do something else if you want something reproducible) 
 
  double f[nfreq]; 
  double A[nfreq]; 
  double ph[nfreq]; 
  double fnyq = 1./(2*dt); 

  for (int i = 0; i < nfreq; i++) 
  {
    {
      f[i] = gRandom->Uniform( min_freq, max_freq); 
      A[i] = gRandom->Uniform(min_amp,max_amp); 
      ph[i] = gRandom->Uniform(-TMath::Pi(),TMath::Pi()); 
      printf("True CW %d: f = %f, A = %f, ph = %f\n", i, f[i], A[i], ph[i]); 
    }
  }


  FFTtools::ThermalNoise noise(N*10, min_noise_freq/fnyq, max_noise_freq/fnyq, rms, 2); 
  TGraph * g = noise.makeGaussianDtWf(N, dt, jitter); 
  for (int i = 0; i < N; i++) 
  {
    for (int j = 0; j < nfreq; j++) 
    {
         g->GetY()[i] += A[j] * TMath::Sin( 2*TMath::Pi()*f[j]* (g->GetX()[i] + ph[j])); 
    }
  }

  TGraph * p = FFTtools::lombScarglePeriodogram(g,oversample_factor,high_factor); 

  p->Draw(); 
  FFTtools::saveWisdom("wisdom.dat"); 
  return p; 
}



 

