#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class FFTWComplex+;
#pragma link C++ class RFSignal+;
#pragma link C++ class RFFilter+;


#ifdef FFTTOOLS_COMPAT_MODE
#pragma link C++ class FFTtools+;
#else
#pragma link C++ namespace FFTtools;


//window stuff
#pragma link C++ class FFTtools::FFTWindow;
#pragma link C++ class FFTtools::FFTWindowType;
#pragma link C++ class FFTtools::RectangularWindow; 
#pragma link C++ class FFTtools::TriangularWindow; 
#pragma link C++ class FFTtools::HannWindow; 
#pragma link C++ class FFTtools::HammingWindow; 
#pragma link C++ class FFTtools::BlackmanWindow; 
#pragma link C++ class FFTtools::BlackmanHarrisWindow; 
#pragma link C++ class FFTtools::KaiserWindow; 
#pragma link C++ class FFTtools::GaussianWindow; 
#pragma link C++ class FFTtools::TukeyWindow; 


#ifndef DONT_HAVE_MINUIT2
//sin subtract stuff
#pragma link C++ class FFTtools::SineSubtract;
#pragma link C++ class FFTtools::SineFitter;
#pragma link C++ class FFTtools::SineSubtractResult+;
#pragma link C++ class FFTtools::SineFitterLimits;
#endif

// filter stuff
#pragma link C++ class FFTtools::DigitalFilter; 
#pragma link C++ class FFTtools::SavitzkyGolayFilter; 
#pragma link C++ class FFTtools::GaussianFilter; 
#pragma link C++ class FFTtools::BoxFilter; 
#pragma link C++ class FFTtools::DifferenceFilter; 
#pragma link C++ class FFTtools::DigitalFilterSeries; 
#pragma link C++ class FFTtools::IIRFilter; 
#pragma link C++ class FFTtools::FIRFilter; 
#pragma link C++ class FFTtools::SincFilter; 
#pragma link C++ class FFTtools::LanczosFilter; 
#pragma link C++ class FFTtools::TransformedZPKFilter; 
#pragma link C++ class FFTtools::RCFilter; 
#pragma link C++ class FFTtools::ButterworthFilter; 
#pragma link C++ class FFTtools::ChebyshevIFilter; 
#pragma link C++ class FFTtools::ThiranFilter; 



// cwt 
#pragma link C++ class FFTtools::CWT; 
#pragma link C++ class FFTtools::CWT::MotherWavelet; 
#pragma link C++ class FFTtools::CWT::Ricker; 
#pragma link C++ class FFTtools::CWT::Ridger; 
#pragma link C++ class FFTtools::CWT::Morlet; 

// signal 
#pragma link C++ class FFTtools::AnalyticSignal; 
#pragma link C++ class FFTtools::CompositeSignal; 
#pragma link C++ class FFTtools::BandlimitedSampledSignal; 
#pragma link C++ class FFTtools::ThermalNoise; 


#pragma link C++ class FFTtools::Averager; 
#endif

#endif

