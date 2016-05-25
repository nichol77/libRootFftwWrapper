#ifndef FFTTOOLS_FFTWINDOW_H
#define FFTTOOLS_FFTWINDOW_H

class TGraph; 
#include <stddef.h>

namespace FFTtools
{


  class FFTWindowType 
  {
    
    public: 
      virtual double value(double, size_t) const = 0; 
      virtual void apply(size_t N, double * x) const; 
      virtual void fill(size_t N, double * x) const; 
      virtual double * make(size_t N) const; 
      virtual ~FFTWindowType() {; }
  }; 

  class FFTWindow 
  {
    public: 
      FFTWindow(int N, const FFTWindowType  *win) { w = win->make(N); sz = N; }
      FFTWindow(int N, const FFTWindowType & win) { w = win.make(N); sz = N; }
      virtual double value (int i) const { return w[i]; } 
      virtual ~FFTWindow() { delete w; } 
      size_t size() const { return sz; } 
    private: 
      double * w; 
      size_t sz; 
  }; 



  class RectangularWindow : public FFTWindowType
  {

    public: 
      RectangularWindow() {;} 
      virtual double value(double, size_t) const{ return 1; } 
      virtual void apply(size_t, double *) const {;} // noop 
  };

  static const RectangularWindow RECTANGULAR_WINDOW; 

  class TriangularWindow : public FFTWindowType
  {
    public: 
      TriangularWindow() {;} 
      virtual double value(double i, size_t N) const;
  };

  static const TriangularWindow TRIANGULAR_WINDOW; 

  class HannWindow  : public FFTWindowType
  {
    public: 
      HannWindow() {;} 
      virtual double value(double i, size_t N)const;
  };

  static const HannWindow HANN_WINDOW; 

  class HammingWindow : public FFTWindowType
  {
    public: 
      HammingWindow() {;} 
      virtual double value(double i, size_t N)const;
  };

  static const HammingWindow HAMMING_WINDOW; 

  class BlackmanWindow  : public FFTWindowType
  {
    public: 
      BlackmanWindow() {;} 
      virtual double value(double i, size_t N)const;
  };

  static const BlackmanWindow BLACKMAN_WINDOW; 

  class BlackmanHarrisWindow  : public FFTWindowType
  {
    public: 
      BlackmanHarrisWindow() {;} 
      virtual double value(double i, size_t N)const;
  };

  static const BlackmanHarrisWindow BLACKMAN_HARRIS_WINDOW; 

  class KaiserWindow : public FFTWindowType
  {
    public: 
      KaiserWindow(double defaultAlpha = 3) : alpha(defaultAlpha) {} 
      virtual double value(double i, size_t N)const;
    private: 
      double alpha;
  };

  static const KaiserWindow KAISER_WINDOW; 


  class GaussianWindow : public FFTWindowType
  {
    public:
      GaussianWindow(double defaultAlpha= 2.5) : alpha(defaultAlpha) {} 
      virtual double value(double i, size_t N) const; 
    private: 
      double alpha;

  }; 


  static const GaussianWindow GAUSSIAN_WINDOW; 


}

#endif
