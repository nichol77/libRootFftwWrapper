#ifndef FFTWCOMPLEX_H

#define FFTWCOMPLEX_H
#include "TMath.h"
#include <ostream>
#include <complex>

//!  This is a wrapper class for a complex number
/*!
  And that's it.
*/
class FFTWComplex{
public:
  FFTWComplex():re(0),im(0){}; ///<Default constructor
  ~FFTWComplex(){}; ///<Destructor

  double re; ///<The real part
  double im; ///<The imaginary part
  
  FFTWComplex(double real, double imag):re(real),im(imag){}

  inline FFTWComplex operator*(const FFTWComplex &rhs){
    return FFTWComplex(*this) *= rhs;
  }
  
  inline FFTWComplex& operator*=(const FFTWComplex &rhs){
    Double_t newRe = re*rhs.re - im*rhs.im;
    Double_t newIm = im*rhs.re + re*rhs.im;
    re = newRe;
    im = newIm;
    return *this;
  }
  

  inline FFTWComplex operator+(const FFTWComplex &rhs){
    return FFTWComplex(*this) += rhs;
  }
  
  inline FFTWComplex& operator+=(const FFTWComplex &rhs){
    re += rhs.re;
    im += rhs.im;
    return *this;
  }

  inline FFTWComplex operator-(const FFTWComplex &rhs){
    return FFTWComplex(*this) -= rhs;
  }
  
  inline FFTWComplex& operator-=(const FFTWComplex &rhs){
    re -= rhs.re;
    im -= rhs.im;
    return *this;
  }    

  inline FFTWComplex operator/(const FFTWComplex &rhs){
    return FFTWComplex(*this) /= rhs;
  }
  
  inline FFTWComplex& operator/=(const FFTWComplex &rhs){
    Double_t norm = rhs.getAbsSq();
    (*this) *= rhs.conj();
    re /= norm;
    im /= norm;    
    return *this;
  }    

  inline void setMagPhase(double mag, double phase){
    re = mag*TMath::Cos(phase);
    im = mag*TMath::Sin(phase);
  }
  
  inline double getAbs() const{
    return TMath::Sqrt(re*re+im*im);
  }
  
  inline double getAbsSq() const{
    return (re*re+im*im);
  }
  
  inline double getPhase() const{
    return TMath::ATan2(im,re);
  }

  operator std::complex<double>() 
  {
    return std::complex<double>(re,im); 

  }

  inline FFTWComplex conj() const
  {
    FFTWComplex copy(*this); 
    copy.im = -copy.im; 
    return copy; 
  }

};

// overload string stream operator... just for fun
std::ostream& operator<<(std::ostream& os, const FFTWComplex& val);


#endif // FFTWCOMPLEX_H
