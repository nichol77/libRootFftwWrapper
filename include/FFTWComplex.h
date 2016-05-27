#ifndef FFTWCOMPLEX_H
#define FFTWCOMPLEX_H
#include "TMath.h"
#include <ostream>

//!  This is a wrapper class for a complex number
/*!
  And that's it.
*/
class FFTWComplex{
public:
  FFTWComplex(); ///<Default constructor
  ~FFTWComplex(); ///<Destructor
  FFTWComplex(double real, double imag); ///<Assignment constructor
  FFTWComplex& operator=(const FFTWComplex &rhs);
  FFTWComplex operator*(const FFTWComplex &rhs);
  FFTWComplex& operator*=(const FFTWComplex &rhs);
  FFTWComplex operator+(const FFTWComplex &rhs);
  FFTWComplex& operator+=(const FFTWComplex &rhs);
  FFTWComplex operator-(const FFTWComplex &rhs);
  FFTWComplex& operator-=(const FFTWComplex &rhs);
  FFTWComplex operator/(const FFTWComplex &rhs);
  FFTWComplex& operator/=(const FFTWComplex &rhs);
  void setMagPhase(double mag, double phase);
  double re; ///<The real part
  double im; ///<The imaginary part
  double getAbs() const{
    return TMath::Sqrt(re*re+im*im);
  }
  double getAbsSq() const{
    return (re*re+im*im);
  }
  double getPhase() const{
    return TMath::ATan2(im,re);
  }

  
};

// overload string stream operator... just for fun
std::ostream& operator<<(std::ostream& os, const FFTWComplex& val);


#endif // FFTWCOMPLEX_H
