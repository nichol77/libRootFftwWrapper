#include "FFTWComplex.h"
#include <iostream>


// FFTWComplex::FFTWComplex() 

// {
//   //Default Constructor
  
// }


// FFTWComplex::~FFTWComplex() 
// {
//   //Default destructor
// }


// FFTWComplex::FFTWComplex(double real, double imag) 
//   :re(real),im(imag)
// {
//   //Assignment Constructor
// }



// FFTWComplex& FFTWComplex::operator=(const FFTWComplex &rhs) 
// {
//   re=rhs.re;
//   im=rhs.im;
//   return *this;
// }



// FFTWComplex FFTWComplex::operator*(const FFTWComplex &rhs)
// {
//   return FFTWComplex(*this) *= rhs;
 
// }

// FFTWComplex& FFTWComplex::operator*=(const FFTWComplex &rhs)
// {
//   // Double_t mag = this->getAbs()*rhs.getAbs();
//   // Double_t phase = this->getPhase()+rhs.getPhase();

//   Double_t newRe = re*rhs.re - im*rhs.im;
//   Double_t newIm = im*rhs.re + re*rhs.im;  
//   re = newRe;
//   im = newIm;

//   return *this;
// }


// FFTWComplex FFTWComplex::operator+(const FFTWComplex &rhs)
// {
//   return FFTWComplex(*this) += rhs;
// }

// FFTWComplex& FFTWComplex::operator+=(const FFTWComplex &rhs)
// {
//   // FFTWComplex other(rhs);
//   this->re+=rhs.re;
//   this->im+=rhs.im;
//   return *this;
// }


// FFTWComplex FFTWComplex::operator-(const FFTWComplex &rhs)
// {
//   return FFTWComplex(*this) -= rhs;
// }

// FFTWComplex& FFTWComplex::operator-=(const FFTWComplex &rhs)
// {
  
//   FFTWComplex other(rhs);
//   this->re-=other.re;
//   this->im-=other.im;
//   return *this;
// }


// FFTWComplex FFTWComplex::operator/(const FFTWComplex &rhs)
// {
//   return FFTWComplex(*this) /= rhs;
// }

// FFTWComplex& FFTWComplex::operator/=(const FFTWComplex &rhs)
// {
  
//   // FFTWComplex other(rhs);
//   Double_t norm = rhs.re*rhs.re + rhs.im*rhs.im;
//   Double_t newRe = re*rhs.re - im*rhs.im;
//   Double_t newIm = im*rhs.re + re*rhs.im;
//   re = newRe/norm;
//   im = newIm/norm;  
  
//   return *this;
// }


std::ostream& operator<<(std::ostream& os, const FFTWComplex& val){  
  return (os << "(" << val.re << "," << val.im << "i)");
}
