#include "FFTWComplex.h"
#include <iostream>


FFTWComplex::FFTWComplex() 
  :re(0),im(0)
{
  //Default Constructor
  
}


FFTWComplex::~FFTWComplex() 
{
  //Default destructor
}


FFTWComplex::FFTWComplex(double real, double imag) 
  :re(real),im(imag)
{
  //Assignment Constructor
}



FFTWComplex& FFTWComplex::operator=(const FFTWComplex &rhs) 
{
  re=rhs.re;
  im=rhs.im;
  return *this;
}



FFTWComplex FFTWComplex::operator*(const FFTWComplex &rhs)
{
  return FFTWComplex(*this) *= rhs;
 
}

FFTWComplex& FFTWComplex::operator*=(const FFTWComplex &rhs)
{
  Double_t mag = this->getAbs()*rhs.getAbs();
  Double_t phase = this->getPhase()+rhs.getPhase();
  
  this->setMagPhase(mag,phase);
  return *this;
}


FFTWComplex FFTWComplex::operator+(const FFTWComplex &rhs)
{
  return FFTWComplex(*this) += rhs;
}

FFTWComplex& FFTWComplex::operator+=(const FFTWComplex &rhs)
{
  // FFTWComplex other(rhs);
  this->re+=rhs.re;
  this->im+=rhs.im;
  return *this;
}


FFTWComplex FFTWComplex::operator-(const FFTWComplex &rhs)
{
  return FFTWComplex(*this) -= rhs;
}

FFTWComplex& FFTWComplex::operator-=(const FFTWComplex &rhs)
{
  
  FFTWComplex other(rhs);
  this->re-=other.re;
  this->im-=other.im;
  return *this;
}


FFTWComplex FFTWComplex::operator/(const FFTWComplex &rhs)
{
  return FFTWComplex(*this) /= rhs;
}

FFTWComplex& FFTWComplex::operator/=(const FFTWComplex &rhs)
{
  
  // FFTWComplex other(rhs);
  Double_t mag=this->getAbs()/rhs.getAbs();
  Double_t phase=this->getPhase()-rhs.getPhase();
  this->setMagPhase(mag,phase);
  return *this;
}


void FFTWComplex::setMagPhase(Double_t mag, Double_t phase) 
{
  this->re=mag*TMath::Cos(phase);
  this->im=mag*TMath::Sin(phase);
}


std::ostream& operator<<(std::ostream& os, const FFTWComplex& val){  
  return (os << "(" << val.re << "," << val.im << "i)");
}
