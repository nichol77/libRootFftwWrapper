
#include "FFTWComplex.h"

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



FFTWComplex& FFTWComplex::operator*(const FFTWComplex &rhs)
{
  return FFTWComplex(*this) *= rhs;
}

FFTWComplex& FFTWComplex::operator*=(const FFTWComplex &rhs)
{
  FFTWComplex other(rhs);
  Double_t mag=this->getAbs()*other.getAbs();
  Double_t phase=this->getPhase()+other.getPhase();
  this->setMagPhase(mag,phase);
  return *this;
}

void FFTWComplex::setMagPhase(Double_t mag, Double_t phase) 
{
  this->re=mag*TMath::Cos(phase);
  this->im=mag*TMath::Sin(phase);
}
