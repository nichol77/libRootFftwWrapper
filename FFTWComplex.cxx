
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
