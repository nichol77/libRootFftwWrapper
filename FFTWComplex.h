#ifndef FFTWCOMPLEX_H
#define FFTWCOMPLEX_H


//!  This is a wrapper class for a complex number
/*!
  And that's it.
*/
class FFTWComplex{
public:
   FFTWComplex(); ///<Default constructor
   ~FFTWComplex(); ///<Destructor
   FFTWComplex(double real, double imag); ///<Assignment constructor
   double re; ///<The real part
   double im; ///<The imaginary part
};

#endif // FFTWCOMPLEX_H
