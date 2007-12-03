#ifndef FFTWCOMPLEX_H
#define FFTWCOMPLEX_H


class FFTWComplex{
public:
    FFTWComplex();
    ~FFTWComplex();
    FFTWComplex(double real, double imag);
    double re;
    double im;
};

#endif // FFTWCOMPLEX_H
