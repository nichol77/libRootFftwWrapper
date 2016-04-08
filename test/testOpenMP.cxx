#include "FFTtools.h" 
#include <omp.h> 
#include <stdlib.h>
#include "TRandom3.h"
#include "TStopwatch.h" 
#include <stdio.h>


int sizes[] = { 16,128,512,1024,4096,1025,4093,123,60,41}; 
const int nsizes = sizeof(sizes)/sizeof(*sizes);


static __thread TRandom3 * rng = 0; 

int main(int nargs, char ** args) 
{
  int howmany = nargs > 1 ? atoi(args[1]) : 100000; 

  FFTtools::loadWisdom("wis.dom"); 

  

  #pragma omp parallel for
  for (int i = 0; i < howmany; i++) 
  {
    TStopwatch sw; 
    if (rng == 0) rng = new TRandom3(0); 

    int size = sizes[rng->Integer(nsizes)]; 

    double * x = new double[size]; 
    for (int j = 0; j < size; j++) 
    {
      x[j] = rng->Gaus(); 
    }

    FFTWComplex * X = FFTtools::doFFT(size,x); 
    double * x2 = FFTtools::doInvFFT(size,X); 

    double resid = 0; 

    for (int j = 0; j < size; j++)
    {
      resid +=  (x[j] - x2[j]) * (x[j]-x2[j]); 
    }
    
    printf("i: %d thread: %d size: %d resid: %f\n", i, omp_get_thread_num(), size, resid); 
    sw.Print("u"); 

    printf("\n"); 

    delete  [] x; 
    delete  [] X; 
    delete  [] x2; 
  }


  FFTtools::saveWisdom("wis.dom"); 
  return 0; 
}


