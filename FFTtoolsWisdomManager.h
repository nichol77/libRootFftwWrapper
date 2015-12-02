/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             What I really want is a FFTtools constructor or destructor.
	     But since I made it a static class I need another class which has a constructor.
	     Then I can make a static instance of that class inside the FFTtools class.
*************************************************************************************************************** */


#ifndef FFTTOOLS_WISDOM_MANAGER_H
#define FFTTOOLS_WISDOM_MANAGER_H

#include <iostream>
#include <complex>
#include <fftw3.h>
#include <map>
#include <stdlib.h>
#include <cstring>

class FFTtoolsWisdomManager{
public:
  FFTtoolsWisdomManager();
  ~FFTtoolsWisdomManager();
  void importWisdom();
  int exportWisdom();

  int unableToRead;
  int unableToWrite;
  char* wisdomEnv;

 /* Takes care of checking whether a plan exists or not */
  static bool makeNewPlanIfNeeded(int len);

  static std::map<int, fftw_plan> fRealToComplex;
  static std::map<int, fftw_plan> fComplexToReal;
  static std::map<int, double*> fReals;
  static std::map<int, std::complex<double>*> fComplex;

  
private:
  const char* wisdomDir;

  /*
     std::maps which hold all the fftw goodies.
     The length is the key so we have an easy way to check if a plan already exists.
     The values are the input/ouput arrays and the plans, real-to-complex and complex-to-real.
  */

};


#endif //FFTTOOLS_WISDOM_MANAGER_H
