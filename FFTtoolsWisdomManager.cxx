/* -*- C++ -*-.***************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             What I really want is a FFTtools constructor or destructor.
	     But since I made it a static class I need another class which has a constructor.
	     Then I can make a static instance of this class inside the FFTtools class.
*************************************************************************************************************** */

#include "FFTtoolsWisdomManager.h"

/* Define static members */
/* https://stackoverflow.com/questions/18433752/c-access-private-static-member-from-public-static-method */
std::map<int, fftw_plan> FFTtoolsWisdomManager::fRealToComplex;
std::map<int, fftw_plan> FFTtoolsWisdomManager::fComplexToReal;
std::map<int, double*> FFTtoolsWisdomManager::fReals;
std::map<int, std::complex<double>*> FFTtoolsWisdomManager::fComplex;

FFTtoolsWisdomManager::FFTtoolsWisdomManager(){
  unableToRead=0;
  unableToWrite=0;
  importWisdom();
}

FFTtoolsWisdomManager::~FFTtoolsWisdomManager(){
  exportWisdom();
}

void FFTtoolsWisdomManager::importWisdom(){

  const char* fftwWisdomEnv = "FFTW_WISDOM_DIR";
  const char* anitaUtilEnv = "ANITA_UTIL_INSTALL_DIR";

  const int numWisdomEnvs = 2;
  const char* envsToTry[numWisdomEnvs] = {fftwWisdomEnv, anitaUtilEnv};
  FILE* wisdomFile = NULL;
  wisdomDir = NULL;
  wisdomEnv = NULL;

  for(int envInd=0; envInd<numWisdomEnvs; envInd++){

    const char* thisWisdomDir = getenv(envsToTry[envInd]);

    // fprintf(stderr, "%s\n", envsToTry[envInd]);
    // fprintf(stderr, "%s\n", thisWisdomDir);
  
    if(thisWisdomDir==NULL){// no environment variable
      // std::cerr << "Warning in FFTtoolsWisdomManager! Please set environment variable " 
      // 		<<  anitaUtilEnv << " so I know where to save fftw3 wisdom." << std::endl;
    }
    else{
      wisdomDir = thisWisdomDir;
      wisdomEnv = (char*) envsToTry[envInd];
      char wisdomFileName[FILENAME_MAX];
      sprintf(wisdomFileName, "%s/fftw.wisdom", wisdomDir);
      wisdomFile = fopen(wisdomFileName, "r");
      if(wisdomFile!=NULL){
	fftw_import_wisdom_from_file(wisdomFile);
	fclose(wisdomFile);
      }

      // We have a set environment variable so stop here.
      break;
    }
  }
  if(wisdomFile==NULL){//then no file
    unableToRead = 1;
    unableToWrite = exportWisdom(); // Try to export the wisdom to see if we can write to the file...
    if(unableToRead!=0 && unableToWrite!=0){
      std::cerr << "******************************************************************************************\n";
      std::cerr << "Warning in " << __FILE__ << "! Unable to read or write to a fftw.wisdom file!" << std::endl;
      if(unableToRead==1){ // then got wisdom dir but couldn't write
	std::cerr << "Tried to write to " << wisdomEnv << " = " << wisdomDir 
		  << ", but either it doesn't exist or this user doesn't have permission." << std::endl;
	if(strcmp(wisdomEnv, anitaUtilEnv)==0){
	  std::cerr << "If you don't want to completely change your setup for " << anitaUtilEnv << ", try setting the environment variable " << fftwWisdomEnv << " to point to somewhere this user has permission to read and write." << std::endl;
	}
	else{
	  std::cerr << "Try setting the environment variable " << fftwWisdomEnv << " to point to somewhere this user has permission to read and write." << std::endl;
	}
      }
      else{
	std::cerr << "Couldn't find either environment variables " << fftwWisdomEnv 
		  << " or " << anitaUtilEnv << std::endl;
      }
      std::cerr << "******************************************************************************************\n";
    }
  }
}

int FFTtoolsWisdomManager::exportWisdom(){
  if(wisdomDir != NULL){
    char wisdomFileName[FILENAME_MAX];
    sprintf(wisdomFileName, "%s/fftw.wisdom", wisdomDir);

    FILE* wisdomFile = fopen(wisdomFileName, "w");
    if(wisdomFile!=NULL){
      fftw_export_wisdom_to_file(wisdomFile);
      fclose(wisdomFile);
      return 0; // success
    }
    else{
      return 1; // got wisdomDir but can't write
    }
  }
  else{
    return 2; // don't have wisdomDir
  }
}


bool FFTtoolsWisdomManager::makeNewPlanIfNeeded(int len){
  /* 
     Function which checks whether we've encountered a request to do an FFT of this length before.
     If we haven't then we need a new plan!
  */
  std::map<int,fftw_plan>::iterator it = fRealToComplex.find(len);
  if(it==fRealToComplex.end()){
    fReals[len] = (double*) fftw_malloc(sizeof(double)*len);
    fComplex[len] = (std::complex<double>*) fftw_malloc(sizeof(fftw_complex)*len);
    fRealToComplex[len] = fftw_plan_dft_r2c_1d(len,fReals[len],(fftw_complex*)fComplex[len],FFTW_MEASURE);
    fComplexToReal[len] = fftw_plan_dft_c2r_1d(len,(fftw_complex*)fComplex[len],fReals[len],FFTW_MEASURE);
    return true;
  }
  else{
    return false;
  }
}


