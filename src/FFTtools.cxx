#include "FFTtools.h"

#include <fftw3.h>
#include "FFTWindow.h"
#include "TRandom.h" 
#include <assert.h>
#include "TF1.h" 
#include <algorithm>

#ifndef __APPLE__
#define SINCOS sincos 
#else
#define SINCOS __sincos
#endif




using namespace std;


Double_t FFTtools::bartlettWindow(Int_t j, Int_t n)
{ 
return 1. - TMath::Abs(Double_t(2*j -n)/n);
}

Double_t FFTtools::welchWindow(Int_t j, Int_t n)
{ 
  return 1. -TMath::Power(Double_t(2*j -n)/n,2);
}

TGraph *FFTtools::getInterpolatedGraph(TGraph *grIn, Double_t deltaT)
{
  //Will use the ROOT::Math::Interpolator function to do this.
  std::vector<double> tVec;
  std::vector<double> vVec;
   
  Int_t numIn=grIn->GetN();
  Double_t tIn,vIn;

  Double_t startTime=0;
  Double_t lastTime=0;
  for (int samp=0;samp<numIn;samp++) {
    grIn->GetPoint(samp,tIn,vIn);
    tVec.push_back(tIn);
    vVec.push_back(vIn);
    //std::cout << "samp " << samp << " t " << tIn << " v " << vIn << " this-last " << tIn-tVec[tVec.size()-2] << std::endl;
    if(samp==0)
      startTime=tIn;
    lastTime=tIn;
   }
   if(tVec.size()<1) {
     std::cout << "Insufficent points for interpolation\n";
     return NULL;
   }

   //Bastards
   ROOT::Math::Interpolator chanInterp(tVec,vVec,ROOT::Math::Interpolation::kAKIMA);
   
   Int_t roughPoints=Int_t((lastTime-startTime)/deltaT);
   

   Double_t *newTimes = new Double_t[roughPoints+100]; //Will change this at some point, but for now
   Double_t *newVolts = new Double_t[roughPoints+100]; //Will change this at some point, but for now
   Int_t numPoints=0;
   for(Double_t time=startTime;time<=lastTime;time+=deltaT) {
      newTimes[numPoints]=time;
      newVolts[numPoints]=chanInterp.Eval(time);
      //      std::cout << numPoints << "\t" << newTimes[numPoints]
      //      		<< "\t" << newVolts[numPoints] << std::endl;
	       
      numPoints++;
   }

   TGraph *grInt = new TGraph(numPoints,newTimes,newVolts);
   delete [] newTimes;
   delete [] newVolts;
   return grInt;

}



TGraph *FFTtools::getInterpolatedGraphFreqDom(TGraph *grIn, Double_t deltaT)
{
   
  Int_t numIn=grIn->GetN();
  Double_t *tIn=grIn->GetX();
  Double_t *vIn=grIn->GetY();
  Double_t oldDt=tIn[1]-tIn[0];
  if(deltaT>oldDt)
    return getInterpolatedGraph(grIn,deltaT);
  
  FFTWComplex *theFFt = doFFT(numIn,vIn);
  Int_t fftLength=(numIn/2)+1;
  Int_t newFFTLength=(oldDt/deltaT)*fftLength;
  FFTWComplex *thePaddedFft = new FFTWComplex[newFFTLength];
  Int_t numPoints=(newFFTLength-1)*2;
  //  std::cerr << numIn << "\t" << fftLength << "\t" << newFFTLength << "\n";
  Double_t scaleFactor=Double_t(numPoints)/Double_t(numIn);
  for(int i=0;i<newFFTLength;i++) {
    if(i<fftLength) {
      thePaddedFft[i]=theFFt[i];
      thePaddedFft[i]*=FFTWComplex(scaleFactor,0);
    }
    else {
      thePaddedFft[i].re=0;
      thePaddedFft[i].im=0;
    }
  }
  
  Double_t *newTimes = new Double_t[numPoints]; //Will change this at some point, but for now
  Double_t *newVolts = doInvFFT(numPoints,thePaddedFft);
  for(Int_t i=0;i<numPoints;i++) {
    newTimes[i]=tIn[0]+deltaT*(i-1);
    //    std::cout << i<< "\t" << newTimes[i]
    //	      << "\t" << newVolts[i] << std::endl;	       
  }

   TGraph *grInt = new TGraph(numPoints,newTimes,newVolts);
   delete [] newTimes;
   delete [] newVolts;
   delete [] theFFt;
   delete [] thePaddedFft;
   return grInt;

}




 /***** WARNING: You are about to enter a zone that is quite complicated!!!
 *
 *   The basic idea is this: FFTW3 figures out a super fast way to do a
 *   particular FFT, based on its size (and memory alignment, but in practice
 *   we always want things aligned optimally). This operation (planning) takes
 *   a while, and while FFTW3 will "remember" plans it already figured out, its
 *   hash function for storing them is needlesly slow (it does an md5sum of
 *   some properties of the transform, which is slow, partially, at least on
 *   Linux, because glibc sincos is so slow... but I digress)  
 *
 *   So the solution here is to use a std::map to store the FFtW3 goodies. This
 *   will use a tree structure, ensuring O(logN) lookup. 
 *
 *   But because arrays you pass are not necessarily going to be aligned
 *   properly for SIMD, we allocate memory here and will copy data in and out
 *   of the arrays. This might seem like it's slower, and it might be in some
 *   cases, but this is what we do. 
 *
 *   If FFTTOOLS_ALLOCATE_CONTIGUOUS is true, the input and output arrays will
 *   be (almost) contiguous, which theoretically might improve cache locality,
 *   but I haven't done any rigorous profiling. 
 *
 *   Finally, you might notice a bunch of threading-related crap littering the
 *   code.  This is because the FFTW3 planning and the static arrays are not
 *   thread-safe.  There are two threading models implemented, both
 *   experimental until such time as someone decides they're not (and
 *   therefore, neither enabled by default). 
 *
 *   The current strategy for thread safety is to lock the planning stage and
 *   then allocate different arrays for each thread. The FFT's are then
 *   executed using the new-array variants of  fft_execute. It would be nice to
 *   just use TSD for this and just be able to mark the maps as thread_local,
 *   but that is not well implemented enough in compilers yet (afaik?) and I'm
 *   not sure if it works with OpenMP or not (but maybe it does?). 
 *
 *
 *   FFTTOOLS_THREAD_SAFE adds the necessary machinery to use ROOT's built-in
 *   threading facilities (i.e. pthreads on Unix, and the Windows thread
 *   thingie on on Windows).
 *
 *   FFTTOOLS_USE_OMP uses the OpenMP macros, which are somewhat easier to use. 
 *
 *   Trying to use both will crash and burn. don't
 *
 *
 *   There is one fairly big drawback to the way things are done; the per-thread arrys
 *   are never deleted. As long as you reuse threads instead of spawning countless
 *   numbers of them, this should be ok. In the future, maybe I'll properly use 
 *   thread_local and shit to make it work, but I don't feel like it's right to make
 *   people use a C++11 compiler until at least when people stop using Fortran77 
 *
 *******************************************************************************/


/** This caches both forward and backwards plans **/
static std::map<int, std::pair<fftw_plan, fftw_plan> > cached_plans; //this caches both plans for a given length



/* error if you try to compile with both */ 
#ifdef FFTTOOLS_THREAD_SAFE
#ifdef FFTTOOLS_USE_OMP
 Sorry, you cannot compile using both FFTTOOLS_THREAD_SAFE and FFTTOOLS_USE_OMP. 
#endif 
#endif


/** ROOT thread stuff **/
#ifdef FFTTOOLS_THREAD_SAFE
#define USE_PER_THREAD_MEMORY
#include "TMutex.h" 
#include "TThread.h" 
static TMutex plan_mutex; 
#define GET_TID() TThread::SelfId()
#endif 


/** OpenMP stuff **/ 
#ifdef FFTTOOLS_USE_OMP
#define USE_PER_THREAD_MEMORY
#include "omp.h"
#define GET_TID() omp_get_thread_num() 
#endif


#ifdef  USE_PER_THREAD_MEMORY
/*define memory for each thread, this will use the thread-id as an index 
 *
 *  For OMP, the thread numbers increase sequentially, so I use a vector to 
 *  increase lookup speed. For pthreads, the thread number is an opaque handle,
 *  so I use a map instead 
 * */

#ifdef FFTTOOLS_USE_OMP
static  std::map<int, std::vector<double *> > cached_xs;
static  std::map<int, std::vector<fftw_complex *> > cached_Xs; 
#else
static  std::map<int, std::map<long,double *> > cached_xs;
static  std::map<int, std::map<long,fftw_complex *> > cached_Xs; 
#endif

#else
static std::map<int, double*> cached_xs;
static std::map<int, fftw_complex*> cached_Xs;
#endif

static fftw_plan getPlan(int len, bool forward = true)
{
  /* 
     Function which checks whether we've encountered a request to do an FFT of this length before.
     If we haven't then we need a new plan!
  */

#ifdef FFTOOLS_THREAD_SAFE
  plan_mutex.Lock(); //root lock 
#endif 



  std::map<int,std::pair<fftw_plan, fftw_plan> >::iterator it = cached_plans.find(len);


  bool must_plan = it == cached_plans.end(); 
  bool must_allocate = must_plan; 

#ifdef USE_PER_THREAD_MEMORY
  unsigned long tid = GET_TID(); 

#ifdef FFTTOOLS_USE_OMP
  must_allocate = must_allocate || cached_xs[len].size() <= tid || cached_xs[len][tid] == 0 ; 
#else
  must_allocate = must_allocate || cached_xs[len].count(tid) == 0 || cached_xs[len][tid] == 0 ; 
#endif

#endif 

  double * mem_x = 0; 
  fftw_complex * mem_X = 0; 

  if(must_allocate)
  {
#ifdef FFTTOOLS_ALLOCATE_CONTIGUOUS
    /* Allocate contiguous memory for both the real and complex arrays 
     *
     *  For proper alignment, we ought to pad sizeof(double) * len to whatever the largest SIMD alignment is. 
     *  Currently that's sizeof(double) * 4, but in the future it might be sizeof(double) * 8; 
     **/ 
 
    void * mem = fftw_malloc(sizeof(double) * (len + len % 8) +  (1 + len/2) * sizeof(fftw_complex)); 
    mem_x = (double*) mem;  //pointer to real array 
    mem_X = (fftw_complex*) (mem_x + len + (len % 8));  //pointer to complex array

//    printf ("%lu %lu\n", ((size_t)mem_x) % 32, ((size_t)mem_X) % 32);  // check that aligned
#else
    mem_x= fftw_alloc_real(len); 
    mem_X= fftw_alloc_complex(len/2+1); 
#endif

#ifdef USE_PER_THREAD_MEMORY


#ifdef FFTTOOLS_USE_OMP
    if (must_plan)
    {
      cached_xs[len] = std::vector<double*>(tid+1,0); 
      cached_Xs[len] = std::vector<fftw_complex*>(tid+1,0); 
    }

    if (cached_xs[len].size() <= tid)
    {
      cached_xs[len].resize(tid+1,0); 
      cached_Xs[len].resize(tid+1,0); 
    }

#endif

    cached_xs[len][tid] = mem_x; 
    cached_Xs[len][tid] = mem_X; 

//    printf("FFTtools:: Allocated memory for thread %lu, length %d\n", tid, len); 
#else
    cached_xs[len] = mem_x;  //insert into maps 
    cached_Xs[len] = mem_X; 
#endif

  }
   

  if (must_plan)
  {

    //create plans
    std::pair<fftw_plan,fftw_plan> plans; 
#ifdef FFTW_USE_PATIENT
    plans.first = fftw_plan_dft_r2c_1d(len,mem_x, mem_X, FFTW_PATIENT | FFTW_PRESERVE_INPUT);
    plans.second = fftw_plan_dft_c2r_1d(len,mem_X, mem_x, FFTW_PATIENT); 
#else
    plans.first = fftw_plan_dft_r2c_1d(len,mem_x, mem_X, FFTW_MEASURE | FFTW_PRESERVE_INPUT);
    plans.second = fftw_plan_dft_c2r_1d(len,mem_X, mem_x, FFTW_MEASURE); 
#endif
    cached_plans[len] = plans; 
  }

#ifdef FFTOOLS_THREAD_SAFE
    plan_mutex.UnLock(); 
#endif 

 if (!must_plan) 
 {
   return forward ? (*it).second.first : (*it).second.second; 
 }
 else
 {
   return forward ? cached_plans[len].first : cached_plans[len].second; 
 }
  
}


void FFTtools::doFFT(int length, const double * in, FFTWComplex * out)
{

  fftw_plan plan;
#ifdef FFTTOOLS_USE_OMP
#pragma omp critical (fft_tools) 
  {
    plan = getPlan(length,true); 
  }
#else
  plan = getPlan(length,true); 
#endif

  fftw_execute_dft_r2c(plan, (double*) in, (fftw_complex*) out);
}






void FFTtools::doInvFFTNoClobber(int length, const FFTWComplex * in, double * out)
{
  fftw_plan plan;
#ifdef FFTTOOLS_USE_OMP
#pragma omp critical (fft_tools) 
  {
    plan = getPlan(length,false); 
  }
#else
  plan = getPlan(length,false); 
#endif

  fftw_complex new_in[length/2+1] __attribute__((aligned(32)));

#if ( __clang__major__ >3 || (__clang_major__==3 && __clang_minor__ >=6)  || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7) || (__GNUC__ > 4))
  FFTWComplex * ain = (FFTWComplex*) __builtin_assume_aligned(in,32);
#else
  FFTWComplex * ain = (FFTWComplex*) in;
#endif
  memcpy(new_in, ain, (length/2+1) * sizeof(FFTWComplex)); 
  fftw_execute_dft_c2r(plan,new_in, out);

#if ( __clang__major__ >3 || (__clang_major__==3 && __clang_minor__ >=6)  || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7) || (__GNUC__ > 4))
  double * aout = (double*) __builtin_assume_aligned(out,32);   
#else
  double * aout = (double*) out;
#endif  


  for (int i = 0; i < length; i++) 
  {
    aout[i] /= length; 
  }
}


void FFTtools::doInvFFTClobber(int length, FFTWComplex * in, double * out)
{
  fftw_plan plan;
#ifdef FFTTOOLS_USE_OMP
#pragma omp critical (fft_tools) 
  {
    plan = getPlan(length,false); 
  }
#else
  plan = getPlan(length,false); 
#endif


#if ( __clang__major__ >3 || (__clang_major__==3 && __clang_minor__ >=6)  || (__GNUC__ == 4 && __GNUC_MINOR__ >= 7) || (__GNUC__ > 4))
  fftw_execute_dft_c2r(plan, (fftw_complex*) in, out);

  double * aout = (double*) __builtin_assume_aligned(out,32);   
#else
  double * aout = (double*) out;
#endif  
  
  for (int i = 0; i < length; i++) 
  {
    aout[i] /= length; 
  }
}


FFTWComplex *FFTtools::doFFT(int length, const double *theInput) {
  //Here is what the sillyFFT program should be doing;    

  const int numFreqs= (length/2)+1;
#ifdef FFTTOOLS_USE_OMP
  fftw_plan plan; 
  fftw_complex * mem_X; 
  double * mem_x; 
  unsigned long tid = GET_TID();
#pragma omp critical (fft_tools)
  {
   plan = getPlan(length, true); 
   mem_X = cached_Xs[length][tid]; 
   mem_x = cached_xs[length][tid]; 
  }
#else
  fftw_plan plan = getPlan(length, true); 

#if FFTTOOLS_THREAD_SAFE 
  unsigned long tid = GET_TID();
  plan_mutex.Lock(); 
  fftw_complex * mem_X = cached_Xs[length][tid]; 
  double * mem_x = cached_xs[length][tid]; 
  plan_mutex.UnLock(); 
#endif

#endif

  FFTWComplex *myOutput = new FFTWComplex [numFreqs];

#ifdef USE_PER_THREAD_MEMORY
//  printf("%x\n", (size_t) mem_x); 
  memcpy(mem_x, theInput, sizeof(double)*length);
  fftw_execute_dft_r2c(plan, mem_x, mem_X);
  memcpy(myOutput, mem_X, sizeof(fftw_complex)*numFreqs);
#else
  memcpy(cached_xs[length], theInput, sizeof(double)*length);
  fftw_execute(plan);
  memcpy(myOutput, cached_Xs[length], sizeof(fftw_complex)*numFreqs);
#endif

  return myOutput;
}



double *FFTtools::doInvFFT(int length, const FFTWComplex *theInput) {
  // This is what sillyFFT should be doing
  //    //Takes account of normailisation 
  // Although note that fftw_plan_dft_c2r_1d assumes that the frequency array is only the positive half, so it gets scaled by sqrt(2) to account for symmetry


#ifdef FFTTOOLS_USE_OMP
  fftw_complex * mem_X; 
  double * mem_x; 
  fftw_plan plan; 
  unsigned long tid = GET_TID(); 
#pragma omp critical (fft_tools)
  {
   plan = getPlan(length, false); 
   mem_X = cached_Xs[length][tid]; 
   mem_x = cached_xs[length][tid]; 
  }
#else 
  fftw_plan  plan = getPlan(length, false); 

#if FFTTOOLS_THREAD_SAFE
  plan_mutex.Lock(); 
  unsigned long tid = GET_TID();
  fftw_complex * mem_X = cached_Xs[length][tid]; 
  double *mem_x = cached_xs[length][tid]; 
  plan_mutex.UnLock(); 
#endif

#endif
 
  double *theOutput = new double [length];

#ifdef USE_PER_THREAD_MEMORY
  memcpy(mem_X, theInput, sizeof(fftw_complex) * (length/2 + 1)); 
  fftw_execute_dft_c2r(plan, mem_X, mem_x);
#else
  memcpy (cached_Xs[length], theInput, sizeof(fftw_complex) * (length/2 +1)); 
  fftw_execute(plan);
  double * mem_x = cached_xs[length]; 
#endif

  /* Normalization needed on the inverse transform */
  for(int i=0; i<length; i++){
    mem_x[i]/=length;
  }

  // Copy output array
  memcpy(theOutput, mem_x, sizeof(double)*length);
  return theOutput;
}



/* done with complicated code */ 


Double_t *FFTtools::combineValuesUsingFFTs(Int_t numArrays, Double_t **thePtrPtr, Int_t eachLength) {
    FFTWComplex **theFFTs = new FFTWComplex* [numArrays];
    for(int i=0;i<numArrays;i++) {
	theFFTs[i]=doFFT(eachLength,thePtrPtr[i]);
    }

    int fftLength=(eachLength/2)+1;
    FFTWComplex *combinedFFT = new FFTWComplex [fftLength];


    for(int i=0;i<fftLength;i++) {
	double tempAbs0=getAbs(theFFTs[0][i]);
	double tempTotAbs=tempAbs0;
	for(int arNum=1;arNum<numArrays;arNum++) {
	    tempTotAbs+=getAbs(theFFTs[arNum][i]);
	}

	combinedFFT[i].re=theFFTs[0][i].re*(tempTotAbs/(tempAbs0*double(numArrays)));
	combinedFFT[i].im=theFFTs[0][i].im*(tempTotAbs/(tempAbs0*double(numArrays)));
    }

    for(int i=0;i<numArrays;i++) {
	delete [] theFFTs[i];
    }
    delete [] theFFTs;
    double *newValues=doInvFFT(eachLength,combinedFFT);
    delete [] combinedFFT;
    return newValues;
    
}


TGraph *FFTtools::combineGraphsUsingFFTs(Int_t numGraphs, TGraph **grPtr,double *theWeights) {

    double totalWeight=0;
    if(theWeights) {
	for(int i=0;i<numGraphs;i++) {
	    totalWeight+=theWeights[i];
//	    cout << "Weight " << i << "\t" << theWeights[i] << endl;
	}
    }

    FFTWComplex **theFFTs = new FFTWComplex* [numGraphs];
    for(int i=0;i<numGraphs;i++) {
	double *oldY=grPtr[i]->GetY();
	int oldLength=grPtr[i]->GetN();
	theFFTs[i]=doFFT(oldLength,oldY);
    }

    int fftLength=((grPtr[0]->GetN())/2)+1;
    FFTWComplex *combinedFFT = new FFTWComplex [fftLength];

    for(int i=0;i<fftLength;i++) {
	if(theWeights) {
	    double tempAbs0=getAbs(theFFTs[0][i]);
	    double tempTotAbs=tempAbs0*theWeights[0];
	    for(int grNum=1;grNum<numGraphs;grNum++) {
		tempTotAbs+=getAbs(theFFTs[grNum][i])*theWeights[grNum];
	    }
	    
	    combinedFFT[i].re=theFFTs[0][i].re*(tempTotAbs/(tempAbs0*double(totalWeight)));
	    combinedFFT[i].im=theFFTs[0][i].im*(tempTotAbs/(tempAbs0*double(totalWeight)));
	}
	else {
	    double tempAbs0=getAbs(theFFTs[0][i]);
	    double tempTotAbs=tempAbs0;
	    for(int grNum=1;grNum<numGraphs;grNum++) {
		tempTotAbs+=getAbs(theFFTs[grNum][i]);
	    }
	    
	    combinedFFT[i].re=theFFTs[0][i].re*(tempTotAbs/(tempAbs0*double(numGraphs)));
	    combinedFFT[i].im=theFFTs[0][i].im*(tempTotAbs/(tempAbs0*double(numGraphs)));
	}
    }

    for(int i=0;i<numGraphs;i++) {
	delete [] theFFTs[i];
    }
    delete [] theFFTs;

    double *newX=grPtr[0]->GetX();
    int newLength=grPtr[0]->GetN();
    double *newY=doInvFFT(newLength,combinedFFT);
    TGraph *grOut = new TGraph(newLength,newX,newY);
    delete [] combinedFFT;
    return grOut;
    
}

TGraph *FFTtools::getNormalisedCorrelationGraphTimeDomain(TGraph *gr1, TGraph *gr2,  Int_t *zeroOffset, Int_t useDtRange, Double_t dtMin, Double_t dtMax) {
  //Will also assume these graphs are zero meaned... may fix this assumption
   //Now we'll extend this up to a power of 2
  int length=gr1->GetN();
  Double_t *y1=gr1->GetY();
  int length2=gr2->GetN();
  if(length2<length) length=length2;
  Double_t *y2=gr2->GetY();
  Double_t denom=gr1->GetRMS(2)*gr2->GetRMS(2);
  
  Double_t *x1=gr1->GetX();
  Double_t *x2=gr2->GetX();
    
  double deltaT=x1[1]-x1[0];
  double waveOffset=x1[0]-x2[0];

  
  int N=2*length-1;
  

  if(zeroOffset) {
    *zeroOffset=N/2;
    (*zeroOffset)+=Int_t(waveOffset/deltaT);
  }

  //Will really assume that N's are equal for now
  //  int firstRealSamp=1+(N-2*length)/2;
  //  int lastRealSamp=firstRealSamp+2*(length-1);
  int firstRealSamp=0;
  int lastRealSamp=N-1;
  int minDtIndex=0;
  int maxDtIndex=N-1;
  if(useDtRange) {
    minDtIndex=TMath::Floor((dtMin-waveOffset)/deltaT)+(N/2);
    if(minDtIndex<0) minDtIndex=0;
    maxDtIndex=TMath::Ceil((dtMax-waveOffset)/deltaT)+(N/2);
    if(maxDtIndex<0) maxDtIndex=0;
    //    std::cout << minDtIndex << "\t" << maxDtIndex << "\t" << waveOffset << "\t" << deltaT << "\t" << dtMin << "\t" << dtMax << "\t" << N << "\t" << TMath::Floor((dtMin-waveOffset)/deltaT) << "\n";
  }
    


  
  double *xVals = new double [N];
  double *corVals= new double [N];
  for(int i=minDtIndex;i<=maxDtIndex;i++) {
    //    if(i<minDtIndex || i>maxDtIndex) continue;
    int dtIndex=(i-minDtIndex);


    xVals[dtIndex]=((i-N/2)*deltaT)+waveOffset;
    corVals[dtIndex]=0;
    if(i>=firstRealSamp && i<=lastRealSamp) {
      //if (i-firstRealSamp)==0 only one entry in correlation
      Int_t firstIndex=(i-firstRealSamp);
      Int_t secondIndex=length-1;
      if(firstIndex>length-1) {
	int offset=firstIndex-(length-1);
	firstIndex=length-1;
	secondIndex-=offset;
      }

      Int_t numSamples=0;
      //firstIndex+1;
      //if(secondIndex<firstIndex)
      //	numSamples=secondIndex+1;
      for(;firstIndex>=0 && secondIndex>=0;firstIndex--) {
	//	std::cout << i << "\t"  << firstIndex << "\t" << secondIndex << "\n";
	corVals[dtIndex]+=y1[firstIndex]*y2[secondIndex];
	numSamples++;
	secondIndex--;
      }
      corVals[dtIndex]/=denom*sqrt(numSamples);
    }      
  }

  TGraph *grCor = new TGraph((maxDtIndex-minDtIndex)+1,xVals,corVals);
  delete [] xVals;
  delete [] corVals;
  return grCor;
}

TGraph *FFTtools::getNormalisedCorrelationGraph(TGraph *gr1, TGraph *gr2, Int_t *zeroOffset) {
  //Will also assume these graphs are zero meaned... may fix this assumption
   //Now we'll extend this up to a power of 2
  int length=gr1->GetN();
  Double_t *y1=gr1->GetY();
  int length2=gr2->GetN();
  Double_t *y2=gr2->GetY();
  Double_t denom=gr1->GetRMS(2)*gr2->GetRMS(2);
  
  int N=int(TMath::Power(2,int(TMath::Log2(length))+2));
  if(N<length2)
    N=int(TMath::Power(2,int(TMath::Log2(length2))+2));
  
  //Will really assume that N's are equal for now
  int firstRealSamp=1+(N-2*length)/2;
  int lastRealSamp=firstRealSamp+2*(length-1);
  TGraph *grCor = getCorrelationGraph(gr1,gr2,zeroOffset);
  Double_t *corVal=grCor->GetY();
  Double_t norm1=0;
  Double_t norm2=0;
  
  for(int i=0;i<N;i++) {
    if(i>=firstRealSamp && i<=lastRealSamp) {
      if(i<=N/2) {
	norm1+=(y1[i-firstRealSamp]*y1[i-firstRealSamp]);
	norm2+=(y2[length-1-(i-firstRealSamp)]*y2[length-1-(i-firstRealSamp)]);
	int effN=1+(i-firstRealSamp);
	corVal[i]/=(sqrt(effN)*denom);
      }
      else if(i<N-1) {
	norm1-=(y1[i-1-(N/2)]*y1[i-1-(N/2)]);
	norm2-=(y2[length-(i-(N/2))]*y2[length-(i-(N/2))]);
	int effN=(1+lastRealSamp-i);
	corVal[i]/=(sqrt(effN)*denom);

      }
    }
  }

  return grCor;
}




TGraph *FFTtools::getCorrelationGraph(TGraph *gr1, TGraph *gr2, Int_t *zeroOffset) {
   //Now we'll extend this up to a power of 2
    int length=gr1->GetN();
    int length2=gr2->GetN();

    int N=int(TMath::Power(2,int(TMath::Log2(length))+2));
    if(N<length2)
       N=int(TMath::Power(2,int(TMath::Log2(length2))+2));

    //Will really assume that N's are equal for now
    int firstRealSamp=(N-length)/2;

    double *oldY1 = new double [N];
    double *oldY2 = new double [N];
    
    double x,y;
    Double_t x2,y2;
    gr1->GetPoint(1,x2,y2);
    gr1->GetPoint(0,x,y);
    double deltaT=x2-x;
    double firstX=x;

    gr2->GetPoint(0,x2,y2);
    double waveOffset=firstX-x2;
    

    //    gr1->GetPoint(N/2,x2,y2);
    //    double offset=x-x2;
    //    std::cout << length << "\t" << length2 << "\n";

    for(int i=0;i<N;i++) {
       
       if(i<firstRealSamp || i>=firstRealSamp+length)
	  y=0;
       else {
	  gr1->GetPoint(i-firstRealSamp,x,y);
       }
       oldY1[i]=y;
	  
       if(i<firstRealSamp || i>=firstRealSamp+length2)
	  y=0;
       else {
	  gr2->GetPoint(i-firstRealSamp,x,y);
       }
       oldY2[i]=y;
              
    }


    //    offset+=waveOffset;
    if(zeroOffset) {
       *zeroOffset=N/2;
       (*zeroOffset)+=Int_t(waveOffset/deltaT);
    }

    double *xVals = new double [N];
    double *yVals = new double [N];
    double *corVals=getCorrelation(N,oldY1,oldY2);
    for(int i=0;i<N;i++) {
       if(i<N/2) {
	  //Positive
	  xVals[i+(N/2)]=(i*deltaT)+waveOffset;
	  yVals[i+(N/2)]=corVals[i];
       }
       else {
	  //Negative
	  xVals[i-(N/2)]=((i-N)*deltaT)+waveOffset;
	  yVals[i-(N/2)]=corVals[i];	  
       }
    }


    TGraph *grCor = new TGraph(N,xVals,yVals);
    delete [] oldY1;
    delete [] oldY2;
    delete [] xVals;
    delete [] yVals;
    delete [] corVals;
    
    return grCor;
}



TGraph *FFTtools::getInterpolatedCorrelationGraph(TGraph *grIn1, TGraph *grIn2, Double_t deltaT)
{
  TGraph *gr1 = getInterpolatedGraph(grIn1,deltaT);
  TGraph *gr2 = getInterpolatedGraph(grIn2,deltaT);
  //  std::cout << gr1 << "\t" << gr2 << "\n"; 
  //  std::cout << gr1->GetN() << "\t" << gr2->GetN() << "\n"; 
  TGraph *grCor = getCorrelationGraph(gr1,gr2);
  //  std::cout << grCor << "\n";
  ///  std::cout << grCor->GetN() << "\n";
  delete gr1;
  delete gr2;
  return grCor;
}

double *FFTtools::getCorrelation(int length,float *oldY1, float *oldY2) 
{
    double *newY1 = new double [length];
    double *newY2 = new double [length];
    for(int i=0;i<length;i++) {
	newY1[i]=(double)oldY1[i];
	newY2[i]=(double)oldY2[i];
    }

    double *theCorr=getCorrelation(length,newY1,newY2);
    delete [] newY1;
    delete [] newY2;
    return theCorr;
}


double *FFTtools::getCorrelation(int length,double *oldY1, double *oldY2) 
{

//    cout << "Here in getCorrelation" << endl;
    FFTWComplex *theFFT1=doFFT(length,oldY1);
    FFTWComplex *theFFT2=doFFT(length,oldY2);
    

    int newLength=(length/2)+1;
//     cout << "newLength " << newLength << endl;
    FFTWComplex *tempStep = new FFTWComplex [newLength];
    int no2=length>>1;
    for(int i=0;i<newLength;i++) {
	double reFFT1=theFFT1[i].re;
	double imFFT1=theFFT1[i].im;
	double reFFT2=theFFT2[i].re;
	double imFFT2=theFFT2[i].im;

	//Real part of output 
	tempStep[i].re=(reFFT1*reFFT2+imFFT1*imFFT2)/double(no2/2);
	//Imaginary part of output 
	tempStep[i].im=(imFFT1*reFFT2-reFFT1*imFFT2)/double(no2/2);
    }
//    cout << "finished messing around" << endl;
    double *theOutput=doInvFFT(length,tempStep);
//    cout << "got inverse" << endl;
    delete [] theFFT1;
    delete [] theFFT2;
    delete [] tempStep;
    return theOutput;

}


double *FFTtools::getCorrelation(TGraph *gr1, TGraph *gr2,int firstIndex,int lastIndex) {
    int tempLength=gr1->GetN();
    if(firstIndex<0 || lastIndex>tempLength) return 0;
    
    int length=lastIndex-firstIndex;
//    double *x1 = gr1->GetX();
    double *y1 = gr1->GetY();
//    double *x2 = gr2->GetX();
    double *y2 = gr2->GetY();
//    TGraph newGr1(length,&x1[firstIndex],&y1[firstIndex]);
//    TGraph newGr2(length,&x2[firstIndex],&y2[firstIndex]);
    
//    return getCorrelation(&newGr1,&newGr2);
    return getCorrelation(length,&y1[firstIndex],&y2[firstIndex]);
}


TGraph *FFTtools::makeInverseInverseSpectrum(TGraph *grWave) {

    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
//    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT=doFFT(length,oldY);
    double *invInvSpectrum = doInvFFT(length,theFFT);

    TGraph *grInvInv = new TGraph(length,oldX,invInvSpectrum);
//     for(int i=0;i<length;i++) {
// 	cout << oldX[i] << "\t" << invInvSpectrum[i] << endl;
//     }
    return grInvInv;

}

double* FFTtools::getHilbertTransform(int length, const double * y)
{

  FFTWComplex *theFFT=doFFT(length,y);  
  int newLength=(length/2)+1;
  for(int i=0;i<newLength;i++) {
    double tempIm=theFFT[i].im;
    theFFT[i].im=theFFT[i].re;
    theFFT[i].re=-1*tempIm;
  }
  double *hilbert = doInvFFT(length,theFFT);
  
  delete [] theFFT;
  return hilbert; 
}
TGraph *FFTtools::getHilbertTransform(TGraph *grWave)
{
  double *oldY = grWave->GetY();
  double *oldX = grWave->GetX();
  int length=grWave->GetN();
  double *hilbert = getHilbertTransform(length, oldY); 
  TGraph *grHilbert = new TGraph(length,oldX,hilbert);
  delete [] hilbert;

  return grHilbert;
}


TGraph *FFTtools::getHilbertEnvelope(TGraph *grWave)
{
  //  cout << "Here" << endl;
  double *realY = grWave->GetY();
  double *x = grWave->GetX();
  TGraph *grHilbert = FFTtools::getHilbertTransform(grWave);
  double *hilY=grHilbert->GetY();
  int length=grWave->GetN();
  Double_t *envY= new Double_t[length];
  for(int i=0;i<length;i++) {

    envY[i]=TMath::Sqrt(realY[i]*realY[i] + hilY[i]*hilY[i]);
    //    cout << i << "\t" << envY[i] << endl;
  }
  TGraph *grEnvelope = new TGraph(length,x,envY);
  delete [] envY;
  delete grHilbert;
  return grEnvelope;
}

TGraph *FFTtools::getBoxCar(TGraph *grWave, Int_t halfWidth) 
{
  //Just do this the lazy way for now
  Double_t *inY = grWave->GetY();
  Double_t *inX = grWave->GetX();
  Int_t length=grWave->GetN();
  Double_t *smoothY = new Double_t[length];
  for(int i=0;i<length;i++) {
    smoothY[i]=0;
    if(i<halfWidth || length-i<=halfWidth) {
      int countVals=0;
      for(int j=i-halfWidth;j<=i+halfWidth;j++) {
	if(j>=0 && j<length) {
	  smoothY[i]+=inY[j];
	  countVals++;
	}
      }
      //      cout << i << "\t" << countVals << endl;
      smoothY[i]/=countVals;
    }
    else {
      for(int j=i-halfWidth;j<=i+halfWidth;j++) {
	smoothY[i]+=inY[j];
      }
      smoothY[i]/=1+2*halfWidth;
    }      
  }
  TGraph *grSmooth = new TGraph(length,inX,smoothY);
  delete [] smoothY;
  return grSmooth;
  
}

TGraph *FFTtools::makePowerSpectrumVoltsSecondsBartlett(TGraph *grWave) {
  Double_t *oldY=grWave->GetY();
  Double_t *t = grWave->GetX();
  Int_t numPoints = grWave->GetN();
  Double_t *newY = new Double_t [numPoints];
  for(int i=0;i<numPoints;i++) {
    newY[i]=oldY[i]*bartlettWindow(i,numPoints);
  }
  TGraph *grNew = new TGraph(numPoints,t,newY);
  delete [] newY;
  return makePowerSpectrumVoltsSeconds(grNew);
}


TGraph *FFTtools::makePowerSpectrum(TGraph *grWave) {

    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT=doFFT(length,oldY);

    int newLength=(length/2)+1;
    double *newY = new double [newLength];
    double *newX = new double [newLength];

    //    double fMax = 1/(2*deltaT);  // In GHz
    double deltaF=1/(deltaT*length);

    double tempF=0;
    for(int i=0;i<newLength;i++) {
       float power=pow(getAbs(theFFT[i]),2);
	if(i>0 && i<newLength-1) power*=2; //account for symmetry
	power/=length;
	//	if (power>0 ) power=10*TMath::Log10(power);
	//	else power=-100; //no reason
	newX[i]=tempF;
	newY[i]=power;
	tempF+=deltaF;
    }



    TGraph *grPower = new TGraph(newLength,newX,newY);
    delete [] theFFT;
    delete [] newY;
    delete [] newX;
    return grPower;

}



TGraph *FFTtools::makePowerSpectrumPeriodogram(TGraph *grWave) {

    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT=doFFT(length,oldY);
    
    int newLength=(length/2)+1;

    double *newY = new double [newLength];
    double *newX = new double [newLength];

    //    double fMax = 1/(2*deltaT);  // In Hz
    double deltaF=1/(deltaT*length);

    double tempF=0;
    for(int i=0;i<newLength;i++) {
       float power=pow(getAbs(theFFT[i]),2);
	if(i>0 && i<newLength-1) power*=2; //account for symmetry
	power/=double(length*length);
	//	if (power>0 ) power=10*TMath::Log10(power);
	//	else power=-100; //no reason
	newX[i]=tempF;
	newY[i]=power;
	tempF+=deltaF;
    }

    TGraph *grPower = new TGraph(newLength,newX,newY);
    delete [] theFFT;
    delete [] newY;
    delete [] newX;
    return grPower;

}

TGraph *FFTtools::makePowerSpectrumVoltsSeconds(TGraph *grWave) {

    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT=doFFT(length,oldY);

    int newLength=(length/2)+1;

    double *newY = new double [newLength];
    double *newX = new double [newLength];

    //    double fMax = 1/(2*deltaT);  // In Hz
    double deltaF=1/(deltaT*length); //Hz
    deltaF*=1e-6; //MHz


    double tempF=0;
    for(int i=0;i<newLength;i++) {
      float power=pow(getAbs(theFFT[i]),2);
      	if(i>0 && i<newLength-1) power*=2; //account for symmetry
	power*=deltaT/(length); //For time-integral squared amplitude
	power/=deltaF;//Just to normalise bin-widths
	//Ends up the same as dt^2, need to integrate the power (multiply by df)
	//to get a meaningful number out.

	newX[i]=tempF;
	newY[i]=power;
	tempF+=deltaF;
    }


    TGraph *grPower = new TGraph(newLength,newX,newY);
    delete [] theFFT;
    delete [] newY;
    delete [] newX;
    return grPower;

}

TGraph *FFTtools::makePowerSpectrumMilliVoltsNanoSeconds(TGraph *grWave) {

//   double *oldY = grWave->GetY(); //in millivolts
//   double *oldX = grWave->GetX(); //in nanoseconds
//   int length=grWave->GetN();

//   double milli = 0.001;
//   double nano = 0.000000001;

//   double *oldY2 = new double [length];
//   double *oldX2 = new double [length];

//   for(int i=0;i<length;i++){
//     oldY2[i]=oldY[i]*milli; //in volts
//     oldX2[i]=oldX[i]*nano; // in seconds
//   }

//   double deltaT=(oldX2[1]-oldX2[0]);
//   FFTWComplex *theFFT=doFFT(length,oldY2);

//   int newLength=(length/2)+1;

//   double *newY = new double [newLength];
//   double *newX = new double [newLength];

//   //    double fMax = 1/(2*deltaT);  // In Hz
//   double deltaF=1/(deltaT*length); //Hz
//   deltaF*=1e-6; //MHz


//   double tempF=0;
//   for(int i=0;i<newLength;i++) {
//     float power=pow(getAbs(theFFT[i]),2);
//     if(i>0 && i<newLength-1) power*=2; //account for symmetry
//     power*=deltaT/(length); //For time-integral squared amplitude
//     power/=deltaF;//Just to normalise bin-widths
//     //Ends up the same as dt^2, need to integrate the power (multiply by df)
//     //to get a meaningful number out.

//     newX[i]=tempF;
//     newY[i]=power;
//     tempF+=deltaF;
//   }

//   TGraph *grPower = new TGraph(newLength,newX,newY);
//   delete [] newY;
//   delete [] newX;
//   delete [] oldY2;
//   delete [] oldX2;
//   ///////THIS BIT COULD DELETE THE POWERSPEC?????????
//   delete [] theFFT;
//   return grPower;
  double *oldY = grWave->GetY();
  double *oldX = grWave->GetX();
  double deltaT=oldX[1]-oldX[0];
  int length=grWave->GetN();
  FFTWComplex *theFFT=doFFT(length,oldY);
  
  int newLength=(length/2)+1;
  
  double *newY = new double [newLength];
  double *newX = new double [newLength];
  
  //    double fMax = 1/(2*deltaT);  // In Hz
  double deltaF=1/(deltaT*length); //Hz
  deltaF*=1e3; //MHz
  //    std::cout << fMax << "\t" << deltaF << "\t" << deltaT << "\t" << length << std::endl;
  
  double tempF=0;
  for(int i=0;i<newLength;i++) {
    float power=pow(getAbs(theFFT[i]),2);
    if(i>0 && i<newLength-1) power*=2; //account for symmetry
    power*=deltaT/(length); //For time-integral squared amplitude
    //    power*=(1e3*1e3)/1e9;
    power/=deltaF;//Just to normalise bin-widths
    //Ends up the same as dt^2, need to integrate the power (multiply by df)
    //to get a meaningful number out.	
    
    //if (power>0 ) power=10*TMath::Log10(power);
    //else power=-1000; //no reason
    newX[i]=tempF;
    newY[i]=power;
    tempF+=deltaF;
  }

  TGraph *grPower = new TGraph(newLength,newX,newY);
  delete [] theFFT;
  delete [] newY;
  delete [] newX;
  return grPower;

}


TGraph *FFTtools::makePowerSpectrumMilliVoltsNanoSecondsdB(TGraph *grWave)
{
    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT=doFFT(length,oldY);

    int newLength=(length/2)+1;

    double *newY = new double [newLength];
    double *newX = new double [newLength];

    //    double fMax = 1/(2*deltaT);  // In Hz
    double deltaF=1/(deltaT*length); //Hz
    deltaF*=1e3; //MHz
    //    std::cout  << deltaF << "\t" << deltaT << "\t" << length << std::endl;

    double tempF=0;
    for(int i=0;i<newLength;i++) {
       float power=pow(getAbs(theFFT[i]),2);
       if(i>0 && i<newLength-1) power*=2; //Changed form 2 by RJN 29/01/10 //account for symmetry
	power*=deltaT/(length); //For time-integral squared amplitude
	//	power/=(1e3*1e3*1e9); //This bit converts from mv*mv*ns to v*v*s
	power/=deltaF;//Just to normalise bin-widths
	//Ends up the same as dt^2, need to integrate the power (multiply by df)
	//to get a meaningful number out.	
	
	if (power>0 ) power=10*TMath::Log10(power);
	else power=-1000; //no reason
	newX[i]=tempF;
	newY[i]=power;
	tempF+=deltaF;
    }


    TGraph *grPower = new TGraph(newLength,newX,newY);
    delete [] theFFT;
    delete [] newY;
    delete [] newX;
    return grPower;


}

TGraph *FFTtools::makePowerSpectrumVoltsSecondsdB(TGraph *grWave) {

    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT=doFFT(length,oldY);

    int newLength=(length/2)+1;

    double *newY = new double [newLength];
    double *newX = new double [newLength];

    //    double fMax = 1/(2*deltaT);  // In Hz
    double deltaF=1/(deltaT*length); //Hz
    deltaF*=1e-6; //MHz


    double tempF=0;
    for(int i=0;i<newLength;i++) {
       double logpower;
       double power=pow(getAbs(theFFT[i]),2);
	if(i>0 && i<newLength-1) power*=2; //account for symmetry
	power*=deltaT/(length); //For time-integral squared amplitude
	power/=deltaF;//Just to normalise bin-widths
	//Ends up the same as dt^2, need to integrate the power (multiply by df)
	//to get a meaningful number out.	
	
	if (power>0 ){
	  logpower=10*TMath::Log10(power);
	}
	else{
          logpower=-1000; //no reason
	}	

	newX[i]=tempF;
	newY[i]=logpower;
	tempF+=deltaF;
    }


    TGraph *grPower = new TGraph(newLength,newX,newY);
    delete [] theFFT;
    delete [] newY;
    delete [] newX;
    return grPower;


}

TGraph *FFTtools::makePowerSpectrumVoltsSecondsPadded(TGraph *grWave, Int_t padFactor) {

   TGraph *grPad=padWave(grWave,padFactor);
   TGraph *grPower=makePowerSpectrumVoltsSeconds(grPad);
   delete grPad;
   return grPower;
   
}


TGraph *FFTtools::makePowerSpectrumVoltsSecondsPaddeddB(TGraph *grWave, Int_t padFactor) {
   TGraph *grPad=padWave(grWave,padFactor);
   TGraph *grPower=makePowerSpectrumVoltsSecondsdB(grPad);
   delete grPad;
   return grPower;
}


TGraph *FFTtools::makeRawPowerSpectrum(TGraph *grWave) {

    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT=doFFT(length,oldY);

    int newLength=(length/2)+1;
    double *newY = new double [newLength];
    double *newX = new double [newLength];

    double deltaF=1/(deltaT*length);
    //    double fMax = 1/(2*deltaT);  // In GHz

    double tempF=0;
    for(int i=0;i<newLength;i++) {
      float power=pow(getAbs(theFFT[i]),2);
      if(i>0 && i<newLength-1) power*=2; //account for symmetry
      newX[i]=tempF;
      newY[i]=power;
      tempF+=deltaF;
    }


    TGraph *grPower = new TGraph(newLength,newX,newY);
    delete [] theFFT;
    delete [] newY;
    delete [] newX;
    return grPower;

}


double FFTtools::getAbs(FFTWComplex &theNum) {
    return sqrt(theNum.re*theNum.re+theNum.im*theNum.im);
}


Double_t FFTtools::sumPower(TGraph *gr,Int_t firstBin,Int_t lastBin)
{
  Double_t integral=0;
  Double_t freq,power;
  if(firstBin<0) firstBin=0;
  if(lastBin<0) lastBin=gr->GetN()-1;
  for(int i=firstBin;i<=lastBin;i++) {
    gr->GetPoint(i,freq,power);
    integral+=power;
  }
  return integral;
}

Double_t FFTtools::integratePower(TGraph *gr,Int_t firstBin,Int_t lastBin)
{
  Double_t integral=0;
  Double_t freq,power;
  gr->GetPoint(1,freq,power);
  Double_t df=freq;
  gr->GetPoint(0,freq,power);
  df-=freq;
  if(firstBin<0) firstBin=0;
  if(lastBin<0) lastBin=gr->GetN()-1;
  for(int i=firstBin;i<=lastBin;i++) {
    gr->GetPoint(i,freq,power);
    integral+=power*df;
  }
  return integral;
}

Double_t FFTtools::sumVoltageSquared(TGraph *gr,Int_t firstBin,Int_t lastBin)
{
  Double_t integral=0;
  Double_t time,volts;
  if(firstBin<0) firstBin=0;
  if(lastBin<0) lastBin=gr->GetN()-1;
  for(int i=firstBin;i<=lastBin;i++) {
    gr->GetPoint(i,time,volts);
    integral+=volts*volts;
  }
  return integral;
}

Double_t FFTtools::integrateVoltageSquared(TGraph *gr,Int_t firstBin,Int_t lastBin)
{
  Double_t integral=0;
  Double_t this_time=0, next_time=0, this_v=0, next_v=0, dt=0;
  
  if(firstBin<0) firstBin=0;
  if(lastBin<0) lastBin = gr->GetN()-1;
  for(int samp=firstBin; samp<lastBin; samp++){
    gr->GetPoint(samp, this_time, this_v);
    gr->GetPoint(samp+1, next_time, next_v);
    dt=next_time-this_time;
    integral+=this_v*this_v*dt;
  }
  //Now account for the last bin -- just use the last dt calculated
  integral+=next_v*next_v*dt;
  
  return integral;

}

Int_t FFTtools::getPeakBin(TGraph *gr) 
{
  Double_t x,y;
  gr->GetPoint(0,x,y);
  Double_t peakVal=y;
  Int_t peakBin=0;
  for(int i=1;i<gr->GetN();i++) {
    gr->GetPoint(i,x,y);
    if(peakVal<y) {
      peakVal=y;
      peakBin=i;
    }      
  }
  return peakBin;
}


Int_t FFTtools::getPeakBin(TGraph *gr, Int_t firstBin, Int_t lastBin) 
{
  if(firstBin<0 || lastBin<0 || firstBin>gr->GetN() || lastBin>gr->GetN() ||
     (lastBin<firstBin))
    return -1;
  Double_t x,y;
  gr->GetPoint(firstBin,x,y);
  Double_t peakVal=y;
  Int_t peakBin=0;
  for(int i=firstBin;i<lastBin;i++) {
    gr->GetPoint(i,x,y);
    if(peakVal<y) {
      peakVal=y;
      peakBin=i;
    }      
  }
  return peakBin;
}


Double_t FFTtools::getPeakVal(TGraph *gr, int *index) 
{
   Double_t x,y;
   gr->GetPoint(0,x,y);
   Double_t peakVal=y;
   Int_t peakBin=0;
   Int_t numPoints=gr->GetN();
   for(int i=1;i<numPoints;i++) {
      gr->GetPoint(i,x,y);
      if(peakVal<y) {
	 peakVal=y;
	 peakBin=i;
      }      
   }
   if(index) 
      *index=peakBin;
   return peakVal;
}

Double_t FFTtools::getPeakSqVal(TGraph *gr, int *index) 
{
   Double_t x,y;
   gr->GetPoint(0,x,y);
   Double_t peakVal=y*y;
   Int_t peakBin=0;
   Int_t numPoints=gr->GetN();
   for(int i=1;i<numPoints;i++) {
      gr->GetPoint(i,x,y);
      if(peakVal<y*y) {
	 peakVal=y*y;
	 peakBin=i;
      }      
   }
   if(index) 
      *index=peakBin;
   return peakVal;

}

void FFTtools::getPeakRmsSqVal(TGraph *gr, Double_t &peak, Double_t &rms, Int_t *index)
{
  Int_t numPoints=gr->GetN();
  Double_t *y = gr->GetY();
  Double_t *ySq = new Double_t [numPoints];
  for(int i=0;i<numPoints;i++) {
    ySq[i]=y[i]*y[i];
  }
  Int_t peakInd=TMath::LocMax(numPoints,ySq);
  peak=ySq[peakInd];
  rms=TMath::RMS(numPoints,ySq);
  if(index)
    *index=peakInd;
  delete [] ySq;
}

void FFTtools::getPeakRmsRectified(TGraph *gr, Double_t &peak, Double_t &rms, Int_t *index)
{
  TGraph *grRec = rectifyWave(gr);
  Int_t numPoints=grRec->GetN();
  Double_t *y = grRec->GetY();
  Int_t peakInd=TMath::LocMax(numPoints,y);
  peak=y[peakInd];
  rms=TMath::RMS(numPoints,y);
  if(index)
    *index=peakInd;
  delete grRec;
}

TGraph *FFTtools::subtractGraphs(TGraph *grA, TGraph *grB) 
{
  Int_t N1=grA->GetN();
  Int_t N2=grB->GetN();
  if(N1!=N2) return NULL;

  Double_t *newY = new Double_t [N1];
  Double_t *xVals=grA->GetX();
  Double_t x,yA,yB;
  for(int i=0;i<N1;i++) {
    grA->GetPoint(i,x,yA);
    grB->GetPoint(i,x,yB);
    newY[i]=yA-yB;
  }
  TGraph *grDiff = new TGraph(N1,xVals,newY);
  delete [] newY;
  return grDiff;
}


TGraph *FFTtools::translateGraph(TGraph *grWave, Double_t deltaT)
{
  Int_t N = grWave->GetN();
  Double_t *X=grWave->GetX();
  Double_t *Y=grWave->GetY();
  Double_t *newX = new Double_t[N];
  for(int i=0;i<N;i++)
    newX[i]=X[i]+deltaT;
  TGraph *grOut = new TGraph(N,newX,Y);
  delete [] newX;
  return grOut;
}

TGraph *FFTtools::divideGraphs(TGraph *grA, TGraph *grB) 
{
  Int_t N1=grA->GetN();
  Int_t N2=grB->GetN();
  if(N1!=N2) return NULL;

  Double_t *newY = new Double_t [N1];
  Double_t *xVals=grA->GetX();
  Double_t x,yA,yB;
  for(int i=0;i<N1;i++) {
    grA->GetPoint(i,x,yA);
    grB->GetPoint(i,x,yB);
    newY[i]=yA/yB;
  }
  TGraph *grRat = new TGraph(N1,xVals,newY);
  delete [] newY;
  return grRat;
}


TGraph *FFTtools::ratioSubtractOneGraphs(TGraph *grA, TGraph *grB) 
{
   Int_t N1=grA->GetN();
  Int_t N2=grB->GetN();
  //  if(N1!=N2) return NULL;
  
  Int_t newN=N1;
  if(N2<N1) {
     newN=N2;
     //return NULL;
  }
  Double_t *xVals=grA->GetX();
  Double_t *xBVals=grB->GetX();
  Double_t deltaF=xVals[1]-xVals[0];
  Double_t deltaFB=xBVals[1]-xBVals[0];
  
  if(TMath::Abs(deltaFB-deltaF)>1) return NULL;
  //  cout << newN << endl;
  Double_t *newY = new Double_t [newN];
  Double_t x,yA,yB;
  for(int i=0;i<newN;i++) {
    grA->GetPoint(i,x,yA);
    grB->GetPoint(i,x,yB);
    newY[i]=1-yA/yB;
  }
  TGraph *grRat = new TGraph(newN,xVals,newY);
  delete [] newY;
  return grRat;
}

TGraph *FFTtools::dbGraphs(TGraph *grA, TGraph *grB) 
{
  Int_t N1=grA->GetN();
  Int_t N2=grB->GetN();
  //  if(N1!=N2) return NULL;
  
  Int_t newN=N1;
  if(N2<N1) {
     newN=N2;
     //return NULL;
  }
  Double_t *xVals=grA->GetX();
  Double_t *xBVals=grB->GetX();
  Double_t deltaF=xVals[1]-xVals[0];
  Double_t deltaFB=xBVals[1]-xBVals[0];
  //  cout << N1 << "\t" << N2 << "\t" << deltaF << "\t" << deltaFB << "\n";


  if(TMath::Abs(deltaFB-deltaF)>1) return NULL;
  //  cout << newN << endl;
  Double_t *newY = new Double_t [newN];
  Double_t x,yA,yB;
  for(int i=0;i<newN;i++) {
    grA->GetPoint(i,x,yA);
    grB->GetPoint(i,x,yB);
    newY[i]=10*TMath::Log10(yA/yB);
  }
  TGraph *grRat = new TGraph(newN,xVals,newY);
  delete [] newY;
  return grRat;
}


TGraph *FFTtools::smoothFFT(TGraph *gr,Int_t factor) 
{
  Int_t N=gr->GetN();
  Int_t newN=N/factor;
  Double_t *xVals=gr->GetX();
  Double_t *yVals=gr->GetY();

  Double_t *newX = new Double_t [newN];
  Double_t *newY = new Double_t [newN];
  Double_t sumX=0;
  Double_t sumY=0;
  //  cerr << N << "\t" << factor << "\t" << newN << "\n";
  for(int i=0;i<N;i++) {
     sumX+=xVals[i];
     sumY+=yVals[i];
     if((i+1)%factor==0) {
	//	cerr << i << "\t" << sumX << "\t" << sumY << "\n";
	//New Point
	newX[(i+1)/factor-1]=sumX/factor;
	newY[(i+1)/factor-1]=sumY/factor;
	sumX=0;
	sumY=0;
     }
  }
  TGraph *grSmooth = new TGraph(newN,newX,newY);
//   delete [] newX;
//   delete [] newY;
  return grSmooth;
}

TGraph *FFTtools::padWave(TGraph *grWave, Int_t padFactor) {
   double *oldY = grWave->GetY();
   double *oldX = grWave->GetX();
   double deltaT=oldX[1]-oldX[0];
   int realLength = grWave->GetN();
   int length = grWave->GetN()*padFactor;
   double *paddedY = new double [length];
   double *paddedX = new double [length];
   int newStart=(realLength*(padFactor-1))/2;
   for(int i=0;i<length;i++) {
      int waveIndex=i-newStart;
      paddedY[i]=0;
      paddedX[i]=(waveIndex*deltaT)+oldX[0];
      if(waveIndex>=0 && waveIndex<realLength) {
	 paddedY[i]=oldY[waveIndex];
      }
   }
   TGraph *grPadded = new TGraph(length,paddedX,paddedY);
   delete [] paddedX;
   delete [] paddedY;
   return grPadded;
}


TGraph *FFTtools::padWaveToLength(TGraph *grWave, Int_t newLength) {
   double *oldY = grWave->GetY();
   double *oldX = grWave->GetX();
   double deltaT=oldX[1]-oldX[0];
   int realLength = grWave->GetN();
   if(newLength<realLength) {
     std::cerr << "Can't pad waveform of length " << realLength 
	       << " to " << newLength << " you need to crop not pad\n";
     return NULL;
   }
   int length = newLength;
   double *paddedY = new double [length];
   double *paddedX = new double [length];
   int newStart=(newLength-realLength)/2;   
   for(int i=0;i<length;i++) {
      int waveIndex=i-newStart;
      paddedY[i]=0;
      paddedX[i]=(waveIndex*deltaT)+oldX[0];
      if(waveIndex>=0 && waveIndex<realLength) {
	 paddedY[i]=oldY[waveIndex];
      }
   }
   TGraph *grPadded = new TGraph(length,paddedX,paddedY);
   delete [] paddedX;
   delete [] paddedY;
   return grPadded;
}



TGraph *FFTtools::rectifyWave(TGraph *gr, Int_t isNeg) {
  Int_t sign=1;
  if(isNeg) 
    sign=-1;

  Int_t numPoints = gr->GetN();
  Double_t *x = gr->GetX();
  Double_t *y = gr->GetY();
  Double_t *yRec = new Double_t [numPoints];
  for(int i=0;i<numPoints;i++) {
    yRec[i]=sign*TMath::Abs(y[i]);
  }
  TGraph *grRec = new TGraph(numPoints,x,yRec);
  delete [] yRec;
  return grRec;
}


TGraph *FFTtools::getSimplePowerEnvelopeGraph(TGraph *gr) {
  Double_t *ySq = new Double_t [gr->GetN()];
  Double_t *xOrig = new Double_t [gr->GetN()];
  Double_t *yEnvelope = new Double_t[gr->GetN()];
  Double_t *xEnvelope = new Double_t[gr->GetN()];
  
  Double_t x,y;
  Int_t numPoints=0;
  
  for(int i=0;i<gr->GetN();i++) {
    gr->GetPoint(i,x,y);
    ySq[i]=y*y;
    xOrig[i]=x;
    if(i==1) {
      if(ySq[0]>ySq[i]) {
        yEnvelope[numPoints]=ySq[0];
        xEnvelope[numPoints]=xOrig[0];
        numPoints++;
      }
    }
    else if(i==gr->GetN()-1 && ySq[i]>ySq[i-1]) {
      yEnvelope[numPoints]=ySq[i];
      xEnvelope[numPoints]=xOrig[i];
      numPoints++;
    }
    else if(ySq[i-1]>ySq[i-2] && ySq[i-1]>ySq[i]) {
      yEnvelope[numPoints]=ySq[i-1];
      xEnvelope[numPoints]=xOrig[i-1];
      numPoints++;
    }
  }                                  
  TGraph *grEnvelope= new TGraph(numPoints,xEnvelope,yEnvelope);
  delete [] ySq;
  delete [] xOrig;
  delete [] yEnvelope;
  delete [] xEnvelope;
  return grEnvelope;
}

TGraph *FFTtools::makePSVSBartlettPaddedOverlap(TGraph *grWave, Int_t padFactor, Int_t numFreqs)
{
  TGraph *grPad=padWave(grWave,padFactor);  
  Int_t numTot=grPad->GetN();
  Double_t *tVals = grPad->GetX();
  Double_t *vVals = grPad->GetY();
  Int_t m2=numFreqs*2;
  Double_t *inTimes = new Double_t[m2];
  Double_t *inVolts = new Double_t[m2];
  Double_t *outFreqs= new Double_t[m2];
  Double_t *outPower= new Double_t[m2];
  Int_t numSegs=(numTot-1)/numFreqs;
  Double_t delta= numSegs>1 ? Double_t(numTot-m2)/(numSegs-1.) : 0.;
  if(numTot < m2) {
    cerr << "Not enough data points for this many frequency bins" << endl;
    return 0;
  }  
  Int_t numOut=0;
  for(int i=0;i<m2;i++) {
    outPower[i]=0;
  }

  for(int k=0;k<numSegs;k++) {
    Int_t noff = (Int_t)(k*delta + 0.5);
    for(int i=0;i<m2;i++) {
      inTimes[i]=tVals[i+noff];
      inVolts[i]=vVals[i+noff]*bartlettWindow(i,m2);
    }
    TGraph grTemp(m2,inTimes,inVolts);
    TGraph *grPowTemp = makePowerSpectrumVoltsSeconds(&grTemp);
    Double_t *tempFreq = grPowTemp->GetX();
    Double_t *tempPower = grPowTemp->GetY();
    numOut=grPowTemp->GetN();
    for(int i=0;i<numOut;i++) {
      outFreqs[i]=tempFreq[i];
      outPower[i]+=tempPower[i]; //To divide by numSegs or not that is the question
    }
    delete grPowTemp;
  }  
  //  std::cout << m2 << "\t" << numOut << "\t" << numSegs << endl;

  TGraph *grRet = new TGraph (numOut,outFreqs,outPower);
  delete [] inTimes;
  delete [] inVolts;
  delete [] outFreqs;
  delete [] outPower;
  return grRet;
}

void FFTtools::takeDerivative(Int_t numPoints, Double_t *inputX, Double_t *inputY, Double_t *outputX, Double_t *outputY) {
  int count=0;
  for(int samp=1;samp<numPoints;samp++) {
    Double_t deltaT=inputX[samp]-inputX[samp-1];
    Double_t deltaV=inputY[samp]-inputY[samp-1];
    outputX[count]=inputX[samp];
    outputY[count]=deltaV/deltaT;
    count++;
  }
}

//Legacy misspelling
TGraph *FFTtools::getDerviative(TGraph *grIn) 
{
  return getDerivative(grIn);
}

TGraph *FFTtools::getDerivative(TGraph *grIn)
{
  Int_t numPoints=grIn->GetN();
  if(numPoints<2) return NULL;
  Double_t *xVals=grIn->GetX();
  Double_t *yVals=grIn->GetY();
  Double_t *newX = new Double_t [numPoints];
  Double_t *newY = new Double_t [numPoints];
  takeDerivative(numPoints,xVals,yVals,newX,newY);
  TGraph *grDeriv = new TGraph(numPoints-1,newX,newY);
  delete [] newX;
  delete [] newY;
  return grDeriv;

}

TGraph *FFTtools::simplePassBandFilter(TGraph *grWave, Double_t minFreq, Double_t maxFreq)
{

    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT=doFFT(length,oldY);

    int newLength=(length/2)+1;

    //    double fMax = 1/(2*deltaT);  // In Hz
    double deltaF=1/(deltaT*length); //Hz
    deltaF*=1e3; //MHz
    //    std::cout << fMax << "\t" << deltaF << "\t" << deltaT << "\t" << length << std::endl;

    double tempF=0;
    for(int i=0;i<newLength;i++) {
      //      std::cout << tempF << "\t" << theFFT[i].re << "\t" << theFFT[i].im << "\n";
      if(tempF<minFreq || tempF>maxFreq) {
	theFFT[i].re=0;
	theFFT[i].im=0;
      }      
      //      std::cout << tempF << "\t" << theFFT[i].re << "\t" << theFFT[i].im << "\n";
      tempF+=deltaF;
    }

    double *filteredVals = doInvFFT(length,theFFT);


    TGraph *grFiltered = new TGraph(length,oldX,filteredVals);
    delete [] theFFT;
    delete [] filteredVals;
    return grFiltered;

}

TGraph *FFTtools::simpleNotchFilter(TGraph *grWave, Double_t minFreq, Double_t maxFreq)
{

    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT=doFFT(length,oldY);

    int newLength=(length/2)+1;

    //    double fMax = 1/(2*deltaT);  // In Hz
    double deltaF=1/(deltaT*length); //Hz
    deltaF*=1e3; //MHz
    //    std::cout << fMax << "\t" << deltaF << "\t" << deltaT << "\t" << length << std::endl;

    double tempF=0;
    for(int i=0;i<newLength;i++) {
      //      std::cout << tempF << "\t" << theFFT[i].re << "\t" << theFFT[i].im << "\n";
      if(tempF>minFreq && tempF<maxFreq) {
	theFFT[i].re=0;
	theFFT[i].im=0;
      }      
      //      std::cout << tempF << "\t" << theFFT[i].re << "\t" << theFFT[i].im << "\n";
      tempF+=deltaF;
    }

    double *filteredVals = doInvFFT(length,theFFT);


    TGraph *grFiltered = new TGraph(length,oldX,filteredVals);
    delete [] theFFT;
    delete [] filteredVals;
    return grFiltered;

}

TGraph *FFTtools::cropWave(TGraph *grWave, Double_t minTime, Double_t maxTime)
{
  Int_t numPoints=grWave->GetN();
  if(numPoints<1) return NULL;
  Double_t *xVals=grWave->GetX();
  Double_t *yVals=grWave->GetY();
  
  Double_t *outX = new Double_t[numPoints];
  Double_t *outY = new Double_t[numPoints];
  Int_t outPoints=0;

  for(int i=0;i<numPoints;i++) {
    if(xVals[i]>=minTime && xVals[i]<=maxTime) {
      outX[outPoints]=xVals[i];
      outY[outPoints]=yVals[i];
      outPoints++;
    }
  }
  
  TGraph *grCrop = new TGraph(outPoints,outX,outY);
  delete [] outX;
  delete [] outY;
  return grCrop;
}


TGraph *FFTtools::multipleSimpleNotchFilters(TGraph *grWave, Int_t numNotches, Double_t minFreq[], Double_t maxFreq[])
{

    double *oldY = grWave->GetY();
    double *oldX = grWave->GetX();
    double deltaT=oldX[1]-oldX[0];
    int length=grWave->GetN();
    FFTWComplex *theFFT=doFFT(length,oldY);

    int newLength=(length/2)+1;

    //    double fMax = 1/(2*deltaT);  // In Hz
    double deltaF=1/(deltaT*length); //Hz
    deltaF*=1e3; //MHz
    //    std::cout << fMax << "\t" << deltaF << "\t" << deltaT << "\t" << length << std::endl;

    double tempF=0;
    for(int i=0;i<newLength;i++) {
      //      std::cout << tempF << "\t" << theFFT[i].re << "\t" << theFFT[i].im << "\n";
      for(int notch=0;notch<numNotches;notch++) {
	if(tempF>minFreq[notch] && tempF<maxFreq[notch]) {
	  std::cout << "notching " << tempF << " notch num " << notch << " min " << minFreq[notch] << " max " << maxFreq[notch] << std::endl;
	  theFFT[i].re=0;
	  theFFT[i].im=0;
	}      
      }
      //      std::cout << tempF << "\t" << theFFT[i].re << "\t" << theFFT[i].im << "\n";
      tempF+=deltaF;
    }

    double *filteredVals = doInvFFT(length,theFFT);


    TGraph *grFiltered = new TGraph(length,oldX,filteredVals);
    delete [] theFFT;
    delete [] filteredVals;
    return grFiltered;



}

//______________________________________________________________
Double_t FFTtools::getWaveformSNR(TGraph *gr){
  Double_t dummyPeak;
  Double_t dummyRms;
  Double_t snr = getWaveformSNR(gr,dummyPeak,dummyRms);
  return snr;
}




//______________________________________________________________
Double_t FFTtools::getWaveformSNR(TGraph *gr,Double_t &peakToPeak,Double_t &rms)
{
  Int_t nRMS=25;

  Int_t nBins = gr->GetN();
  //  Double_t *xVals = gr->GetX();
  Double_t *yVals = gr->GetY();

  Double_t mean=0.;
  Double_t meanSq=0.;

  for(int i=0;i<nRMS;i++){
    mean+=yVals[i];
    meanSq+=yVals[i]*yVals[i];
  }
  mean/=static_cast<double>(nRMS);
  meanSq/=static_cast<double>(nRMS);

  Int_t trending=3;
  Double_t p2p=0;
  Int_t firstBin=0;
  Double_t y;

  for(int i=0;i<nBins;i++){
    y=yVals[i];
    if(i>0){
      if(y<yVals[i-1] && trending==0){
	if(TMath::Abs(y-yVals[firstBin]>p2p)){
	  p2p=TMath::Abs(y-yVals[firstBin]);
	}
      }
      else if(y<yVals[i-1] && (trending==1 || trending==2)){
	trending=0;
	firstBin=i-1;
	if(TMath::Abs(y-yVals[firstBin]>p2p)){
	  p2p=TMath::Abs(y-yVals[firstBin]);
	}
      }
      else if(y>yVals[i-1] && (trending==0 || trending==2)){
	trending=1;
	firstBin=i-1;
	if(TMath::Abs(y-yVals[firstBin]>p2p)){
	  p2p=TMath::Abs(y-yVals[firstBin]);
	}
      }
      else if(y>yVals[i-1] && trending==1){
	if(TMath::Abs(y-yVals[firstBin]>p2p)){
	  p2p=TMath::Abs(y-yVals[firstBin]);
	}
      }
      else if(y==yVals[i-1]){
	trending=2;
      }
      else if(trending==3){
	if(y<yVals[i-1]){
	  trending=0;
	  firstBin=0;
	}
	if(y>yVals[i-1]){
	  trending=1;
	  firstBin=0;
	}
      }
      else{
	std::cout << "trending cock up!" << std::endl;
	std::cout << "y " << y << " yVals[i] " << yVals[i] << " yVals[i-1] " << yVals[i-1] << std::endl;
	return -1;
      }
    }
  }

  p2p/=2.;
  
  rms=sqrt(meanSq-mean*mean);
  peakToPeak = p2p;

  return p2p/rms;

}




//______________________________________________________________
Double_t FFTtools::getWaveformPeak(TGraph *gr)
{

  Int_t nBins = gr->GetN();
  Double_t yMax = -9999;
  Double_t y,x;

  for(int i=0;i<nBins;i++){
    gr->GetPoint(i,x,y);
    if(y>yMax)
      yMax = y;
  }

  return yMax;
}



//______________________________________________________________
Double_t FFTtools::getEnvelopeSNR(TGraph *gr){
  Double_t dummyPeak;
  Double_t dummyRms;
  Double_t dummyTime;
  Double_t snr = getEnvelopeSNR(gr,dummyPeak,dummyRms,dummyTime);
  return snr;
}




//______________________________________________________________
Double_t FFTtools::getEnvelopeSNR(TGraph *gr,Double_t &peak,Double_t &rms,Double_t &timeOfPeak)
{

  Int_t nRMS=25;

  Int_t nBins = gr->GetN();
  Double_t *xVals = gr->GetX();
  Double_t *yVals = gr->GetY();

  Double_t mean=0.;
  Double_t meanSq=0.;

  for(int i=0;i<nRMS;i++){
    mean+=yVals[i];
    meanSq+=yVals[i]*yVals[i];
  }
  mean/=static_cast<double>(nRMS);
  meanSq/=static_cast<double>(nRMS);

  Double_t p=0;
  Double_t t=0;
  Double_t y;
  Double_t x;

  for(int i=0;i<nBins;i++){
    if(yVals[i]<0){
      std::cout << "this isn't an envelope!!!" << std::endl;
      return -1;
    }
    y=yVals[i];
    x=xVals[i];
    if(i>0){
      if(y>p){
	p=y;
	t=x;
      }
    }
  }

  rms=sqrt(meanSq-mean*mean);
  peak = p;
  timeOfPeak = t;

  return peak/rms;

}



TGraph *FFTtools::convertMagnitudeToTimeDomain(TGraph *inputMag)
{
  Double_t *freqs=inputMag->GetX();
  Double_t *newVmmhz=inputMag->GetY();
  Int_t numFreqs=inputMag->GetN();
  Double_t *vmhz = new Double_t [numFreqs];
  Int_t numPoints=2*(numFreqs-1);
  FFTWComplex *freqDom = new FFTWComplex [numFreqs];
  Double_t df=freqs[1]-freqs[0];
  for(int i=0;i<numFreqs;i++) {
    vmhz[i]=newVmmhz[i]/1e6;
    //    freqDom[i].re=vmhz[i];
    //    freqDom[i].im=0;
    freqDom[i].im=vmhz[i];
    freqDom[i].re=0;
  }
  Double_t *tempV=FFTtools::doInvFFT(numPoints,freqDom);
  Double_t *newT = new Double_t [numPoints];
  Double_t *newV = new Double_t [numPoints];
  Double_t dt=1./(numPoints*df);
  for(int i=0;i<numPoints;i++) {
    if(i<numPoints/2) {
      Int_t tempInd=(numPoints/2)-(i+1);
      //      std::cout << "First: " << i << "\t" << tempInd << "\n";
      newT[tempInd]=-i*dt;
      newV[tempInd]=tempV[i]*df*numPoints/sqrt(2);
      //The sqrt(2) is for the positive and negative frequencies
    }
    else {
      Int_t tempInd=numPoints+(numPoints/2)-(i+1);
      //      std::cout << "Second: " << i << "\t" << tempInd << "\n";
      newT[tempInd]=(numPoints-i)*dt;
      newV[tempInd]=tempV[i]*df*numPoints/sqrt(2);
      //The sqrt(2) is for the positive and negative frequencies
    }
  }
 //  for(int i=0;i<numPoints;i++) {
//     std::cout << i << "\t" << newT[i] << "\n";
//   }
  TGraph *grWave = new TGraph(numPoints,newT,newV);
  delete [] vmhz;
  delete [] newT;
  delete [] newV;
  delete [] freqDom;
  return grWave;

}

Double_t FFTtools::simpleInterploate(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t x)
{
  return ((y2 - y1)* ((x - x1) / (x2-x1)) + y1);
}

TGraph *FFTtools::getConvolution(TGraph *grA, TGraph *grB)
{
  Int_t numPointsA=grA->GetN();
  Int_t numPointsB=grB->GetN();

  Double_t *tA=grA->GetX();
  Double_t *inA=grA->GetY();
  Double_t *inB=grB->GetY();

  Double_t t0=tA[0];
  Double_t deltaT=tA[1]-tA[0];

  //First up make them the same size
  Int_t numPoints=numPointsA;
  if(numPointsB>numPoints) { 
    numPoints=numPointsB;
  }
  Double_t *A= new Double_t[numPoints];
  Double_t *B= new Double_t[numPoints];
  Double_t *T= new Double_t[numPoints];

  for(int i=0;i<numPoints;i++) {
    Int_t indA=(i-(numPoints-numPointsA)/2);
    T[i]=t0+indA*deltaT;
    if(indA<0 || indA>=numPointsA)
	A[i]=0;
    else
      A[i]=inA[indA];

    Int_t indB=(i-(numPoints-numPointsB)/2);
    if(indB<0 || indB>=numPointsB)
	B[i]=0;
    else
      B[i]=inB[indB];

    //    std::cout << i << "\t" << indA << "\t" << indB <<  "\t" << A[i] << "\t" << B[i] << "\n";
  }
  
  Int_t numFreqs=(numPoints/2)+1;
  FFTWComplex *fftA=doFFT(numPoints,A);
  FFTWComplex *fftB=doFFT(numPoints,B);
  // FFTWComplex *fftAB= new FFTWComplex [numFreqs];
  Double_t freq=0;
  Double_t deltaF=1./(numPoints*deltaT);
  for(int i=0;i<numFreqs;i++) {
    // fftAB[i]=(fftA[i]*fftB[i]);
    fftA[i]*=fftB[i];

    // std::cout << "vals\t" << fftAB[i] << "\t" << fftA[i] << "\t" << fftB[i] << std::endl;
    // std::cout << "diffs\t" << (fftAB[i] - fftA[i]) << "..\t" << (fftAB[i] - fftB[i]) << std::endl;
    // std::cout << "elemets\t"
    // 	      << fftAB[i].re - fftA[i].re << ", " << fftAB[i].im - fftA[i].im << "...\t"
    // 	      << fftAB[i].re - fftB[i].re << "\t" << fftAB[i].im - fftB[i].im << std::endl;
    // std::cout << std::endl;
    
    //    std::cout << freq << "\t" << fftAB[i].getAbs() << "\t" << fftA[i].getAbs() << "\t" << fftB[i].getAbs()
    //	      << "\t" << fftA[i].getAbs()*fftB[i].getAbs() << "\n";
    //    std::cout << freq << "\t" << fftAB[i].getPhase() << "\t" << fftA[i].getPhase() << "\t" << fftB[i].getPhase()
    //	      << "\t" << fftA[i].getPhase()+fftB[i].getPhase() << "\n";
    freq+=deltaF;
  }
  
  // Double_t *AB=doInvFFT(numPoints,fftAB);
  Double_t *AB=doInvFFT(numPoints,fftA);  
  Double_t *newAB = new Double_t[numPoints];
  for(int i=0;i<numPoints;i++) {
    if(i<numPoints/2) {
      newAB[i]=AB[(numPoints/2)+i];
      //newAB[i]=AB[i];
    }
    else {
      newAB[i]=AB[i-(numPoints/2)];
      //newAB[i]=AB[i];
    }
  }
  TGraph *grConv = new TGraph(numPoints,T,newAB);
  // delete [] fftAB;
  delete [] fftA;
  delete [] fftB;				
  delete [] A;
  delete [] B;
  delete [] T;
  delete [] AB;
  delete [] newAB;
  
  return grConv;
}



RFSignal *FFTtools::getConvolution(RFSignal *grA, RFSignal *grB)
{
  Int_t numPointsA=grA->GetN();
  Int_t numPointsB=grB->GetN();
  if(numPointsA!=numPointsB) {
    std::cout << "gr method " << numPointsA << " " << numPointsB << "\n";
    TGraph *grRet =getConvolution((TGraph*)grA,(TGraph*)grB);
    RFSignal *rfRet = new RFSignal(grRet);
    delete grRet;
    return rfRet;
  }

//   std::cout << "rf method " << numPointsA << " " << numPointsB << "\n";
  Double_t *tA=grA->GetX();
  Double_t deltaT=tA[1]-tA[0];

  Int_t numPoints=numPointsA;
  Int_t numFreqs=grA->getNumFreqs();
  FFTWComplex *fftA=grA->getComplexNums();
  FFTWComplex *fftB=grB->getComplexNums();
  FFTWComplex *fftAB= new FFTWComplex [numFreqs];
  Double_t freq=0;
  Double_t deltaF=1./(numPoints*deltaT);
  for(int i=0;i<numFreqs;i++) {
    fftAB[i]=fftA[i];
    fftAB[i]*=fftB[i];
    freq+=deltaF;
  }
  
  Double_t *AB=doInvFFT(numPoints,fftAB);
  Double_t *newAB = new Double_t[numPoints];
  for(int i=0;i<numPoints;i++) {
    if(i<numPoints/2) {
      newAB[i]=AB[(numPoints/2)+i];
      //newAB[i]=AB[i];
    }
    else {
      newAB[i]=AB[i-(numPoints/2)];
      //newAB[i]=AB[i];
    }
  }
  RFSignal *grConv = new RFSignal(numPoints,tA,newAB);
  delete [] fftAB;				
  delete [] AB;
  delete [] newAB;  
  return grConv;
}

TGraph *FFTtools::interpolateCorrelateAndAverage(Double_t deltaTInt,Int_t numGraphs, TGraph **grPtrPtr)
{
  TGraph **grInt = new TGraph* [numGraphs];
  for(int i=0;i<numGraphs;i++)
    grInt[i]=getInterpolatedGraph(grPtrPtr[i],deltaTInt);
  TGraph *grAvg=correlateAndAverage(numGraphs,grInt);
  for(int i=0;i<numGraphs;i++)
    delete grInt[i];
  delete [] grInt;
  return grAvg;

}

//___________________________________________//
TGraph *FFTtools::correlateAndAverage(Int_t numGraphs, TGraph **grPtrPtr)
{
  //Assume they are all at same sampling rate

  // Can't correlate and average if there's only one graph.
  // So return 0
  if(numGraphs<2) return NULL;

  // TGraph *grA = grPtrPtr[0];
  TGraph *grA = (TGraph*) grPtrPtr[0]->Clone(); // Make copy of graph rather than using first graph.

  // Copy times from grA into new array.
  Int_t numPoints=grA->GetN();  
  Double_t *timeVals= grA->GetX();
  Double_t *safeTimeVals = new Double_t[numPoints];
  Double_t *sumVolts = new Double_t [numPoints];
  for(int i=0;i<numPoints;i++) 
    safeTimeVals[i]=timeVals[i];  
  

  // Loop over graph array.
  int countWaves=1;
  for(int graphNum=1;graphNum<numGraphs;graphNum++) {
    TGraph *grB = grPtrPtr[graphNum];
    if(grB->GetN()<numPoints)
      numPoints=grB->GetN();
    TGraph *grCorAB = FFTtools::getCorrelationGraph(grA,grB);

    Int_t peakBin = FFTtools::getPeakBin(grCorAB);
    //    Double_t *deltaTVals=grCorAB->GetX();
    //    cout << peakBin << "\t" << grCorAB->GetN() << endl;
    Int_t offset=peakBin-(grCorAB->GetN()/2);
    //    cout << deltaTVals[peakBin] << "\t" << safeTimeVals[offset] << endl;
 
    Double_t *aVolts = grA->GetY();
    Double_t *bVolts = grB->GetY();

    for(int ind=0;ind<numPoints;ind++) {
      int aIndex=ind;
      int bIndex=ind-offset;
      
      if(bIndex>=0 && bIndex<numPoints) {
	sumVolts[ind]=(aVolts[aIndex]+bVolts[bIndex]);
      }
      else {
	sumVolts[ind]=aVolts[aIndex];
      }
    }
    

    TGraph *grComAB = new TGraph(numPoints,safeTimeVals,sumVolts);

    //    delete grB;
    delete grCorAB;
    if(graphNum>1)
      delete grA;
    grA=grComAB;
    countWaves++;

  }
  for(int i=0;i<numPoints;i++) {
    sumVolts[i]/=countWaves;
  }
  Double_t meanVal=TMath::Mean(numPoints,sumVolts);
  for(int i=0;i<numPoints;i++) {
    sumVolts[i]-=meanVal;
  }
  delete grA;
  TGraph *grRet = new TGraph(numPoints,safeTimeVals,sumVolts);
  delete [] safeTimeVals;
  delete [] sumVolts;
  return grRet;
}








//___________________________________________//
double FFTtools::sinc(double x, double eps) 
{
  if (fabs(x)<=eps)
  {
    return 1.; 
  }
  else
  {
    return sin(TMath::Pi()*x) / (TMath::Pi()* x); 
  }

}





//___________________________________________//
void FFTtools::rotate(TGraph *g, int rot) 
{
  rot *=-1; //because of way trigger cell is defined 
  int N = g->GetN(); 
  if (rot < 0) 
  {
    rot+= N*(1+(-rot-1)/N);; 
  }
  rot %= N; 
  std::rotate(g->GetY(), g->GetY() + rot, g->GetY() + N); 
}


//___________________________________________//
double * FFTtools::FFTCorrelation(int length, const FFTWComplex * A, const FFTWComplex * B, FFTWComplex * work, int min_i, int max_i, int order)
{

  int fftlen = length/2 +1; 
  if (max_i <= 0 || max_i >= fftlen) max_i = fftlen-1; 
  if (min_i < 0) min_i = 0; 



  bool work_given = work; 
  if (!work)
  {
    work = new FFTWComplex[fftlen]; 
  }
  else
  {
    memset(work,0,sizeof(FFTWComplex) * fftlen); 
  }

  double rmsA = 0; 
  double rmsB = 0; 
  for(int i=0;i<fftlen;i++)
  {
    double reFFT1=A[i].re;
    double imFFT1=A[i].im;
    double reFFT2=B[i].re; 
    double imFFT2=B[i].im; 



    double weight = 1; 

    if (min_i > 0) 
    {
      weight /= (1 + TMath::Power(double(min_i)/i,2*order)); 
    }

    if (max_i < fftlen - 1) 
    {
      weight /=( 1 +TMath::Power(double(i)/max_i,2*order));     
    }


    if (i > 0)  //don't add DC component to RMS! 
    {
      rmsA += weight*(reFFT1 * reFFT1 + imFFT1 * imFFT1) * 2 / (length*length); 
      rmsB += weight*(reFFT2 * reFFT2 + imFFT2 * imFFT2) * 2 / (length*length); 
    }
    work[i].re=weight*(reFFT1*reFFT2+imFFT1*imFFT2)/length;
    work[i].im=weight*(imFFT1*reFFT2-reFFT1*imFFT2)/length; 
  }

  double *answer=FFTtools::doInvFFT(length,work);

  double norm = (rmsA && rmsB) ? 1./(sqrt(rmsA*rmsB)) : 1; 
  for (int i = 0; i< length; i++) 
  {
    answer[i] *=norm; 
  }

  if (!work_given)
  {
    delete [] work; 
  }

  return answer; 

}


//___________________________________________//
void FFTtools::applyWindow(TGraph *g, const FFTWindowType* win)
{
  win->apply(g->GetN(), g->GetY()); 
}


//___________________________________________//
void FFTtools::inPlaceShift(int N, double *x)
{

  for (int i = 0; i < N/2; i++) 
  {
    double tmp = x[i]; 
    x[i] = x[i+N/2];
    x[i+N/2] = tmp; 
  }
}


//___________________________________________//
double FFTtools::getDt(const TGraph * g, int realN) 
{

  int n = realN ? realN : g->GetN(); 
  double x0 = g->GetX()[0]; 
  double x1 = g->GetX()[g->GetN()-1]; 
  return (x1-x0)/(n-0.5); 
}


//___________________________________________//
void FFTtools::polySubtract(TGraph *g, int order) 
{
  TF1 f("polysub",TString::Format("pol%d",order)); 
  g->Fit(&f,"NQ"); 
  for (int i = 0; i < g->GetN(); i++)
  { 
    g->GetY()[i]-= f.Eval(g->GetX()[i]); 
  } 
}


//___________________________________________//
double * FFTtools::directConvolve(int N, const double * x, int M, const double *h,  double *y, int delay, DirectConvolveEdgeBehavior edge) 
{


  if (!y) y = new double[N]; 
  double start_val = edge == ZEROES_OUTSIDE ? 0 : x[0]; 
  double end_val = edge == ZEROES_OUTSIDE ? 0 : x[N-1]; 
  for (int i = 0; i < N; i++) 
  {
    y[i] = 0; 
    for (int j = delay-M/2; j < delay +(M+1)/2; j++) 
    {
      double X = 0;
      if (i + j < 0) 
      {
        X = start_val ; 
      }
      else if ( i + j >= N) 
      {
        X = end_val; 
      }
      else 
      {
        X = x[i+j]; 
      }

      y[i] += X * h[delay + M/2+j]; 
    }
  }

  return y; 
}


//___________________________________________//
double FFTtools::randomRayleigh(double sigma, TRandom * rng)
{
  double u = rng ? rng->Uniform(0,1) : gRandom->Uniform(0,1); 
  return  sigma * sqrt(-2*log(u)); 
}

//___________________________________________//
template<typename T> 
static void unwrap(size_t N, T * vals, T period) 
{
  T adjust = 0; 
  for (size_t i = 1; i < N; i++) 
  {
    if (vals[i] - vals[i-1] + adjust > period/2)
    {
      adjust -= period; 
    }
    else if (vals[i] - vals[i-1]  + adjust < -period/2)
    {
      adjust += period; 
    }

    vals[i] += adjust; 
  }
}

//___________________________________________//
template <typename T> 
static void wrap(size_t N, T * vals, T period, T center) 
{
  T start = center - period/2; 
  for (size_t i = 0; i < N; i++) 
  {
    int nsub = FFTtools::fast_floor((vals[i]-start) / period); 
    vals[i] -= nsub * period; 
  }
}

//___________________________________________//
template <typename T> 
static void wrap(size_t N, T * vals, T period) 
{
  wrap<T>(N,vals,period,period/2); 

}

//___________________________________________//
void FFTtools::wrap(size_t N, float * vals, float period) 
{
  ::wrap<float>(N, vals, period); 
}

//___________________________________________//
void FFTtools::unwrap(size_t N, float * vals, float period) 
{
  ::unwrap<float>(N, vals, period); 
}

//___________________________________________//
void FFTtools::wrap(size_t N, double * vals, double period) 
{
  ::wrap<double>(N, vals, period); 
}

//___________________________________________//
void FFTtools::unwrap(size_t N, double * vals, double period) 
{
  ::unwrap<double>(N, vals, period); 
}


//___________________________________________//
void FFTtools::wrap(size_t N, float * vals, float period, float center) 
{
  ::wrap<float>(N, vals, period, center); 
}



//___________________________________________//
void FFTtools::wrap(size_t N, double * vals, double period, double center) 
{
  ::wrap<double>(N, vals, period, center); 
}




//___________________________________________//
double FFTtools::evalEvenGraph(const TGraph * g, double t)
{
  double t0 = g->GetX()[0]; 
  if (t < t0) return 0;  //return 0 if before first sample
  double dt = g->GetX()[1] - t0; 

  int bin_low = int ((t-t0)/dt); 

  if (bin_low >= g->GetN()) return 0; 
  if (bin_low ==  g->GetN()-1) return g->GetY()[g->GetN()-1]; 

  int bin_high = bin_low + 1; 
  double frac = (t - (t0 + dt * bin_low)) / dt; 

  double val_low = g->GetY()[bin_low]; 
  double val_high = g->GetY()[bin_high]; 

  return frac * val_high + (1-frac) * val_low; 

}



//___________________________________________//
int FFTtools::saveWisdom(const char * file)
{
  // return fftw_export_wisdom_to_filename(file);

  // No fftw_export_wisdom_to_filename(const char*) in ancient fftw :(
  // Apparently it returns non-zero on success, I will wrap that too.
  FILE* filePtr = fopen(file, "w");
  int retVal = 0;
  if(filePtr){
    fftw_export_wisdom_to_file(filePtr);
    retVal = 1;
    fclose(filePtr);
  }
  else{
    std::cerr << "Warning! " << __PRETTY_FUNCTION__ << " could not open " << file << " in write mode."
	      << std::endl;
    retVal = 0;
  }
  return retVal;
}

//___________________________________________//
int FFTtools::loadWisdom(const char * file)
{

  //  return fftw_import_wisdom_from_filename(file);
  
  // No fftw_export_wisdom_to_filename(const char*) in ancient fftw :(
  // Apparently it returns non-zero on success, I will wrap that too.
  FILE* filePtr = fopen(file, "r");
  int retVal = 0;
  if(filePtr){
    fftw_import_wisdom_from_file(filePtr);
    retVal = 1;
    fclose(filePtr);
  }
  else{
    std::cerr << "Warning! " << __PRETTY_FUNCTION__ << " could not open " << file << " in read mode."
	      << std::endl;
    retVal = 0;
  }
  return retVal;
  
}


//___________________________________________//

#ifdef ENABLE_VECTORIZE
#include "vectorclass.h"
#include "vectormath_trig.h"
#define VEC Vec4d 
#define VEC_N 4
#define VEC_T double
#endif



void FFTtools::stokesParameters(int N, const double * __restrict x, const double *  __restrict xh,
                             const double *  __restrict y, const double *  __restrict yh, 
                             double *I, double * Q, double * U, double * V) 
{
  double sum_x2 = 0; 
  double sum_xh2 = 0;
  double sum_y2 = 0; 
  double sum_yh2 = 0; 
  double sum_xy = 0; 
  double sum_xh_yh =0; 
  double sum_xh_y = 0; 
  double sum_x_yh = 0; 

#ifdef ENABLE_VECTORIZE
  int leftover = N % VEC_N; 
  int nit = N / VEC_N + leftover ? 1 : 0; 


  VEC vx; 
  VEC vxh; 
  VEC vy; 
  VEC vyh; 
  for (int i = 0; i < nit; i++)
  {
    vx.load(x + VEC_N * i);
    vy.load(y + VEC_N * i);
    vxh.load(xh + VEC_N * i);
    vyh.load(yh + VEC_N * i);

    if (i == nit-1 && leftover) 
    {
      vx.cutoff(leftover); 
      vy.cutoff(leftover); 
      vxh.cutoff(leftover); 
      vyh.cutoff(leftover); 
    }

    if (I || Q) 
    {
      sum_x2 += horizontal_add(vx*vx); 
      sum_y2 += horizontal_add(vy*vy); 
      sum_yh2 += horizontal_add(vyh*vyh); 
      sum_xh2 += horizontal_add(vxh*vxh); 
    }

    if (U)
    {
      sum_xy += horizontal_add(vx*vy); 
      sum_xh_yh += horizontal_add(vxh*vyh); 
    }

    if (V)
    {
      sum_x_yh += horizontal_add(vx*vyh); 
      sum_xh_y += horizontal_add(vxh*vy); 
    }
  }

#else
  for (int i = 0; i < N; i++)
  {
    if (I || Q)
    {
      sum_x2 += x[i]*x[i]; 
      sum_xh2 += xh[i]*xh[i]; 
      sum_y2  += y[i]*y[i]; 
      sum_yh2 += yh[i]*yh[i]; 
    }

    if (U) 
    {
      sum_xy += x[i]*y[i]; 
      sum_xh_yh += xh[i]*yh[i]; 
    }
    if (V)
    {
      sum_x_yh += x[i]*yh[i]; 
      sum_xh_y += xh[i]*y[i]; 
    }
  }
#endif

  if (I) *I = (sum_x2 + sum_xh2 + sum_y2 + sum_yh2) / N; 
  if (Q) *Q = (sum_x2 + sum_xh2 - sum_y2 - sum_yh2) / N; 
  if (U) *U = (sum_xy  + sum_xh_yh)*2/ N; 
  if (V) *V = (sum_xh_y  - sum_x_yh)*2/ N; 
}


void FFTtools::dftAtFreqAndMultiples(const TGraph * g, double f, int nmultiples, double * phase, double *amp, double * real, double * imag)
{

  if (nmultiples < 1) return; 
  if (nmultiples == 1) return dftAtFreq(g,f,phase,amp,real,imag); 

  double w = 2 * TMath::Pi() * f; 

  const double * t = g->GetX();
  const double * y = g->GetX();
  int N = g->GetN(); 

  //build up sin table 
  double sin_f[N], cos_f[N]; 
  double sin_nf[N], cos_nf[N]; 

  double vcos = 0;
  double vsin = 0; 

#ifdef ENABLE_VECTORIZE
  VEC vecy; 
  VEC vect; 
  VEC vw = w; 

  int leftover = N % VEC_N; 
  int nit = N / VEC_N + (leftover ? 1 : 0); 

  // do first pass and compute base frequency; 

  for (int i = 0; i < nit; i++) 
  {
    if (i < nit -1 || !leftover)
    {
       vect.load(t+VEC_N*i); 
       vecy.load(y+VEC_N*i); 
    }
    else
    {
       vect.load_partial(leftover, t+VEC_N*i); 
       vecy.load_partial(leftover, y+VEC_N*i); 
    }


    VEC ang = vect * vw; 
    VEC vec_sin, vec_cos; 
    vec_sin = sincos(&vec_cos, ang); 



    if ( i < nit -1 || !leftover) 
    {
      vec_sin.store(sin_f + i * VEC_N); 
      vec_cos.store(cos_f + i * VEC_N); 
    }
    else
    {

      vec_sin.store_partial(leftover, sin_f + i * VEC_N); 
      vec_cos.store_partial(leftover, cos_f + i * VEC_N); 
    }


    vec_sin *= vecy; 
    vec_cos *= vecy; 
    vsin += horizontal_add(vec_sin); 
    vcos += horizontal_add(vec_cos); 
  }


#else
  
  //do first pass and compute base frequency
 

  for (int i = 0; i < N; i++)
  {
    double c,s; 
    SINCOS(w* t[i], &s,&c); 
    sin_f[i] = s; 
    cos_f[i] = c; 
    double v = y[i]; 
    vcos += v*c; 
    vsin += v*s; 
  }
#endif

  if (phase) phase[0] = atan2(vsin,vcos); 
  if (amp) amp[0] = sqrt(vsin*vsin + vcos*vcos); 
  if (real) real[0] = vcos;
  if (imag) imag[0] = vsin; 

  memcpy(sin_nf, sin_f, sizeof(sin_f)); 
  memcpy(cos_nf, cos_f, sizeof(cos_f)); 



  for (int j = 1; j < nmultiples; j++)
  {
    vcos = 0; 
    vsin = 0; 

#ifdef ENABLE_VECTORIZE
    VEC vsin_f; 
    VEC vcos_f; 
    VEC vsin_nf; 
    VEC vcos_nf; 
    for (int i = 0; i < nit; i++)
    {
      if (i < nit -1 || !leftover)
      {
         vect.load(t+VEC_N*i); 
         vecy.load(y+VEC_N*i); 
         vsin_f.load(sin_f + VEC_N * i); 
         vcos_f.load(cos_f + VEC_N * i); 
         vsin_nf.load(sin_nf + VEC_N * i); 
         vcos_nf.load(cos_nf + VEC_N * i); 
      }
      else
      {
         vect.load_partial(leftover, t+VEC_N*i); 
         vecy.load_partial(leftover, y+VEC_N*i); 
         vsin_f.load_partial(leftover, sin_f + VEC_N * i); 
         vcos_f.load_partial(leftover, cos_f + VEC_N * i); 
         vsin_nf.load_partial(leftover, sin_nf + VEC_N * i); 
         vcos_nf.load_partial(leftover, cos_nf + VEC_N * i); 
      }

      VEC vec_sin = vcos_f * vcos_nf - vsin_f * vsin_nf; 
      VEC vec_cos = vcos_f * vsin_nf + vsin_f * vcos_nf; 

      if ( i < nit -1 || !leftover) 
      {
        vec_sin.store(sin_nf + i * VEC_N); 
        vec_cos.store(cos_nf + i * VEC_N); 
      }
      else
      {

        vec_sin.store_partial(leftover, sin_nf + i * VEC_N); 
        vec_cos.store_partial(leftover, cos_nf + i * VEC_N); 
      }

      vec_sin *= vecy; 
      vec_cos *= vecy; 

      vsin += horizontal_add(vec_sin); 
      vcos += horizontal_add(vec_cos); 
    }

#else
    for (int i = 0; i < N; i++)
    {
      double tempcos = cos_f[i] * cos_nf[i] - sin_f[i] * sin_nf[i]; 
      double tempsin = cos_f[i] * sin_nf[i] + sin_f[i] *cos_nf[i]; 
      cos_nf[i] = tempcos; 
      sin_nf[i] = tempsin; 

      vcos += tempcos * y[i]; 
      vsin += tempsin * y[i]; 
    }
#endif

    if (phase) phase[j] = atan2(vsin,vcos); 
    if (amp) amp[j] = sqrt(vsin*vsin + vcos*vcos); 
    if (real) real[j] = vcos; 
    if (imag) imag[j] = vsin; 

  }

}


void FFTtools::dftAtFreq(const TGraph * g, double f, double * phase , double * amp, double * real, double * imag) 
{

  double w = 2 * TMath::Pi() * f;
  double vcos = 0;
  double vsin = 0;

  const double * t = g->GetX(); 
  const double * y = g->GetY(); 
  int N = g->GetN(); 

#ifdef ENABLE_VECTORIZE

  VEC vw= w; 
  VEC vecy; 
  VEC vect ;
  
  int leftover = N % VEC_N;
  int nit = N/VEC_N + (leftover ? 1 : 0); 

  for (int i = 0; i < nit; i++)
  {
    if (i < nit -1 || !leftover)
    {
       vect.load(t+VEC_N*i); 
       vecy.load(y+VEC_N*i); 
    }
    else
    {
       vect.load_partial(leftover, t+VEC_N*i); 
       vecy.load_partial(leftover, y+VEC_N*i); 
    }

    VEC ang = vect * vw; 
    VEC vec_sin, vec_cos; 
    vec_sin = sincos(&vec_cos, ang); 

    vec_sin *= vecy; 
    vec_cos *= vecy; 
    vsin += horizontal_add(vec_sin); 
    vcos += horizontal_add(vec_cos); 
  }
#else
  for (int i = 0; i < N; i++)
	{
    double c,s;
    SINCOS(w*t[i], &s,&c);
    double v = y[i];
    vcos +=c*v;
    vsin +=s*v;
  }
#endif
  if (real)
    *real = vcos; 

  if (imag)
    *imag = vsin; 

  if (phase) 
    *phase=atan2(vsin,vcos);

  if (amp) 
    *amp = sqrt(vsin*vsin + vcos*vcos);  
}




//////////////////////////////////////



FFTWComplex * FFTtools::makeMinimumPhase(int N, const double * G, double mindb)
{
  double minval =  TMath::Power(10,mindb/20); 
  double min = log(minval); 


  double * logmag = (double*) fftw_malloc(sizeof(double) * N); 
  FFTWComplex * output = (FFTWComplex*) fftw_malloc(N * sizeof(FFTWComplex)); 

  for (int i = 0; i < N; i++)
  {
    logmag[i]=  G[i] < minval ? min : log(G[i]); 
  }

  double * phase = getHilbertTransform(N, logmag); 
  
  for (int i = 0; i < N; i++) 
  {
    double mag = G[i] < minval ? minval : G[i]; 

    if ( i == 0 || i == N-1) 
    {
      output[i].re = mag; 
      output[i].im = 0; 

    }
    else
    {
      output[i].re = mag/sqrt(2) * cos(phase[i]); 
      output[i].im = mag/sqrt(2) * sin(phase[i]); //TODO check if this negative sign should be here!!! Might be a convention thing. 
    }
  }

  delete logmag; 
  delete phase; 
  return output; 
}



double FFTtools::checkCausality(int N, const FFTWComplex * signal)
{
  double * real = new double[N]; 
  double * imag = new double[N]; 

  for (int i = 0; i < N; i++) 
  {
    real[i] = signal[i].re; 
    imag[i] = signal[i].im; 
  }

  double * hilbert = getHilbertTransform(N, real); 

  double sum2=0; 
  for (int i = 0; i < N; i++) 
  {
    double diff=hilbert[i]-imag[i]; 
    sum2+= diff*diff; 
  }


  delete real;
  delete hilbert; 
  delete imag; 

  return sqrt(sum2/N); 
}


