#include "FFTtools.h"

using namespace std;

FFTtools::FFTtools()
{
}

FFTtools::~FFTtools()
{
}

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





// FFTWComplex *FFTtools::doFFT(int length, double *theInput) {
//   //Here is what the sillyFFT program should be doing;    
//     fftw_complex *theOutput = new fftw_complex [(length/2)+1];
//     double *newInput = new double [length]; 
    
//     //cout << length  << " " << input[0] << " " << theFFT[0] << endl;
//     fftw_plan thePlan = fftw_plan_dft_r2c_1d(length,newInput,theOutput,FFTW_MEASURE);
//     if(!thePlan) {
// 	cout << "Bollocks" << endl;
//     }
        
//     for(int i=0;i<length;i++) {
// 	newInput[i]=theInput[i];
//     }

//     for (int i=0;i<(length/2)+1;i++) {
// 	theOutput[i][0]=0.;
// 	theOutput[i][1]=0.;
//     }
//     fftw_execute(thePlan);
//     delete [] newInput; 
//     fftw_destroy_plan(thePlan);
    

//     FFTWComplex *myOutput = new FFTWComplex [(length/2)+1];
//     for (int i=0;i<(length/2)+1;i++) {
// 	myOutput[i].re=theOutput[i][0];
// 	myOutput[i].im=theOutput[i][1];
//     }
//     delete [] theOutput;
//     return myOutput;
// }


FFTWComplex *FFTtools::doFFT(int length, double *theInput) {
  //Here is what the sillyFFT program should be doing;    

  makeNewPlanIfNeeded(length);
  memcpy(FFTtools::fReals[length], theInput, sizeof(double)*length);
  fftw_execute(FFTtools::fRealToComplex[length]);

  const int numFreqs= (length/2)+1;
  FFTWComplex *myOutput = new FFTWComplex [numFreqs];
  memcpy(myOutput, FFTtools::fComplex[length], sizeof(fftw_complex)*numFreqs);
  return myOutput;
}



// double *FFTtools::doInvFFT(int length, FFTWComplex *theInput) {
//   // This is what sillyFFT should be doing
//   //    //Takes account of normailisation 
//   // Although note that fftw_plan_dft_c2r_1d assumes that the frequency array is only the positive half, so it gets scaled by sqrt(2) to account for symmetry
//     fftw_complex *newInput = new fftw_complex [(length/2)+1];
//     double *theOutput = new double [length]; 
//     fftw_plan thePlan = fftw_plan_dft_c2r_1d(length,newInput,theOutput,FFTW_MEASURE);
           
//     for(int i=0;i<((length/2)+1);i++) {
// 	newInput[i][0]=theInput[i].re;
// 	newInput[i][1]=theInput[i].im;
//     }
//     for(int i=0;i<length;i++) {
// 	theOutput[i]=0;
//     }
    
//     fftw_execute(thePlan);
//     delete [] newInput; 
//     fftw_destroy_plan(thePlan);
//     for(int i=0;i<length;i++) {
// 	theOutput[i]/=length;
//     }
//     return theOutput;
// }

double *FFTtools::doInvFFT(int length, FFTWComplex *theInput) {
  // This is what sillyFFT should be doing
  //    //Takes account of normailisation 
  // Although note that fftw_plan_dft_c2r_1d assumes that the frequency array is only the positive half, so it gets scaled by sqrt(2) to account for symmetry


  makeNewPlanIfNeeded(length);
  std::complex<double>* tempVals = (std::complex<double>*) FFTtools::fComplex[length];
  const int numFreqs= (length/2)+1;
  for(int i=0; i<numFreqs; i++){
    tempVals[i].real(theInput[i].re);
    tempVals[i].imag(theInput[i].im);    
  }

  // Do inverse FFT.
  fftw_execute(FFTtools::fComplexToReal[length]);

  /* Normalization needed on the inverse transform */
  double* invFftOutPtr = FFTtools::fReals[length];
  for(int i=0; i<length; i++){
    invFftOutPtr[i]/=length;
  }

  // Copy output array
  double *theOutput = new double [length];
  memcpy(theOutput, invFftOutPtr, sizeof(double)*length);
  return theOutput;
}





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

TGraph *FFTtools::getHilbertTransform(TGraph *grWave)
{
  double *oldY = grWave->GetY();
  double *oldX = grWave->GetX();
  //    double deltaT=oldX[1]-oldX[0];
  int length=grWave->GetN();
  FFTWComplex *theFFT=doFFT(length,oldY);  
  int newLength=(length/2)+1;
  for(int i=0;i<newLength;i++) {
    double tempIm=theFFT[i].im;
    theFFT[i].im=theFFT[i].re;
    theFFT[i].re=-1*tempIm;
  }
  double *hilbert = doInvFFT(length,theFFT);
  
  TGraph *grHilbert = new TGraph(length,oldX,hilbert);
  delete [] hilbert;
  delete [] theFFT;

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
  FFTWComplex *fftAB= new FFTWComplex [numFreqs];
  Double_t freq=0;
  Double_t deltaF=1./(numPoints*deltaT);
  for(int i=0;i<numFreqs;i++) {
    fftAB[i]=fftA[i]*fftB[i];
    //    std::cout << freq << "\t" << fftAB[i].getAbs() << "\t" << fftA[i].getAbs() << "\t" << fftB[i].getAbs()
    //	      << "\t" << fftA[i].getAbs()*fftB[i].getAbs() << "\n";
    //    std::cout << freq << "\t" << fftAB[i].getPhase() << "\t" << fftA[i].getPhase() << "\t" << fftB[i].getPhase()
    //	      << "\t" << fftA[i].getPhase()+fftB[i].getPhase() << "\n";
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
  TGraph *grConv = new TGraph(numPoints,T,newAB);
  delete [] fftAB;
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
    fftAB[i]=fftA[i]*fftB[i];
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
    grInt[i]=getInterpolatedGraph(grPtrPtr[0],deltaTInt);
  TGraph *grAvg=correlateAndAverage(numGraphs,grInt);
  for(int i=0;i<numGraphs;i++)
    delete grInt[i];
  delete [] grInt;
  return grAvg;

}

TGraph *FFTtools::correlateAndAverage(Int_t numGraphs, TGraph **grPtrPtr)
{
  //Assume they are all at same sampling rate
  if(numGraphs<2) return NULL;
  TGraph *grA = grPtrPtr[0];
  Int_t numPoints=grA->GetN();  
  Double_t *timeVals= grA->GetX();
  Double_t *safeTimeVals = new Double_t[numPoints];
  Double_t *sumVolts = new Double_t [numPoints];
  for(int i=0;i<numPoints;i++) 
    safeTimeVals[i]=timeVals[i];  
  
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











////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PRIVATE THINGS
////////////////////////////////////////////////////////////////////////////////////////////////////////////



/* Define static members */
/* https://stackoverflow.com/questions/18433752/c-access-private-static-member-from-public-static-method */
std::map<int, fftw_plan> FFTtools::fRealToComplex;
std::map<int, fftw_plan> FFTtools::fComplexToReal;
std::map<int, double*> FFTtools::fReals;
std::map<int, std::complex<double>*> FFTtools::fComplex;

bool FFTtools::makeNewPlanIfNeeded(int len){
  /* 
     Function which checks whether we've encountered a request to do an FFT of this length before.
     If we haven't then we need a new plan!
  */

  std::map<int,fftw_plan>::iterator it = FFTtools::fRealToComplex.find(len);
  if(it==FFTtools::fRealToComplex.end()){
    // std::cout << len << "\t" << threadInd << std::endl;
    FFTtools::fReals[len] = (double*) fftw_malloc(sizeof(double)*len);
    FFTtools::fComplex[len] = (std::complex<double>*) fftw_malloc(sizeof(fftw_complex)*len);
    FFTtools::fRealToComplex[len] = fftw_plan_dft_r2c_1d(len,
							 FFTtools::fReals[len],
							 (fftw_complex*)FFTtools::fComplex[len],
							 FFTW_MEASURE);
    FFTtools::fComplexToReal[len] = fftw_plan_dft_c2r_1d(len,
							 (fftw_complex*)FFTtools::fComplex[len],
							 FFTtools::fReals[len],
							 FFTW_MEASURE);
    return true;
  }
  else{
    return false;
  }
}
