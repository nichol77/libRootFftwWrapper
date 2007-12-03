
#ifndef FFTTOOLS_H
#define FFTTOOLS_H

//#include <fftw3.h>
#include "TGraph.h"
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"

//My includes
#include "FFTWComplex.h"


class FFTtools
{
public:
    FFTtools();
    ~FFTtools();    
    static double getAbs(FFTWComplex &theNum);
    static double *doInvFFT(int length, FFTWComplex *theInput);
    static FFTWComplex *doFFT(int length,double *theInput);
    

    static TGraph *getBoxCar(TGraph *grWave, Int_t halfWidth);
    static TGraph *getHilbertTransform(TGraph *grWave);
    static TGraph *getHilbertEnvelope(TGraph *grWave);
    static TGraph *makePowerSpectrum(TGraph *grWave);
    static TGraph *makePowerSpectrumPeriodogram(TGraph *grWave);
    static TGraph *makePowerSpectrumVoltsSeconds(TGraph *grWave);
    static TGraph *makePowerSpectrumVoltsSecondsdB(TGraph *grWave);
    static TGraph *makePowerSpectrumVoltsSecondsPadded(TGraph *grWave, Int_t padFactor=4);
    static TGraph *makePowerSpectrumVoltsSecondsPaddeddB(TGraph *grWave, Int_t padFactor=4);
    static TGraph *makeRawPowerSpectrum(TGraph *grWave);
    static TGraph *getCorrelationGraph(TGraph *gr1, TGraph *gr2);
    static double *getCorrelation(TGraph *gr1, TGraph *gr2,int firstIndex,int lastIndex);
    static double *getCorrelation(int length,float *oldY1, float *oldY2);
    static double *getCorrelation(int length,double *oldY1, double *oldY2);
    
    static TGraph *makeInverseInverseSpectrum(TGraph *grWave);
    static TGraph *combineGraphsUsingFFTs(Int_t numGraphs, TGraph **grPtr,double *theWeights=0);
    static Double_t *combineValuesUsingFFTs(Int_t numArrays, Double_t **thePtrPtr, Int_t eachLength);


    //Utility functions
    static Double_t sumPower(TGraph *gr,Int_t firstBin,Int_t lastBin);
    static Double_t integratePower(TGraph *gr,Int_t firstBin,Int_t lastBin);
    static Double_t sumVoltageSquared(TGraph *gr,Int_t firstBin,Int_t lastBin);
    static Double_t integrateVoltageSquared(TGraph *gr,Int_t firstBin,Int_t lastBin);
    static Int_t getPeakBin(TGraph *gr); 
    
    //Graph returning utility funcs
    static TGraph *getSimplePowerEnvelopeGraph(TGraph *gr);
    static TGraph *smoothFFT(TGraph *gr,Int_t factor) ;
    static TGraph *subtractGraphs(TGraph *grA, TGraph *grB);
    static TGraph *divideGraphs(TGraph *grA, TGraph *grB);
    static TGraph *ratioSubtractOneGraphs(TGraph *grA, TGraph *grB) ;
    static TGraph *dbGraphs(TGraph *grA, TGraph *grB);
    static TGraph *padWave(TGraph *grA, Int_t padFactor);

};
   
#endif //FFTTOOLS_H
