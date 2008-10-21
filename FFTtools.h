
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

    
    //Interpolation Routines that use ROOT::Math::Interpolator
    static TGraph *getInterpolatedGraph(TGraph *grIn, Double_t deltaT);

    //Worker functions that use FFTW
    static double getAbs(FFTWComplex &theNum);
    static double *doInvFFT(int length, FFTWComplex *theInput);
    static FFTWComplex *doFFT(int length,double *theInput);    
    static double *getCorrelation(TGraph *gr1, TGraph *gr2,int firstIndex,int lastIndex);
    static double *getCorrelation(int length,float *oldY1, float *oldY2);
    static double *getCorrelation(int length,double *oldY1, double *oldY2);
    static Double_t *combineValuesUsingFFTs(Int_t numArrays, Double_t **thePtrPtr, Int_t eachLength);

    //Higher level functions that take and return TGraphs
    static TGraph *makePowerSpectrum(TGraph *grWave);
    static TGraph *makePowerSpectrumPeriodogram(TGraph *grWave);
    static TGraph *makePowerSpectrumVoltsSeconds(TGraph *grWave);
    static TGraph *makePowerSpectrumVoltsSecondsBartlett(TGraph *grWave);
    static TGraph *makePSVSBartlettPaddedOverlap(TGraph *grWave, Int_t padFactor=4, Int_t numFreqs=64);
    static TGraph *makePowerSpectrumVoltsSecondsdB(TGraph *grWave);
    static TGraph *makePowerSpectrumMilliVoltsNanoSecondsdB(TGraph *grWave);
    static TGraph *makePowerSpectrumVoltsSecondsPadded(TGraph *grWave, Int_t padFactor=4);
    static TGraph *makePowerSpectrumVoltsSecondsPaddeddB(TGraph *grWave, Int_t padFactor=4);
    static TGraph *makeRawPowerSpectrum(TGraph *grWave);
    static TGraph *getCorrelationGraph(TGraph *gr1, TGraph *gr2);
    static TGraph *getInterpolatedCorrelationGraph(TGraph *grIn1, TGraph *grIn2, Double_t deltaT);
    static TGraph *makeInverseInverseSpectrum(TGraph *grWave);
    static TGraph *combineGraphsUsingFFTs(Int_t numGraphs, TGraph **grPtr,double *theWeights=0);
    static TGraph *getBoxCar(TGraph *grWave, Int_t halfWidth);
    static TGraph *getHilbertTransform(TGraph *grWave);
    static TGraph *getHilbertEnvelope(TGraph *grWave);

    //Utility functions (not necessarily FFT related but they'll live here for now
    static Double_t sumPower(TGraph *gr,Int_t firstBin=-1,Int_t lastBin=-1);
    static Double_t integratePower(TGraph *gr,Int_t firstBin=-1,Int_t lastBin=-1);
    static Double_t sumVoltageSquared(TGraph *gr,Int_t firstBin,Int_t lastBin);
    static Double_t integrateVoltageSquared(TGraph *gr,Int_t firstBin=-1,Int_t lastBin=-1);
    static Int_t getPeakBin(TGraph *gr); 
    static Double_t getPeakVal(TGraph *gr, int *index=0);
    static Double_t getPeakSqVal(TGraph *gr, int *index=0);
    static void getPeakRmsSqVal(TGraph *gr, Double_t &peak, Double_t &rms, int *index=0);
    static void getPeakRmsRectified(TGraph *gr, Double_t &peak, Double_t &rms, int *index=0);

    //Graph returning utility funcs
    static TGraph *getSimplePowerEnvelopeGraph(TGraph *gr);
    static TGraph *smoothFFT(TGraph *gr,Int_t factor) ;
    static TGraph *subtractGraphs(TGraph *grA, TGraph *grB);
    static TGraph *divideGraphs(TGraph *grA, TGraph *grB);
    static TGraph *ratioSubtractOneGraphs(TGraph *grA, TGraph *grB) ;
    static TGraph *dbGraphs(TGraph *grA, TGraph *grB);
    static TGraph *padWave(TGraph *grA, Int_t padFactor);
    static TGraph *rectifyWave(TGraph *gr, Int_t makeNeg=0);

    //Window functions
    static Double_t  bartlettWindow(Int_t j, Int_t n);
    static Double_t welchWindow(Int_t j, Int_t n);



};
   
#endif //FFTTOOLS_H