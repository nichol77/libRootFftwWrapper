#include <cmath>
#include "FFTWComplex.h"
#include "FFTtools.h"
#include "TMath.h"
#include "FFTtoolsRev.h"


double * FFTtools::getCrossCov(int length, const double * oldY1, const double * oldY2) {

	//  Do FFTs on input arrays.
	FFTWComplex * theFFT1 = FFTtools::doFFT(length, oldY1);
	FFTWComplex * theFFT2 = FFTtools::doFFT(length, oldY2);

	int newLength = length / 2 + 1;  //  Following the same convention in the initial FFTtools.cxx. Follows from FFTW3 using unused padding when "newLength" is odd numbered.

	FFTWComplex * tempStep = new FFTWComplex[newLength];

	for (int i = 0; i < newLength; ++i) {

        	//  Evaluate the real and imaginary parts of tempStep independently.
    		tempStep[i].re = theFFT1[i].re * theFFT2[i].re + theFFT1[i].im * theFFT2[i].im;
    		tempStep[i].im = theFFT1[i].im * theFFT2[i].re - theFFT1[i].re * theFFT2[i].im;
    	}

	double * theOutput = FFTtools::doInvFFT(length, tempStep);

	delete [] tempStep;
	delete [] theFFT1, delete [] theFFT2;

	return theOutput;
}


double * FFTtools::getCrossCorr(int length, const double * oldY1, const double * oldY2) {

	//  Calculating the normalization. Assuming the arrays have a zero mean.
	double oldY1SqSum = 0, oldY2SqSum = 0;

	for (int i = 0; i < length; ++i) {

		oldY1SqSum += oldY1[i] * oldY1[i];
		oldY2SqSum += oldY2[i] * oldY2[i];
	}

	double norm = sqrt(oldY1SqSum * oldY2SqSum);

	double * theOutput = FFTtools::getCrossCov(length, oldY1, oldY2);

	for (int i = 0; i < length; ++i) theOutput[i] /= norm;

	return theOutput;
}


TGraph * FFTtools::getCovGraph(const TGraph * gr1, const TGraph * gr2, int * zeroOffset) {

	//  Double graph's length from next power of 2.
	int length1 = gr1 -> GetN();
	int length2 = gr2 -> GetN();
	int length = (length1 >= length2) ? length1 : length2;
	int N = pow(2, ceil(log2(length)) + 1);

	//  Sample index relative to gr1.
	int firstRealSamp = (N - length1) / 2;

	//  Zero padding the waveforms.
	double * padY1 = new double[N];
	double * padY2 = new double[N];
	for (int i = 0; i < N; ++i) {

		padY1[i] = (i < firstRealSamp || i >= firstRealSamp + length1) ? 0 : gr1 -> GetY()[i - firstRealSamp];
		padY2[i] = (i < firstRealSamp || i >= firstRealSamp + length2) ? 0 : gr2 -> GetY()[i - firstRealSamp];
	}

	double deltaT = (gr1 -> GetX()[length1 - 1] - gr1 -> GetX()[0]) / (length1 - 1);
	double waveOffset = gr2 -> GetX()[0] - gr1 -> GetX()[0];

	//  Optional return of zero offset index.
	if (zeroOffset) * zeroOffset = N / 2 - int(waveOffset / deltaT);

	//  Generating the TGraph to be returned.
	double * covVals = FFTtools::getCrossCov(N, padY1, padY2);

	TGraph * grCov = new TGraph(N);
	double * XCov = grCov -> GetX();
	double * YCov = grCov -> GetY();
	for (int i = 0; i < N; ++i) {

		XCov[i] = (i - N / 2) * deltaT - waveOffset;
		//  Weird wrapping follows. Perhaps from the zero padding at the beginning of the arrays?
		if (i < N / 2) YCov[i + N / 2] = covVals[i];  //  Positive lag time portion.
		else YCov[i - N / 2] = covVals[i];  //  Negative lag time portion.
	}

	delete [] covVals;
	delete [] padY1, delete [] padY2;

	return grCov;
}


TGraph * FFTtools::getCorrGraph(const TGraph * gr1, const TGraph * gr2, int * zeroOffset) {

	TGraph * grCorr = FFTtools::getCovGraph(gr1, gr2, zeroOffset);

	//  Calculating normalization.
	double * Y1 = gr1 -> GetY();
	double totPow1 = 0;
	for (int i = 0; i < gr1 -> GetN(); ++i) totPow1 += Y1[i] * Y1[i];

	double * Y2 = gr2 -> GetY();
	double totPow2 = 0;
	for (int i = 0; i < gr2 -> GetN(); ++i) totPow2 += Y2[i] * Y2[i];

	double norm = sqrt(totPow1 * totPow2);

	//  Applying the normalization.
	double * YCorr = grCorr -> GetY();
	for (int i = 0; i < grCorr -> GetN(); ++i) YCorr[i] /= norm;

	return grCorr;
}


TGraph * FFTtools::getInterpCovGraph(const TGraph * gr1, const TGraph * gr2, double deltaT) {

	TGraph * gr1Interp = FFTtools::getInterpolatedGraph(gr1, deltaT);
	TGraph * gr2Interp = FFTtools::getInterpolatedGraph(gr2, deltaT);
	TGraph * grOut = FFTtools::getCovGraph(gr1Interp, gr2Interp);

	delete gr1Interp, delete gr2Interp;

	return grOut;
}


TGraph * FFTtools::getInterpCorrGraph(const TGraph * gr1, const TGraph * gr2, double deltaT) {

	TGraph * gr1Interp = FFTtools::getInterpolatedGraph(gr1, deltaT);
	TGraph * gr2Interp = FFTtools::getInterpolatedGraph(gr2, deltaT);
	TGraph * grOut = FFTtools::getCorrGraph(gr1Interp, gr2Interp);

	delete gr1Interp, delete gr2Interp;

	return grOut;
}
