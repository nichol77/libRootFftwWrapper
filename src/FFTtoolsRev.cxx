#include <cmath>
#include "FFTWComplex.h"
#include "FFTtools.h"
#include "TMath.h"
#include "FFTtoolsRev.h"


double * FFTtools::getCrossCov(int length, const double * oldY1, const double * oldY2) {

	//  Assure FFTs are done for next power of 2 at or above length.
	int adjLength = pow(2, ceil(log2(length)));
	int newLength = adjLength / 2;  //  Half the adjusted length, for summation purposes and output.

	//  Do FFTs on input arrays.
	FFTWComplex * theFFT1 = FFTtools::doFFT(adjLength, oldY1);
	FFTWComplex * theFFT2 = FFTtools::doFFT(adjLength, oldY2);

	FFTWComplex * tempStep = new FFTWComplex[newLength];
	for (int i = 0; i < newLength; ++i) {

		//  Evaluate the real and imaginary parts of tempStep independently.
    		tempStep[i].re = theFFT1[i].re * theFFT2[i].re + theFFT1[i].im * theFFT2[i].im;
    		tempStep[i].im = theFFT1[i].im * theFFT2[i].re - theFFT1[i].re * theFFT2[i].im;
	}
	double * theOutput = FFTtools::doInvFFT(adjLength, tempStep);

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
	int adjLength = pow(2, ceil(log2(length)));  //  Adjusted length of theOutput.
	for (int i = 0; i < adjLength; ++i) theOutput[i] /= norm;

	return theOutput;
}


TGraph * FFTtools::getCovGraph(const TGraph * gr1, const TGraph * gr2, int * zeroOffset) {

	//  Double graph's length from next power of 2.
	int length1 = gr1 -> GetN();
	int length2 = gr2 -> GetN();
	int length = (length1 > length2) ? length1 : length2;
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

	//  Calculating the normalization. Assuming GetY() arrays for gr1 and gr2 have zero mean.
	double norm = sqrt(gr1 -> GetN() * gr2 -> GetN()) * gr1 -> GetRMS(2) * gr2 -> GetRMS(2);
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


double * FFTtools::getBrickWalledCrossCov(int length, const double * oldY1, const double * oldY2, unsigned int idxStart, unsigned int idxEnd) {

	//  If idxStart > idxEnd, flip them.
	if (idxStart > idxEnd) {

		int oldIdxStart = idxStart, oldIdxEnd = idxEnd;
		idxStart = oldIdxEnd;
		idxEnd = oldIdxStart;
	}

	//  Assure FFTs are done for next power of 2 at or above length.
	int adjLength = pow(2, ceil(log2(length)));
	int newLength = adjLength / 2;  //  Half the adjusted length, for summation purposes and output.

	//  Check the input indices are reasonable.
	if (idxEnd > newLength) {

		fprintf(stderr, "Chosen indices for idxStart and idxEnd outside range of Fourier transform.\n");
		exit(EXIT_FAILURE);
	}

	//  Do FFTs on input arrays.
	FFTWComplex * theFFT1 = FFTtools::doFFT(adjLength, oldY1);
	FFTWComplex * theFFT2 = FFTtools::doFFT(adjLength, oldY2);

	FFTWComplex * tempStep = new FFTWComplex[newLength];
	for (int i = 0; i < newLength; ++i) {

		//  Evaluate the real and imaginary parts of tempStep independently.
		if (i <= idxStart && i >= idxEnd) {

	    		tempStep[i].re = theFFT1[i].re * theFFT2[i].re + theFFT1[i].im * theFFT2[i].im;
    			tempStep[i].im = theFFT1[i].im * theFFT2[i].re - theFFT1[i].re * theFFT2[i].im;
		} else {

			tempStep[i] = 0;
		}
	}
	double * theOutput = FFTtools::doInvFFT(adjLength, tempStep);

	delete [] tempStep;
	delete [] theFFT1, delete [] theFFT2;

	return theOutput;
}


double * FFTtools::getBrickWalledCrossCorr(int length, const double * oldY1, const double * oldY2, int idxStart, int idxEnd) {

	//  For normalization purposes, calculate brick-walled variance of oldY1 and oldY2. Assumes oldY1 and oldY2 are zero meaned.
	double brickWalledY1Var = FFTtools::getBrickWalledCrossCov(length, oldY1, oldY1, idxStart, idxEnd)[0];
	double brickWalledY2Var = FFTtools::getBrickWalledCrossCov(length, oldY2, oldY2, idxStart, idxEnd)[0];

	//  Calculating the normalization.
	double norm = sqrt(brickWalledY1Var * brickWalledY2Var);

	double * theOutput = FFTtools::getBrickWalledCrossCov(length, oldY1, oldY2, idxStart, idxEnd);
	int adjLength = pow(2, ceil(log2(length)));  //  Adjusted length of theOutput.
	for (int i = 0; i < adjLength; ++i) theOutput[i] /= norm;

	return theOutput;
}


TGraph * FFTtools::getBrickWalledCovGraph(const TGraph * gr1, const TGraph * gr2, int * zeroOffset, double fStart, double fEnd) {

	//  Double graph's length from next power of 2.
	int length1 = gr1 -> GetN();
	int length2 = gr2 -> GetN();
	int length = (length1 > length2) ? length1 : length2;
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
	int idxStart = N * fStart * deltaT;  //  fStart = idxStart * deltaF, deltaF * deltaT = 1 / N.
	int idxEnd = N * fEnd * deltaT;
	double * covVals = FFTtools::getBrickWalledCrossCov(N, padY1, padY2, idxStart, idxEnd);

	TGraph * grCov = new TGraph(N);
	double * XCov = grCov -> GetX();
	double * YCov = grCov -> GetY();
	for (int i = 0; i < N; ++i) {

		XCov[i] = (i - N / 2) * deltaT - waveOffset;
		//  Unwrapping so that negative lag time portion precedes positive lag time portion.
		if (i < N / 2) YCov[i + N / 2] = covVals[i];  //  Positive lag time portion.
		else YCov[i - N / 2] = covVals[i];  //  Negative lag time portion.
	}

	delete [] covVals;
	delete [] padY1, delete [] padY2;

	return grCov;
}


TGraph * FFTtools::getBrickWalledCorrGraph(const TGraph * gr1, const TGraph * gr2, int * zeroOffset, double fStart, double fEnd) {

	TGraph * grCorr = FFTtools::getBrickWalledCovGraph(gr1, gr2, zeroOffset, fStart, fEnd);

	//  For normalization purposes, calculate variance for brick-walled gr1 and gr2. Assumes gr1 and gr2 have zero mean.
	double brickWalledGr1Var = FFTtools::getBrickWalledCovGraph(gr1, gr1, zeroOffset, fStart, fEnd) -> Eval(0);
	double brickWalledGr2Var = FFTtools::getBrickWalledCovGraph(gr2, gr2, zeroOffset, fStart, fEnd) -> Eval(0);

	//  Calculating the normalization.
	double norm = sqrt(brickWalledGr1Var * brickWalledGr2Var);
	//  Applying the normalization.
	double * YCorr = grCorr -> GetY();
	for (int i = 0; i < grCorr -> GetN(); ++i) YCorr[i] /= norm;

	return grCorr;
}


TGraph * FFTtools::getInterpBrickWalledCovGraph(const TGraph * gr1, const TGraph * gr2, double deltaT, double fStart, double fEnd) {

	TGraph * gr1Interp = FFTtools::getInterpolatedGraph(gr1, deltaT);
	TGraph * gr2Interp = FFTtools::getInterpolatedGraph(gr2, deltaT);
	TGraph * grOut = FFTtools::getBrickWalledCorrGraph(gr1Interp, gr2Interp);

	delete gr1Interp, delete gr2Interp;

	return grOut;
}


TGraph * FFTtools::getInterpBrickWalledCorrGraph(const TGraph * gr1, const TGraph * gr2, double deltaT, double fStart, double fEnd) {

	TGraph * gr1Interp = FFTtools::getInterpolatedGraph(gr1, deltaT);
	TGraph * gr2Interp = FFTtools::getInterpolatedGraph(gr2, deltaT);
	TGraph * grOut = FFTtools::getBrickWalledCorrGraph(gr1Interp, gr2Interp);

	delete gr1Interp, delete gr2Interp;

	return grOut;
}
