#ifndef FFTTOOLSREV_H
#define FFTTOOLSREV_H

/*  Tools written by John Russell <jwruss@hawaii.edu> which make revisions to existing FFTtools
 *  functions such that they exhibit conventional behavior.
 */


#include "TGraph.h"


namespace FFTtools {

	/*  It appears that FFTtools::getCorrelation(int length, const double * oldY1, const double * oldY2)
	 *  has been normalized inappropriately by fNpoints / 4, so here we replace it and instead call it
	 *  FFTtools::crossCov(), for "cross-covariance". We call its normalization the cross-correlation.
	 *  We assume oldY1 and oldY2 are of the same length.
	 */
	double * getCrossCov(int length, const double * oldY1, const double * oldY2);

	/*  This does what FFTtools::getCorrelation() tries to do. By dividing FFTtools::getCrossCov() by the
	 *  geometric mean of sum of the squares of oldY1 and oldY2, we then have the cross-correlation.
	 *  We assume oldY1 and oldY2 are of the same length.
	 */
	double * getCrossCorr(int length, const double * oldY1, const double * oldY2);

	/*  Since FFTtools::getCorrelation() is replaced above with FFTtools::getCrossCov(), it stands to reason
	 *  that method changes means replacing FFTtools::getCorrelationGraph() as well. We do that here with
	 *  FFTtools::getCovGraph().
	 */
	TGraph * getCovGraph(const TGraph * gr1, const TGraph * gr2, int * zeroOffset = 0);
	
	/*  This is to replace FFTtools::getCorrelation(), or FFTtools::getNormalisedCorrelation() since this
	 *  returns the cross-correlation graph as opposed to the cross-covariance graph.
	 */
	TGraph * getCorrGraph(const TGraph * gr1, const TGraph * gr2, int * zeroOffset = 0);

	/*  Prior to calling FFTtools::getCovGraph(), the input gr1 and gr2 are Akima spline interpolated to a
	 *  uniform time step size of deltaT nanoseconds.
	 */
	TGraph * getInterpCovGraph(const TGraph * gr1, const TGraph * gr2, double deltaT);

	/*  Prior to calling FFTtools::getCorrGraph(), the input gr1 and gr2 are Akima spline interpolated to a
	 *  uniform time step size of deltaT nanoseconds.
	 */
	TGraph * getInterpCorrGraph(const TGraph * gr1, const TGraph * gr2, double deltaT);

	/*  This does the same as FFTtools::getCrossCov() above, except it applies a brick-wall filter in the
	 *  Fourier domain over the input index range before inverting back in the original domain.
	 *  Disclaimer: You have to know the range of bin indices in the Fourier domain you want to consider!
	 */
	double * getBrickWalledCrossCov(int length, const double * oldY1, const double * oldY2, unsigned int idxStart, unsigned int idxEnd);

	/*  What FFTtools::getCrossCorr() does for FFTtools::getCrossCov(), this does for FFTtools::getWeightedCrossCov().
	 *  Disclaimer: You have to know the range of bin indices in the Fourier domain you want to consider!
	 */
	double * getBrickWalledCrossCorr(int length, const double * oldY1, const double * oldY2, unsigned int idxStart, unsigned int idxEnd);

	/*  Basically does FFTtools::getCovGraph(), but uses FFTtools::getBrickWalledCrossCov() instead.
	 *  fStart and fEnd are the frequencies corresponding to the beginning and end of the brick wall filter.
	 */
	TGraph * getBrickWalledCovGraph(const TGraph * gr1, const TGraph * gr2, int * zeroOffset = 0, double fStart, double fEnd);

	/*  Basically does FFTtools::getCorrGraph(), but uses FFTtools::getBrickWalledCrossCov() instead.
	 *  fStart and fEnd are the GHz frequencies corresponding to the beginning and end of the brick wall filter.
	 */
	TGraph * getBrickWalledCorrGraph(const TGraph * gr1, const TGraph * gr2, int * zeroOffset = 0, double fStart, double fEnd);

	/*  Prior to calling FFTtools::getBrickWalledCovGraph(), the input gr1 and gr2 are Akima spline interpolated to a
	 *  uniform time step size of deltaT nanoseconds.
	 *  fStart and fEnd are the GHz frequencies corresponding to the beginning and end of the brick wall filter.
	 */
	TGraph * getBrickWalledWeightedCovGraph(const TGraph * gr1, const TGraph * gr2, double deltaT, double fStart, double fEnd);

	/*  Prior to calling FFTtools::getBrickWalledCorrGraph(), the input gr1 and gr2 are Akima spline interpolated to a
	 *  uniform time step size of deltaT nanoseconds.
	 *  fStart and fEnd are the GHz frequencies corresponding to the beginning and end of the brick wall filter.
	 */
	TGraph * getInterpBrickWalledCorrGraph(const TGraph * gr1, const TGraph * gr2, double deltaT, double fStart, double fEnd);
}

#endif
