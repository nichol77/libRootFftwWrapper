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
}

#endif
