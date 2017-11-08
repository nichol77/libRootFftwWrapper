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

	/*  This does what FFTtools::getCorrelation() tries to do. By dividing getCrossCov() by the geometric mean
	 *  of sum of the squares of oldY1 and oldY2, we then have the cross-correlation. We assume oldY1 and oldY2
	 *  are of the same length.
	 */
	double * getCrossCorr(int length, const double * oldY1, const double * oldY2);

	/*  Since FFTtools::getCorrelation() is replaced above with crossCov(), it stands to reason
	 *  that method changes means replacing getCorrelationGraph() as well. We do that here with covGraph().
	 */
	TGraph * getCrossCovGraph(const TGraph * gr1, const TGraph * gr2, int & zeroOffset);

	/*  This is to replace FFTtools::getCorrelation(), or FFTtools::getNormalisedCorrelation() since this
	 *  returns the cross-correlation graph as opposed to the cross-covariance graph.
	 */
	TGraph * getCrossCorrGraph(const TGraph * gr1, const TGraph * gr2, int & zeroOffset);
}

#endif
