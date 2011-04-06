/**
  @namespace GoodnessOfFit

  @author Greig A Cowan greig.cowan@cern.ch
  @date 2011-04-06
 */
#ifndef RAPIDFIT_GOODNESSOFFIT
#define RAPIDFIT_GOODNESSOFFIT
#include "IDataSet.h"
#include "IPDF.h"
#include "PhaseSpaceBoundary.h"
#include "TH1D.h"

namespace GoodnessOfFit 
{
	double pValueFromPoint2PointDissimilarity( IDataSet * data, IDataSet * mc );
	double calculateTstatistic( IDataSet * data, IDataSet * mc);
	vector<double> permutation( IDataSet * data, IDataSet * mc, int nPerm );
	double permutationCore( IDataSet * data, IDataSet * mc, int iteration );
	double sumEvents( IDataSet * data );
	double sumDataMCEvents( IDataSet * data, IDataSet * mcData );
	void plotUstatistic( IPDF * pdf, IDataSet * data, PhaseSpaceBoundary * phase, string plot );
	void calculateUstatistic( IPDF * pdf, IDataSet * data, PhaseSpaceBoundary * phase, TH1D * distances );
	double getDistance( DataPoint * x, DataPoint * y );
}

#endif
