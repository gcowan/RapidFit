
#ifndef MULTIDIMCHI2_H
#define MULTIDIMCHI2_H

#include "THn.h"
#include <vector>
#include "ObservableRef.h"
#include "PhaseSpaceBoundary.h"
#include "IDataSet.h"
#include "PDFWithData.h"
#include <string>
#include "IPDF.h"

using namespace::std;

struct ThisObsBinning
{
	string ObservableName;
	vector<double> binCenters;
	double thisMin;
	double thisMax;
	double thisStepSize;
	unsigned int theseBins;
};

class MultiDimChi2
{
	public:
		MultiDimChi2( vector<PDFWithData*> allObjects, PhaseSpaceBoundary* thisBound, vector<string> wantedObservables );

		void PerformMuiltDimTest();

	private:

		vector<double> x_min;
		vector<double> x_max;
		vector<ObservableRef> goodObservables;
		vector<int> x_bins;

		THnD* internalHisto;

		vector<ThisObsBinning*> theseDimensions;

		void ConstructBinCenters();

		void ConstructInternalHisto( vector<string> wantedObservables, PhaseSpaceBoundary* thisBound );

		void ConstructAllCoordinates();

		unsigned int nDim;

		vector<vector<double> >* allBinCenters;

		void AddCoordinates( unsigned int thisDim );

		vector<IPDF*> allPDFs;
		vector<IDataSet*> allDataSets;

		vector<PDFWithData*> allPDFData;

		vector<PhaseSpaceBoundary*> allBoundaries;

		double CalcChi2( vector<double> expected_events, vector<double> observed_events );

		double CalculateTotalExpected( vector<double> thisBinCenter );

		void ConstructBoundaries( PhaseSpaceBoundary* totalPhaseSpace, vector<string> );

		double PDF2DataNormalisation( unsigned int PDFNum, const unsigned int combinationIndex, DataPoint* thisDataPoint );

		unsigned int data_binning;

		vector<vector<double> > ratioOfIntegrals;
		vector<vector<double> > combinationIntegrals;
		vector<double> weightNorms;

		void ConstructIntegralsRatios( vector<string> wantedObservables );

		double CalculateRange( PhaseSpaceBoundary* thisBound );

		void populateAllObjects( vector<PDFWithData*> allObjects );

		double CorrectIntegral( double input_Integral, DataPoint* thisPoint, PhaseSpaceBoundary* thisPhaseSpace, RapidFitIntegrator* thisPDFIntegrator );

		double CorrectYield( IDataSet* thisSet, DataPoint* thisPoint );
};

#endif


