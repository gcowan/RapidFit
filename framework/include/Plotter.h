/**
	@class Plotter

	A class for plotting PDF projections onto histograms

	@author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-13
 */

#ifndef PLOTTER_H
#define PLOTTER_H

#include "IPDF.h"
#include "TH1F.h"
#include "IDataSet.h"
#include "RapidFitIntegrator.h"
#include "TFile.h"

class Plotter
{
	public:
		Plotter();
		Plotter( IPDF*, IDataSet* );
		~Plotter();

		void SetWeightsWereUsed( string ) ;
		void PlotAllObservables( string );
		void PlotObservables( string, vector<string> );
		vector<double> ProjectObservable( DataPoint, string );

	private:
		vector<double> ProjectObservable( DataPoint, string, double, double, int, double& );
		vector<DataPoint> GetDiscreteCombinations( vector<double>&, vector<string>& );
		vector<double> GetStatistics( string, double&, double&, int& );
		void MakeObservablePlots( string, vector<DataPoint>, vector<double>, vector<string>, TFile* );
		void MakePlotCanvas( string, string, TH1F*, double*, double*, int );
	
		bool weightsWereUsed ;
		string weightName ;

		RapidFitIntegrator * pdfIntegrator;
		IPDF * plotPDF;
		IDataSet * plotData;
};

#endif
