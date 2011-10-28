/**
        @class DataPoint

        Holds all observables for a given event

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/


#ifndef DATA_POINT_H
#define DATA_POINT_H

//	RapidFit Headers
#include "Observable.h"
#include "ObservableRef.h"
//	System Headers
#include <vector>
#include <string>

using namespace std;

class DataPoint
{
	public:
		DataPoint();
		DataPoint( vector<string> );
		~DataPoint();

		vector<string> GetAllNames();
		Observable * GetObservable( pair<string,int>* );	//  Return wanted parameter & cache lookup ref
		Observable * GetObservable( string const );
		Observable * GetObservable( ObservableRef& );

		//	Pseudo-Observable functions which return an Observable or create one using the given dependency and relation
		//	Useful as once th observable is created it's permenantly stored in the DataPoint and the result is cached
		//	This is intended to reduce the number of angular calculations by effectivley pre-processing the dataset
		Observable* GetPseudoObservable( ObservableRef&, string, double (*pseudoRelation)(vector<double>) );
		Observable* GetPseudoObservable( ObservableRef&, string, pair<double,double> (*pseudoRelation)(vector<double>) );

		bool SetObservable( string, Observable* );
		bool SetObservable( const string, const double, const double, const string, const bool=false, const int=-1);

		//	Wanted for sorting datapoints
		bool operator() ( pair<DataPoint , pair<string,int> >, pair<DataPoint , pair<string,int> > );

	private:
		vector<Observable> allObservables;
		vector<string> allNames;

		vector<string> allPseudoNames;
		vector<Observable> allPseudoObservables; 
};

#endif
