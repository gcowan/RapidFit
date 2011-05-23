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
		Observable const* GetSafeObservable( string const ) const;
		bool SetObservable( string, Observable* );
		bool SetObservable( const string, const double, const double, const string, const bool=false, const int=-1);

		//	Wanted for sorting datapoints
		bool operator() ( pair<DataPoint , pair<string,int> >, pair<DataPoint , pair<string,int> > );

	private:
		vector<Observable> allObservables;
		vector<string> allNames;
};

#endif
