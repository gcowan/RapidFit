/**
        @class DataPoint

        Holds all observables for a given event

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/


#ifndef DATA_POINT_H
#define DATA_POINT_H

#include <vector>
#include <string>
#include "Observable.h"

using namespace std;

class DataPoint
{
	public:
		DataPoint();
		DataPoint( vector<string> );
		~DataPoint();

		vector<string> GetAllNames();
		Observable * GetObservable(string);
		bool SetObservable( string, Observable* );
		bool SetObservable( string, double, double, string );

	private:
		vector<Observable> allObservables;
		vector<string> allNames;
};

#endif
