/**
        @class PhysicsParameter

	Header to class for Scan Parameters
	This class constructs a Scan Parameter which contains the information on:
	the Parameter 'Name' of the Scan
	the Scan 'Type'
	the Scan Limits 'Max, Min' or number of 'Sigma' from the minima
	the Scan Resolution, or the number of 'Points' to plot

        @author Rob Currie rcurrie@cern.ch
	@date 2011-02
*/

#ifndef SCANPARAM_H
#define SCANPARAM_H

#include <string>
#include <vector>

using namespace std;

class ScanParam
{
	public:
		ScanParam();
		ScanParam( vector<string>, vector<string>, vector<double>, vector<double>, vector<int>, vector<int> );
		ScanParam( string, string, double, double, int );
		ScanParam( string, string, int, int );
		ScanParam( string, string, int );
		ScanParam( string );
		~ScanParam();

		//  These should probably be implemented in the .cpp but I'm being VERY lazy

		bool HasName(){   return !name.empty();                                                           };
		string GetName() { if( !name.empty() ){   return name[0]; }    else return "";                    };
		void SetName(string new_val) { while( !name.empty() ){   name.pop_back(); };    name.push_back(    new_val ); };

		bool HasType(){   return !type.empty();                                                           };
		string GetType() { if( !type.empty() ){   return type[0]; }    else return "";                    };
		void SetType(string new_val) { while( !type.empty() ){   type.pop_back(); };    type.push_back(    new_val ); };

		bool HasMax(){    return !maximum.empty();                                                        };
		double GetMax() { if( !maximum.empty() ){ return maximum[0]; } else return 0;                     };
		void SetMax(double new_val) { while( !maximum.empty() ){ maximum.pop_back(); }; maximum.push_back( new_val ); };

		bool HasMin(){    return !minimum.empty();                                                        };
		double GetMin() { if( !minimum.empty() ){ return minimum[0]; } else return 0;                     };
		void SetMin(double new_val) { while( !minimum.empty() ){ minimum.pop_back(); }; minimum.push_back( new_val ); };

		bool HasSigma(){  return !sigma.empty();                                                          };
		int GetSigma() { if( !sigma.empty() ) {   return sigma[0]; }   else return 5;                     };
		void SetSigma(int new_val) { while( !sigma.empty() ){    sigma.pop_back(); };   sigma.push_back(   new_val ); };

		bool HasPoints(){ return !points.empty();                                                         };
		int GetPoints() { if( !points.empty() ){  return points[0]; }   else return 10;                    };
		void SetPoints(int new_val) { while( !points.empty() ){  points.pop_back(); };  points.push_back(  new_val ); };

		void print();

	private:
		vector<string> name;
		vector<string> type;
		vector<double> minimum;
		vector<double> maximum;
		vector<int> sigma;
		vector<int> points;

};

#endif
