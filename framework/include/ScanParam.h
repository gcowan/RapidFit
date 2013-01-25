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

#pragma once
#ifndef SCANPARAM_H
#define SCANPARAM_H

//	System Headers
#include <string>
#include <vector>

using namespace::std;

class ScanParam
{
	public:
		ScanParam( vector<string>, vector<double>, vector<double>, vector<int>, vector<int> );
		ScanParam( string, double, double, int );
		ScanParam( string, int, int );
		ScanParam( string, int );
		ScanParam( string );

		/*!
		 * @brief Destructor
		 */
		~ScanParam();

		//  These should probably be implemented in the .cpp but I'm being VERY lazy

		bool HasName();
		string GetName();
		void SetName(string new_val);

		bool HasMax();
		double GetMax();
		void SetMax(double new_val);

		bool HasMin();
		double GetMin();
		void SetMin(double new_val);

		bool HasSigma();
		int GetSigma();
		void SetSigma(int new_val);

		bool HasPoints();
		int GetPoints();
		void SetPoints(int new_val);

		void print();

	private:
		vector<string> name;
		vector<double> minimum;
		vector<double> maximum;
		vector<int> sigma;
		vector<int> points;

};

#endif

