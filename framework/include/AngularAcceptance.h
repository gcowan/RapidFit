/*!
 * @class AngularAcceptance
 *
 * @brief A class for holding angular acceptance information.
 *
 * @author Pete Clarke
 * @data 2012-03-29
 */

#pragma once
#ifndef ANGULAR_ACCEPTANCE_H
#define ANGULAR_ACCEPTANCE_H

//	RapidFit Headers
//	System Headers
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include "TMath.h"
#include "TAxis.h"
#include "TFile.h"
#include "TH3D.h"

#include "DataPoint.h"
#include "Observable.h"

using namespace::std;


//=======================================
class AngularAcceptance
{
	public:

		/*!
		  * @class AngularAcceptance
		  *
		  * @breif Constructor Function
		  *
		  * @param fileName Name of the file that contains the histogram to be used for angular Acceptance and the weights in a vector<double> TBranch in a TTree
		  *
		  * @param useHelicityBasis True: Look for the helicity Histogram object
		  *                         False: Look for the transversity Histogram object
		  *
		  * @param IgnoreAcceptanceHisto if True the class will not attempt to open the acceptance Histogram for use in the numerator of the PDF
		  *
		  * @param quiet True: reduce the amount of cout stataments
		  *              False: dump lots of info to cout
		  */
		AngularAcceptance( string fileName, bool useHelicityBasis, bool IgnoreAcceptanceHisto=false, bool quiet=false ) ;
		~AngularAcceptance();
		AngularAcceptance( const AngularAcceptance& );

		// Methods for numerator of PDF to return acceptance factor
		inline double af1() const { return _af1 ; } 
		inline double af2() const { return _af2 ; }
		inline double af3() const { return _af3 ; }
		inline double af4() const { return _af4 ; }
		inline double af5() const { return _af5 ; }
		inline double af6() const { return _af6 ; }
		inline double af7() const { return _af7 ; }
		inline double af8() const { return _af8 ; }
		inline double af9() const { return _af9 ; }
		inline double af10() const { return _af10 ; }

		// To get the acceptance for a given angular bin
		double getValue( double cosPsi, double cosTheta, double phi ) const;
		double getValue( Observable* cosPsi, Observable* cosTheta, Observable* phi ) const;
		double getValue( DataPoint* input ) const;

		void Print() const;

		double GetAvgBinContent() const { return average_bin_content; };
		double GetNumEntries() const { return total_num_entries; };
		double GetXmin() const { return xmin; };
		double GetXmax() const { return xmax; };
		double GetYmin() const { return ymin; };
		double GetYmax() const { return ymax; };
		double GetZmin() const { return zmin; };
		double GetZmax() const { return zmax; };

		double GetZeroBins() const { return zeroBins; };
	private:
		AngularAcceptance& operator= ( const AngularAcceptance& );

		//	double stream(ifstream& stream) ;

		double _af1, _af2, _af3, _af4, _af5, _af6, _af7, _af8, _af9, _af10 ; 
		bool useFlatAngularAcceptance ;

		TH3D *histo;
		TAxis *xaxis, *yaxis, *zaxis;
		int nxbins, nybins, nzbins;
		double xmin, xmax, ymin, ymax, zmin, zmax, deltax, deltay, deltaz;
		double total_num_entries;
		double average_bin_content;

		string openFile( string fileName, bool quiet=false ) ;
		double processHistogram( bool quiet=false ) ;

		ObservableRef cosThetaName, cosPsiName, phiName, helcosthetaKName, helcosthetaLName, helphiName;
		bool _useHelicityBasis;

		double zeroBins;

		mutable int psi_num;
		mutable int globalbin;
		mutable double num_entries_bin;
		mutable double _acc;
		mutable int xbin;
		mutable int ybin;
		mutable int zbin;
		mutable Observable* ThetaObs;
		mutable Observable* PsiObs;
		mutable Observable* PhiObs;
};

#endif

