/**
        @class LLscanResult

        Container for all results from a minimisation

        @author Pete Clarke
	@date 2010-11-22
*/

//	ROOT Headers
#include "TGraph.h"
#include "TFrame.h"
#include "TAxis.h"
//	RapidFit Headers
#include "LLscanResult.h"
//	System Headers
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <math.h>

#define DOUBLE_TOLERANCE 1E-6

//Default constructor
LLscanResult::LLscanResult() : parameterName(), llmin(), parameterValues(), llvalues(), llvalues_offset()
{
}

//Main constructor
LLscanResult::LLscanResult( string _parameterName, vector<double> _parameterValues, vector<double> _llvalues  ):
	parameterName(_parameterName) , llmin(),
	parameterValues(_parameterValues),
	llvalues(_llvalues), llvalues_offset()
{
	// Reality check
	if( llvalues.size() != parameterValues.size() ) 
	{
		cout << endl << "In LLscanResult.print() llvalues.size()" << llvalues.size() << "  != parameterValues.size() " <<  parameterValues.size() << endl ;
		exit(1);
	}
		
	// Initialise llmin to a "safe"  large value. This is the only way to do it 
	llmin = 0;
	for(unsigned int i=0; i < llvalues.size() ;++i ) { if( llvalues[i] > llmin ) llmin = llvalues[i] ; }

	// Now find the minimum value in llvalues and set llmin to hold it.
	for(unsigned int i=0; i < llvalues.size() ;++i )
	{
		if( (llvalues[i] < llmin) && ((llvalues[i] - LLSCAN_FIT_FAILURE_VALUE ) > DOUBLE_TOLERANCE ) ) {
			llmin = llvalues[i] ;
		}
	}
	
	//Create llvalues_offset (i.e. offset to the minimum )
	for(unsigned int i=0; i < llvalues.size(); ++i )
	{		
		if( fabs( llvalues[i] - LLSCAN_FIT_FAILURE_VALUE ) < DOUBLE_TOLERANCE ) {
			llvalues_offset.push_back( 0. ) ;
		}
		else
		{
			llvalues_offset.push_back( llvalues[i] - llmin ) ;
		}
	}
	
}

//Destructor
LLscanResult::~LLscanResult()
{
}

//General Print Method
void LLscanResult::print()
{
	cout << endl << "LL scan results for parameter: " << parameterName << endl ;

	for(unsigned int i =0; i < llvalues.size() ; ++i ) {
		cout << setprecision(3) << "  " << parameterValues[i] << "      "<< llvalues_offset[i] << endl ;
	}
}

//Return scanned parameter values
vector<double> LLscanResult::GetParameterValues() { return parameterValues ; }

//Return raw LLscan values (i.e. without offset central value
vector<double> LLscanResult::GetRawLLvalues() { return llvalues ; }

//Return LL values offset to zero at the central value point
vector<double> LLscanResult::GetLLvalues() { return llvalues_offset ; } 

//Return a graph of the llscan
TGraph * LLscanResult::GetGraph() 
{
	double*  pvs = new double[parameterValues.size()] ;
	double* llvs = new double[parameterValues.size()] ;
	double llmax = 0 ;	
	for(unsigned int i=0; i< parameterValues.size() ; ++i ){
		pvs[i] = parameterValues[i] ;
		llvs[i] = llvalues_offset[i] ;
		if( llvs[i] > llmax ) llmax = llvs[i] ;
	}	

	TGraph* gr = new TGraph( Int_t(parameterValues.size()), pvs, llvs ) ;
	//gr->SetTitle("LL Scan for Parameter xxx");	
	gr->SetMarkerStyle(1);
	gr->SetLineWidth(2);
	gr->SetMarkerColor(4);
	gr->SetLineColor(4);
	gr->GetYaxis()->SetRangeUser( 0., llmax*1.2 );
	gr->GetYaxis()->SetLimits( 0., llmax*1.2 );
	gr->SetMinimum( 0.0 );
	gr->SetMaximum( llmax*1.2 );
	gr->Draw("ALP");
	string title("LL Scan for Parameter ") ;
	title+=parameterName.c_str();
	gr->SetTitle(title.c_str());	
	gr->GetXaxis()->SetTitle(parameterName.c_str());

	return gr ;
}
