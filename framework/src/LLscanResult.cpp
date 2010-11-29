/**
        @class LLscanResult

        Container for all results from a minimisation

        @author Pete Clarke
	@date 2010-11-22
*/

#include <iostream>
#include <iomanip>

#include "LLscanResult.h"
#include <TGraph.h>
#include <TFrame.h>
#include <TAxis.h>


//Default constructor
LLscanResult::LLscanResult()
{
}

//Main constructor
LLscanResult::LLscanResult( string _parameterName, double _centralParameterValue, double _parameterError, double _llmin, vector<double> _parameterValues, vector<double> _llvalues  ):
	parameterName(_parameterName) ,
	centralParameterValue( _centralParameterValue ),
	parameterError( _parameterError ),
	llmin(_llmin ),
	parameterValues(_parameterValues),
	llvalues(_llvalues)
{
}

//Destructor
LLscanResult::~LLscanResult()
{
}


//General Print Method
void LLscanResult::print()
{
	cout << endl << "LL scan results for parameter: " << parameterName <<  "   Central param value = " << centralParameterValue <<  endl ;
	for( int i =0; i < llvalues.size() ; i++ ) {
		cout << setprecision(3) << "  " << parameterValues[i] << "      "<< llvalues[i]-llmin << endl ;
	}
}

//Return scanned parameter values
vector<double> LLscanResult::GetParameterValues() { return parameterValues ; }

//Return LL values offset to zero
vector<double> LLscanResult::GetLLvalues() { 
	vector<double> offsetLL ;
	for(int ii=0; ii < llvalues.size(); ii++ )
	{
		offsetLL.push_back( llvalues[ii] - llmin ) ;
	}
	return offsetLL ; 
}



TGraph * LLscanResult::GetGraph() 
{
	
	double pvs[parameterValues.size()] ;
	double llvs[parameterValues.size()] ;
	double llmax = 0 ;	
	for( int i=0; i< parameterValues.size() ; i++ ){
		pvs[i] = parameterValues[i] ;
		llvs[i] =  llvalues[i] - llmin ;
		if( llvs[i] > llmax ) llmax = llvs[i] ;
	}	

	TGraph* gr = new TGraph( parameterValues.size(), pvs, llvs ) ;
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
