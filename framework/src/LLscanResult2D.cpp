/**
        @class LLscanResult2D

        @author Pete Clarke
	@date 2010-11-22
*/

//	ROOT Headers
#include "TH2.h"
#include "TFrame.h"
#include "TAxis.h"
#include "TCanvas.h"
//	RapidFit Headers
#include "LLscanResult2D.h"
//	System Headers
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <math.h>

#define DOUBLE_TOLERANCE 1E-6

//Default constructor
LLscanResult2D::LLscanResult2D() : parameterName (), parameterValues(), parameterName2 (), parameterValues2(), listOfLLscans()
{
}

//Main constructor
LLscanResult2D::LLscanResult2D( string _parameterName, vector<double> _parameterValues, string _parameterName2, vector<double> _parameterValues2, vector<LLscanResult*> _llscans  ):
	parameterName(_parameterName) ,
	parameterValues(_parameterValues),
	parameterName2(_parameterName2) ,
	parameterValues2(_parameterValues2),
	listOfLLscans(_llscans)
{
	if( parameterValues.size() != listOfLLscans.size() ) {
		cout << "Constructing LLscanResult2D: number parameters .ne. number of scans " << endl ;
		exit(1);
	}
	if( parameterValues.size() == 0 ) {
		cout << "Constructing LLscanResult2D: zero first parameter list  " << endl ;
		exit(1);
	}
	if( parameterValues2.size() == 0  ) {
		cout << "Constructing LLscanResult2D: zero size second parameter list " << endl ;
		exit(1);
	}
	
	if( parameterValues2.size() != (listOfLLscans[0]->GetRawLLvalues()).size() ) {
		cout << "Constructing LLscanResult2D: number parameters2 .ne. size of LLscan " << endl ;
		exit(1);
	}
	
/*
 //PELC WORK IN PROGRESS

 //Now determine the minimum ll value, and create the offset llvalues

	//This is a way to set any safe llmin value to start with
	double llmin = 0 ;
	for( i1=0; i1 < parameterValues.size() ; ++i1 ) {
		vector<double> llscan = listOfLLscans[i1]->GetRawLLvalues() ; 
		for(i2=0; i2 < parameterValues2.size() ; ++i2 ) {
			if( llscan[i2] != LLSCAN_FIT_FAILURE_VALUE ) {
				llmin = llscan[i2] ;
				break ;
			}
		}
	}
	
	for( i1=0; i1 < parameterValues.size() ; ++i1 ) {
		vector<double> llscan_raw = listOfLLscans[i1]->GetRawLLvalues() ; 
		vector<double> llscan_offset ;
		for(i2=0; i2 < parameterValues2.size() ; ++i2 ) {
			if( llscan_raw != LLSCAN_FIT_FAILURE_VALUE ) {
				llscan_offset.push_back( llscan_raw - llmin ) ;
			}
			else {
				llscan_offset.push_back( 0 ) ;
			}
		}
		????? listOfLLscans_offset.push_back(llscan_offset) ;
	}
	
	//Now need to 
*/		
	
	
}

//Destructor
LLscanResult2D::~LLscanResult2D()
{
}


//General Print Method
void LLscanResult2D::print()
{
	cout << endl << "2D LL scan results "  <<  endl ;
	for(unsigned int i =0; i < listOfLLscans.size() ; ++i ) {
		cout << setprecision(5) << endl << "Outer loop over parameter " << parameterName << "  fixed to " << parameterValues[i] << endl ;
		listOfLLscans[i]->print() ;
	}
}



// Return a 2D histogram
TH2D * LLscanResult2D::GetTH2D() 
{

	int numberOfPoints = int(parameterValues.size())*int(parameterValues2.size());

	double* pvx = new double[unsigned(numberOfPoints)];
	double* pvy = new double[unsigned(numberOfPoints)];
	double* pvz = new double[unsigned(numberOfPoints)];

	//Extract lists of x,y, z (z is LL value)
	int ind = 0 ;
	for(unsigned int ix=0; ix < parameterValues.size() ; ++ix ) {
		vector<double> LLvalues= listOfLLscans[ix]->GetRawLLvalues() ;
		for(unsigned int iy=0; iy < parameterValues2.size() ; ++iy ) {
			pvx[ind] = parameterValues[ix]  ;  
			pvy[ind] = parameterValues2[iy] ; 
			pvz[ind] = LLvalues[iy];
			++ind;
		}
	}
	cout << endl ;
	cout << setprecision(6) ;
	//Now find the minimum value in the grid, excluding the -9999. values which indicate a failed fit.
	double llminNew = 0 ;
	int ipmin = -1 ;
	for(int ip=0; ip < numberOfPoints ; ++ip ) { if( (pvz[ip] > llminNew) ) llminNew = pvz[ip] ;}  // This is only safe way to initialise llminNew to a large value
	for(int ip=0; ip < numberOfPoints ; ++ip ) { 
		if( (pvz[ip] < llminNew) && (fabs(pvz[ip] + 9999.)<DOUBLE_TOLERANCE) ) 
		{
			llminNew = pvz[ip] ;
			ipmin = ip ;
		}
	} 

	//Now adjust the LL values to be relative to the minimum
	for(int ip=0; ip < numberOfPoints ; ++ip ) {
		
		if( fabs( pvz[ip] + 9999. ) < DOUBLE_TOLERANCE )
		{
			pvz[ip]-=llminNew ;
		}
		else
		{
			pvz[ip] = -1. ;
		}
	}

	// Write out the adjusted values
	for(int ii=0; ii<numberOfPoints; ++ii) {
		if( fabs( ii - ipmin ) < DOUBLE_TOLERANCE ) cout << " >>>>>>>>>>>>>MINIMUM>>>>>>>>>>>>>>" <<endl;
		cout << " i:  " << ii << " x:  "<< pvx[ii] << " y:  "  << pvy[ii]  << " LL: " <<pvz[ii] << endl ;
		if( fabs( ii - ipmin ) < DOUBLE_TOLERANCE ) cout << " >>>>>>>>>>>>>MINIMUM>>>>>>>>>>>>>>" <<endl;
	}
 
	
	// Make some plots
	TCanvas * c0 = new TCanvas;	

	TGraph2D * gr = new TGraph2D(numberOfPoints, pvx, pvy, pvz ) ;	
	gr->SetTitle(" ");
	
	gr->Draw("colz") ;   //cont, //cont1, lego
	c0->Update() ;	
	c0->SaveAs( "llcontour-graph-colz.eps" ) ;

	gr->Draw("cont1") ;   //cont, //cont1, lego
	c0->Update() ;	
	c0->SaveAs( "llcontour-graph-cont1.eps" ) ;

	TH2D* hist = gr->GetHistogram() ;
	hist->SetTitle(" ");
	
	double contours[4] ;
	contours[0]=1.15;    //not 0.5,  -2log(0.317) = 2.2977 / This must = 2* Delta_log(LL) /  Hence Delta_log(LL) = 1.15
	contours[1]=2.30;    //not ?,  -2log(0.1) = 4.605 / This must = 2* Delta_log(LL) /  Hence Delta_log(LL) = 2.30
	contours[2]=3.00;    //not 2,  -2log(0.05) = 5.99 / This must = 2* Delta_log(LL) /  Hence Delta_log(LL) = 3.00
	contours[3]=4.61;    //not 4.5,-2log(0.01) = 9.21 / This must = 2* Delta_log(LL) /  Hence Delta_log(LL) = 4.61 
	hist->SetContour(4, contours );
	hist->Draw("cont1") ;
	c0->Update() ;	
	c0->SaveAs( "llcontour.eps" ) ;
	
	
	//SetTitle( "" )  to get rid of big title
	//
	
	/*
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
	 */
								 
	return hist ;
}

