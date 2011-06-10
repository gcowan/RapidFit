/**
 @class SlicedAcceptance
 
 A class for holding a sliced propertime acceptance
 
 @author Pete Clarke
 @data 2011-06-07
 */


#include "SlicedAcceptance.h"


//............................................
// Constructor for flat acceptance
SlicedAcceptance::SlicedAcceptance( double tl, double th ) :
  nullSlice(new AcceptanceSlice(0.,0.,0.)), tlow(tl), thigh(th), beta(0)
{
	
	//Reality checks
	if( tl > th ) { cout << "SlicedAcceptance::SlicedAcceptance :  tl > th  - exiting" << endl ; exit(1) ; }
	
	//Create single slice with acceptance=1.0
	slices.push_back( new AcceptanceSlice( tlow, thigh, 1.0 ) ) ;
		
	//....done.....
	
}


//............................................
// Constructor for simple upper time acceptance only
SlicedAcceptance::SlicedAcceptance( double tl, double th, double b ) :
 nullSlice(new AcceptanceSlice(0.,0.,0.)), tlow(tl), thigh(th), beta(b)
{

	//Reality checks
	if( tl > th ) { cout << "SlicedAcceptance::SlicedAcceptance :  tl > th  - exiting" << endl ; exit(1) ; }
	if( beta < 0. ) { cout << "SlicedAcceptance::SlicedAcceptance :  beta < 0  - exiting" << endl ; exit(1) ; }
	
	//First step is to create the Time slices.
	//In this constructor it is the default for upper time acceptance
	
	// The upper time acceptance formula is (1 - beta*time) so determine highest and lowest acceptance 
	double accUp = 1.0 - beta*tlow ;
	double accLow = 1.0 - beta*thigh ;
	
	
	//Create N more, equispaced 
	int N = 100 ;
	double dh = (accUp-accLow) / N ;
	double dt = (thigh - tlow) / N ;

	//First slice which goes underneath the whole slope at the accLow level +dh/2
	slices.push_back( new AcceptanceSlice( tlow, thigh, (accLow+(dh/2.0)) ) ) ;

	for( int is=1; is <= (N-1) ; is++ )    //The N-1 is important as the N one is included as the +dh/2 in first slice
	{
		double th = tlow + is*dt  ;
		slices.push_back( new AcceptanceSlice( tlow, th, dh ) ) ;
	}
		
	//....done.....
	
}


//............................................
// Return numerator for evaluate
double SlicedAcceptance::getValue( double t ) const
{
	double returnValue = 0 ;	
	for( unsigned int is = 0; is < slices.size() ; is++ ){
		if( (t>=slices[is]->tlow()) && (t<slices[is]->thigh()) ) returnValue += slices[is]->height() ;
	}	
	return returnValue ;
}

//............................................
// Return the number of slices
int SlicedAcceptance::numberOfSlices() const { return slices.size() ; }

//............................................
// Return a slice
AcceptanceSlice * SlicedAcceptance::getSlice( int s ) const
{
	if( (s<0) || (s>=(int)slices.size()) ) return nullSlice ;
	return slices[s] ;
}
