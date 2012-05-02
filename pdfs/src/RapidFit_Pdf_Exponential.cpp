/** @class RapidFitPdfExponential RapidFitPdfExponential.cxx
 *  
 *
 *  PDf for a simple exponential with a single Gaussian resolution convoluted in analytically
 *
 *  For this PDF the only PhysicsPArameter is "gamma" , the exponential decay width
 *
 *  For this PDF the only relevant Observable is the decay time "time"
 *
 * The additional detector parameter is the resolution width
 *
 *  @author Pete clarke
 *  @date   2009-05-01
 */



#include "RapidFit_Pdf_Exponential.h"
#include <TMath.h>
#include <cmath>
#include <iostream>

PDF_CREATOR( RapidFitPdfExponential );

using namespace std;

//.................................
// Default constructor
RapidFitPdfExponential::RapidFitPdfExponential( PDFConfigurator* configurator ) : gamma(), resolution(0), valid(true), gammaName("gamma"), timeName("time")
{
	(void) configurator;
}

//	NO LONGER ACCESSIBLE 2011-10
//
//.................................
// Constructor with resolution
//RapidFitPdfExponential::RapidFitPdfExponential( double res ) : gamma(), resolution(), valid(), gammaName("gamma"), timeName("time")
//{
//  if( res >= 0.0 ) {
//    resolution = res ; 
//    valid = true ;
//  }
//  else {
//    resolution = 0.0; 
//    valid = false ;
//  }
//}


//....................................
//Copy constructor  
RapidFitPdfExponential::RapidFitPdfExponential( const RapidFitPdfExponential & other  ) :
   BasePDF( (BasePDF) other ), gamma( other.gamma ), resolution(other.resolution), valid(other.valid), gammaName("gamma"), timeName("time")
{  
}
 

//.................................
// Status
bool RapidFitPdfExponential::IsValid()
{ return valid ; }

//.................................
// Set Physics Parameters
//
// This PDF needs only one PhysicsParameter which is the decay width "Gamma"
// It must look for this in the set of PhysicsParameters which have been passed to it.
// If it cant find it it the method returns FALSE
// If it finds it OK it updates its internal copy of "Gamma"
// It only does this if the "Gamma" it has been passed is valid.

bool RapidFitPdfExponential::SetPhysicsParameters( ParameterSet * params )
{
	PhysicsParameter * newGamma = params->GetPhysicsParameter( gammaName );
	if ( newGamma->GetUnit() == "NameNotFoundError" )
	{
		cerr << "Parameter gamma not provided" << endl;
		return false;
	}
	else
	{
		gamma = newGamma;
		return true;
	}
  /*
  if( ! params->has("Gamma") ) { return false ; }
  
  if( params->get("Gamma").isValid() ) {
    gamma.set( params->get("Gamma") ) ;
	return true ;
  }
  else
    return false ; 
    */
}


//..................................
// Evaluate numerator of pdf
//
// This is called externally per DataPoint to evaluate the value of the PDF
// It must check that the relevant Observable "ProperTime"  is present in the DataPoint
// It shold also then check whether the Observable is OK for this PDF (e.g. in range)
// Once it has done all checks it passes to the internal method for calculation

double RapidFitPdfExponential::Evaluate( DataPoint * dataPoint )
{
	Observable * newTime = dataPoint->GetObservable( timeName );
	if ( newTime->GetUnit() == "NameNotFoundError" )
	{
		cerr << "Observable time not provided" << endl;
		return 0.0;
	}
	else
	{
		double t = newTime->GetValue();
		return this->evaluateNumerator( t ) ;  
	}
  /*
  if( ! dataPoint->has("ProperTime") ) { return -1.0 ; }

  double t = dataPoint->get("ProperTime").getValue() ; 
  
     	 
  return this.evaluateNumerator( t ) ;  
  */
}



//..................................
// Evaluate integral of pdf	
//
// This is called to evaluate the interal of the numerator over a specified set of bounds of
// the Observables
//
// It MUST return a negative number if it cannot do the integral.
//
// Typically it will only be called once when iterating over a DataSet because it 
// doesnt change for any given set of PhysicsParameters and DataPointBounds
// It does checks on the input and then passes to an internal method for calculation.
//

double RapidFitPdfExponential::Integral( DataPoint* dataPoint, PhaseSpaceBoundary * bounds )
{
	(void)dataPoint;
	IConstraint * timeBound = bounds->GetConstraint(timeName);
	if ( timeBound->GetUnit() == "NameNotFoundError" )
	{
		cerr << "Bound on time not provided" << endl;
		return 1.0;
	}
	else
	{
		double tlow = timeBound->GetMinimum();
		double thigh = timeBound->GetMaximum();
		return this->evaluateIntegral( tlow, thigh );
	}
  /*
  if( ! bounds->has("ProperTime") ) { return -1.0 ; }

  double tlow = bounds->get("ProperTime").getLow() ; 
  double thigh = bounds->get("ProperTime").getHigh() ; 
  
  if( thigh <  tlow ) 
     return -1.0 ;
  else 
    return this.evaluateIntegral( tlow, thigh ) ;  
  */
}


//...........................
// internal method to calculate numberator of PDF
//
// This is where you put all the maths

double  RapidFitPdfExponential::evaluateNumerator( double t ) const
{ 
   double gam = gamma->GetValue() ;

   if(resolution > 0.0 ) {     
      double theExp = exp( -gam*t + resolution*resolution * gam*gam / 2. ) ;
      double theErfc = TMath::Erfc(  -( t - resolution*resolution*gam ) /sqrt(2.)/resolution )  ;
      return theExp * theErfc  / 2.0 ;
   }	
   else
      if( t < 0.0 ) return 0.0 ;
   else
      return exp( -gam * t ) ;
  
}


//...........................
// internal method to calculate integral of PDF
//
// This is where you put all the maths
// Remember that this must return a negative number if it determines it cannot do the integral.

double RapidFitPdfExponential::evaluateIntegral( double tlow, double thigh) const 
{ 
   double gam = gamma->GetValue() ;

   if( resolution > 0.0 ) {
     // Cant do integral with resolution with finite limits
	 if( ( tlow <  - 5.0*resolution ) && (gam*thigh > 20.0 )  ) 
	   return ( 1.0 / gam) ;
	 else 
	   return -1.0 ;  
   }	  
   else
     if( tlow <= 0 )
       return ( 1.0 / gam  *  ( 1.0 - exp(thigh) )  ) ; 
     else
       return ( 1.0 / gam  *  ( exp(tlow) - exp(thigh) )  ) ; 
}

vector<string> RapidFitPdfExponential::GetPrototypeDataPoint()
{
	vector<string> dataPoint;
	dataPoint.push_back( "time" );
	return dataPoint;
}
vector<string> RapidFitPdfExponential::GetPrototypeParameterSet()
{
	vector<string> parameterSet;
	parameterSet.push_back( "gamma" );
	return parameterSet;
}

