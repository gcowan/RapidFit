// 
#ifndef RAPIDFITPDFEXPONENTIAL_H 
#define RAPIDFITPDFEXPONENTIAL_H 1

// Include files

#ifndef __CINT__
#include "IPDF.h"
#endif
#ifdef __CINT__
#include "framework/include/IPDF.h"
#endif

//#include "DataPoint.h"
//#include "PhaseSpaceBoundary.h"
//#include "PhysicsParameter.h"
//#include "ParameterSet.h"


/** @class RapidFitPdfExponential RapidFitPdfExponential.h
 *  
 *  This is an example for how to write a PDF implementing IPDF
 *  Because it is an example it is very heavily commented.
 *
 *  This PDF is for a simple exponential with a single Gaussian resolution convoluted in analytically
 *
 *  @author Pete clarke
 *  @date   2009-05-01
 */
class RapidFitPdfExponential : public IPDF  {

public: 
   //......................................
   // Default Constructor  which sets resolution to zero
   RapidFitPdfExponential( ) ;
   
   //......................................
   // Constructor with resolution
   RapidFitPdfExponential( double res ) ;
   
   //......................................
   // Copy constructor  
   RapidFitPdfExponential( const RapidFitPdfExponential& other );
      
   //......................................
   //  Destructor   
   inline virtual ~RapidFitPdfExponential() { };
   
   //......................................
   //  This tells you if the object has been constructed correctly with valid arguments 
   bool IsValid();
   
   //...................................... 
   // IPDF methods  - see IPDF documentation
   
   double Evaluate( DataPoint * d );
   
   double Integral( DataPoint* d, PhaseSpaceBoundary * b );
   
   bool SetPhysicsParameters( ParameterSet * p ) ;

   virtual vector<string> GetPrototypeDataPoint();
   virtual vector<string> GetPrototypeParameterSet();

private:

   // Physics Parameters - set by setPhysicsParameters() 
   PhysicsParameter gamma ;
   
   //Detector parameters
   double resolution ;
      
   // Internal Methods to evaluate numerator and denominator of PDF
   double evaluateNumerator( double t ) const;
   double evaluateIntegral( double tlow, double thigh ) const;
   
   // Internal status
   bool valid ;

   ObservableRef gammaName,timeName;   
};

#endif // RAPIDFITPDFEXPONENTIAL_H
