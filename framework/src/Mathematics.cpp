/**
	@namespace Mathematics

	Namespace holding some common mathematical functions that are used
	by many PDFs. For example, gaussian's convolved with trigonometric
	functions.

        @author Greig A Cowan greig.cowan@cern.ch
        @date 2009-12-22
*/

#include <iostream>
#include "math.h"
#include "TMath.h"
#include "RooMath.h"

namespace Mathematics 
{

//--------------------------- exp and exp*sin and exp*cos time functions -------------------------
// time functions for use in PDFs with resolution
	
	//...................................................................
	// All of these functions were taken from Yue Hongs code
	
	// Calculate exp(-u^2) cwerf(swt*c + i(u+c)), taking care of numerical instabilities
	//RooComplex TimeFunctionUtility::evalCerf(double swt, double u, double c) const {
	//  RooComplex z(swt*c,u+c);
	//  return (z.im()>-4.0) ? RooMath::FastComplexErrFunc(z)*exp(-u*u) : evalCerfApprox(swt,u,c) ;
	///}  DIDNT APPEAR TO BE USED
	
	// use the approximation: erf(z) = exp(-z*z)/(sqrt(pi)*z) to explicitly cancel the divergent exp(y*y) behaviour of
	// CWERF for z = x + i y with large negative y
	RooComplex evalCerfApprox(double swt, double u, double c) 
	{
		static double rootpi= sqrt(atan2(0.,-1.));
		RooComplex z(swt*c,u+c);  
		RooComplex zc(u+c,-swt*c);
		RooComplex zsq= z*z;
		RooComplex v= -zsq - u*u;
		return v.exp()*(-zsq.exp()/(zc*rootpi) + 1)*2 ; //why shoule be a 2 here?
		//   return v.exp()*(-zsq.exp()/(zc*rootpi) + 1);
	}


	// Calculate Re(exp(-u^2) cwerf(swt*c + i(u+c))), taking care of numerical instabilities
	double evalCerfRe(double swt, double u, double c)  {
		RooComplex z(swt*c,u+c);
		return (z.im()>-4.0) ? RooMath::FastComplexErrFuncRe(z)*TMath::Exp(-u*u) : evalCerfApprox(swt,u,c).re() ;
	}
	
	// Calculate Im(exp(-u^2) cwerf(swt*c + i(u+c))), taking care of numerical instabilities
	double evalCerfIm(double swt, double u, double c)  {
		RooComplex z(swt*c,u+c);
		return (z.im()>-4.0) ? RooMath::FastComplexErrFuncIm(z)*TMath::Exp(-u*u) : evalCerfApprox(swt,u,c).im() ;
	}
			
	
	//........................................
	//evaluate a simple exponential with single gaussian time resolution
	
	double Exp( double t, double gamma, double resolution )  
	{
		
		if(resolution > 0.) {     
			
			double theExp = TMath::Exp( -t*gamma + resolution*resolution * gamma*gamma / 2. ) ;
			double theErfc = RooMath::erfc(  -( t - resolution*resolution*gamma ) /sqrt(2.)/resolution )  ;
			return theExp * theErfc  / 2.0 ;
			
			// Yue hongs code
			//double c = gamma * resolution /sqrt(2.); 
			//double u = t / resolution / sqrt(2.);
			//return exp( c*c - gamma*t ) * RooMath::erfc(c-u) / 2.;
		}	
		else {
			if( t < 0.0 ) return 0.0 ;
			return TMath::Exp( -gamma * t ) ;
		}
		
	}
	
	
	//.......................................
	// Evaluate integral of a simple exponential with single gaussian time resolution
	
	double ExpInt( double tlow, double thigh, double gamma, double resolution  )  
	{		
		if( thigh < tlow ) {
			std::cerr << " Mathematics::ExpInt: thigh is < tlow " << std::endl ;
			return -1.0 ;				
		}
		
		if( resolution > 0. ) {
			// This is a placeholder as I havnt put the correct code in yet as i dont know it.
			// So it only works if time limits are large and start from < 0
			if( ( tlow > -5.0*resolution ) || ( thigh < 5. ) ) {
				std::cerr << " Mathematics::ExpInt: cannot handle tlow > -"<<5.0*resolution<<" or thigh < 5  with resolution on" << std::endl ;
				return -1. ;				
			}
			return (1/gamma) * ( 1.0 - TMath::Exp(-gamma*thigh) ) ;
			
		}
		
		else {
			if( tlow < 0. ) return (1/gamma) * ( 1.0 - TMath::Exp(-gamma*thigh) ) ;
			else return (1/gamma) * ( TMath::Exp(-gamma*tlow) - TMath::Exp(-gamma*thigh) ) ;
		}
	}
	
	
	//.................................................................
	// Evaluate exponential X cosine with single time resolution
	
	double ExpCos( double t, double gamma, double deltaM, double resolution ) 
	{
		
		if(resolution > 0.) {     
			
			//Yue Hongs code
			double c = gamma * resolution/sqrt(2.);
			double u = t / resolution / sqrt(2.) ;
			double wt = deltaM / gamma ;
			return ( evalCerfRe(wt,-u,c) + evalCerfRe(-wt,-u,c) ) / 4. ;
			
			// My code which didnt work due to numerical instability
			//double theExp = exp( -t*gamma + timeRes*timeRes * ( gamma*gamma - deltaM*deltaM ) / 2. ) ;
			//double theCos = cos( deltaM * ( t - timeRes*timeRes*gamma ) ) ;
			//double theSin = sin( deltaM * ( t - timeRes*timeRes*gamma ) ) ;
			//RooComplex z( -( t - timeRes*timeRes*gamma )/sqrt(2.)/timeRes,  - timeRes*deltaM/sqrt(2.) ) ;
			//double theReErfc = (z.im()>-4.0) ? ( 1.0 - RooMath::FastComplexErrFuncRe(z) ) : ( 1.0 - RooMath::FastComplexErrFuncRe(z) );
			//double theImErfc = (z.im()>-4.0) ? ( 1.0 - RooMath::FastComplexErrFuncIm(z) ) : ( 1.0 - RooMath::FastComplexErrFuncIm(z) );
			//return theExp * ( theCos*theReErfc - theSin*theImErfc ) / 2.0 ;
			
		}	
		else {
			if( t < 0.0 ) return 0.0 ;
			return TMath::Exp( -gamma *t ) * cos( deltaM * t )  ;
		}
		
	}
	
	
	//.................................................................
	// Evaluate integral of exponential X cosine with single time resolution
	
	
	double ExpCosInt( double tlow, double thigh, double gamma, double deltaM, double resolution  )  
	{	
		if( thigh < tlow ) {
			std::cerr << " Mathematics::ExpInt: thigh is < tlow " << std::endl ;
			return -1.0 ;				
		}

		if( resolution > 0. ) {
			// This is a placeholder as I havnt put the correct code in yet as i dont know it.
			// So it only works if time limits are large and start from < 0
			if( ( tlow > -5.0*resolution ) || ( thigh < 5. ) ) {
				std::cerr << " Mathematics::ExpInt: cannot handle tlow > -"<<5.0*resolution<<" or thigh < 5  with resolution on" << std::endl ;
				return -1. ;				
			}
		}
		
		if( tlow < 0. ) tlow = 0. ;
		
		return (1/(gamma*gamma + deltaM*deltaM)) * (
													( TMath::Exp(-gamma*tlow)* (gamma*cos(deltaM*tlow) - deltaM*sin(deltaM*tlow)))
													-( TMath::Exp(-gamma*thigh)* (gamma*cos(deltaM*thigh) - deltaM*sin(deltaM*thigh)))
													);
		
	}
	
	
	//.................................................................
	// Evaluate exponential X sine with single time resolution
	
	double ExpSin( double t, double gamma, double deltaM, double resolution )  
	{
		if(resolution > 0.) {     
			
			//Yue Hongs code
			double c = gamma * resolution/sqrt(2.);
			double u = t / resolution / sqrt(2.) ;
			double wt = deltaM / gamma ;
			return ( evalCerfIm(wt,-u,c) - evalCerfIm(-wt,-u,c) ) /4. ;
			
			// My code which didnt work due to numerical instability
			//double theExp = exp( -t*gamma() + resolution*resolution * ( gamma()*gamma() - delta_ms*delta_ms ) / 2. ) ;
			//double theCos = cos( delta_ms * ( t - resolution*resolution*gamma() ) ) ;
			//double theSin = sin( delta_ms * ( t - resolution*resolution*gamma() ) ) ;
			//RooComplex z( -( t - resolution*resolution*gamma() )/sqrt(2.)/resolution,  - resolution*delta_ms/sqrt(2.) ) ;
			//double theReErfc = (z.im()>-4.0) ? ( 1.0 - RooMath::FastComplexErrFuncRe(z) ) : ( 1.0 - RooMath::FastComplexErrFuncRe(z) );
			//double theImErfc = (z.im()>-4.0) ? ( 1.0 - RooMath::FastComplexErrFuncIm(z) ) : ( 1.0 - RooMath::FastComplexErrFuncIm(z) );
			//return theExp * ( theCos*theImErfc + theSin*theReErfc ) / 2.0 ;
			
			
		}	
		else {
			if( t < 0.0 ) return 0.0 ;
			return TMath::Exp( -gamma *t ) * sin( deltaM * t )  ;
		}
		
	}
	
	
	//.................................................................
	// Evaluate integral of exponential X cosine with single time resolution
	
	
	double ExpSinInt( double tlow, double thigh, double gamma, double deltaM, double resolution  )  
	{
		if( thigh < tlow ) {
			std::cerr << " Mathematics::ExpInt: thigh is < tlow " << std::endl ;
			return -1.0 ;				
		}
		
		if( resolution > 0. ) {
			// This is a placeholder as I havnt put the correct code in yet as i dont know it.
			// So it only works if time limits are large and start from < 0
			if( ( tlow > -5.0*resolution ) || ( thigh < 5. ) ) {
				std::cerr << " Mathematics::ExpInt: cannot handle tlow > -"<<5.0*resolution<<" or thigh < 5  with resolution on" << std::endl ;
				return -1.0 ;				
			}
		}
		
		if( tlow < 0. ) tlow = 0. ;
		
		return (1/(gamma*gamma + deltaM*deltaM)) * (
													( TMath::Exp(-gamma*tlow)* (gamma*sin(deltaM*tlow) + deltaM*cos(deltaM*tlow)))
													-( TMath::Exp(-gamma*thigh)* (gamma*sin(deltaM*thigh) + deltaM*cos(deltaM*thigh)))
													);
		
	}
	
	
	
	
	
}

