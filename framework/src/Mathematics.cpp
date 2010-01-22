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
	// Mathematica integral of the exp * erf
	//Integrate[(1*Exp[-(x/t) + s^2/(2*t^2)]* Erfc[-((x - s^2/t)/(Sqrt[2]*s))])/2, x] ==
	//(t*(Erf[x/(Sqrt[2]*s)] - E^((s^2 - 2*t*x)/(2*t^2))* Erfc[(s^2 - t*x)/(Sqrt[2]*s*t)]))/2
	double expErfInt( double tlimit, double tau, double sigma)
	{
		double val = 0.5 * (tau * ( RooMath::erf( tlimit/(sqrt(2.)*sigma) )
					- exp( (sigma*sigma - 2*tau*tlimit)/(2.*tau*tau) )
					* RooMath::erfc( (sigma*sigma - tau*tlimit)/(sqrt(2.)*sigma*tau) )
					)
				);
		return val;
	}

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
	/* From Mathematica we have:
	   R2[t_] := 1/2 * Exp[-t*gamma + timeRes*timeRes*gamma*gamma/2] * Erfc[-(t - timeRes*timeRes*gamma)/(Sqrt[2]*timeRes)]
	   CForm[FullSimplify[Simplify[Integrate[R2[t], {t, tlow, thigh}]]]]
	   (
	   -2*exp((gamma*gamma*timeRes*timeRes)/2.) - 
	   exp( gamma*thigh)*Erfc(thigh/(sqrt(2)*timeRes)) + 
	   exp((gamma*gamma*timeRes*timeRes)/2.) 		      * Erfc((thigh - gamma*timeRes*timeRes)/(sqrt(2)*timeRes)) + 
	   exp((gamma*(2*thigh + gamma*timeRes*timeRes - 2*tlow))/2.) * Erfc((gamma*timeRes*timeRes - tlow )/(sqrt(2)*timeRes)) + 
	   exp( gamma*thigh)*Erfc(tlow/(sqrt(2)*timeRes))
	   )/(2.*exp(gamma*thigh)*gamma)

	   We should be able to use the above as a direct replacement of the code below
	 */
	double ExpInt( double tlow, double thigh, double gamma, double resolution  )  
	{		
		if( thigh < tlow ) 
		{
			std::cerr << " Mathematics::ExpInt: thigh is < tlow " << std::endl ;
			return -1.0 ;				
		}

		if( resolution > 0. ) 
		{
			// This is a placeholder as I havnt put the correct code in yet as i dont know it.
			// So it only works if time limits are large and start from < 0
			/*
			   if( ( tlow > -5.0*resolution ) || ( thigh < 5. ) ) {
			   std::cerr << " Mathematics::ExpInt: cannot handle tlow > -"<<5.0*resolution<<" or thigh < 5  with resolution on" << std::endl ;
			   return -1. ;				
			   }
			   return (1/gamma) * ( 1.0 - TMath::Exp(-gamma*thigh) ) ;
			 */
			return expErfInt(thigh, 1./gamma, resolution) - expErfInt(tlow, 1./gamma, resolution);

		}
		else
		{
			if( tlow < 0. ) return (1/gamma) * ( 1.0 - TMath::Exp(-gamma*thigh) ) ;
			else return (1/gamma) * ( TMath::Exp(-gamma*tlow) - TMath::Exp(-gamma*thigh) ) ;
		}
	}

	//........................................
	//evaluate a simple exponential X cosh with single gaussian time resolution
	//When you express the cosh as a sum of exp and then multiply out, you are
	//left with just a sum of exponentials.
	double ExpCosh( double t, double gamma, double deltaGamma, double resolution )
	{
		double gammaL = gamma + deltaGamma/2.;
		double gammaH = gamma - deltaGamma/2.;
		return ( Exp( t, gammaH, resolution ) + Exp( t, gammaL, resolution ) ) / 2.;
	}

	//.......................................
	// Evaluate integral of a simple exponential X cosh with single gaussian time resolution
	double ExpCoshInt( double tlow, double thigh, double gamma, double deltaGamma, double resolution  )
	{
		double gammaL = gamma + deltaGamma/2.;
		double gammaH = gamma - deltaGamma/2.;
		return ( ExpInt( tlow, thigh, gammaH, resolution ) + ExpInt( tlow, thigh, gammaL, resolution ) )/2.;
	}

	//........................................
	//evaluate a simple exponential with single gaussian time resolution
	//When you express the sinh as a sum of exp and then multiply out, you are
	//left with just a sum of exponentials.
	double ExpSinh( double t, double gamma, double deltaGamma, double resolution )
	{
		double gammaL = gamma + deltaGamma/2.;
		double gammaH = gamma - deltaGamma/2.;
		return ( Exp( t, gammaH, resolution ) - Exp( t, gammaL, resolution ) )/2.;
	}

	//.......................................
	// Evaluate integral of a simple exponential X sinh with single gaussian time resolution
	double ExpSinhInt( double tlow, double thigh, double gamma, double deltaGamma, double resolution  )
	{
		double gammaL = gamma + deltaGamma/2.;
		double gammaH = gamma - deltaGamma/2.;
		return ( ExpInt( tlow, thigh, gammaH, resolution ) - ExpInt( tlow, thigh, gammaL, resolution ))/2.;
	}

	// Mathematica integral of the exp * cos * erf
	//Integrate[(1*Exp[-(x/t) + s^2/(2*t^2)]* Erfc[-((x - s^2/t)/(Sqrt[2]*s))])/2, x] ==

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
				std::cerr << " Mathematics::ExpCosInt: cannot handle tlow > -"<<5.0*resolution<<" or thigh < 5  with resolution on" << std::endl ;
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
				std::cerr << " Mathematics::ExpSinInt: cannot handle tlow > -"<<5.0*resolution<<" or thigh < 5  with resolution on" << std::endl ;
				return -1.0 ;				
			}
		}

		if( tlow < 0. ) tlow = 0. ;

		return (1/(gamma*gamma + deltaM*deltaM)) * (
				( TMath::Exp(-gamma*tlow)* (gamma*sin(deltaM*tlow) + deltaM*cos(deltaM*tlow)))
				-( TMath::Exp(-gamma*thigh)* (gamma*sin(deltaM*thigh) + deltaM*cos(deltaM*thigh)))
				);
	}


	void getBs2JpsiPhiAngularFunctions( double & f1
			, double & f2
			, double & f3
			, double & f4
			, double & f5
			, double & f6
			, double cosTheta
			, double phi
			, double cosPsi)
	{
		double sinTheta  = sqrt(1. - cosTheta*cosTheta);
		double sinPsi    = sqrt(1. - cosPsi*cosPsi);

		double cosPhi    = cos(phi);
		double sinPhi    = sin(phi);

		double sin2Theta = 2.*sinTheta*cosTheta;
		double sin2Psi   = 2.*sinPsi*cosPsi;
		double sin2Phi   = 2.*sinPhi*cosPhi;

		//double norm = 9./32./TMath::Pi();
		//this is the factor that drops out when you integrate over the angles
		// i.e., int(2*cospsi*cospsi*(1-(1-costh*costh)*cos(phi)*cos(phi)),cospsi=-1..1, costh=-1..1,phi=-Pi..Pi);
		// same factor for f1, f2, f3. The remaining terms f4, f5, f6 give 0
		// int(-(1-cospsi*cospsi)*2*sqrt(1-costh*costh)*costh*sin(phi),cospsi=-1..1, costh=-1..1,phi=-Pi..Pi); 
		f1 =  2.* cosPsi*cosPsi * ( 1. - sinTheta*sinTheta * cosPhi*cosPhi );
		f2 =      sinPsi*sinPsi * ( 1. - sinTheta*sinTheta * sinPhi*sinPhi );
		f3 =      sinPsi*sinPsi * sinTheta*sinTheta;
		f4 = -1.* sinPsi*sinPsi * sin2Theta * sinPhi;
		f5 = sin2Psi * sinTheta*sinTheta * sin2Phi/sqrt(2.);
		f6 = sin2Psi * sin2Theta * cosPhi/sqrt(2.);
		return;
	}

}

