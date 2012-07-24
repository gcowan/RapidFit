/**
  @namespace Mathematics

  Namespace holding some common mathematical functions that are used
  by many PDFs. For example, gaussian's convolved with trigonometric
  functions.

  @author Greig A Cowan greig.cowan@cern.ch
  @date 2009-12-22
 */

///	ROOT Headers
#include "TMath.h"
#include "RooMath.h"
///	RapidFit Headers
#include "main.h"
#include "Mathematics.h"
#include "Bd2JpsiKstar_sWave.h"
///	System Headers
#include <vector>
#include <cmath>
#include <iostream>
#include <pthread.h>
#include <iomanip>

pthread_mutex_t ROOT_Lock;

using namespace::std;

namespace Mathematics
{

	// Mathematica integral of the exp * erf
	//Integrate[(1*Exp[-(x/t) + s^2/(2*t^2)]* erfc[-((x - s^2/t)/(Sqrt[2]*s))])/2, x] ==
	//(t*(erf[x/(Sqrt[2]*s)] - E^((s^2 - 2*t*x)/(2*t^2))* erfc[(s^2 - t*x)/(Sqrt[2]*s*t)]))/2
	double expErfInt( double tlimit, double tau, double sigma)
	{
		const double sigma_2=sigma*sigma;
		const double inv_r2_sigma = 1./(sqrt_2*sigma);
		const double inv_r2_sigma_tau = inv_r2_sigma/tau;
		const double tau_tlimit = tau*tlimit;
		return 0.5 * (tau * ( erf( tlimit*inv_r2_sigma )
					- exp( (sigma_2*0.5 - tau_tlimit)/(tau*tau) )
					* erfc( (sigma_2 - tau_tlimit)*inv_r2_sigma_tau )
				    )
			     );
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
	RooComplex evalCerfApprox( double swt, double u, double c)
	{
		const double swt_c=swt*c;
		RooComplex z(swt_c,u+c);
		RooComplex zc(u+c,-swt_c);
		RooComplex zsq= z*z;
		RooComplex v= -zsq - u*u;
		return v.exp()*(-zsq.exp()/(zc*rootpi) + 1)*2 ; //why shoule be a 2 here?
		//   return v.exp()*(-zsq.exp()/(zc*rootpi) + 1);
	}

	// Calculate Re(exp(-u^2) cwerf(swt*c + i(u+c))), taking care of numerical instabilities
	double evalCerfRe( double swt, double u, double c)  {
		RooComplex* z = new RooComplex( swt*c, u+c );
		double returnable = 0.;
		if( z->im() > -4.0 )
		{
			pthread_mutex_lock( &ROOT_Lock );
			returnable = math_object->ITPComplexErrFuncRe( *z, 10 )*exp( -u*u );
			pthread_mutex_unlock( &ROOT_Lock );
		}
		else
		{
			returnable = evalCerfApprox( swt, u, c ).re();
		}
		delete z;
		return returnable;
	}

	// Calculate Im(exp(-u^2) cwerf(swt*c + i(u+c))), taking care of numerical instabilities
	double evalCerfIm( double swt, double u, double c)  {
		RooComplex* z = new RooComplex( swt*c, u+c );
		double returnable = 0.;
		if( z->im() > -4.0 )
		{
			pthread_mutex_lock( &ROOT_Lock );
			returnable = math_object->ITPComplexErrFuncIm( *z, 10 )*exp( -u*u );
			pthread_mutex_unlock( &ROOT_Lock );
		}
		else
		{
			returnable = evalCerfApprox( swt, u, c ).im();
		}
		delete z;
		return returnable;
	}

	//----------------------------------------------------------------------------------------------
	//........................................
	//evaluate a simple exponential with single gaussian time resolution
	double Exp( double t, double gamma, double resolution )
	{
		if(resolution > 0.) {

			const double t_gamma=t*gamma;
			const double resolution_2_gamma=resolution*resolution*gamma;
			const double theExp = exp( -t_gamma + resolution_2_gamma*gamma *0.5 ) ;
			const double theErfc = erfc(  -( t - resolution_2_gamma ) *_over_sqrt_2 /resolution )  ;
			return theExp * theErfc  *0.5 ;

			// Yue hongs code
			//double c = gamma * resolution *_over_sqrt_2;
			//double u = t / resolution *_over_sqrt_2;
			//return exp( c*c - gamma*t ) * erfc(c-u) / 2.;
		}
		else {
			if( t < 0.0 ) return 0.0 ;
			return exp( -t*gamma ) ;
		}

	}

	//.......................................
	// Evaluate integral of a simple exponential with single gaussian time resolution
	/* From Mathematica we have:
	   R2[t_] := 0.5 * Exp[-t*gamma + timeRes*timeRes*gamma*gamma*0.5] * Erfc[-(t - timeRes*timeRes*gamma)/(Sqrt[2]*timeRes)]
	   CForm[FullSimplify[Simplify[Integrate[R2[t], {t, tlow, thigh}]]]]
	   (
	   -2*exp((gamma*gamma*timeRes*timeRes)*0.5.) -
	   exp( gamma*thigh)*Erfc(thigh/(sqrt(2)*timeRes)) +
	   exp((gamma*gamma*timeRes*timeRes)*0.5) 		      * Erfc((thigh - gamma*timeRes*timeRes)/(sqrt(2)*timeRes)) +
	   exp((gamma*(2*thigh + gamma*timeRes*timeRes - 2*tlow))*0.5) * Erfc((gamma*timeRes*timeRes - tlow )/(sqrt(2)*timeRes)) +
	   exp( gamma*thigh)*Erfc(tlow/(sqrt(2)*timeRes))
	   )/(2.*exp(gamma*thigh)*gamma)
	 */
	double ExpInt( double tlow, double thigh, double gamma, double resolution  )
	{
		if( thigh < tlow )
		{
			std::cerr << " Mathematics::ExpInt: thigh is < tlow " << std::endl ;
			return -1.0 ;
		}

		const double invgamma=1./gamma;
		if( resolution > 0. )
		{
			return expErfInt(thigh, invgamma, resolution) - expErfInt(tlow, invgamma, resolution);
		}
		else
		{
			const double exp_gamma_thigh=exp(-gamma*thigh);
			if( tlow < 0. ) return invgamma * ( 1.0 - exp_gamma_thigh ) ;
			else return invgamma * ( exp(-gamma*tlow) - exp_gamma_thigh ) ;
		}
	}

	//------------------------------------------------------------------------------------------
	//........................................
	//evaluate a simple exponential with single gaussian time resolution, allowing for an acceptance (1.0 - b*t) - formula from wolfram
	double Exp_betaAcceptance( double t, double gamma, double resolution, double acceptanceParameter  )
	{
		if(resolution > 0.)
		{
			// At present we dont know how to do this with resolution included
			cout <<" Mathematics::Exp_betaAcceptance - with (1-bt) acceptance : This doesnt work when resolution .ne. 0 yet " << endl ;
			exit(1) ;
		}
		else {
			if( t < 0.0 ) return 0.0 ;
			return exp( -gamma * t ) * (1. - acceptanceParameter * t ) ;
		}
	}


	//.....................................................
	// Evaluate integral of a simple exponential with single gaussian time resolution, allowing for an acceptance (1.0 - b*t) - formula from wolfram
	// I = -1/G * exp (-tG) (1 - b(1/G +t ))  instead of I = -1/G * exp (-tG)
	double ExpInt_betaAcceptance( double tlow, double thigh, double gamma, double resolution, double acceptanceParameter  )
	{
		if( thigh < tlow )
		{
			std::cerr << " Mathematics::ExpInt: thigh is < tlow " << std::endl ;
			exit(1) ;
		}

		if( resolution > 0. )
		{
			// At present we dont know how to do this with resolution included
			cout <<" Mathematics::ExpInt_betaAcceptance - with (1-bt) acceptance : This doesnt work when resolution .ne. 0 yet " << endl ;
			exit(1) ;
		}
		else
		{
			double LoFactor=0, UpFactor=0 ;

			const double invgamma = 1./gamma;
			if( tlow < 0. ) {
				LoFactor = (1. - acceptanceParameter*invgamma) * (-invgamma) ;
			}
			else {
				LoFactor = (1. - acceptanceParameter*(invgamma+tlow)) * (-invgamma) * exp(-gamma*tlow) ;
			}
			UpFactor = (1. - acceptanceParameter*(invgamma+thigh)) * (-invgamma) * exp(-gamma*thigh) ;

			return UpFactor - LoFactor ;

		}
	}

	//----------------------------------------------------------------------------------------------------------
	//........................................
	//evaluate a simple exponential X cosh with single gaussian time resolution
	//When you express the cosh as a sum of exp and then multiply out, you are
	//left with just a sum of exponentials.
	double ExpCosh( double t, double gamma, double deltaGamma, double resolution )
	{
		const double dg_2 = deltaGamma*0.5;
		const double gammaL = gamma + dg_2;
		const double gammaH = gamma - dg_2;
		return ( Exp( t, gammaH, resolution ) + Exp( t, gammaL, resolution ) ) *0.5;
	}

	//.......................................
	// Evaluate integral of a simple exponential X cosh with single gaussian time resolution
	double ExpCoshInt( double tlow, double thigh, double gamma, double deltaGamma, double resolution  )
	{
		const double dg_2 = deltaGamma*0.5;
		const double gammaL = gamma + dg_2;
		const double gammaH = gamma - dg_2;
		return ( ExpInt( tlow, thigh, gammaH, resolution ) + ExpInt( tlow, thigh, gammaL, resolution ) )*0.5;
	}

	//........................................
	//evaluate a simple exponential with single gaussian time resolution
	//When you express the sinh as a sum of exp and then multiply out, you are
	//left with just a sum of exponentials.
	double ExpSinh( double t, double gamma, double deltaGamma, double resolution )
	{
		const double dg_2 = deltaGamma*0.5;
		const double gammaL = gamma + dg_2;
		const double gammaH = gamma - dg_2;
		return ( Exp( t, gammaH, resolution ) - Exp( t, gammaL, resolution ) )*0.5;
	}

	//.......................................
	// Evaluate integral of a simple exponential X sinh with single gaussian time resolution
	double ExpSinhInt( double tlow, double thigh, double gamma, double deltaGamma, double resolution  )
	{
		const double dg_2 = deltaGamma*0.5;
		const double gammaL = gamma + dg_2;
		const double gammaH = gamma - dg_2;
		return ( ExpInt( tlow, thigh, gammaH, resolution ) - ExpInt( tlow, thigh, gammaL, resolution ))*0.5;
	}

	// Mathematica integral of the exp * cos * erf
	//Integrate[(1*Exp[-(x/t) + s^2/(2*t^2)]* Erfc[-((x - s^2/t)/(Sqrt[2]*s))])/2, x] ==

	//.................................................................
	// Evaluate exponential X cosine with single time resolution
	double ExpCos( double t, double gamma, double deltaM, double resolution )
	{

		if(resolution > 0.) {

			//Yue Hongs code
			const double c = gamma * resolution*_over_sqrt_2;
			const double u = (t / resolution) *_over_sqrt_2 ;
			const double wt = deltaM / gamma ;
			return ( evalCerfRe(wt,-u,c) + evalCerfRe(-wt,-u,c) ) *0.25 ;

			// My code which didnt work due to numerical instability
			//double theExp = exp( -t*gamma + timeRes*timeRes * ( gamma*gamma - deltaM*deltaM ) / 2. ) ;
			//double theCos = cos( deltaM * ( t - timeRes*timeRes*gamma ) ) ;
			//double theSin = sin( deltaM * ( t - timeRes*timeRes*gamma ) ) ;
			//RooComplex z( -( t - timeRes*timeRes*gamma )*_over_sqrt_2/timeRes,  - timeRes*deltaM*_over_sqrt_2 ) ;
			//double theReErfc = (z.im()>-4.0) ? ( 1.0 - RooMath::FastComplexErrFuncRe(z) ) : ( 1.0 - RooMath::FastComplexErrFuncRe(z) );
			//double theImErfc = (z.im()>-4.0) ? ( 1.0 - RooMath::FastComplexErrFuncIm(z) ) : ( 1.0 - RooMath::FastComplexErrFuncIm(z) );
			//return theExp * ( theCos*theReErfc - theSin*theImErfc ) / 2.0 ;

		}
		else {
			if( t < 0.0 ) return 0.0 ;
			return exp( -gamma *t ) * cos( deltaM * t )  ;
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

		if( resolution > 0.) {
			//Added by Pete after getting code from Yuehong 120118
			const double c = gamma * resolution * _over_sqrt_2;
			const double umax = (thigh / resolution) *_over_sqrt_2 ;	
			const double umin = (tlow / resolution) *_over_sqrt_2 ;	
			const double wt = deltaM / gamma ;
			RooComplex evalDif(evalCerfApprox(-wt,-umax,c) - evalCerfApprox(-wt,-umin,c)) ;			
			double deltaCos = -0.5/gamma/(1+wt*wt) * ( evalDif.re() + wt*evalDif.im() + RooMath::erf(-umax) - RooMath::erf(-umin) ) ;
			return deltaCos ;
		}

		double real_tlow=tlow;
		if( tlow < 0. ) real_tlow = 0. ;

		const double deltaM_tlow=deltaM*real_tlow;
		const double deltaM_thigh=deltaM*thigh;
		return (
				( exp(-gamma*real_tlow)* (gamma*cos(deltaM_tlow) - deltaM*sin(deltaM_tlow)))
				-( exp(-gamma*thigh)* (gamma*cos(deltaM_thigh) - deltaM*sin(deltaM_thigh)))
		       )/(gamma*gamma + deltaM*deltaM);

	}


	//.................................................................
	// Evaluate exponential X sine with single time resolution
	double ExpSin( double t, double gamma, double deltaM, double resolution )
	{
		if(resolution > 0.) {

			//Yue Hongs code
			const double c = gamma * resolution*_over_sqrt_2;
			const double u = (t / resolution) *_over_sqrt_2 ;
			const double wt = deltaM / gamma ;
			const double val1 = evalCerfIm(wt,-u,c);
			const double val2 = evalCerfIm(-wt,-u,c);
			return ( val1 - val2 ) * 0.25 ;

			// My code which didnt work due to numerical instability
			//double theExp = exp( -t*gamma() + resolution*resolution * ( gamma()*gamma() - delta_ms*delta_ms ) / 2. ) ;
			//double theCos = cos( delta_ms * ( t - resolution*resolution*gamma() ) ) ;
			//double theSin = sin( delta_ms * ( t - resolution*resolution*gamma() ) ) ;
			//RooComplex z( -( t - resolution*resolution*gamma() )*_over_sqrt_2/resolution,  - resolution*delta_ms*_over_sqrt_2) ;
			//double theReErfc = (z.im()>-4.0) ? ( 1.0 - RooMath::FastComplexErrFuncRe(z) ) : ( 1.0 - RooMath::FastComplexErrFuncRe(z) );
			//double theImErfc = (z.im()>-4.0) ? ( 1.0 - RooMath::FastComplexErrFuncIm(z) ) : ( 1.0 - RooMath::FastComplexErrFuncIm(z) );
			//return theExp * ( theCos*theImErfc + theSin*theReErfc ) / 2.0 ;


		}
		else {
			if( t < 0.0 ) return 0.0 ;
			return exp( -gamma *t ) * sin( deltaM * t )  ;
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
			//Added by Pete after getting code from Yuehong 120118
			const double c = gamma * resolution * _over_sqrt_2;
			const double umax = (thigh / resolution) *_over_sqrt_2 ;	
			const double umin = (tlow / resolution) *_over_sqrt_2 ;	
			const double wt = deltaM / gamma ;
			RooComplex evalDif(evalCerfApprox(-wt,-umax,c) - evalCerfApprox(-wt,-umin,c)) ;			
			double deltaSin = -0.5/gamma/(1.+wt*wt) * ( -evalDif.im() +   wt*evalDif.re() -   -wt*(RooMath::erf(-umax) - RooMath::erf(-umin)) ) ;
			return deltaSin ;
		}

		double real_tlow=tlow;
		if( tlow < 0. ) real_tlow = 0. ;

		double deltaM_tlow=deltaM*real_tlow;
		double deltaM_thigh=deltaM*thigh;
		return (	( exp(-gamma*real_tlow)* (gamma*sin(deltaM_tlow) + deltaM*cos(deltaM_tlow)))
				-( exp(-gamma*thigh)* (gamma*sin(deltaM_thigh) + deltaM*cos(deltaM_thigh)))
		       )/(gamma*gamma + deltaM*deltaM);
	}

	//......................................................................
	void getBs2JpsiPhiAngularFunctions( double & f1
			, double & f2
			, double & f3
			, double & f4
			, double & f5
			, double & f6
			, const double cosTheta
			, const double phi
			, const double cosPsi)
	{
		const double sinTheta  = sqrt(1. - cosTheta*cosTheta);
		const double sinPsi    = sqrt(1. - cosPsi*cosPsi);

		const double cosPhi    = cos(phi);
		const double sinPhi    = sin(phi);

		const double sin2Theta = 2.*sinTheta*cosTheta;
		const double sin2Psi   = 2.*sinPsi*cosPsi;
		const double sin2Phi   = 2.*sinPhi*cosPhi;

		const double norm = global_frac;//9./32./TMath::Pi();
		const double sinTheta_sinTheta=sinTheta*sinTheta;
		const double sinPsi_sinPsi=sinPsi*sinPsi;

		// This is the factor that drops out when you integrate over the angles
		// Multiply by it here to get overall normalisation of 1.
		// i.e., int(2*cospsi*cospsi*(1-(1-costh*costh)*cos(phi)*cos(phi)),cospsi=-1..1, costh=-1..1,phi=-Pi..Pi);
		// same factor for f1, f2, f3. The remaining terms f4, f5, f6 give 0
		// int(-(1-cospsi*cospsi)*2*sqrt(1-costh*costh)*costh*sin(phi),cospsi=-1..1, costh=-1..1,phi=-Pi..Pi);
		f1 =  2.* cosPsi*cosPsi * ( 1. - sinTheta_sinTheta * cosPhi*cosPhi ) * norm;
		f2 =      sinPsi_sinPsi * ( 1. - sinTheta_sinTheta * sinPhi*sinPhi ) * norm;
		f3 =      sinPsi_sinPsi * sinTheta_sinTheta * norm;
		f4 = -1.* sinPsi_sinPsi * sin2Theta * sinPhi * norm;
		f5 = sin2Psi * sinTheta_sinTheta * sin2Phi*_over_sqrt_2 * norm;
		f6 = sin2Psi * sin2Theta * cosPhi*_over_sqrt_2 * norm;
		return;
	}

	void getBs2JpsiPhiAngularFunctionsWithSwave( double & f1
			, double & f2
			, double & f3
			, double & f4
			, double & f5
			, double & f6
			, double & f7
			, double & f8
			, double & f9
			, double & f10
			, const double cosTheta
			, const double phi
			, const double cosPsi)
	{
		const double sinTheta  = sqrt(1. - cosTheta*cosTheta);
		const double sinPsi    = sqrt(1. - cosPsi*cosPsi);

		const double cosPhi    = cos(phi);
		const double sinPhi    = sin(phi);

		const double sin2Theta = 2.*sinTheta*cosTheta;
		const double sin2Psi   = 2.*sinPsi*cosPsi;
		const double sin2Phi   = 2.*sinPhi*cosPhi;

		const double norm = global_frac;//9./32./TMath::Pi();
		const double sinTheta_sinTheta=sinTheta*sinTheta;
		const double sinPsi_sinPsi=sinPsi*sinPsi;
		const double cosPhi_cosPhi=cosPhi*cosPhi;

		f1 =  2.* cosPsi*cosPsi * ( 1. - sinTheta_sinTheta * cosPhi_cosPhi ) * norm;
		f2 =      sinPsi_sinPsi * ( 1. - sinTheta_sinTheta * sinPhi*sinPhi ) * norm;
		f3 =      sinPsi_sinPsi * sinTheta_sinTheta * norm;
		f4 = -1.* sinPsi_sinPsi * sin2Theta * sinPhi * norm;
		f5 = sin2Psi * sinTheta_sinTheta * sin2Phi*_over_sqrt_2 * norm;
		f6 = sin2Psi * sin2Theta * cosPhi*_over_sqrt_2 * norm;
		// Need to make sure that we deal with the normalisation of the following terms properly in the code that uses them
		f7 =  2*third * (1. - sinTheta_sinTheta * cosPhi_cosPhi) * norm ;			//Check: norm = (9./64./TMath::Pi())
		f8 =  third * root_6 * sinTheta_sinTheta * sin2Phi * sinPsi * norm;        		//Check: norm = 0
		f9 =  third * root_6 * sin2Theta * cosPhi * sinPsi * norm;         			//Check: norm = 0
		f10= (1+third) * root_3 * (1. - sinTheta_sinTheta * cosPhi_cosPhi) * cosPsi * norm;   	//Check: norm = 0
		return;
	}

	vector<double> calculateAcceptanceWeights( IDataSet * dataSet, IPDF * PDF )
	{
		vector<double> weights;
		// This tries to implement the NIKHEF method for calculating the
		// average acceptance weights for Bs2JpsiPhi.
		RapidFitIntegrator * rapidInt = new RapidFitIntegrator( PDF, true);
		PhaseSpaceBoundary * boundary = dataSet->GetBoundary();

		const int numAngularTerms = 10;//6;
		double*  f = new double[numAngularTerms]; // the angular functions
		double xi[numAngularTerms]; // the angular weights
		double cosTheta, phi, cosPsi, time; (void) time;
		double evalPDFraw, evalPDFnorm, val;
		int numEvents = dataSet->GetDataNumber();
		for (int i = 0; i < numAngularTerms; i++) xi[i] = 0.0;

		double Sum[numAngularTerms];
		double Sum_sq[numAngularTerms][numAngularTerms];
		double cov[numAngularTerms][numAngularTerms];
		double cor[numAngularTerms][numAngularTerms];

		for (int i = 0; i < numAngularTerms; i++) {
			Sum[i] = 0.0;
			for (int k = 0; k < numAngularTerms; k++) {
				cov[i][k] = 0.0;
				Sum_sq[i][k] = 0.0;
			}
		}

		double testEval = 0.;
		try
		{
			testEval = PDF->EvaluateTimeOnly( dataSet->GetDataPoint(0) );
		}
		catch(...)
		{
			testEval=-1.;
		}
		bool usePDF = testEval>=0. && !(isnan(testEval));

		for (int e = 0; e < numEvents; e++)
		{
			if (e % 1000 == 0) cout << "Event # " << e << "\t\t" << setprecision(4) << 100.*(double)e/(double)numEvents << "\% Complete\b\b\b\b\b\b\b\r\r\r\r\r\r\r\r\r\r\r";
			DataPoint * event = dataSet->GetDataPoint(e);
			cosTheta = event->GetObservable("cosTheta")->GetValue();
			phi      = event->GetObservable("phi")->GetValue();
			cosPsi   = event->GetObservable("cosPsi")->GetValue();
			time     = event->GetObservable("time")->GetValue();
			getBs2JpsiPhiAngularFunctionsWithSwave( f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9], cosTheta, phi, cosPsi);

			// The method sums the angular factor f_i divided by the sum_i(A_i*f_i)
			// for each accepted event. I'm not sure if dividing by Evaluate is exactly the
			// same here, particularly if you look at untagged events.
			evalPDFraw = PDF->Evaluate( event );
			// Now need to calculate the normalisation when integrated over the 3 angles
			vector<string> dontIntegrate = PDF->GetDoNotIntegrateList();
			dontIntegrate.push_back("time");
			//dontIntegrate.push_back("tag");
			dontIntegrate.push_back("tag");
			//evalPDFnorm = rapidInt->DoNumericalIntegral( event, boundary, dontIntegrate );
			if( usePDF )
			{
				evalPDFnorm = PDF->EvaluateTimeOnly( event );
			}
			else
			{
				evalPDFnorm = rapidInt->NumericallyIntegrateDataPoint( event, boundary, dontIntegrate );
			}
			val = evalPDFraw/evalPDFnorm;
			//cout << f[0] << " " << f[1]<< " " <<  f[2]<< " " <<  f[3]<< " " <<  f[4]<< " " <<  f[5]<< " " << f[6]<< " " <<  f[7]<< " " <<  f[8]<< " " <<  f[9]<< endl;
			//cout << time << " " << cosTheta  << " " << phi << " " << cosPsi << " " << evalPDFraw << " " << evalPDFnorm << " " << val << endl;
			for (int i = 0; i < numAngularTerms; i++)
			{
				Sum[i] = Sum[i] + f[i]/val;
				xi[i] += f[i]/val;

				for (int k = 0; k < numAngularTerms; k++)
				{
					Sum_sq[i][k] += f[i]/val*f[k]/val;
				}
			}
		}

		cout << endl;

		cout << "Covariance matrix " << endl;
		for (int i = 0; i < numAngularTerms; i++)
		{
			for (int k = 0; k < numAngularTerms; k++)
			{
				cov[i][k] = 1./numEvents/numEvents * ( Sum_sq[i][k] - Sum[i]*Sum[k]/numEvents);
				cout << cov[i][k] << "\t";
			}
			cout << endl;
		}

		cout << "Correlation matrix " << endl;
		for (int i = 0; i < numAngularTerms; i++)
		{
			for (int k = 0; k < numAngularTerms; k++)
			{
				cor[i][k] = cov[i][k]/(sqrt(cov[i][i])*sqrt(cov[k][k]));
				cout << cor[i][k] << "\t";
			}
			cout << endl;
		}

		cout << "Weight +- error " << endl;
		for (int i = 0; i < numAngularTerms; i++)
		{
			cout << fixed << setprecision(4) << xi[i]/numEvents << " \\pm " << sqrt(cov[i][i]) << endl;
			weights.push_back(xi[i]/numEvents);
		}	
		return weights;
	}

	void calculateAcceptanceWeightsWithSwave( IDataSet * dataSet, IPDF * PDF )
	{
		// This tries to implement the NIKHEF method for calculating the
		// average acceptance weights for Bs2JpsiPhi.
		RapidFitIntegrator * rapidInt = new RapidFitIntegrator( PDF, true);	(void)	rapidInt;	//	shutup gcc for unused param
		PhaseSpaceBoundary * boundary = dataSet->GetBoundary();
		int numAngularTerms = 10;
		double*  f = new double[(unsigned)numAngularTerms]; // the angular functions
		double* xi = new double[(unsigned)numAngularTerms]; // the angular weights
		double cosTheta, phi, cosPsi, time; (void) time;
		double evalPDFraw, evalPDFnorm, evalPDFnorm2, val;		(void) evalPDFnorm2;	//	shutup gcc for unused param!
		int numEvents = dataSet->GetDataNumber();
		for (int i = 0; i < numAngularTerms; i++) xi[i] = 0.0;

		for (int e = 0; e < numEvents; e++)
		{
			if (e % 10000 == 0) cout << "Event # " << e << endl;
			DataPoint * event = dataSet->GetDataPoint(e);
			cosTheta = event->GetObservable("cosTheta")->GetValue();
			phi      = event->GetObservable("phi")->GetValue();
			cosPsi   = event->GetObservable("cosPsi")->GetValue();
			time     = event->GetObservable("time")->GetValue();
			double weight = 1.;//event->GetObservable("NSig_sw")->GetValue();

			//cout << weight << endl;
			getBs2JpsiPhiAngularFunctionsWithSwave( f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9], cosTheta, phi, cosPsi);
			// The method sums the angular factor f_i divided by the sum_i(A_i*f_i)
			// for each accepted event. I'm not sure if dividing by Evaluate is exactly the
			// same here, particularly if you look at untagged events.

			evalPDFraw = PDF->Evaluate( event );
			// Now need to calculate the normalisation when integrated over the 3 angles
			vector<string> dontIntegrate = PDF->GetDoNotIntegrateList();
			dontIntegrate.push_back("time");
			dontIntegrate.push_back("tag");

			Bd2JpsiKstar_sWave * bpdf = (Bd2JpsiKstar_sWave *) PDF; //casting instance of IPDF as bpdf
			evalPDFnorm = bpdf->NormAnglesOnlyForAcceptanceWeights(event, boundary);
			//evalPDFnorm2 = rapidInt->DoNumericalIntegral( event, boundary, dontIntegrate );

			//	cout << "Normalisation using PDF method = " << evalPDFnorm << " : Normalisation using Numerical Integration = " << evalPDFnorm2 << endl;

			val = evalPDFraw  / evalPDFnorm;

			for (int i = 0; i < numAngularTerms; i++) //For each event work out 10 f values, and 10 xi values and find the average
			{
				xi[i] +=  weight * f[i] / val; ///1501544.3013043751;
			}
			//cout << f[0]<<" " << cosTheta << " " << phi << " " << cosPsi << " " << evalPDFraw  << " " << xi[0] << " " << evalPDFraw*weight << " " << val << endl;
		}


		float xi_scale = float( 3./ (xi[0] + xi[1] + xi[2]) );
		//cout << "[" << xi[0] << ", " << xi[1] << ", " << xi[2] <<  ", " << xi[3] << ", " << xi[4] << ", " << xi[5] <<  "," << xi[7] << "," << xi[8] << ","  << xi[9] <<"]" <<  endl;
		cout << "[" << xi[0]*xi_scale << ", " << xi[1]*xi_scale << ", " << xi[2]*xi_scale <<  ", " << xi[3]*xi_scale << ", " << xi[4]*xi_scale << ", " << xi[5]*xi_scale <<  ", " << xi[6]*xi_scale << "," << xi[7]*xi_scale << "," << xi[8]*xi_scale << ","  << xi[9]*xi_scale <<"]" <<  endl;
		//cout << "[" << xi[0]/numEvents << ", " << xi[1]/numEvents << ", " << xi[2]/numEvents <<  ", " << xi[3]/numEvents << ", " << xi[4]/numEvents << ", " << xi[5]/numEvents  << "," <<  xi[6]/numEvents << "," << xi[7]/numEvents << "," << xi[8]/numEvents << ","  << xi[9]/numEvents << "]" <<  endl;
		return;
	}

	double Exp_Wrapper( vector<double> input )
	{
		return Exp( input[0], input[1], input[2] );
	}

	double ExpInt_Wrapper( vector<double> input )
	{
		return ExpInt( input[0], input[1], input[2], input[3] );
	}

	double Exp_betaAcceptance_Wrapper( vector<double> input )
	{
		return Exp_betaAcceptance( input[0], input[1], input[2], input[3] );
	}

	double ExpInt_betaAcceptance_Wrapper( vector<double> input )
	{
		return ExpInt_betaAcceptance( input[0], input[1], input[2], input[3], input[4] );
	}

	double ExpCosh_Wrapper( vector<double> input )
	{
		return ExpCosh( input[0], input[1], input[2], input[3] );
	}

	double ExpCoshInt_Wrapper( vector<double> input )
	{
		return ExpCoshInt( input[0], input[1], input[2], input[3], input[4] );
	}

	double ExpSinh_Wrapper( vector<double> input )
	{
		return ExpSinh( input[0], input[1], input[2], input[3] );
	}

	double ExpSinhInt_Wrapper( vector<double> input )
	{
		return ExpSinhInt( input[0], input[1], input[2], input[3], input[4] );
	}

	double ExpCos_Wrapper( vector<double> input )
	{
		return ExpCos( input[0], input[1], input[2], input[3] );
	}

	double ExpCosInt_Wrapper( vector<double> input )
	{
		return ExpCosInt( input[0], input[1], input[2], input[3], input[4] );
	}

	double ExpSin_Wrapper( vector<double> input )
	{
		return ExpSin( input[0], input[1], input[2], input[3] );
	}

	double ExpSinInt_Wrapper( vector<double> input )
	{
		return ExpSinInt( input[0], input[1], input[2], input[3], input[4] );
	}

	double expErfInt_Wrapper( vector<double> input )
	{
		return expErfInt( input[0], input[1], input[2] );
	}

}


