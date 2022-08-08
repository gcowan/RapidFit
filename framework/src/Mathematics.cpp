/**
  @namespace Mathematics

  Namespace holding some common mathematical functions that are used
  by many PDFs. For example, gaussian's convolved with trigonometric
  functions.

  @author Greig A Cowan greig.cowan@cern.ch
  @date 2009-12-22
 */

///	ROOT Headers
#include "TNtuple.h"
#include "TMath.h"
#include "TH1D.h"
#include "TF3.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "RooMath.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/GSLMCIntegrator.h"
///	RapidFit Headers
#include "main.h"
#include "Mathematics.h"
#include "Bd2JpsiKstar_sWave.h"
#include "Bd2JpsiKstar_sWave_Fscopy.h"
#include "Bd2JpsiKstar_sWave_Fs_withAcc.h"
#include "DPHelpers.hh"
#include "LegendreMomentShape.h"
///	System Headers
#include <vector>
#include <cmath>
#include <iostream>
#include <pthread.h>
#include <iomanip>
#include <complex>

pthread_mutex_t ROOT_Lock = pthread_mutex_t();

using namespace::std;

bool RooMathinit=false;

#ifdef __RAPIDFIT_USE_GSL
///	GSL Headers
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_qrng.h>

const gsl_complex prefix = gsl_complex_rect( (2./sqrt( Mathematics::Pi())), 0. );

gsl_complex gsl_erf( gsl_complex z )
{
	gsl_complex result = gsl_complex_rect( 0., 0. );
	gsl_complex numerator  = gsl_complex_rect( 0., 0.);
	gsl_complex denominator  = gsl_complex_rect( 0., 0. );

	//	Common factor which changes term to term
	int factor=-1;
	//gsl_complex gsl_factor = gsl_complex_rect( 0., 0. );

	//	Running numerator turn in the expansion
	gsl_complex z_raised = z;

	//	Common Factor in the expansion
	gsl_complex z_square = gsl_complex_mul( z, z );

	//	Running factor to govern the size of each step
	unsigned int fact =1;

	numerator = z;
	GSL_SET_REAL( &denominator, fact );
	result  = gsl_complex_add( result, gsl_complex_div( numerator, denominator ) );

	//	Can't get better than double precision
	//	If the numerator varies by less than this we should assume we've reached as good as we can get
	for( unsigned int i = 1;
			(GSL_REAL(z_raised)>(1E-25)*GSL_REAL(z)) && (GSL_IMAG(z_raised)>(1E-25)*GSL_IMAG(z));
			++i )
	{
		factor*=-1;
		//	Assumed Quicker than constructing a complex number and doing a complex multiplication
		GSL_SET_REAL( &z, (double)factor*GSL_REAL(z) );

		z_raised = gsl_complex_mul( z_raised, z_square );
		fact*=i;

		//numerator  = gsl_complex_pow_real( z, 2*i+1 );
		numerator = z_raised;

		//	Only a real number
		GSL_SET_REAL( &denominator, (2*i + 1)*fact );

		//	Add this term in the series to the result
		result  = gsl_complex_add( result, gsl_complex_div( numerator, denominator ) );
	}

	result = gsl_complex_mul( prefix, result );

	return result;
}

#endif

namespace Mathematics
{

	// Mathematica integral of the exp * erf
	//Integrate[(1*Exp[-(x/t) + s^2/(2*t^2)]* erfc[-((x - s^2/t)/(Sqrt[2]*s))])/2, x] ==
	//(t*(erf[x/(Sqrt[2]*s)] - E^((s^2 - 2*t*x)/(2*t^2))* erfc[(s^2 - t*x)/(Sqrt[2]*s*t)]))/2
	double expErfInt( double tlimit, double tau, double sigma)
	{
		const double sigma_2=sigma*sigma;
		const double tau_2 = tau*tau;
		const double inv_tau2 = 1./tau_2;
		const double inv_r2_sigma = 1./(sqrt_2*sigma);
		const double inv_r2_sigma_tau = inv_r2_sigma/tau;
		const double tau_tlimit = tau*tlimit;
		return 0.5 * (tau * ( erf( tlimit*inv_r2_sigma )
					- exp( (sigma_2*0.5 - tau_tlimit)*inv_tau2 )
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
	complex<double> evalCerfApprox( double swt, double u, double c)
	{
		const double swt_c=swt*c;
		complex<double> z(swt_c,u+c);
		complex<double> zc(u+c,-swt_c);
		complex<double> zsq= z*z;
		complex<double> v= -zsq - u*u;
		double v_mag = exp( v.real() );
		complex<double> v_exp( v_mag*cos(v.imag()), v_mag*sin(v.imag()) );
		double zsq_mag = exp( zsq.real() );
		complex<double> zsq_exp( zsq_mag*cos(zsq.imag()), zsq_mag*sin(zsq.imag()) );
		return v_exp*(-zsq_exp/(zc*rootpi) + 1.)*2.; //why shoule be a 2 here?
		//   return v.exp()*(-zsq.exp()/(zc*rootpi) + 1);
	}

	complex<double> evalCerf( double swt, double u, double c )
	{
#if ROOT_VERSION_CODE < ROOT_VERSION(5,34,10)
		RooComplex z( swt*c, u+c );
		complex<double> returnable = 0.;
		if( (u+c) > -4.0 )
		{
			if( !RooMathinit )
			{
				pthread_mutex_lock( &ROOT_Lock );
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,34,0)
				if( !RooMathinit ) RooMath::initFastCERF( 800, -4.0, 4.0, 1000, -4.0, 6.0 );
#else
				if( !RooMathinit ) RooComplex RooReturnable = RooMath::ITPComplexErrFunc( z, 13 )*exp( -u*u );
#endif
				RooMathinit = true;
				pthread_mutex_unlock( &ROOT_Lock );
			}
			RooComplex thisRes = RooMath::ITPComplexErrFunc( z, 13 )*exp( -u*u );
			returnable = complex<double>( thisRes.re(), thisRes.im() );
		}
		else
		{
			returnable = evalCerfApprox( swt, u, c );
		}
		return returnable;
#else
		complex<double> z_stl( swt*c, u+c );
		complex<double> returnable = 0.;
		if( (u+c) > -4.0 )
		{
			if( !RooMathinit )
			{
				pthread_mutex_lock( &ROOT_Lock );
				if( !RooMathinit ) returnable = RooMath::faddeeva( z_stl )*exp( -u*u );
				RooMathinit = true;
				pthread_mutex_unlock( &ROOT_Lock );
			}
			returnable = RooMath::faddeeva( z_stl )*exp( -u*u );
		}
		else
		{
			returnable = evalCerfApprox( swt, u, c );
		}
		return returnable;
#endif
	}

	double evalCerfRe( double swt, double u, double c )
	{
#if ROOT_VERSION_CODE < ROOT_VERSION(5,34,10)
		RooComplex z( swt*c, u+c );
		double returnable = 0.;
		if( (u+c) > -4.0 )
		{
			if( !RooMathinit )
			{
				pthread_mutex_lock( &ROOT_Lock );
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,34,0)
				if( !RooMathinit ) RooMath::initFastCERF( 800, -4.0, 4.0, 1000, -4.0, 6.0 );
#else
				if( !RooMathinit ) RooComplex RooReturnable = RooMath::ITPComplexErrFuncRe( z, 13 )*exp( -u*u );
#endif
				RooMathinit = true;
				pthread_mutex_unlock( &ROOT_Lock );
			}
			returnable = RooMath::ITPComplexErrFuncRe( z, 13 )*exp( -u*u );
		}
		else
		{
			returnable = evalCerfApprox( swt, u, c ).real();
		}
		return returnable;
#else
		complex<double> z_stl( swt*c, u+c );
		double returnable = 0.;
		if( (u+c) > -4.0 )
		{
			if( !RooMathinit )
			{
				pthread_mutex_lock( &ROOT_Lock );
				if( !RooMathinit ) returnable = RooMath::faddeeva( z_stl ).real()*exp( -u*u );
				RooMathinit = true;
				pthread_mutex_unlock( &ROOT_Lock );
			}
			returnable = RooMath::faddeeva( z_stl ).real()*exp( -u*u );
		}
		else
		{
			returnable = evalCerfApprox( swt, u, c ).real();
		}
		return returnable;
#endif
	}

	double evalCerfIm( double swt, double u, double c)
	{
#if ROOT_VERSION_CODE < ROOT_VERSION(5,34,10)
		RooComplex z( swt*c, u+c );
		double returnable = 0.;
		if( (u+c) > -4.0 )
		{
			if( !RooMathinit )
			{
				pthread_mutex_lock( &ROOT_Lock );
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,34,0)
				if( !RooMathinit ) RooMath::initFastCERF( 800, -4.0, 4.0, 1000, -4.0, 6.0 );
#else
				if( !RooMathinit ) returnable = RooMath::ITPComplexErrFuncIm( z, 13 )*exp( -u*u );
#endif
				RooMathinit = true;
				pthread_mutex_unlock( &ROOT_Lock );
			}
		}
		else
		{
			returnable = evalCerfApprox( swt, u, c ).imag();
		}
		return returnable;
#else
		complex<double> z_stl( swt*c, u+c );
		double returnable = 0.;
		if( (u+c) > -4.0 )
		{
			if( !RooMathinit )
			{
				pthread_mutex_lock( &ROOT_Lock );
				if( !RooMathinit ) returnable = RooMath::faddeeva( z_stl ).imag()*exp( -u*u );
				RooMathinit = true;
				pthread_mutex_unlock( &ROOT_Lock );
			}
			returnable = RooMath::faddeeva( z_stl ).imag()*exp( -u*u );
		}
		else
		{
			returnable = evalCerfApprox( swt, u, c ).imag();
		}
		return returnable;
#endif
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
			const double sigma_2=resolution*resolution;
			const double inv_tau2 = gamma*gamma;
			const double inv_r2_sigma = 1./(sqrt_2*resolution);
			const double inv_r2_sigma_tau = inv_r2_sigma*gamma;
			const double tau_thigh = invgamma*thigh;
			const double tau_tlow = invgamma*tlow;
			const double half_invgamma = 0.5*invgamma;
			const double half_sigma_2 = sigma_2*0.5;
			const double int_tHigh = erf( thigh*inv_r2_sigma )
				               - exp( (half_sigma_2 - tau_thigh )*inv_tau2 )
				              * erfc( (sigma_2 - tau_thigh )*inv_r2_sigma_tau ) ;

			const double int_tLow = erf( tlow*inv_r2_sigma )
                                              - exp( (half_sigma_2 - tau_tlow )*inv_tau2 )
					     * erfc( (sigma_2 - tau_tlow )*inv_r2_sigma_tau ) ;

			return half_invgamma * ( int_tHigh - int_tLow );
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

	pair<double, double> ExpCosSin( double t, double gamma, double deltaM, double resolution )
	{
		if(resolution > 0.) {
			const double c = gamma * resolution*_over_sqrt_2;
			const double u = (t / resolution) *_over_sqrt_2 ;
			const double wt = deltaM / gamma ;
			const complex<double> _cerf_plus = evalCerf( wt,-u,c );
			const complex<double> _cerf_minus = evalCerf( -wt,-u,c );

			const double ExpCos =  ( _cerf_plus.real() + _cerf_minus.real() ) * 0.25;
			const double ExpSin =  ( _cerf_plus.imag() - _cerf_minus.imag() ) * 0.25;
			return make_pair( ExpCos, ExpSin );
		}
		else {
			const double deltaM_t = deltaM * t;
			const double exp_val = exp( -gamma *t );
			const double ExpCos = exp_val * cos( deltaM_t );
			const double ExpSin = exp_val * sin( deltaM_t );
			return make_pair( ExpCos, ExpSin );
		}
	}


	pair<double,double> ExpCosSinInt( double tlow, double thigh, double gamma, double deltaM, double resolution  )
	{
		if( thigh < tlow ) {
			std::cerr << " Mathematics::ExpInt: thigh is < tlow " << std::endl ;
			return make_pair(-1.0,-1.0);
		}

		if( resolution > 0.) {
			//Added by Pete after getting code from Yuehong 120118
			const double c = gamma * resolution * _over_sqrt_2;
			const double inv_res = 1./resolution;
			const double umax = (thigh * inv_res ) *_over_sqrt_2 ;
			const double umin = (tlow * inv_res ) *_over_sqrt_2 ;
			const double wt = deltaM / gamma ;
			complex<double> evalDif( evalCerf(-wt,-umax,c) - evalCerf(-wt,-umin,c) );

			const double evalDifRe = evalDif.real();
			const double evalDifIm = evalDif.imag();

			const double factor = -0.5/gamma/(1+wt*wt);

			const double erf_max = RooMath::erf(-umax);
			const double erf_min = RooMath::erf(-umin);
			const double deltaCos = factor * ( evalDifRe + wt*evalDifIm + erf_max - erf_min );
			const double deltaSin = factor * ( -evalDifIm + wt*evalDifRe - -wt*(erf_max - erf_min) );

			return make_pair( deltaCos, deltaSin );
		}
		else
		{
			double real_tlow=tlow;
			if( tlow < 0. ) real_tlow = 0. ;

			const double deltaM_tlow=deltaM*real_tlow;
			const double gamma_deltaM = gamma*deltaM;
			const double deltaM_thigh=deltaM*thigh;

			const double factor = 1./( gamma_deltaM*gamma_deltaM );
			const double exp_lo = exp(-gamma*real_tlow);
			const double exp_hi = exp(-gamma*thigh);

			const double sin_hi = sin(deltaM_thigh);
			const double cos_hi = cos(deltaM_thigh);
			const double sin_lo = sin(deltaM_tlow);
			const double cos_lo = cos(deltaM_tlow);

			const double expCos = ( ( exp_lo * ( gamma*cos_lo - deltaM*sin_lo ))
					-( exp_hi * ( gamma*cos_hi - deltaM*sin_hi ))
					) * factor;

			const double expSin = ( ( exp_lo * ( gamma*sin_lo + deltaM*cos_lo ) )
					-( exp_hi * ( gamma*sin_hi + deltaM*cos_hi ) )
					) * factor;

			return make_pair( expCos, expSin );
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
			const double inv_res = 1./resolution;
			const double umax = (thigh * inv_res ) *_over_sqrt_2 ;
			const double umin = (tlow * inv_res ) *_over_sqrt_2 ;
			const double wt = deltaM / gamma ;
			complex<double> evalDif( evalCerf(-wt,-umax,c) - evalCerf(-wt,-umin,c) );

			const double evalDifRe = evalDif.real();
			const double evalDifIm = evalDif.imag();

			double deltaCos = -0.5/gamma/(1+wt*wt) * ( evalDifRe + wt*evalDifIm + RooMath::erf(-umax) - RooMath::erf(-umin) );
			return deltaCos;
		}
		else
		{
			double real_tlow=tlow;
			if( tlow < 0. ) real_tlow = 0. ;

			const double deltaM_tlow=deltaM*real_tlow;
			const double gamma_deltaM = gamma*deltaM;
			const double deltaM_thigh=deltaM*thigh;
			return (
					( exp(-gamma*real_tlow)* (gamma*cos(deltaM_tlow) - deltaM*sin(deltaM_tlow)))
					-( exp(-gamma*thigh)* (gamma*cos(deltaM_thigh) - deltaM*sin(deltaM_thigh)))
			       )/( gamma_deltaM*gamma_deltaM );
		}
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
			const double inv_res = 1. / resolution;
			const double umax = (thigh * inv_res ) *_over_sqrt_2;
			const double umin = (tlow * inv_res ) *_over_sqrt_2;
			const double wt = deltaM / gamma;
			complex<double> evalDif( evalCerf(-wt,-umax,c) - evalCerf(-wt,-umin,c) );

			const double evalDifRe = evalDif.real();
			const double evalDifIm = evalDif.imag();

			double deltaSin = -0.5/gamma/(1.+wt*wt) * ( -evalDifIm +   wt*evalDifRe -   -wt*(RooMath::erf(-umax) - RooMath::erf(-umin)) ) ;
			return deltaSin;
		}
		else
		{
			double real_tlow=tlow;
			if( tlow < 0. ) real_tlow = 0. ;

			const double deltaM_tlow=deltaM*real_tlow;
			const double gamma_deltaM = gamma*deltaM;
			const double deltaM_thigh=deltaM*thigh;
			return (	 ( exp(-gamma*real_tlow) * (gamma*sin(deltaM_tlow) + deltaM*cos(deltaM_tlow)) )
					-( exp(-gamma*thigh) * (gamma*sin(deltaM_thigh) + deltaM*cos(deltaM_thigh)) )
			       )/ ( gamma_deltaM * gamma_deltaM );
		}
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
		for (int i = 0; i < numAngularTerms; ++i) xi[i] = 0.0;

		double Sum[numAngularTerms];
		double Sum_sq[numAngularTerms][numAngularTerms];
		double cov[numAngularTerms][numAngularTerms];
		double cor[numAngularTerms][numAngularTerms];

		for (int i = 0; i < numAngularTerms; ++i) {
			Sum[i] = 0.0;
			for (int k = 0; k < numAngularTerms; ++k) {
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
		bool usePDF = testEval>=0. && !(std::isnan(testEval));

		for (int e = 0; e < numEvents; ++e)
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

			//val = 1./(8*TMath::Pi());
			//cout << f[0] << " " << f[1]<< " " <<  f[2]<< " " <<  f[3]<< " " <<  f[4]<< " " <<  f[5]<< " " << f[6]<< " " <<  f[7]<< " " <<  f[8]<< " " <<  f[9]<< endl;
			//cout << time << " " << cosTheta  << " " << phi << " " << cosPsi << " " << evalPDFraw << " " << evalPDFnorm << " " << val << endl;

			for (int i = 0; i < numAngularTerms; ++i)
			{
				Sum[i] = Sum[i] + f[i]/val;
				xi[i] += f[i]/val;

				for (int k = 0; k < numAngularTerms; ++k)
				{
					Sum_sq[i][k] += f[i]/val*f[k]/val;
				}
			}
		}

		cout << endl;

		cout << "Covariance matrix " << endl;
		for (int i = 0; i < numAngularTerms; ++i)
		{
			for (int k = 0; k < numAngularTerms; ++k)
			{
				cov[i][k] = 1./numEvents/numEvents * ( Sum_sq[i][k] - Sum[i]*Sum[k]/numEvents);
				cout << cov[i][k] << "\t";
			}
			cout << endl;
		}

		cout << "Correlation matrix " << endl;
		for (int i = 0; i < numAngularTerms; ++i)
		{
			for (int k = 0; k < numAngularTerms; ++k)
			{
				cor[i][k] = cov[i][k]/(sqrt(cov[i][i])*sqrt(cov[k][k]));
				cout << cor[i][k] << "\t";
			}
			cout << endl;
		}

		cout << "Weight +- error " << endl;
		for (int i = 0; i < numAngularTerms; ++i)
		{
			cout << fixed << setprecision(5) << xi[i]/numEvents << " \\pm " << sqrt(cov[i][i]) << endl;
			weights.push_back(xi[i]/numEvents);
		}
		return weights;
	}


	struct AccParam
	{
		AccParam(double * coefficients, int i_MAX, int j_MAX, int k_MAX, int l_MAX, int numEvents):
			pcoeff(coefficients)
			, l_max(l_MAX)
			, i_max(i_MAX)
			, k_max(k_MAX)
			, j_max(j_MAX)
			, num(numEvents)
		{}
		double parameterisation( double * omega, double * p )
		{
			(void) p;
			double returnValue(0.);
			double cosTheta(omega[0]);
			double phi     (omega[1]);
			double cosPsi  (omega[2]);
			double mKpi    (omega[3]);
			double Q_l(0.);  // mKpi
			double P_i(0.);  // cosPsi
			double Y_jk(0.); // cosTheta, phi
#ifdef __RAPIDFIT_USE_GSL
			for ( int l = 0; l < l_max+1; ++l )
			{
				for ( int i = 0; i < i_max+1; ++i )
				{
					for ( int k = 0; k < k_max+1; ++k )
					{
						for ( int j = 0; j < j_max+1; ++j ) // must have l >= k
						{
							if (j < k) continue;
							Q_l  = gsl_sf_legendre_Pl     (l,    mKpi);
							P_i  = gsl_sf_legendre_Pl     (i,    cosPsi);
							// only consider case where k >= 0
							// these are the real valued spherical harmonics
							if ( k == 0 ) Y_jk =           gsl_sf_legendre_sphPlm (j, k, cosTheta);
							else          Y_jk = sqrt(2) * gsl_sf_legendre_sphPlm (j, k, cosTheta) * cos(k*phi);
							//returnValue +=                pcoeff[             (i  * (k_max+1) + k) * (j_max+1) + j]*(      P_i * Y_jk);
							returnValue +=                pcoeff[((l*(i_max+1) + i) * (k_max+1) + k) * (j_max+1) + j]*(Q_l * P_i * Y_jk);
							//cout << i << k << j << " " << (i*(k_max+1)+k)*(j_max+1)+j << " " << pcoeff[(i*(k_max+1)+k)*(j_max+1)+j]/num << endl;
						}
					}
				}
			}
#else
			cerr << "Can't do this without GSL" << endl;
			exit(0);
#endif
			return returnValue/num;
		}
		double * pcoeff;
		int l_max;
		int i_max;
		int k_max;
		int j_max;
		int num;
	};


	int calculateAcceptanceCoefficients( IDataSet * dataSet, IPDF * PDF)
	{
		(void)PDF; // TODO hopefully the config is stored here
		// Implementation of the NIKHEF method for calculating the acceptance corefficients using Legendre polynomials and real valued spherical harmonics.
		cout << "Calculating acceptance coefficients" << endl;
		PhaseSpaceBoundary* boundary = dataSet->GetBoundary();
		const int l_max(6); // mKK
		const int i_max(6); // cos(theta_1)
		const int k_max(6); // phi
		const int j_max(6); // cos(theta_2)
		LegendreMomentShape lms; // All of the legwork is done in this class. Avoids duplicating lines of code.
		lms.SetMax(l_max+1, i_max+1, k_max+1, j_max+1);
		lms.Generate(dataSet, boundary, "mKK", "phi", "ctheta_1", "ctheta_2");
		lms.Save("LegendreMoments.root");
		TNtuple * dataacctree = new TNtuple("dataacctuple", "", "mKK:phi:ctheta_1:ctheta_2:weight");
		const double mK  = 0.493677; // TODO: read these from config somehow
		const double mBs  = 5.36677;
		const double mPhi= 1.019461;
		int numEvents = dataSet->GetDataNumber();
		for (int e = 0; e < numEvents; ++e)
		{
			// Retrieve the data point
			DataPoint * event = dataSet->GetDataPoint(e);
			double phi        = event->GetObservable("phi")->GetValue();
			double ctheta_1   = event->GetObservable("ctheta_1")->GetValue();
			double ctheta_2   = event->GetObservable("ctheta_2")->GetValue();
			double mKK        = event->GetObservable("mKK")->GetValue();
			// Calculate phase space element
			double p1_st  = DPHelpers::daughterMomentum(mKK, mK, mK);
			double p3     = DPHelpers::daughterMomentum(mBs,mKK,mPhi);
			double val    = p1_st*p3;
			dataacctree->Fill((float)mKK, (float)phi, (float)ctheta_1, (float)ctheta_2, (float)(1./val));
		}
		// Now sample the acceptance surface so that we can make projections to check that it looks sensible
		TNtuple * sampledtree = new TNtuple("sampledtuple", "", "mKK:phi:ctheta_1:ctheta_2:weight");
		gsl_qrng * q = NULL;
		try
		{
			q = gsl_qrng_alloc( gsl_qrng_sobol, 4);
		}
		catch(...)
		{
			cerr << "can't allocate random numbers for integration" << endl;
			exit(123);
		}
		unsigned int nSample(500000);
		vector<double> minima = {lms.mKK_min,-M_PI,-1,-1};
		vector<double> maxima = {lms.mKK_max,+M_PI,+1,+1};
		for ( int i = 0; i < (int)nSample; ++i )
		{
			double* point = new double[4];
			double mKK, phi, ctheta_1, ctheta_2, weight;
			vector<double*> point_mapped = {&mKK, &phi, &ctheta_1, &ctheta_2}; // So these can be looped over
			gsl_qrng_get( q, point );
			for ( unsigned int j = 0; j < 4; ++j )
			{
				*point_mapped[j] = point[j] * (maxima[j] - minima[j]) + minima[j];
			}
			weight = std::erf(186*(mKK)-2*mK)*lms.Evaluate(mKK, phi, ctheta_1, ctheta_2);
			sampledtree->Fill((float)mKK, (float)phi, (float)ctheta_1, (float)ctheta_2, (float)weight);
			delete[] point;
		}
		TFile * acceptance_file = TFile::Open("sampled_LegendreMomentShape.root","RECREATE");
		sampledtree->Write();
		dataacctree->Write();
		acceptance_file->Close();
		delete sampledtree;
		delete dataacctree;
		delete acceptance_file;
		return 1.;
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
		for (int i = 0; i < numAngularTerms; ++i) xi[i] = 0.0;

		for (int e = 0; e < numEvents; ++e)
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

			for (int i = 0; i < numAngularTerms; ++i) //For each event work out 10 f values, and 10 xi values and find the average
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


