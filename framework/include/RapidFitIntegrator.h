/*!
 * @class RapidFitIntegrator
 *
 * A numerical integrator to be used for projecting the PDF, and as a fallback if the PDF does not provide its own normalisation
 * This class uses two integrator classes provided by root: AdaptiveIntegratorMultiDim and GaussLegendreIntegrator.
 * Both of these assume the function to be integrated is well behaved.
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
 */

#pragma once
#ifndef RAPIDFIT_INTEGRATOR_H
#define RAPIDFIT_INTEGRATOR_H

///	ROOT Headers
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "Math/Integrator.h"
///	RapidFit Headers
#include "RapidFitIntegratorConfig.h"
#include "IntegratorFunction.h"
#include "FoamIntegrator.h"
#include "IDataSet.h"
#include "DebugClass.h"
///	System Headers
#include <string>
#include <vector>

class IPDF;
class FoamIntegrator;
class IntegratorFunction;

using namespace ROOT::Math;
using namespace::std;

static vector<DataPoint*> _global_doEval_points;
static vector<double> _global_range_minima;
static vector<double> _global_range_maxima;
static vector<string> _global_observable_names;

class RapidFitIntegrator
{
	public:
		/*!
		 * @brief Constructor with required objects
		 *
		 * This class acts as an abstraction layer between the PDF and RapidFit
		 *
		 * This class will return the appropriate Integral for a PDF which is either Analytically of Numerically derived
		 *
		 * This is intended to be a 'transparent' object in the way that it will just return the Integral of a given PDF as a double once it's been setup correctly
		 *
		 * @param InputPDF            This is the PDF which we are intending to return the Integral of
		 *                            This class DOES NOT make a copy of this PDF and unless changed with SetPDF the class will only ever use this PDF for all calculations
		 *
		 * @param ForceNumerical      default false: This class will determine if the PDF can analytically integrate and return this value if it can
		 *                                           will safely fall back to Numerical Integration
		 *
		 *                                    true:  This will force the class to
		 *                                           When this is true the class will not perform a comparison between the Analytical and Numerical as it's assumed to be not required
		 *
		 *                            The value of this is stored in the object RapidFitIntegratorNumerical
		 */
		RapidFitIntegrator( IPDF* InputPDF, bool ForceNumerical = false, bool UsePseudoRandomIntegration = false );

		/*!
		 * @breif Copy Constructor
		 *
		 * @warning         This Class DOES NOT take ownership of a given PDF and will not attempt to duplicate the PDF
		 *                  if you call multiple objects in parallel (i.e. multithreaded) be sure that you have seperate IPDF instances to avoid thread-unsafe operations taking place
		 *                  This can be achieved with SetPDF( IPDF* )
		 *
		 * @param Input     RapidFitIntegrator instance to duplicate internal state of
		 */
		RapidFitIntegrator ( const RapidFitIntegrator& );

		/*!
		 * Destructor Function
		 */
		~RapidFitIntegrator();

		void SetUpIntegrator( const RapidFitIntegratorConfig* config );

		void SetMaxIntegrationSteps( const unsigned int );

		void SetIntegrationAbsTolerance( const double input );

		void SetIntegrationRelTolerance( const double input );

		bool GetUseGSLIntegrator() const;
		void SetUseGSLIntegrator( const bool input );

		void SetFixedIntegralPoints( const unsigned int input );

		void SetNumThreads( const unsigned int input );

		/*!
		 *
		 * @brief Setup the Integrator for Projections (potentially speeds up the process slightly)
		 *
		 * @warning This will have no effect if compiled against ROOT < 5.28 as it relies on new features added here
		 *
		 * This actually changes the internal settings of the Numerical Integrator being used such that the Integral is 'good enough' for plotting
		 *
		 * This relies on the fact that you cannot see differences at say the <0.1% on a plotted projection.
		 */
		void ProjectionSettings();

		/*!
		 * @breif Get the Integral over this PhaseSpace for a given DataPoint template
		 *
		 * @param InputDataPoint      This is the DataPoint to be Integrated over
		 *
		 * @param IntegratePhaseSpace This is the PhaseSpaceBoundary to Integrate within
		 *
		 *
		 *
		 * This will Look to the IPDF* object which the class is currently using to determine if it should perform a Numerical Integral or Not
		 *
		 * If this PDF requires event by event Normalisation then only the Discrete Combination in the PhaseSpace that matches the DataPoint is used in calculation
		 *
		 *
		 * To get Analytical Integration:
		 *
		 * If the PDF can Analytically Integrate and RapidFitIntegratorNumerical is false  This class returns the Analytical Integral.
		 *
		 * Only under this set of Conditions is the PDF Analytical Integral compared to the Numerical Integral
		 *
		 *
		 * To get Numerical Integration:
		 *
		 * If the PDF cannot Integrate Analytically _OR_ RapidFitIntegratorNumerical is true This class will calculate and return the Numerical Integral
		 *
		 *
		 * If the PDF requests Event by Event Normalisation the Calculation is repeated every time this method is called. THIS IS EXTREMELY EXPENSIVE IN CPU TERMS!
		 *      In this configuration the Method only ever calculates the Integral corresponding to the Discrete Combination matching the provided DataPoint
		 *      This Uses the Method: NumericallyIntegrateDataPoint
		 *
		 * If the PDF can cache the Normalisation the Calculation is Perfomed Integrating over the whole PhaseSpace.
		 *      The resulting Normalisation is then Stored in the Normalisation Cache in the BasePDF.
		 *      This Uses the Method: NumericallyIntegratePhaseSpace
		 *
		 *
		 * @return This will always return a double > 0 unless an error occurs
		 */
		double Integral( DataPoint* InputDataPoint, PhaseSpaceBoundary* InputPhaseSpace );

		/*!
		 * @brief Test the Numerical vs Analytical Integral over a PhaseSpace for a given DataPoint template
		 *
		 * This is a Method to print out some information to the screen which allows for the User to determine if the Analytical and Numerical Integrals Agree (They Should!)
		 *
		 * This should only ever be performed once and after it has been performed the boolean haveTestedIntegral is set to true
		 *
		 * This can be disabled from being run as part of the call to Integral by using ForceTestStatus( true ), which sets the test as having completed
		 *
		 * @param InputDataPoint      This is required to be passed to the PDF being tested
		 *
		 * @param InputPhaseSpace     This is the PhaseSpace that will be Integrated over Analytically and Numerically
		 *
		 * @return Returns the correct Integral with the Analytical Integral being preferred
		 *
		 * @warning because it may return an Analytical Integral which may not be valid this is Never called if RapidFitIntegratorNumerical is true
		 *
		 * @warning If you compare the Numerical Integral to an Analytical for a PDF with per event Normalisation enabled expect some discrepancy which cannot be resolved!
		 */
		double TestIntegral( DataPoint* InputDataPoint, PhaseSpaceBoundary* InputPhaseSpace );

		/*!
		 * @brief Force the status of the Integration test
		 *
		 * This enables/disables the test of the Integral during a call to the Intergral Function.
		 *
		 * @warning   Note if it is your intent to enable the Test it will only be enabled for the First call to Integral!
		 *
		 * @param Input        true: The PDF will Not be tested on a call to Integral( DataPoint*, PhaseSpaceBoundary* )
		 *                     false: Forces the PDF to be tested only on the next call to Integral( DataPoint*, PhaseSpaceBoundary* )
		 * @return Void
		 */
		void ForceTestStatus( bool Input );

		/*!
		 * @brief Project a given observable over a given PhaseSpace for a template DataPoint for a given Parameter and a given component ID
		 *
		 * @param InputDataPoint  This is the DataPoint corresponding to the single Discrete Combination inthe PhaseSpace which is to be intergated over
		 *
		 * @param InputPhaseSpace This is the PhaseSpace that the DataPoint is to be integrated over.
		 *
		 * @param Projectile      This is the name of the Observable being projected over and as such the name needs to be added to the doNotIntegrate list
		 *
		 * @return the Integral of this Projectile
		 */
		double ProjectObservable( DataPoint* InputDataPoint, PhaseSpaceBoundary* InputPhaseSpace, string Projectile, ComponentRef* InputRef = NULL );

		/*!
		 * @brief Get the ratio between the 2 integrals
		 *
		 * This returns the Ratio of Analytical to Numerical Integrals
		 *
		 * It is used as a minor correction factor in Projections to correct for the fact that the Numerical Integral might under/over estimate the total integral by some small amount
		 *
		 * @warning This is a 'fudge factor'. There is no way of Fixing this problem completely,
		 *          When this Number is approximately you can trust your PDF and the Projections
		 *          When this Number is far from 1 you have a SERIOUS problem with your PDF and probably shouldn't be making projections!
		 *
		 * @return This returns a number which ideally should be exactly 1 and is a fudge factor to bring Numerical and Analytical Integrations inline
		 *
		 *         If you are relying on a Numerical Integral in your PDF this will return exactly 1. as there's nothing analytically to compare this to
		 */
		double GetRatioOfIntegrals();

		/*!
		 * @brief Get the internal PDF used by this integrator
		 *
		 * @warning This doesn't make any test that the PDF is valid or exists as an object
		 *
		 * @return Returns a Pointer to the PDF being used Internally
		 */
		IPDF * GetPDF() const;

		/*!
		 * @brief Set the internal PDF used by this integrator
		 *
		 * @param InputPDF  This replaces the internal PDF reference with this input
		 *
		 * @return Void
		 */
		void SetPDF( IPDF* InputPDF );

		/*!
		 * @brief Integrate over all possible discrete combinations within the Phase-Space
		 *
		 * This Will calculate the Integral over the whole PhaseSpace. This includes all Discrete Combinations
		 *
		 * The average value from the separate Discrete Combinations is returned
		 *
		 * @param InputPhaseSpace  This is the PhaseSpace that should be evaluated for all Discrete Combinations
		 *
		 * @param InputRef         This is the reference to the Component which should be integrated over the whole PhaseSpace
		 *
		 * @param DoNotIntegrate   These are the parameters Not to be Integrated over
		 *
		 * @return This Should return a double > 0 unless there has been an error
		 */
		double NumericallyIntegratePhaseSpace( PhaseSpaceBoundary* InputPhaseSpace, vector<string> DoNotIntegrate, ComponentRef* InputRef = NULL );

		/*!
		 * @brief Integrate over the given phase space using this given DataPoint. i.e. ONLY 1 UNIQUE discrete combination
		 *
		 * @param InputDataPoint   This is the DataPointe corresponding to the Discrete Combination in the PhaseSpace which should be integrated
		 *
		 * @prarm InputPhaseSpace  This is the PhaseSpace that the DataPoint should be integrated over
		 *
		 * @param DoNotIntegrate    These are the parameters Not to be Integrated over
		 *
		 * @param InputRef         This is the reference to the Component which should be integrated over the Discrete Combination
		 *
		 * @return This Should return a double > 0 unless there has been an error
		 */
		double NumericallyIntegrateDataPoint( DataPoint* InputDataPoint, PhaseSpaceBoundary* InputPhaseSpace, vector<string> DoNotIntegrate, ComponentRef* InputRef = NULL );


		/*!
		 * @brief Return a list of unclaimed and un-Integrable observables within this phase-space
		 *
		 * If an observable exists in the Phase-Space don't pass it to the Integrator, there will be errors.
		 * To avoid these slip ups on behalf of the user request that we check the typical DataPoint and add any unlcaimed observables to the list
		 *
		 * @return a vector of all non-Integrable Observables found in this DataPoint.
		 */
		vector<string> DontNumericallyIntegrateList( const DataPoint*, vector<string> = vector<string>() );

		void SetDebug( DebugClass* debug );

		static void clearGSLIntegrationPoints();

		static vector<DataPoint*> getGSLIntegrationPoints( unsigned int number, vector<double> maxima, vector<double> minima, DataPoint* templateDataPoint, vector<string> doIntegrate );

	private:
		/*!
		 * Don't Copy the class this way!
		 */
		RapidFitIntegrator& operator = ( const RapidFitIntegrator& );

		/*!
		 * @brief Calculate the Numerical Integral from the provided function
		 *
		 * Both NumericallyIntegrateDataPoint and NumericallyIntegratePhaseSpace reference this Function,
		 *
		 * However, it's much clear to work out what is going on with their function names and this saves duplicating code
		 *
		 * @param InputDataPoint    This is the DataPoint corresponding to the Discrete Combination in the PhaseSpace which should be integrated, if any in particular
		 *
		 * @param InputPhaseSpace   This is always the PhaseSpace that should be integrated over
		 *
		 * @param DoNotIntegrate    These are the parameters Not to be Integrated over
		 *
		 * @param InputRef          This is the reference to the Component which should be integrated over the Discrete Combination
		 *
		 * @param decision          This decides should we Integrate This PhaseSpace or DataPoint?
		 *
		 * @return This Should return a double > 0 unless there has been an error
		 */
		double DoNumericalIntegral( const DataPoint* InputDataPoint, PhaseSpaceBoundary* InputPhaseSpace, const vector<string> DoNotIntegrate, ComponentRef* InputRef = NULL, const bool decision = false );


		/*!
		 * @brief This is the Interface which allows us to perform psuedo-random number integration using GSL
		 *
		 * @return This Should return a double > 0 unless there has been an error
		 */
		static double PseudoRandomNumberIntegral( IPDF* functionToWrap, const DataPoint * NewDataPoint, const PhaseSpaceBoundary * NewBoundary, ComponentRef* componentIndex,
				vector<string> doIntegrate, vector<string> doNotIntegrate, unsigned int GSLFixedPoints=10000 );

		static double PseudoRandomNumberIntegralThreaded( IPDF* functionToWrap, const DataPoint * NewDataPoint, const PhaseSpaceBoundary * NewBoundary, ComponentRef* componentIndex,
				vector<string> doIntegrate, vector<string> doNotIntegrate, unsigned int num_threads=4, unsigned int GSLFixedPoints=10000, DebugClass* debug=new DebugClass(false) );

		/*!
		 * @brief This is the Interface to The MuliDimentional Integral class within ROOT
		 *
		 * This has been moved to be a static double in order to test which object weren't thread safe in an easier way
		 *
		 * @return This Should return a double > 0 unless there has been an error
		 */
		static double MultiDimentionIntegral( IPDF* functionToWrap, AdaptiveIntegratorMultiDim* thisIntegrator, const DataPoint * NewDataPoint, const PhaseSpaceBoundary * NewBoundary,
				ComponentRef* componentIndex, vector<string> doIntegrate, vector<string> dontIntegrate, bool, DebugClass* debug=new DebugClass(false) );

		/*!
		 * @brief This is the Interface to the 1D Integrator class within ROOT
		 *
		 * This has been moved into a separate function by itself to make code clearer, although a lot of it looks like the MultiDim case
		 *
		 * @return This Should return a double > 0 unless there has been an error
		 */
		static double OneDimentionIntegral( IPDF* functionToWrap, IntegratorOneDim * oneDimensionIntegrator, const DataPoint * NewDataPoint, const PhaseSpaceBoundary * NewBoundary,
				ComponentRef* componentIndex, vector<string> doIntegrate, vector<string> dontIntegrate, bool, DebugClass* debug=new DebugClass(false) );

		/*!
		 * This stores the Ratio of the Integrals after the Comparison between Analytical and Numerical
		 *
		 * This value is 1 if the PDF cannot Integrate or this class has been set to rely on the Numerical Integrator
		 *
		 * @return This should be a value close to or exactly 1
		 */
		double ratioOfIntegrals;

		/*!
		 * @brief Pointer to the FoamIntegrator object.
		 *
		 * This appears to be No longer Accessible?     Can someone think about this?
		 */
		FoamIntegrator * fastIntegrator;

		/*!
		 * @brief Pointer to the PDF that is being Integrated
		 */
		IPDF * functionToWrap;

		/*!
		 * @brief Internal MultiDimentional Integrator Instance
		 */
		AdaptiveIntegratorMultiDim * multiDimensionIntegrator;

		/*!
		 * @brief Internal 1D Integrator Instance
		 */
		IntegratorOneDim * oneDimensionIntegrator;

		/*!
		 * @brief Has this class been constructed to force the Integral to always be calculated Numerically
		 */
		bool RapidFitIntegratorNumerical;

		/*!
		 * @brief This contains the decision if the functionToWrap can Integrate   true: the PDF can Integrate (somehow)   false: the PDF can't Integrate
		 */
		bool functionCanIntegrate;

		/*!
		 * @brief Boolean that controls if the Analytical vs Numerical test has been performed
		 */
		bool haveTestedIntegral;

		/*!
		 * @brief This contains a list of Observables not claimed by the PDF present in the DataSet and so shouldn't be integrated
		 */
		vector<string> checked_list;

		/*!
		 * @brief A boolean to make sure that the Observables DataSet is only checked to be Integrable only once
		 */
		bool obs_check;

		/*!
		 * @brief to let the user choose if they want to use MC integration using pseudo-random numbers.
		 */
		bool pseudoRandomIntegration;

		DebugClass* debug;

		unsigned int num_threads;

		unsigned int GSLFixedPoints;


		//static vector<double*> initGSLDataPoints( unsigned int number, vector<double> maxima, vector<double> minima );

		static vector<DataPoint*> initGSLDataPoints( unsigned int number, vector<double> maxima, vector<double> minima, DataPoint* thisDataPoint, vector<string> doIntegrate );

};

#endif

