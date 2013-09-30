/*!
 * @class IntegratorFunction
 *
 * @brief A wrapper to make IPDF usable by the numerical integrator
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
 */

#pragma once
#ifndef INTEGRATOR_FUNCTION_H
#define INTEGRATOR_FUNCTION_H

///	ROOT Headers
#include "Math/IFunction.h"
#include "TFoamIntegrand.h"
///	RapidFit Headers
#include "ObservableRef.h"
#include "DataPoint.h"
#include "PhaseSpaceBoundary.h"
#include "ComponentRef.h"
#include "DebugClass.h"
///	System Headers
#include <vector>

class IPDF;

using namespace::std;
///	This is required to use the ROOT::Math::IBaseFunctionMultiDim::Clone Copy function
using namespace ROOT::Math;

class IntegratorFunction : public IBaseFunctionMultiDim, public IBaseFunctionOneDim, public TFoamIntegrand
{
	public:
		/*!
		 * @brief Constructor Used for interfacing with ROOT Integrator classes
		 *
		 * This Constructor is expected to be called from within RapidFitIntegrate to supply an interface between the PDF and ROOT's Numerical Integral classes
		 *
		 * @param InputPDF      This is the PDF that is used for (in this constructor) Integration   i.e. this is the PDF which is called for evaluate
		 *
		 * @param InputPoint    This is the DataPoint instance which governs which Discrete Combination this PDF is to be evaluated at. Can only ever do one at a time by definition
		 *
		 * @param doIntegrate   This is the names of all of the Observables which are being integrated over. The size of this vector is the number of dimentions that are being integrated over
		 *
		 * @param dontIntegrate This is the names of all of the Observables which are NOT being integrated over. The DataPoints will always contain the Central Value if this Observable was constrained over a Range
		 *
		 * @param InputBoundary This is the PhaseSpace which is being integrated over. If a DataPoint is generated outside of the Range (shouldn't be but I simply don't trust ROOT) this forces Eval to return 0.0;
		 *
		 * @param WantedComponent  This is the ComponentIndex which corresponds to the Component Being Integrated over. By default this will give the same effect as Component '0' i.e. whole PDF
		 *
		 * @param min           This is the list of minimum values extracted from the PhaseSpace (This should probably be integrally generated, but it's redily available in the programs which call this constructor
		 *
		 * @param max           This is the list of maximum values extracted from the PhaseSpace (This should probably be integrally generated, but it's redily available in the programs which call this constructor
		 */
		IntegratorFunction( IPDF* InputPDF, const DataPoint* InputPoint, vector<string> doIntegrate, vector<string> dontIntegrate, const PhaseSpaceBoundary* InputBoundary, ComponentRef* WantedComponent= NULL,
				vector<double> min=vector<double>(), vector<double> max =vector<double>() );

		/*!
		 * @brief Constructor Used for interfacing with TFoam class
		 *
		 * This Constructor configures the class to be an interface between TFoam and the PDF
		 *
		 * @param InputPDF      This is the PDF that is used for (in this constructor) Integration   i.e. this is the PDF which is called for evaluate
                 *                                                         
                 * @param InputPoint    This is the DataPoint instance which governs which Discrete Combination this PDF is to be evaluated at. Can only ever do one at a time by definition
                 *                                                         
                 * @param doIntegrate   This is the names of all of the Observables which are being integrated over. The size of this vector is the number of dimentions that are being integrated over
                 *                                                         
                 * @param dontIntegrate This is the names of all of the Observables which are NOT being integrated over. The DataPoints will always contain the Central Value if this Observable was constrained over a Range
                 *  
		 * @param min           This is the list of minimum values extracted from the PhaseSpace (This should probably be integrally generated, but it's redily available in the programs which call this constructor
		 *
		 * @param max           This is the list of maximum values extracted from the PhaseSpace (This should probably be integrally generated, but it's redily available in the programs which call this constructor
		 *
		 * @param InputBoundary This is the PhaseSpace which is being integrated over. If a DataPoint is generated outside of the Range (shouldn't be but I simply don't trust ROOT) this forces Eval to return 0.0;
		 */
		IntegratorFunction( IPDF* InputPDF, const DataPoint* InputPoint, vector<string> doIntegrate, vector<string> dontIntegrate, vector<double> min, vector<double> range, const PhaseSpaceBoundary* InputBoundary );

		/*!
		 * @brief Return a pointer to the PDF being wrapped
		 *
		 * @return Returns a pointer to the IPDF object being used in the interface
		 */
		IPDF * GetWrappedFunction() const;

		/*!
		 * @brief Destructor
		 */
		~IntegratorFunction();

		/*!
		 * @brief ROOT copy functions
		 */
		using ROOT::Math::IBaseFunctionMultiDim::Clone;

		/*!
		 * @brief ROOT copy functions
		 */
		IntegratorFunction * Clone() const;

		/*!
		 * @brief ROOT copy functions
		 */
		IntegratorFunction * Clone(const char *newname="") const;

		/*!
		 * @brief Number of Free dimensions to integrate over
		 *
		 * @return Returns the size of the doIntegrate vector
		 */
		unsigned int NDim() const;

		/*!
		 * @brief Interface wrapper to DoEval for 1D
		 *
		 * This is a wrapper back to DoEval( double* )
		 *
		 * @param x  This is the value the ROOT method wishes us to evaluate
		 *
		 * @return Returns the value from the PDF for this value
		 */
		double operator()( Double_t x ) const;

		/*!
		 * @brief Interface wrapper to DoEval for n-D
		 *
		 * This is a wrapper back to DoEval( double* )
		 *
		 * @param x   This is a pointer to the array of values that the ROOT method wants us to evaluate
		 *
		 * @return Returns the value from the PDF for these values
		 */
		double operator()( const Double_t * x ) const;

		/*!
		 * @brief Return the density to Foam
		 * 
		 * From the code the method has to translate the array provided here to be compatible with the array to pass to the operator() or DoEval methods
		 *
		 * The array here is an array of multiples which are defined between 0 and 1 and are the relative position in each dimension from the minimum(0) and maximum(1) of the Observable Range
		 *
		 * (Difficult to describe but the code is perfectly clear)
		 *
		 * @warning not defined const in the ROOT API!
		 *
		 * @param ndim  This is the number of Dimensions (size of the array)
		 *
		 * @param x     This is the array of multiples that are to be used to calculate the Density at this point
		 *
		 * @return Returns the calculated density
		 */
		Double_t Density(Int_t ndim, Double_t * x );

		/*!
		 * @brief Interface function for Foam
		 *
		 * This is another Wrapper to DoEval( double* )
		 *
		 * @param x This is the value that the PDF is to be evaluated for
		 *
		 * @return Returns the value of the PDF for the given value
		 */
		double DoEval( double x ) const;

		/*!
		 * @brief Interface function for Foam
		 *
		 * This is the Function which most other calls wrap back to to save a large amount of duplicated code
		 *
		 * @param x This is the array of values that are to be passed to the PDF
		 *
		 * @return Returns the value of the PDF after being evaluated
		 */
		double DoEval( const double * x ) const;

		/*!
		 * @brief Copy Constructor
		 */
		IntegratorFunction ( const IntegratorFunction& );

		/*!
		 * @brief Assignment Operator
		 */
		IntegratorFunction& operator = ( const IntegratorFunction& );

		void SetDebug( DebugClass* debug );

		DataPoint* GetCurrentDataPoint() const;
	private:

		IPDF * wrappedFunction;					/*!	Pointer to the PDF being interfaced with ROOT objects		*/
		DataPoint * currentPoint;				/*!	Current DataPoint used to define which Discrete Combination being explored	*/
		vector<string> doIntegrate, dontIntegrate;		/*!	do and dont Integrable list of Parameters			*/
		vector<double> minima, ranges;				/*!	Minima and Range of allowed parameters in the PhaseSpace	*/
		vector<int> cache_positions;				/*!	Cached Results of Name Lookups for the do/dont Integrate lists	*/
		ComponentRef* componentIndex;				/*!	ComponentRef used when Interfacing with a Component to a PDF			*/
		mutable DataPoint * newDataPoint;			/*!	DataPoint that is passed to the PDF, mutable as we always change it befor call to PDF	*/
		vector<int> cache_lookup;				/*!	Cache the lookup for Parameter Names, reduce number of lookup operations	*/
		vector<double> lower_limit, upper_limit;		/*!	Upper and Lower Limits of the parameters being explored	*/

		/*!
		 * If the Numerical Normalisation Constructor was called, generateFunc = false integrateFunc = true
		 *
		 * If the Foam or Generator Constructor was called, generateFunc = true, integrateFunc = false;
		 *
		 * These were extremely useful in debugging, however they may no longer be required, but I would prefer to leave them in
		 */
		bool generateFunc, integrateFunc;

		PhaseSpaceBoundary* myPhaseSpaceBoundary;		/*!	Pointer to internally explorable PhaseSpace */

		DebugClass* debug;
};

#endif

