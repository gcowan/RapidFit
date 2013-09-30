/*!
 * @class BasePDF
 *
 * @brief   Class that provides a general implementation of most functions accessible through IPDF.
 *          Physics PDF can inherit from this to make a PDF without worrying about the details.
 *
 * It is intended that All Physics PDFs should inherit from this class in order to exploit a lot of pre-programmed features
 *
 * All new PDFs need to provide the following:
 *
 * REQUIRED:
 *
 * Overload:
 *          Evaluate( DataPoint* )
 *             (Required or there is not point!)
 *
 *          Normalisation( DataPoint*, PhaseSpaceBoundary* ) and/or Normalisation( PhaseSpaceBoundary* ) or set 'forceNumerical=true;'
 *             (Required some way of Normalising the Evaluate Function, return 1 if it already is)
 *
 *
 * BEST PRACTICE:
 *
 * Overload:
 *          All of the above +
 *
 *          SetParameterSet( ParameterSet* )
 *             (This is only called once per IMinimiser call so pre-calculate as much as possible here to save repeating it for every event!)
 *
 * Re-Use the variables:
 *                      bool forceNumerical;
 *                      bool e_by_e_Normalisation;
 *                      vector<string> allObservables;
 *                      ParamaterSet allParameters;
 *                      vector<string> doNotIntegrateList;
 *                      bool cacheValid;
 *
 * The PDF can reimplement these, but if we can all stick to the same coding style then it makes life easier.
 *
 *
 * ADVANCED PDFs:
 *
 * Overload:
 *
 *          EvaluateComponent( DataPoint, ComponentRef* )
 *             (This only needs to be done once you want to look at things like CP eigenstates)
 *             (You can use the ComponentRef object to store the result of lookups to reduce the length of time doing this
 *
 *          Can Overload PDFComponents or simply store all PDF components in: vector<string> component_list;
 *
 *          GetPhysicsParameters()
 *             (Overloading this allows you to mix/match which parameters you want based on PDFConfigurator conditions)
 *
 *          GetPrototypeDataPoint()
 *             (Overloading this allows you to mix/match which observables you want based on the PDFConfigurator conditions)
 *
 *          GetDoNotIntegrateList()
 *             (Overloading this allows you to mix/match which observables you can't model based on the PDFConfigurator conditions)
 *
 *
 * EXPERT/DEVELOPER:
 *
 * Overload:
 *          All of the above +
 *
 *          Integral(DataPoint*,PhaseSpaceBoundary*)
 *             (This has interfaces directly with the RapidFitIntegrator)
 *
 *          EvaluateForNumericIntegral( DataPoint* )
 *             (This interfaces with the NumericalIntegrator class)
 *
 *          EvaluateForNumericGeneration( DataPoint* )
 *             (This interfaces with the Toy DataSet Generators like Foam)
 *
 *
 *
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
 *
*/

#pragma once
#ifndef BASE_PDF_H
#define BASE_PDF_H

///	ROOT Headers
#include "TRandom.h"
#include "TRandom3.h"
///	RapidFit Headers
#include "IPDF.h"
#include "BasePDF_Framework.h"
#include "BasePDF_MCCaching.h"
#include "ObservableRef.h"
#include "PDFConfigurator.h"
#include "ComponentRef.h"
#include "PhaseSpaceBoundary.h"
#include "ParameterSet.h"
///	System Headers
#include <vector>
#include <cmath>
#include <pthread.h>

using namespace::std;

class BasePDF : public BasePDF_Framework, public BasePDF_MCCaching, public IPDF
{
	public:
		/*!
		 * Default  Constructor
		 */
		BasePDF();

		/*!
		 * @brief Explicitly use this Copy Constructor
		 *
		 *
		 * Derived PDF SHOULD write their own copy constructor if they include pointers to objects
		 *
		 * @param Input   Input BasePDF cast of input IPDF
		 */
		BasePDF(const BasePDF& Input);

		/*!
		 * @brief   Interface Function: Default Destructor
		 *
		 * Derived class should only write their own if they contain objects initialized with 'new'
		 */
		virtual ~BasePDF();

		/*!
		 * @brief Explicitly request caching to be turned OFF
		 *
		 * @return Void
		 */
		virtual void TurnCachingOff();

		/*!
		 * @brief   Interface Function: Is the Normalisation Cache Valid?
		 *
		 * The derived PDF MUST, either:   1) Set the internal object 'cacheValid' with Normalisation cache validity
		 *                                 2) Overload this function to return the Validity of the cache
		 *
		 * @param InputPoint     Input DataPoint corresponding to the Discete Combination being checked
		 *
		 * @param InputBoundary  Input PhaseSpaceBoundary containing all the Discrete Combinations
		 *
		 * @return        true = cache is valid , false = cache is not valid
		 */
		virtual bool CacheValid( DataPoint* InputPoint, PhaseSpaceBoundary* InputBoundary );

		/*!
		 * @brief   Interface Function: Is it valid to Cache the Numerical Integral
		 *
		 * Only specialized PDFs should overload this function!
		 *
		 * @return        true = Caching has been Enabled, false = Caching has NOT been Enabled
		 */
		virtual bool GetCachingEnabled() const;

		/*!
		 * @brief   Interface Function: Change the behaviour of the Caching in this PDF
		 *
		 * Only specialized PDFs should overload this function!
		 *
		 * @param Input   true = Enable Caching,  false Disable Caching
		 *
		 * @return Void
		 */
		virtual void SetCachingEnabled( bool );

		virtual void SetComponentStatus( const bool input );

		virtual bool GetComponentStatus() const;

		/*!
		 * @brief   Interface Function: Does the PDF want to use Numerical Normalisation
		 *
		 * @return        true = use Numerical normalisation , false = use Analytical normalisation
		 */
		bool GetNumericalNormalisation() const;

		/*!
		 * @brief Externally Set whether the PDF wants to use Numerical Normalisation
		 *
		 * @param Input   true = Use Numerical Normalisation , false = Use Analytical Normalisation
		 *
		 * @return Void
		 */
		void SetNumericalNormalisation( bool Input );

		/*!
		 * @brief Force the Normalisation Cache for this Discrete Combination given by the DataPoint to be defined as the Input
		 *
		 * @param Input   Value of the Normalisation as calculated externally to the PDF
		 *
		 * @param InputPoint  This is the DataPoint which corresponds to some Discrete Combination
		 *
		 * @param InputBound  This is the PhaseSpace which defines all possible Discrete Combinations
		 *
		 * @return Void
		 */
		void SetCache( double Input, DataPoint* InputPoint, PhaseSpaceBoundary* InputBound );

		/*!
		 * @brief Externally invalidate the Normalisation Cache
		 *
		 * This sets the cacheValid = true and the value of the stored normalisation cache to -1.
		 *
		 * @return Void
		 */
		void UnsetCache();

		/*!
		 * @brief   Interface Function:  Return the function value at the given point
		 *
		 * This MUST BE OVERLOADED BY derived PDF
		 *
		 * @param Input   DataPoint that should be Evaluated
		 *
		 * @return        (+ve) Returns the Likelihood of this DataPoint for this PDF (not required to be normalised)
		 */
		virtual double Evaluate( DataPoint* Input );

		/*!
		 * @brief   Interface Function:  This is a new wrapper between this PDF and the framework in RapidFit
		 *
		 * @param Input   This is the Input ParameterSet which is to be updated in the PDF
		 *
		 * @return Void
		 */
		void UpdatePhysicsParameters( ParameterSet* Input );

		/*!
		 * @brief Return the integral of the function over the given boundary
		 *
		 * @param InputPoint       DataPoint to be used as a guide for the Integration, i.e. tag +/-
		 * @param InputPhaseSpace  Whole PhaseSpace to be integrated, unless the DataPoint is relied upon by the PDF
		 *
		 * @warning       If both Normalise( PhaseSpaceBoundary* ) and Normalise( DataPoint*, PhaseSpaceBoundary* ) are defined
		 *                Normalise( PhaseSpaceBoundary* ) will take precedent!!!
		 *
		 * @return        Value of the PDF normalisation either from a valid cache of from an Analytical calculation
		 */
		double Integral( DataPoint* InputPoint, PhaseSpaceBoundary* InputPhaseSpace );

		/*!
		 * @brief Return the function value at the given point for generation
		 *
		 * In BasePDF this internally wraps around to Evaluate
		 *
		 * This can be, but isn't required to be overloaded by derived PDF
		 *
		 * @param Input   DataPoint to be evaluated for the Toy Dataset Generation tools in RapidFit
		 *
		 * @return        Result of Evaluating the Input DataPoint
		 */
		virtual double EvaluateForNumericGeneration( DataPoint* Input );

		/*!
		 * @brief Return the function value at the given point for use in numeric integral
		 *
		 * In BasePDF this internally wraps around to Evaluate
		 *
		 * This can be, but isn't required to be overloaded by derived PDF
		 *
		 * @param Input   DataPoint to be evaluated for the Integration tools in RapidFit
		 *
		 * @return        Result of Evaluating the Input DataPoint
		 */
		virtual double EvaluateForNumericIntegral( DataPoint* Input );

		virtual double EvaluateTimeOnly( DataPoint* Input );

		/*!
		 * @brief Interface Function: Return a prototype data point
		 *
		 * This returns the list of Observables which this PDF requires which can be passed to a DataPoint Constructor
		 *
		 * The derived PDF MUST, either:   1)   Put the names of all required Observables in member object allObservables
		 *                                 2)   Overload This Function to provide a vector of all required Observable Names
		 *
		 * @return        Vector of individual Names of the Observables required by the PDF
		 */
		virtual vector<string> GetPrototypeDataPoint();

		/*!
		 * @brief Interface Function: Get the function parameters
		 *
		 * Returns a Pointer to the ParameterSet being used internally by the PDF
		 *
		 * The derived PDF MUST, either:   1)   Use the member object allParameters for it's internal Physics ParameterSet
		 *                                 2)   Overload This Function to provide a link to it's internal ParameterSet object
		 *
		 * @return        Returns a pointer to the actual ParameterSet within a PDF, only modify it if you know what you're doing!
		 */
		virtual ParameterSet* GetPhysicsParameters();

		/*!
		 * @brief Interface Function: Return a prototype set of physics parameters
		 *
		 * Returns a list of the Physics Parameters this PDF requires
		 *
		 * The derived PDF MUST, either:   1)   Use allParameters for the internal ParameterSet (wraps around to ParameterSet::AllNames() )
		 *                                 2)   Overload this function to provide the list of all PhysicsParameters it requires
		 *
		 * @return        Vector of individual PhysicsParameter Names which can be passed to a ParameterSet constructor
		 */
		virtual vector<string> GetPrototypeParameterSet();

		/*!
		 * @brief Interface Function: Return a list of parameters not to be integrated by this PDF
		 *
		 * Returns a list of the Observables the PDF cannot correctly model
		 *
		 * The derived PDF MUST, either:   1)   Put the name of Observables it requires, but can't correctly integrate over in the member object doNotIntegrateList
		 *                                 2)   Overload This function to provide the list of Observables it can't integrate over
		 *
		 * @return        Vector of Observable Names that the PDF requires but doesn't fully model, e.g. mistag/time resolution (most detector influenced parameters)
		 */
		virtual vector<string> GetDoNotIntegrateList();

		/*!
		 * @brief Interface Function: Return a list of PDF components addresses in string format
		 *
		 * Returns a vector of available components in this PDF
		 *
		 * At a higher level is 0 is not a Component of this PDF it will be interrogated for a 0'th component regardless of wether the PDF claims to provide it
		 * Be sure if this or a safe default for it is defined!
		 *
		 * The derived PDF MUST, either:   1)   Do nothing and not expect to provide components
		 *                                 2)   Use the internal object component_list to provide a list of the components that the PDF can provide an Evaluate Statement for
		 *                                 3)   Provide a list of the Names of all of the Components it Hopes to plot
		 *
		 * Strings are used as they're much more adaptable to identify all of the components in a PDF even if it is nested within Sum/Product PDF objects
		 * Strings will be taken from the most base class and prepended until the final call to the whole derived object will give the full list of appended objects and names
		 *
		 * @return        Vector of strings which uniquely identify all of the unique components the PDF can provide
		 */
		virtual vector<string> PDFComponents();

		virtual string GetComponentName( ComponentRef* = NULL );

		/*!
		 * @brief Interface Function: Return the value of the component given by ComponentRef, component 0 by default
		 *
		 * When no ComponentRef is provided this wraps around to the Evaluate method
		 *
		 * When the Name of the ComponentRef the PDF should return the same value as the Evaluate Method.
		 * (This is not tested, but a 0'th component is ALWAYS added to PDFComponents when missing for design reasons)
		 *
		 * The derived PDF MUST, either:   1)   Do nothing and not expect to provide components
		 *                                 2)   Provide some method to evaluate the PDF for the given DataPoint(Unique Discrete Combination) for the requested CompoentRef
		 *
		 * @param InputDataPoint   This is the DataPoint (a single Unique Discrete Combination from the PhaseSpaceBoundary) being evaluated
		 * @param InputRef         This is the ComponentRef object being interrogated  (Optional)
		 *                         I make use of this and not a simple string as at a higher level I want to avoid a compound PDF performing many lookups per single call (of potentially hundreds)
		 *
		 * @return        Value of this single component evaluated by the PDF
		 */
		virtual double EvaluateComponent( DataPoint* InputDataPoint, ComponentRef* InputRef = NULL );

	protected:

		/*!
		 * @brief   Interface Function:  This function is called ONCE per call from Minuit
		 *
		 *
		 * Using this, complex variables from the ParameterSet can be calculated/cached.
		 *
		 * The derived PDF MUST, either:   1) Use the Internal ParameterSet allParameters
		 *                                 2) Overload this function but set the internal ParameterSet as appropriate
		 *
		 * @param Input    ParameterSet Normally as defined by the IMinimiser
		 *
		 * @return        true = success , false = fail... (this can probably be replaced with a Void return)
		 *
		 */
		virtual bool SetPhysicsParameters( ParameterSet* Input );


		/*!
		 * @brief Protected Function to Get the correct Cache for this DataPoint
		 *
		 * @return Returns the Cache for this DataPoint or calculates it if it needs re-calculating
		 */
		double GetCache( DataPoint*, PhaseSpaceBoundary* );

		/*!
		 * @brief This Checks the contents of the vector of caches.
		 *
		 * If it finds that this PhaseSpace requires a different number of caches it invalidates the current cache and resizes it
		 *
		 * @param InputBoundary This is the PhaseSpaceBoundary which contains all possible Data Combinations we want to cache a Normalisation for
		 *
		 * @return Void
		 */
		void CheckCaches( PhaseSpaceBoundary* InputBoundary );

		/*!
		 * @brief This is what actually turns the caching off becasue 'TurnCachingOff' has to be overloaded for Specialist PDFs
		 */
		void ReallyTurnCachingOff();

		bool CheckFixed( PhaseSpaceBoundary* NewBoundary );

		/*!
		 * @brief Protected Function for each PDF which provides a method for the PDF to analytically integrate over the whole phase space
		 *
		 * @param Input This is the PhaseSpace to be Integrated
		 *
		 * @return Returns the value of the Normalisation for the whole PhaseSpace
		 */
		virtual double Normalisation( PhaseSpaceBoundary* Input );

		/*!
		 * @brief Protected Function for each PDF which provides a method for the PDF to analytically integrate over the whole phase space using this datapoint
		 *
		 * @param InputDataPoint this is the DataPoint that is to be Integrated
		 *
		 * @param InputPhaseSpace This is the PhaseSpace to be Integrated
		 *
		 * @return Returns the value of the Normalisation for the whole PhaseSpace for this DataPoint
		 */
		virtual double Normalisation( DataPoint* InputDataPoint, PhaseSpaceBoundary* InputPhaseSpace );

		ParameterSet allParameters;		/*!	The internal ParameterSet object which contains all of the PhysicsParameters required for the PDF	*/

		vector<string> allObservables;		/*!	A list of all the Observable Names this PDF requires			*/
		vector<string> doNotIntegrateList;	/*!	A list of Observables that this PDF cannot correctly integrate over	*/

		vector<string> observableDistNames;	/*!	Undocumented	*/
		vector<IPDF*> observableDistributions;	/*!	Undocumented	*/

		vector<string> component_list;	/*!	This is the list of components that this PDF can provide	*/

		pthread_mutex_t* debug_mutex;	/*!	This is the pthread mutex object which is to be used for thread-locking sections of code in the PDFs that is not thread safe, e.g. using streamers*/

		bool can_remove_mutex;		/*!	This is internal to let the PDF know if it's safe to remove the shared mutex	*/

		DebugClass* debug;		/*!	This is the DebugClass which will allow you to turn on selective debugging at runtime	*/

		void ReallySetComponentStatus( const bool input );

		bool ReallyGetComponentStatus() const;

	private:

		bool numericalNormalisation;                  /*!     Does this PDF require Numerical Integration, or has Numerical Integration been requested?       */

		bool cachingEnabled;	/*!	Does this PDF require event by event normalisation	*/

		bool discrete_Normalisation;	/*!	Does this PDF require a unique normalisation to be calculated/applied depending on each unique discrete combination	*/

		mutable vector<double>* DiscreteCaches;	/*!	Normalisation Caches, each one corresponding to a unique Discrete Combination in the PhaseSpace		*/

		/*!
		 * Don't Copy the PDF this way!
		 */
		BasePDF& operator=(const BasePDF&);

		bool haveTestedIntegral;	/*!	Has the Integral Been tested?	*/
		bool requiresBoundary;		/*!	Does the PDF require the PhaseSpaceBoundary to calculate the normalisation?	*/

		bool fixed_checked, isFixed;
		size_t fixedID;

		bool _basePDFComponentStatus;

};

#endif


