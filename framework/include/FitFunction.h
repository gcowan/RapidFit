/*!
 * @class FitFunction
 *
 * @brief Parent class for the function to minimise
 * 
 * Overload the evaluate methods and UP value for Chi2, NLL, etc.
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
*/

#pragma once
#ifndef FIT_FUNCTION_H
#define FIT_FUNCTION_H

///	ROOT Headers
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
///	RapidFit Headers
#include "IFitFunction.h"
#include "PhysicsBottle.h"
#include "RapidFitIntegrator.h"
#include "RapidFitIntegratorConfig.h"
#include "ObservableRef.h"
#include "DebugClass.h"

#include <vector>
#include <string>

using namespace::std;

class FitFunction : public IFitFunction
{
	public:
		/*!
		 * @breif Default Constructor
		 */
		FitFunction();

		/*!
		 * @brief Default Detructor
		 */
		virtual ~FitFunction();

		/*!
		 * @brief Setup the Trace to record all of the output from the Minimiser
		 */
		void SetupTrace( const TString FileName, const int traceNum );

		/*!
		 * @brief Set the name of the numerical integration method to use 
		 */
		void SetIntegratorConfig( const RapidFitIntegratorConfig* gsl );

		/*!
		 * @brief Set the Physics Bottle to be used
		 *
		 * @param Input This is the PhysicsBottle that the FitFunction will use to calculate the result in Evalute
		 *
		 * @return Void
		 */
		void SetPhysicsBottle( const PhysicsBottle* Input );

		/*!
		 * @brief Get the Pointer to the Internal Physics Bottle
		 *
		 * @return Returns a pointer to the PhysicsBottle the FitFunction is using to Evaluate
		 */
		PhysicsBottle * GetPhysicsBottle() const;

		/*!
		 * @brief Change/Update the ParameterSet
		 *
		 * @param Input  This is the ParameterSet that the internal ParameterSet(s) should be chaned to
		 *
		 * @return true Changed the ParameterSet with no problems, false there was an error (I don't think this can be trusted and this should be made void)
		 */
		void SetParameterSet( const ParameterSet* );

		/*!
		 * @brief Get a Pointer to the Internal Parameter Set
		 *
		 * @return Returns a pointer to the Internal ParameterSet
		 */
		ParameterSet * GetParameterSet() const;

		/*!
		 * @brief Evaluate the FitFunction
		 *
		 * @return Returns a final value of the whole ParameterSet as a single double
		 */
		virtual double Evaluate();

		/*!
		 * @brief Set the Name of the Weights to use and the fact that Weights were used in the fit
		 *
		 * @param Name    This sets the name of the Weights to be used when Evaluating the DataSet
		 *
		 * @return Void
		 */
		void UseEventWeights( const string Name );

		/*!
		 * @brief Return if Weights were used as part of the analysis
		 *
		 * @return bool  true = Weights were used, false = Weights were NOT used
		 */
		bool GetWeightsWereUsed() const;

		string GetWeightName() const;
		
		/*!
		 * @brief Set the FitFunction to use Weights squared
		 *
		 * @param Input   true this causes Weights squared to be used in the fit, false don't use Weights squared
		 *
		 * @return Void
		 */
		void SetUseWeightsSquared( const bool Input );

		/*!
		 * @brief Set the Number of threads to be used by the FitFunction
		 *
		 * @param Input   This sets the number of threads that should be used by this FitFunction during Evaluate
		 *
		 * @return Void
		 */
		void SetThreads( const int Input );

		/*!
		 * @brief Get the Number of threads this FitFunction is using
		 *
		 * @return Returns the number of Threads this FitFunction is attempting to use
		 */
		int GetThreads() const;

		/*!
		 * @brief Set whether any RapidFitIntegrator Objects created internally should check the PDF/Numerical Integral
		 *
		 * @param Input  Should The Integrals be tested?  true = yes  false = no
		 *
		 * @return Void
		 */
		void SetIntegratorTest( const bool Input );

		/*!
		 * @brief What size of step in the function defines the error, 0.5 for NLL, 1 for chi2... etc.
		 *
		 * @input n  This is the nsigma Error you wish to calculate
		 *
		 * @return Returns the equivalent rise in the Function value that should be used to calculate the Error Value
		 */
		virtual double UpErrorValue( const int n );

		virtual vector<string> ConstrainedParameter() const;

		void SetDebug( DebugClass* debug );

		vector<string> GetNotConstrainedList() const;

		unsigned int GetCallNum();
	protected:
		/*!
		 * Don't Copy the class this way!
		 */
		FitFunction ( const FitFunction& );

		/*!
		 * Don't Copy the class this way!
		 */
		FitFunction& operator = ( const FitFunction& );


		/*!
		 * Overload these functions in child classes
		 */
		virtual double EvaluateDataSet( IPDF*, IDataSet*, int );

		PhysicsBottle * allData;			/*!	Undocumented	*/
		double testDouble;			/*!	Undocumented	*/
		bool useWeights;			/*!	Undocumented	*/
		string weightObservableName;		/*!	Undocumented	*/
		RapidFitIntegratorConfig* integrationConfig;

		//	This is for training Minuit
		//	(this could give some VERY cool graphs in ResultSpace :D )

		TFile* Fit_File;			/*!	Undocumented	*/
		TTree* Fit_Tree;			/*!	Undocumented	*/
		vector<Double_t> branch_objects;	/*!	Undocumented	*/
		vector<ObservableRef> branch_names;	/*!	Undocumented	*/
		Double_t fit_calls;			/*!	Undocumented	*/
		int traceNum;				/*!	Undocumented	*/
		void SetupTraceTree();

		int Threads;				/*!	Undocumented	*/
		vector<IPDF*> stored_pdfs;		/*!	Undocumented	*/
		vector<PhaseSpaceBoundary*> StoredBoundary;			/*!	Undocumented	*/
		vector< vector<vector<DataPoint*> > > StoredDataSubSet;		/*!	Undocumented	*/
		vector<RapidFitIntegrator*> StoredIntegrals;			/*	Undocumented	*/
		bool finalised;				/*!	Undocumented	*/
		struct Fitting_Thread* fit_thread_data;	/*!	Undocumented	*/

		bool testIntegrator;			/*!	Undocumented	*/

		bool weightsSquared;			/*!	Should Weights squared be used in evaluating the function	*/

		DebugClass* debug;				/*!	Flag controlling the turning on/off of extra debug information	*/

		string Name;

		double step_time;

		unsigned int callNum;

		vector<vector<IDataSet*> > stored_datasets;

};

#endif

