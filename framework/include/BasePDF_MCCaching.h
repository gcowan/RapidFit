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
#ifndef BASE_MCCACHING_PDF_H
#define BASE_MCCACHING_PDF_H

	///	ROOT Headers
#include "TRandom.h"
#include "TRandom3.h"
	///	RapidFit Headers
#include "IPDF_MCCaching.h"
#include "PhaseSpaceBoundary.h"
	/*#include "ObservableRef.h"
#include "PDFConfigurator.h"
#include "ComponentRef.h"
#include "PhaseSpaceBoundary.h"
#include "ParameterSet.h"*/
	///	System Headers
#include <vector>
#include <cmath>
#include <pthread.h>

	using namespace::std;

	class IPDF;

	class BasePDF_MCCaching : public virtual IPDF_MCCaching
{
	public:

		/*!
		 * @brief Explicitly use this Copy Constructor
		 *
		 *
		 * Derived PDF SHOULD write their own copy constructor if they include pointers to objects
		 *
		 * @param Input   Input BasePDF cast of input IPDF
		 */
		BasePDF_MCCaching(const BasePDF_MCCaching& Input);

		/*!
		 * @brief   Interface Function: Default Destructor
		 *
		 * Derived class should only write their own if they contain objects initialized with 'new'
		 */
		~BasePDF_MCCaching();

		/*!
		 * @brief Set the MC cache status
		 *
		 * @param Input    true = This PDF has an MC Cache Status on Disk , False This PDF does NOT have an MC Cache Status on Disk
		 *
		 * @return Void
		 */
		void SetMCCacheStatus( bool Input );

		/*!
		 * @brief Get the cache status
		 *
		 * @return         true = This PDF has an MC Cache Status on Disk , False This PDF does NOT have an MC Cache Status on Disk
		 */
		bool GetMCCacheStatus() const;

		/*!
		 * @brief Get the names of the MC cache files (without the .root extensions)
		 *
		 * This is a vector of filenames as there may be multiple Discrete Combinations in the PhaseSpace requested
		 * These each call the PDF separately to be sure that eg a tag+1 event is not modeled as being a tag0 event.
		 *
		 * @return	Vector of the Files containing the MC Caches on disk without the .root extentions on the filenames
		 */
		vector<string> GetMCCacheNames() const;

		/*!
		 * @brief Start a new TRandom3 instance with a given seed value
		 *
		 * @param num   This is the new number to use as a seed for creating a new TRandom3 instance and delete the existing one
		 *
		 * @return Void
		 */
		void SetRandomFunction( int num );

		/*!
		 * @brief Replace TRandom3 instance with a new function
		 *
		 * @param Input This is the new TRandom3 function we want to associate with this PDF
		 *
		 * @return Void
		 */
		void SetRandomFunction( TRandom3* Input );

		/*!
		 * @brief Get the seed number used to initialise the random number gen
		 *
		 * @return Returns the number used as the seed of the TRandom3 instance
		 */
		int GetSeedNum() const;

		/*!
		 * @brief Remove the cached files
		 *
		 * This also Sets the cache validity as false
		 *
		 * @return Void
		 */
		void Remove_Cache();

		/*!
		 * @brief Set whether this PDF can remove the MC cache object, i.e. is this instance the owner
		 *
		 * @param Input   false means that the PDF will not attempt to remove the TFoam disk instance on disk true means it will
		 *
		 * @return Void
		 */
		void Can_Remove_Cache( bool Input );

		/*!
		 * @brief Add a cache object
		 *
		 * @pram Name   Add a new ROOT file name (minus the .root) which is a cache of the TFoam generator
		 *
		 * @return Void
		 */
		void AddCacheObject( string Name );

		/*!
		 * @brief Get the Random function stored in this PDF
		 *
		 * @return pointer to the TRandom3 instance inside the PDF
		 */
		TRandom3* GetRandomFunction() const;

		void ChangePhaseSpace( PhaseSpaceBoundary * InputBoundary );

	protected:

		/*!
		 * Default  Constructor
		 */
		BasePDF_MCCaching();

	private:

		/*!
		 * Don't Copy the PDF this way!
		 */
		BasePDF_MCCaching& operator=(const BasePDF_MCCaching&);

		/*!
		 * List of files on disk which contain TFoam caches
		 * These should all be removed on exit or on call of Remove_Cache
		 */
		vector<string> cached_files;

		/*!
		 * Does this PDF have an MC cache generator stored on disk, i.e. is there a TFoam object on disk
		 */
		bool hasCachedMCGenerator;

		/*!
		 * vector of a single TRandom3 object for this PDF, this allows us to have a reproducable result for a defined seed
		 * the seed_function.empty() is used to see if this is defined, should probably check for NULL pointer, but oh well
		 */
		mutable TRandom3 * seed_function;

		/*!
		 * vector of single seed value used in the construction of the TRandom3, again should probably check for a null value, but oh well
		 */
		int seed_num;

		/*!
		 * Does this PDF control the on disk TFoam cache
		 * This allows the on disk cach to removed after RapidFit exists cleanly and only removed once
		 */
		bool do_i_control_the_cache;


};

#endif

