/**
  @class PDFWithData

  A class for creating/storing a PDF and its associated data set

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-5
  */

#pragma once
#ifndef PDF_WITH_DATA_H
#define PDF_WITH_DATA_H

//	RapidFit Headers
#include "IPDF.h"
#include "PhaseSpaceBoundary.h"
#include "IPrecalculator.h"
#include "DataSetConfiguration.h"
#include "DebugClass.h"
//	System Headers
#include <string>
#include <vector>

using namespace::std;

class PDFWithData
{
	public:
		/*!
		 * @brief Constructor for PDFWithData
		 *
		 * @param InputPDF          This is the PDF that you wish to use to use to evaluate this dataset
		 *
		 * @param InputPhaseSpace   This is the PhaseSpace that the DataSet is allowed to Occupy
		 *
		 * @param InputConfig       This is the DataSetConfiguration instances which correspond to all of the DataSets you wish to Fit with this single PDF instance
		 */
		PDFWithData( IPDF* InputPDF, PhaseSpaceBoundary* InputPhaseSpace, vector<DataSetConfiguration*> InputConfig );

		/*!
		 * @brief Destructor
		 *
		 * @warning    By Default this will remove the IPDF which has control of the TFoam objects on disk which delete's them
		 *             It also deletes the IDataSet objects from memory.
		 *
		 *             If you use a PDF copied from this internal PDF and/or want to use a DataSet after removing this class be careful in what you do!
		 */
		~PDFWithData();

		/*!
		 * @brief Add a set of datasets to the cache, i.e. generate a load of datasets then call them sequentially
		 *
		 * @warning The destructor will destroy ALL datasets reference in the internal cache!
		 *
		 * @param Input A list of DataSets which are stored in the internal Cache and given to this instance to look after
		 *
		 * @return Void
		 */
		void AddCachedData( vector<IDataSet*> Input );

		/*!
		 * @brief Add a set of datasets to the cache, i.e. generate a load of datasets then call them sequentially
		 *
		 * @warning The destructor will destroy ALL datasets reference in the internal cache!
		 *
		 * @param Input A DataSet which is stored in the internal Cache and given to this instance to look after
		 *
		 * @return Void
		 */
		void AddCachedData( IDataSet* );

		/*!
		 * @brief Set the Physics Parameters from the given ParameterSet
		 *
		 * @param Input    This is the ParameterSet which you wish to use to change the internally stored ParameterSet
		 *
		 * @return true if no errors, false if error occured (this should probably be set to Void)
		 */
		bool SetPhysicsParameters( ParameterSet* Input );

		/*!
		 * @brief Get a pointer to the PDF used for generating/fitting the Data
		 *
		 * @return Returns a pointer to the PDF stored in this class
		 */
		IPDF * GetPDF() const;

		/*!
		 * @brief Get a new/cached dataset.
		 *        When the source is a 'File' it is permenantly cached in memory and by definition doesn't ever change
		 *        When the source is a Foam dataset it is created from the generator
		 *
		 * @return Returns a pointer to the DataSet. It is assumed that the DataSet is NOT destroyed externally to this class
		 */
		IDataSet * GetDataSet() const;

		/*!
		 * @brief Get the dataset Configuration used to generate the data
		 *
		 * @param Num   This is the number DataSetConfiguration instance wanted
		 *
		 * @return Returns a pointer to the DataSetConfiguration requested, first one in the vector by default
		 */
		DataSetConfiguration* GetDataSetConfig(int Num =0);

		/*!
		 * @brief Get all dataset Configurators
		 *
		 * @return Returns a vector with all of the DataSetConfiguration pointers inside
		 */
		vector<DataSetConfiguration*> GetAllDataSetConfigs();

		//Get a copy of the list of datasets that this class knows about
		vector<IDataSet*> GetCacheList();

		/*!
		 * @brief Remove all datasets in the cache
		 *
		 * @warning this will delete all DataSets in the cache from memory/remove access to them in RapidFit
		 *
		 * @return Void
		 */
		void ClearCache();

		/*!
		 * @brief Remove a specific dataset from the cache
		 *
		 * @warning  All other Datsets above this in the vector not drop down one!
		 *
		 * @param num   This is the number of theDataSet to remove from the cache.
		 *
		 * @return Void
		 */
		void RemoveFromCache( int num );

		/*!
		 * @brief Get a specific dataset from the cache
		 *
		 * @param num   This is the number of the DataSet in the cache that is requeted. NULL if out of limits
		 *
		 * @return Returns a Pointer to a DataSet in the cache
		 */
		IDataSet* GetFromCache( int num );

		/*!
		 * @brief Print some information Useful for debugging
		 *
		 * @return Void
		 */
		void Print() const;

		/*!
		 * @brief Returns the Whole section of a ToFit required to regenerate this Object
		 *
		 * @return Returns the required XML in a flat string object
		 */
		string XML() const;

		bool GetUseCache() const;

		void SetUseCache( bool );

		void SetDebug( DebugClass* debug );

		/*!
		 * @brief Set the Name of the Weights to use and the fact that Weights were used in the fit
		 *
		 * @param Name    This sets the name of the Weights to be used when Evaluating the DataSet
		 *
		 * @return Void
		 */
		void UseEventWeights( const string Name );

		bool GetWeightsWereUsed() const;

		string GetWeightName() const;

	private:
		/*!
		 * Don't Copy the class this way!
		 */
		PDFWithData ( const PDFWithData& );
		/*!
		 * Don't Copy the class this way!
		 */
		PDFWithData& operator = ( const PDFWithData& );

		IPDF * fitPDF;					/*!	This if the PDF which is used to Evaluate the DataSet					*/
		PhaseSpaceBoundary * inputBoundary;		/*!	This contains the PhaseSpace that the DataSets have been allowed to Occupy		*/
		bool parametersAreSet;				/*!	Undocumented	*/
		vector< DataSetConfiguration* > dataSetMakers;	/*!	This is the DataSetConfiguration objects which govern the creation of new DataSets	*/
		mutable vector< IDataSet* > cached_data;	/*!	This is the internal Cache of DataSets that this instance 'looks after'			*/

		mutable bool useCache;				/*!	Should PDFWithData return the last cached DataSet that it has (default false)		*/

		string WeightName;
		bool useWeights;

		DebugClass* debug;
};

#endif

