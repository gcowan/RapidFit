/*!
 * @interface IDataSet
 *
 * @brief Interface for all collections of data points
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 */

#pragma once
#ifndef I_DATA_SET_H
#define I_DATA_SET_H

//	RapidFit Headers
#include "DataPoint.h"
#include "PhaseSpaceBoundary.h"
//	System Headers
#include <string>
#include <vector>

using namespace::std;

class IDataSet
{
	public:
		/*!
		 * @brief Interface Function:
		 *        Return DataPoint at given number
		 *
		 * @warning It is Safe to alter the DataPoint object within the DataSet (although normally ill advised) but you CAN NOT alter the pointer itself, to do so WILL cause a malloc double free error
		 *
		 * @param Input   This is the number of the DataPoint in the whole DataSet you wish to get a pointer to
		 *
		 * @return Returns a Pointer to the DataPoint contained in the DataSet.
		 */
		virtual DataPoint * GetDataPoint( int Input ) const = 0;

		/*!
		 * @brief Interface Function:
		 *        Add a DataPoint to this DataSet
		 *
		 * @param Input   Add a DataPoint to This DataSet
		 *
		 * @return true, DataPoint was Added sucessfully false, DataPoint was not added
		 */
		virtual bool AddDataPoint( DataPoint* Input ) = 0;

		/*!
		 * @brief Interface Function:
		 *        Get the total number of points in this DataSet
		 *
		 * @return Returns the size of the DataSet (this could/should be an unsigned int, but this would require fixing a lot more warnings in RapidFit)
		 */
		virtual int GetDataNumber( DataPoint* templateDataPoint =NULL ) const = 0;

		/*!
		 * @brief Interface Function:
		 *        Get the phaseSpaceBoundary used to define the range that this dataset can span
		 *
		 * @return Returns a pointer to the PaseSpace that the DataSet was allowed to fill
		 */
		virtual PhaseSpaceBoundary* GetBoundary() const = 0;

		virtual void SetBoundary( const PhaseSpaceBoundary* ) = 0;

		/*!
		 * @brief Interface Function:
		 *        Sort the DataSet by a chosen parameter
		 *
		 * This allows you to sort the DataSet according to Where each DataPoint lies in a single Observable i.e. sort each event by 'time'
		 *
		 * @param ObsName   This is the Name of the Observable you wish to sort the DataPoints by
		 *
		 * @return Void
		 */
		virtual void SortBy( string ObsName ) = 0;

		/*!
		 * @brief Interface Function:
		 *        Get a discrete subset from within this dataset which matches the conditions passed as the input
		 *
		 * @param Names     These are the Names of the Observables which are Discrete
		 *
		 * @param Values    These are the values of all of the Discete Observables we want to get the corresponding DataPoint for out of the whole DataSet
		 *
		 * @return Returns a vector of DataPoints, This could be a seperate DataSet, but we do NOT duplicate the data or expect it to be changed/deleted this is effectively a wrapper to many calls to GetDataPoint
		 */
		virtual vector<DataPoint*> GetDiscreteSubSet( const vector<string> Names, const vector<double> Values ) const = 0;

		virtual vector<DataPoint*> GetDiscreteSubSet( const vector<ObservableRef> discreteParam, const vector<double> discreteVal ) const = 0;

		virtual IDataSet* GetDiscreteDataSet( const vector<ObservableRef> discreteParam, const vector<double> discreteVal ) const = 0;

		/*!
		 * @brief Interface Function:
		 *        Print out the dataset to provide extremely verbose output
		 *
		 * @warning This gets stupidly verbose for large amounts of data so use with caution
		 */
		virtual void Print() const = 0;

		/*!
		 * @brief Virtual Destructor
		 */
		virtual ~IDataSet() {};

		/*!
		 * @brief Returns an estimate of the total Yield
		 */
		virtual double Yield() = 0;
		virtual double YieldError() = 0;

		virtual void UseEventWeights( const string Name ) = 0;

		virtual bool GetWeightsWereUsed() const = 0;

		virtual string GetWeightName() const = 0;

		virtual double GetSumWeights() const = 0;

		virtual double GetSumWeightsSq() const = 0;

		virtual void ApplyAlpha( const double, const double ) = 0;

		virtual double GetAlpha() const =0;

		virtual void ApplyExternalAlpha( const string alphaName ) = 0;

		virtual void NormaliseWeights() = 0;

		virtual void ClearAllPseudoObservables() = 0;

	protected:
		/*!
		 * @brief Default Destructor
		 */
		IDataSet() {};

	private:

		/*!
		 * The general idea is to have one memory/disk/wherever resident DataSet for each unique set of data and to not duplicate this data.
		 *
		 * This gives us a lower Memory/Disk Footprint. It also means that we can rely on this object existing until the PDFWithData object that likely controls it is Destroyed
		 */

		/*!
		 * Don't Copy the class this way!
		 */
		IDataSet( const IDataSet& );

		/*!
		 * Don't Copy the class this way!
		 */
		IDataSet& operator=(const IDataSet& );
};

#endif


