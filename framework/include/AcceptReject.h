/*!
 * @ingroup Generators
 * @class AcceptReject
 *
 * @brief Class for generating toy data from a PDF.
 *        Can inherit from this to implement preselection for a particular PDF.
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 */

#pragma once
#ifndef ACCEPT_REJECT_H
#define ACCEPT_REJECT_H

///	ROOT Headers
#include "TRandom3.h"
///	RapidFit Headers
#include "IDataGenerator.h"
#include "MemoryDataSet.h"
#include "IPDF.h"
#include "PhaseSpaceBoundary.h"
#include "IDataSet.h"

using namespace::std;

class AcceptReject : public IDataGenerator
{
	public:
		/*!
		 * @brief AcceptReject Constructor
		 *
		 * @param InputPhaseSpace This is the PhaseSpace that is to be filled
		 *
		 * @param InputPDF This is the PDF to be evaluated in the PhaseSpace
		 * 
		 * @warning Neither of the Objects are copied by the class and it only ever refrences them!
		 */
		AcceptReject( PhaseSpaceBoundary* InputPhaseSpace, IPDF* InputPDF );

		/*!
		 * @brief Default Destructor
		 */
		~AcceptReject();

		/*
		 * @brief Interface Function to generate a dataset which is stored internally to this class
		 *
		 * @param Input the number of DataPoints requested
		 *
		 * @return the Number DataPoints created
		 */
		virtual int GenerateData( int Input );

		/*!
		 * @brief Interface Function to get a pointer to this dataset that has been generated
		 * 
		 * @warning This class does not delete the object refrenced by this pointer
		 *
		 * @return This passes the IDataSet to some external function, gives all control
		 */
		virtual IDataSet * GetDataSet() const;

	protected:
		/*!
		 * Don't Copy the class this way!
		 */
		AcceptReject ( const AcceptReject& );

		/*!
		 * Don't Copy the class this way!
		 */
		AcceptReject& operator = ( const AcceptReject& );

		virtual bool Preselection( DataPoint*, double );	/*!	Undocumented!	*/

		IPDF * generationFunction;			/*!	Pointer to the PDF given at construction	*/
		PhaseSpaceBoundary * generationBoundary;	/*!	Pointer to PhaseSpaceBoundary given at constructtion	*/

		int dataNumber;					/*!	Number of DataPoints in newDataSet	*/
		MemoryDataSet * newDataSet;			/*!	Pointer to internally stored dataset	*/
		TRandom3 * rootRandom;		/*!	Pointer to the TRandom3 object stored in generationFunction	*/

		double moreThanMaximum;		/*!	Undocumented!	*/
		int numberAttempts;		/*!	Undocumented!	*/
};

#endif

