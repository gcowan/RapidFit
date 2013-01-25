/*!
 * @ingroup Generators
 * @class Foam
 *
 * @breif Class for generating toy data from a PDF.
 * 
 * Just a wrapper for the Root TFoam generator
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 */

#pragma once
#ifndef FOAM_H
#define FOAM_H

///	ROOT Headers
#include "TRandom3.h"
#include "TFoam.h"
#include "TFile.h"
///	RapidFit Headers
#include "IDataGenerator.h"
#include "MemoryDataSet.h"
#include "IntegratorFunction.h"
#include "ObservableRef.h"
///	System Headers
#include <vector>

using namespace::std;

class Foam : public IDataGenerator
{
	public:
		/*!
		 * @brief Constructor which will initialize TFoam to generate DataPoints over the PhaseSpace using the given PDF
		 *
		 * @param InputBoundary This is the PhaseSpaceBoundary which is used to be populated by the PDF. A copy is taken by this constructor
		 *
		 * @param InputPDF      This is the PDF which is used to create DataPoints which populate the PhaseSpace.
		 *                      This does NOT use an intenal PDF, it uses the PDF passed in this constructor.
		 */
		Foam( PhaseSpaceBoundary* InputBoundary, IPDF* InputPDF );

		/*!
		 * @brief Destruction
		 */
		virtual ~Foam();

		/*!
		 * @brief Interface Function to Generate a DataSet of requested size
		 *
		 * @param Input   This is the total number of DataPoints which should be Generated within this DataSet
		 *
		 * @return Returns the number of DataPoints actually constructed
		 */
		virtual int GenerateData( int Input );

		/*!
		 * @brief Interface Function to Return a pointer to the new DataSet
		 *
		 * @return Returns a pointer to the last DataSet constructed by this class (does NOT remember or alter previous DataSets
		 */
		virtual IDataSet * GetDataSet() const;

	protected:

		/*!
		 * Don't Copy the class this way!
		 */
		Foam ( const Foam& );

		/*!
		 * Don't Copy the class this way!
		 */
		Foam& operator = ( const Foam& );

		/*!
		 * Internal Function to Generate a new TFoam instance, or to attempt to access an existing object from an on disk cache
		 */
		void Init();

		/*!
		 * Internal Function to Remove the internal Generator Instance
		 */
		void RemoveGenerator();

		/*!
		 * A list of the number of files currently Open.
		 * These Files are closed on destruction of this class
		 */
		vector<TFile*> Open_Files;

		/*!
		 * A pointer to the requested PDF from Construction
		 */
		IPDF * InputPDF;


		IntegratorFunction * generationFunction;		/*!	Instance of the Wrapper Class between ROOT and the PDFs in RapidFit, this may be redudant, requires checking	*/
		PhaseSpaceBoundary * generationBoundary;		/*!	Pointer to the PhaseSpace to be filled from Construction		*/
		MemoryDataSet * newDataSet;				/*!	Internal Pointer to the DataSet that was requested			*/
		TRandom3 * rootRandom;					/*!	Pointer to the TRandom3 instance from the PDF provided			*/
		vector< TFoam* > foamGenerators;			/*!	A vector of all TFoam instances for this PDF and PhaseSpace		*/
		vector< IntegratorFunction* > storedIntegrator;		/*!	Vector of Integrator Functions Used as wrappers to this PDF from ROOT	*/
		vector< DataPoint* > storedDatapoint;			/*!	A vector of the DataPoints that have been created before they were committed to a DataSet	*/
		//int dataNumber;
		vector< vector<double> > discreteCombinations;		/*!	A vector containing the values of all of the Discrete Combinations in the Phase-Space		*/
		vector<string> allNames, discreteNames, continuousNames;/*!	The Names of the Doscrete and Continuous Observables in this PhaseSpace				*/
		vector<ObservableRef*> continuousNames_ref, discreteNames_ref;	/*!	vectors of pointers to ObservableRef objects to reduce expensive lookup operations	*/
		vector< vector<double> > discreteValues;		/*!	Undocumented, may be duplicate of discreteCombinations 			*/
		vector<double> minima, ranges;				/*!	Minima and Range of Each Observable in the PhaseSpace, required to be passed to IntegratorFunction as a cross check	*/
};

#endif

