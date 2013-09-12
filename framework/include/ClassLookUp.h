/**
 * @class ClassLookUp
 *
 * @brief Central place to hold the methods for returning an instance of a class with a given name.
 *
 * @author Benjamin M Wynne bwynne@cern.ch
 * @author Robert Currie rcurrie@cern.ch
 */

#pragma once
#ifndef CLASS_LOOK_UP_H
#define CLASS_LOOK_UP_H

///	RapidFit Headers
#include "IPDF.h"
#include "BasePDF.h"
#include "PDFConfigurator.h"
#include "IFitFunction.h"
#include "IMinimiser.h"
#include "IDataGenerator.h"
#include "IPrecalculator.h"
///	System Headers
#include <iostream>
#include <vector>

using namespace::std;

class ClassLookUp
{
	public:
		/*!
		 * @brief Function to return a named PDF by looking for it's matching symbol in the compiled object
		 * 
		 * @param Name   Name of the PDF as defined in the source
		 *
		 * @param config instance of a PDFConfigurator which always has to be accepted by the PDF
		 *
		 * @return returns a pointer to the PDF with the Name provided, NULL if it isn't found in the object itself
		 */
		static IPDF* LookUpPDFName( string Name, PDFConfigurator* config );

		/*!
		 * @brief Copy the input PDF using the correct copy constructor object
		 * 
		 * Copy Constructor cached in IPDF once found
		 *
		 * @param InputPDF  The PDF will be given a pointer to it's copy constructor once found, and that will be used from then on
		 *
		 * @return a pointer to a new IPDF instance which is a copy of the input PDF to the highest derived instance
		 */
		static IPDF* CopyPDF( const IPDF* InputPDF );

		/*!
		 * @brief Function to return a named FitFunction by using a standard string comparison
		 *
		 * @pram Name  This is the name that will be searched for in the function
		 *
		 * @return This returns a pointer to the FitFunction instance
		 */
		static IFitFunction * LookUpFitFunctionName( string Name );

		/*!
		 * @brief Function to return a named Minimiser by using a standard string comparison
		 *
		 * @param Name This is the name of the Minimiser that will be searched for within this function
		 *
		 * @param NumberParameters This is the number of parameters in the fit (free + fixed) required by TMinuit so must be provided
		 *
		 * @return This returns a pointer to the Minimiser instance
		 */
		static IMinimiser * LookUpMinimiserName( string Name , int NumberParameters );

		/*!
		 * @brief Function to return a named DataGenerator by using a standard string comparison
		 *
		 * @param Name This is the name of the DataGenerator that will be searched for within this function
		 *
		 * @param InputPhaseSpace This is the PhaseSpace which will be filled by the DataGenerator
		 *
		 * @return Returns an IDataGenerator instance which is capable of creating Datasets
		 */
		static IDataGenerator * LookUpDataGenerator( string Name, PhaseSpaceBoundary* InputPhaseSpace, IPDF* InputPDF );

		/*!
		 * @brief Function to return a names PreCalculator by using a standard string comparison
		 *
		 * @param Name This is the name of the PreCalculator that will be searched for within this function
		 *
		 * @param WeightName This is the name of the Weight that we wish to create
		 *
		 * @param InputResult This is the FitResult which should be passed to the PreCalculator to interprate
		 *
		 * @param config  This allows the Precalculator to distinguish which component is the Signal and which is Noise
		 *
		 * @return Returns a PreCalculator instance which can weight a Dataset based on a FitResult and certain configuration
		 */
		static IPrecalculator * LookUpPrecalculator( string Name, string WeightName, FitResult* InputResult, unsigned int config );

		/*!
		 * @brief Get the path of the running executable on this system
		 *
		 * Arch independent, and an example of the windows code has not been dreamt up
		 * (and I don't care to unless someone throws their toys out of the pram!)
		 *
		 * @return a char pointer which contains the full running path of the executable object
		 */
		static char* getSelfPath();

		/*!
		 * @brief Resolve the object 'Name' within the binary object containing RapidFit
		 *
		 * @param Name   This is the Name of the object existing in this binary executable
		 *
		 * @return Returns the pointer to the object on disk/in memroy... This has to be a void pointer which is recast
		 */
		static void* getObject( string Name );

};

#endif

