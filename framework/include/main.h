
#pragma once

#ifndef _main_RapidFit_H
#define _main_RapidFit_H

#define _STR(x) #x
#define STR(x) _STR(x)

#include <vector>

using namespace::std;

/*!
 * This is required to be able to compile all of the code for running on the grid.
 *
 * As main represents the standard entry point we need to extract it away
 *
 * This is effectively main but with another name which allows us to call it within main and from elsewhere
 */
int RapidFit( vector<string> );

#endif

/*!
 *
 * @mainpage This is the Documentation for RapidFit
 *
 * @authors Edinburgh PPE Group
 *
 *
 * @section Intoduction
 *
 * RapidFit is an XML driven Minuit Fitting tool written in C++ and makes heavy use of inheritance and baseclasses
 *
 *
 * RapidFit is used in multiple Analyses with many PDFs coded up in the pdfs/ folder of the source code
 *
 * RapidFit supports Mutli-Threaded fitting as well as being compiled to run on the Grid as a ROOT library which allows for complex analyses to be reproduced quickly and efficiently
 *
 *
 * @section PDFs
 *
 * PDFs in RapidFit interface with the main codebase by providing methods to Evaluate( DataPoint* ) and Normalise( DataPoint*, PhaseSpace* )
 *
 *
 * @section Scans
 *
 * RapidFit is capable of performing Complex Scans in 1/2 dimension(s) whilst making full use of all of the underlying technology available to a normal fit which makes performing complex scans straight forward
 *
 * @section FeldManCousins
 *
 * RapidFit is capable of performing a complex FeldmanCousins Analysis to act as a cross-check with results which have been produced using simpler scanning code
 *
 * @section GOF
 *
 * @section XML
 *
 * RapidFit is capable of producing XML templates for various reasons as well as having it's own simple XML interpretor for portability and quick features
 *
 *
 * @section Usage
 *
 *
 *
 * For More information See the Source Code
 *
 */

/*!
 * @defgroup Generators  Generator Classes
 *
 * @defgroup Configurators Configuration Classes
 *
 * @defgroup Minimisers  Minimiser Objects
 *
 * @defgroup FitFunctions  FitFunction Objects
 *
 * @defgroup SpecialPDFs  Special PDFs
 *
 * @defgroup Study  RapidFit Studies
 *
 * @defgroup FitObjects  Normal Fit Objects
 *
 */

