// $Id: EdStyle.h,v 1.1 2009/11/10 10:35:44 gcowan Exp $
/*!
 * @class EdStyle
 *
 * @brief Default style for RapidFit plots
 *
 * @author Greig A Cowan greig.alan.cowan@cern.ch
 * @author Robert Currie rcurrie@cern.ch
*/

#pragma once
#ifndef EDSTYLE_H
#define EDSTYLE_H

///	ROOT Headers
#include "TStyle.h"
#include "TString.h"
#include "TText.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TCanvas.h"
///	System Headers
#include <string.h>

using namespace::std;

class EdStyle {
	public:
		/*!
		 * @brief Constructor
		 */
		EdStyle( );

		/*!
		 * @brief This sets the default and current style to be the same as the LHCb Style
		 */
		static void SetStyle();

		/*!
		 * @brief This will return the Latex formatted Unit for a given Parameter
		 *
		 * @param Input   This is the Unit of a Physics Parameter as used by the PDF
		 *                If we all stick to a standard naming convention of parameters in the XML
		 *                this allows RapidFit to spit out pre-formatted Latex tables :)
		 *
		 * @return TString containing the Latex formatted Unit
		 */
		static TString GetParamLatexUnit( string );

		/*!
		 * @brief This will return the ROOT formatted Unit for a given Parameter
		 *
		 * @param Input   This is the Unit of a Physics Parameter as used by the PDF
		 *                If we all stick to a standard naming convention of parameters in the XML
		 *                this allows RapidFit and RapidPlot to spit out nicer plots :)
		 *
		 * @return TString containing the ROOT formatted Unit
		 */
		static TString GetParamRootUnit( string );

		/*!
		 * @brief This will return the Latex formatted name for a given Parameter
		 *
		 * @param Input   This is the name of a Physics Parameter as used by the PDF
		 *                If we all stick to a standard naming convention of parameters in the XML
		 *                this allows RapidFit to spit out pre-formatted Latex tables :)
		 *
		 * @return TString containing the Latex formatted Name
		 */
		static TString GetParamLatexName( string );

		/*!
		 * @brief This will return the Latex formatted name for a given Parameter
		 *
		 * @param Input   This is the name of a Physics Parameter as used by the PDF
		 *                If we all stick to a standard naming convention of parameters in the XML
		 *                this allows RapidFit and RapidPlot to spit out nicer plots :)
		 *
		 * @return TString containing the ROOT formatted Name
		 */
		static TString GetParamRootName( string );

		/*!
		 * @brief This is a wrapper to GetParamLatexUnit( string )
		 */
		static TString GetParamLatexUnit( TString );

		/*!
		 * @brief This is a wrapper to GetParamROOTUnit( string )
		 */
		static TString GetParamRootUnit( TString );

		/*!
		 * @brief This is a wrapper to GetParamLatexName( string )
		 */
		static TString GetParamLatexName( TString );

		/*!
		 * @brief This is a wrapper to GetParamRootName( string )
		 */
		static TString GetParamRootName( TString );

		/*!
		 * @brief This tests the input to see if there is a defined Suffix and returns the starting character number of the Suffix
		 */
		static int Test_Suffix( string );

		/*!
                 * @brief This is a wrapper to Test_Suffix( string )
                 */
		static int Test_Suffix( TString );

		/*!
		 * @brief This will work out what the suffix on a branch Name is
		 *
		 * @param Input   This is the name of an Ntuple Branch from a RapidFit output
		 *
		 * @return This will return the Suffix from the Input
		 */
		static string Get_Suffix( string );

		/*!
		 * @brief This will work out what the suffix on a branch Name is
		 *
		 * @param Input   This is the name of an Ntuple Branch from a RapidFit output
		 *
		 * @return This will return the Suffix from the Input
		 */
		static TString Get_Suffix( TString );

		/*!
		 * @brief This takes a Branch name and Removes the '_pull'/'_gen'/'_value'/... and returns the Parameter Name
		 *
		 * @param Input   This is the name of an Ntuple Branch from a RapidFit output
		 *
		 * @return Returns the Name of the parameter - any suffix
		 */
		static TString Remove_Suffix( TString Input );

		static void FormatTText(TText*);

		static void FormatTLatex(TLatex*);

		/*!
		 * @brief This takes a Branch name and Removes the '_pull'/'_gen'/'_value'/... and returns the Parameter Name
		 *
		 * @param Input   This is the name of an Ntuple Branch from a RapidFit output
		 *
		 * @return Returns the Name of the parameter - any suffix
		 */
		static string Remove_Suffix( string );

		/*!
		 * @brief This will spit out a TLegend with the formatting required by the LHCb style...
		 *
		 * @return Returns a pointer to the TLegend object for plotting
		 */
		static TLegend* LHCbLegend();
		static TLegend* LHCbLeftLegend();
		static TLegend* LHCbBottomLegend();
		static TLegend* LHCbBottomLeftLegend();

		static TPaveText* LHCbLabel();
		static TPaveText* RightLHCbLabel();

		static Size_t GetLHCbTextSize();
		static Size_t GetLHCbAxisTextSize();
		static Width_t GetLHCbFunctionLineWidth();
		static Width_t GetLHCbAxisLineWidth();
		static Short_t GetLHCbTextAlign();
		static Style_t GetTransparentFillStyle();
		static Style_t GetSolidFillStyle();
		static Font_t GetLHCbFont();

		static TCanvas* RapidFitCanvas( TString Name, TString Title="" );

		/*!
		 * @brief Destructor
		 */
		~EdStyle( );

		static void SetColzPlotStyle();
	private:
		/*!
		 * Don't Copy the class this way!
		 */
		EdStyle ( const EdStyle& );

		/*!
		 * Don't Copy the class this way!
		 */
		EdStyle& operator = ( const EdStyle& );

		TStyle* edStyle;	/*!	Undocumented	*/
		Int_t icol;		/*!	Undocumented	*/
		Int_t font;		/*!	Undocumented	*/
		Double_t tsize;		/*!	Undocumented	*/
};

#endif

