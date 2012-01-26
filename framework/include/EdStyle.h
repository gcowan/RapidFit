// $Id: EdStyle.h,v 1.1 2009/11/10 10:35:44 gcowan Exp $
/**
        @class EdStyle

        Default style for RapidFit plots

        @author Greig A Cowan greig.alan.cowan@cern.ch
	@date 2009-10-02
*/

#ifndef EDSTYLE_H 
#define EDSTYLE_H

// Include files

//	ROOT Headers
#include "TStyle.h"
#include "TString.h"
//	System Headers
#include <string.h>

using namespace std;

/** @class style EdStyle.h
 *  
 *
 *  @author Greig Cowan
 *  @date   2007-12-07
 */
class EdStyle {
	public: 
		/// Standard constructor
		EdStyle( ); 
   
		static void SetStyle();

		static TString GetParamLatexUnit( string );
		static TString GetParamRootUnit( string );
		static TString GetParamLatexName( string );
		static TString GetParamRootName( string );

                static TString GetParamLatexUnit( TString );
                static TString GetParamRootUnit( TString );
                static TString GetParamLatexName( TString );
                static TString GetParamRootName( TString );

		static int Test_Suffix( string );
		static int Test_Suffix( TString );  
		static string Get_Suffix( string );
		static TString Get_Suffix( TString ); 
		static TString Remove_Suffix( TString );

		static string Remove_Suffix( string );

		~EdStyle( ); ///< Destructor

	private:
		//	Uncopyable!
		EdStyle ( const EdStyle& );
		EdStyle& operator = ( const EdStyle& );
		TStyle* edStyle;
		Int_t icol;
		Int_t font;
		Double_t tsize;
};
#endif // EDSTYLE_H
