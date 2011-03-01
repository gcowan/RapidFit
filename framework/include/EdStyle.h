// $Id: EdStyle.h,v 1.1 2009/11/10 10:35:44 gcowan Exp $
/**
        @class EdStyle

        Default style for RapidFit plots

        @author Greig A Cowan greig.alan.cowan@cern.ch
	@date 2009-10-02
*/

#ifndef EDSTYLE_H 
#define EDSTYLE_H 1

// Include files

#include "TStyle.h"
#include <string.h>

using namespace::std;

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
   
   virtual ~EdStyle( ); ///< Destructor

private:

   TStyle* edStyle;
   Int_t icol;
   Int_t font;
   Double_t tsize;
   

};
#endif // EDSTYLE_H
