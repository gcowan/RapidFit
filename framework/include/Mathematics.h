/**
	@namespace Mathematics

	Namespace holding some common mathematical functions that are used
	by many PDFs. For example, gaussian's convolved with trigonometric
	functions.

        @author Greig A Cowan, Pete Clalrke greig.cowan@cern.ch
        @date 2009-12-22
*/


namespace Mathematics 
{

//--------------------------- exp and exp*sin and exp*cos time functions -------------------------
// time functions for use in PDFs with resolution
	
	double Exp( double t, double gamma, double resolution )  ;
	double ExpInt( double tlow, double thigh, double gamma, double resolution  )  ;
	double ExpCos( double t, double gamma, double deltaM, double resolution ) ;
	double ExpCosInt( double tlow, double thigh, double gamma, double deltaM, double resolution  )  ;
	double ExpSin( double t, double gamma, double deltaM, double resolution )  ;
	double ExpSinInt( double tlow, double thigh, double gamma, double deltaM, double resolution  )  ;

}