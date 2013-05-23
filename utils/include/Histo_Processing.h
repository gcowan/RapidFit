#ifndef _HISTO_PROCESS
#define _HISTO_PROCESS

///	ROOT Headers
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TPolyMarker3D.h"
#include "TRandom.h"
#include "TString.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TMultiGraph.h"
///	Utils Headers
#include "Template_Functions.h"
///	System Headers
#include <vector>
#include <limits>
#include <float.h>

using namespace::std;

class Histogram_Processing
{
	public:

		/*! 
		 * @brief This has been adapted from the original code in RapidFit's Statistics code
		 *        It is intended to take a histogram and automatically rebin according to this function
		 *        As such it uses as much inbuilt functionality in root as possible
		 *        Return the ideal number of bins for a histogram of a vector of doubles
		 *        Uses D. Scott's method, published 1979
		 *
		 * This WILL rebin the input histogram according to the result from this method
		 *
		 * @param input_hist    This is the Histogram you wish to rebin in a given axis
		 *
		 * @param axis          This the axis you wish to look at (x=1, y=2, z=3) default is x
		 *
		 * @return Returns the Optimal Number of Bins as a signed integer
		 */
		static int OptimumBinNumber( TH1* input_hist, int axis=1 );

		/*!
		 * @brief Return the optimal number of bins for a given axis
		 *        Return the ideal number of bins for a histogram of a vector of doubles
		 *        Uses D. Scott's method, published 1979
		 *
		 * @param input_hist    This is the Histogram you wish to calculate the number of bins for.
		 *
		 * @return Returns the Optimal Number of bins as an unsigned integer
		 */
		static unsigned int GetOptimalBins( TH1* input_hist, int axis=1 );

		/*!
		 * @brief Optimally rebins the histogram you provide
		 *         By default this operates on the X axis assuming that you have a 1D histo but can rebin each axis independently upto 3D plots
		 *
		 *
		 * Behind the scenes this is a wrapper to the OptimumBinNumber method which stops the compiler complaining about any unused output
		 *
		 * @param input_hist    This is the Histogram you wish to rebin
		 *
		 * @param axis          This is the axis you wish to rebin
		 *
		 * @return Void
		 */
		static void OptimallyRebin( TH1* input_hist, int axis=1 );

		/*!
		 * @brief Get a TH1 from a vector of Doubles (could even write this to be type agnostric but I have the recast function in the template header
		 *
		 * @param input      This is the vector of values you wish to put into a Histogram
		 *
		 * @param rand       This is an optional pointer to a Random number generator to make sure the output TH1 has a unique name
		 *                   Slightly horrible solution to the HORRIBLE ROOT NAMESPACE
		 *                   You can safely pass NULL to this in doubt
		 *
		 * @param bins       This is the number of bins to use for this Histogram
		 *
		 * @param X_min      This is optional minimum of the histogram X axis
		 *
		 * @param X_max      This is the optiomal maximum of the histogram X axis
		 *
		 * @return Returns a pointer to a new TH1D that has been created
		 */
		static TH1* Get_TH1( const vector<Double_t>& input, TRandom* rand=NULL, int bins=100, double X_min=-DBL_MAX, double X_max=DBL_MAX );

		/*!
		 * @brief Get a TH2 from a vector of vectors which should be of outer dimension 2
		 *
		 * If you have a vector of vector which is of inner dimension 2 you can use the rotate() template method to change your object to be used here
		 *
		 * @warning If your input has an outer dimention > 2 only the first 2 are used. There is no warning about this!
		 *
		 * @param rand       This is an optional pointer to a Random number generator to make sure the output TH1 has a unique name                   
		 *                   Slightly horrible solution to the HORRIBLE ROOT NAMESPACE                                                                
		 *                   You can safely pass NULL to this in doubt
		 *
		 * @param bins1      This is the number of bins in X you wish to use
		 *
		 * @param bins2      This is the number of bins in Y you wish to use
		 *
		 * @return Returns a pointer to a new TH2D that has been created
		 */
		static TH2* Get_TH2( const vector<vector<Double_t> >& input, TRandom* rand=NULL, int bins1=100, int bins2=100 );

		/*!
		 * @brief Get a TH3 from a vector of vectors
		 *
		 * @param input 
		 *
		 *
		 * @param rand
		 *
		 * @param bins1
		 *
		 * @param bins2
		 *
		 * @param bins3
		 *
		 * @return Returns a pointer to a new TH3D that has been created.
		 */
		static TH3* Get_TH3( const vector<vector<Double_t> >& input, TRandom* rand=NULL, int bins1=100, int bins2=100, int bins3=100 );

		/*!
		 * @brief get a TGraph2D object from the provided vector vectors
		 *
		 * This uses the Get_TH2 method behind the senes to get a Histogram to create a TGraph2D from
		 *
		 * @return Returns a pointer to the new TGraph2D object which has been created
		 */
		static TGraph2D* Get_TGraph2D( const vector<vector<Double_t> >& input, TRandom* rand=NULL );

		/*!
		 * @brief Common interface to Get_{TH1,TH2,TH3}
		 *
		 * This allows you to be lazy (or generic) and just use a common interface for getting the correctly sized THx object from your input
		 *
		 * @return Returns a pointer to the THxD object that has just been created
		 */
		static TH1* Get_Histo( const vector<vector<Double_t> >& input, TRandom* rand=NULL, int bins1=100, int bins2=100, int bins3=100 );

		/*!
		 * @brief Get a histogram from your draw string, after weighting from the input_tree
		 *
		 * @return Returns a pointer to the new THxD object which has beeen creted
		 */
		static TH1* Get_Histo( TTree* input_tree, TString draw_str, TString weight_str, TRandom* rand=NULL );

		/*!
		 * @brief 
		 *
		 * @return Returns a pointer to the new TGraph object which has just been created
		 */
		static TGraph* Get_TGraph( const vector<vector<Double_t> >& input, TRandom* rand=NULL );

		//	Get a graph from your input draw string, after weighting from the input_tree
		static TGraph* Get_Graph( TTree* input_tree, TString draw_str, TString weight_str, TRandom* rand=NULL );

		//	Return the string which corresponds to the best fit function for the dataset by the best chi2
		//	Evaluates and compares the gaus, gamma & landau functions from ROOT
		static TString Best_Fit_Function( TH1* input, int OutputLevel=-1 );

		//	Fit the TH1 whilst catching a lot of unwanted output
		static void Silent_Fit( TH1* input_histo, TString fit_type, int OutputLevel=-1 );

		//	Draw something whilst catching all of the root output to the standard streams
		static void Silent_Draw( TCanvas* c1, TH1* input_histo, TString options="", int OutputLevel=-1 );

		//	Print whilst catching all of the root output to the standard streams
		static void Silent_Print( TCanvas* c1, TString Print_String, int OutputLevel=-1 );

		//	Addds the LHCb text to a plot
		static TPaveText* addLHCbLabel(TString footer, bool final=false );

		//	Plot a vector to a 1D file and return a pointer to the histogram created when this was done
		static TH1* Plot_1D( const vector<double>& input, TString Filename, TString Options, TRandom* rand=NULL, int bins=100, double min=-DBL_MAX, double max=DBL_MAX );

		//	Plot 2 vectors to a 2D file and return a pointer to the histogram created when this was done
		static TH2* Plot_2D( const vector<double>& X, const vector<double>& Y, TString Filename, TString Option, TRandom* rand=NULL );

		//	Plot 3 vectors to a 3D file and return a pointer to the histogram created when this was done
		static TH3* Plot_3D( const vector<double>& X, const vector<double>& Y, const vector<double>& Z, TString Filename, TString Option, TRandom* rand=NULL );

		/*!
		 * @brief Plot the data in the vector of vectors and return the histogram created when this was done
		 *
		 * @return Returns a pointer to the new Histogram which has been created
		 */
		static TH1* Plot( const vector<vector<double> >& input, TString Filename, TString Option, TRandom* rand=NULL );

		/*!
		 * @brief Get MultiGraph objects which contain the contours wanted from a TH2 histogram :D
		 *        Another short algorithm I'm overly proud of. It makes coding at a higher level SO MUCH EASIER :D :D :D
		 *
		 * This method takes the provided TH2 object and constructs the contours requested by the static levels defined in the contour list
		 * This collects the parts of each contour and constructs 1 TMultiGraph object for each contour requested.
		 *
		 * @param input_th2
		 *
		 * @param contour_list
		 *
		 * @param rand
		 *
		 * @return Returns a vector of pointer to new TMultiGraph objects which you can use normal
		 */
		static vector<TMultiGraph*> GetContoursFromTH2( TH2* input_th2, const vector<double>& contour_list, TRandom* rand );
};

#endif

