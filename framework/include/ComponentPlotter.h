/*!
 * @class ComponentPlotter
 *
 * @brief This class does the work involved in Numerically Projecting the Components of a PDF
 *
 *
 * This class does all of the work and contains all of the functions required to project components from a PDF
 *
 * This requires this class examining the PDF to determine the Components in it
 *
 * Projecting each of these Components over all of the separate Discrete Combinations Present in the PhaseSpace
 *
 * This has to be able to merge the output to give the total Projection for each component onto the dataset
 *
 * This has an output which has been designed such that the total projection from multiple compatible PDFs and multiple data subsets can be merged after this class has finished.
 *
 *
 * This class can also perform a Chi2 calculation of the projection using the PDF as a truth function to the dataset.
 *
 *
 * Some of this Class will seem non obvious or even confusing, I Apologise. However this class has to be able to deal with a lot of varying parameters and keep track of what is going on
 *
 *
 *
 *
 *
 * Methodology:
 *
 *
 * This PDF first loops over all components that the PDF advertises
 *
 * For each of these components the class loops over however many many requested points in the chosen Observable space
 *
 * For each of these points it projects the PDF and renormalises the value to compare it to data
 *
 *
 * Component 0 is defined as the total value from the PDF (the PDF should return the full Evaluate here
 *
 * Discrete Combination 0 is defined as the total of each of the Discrete Combinations within the PhaseSpace
 *
 *
 *
 *
 *
 * Normalisation:
 *
 * The normalisation in this class is tricky you have to compare a PDF which returns a likelihood over the whole Observable Space which varies and never exceeds 1 or falls below 0
 *
 * The data is binned and as such the PDF has to be renormalised to give the predicted amount of data that should be in this bin from a dataset of a given size.
 *
 * This requires knowledge of the yield and such
 *
 *
 *
 *
 *
 * Chi2/N:
 *
 * The Chi2 test requires everything else to be working and has been tested on toys and code gives a correct theory curve and sensible Chi2 values
 *
 *
 * A projection is just a fancy way of saying that the PDF has been integrated over all free parameters while the Observable requested is fixed.
 *
 *
 * @author Robert Currie rcurrie@cern.ch
 */

#pragma once
#ifndef COMPONENT_PLOTTER_H
#define COMPONENT_PLOTTER_H

///	System Headers
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
///	RapidFit Headers
#include "IPDF.h"
#include "IDataSet.h"
#include "RapidFitIntegrator.h"
#include "EdStyle.h"
#include "PhaseSpaceBoundary.h"
#include "ComponentPlotter_config.h"
///	System Headers
#include <vector>
#include <string>

using namespace::std;

class ComponentPlotter
{
	public:
		/*!
		 * @brief Constructor With Correct Number of input Objects
		 * Unlike the previous projection Plotting Code this Class has to exist as a single instance which has knowledge of the actual Component it is running over
		 *
		 * @param InputPDF      This is the PDF which is being projected
		 *
		 * @param InputDataSet  This is the DataSet which is being compared to the projection
		 *
		 * @param PDFStr        The path that will be created in the outputfile to be unique for each PDF
		 *
		 * @param filename	This is the File which contains the ROOT object outputs
		 *
		 * @param Name          This is the Name of the Observable being Projected
		 *
		 * @param config        This is the object which carries configuration information from the XML to this class for projecting the data
		 *
		 * @warning This class will make copies of the PDF and config, however it will NOT make copies of the DataSet and file pointer which should always remain valid
		 */
		ComponentPlotter( IPDF* InputPDF, IDataSet* InputDataSet, TString PDFStr, TFile* filename, string Name, CompPlotter_config* config=NULL, int PDF_Num=0, const DebugClass* debug=NULL );

		/*!
		 * @brief Destructor Function
		 */
		~ComponentPlotter();

		/*!
		 * @brief Main hooks to the rest of RapidFit
		 *
		 * Wrapper function to GenerateProjectionData
		 *
		 * @return Void
		 */
		void ProjectObservable();

		/*!
		 * @brief Let the Class know that sWeights were used.
		 *
		 * This is a separate function as it initialises several internal objects which are set up differently for sWeights and not having sWeights
		 *
		 * @param Input  This is the name of the sWeight
		 *
		 * @return Void
		 */
		void SetWeightsWereUsed( string Input );

		/*!
		 * @brief Function to Extract the pointers to the various Component Projections
		 *
		 * Each Discrete Data SubSet has a Discretely different Set of Projections, hence the nested vector structure
		 *
		 * @warning This class will remove the objects when it is destroyed!
		 *
		 * @return This returns the pointers to the components plotted and stored in TGraphs
		 */
		vector<vector<TGraph*> > GetComponents();

		/*!
		 * @brief Get a vector of the different Data SubSets which correspond to unique Discrete Combinations, with the 0th component being all Data
		 *
		 * @warning This class will remove the objects when it is destoryed!
		 *
		 * @return This returns the Binned Data Subsets corresponding to the different configurations
		 */
		vector<TGraphErrors*> GetBinnedData();

		/*!
		 * @brief Get the Chi2 values after they have been calculated
		 *
		 * @return This returns the pair or make_pair< chi2, N >
		 */
		pair<double,double> GetChi2Numbers();

		/*!
		 * @brief Return Data which is to be used by the Pull Plot
		 *
		 * @return This returns a vector of values for each mass bin
		 *         This has dimension MassBins x 3
		 *         Each nested vector contains the total PDF value at the lower edge, center and upper edge of each bin
		 */
		vector<double> GetFunctionEval();

		static TGraphErrors* PullPlot1D( vector<double> input_bin_theory_data, TGraphErrors* input_data, string observableName, string CombinationDescription,

				TRandom* rand = NULL, CompPlotter_config* conf=NULL, DebugClass* debug = NULL );

		/*!
		 * @brief A Function which can correctly merge seperate datasets into a single TGraphErrors object
		 *
		 * @param input This contains all of the Binned data subsets which you wish to merge
		 *
		 * @param rand This is strictly optional as the class will pick up gRandom, but it's nicer to be explicit
		 *
		 * @returns a single TGraphErrors instance which contains all of the binned data subsets with the correct errors per bin
		 */
		static TGraphErrors* MergeBinnedData( vector<TGraphErrors*> input, TRandom* rand=NULL );

		/*!
		 * @brief A Function which can correctly merge the Components projected over multiple Data SubSets to get the overall Projection
		 *
		 * @param input This contains all of the projections (all components) which you wish to merge
		 *
		 * @param rand This is strictly optional as the class will pick up gRandom, but it's nicer to be explicit
		 *
		 * @return This returns a single set of all of the components merged into a single component projection set
		 */
		static vector<TGraph*> MergeComponents( vector<vector<TGraph*> > input, TRandom* rand=NULL );


		/*!
		 * @brief Helper Function which actually Produces the final plots which are written to an already open ROOT file and images to disk
		 *
		 * @param input_data       This is the data which has been pre-binned into a TGraphErrors instance
		 *
		 * @param input_components These are the various component projections from a PDF
		 *
		 * @param observableName   This is the name of the Observable which has been projected
		 *
		 * @param CombinationDescription This is the String which describes this particular Discrete Dataset in terms of the Discrete Observables in the PhaseSpace
		 *
		 * @param total_boundary   This is the total Phase-Space whcih has been integrated over (except for observableName, of course)
		 *
		 * @param rand             This is strictly optional as the class will pick up gRandom, but it's nicer to be explicit
		 *
		 * @param config           This contains the Configuration from the XML, if none is provided this will go with whatever ROOT wants to do
		 *
		 * @param debug            Satandard Debugging class interface in RapidFit
		 *
		 * @param input_bin_theory_data  This contains the (optional) Data required to produce a plot of the residuals across the fit
		 *
		 * @return Void
		 */

		static void OutputPlot( TGraphErrors* input_data, vector<TGraph*> input_components, string observableName, string CombinationDescription, PhaseSpaceBoundary* total_boundary,

				TRandom* rand=NULL, CompPlotter_config* config=NULL, DebugClass* =NULL, vector<double> input_bin_theory_data=vector<double>() );

		/*!
		 * @breif helper function to write out the numerical data from the Final merged component plot into a TTree to make any post-projection work easier
		 *
		 * @param Total_BinnedData  This is the TGraphErrors used in the projection of the data in a plot
		 *
		 * @param Total_Components  This is the vector of TGraphs corresponding to the various components in the final Projection plot
		 *
		 * @param destination       This is the name of the TTree which should be used for the output
		 *
		 * @param config            This will be used to assign any custom titles to the branches in the TTree which we are using
		 *
		 * @return Void
		 */
		static void WriteData( TGraphErrors* Total_BinnedData, vector<TGraph*> Total_Components, TString destination );

		/*!
		 * @brief Required for wrapping this class in a TF1 for Chi2 calculations
		 *
		 * @param x  this is a pointer to the Oservable Value
		 *
		 * @param p  this is explicitly unused but is there to interface with ROOT
		 *
		 * @return This returns the Projection at the observable value X and allows the PDF to be truth for a Chi2 test
		 */
		double operator() (double *x, double *p);


		void SetDebug( DebugClass* debug );

		static string XML( const int projectionType=1 );

	private:

		/*!
		 * @brief Uncopyable!
		 * No good reason to want to copy these objects so these won't be written
		 */
		ComponentPlotter ( const ComponentPlotter& );

		/*!
		 * @brief Uncopyable!
		 * No good reason to want to copy these objects so these won't be written
		 */
		ComponentPlotter& operator = ( const ComponentPlotter& );



		//	Functions used to move from the Interface functions to perform the actual Projections
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		/*!
		 * @brief Project this component of the observable
		 *
		 * @param InputDataPoint    This datapoint determines which Discrete Combination is to be Projected
		 *
		 * @param ObservableName    This is the Observable being Projected
		 *
		 * @param min               This is the minimum value that the Observable takes
		 *
		 * @param num_of_steps      This is the total number of steps to take in the Observable
		 *
		 * @param step_size         This is the distance between points in the Observable being Projected
		 *
		 * @param compName          This is the name of the Component which is being projected by the PDF
		 *
		 * @return This returns a pointer to a list of the integrals of the PDF for this given component at each of the requested coordinates in the Observable Space
		 */
		vector<double>* ProjectObservableComponent( DataPoint* InputDataPoint, string ObservableName, double min, int num_of_steps, double step_size, string compName );

		/*!
		 * @brief Generate the data for projection plots
		 *
		 * @return Void
		 */
		void GenerateProjectionData();



		//	Functions for Projection Data for Plotting:
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		/*!
		 * @brief Construct the x data points for all conditions for the projections
		 *
		 * This is a simple algorithm that calculates the value of the X axis for each of the projections which are performed and makes a copy for each Discrete Combination
		 *
		 * @param Input  This is the number of points in the Observable which will be Projected
		 *
		 * @return This returns the pointer to the vector of copies of all of the coordinates in X
		 */
		vector<vector<double>* >* MakeXProjectionData( unsigned int Input );

		/*!
		 * @brief Construct the y data points for all conditions for the projections of a requested Component
		 *
		 * @param Name  This is the name of the Component being interrogated
		 *
		 * @return This returns the Y coordinates of the various Projections in the chosen Observable
		 */
		vector<vector<double>* >* MakeYProjectionData( string Name );

		/*!
		 * @brief When we have more than 1 discrete component we need to create component 0 which contains the total PDF result at this coordinate
		 *
		 * This is best written as a separate function to make MakeYProjectionData a bit less cluttered
		 *
		 * @param input     This is a pointer to the projections for multiple discrete combinations of one component from one PDF
		 *
		 * @return This returns the input with a vector prepended to the start of the set which is a total of the vectors given at input
		 */
		vector<vector<double>* >* GenerateComponentZero( vector<vector<double>* >* new_dataset );



		//	Functions for Plotting the output and putting it in a ROOT FILE
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		/*!
		 * @brief Write the Projection Data to a file
		 *
		 * @param X   This contains the X coordinates data of all of the Combinations and Discrete Combinations
		 *
		 * @param Y   This contains the Y coordinate data of all of the Combinations and Discrete Combinations
		 *
		 * @param combinationDescription  This is a vector of descriptions of each of the Discete Combinations from the PhaseSpace
		 *
		 * @return Void
		 */
		void WriteOutput( vector<vector<vector<double>* >* >* X, vector<vector<vector<double>* >* >* Y, vector<string> combinationDescriptions  );

		/*!
		 * @brief Get a TH1 containing this discrete combination
		 *
		 * @param combinationNumber  This is the Discrete Combination number
		 *
		 * @return returns a TH1 contining the Data for the Given Discrete Combination from the PhaseSpace
		 */
		TH1* FormatData( unsigned int combinationNumber );

		/*!
		 * @brief This is a slimmed down version of the AddBranch function I have written elsewhere
		 *
		 * This requires that the branch_data has the same number of points as already in the tuple and will ONLY write the values as doubles.
		 *
		 * @param input_tree   This is the Tree we want to add a branch to
		 *
		 * @param Branch_Name  This is the Name of the Branch that we want to add
		 *
		 * @param branch_data  This is the data that we wan to fill the Branch with
		 *
		 * @return Void
		 */
		static void WriteBranch( TTree* input_tree, TString Branch_Name, vector<double>* branch_data );


		//	Functions used internally within this class
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		/*!
		 * @brief Used for setting up a set of strings to describe each combination
		 *
		 * @return void
		 */
		void SetupCombinationDescriptions();

		/*!
		 * @brief Function to determine the normalisation required to map the PDF to the dataset
		 *
		 * @param combinationIndex  This is the combination number which is being Normalised
		 *
		 * @return Retuns the scaling normalisation factor to project the PDF function over data
		 */
		double PDF2DataNormalisation( unsigned int combinationIndex ) const;

		/*!
		 * @brief helper function for testing that the integrator has actually worked
		 *
		 * This is a check that everything worked.
		 *
		 * If the PDF returns nothing but 0 for the interrogated component for the whole Observable Space it will exit.
		 *
		 * Plotting lines of 0 can cause ROOT problems later and it's likely you have a PDF error!
		 *
		 * @param pointValues   This contains the Y values from the projection of this component for this Discrete Combination
		 *
		 * @param component     This contains the component Name
		 *
		 * @return Void
		 */
		void Sanity_Check( vector<double>* pointValues, TString component );



		//	Private member objects of this class
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		string observableName;					/*!	Observable being projected	*/

		double boundary_min, boundary_max;			/*!	Range of Observable extracted from PhaseSpaceBoundary	*/
		double step_size;


		//	Member variables of the Plotter class
		IPDF * plotPDF;						/*!	Copy of the plotPDF passed to the constructor		*/
		IDataSet * plotData;					/*!	Pointer to the whole dataset				*/
		RapidFitIntegrator * pdfIntegrator;			/*!	Copy of the Integrator class passed to the constructor	*/

		vector<double> combination_integral, ratioOfIntegrals;


		//	sWeights
		bool weightsWereUsed;					/*!	True if Weights have been used and have been setup	*/
		string weightName;					/*!	Name of the weight parameter				*/
		double weight_norm;					/*!	Normalisation Value to correct the effect of weighting	*/
		vector<double> wanted_weights;


		EdStyle* format;					/*!	Pointer to an instance of the style			*/

		TString pdfStr;						/*!	Used for seperating the different PDF components	*/
		TFile* PlotFile;					/*!	Whole file where we write to data			*/
		PhaseSpaceBoundary* full_boundary;			/*!	Copy of the PhaseSpaceBoundary passed to the constructor*/



		vector<string> discreteNames, continuousNames;		/*!	These lists are useful for keeping record of the discrete/continous observables	*/



		/*!
		 * Data Subsets which have been presorted
		 * element 0 is whole dataset
		 * element 1 is discrete combination 1
		 * element 2 is discrete combination 2 ... etc
		 */
		vector<vector<DataPoint*> > data_subsets;


		/*!
		 * A vector of DataPoints, each representative of a Unique combinations for this dataset
		 */
		vector<DataPoint*> allCombinations;
		/*!
		 * Weights of each of the unique combinations
		 */
		vector<double> combinationWeights;

		/*!
		 * Weights of each of the unique combinations
		 */
		vector<string> combinationDescriptions;

		/*!
		 * Vectors containing the unique discrete values of each combination
		 */
		vector<vector<double> > observableValues;


		/*!
		 * These are configured from the XML
		 */
		const int total_points;
		const int data_binning;
		const bool logY;
		const bool logX;

		/*!
		 * These are effectivley the output objects which can be requested as long as this class was still alive
		 */

		vector<TGraphErrors*> binned_data;			/*!	Binned data for the various combinations		*/
		vector<vector<TGraph*> > total_components;		/*!	Pointers to the components for each of the combinations */


		double chi2;				/*!	Chi2 value when calculated	*/
		double N;				/*!	N value when calculated		*/


		/*!	Internal pointer to the original configuration sent to the class constructor		*/
		CompPlotter_config* this_config;

		bool onlyZero;				/*!	Should this class mimic the old Plotter behaviour?	*/

		vector<double> allPullData;	/*!	This will store the value of the total PDF evaluated at each bin, lower edge, center, and upper edge	*/

		DebugClass* debug;

		int PDFNum;
};

#endif

