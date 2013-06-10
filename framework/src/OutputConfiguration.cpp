/**
  @class OutputConfiguration

  A container for all information configuring RapidFit output

  @author Benjamin Wynne bwynne@cern.ch
  */

//	ROOT Headers
#include "TSystem.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
//	RapidFit Headers
#include "OutputConfiguration.h"
#include "ResultFormatter.h"
#include "ComponentPlotter.h"
#include "ScanParam.h"
#include "StringProcessing.h"
//	System Headers
#include <time.h>
#include <stdlib.h>
#include <iomanip>

using namespace::std;

//Constructor with correct arguments
OutputConfiguration::OutputConfiguration( vector< pair< string, string > > InputContours, string PullPlotType, vector<ScanParam*> ScanParameters, vector<pair<ScanParam*, ScanParam*> > _2DScanParameters, vector<CompPlotter_config*> Components ) :
	contours(InputContours),
	LLscanList(),
	pullType(PullPlotType),
	makeAllPlots(false),
	weightName(),
	pullFileName("pullPlots.root"),
	projectionFileName("projectionPlots.root"),
	LLscanFileName("LLscanPlots.root"),
	LLcontourFileName("LLcontourPlots.root"),
	contourFileName("contourPlots.root"),
	weightedEventsWereUsed(false),
	Global_Scan_List( ScanParameters ),
	Global_2DScan_List( _2DScanParameters ),
	Stored_Fit_Results(),
	proj_components( Components ),
	debug(new DebugClass() )
{
}

//Destructor
OutputConfiguration::~OutputConfiguration()
{
}

//Return the requested contour plots
vector< pair< string, string > > OutputConfiguration::GetContourPlots()
{
	return contours;
}

//Return the requested projections
vector<string> OutputConfiguration::GetProjections()
{
	vector<string> projections;
	for( unsigned int i=0; i< proj_components.size(); ++i )
	{
		projections.push_back( proj_components[i]->observableName );
	}
	return projections;
}

ScanParam* OutputConfiguration::GetScanParam( string param_name )
{

	ScanParam* Returnable_Param=NULL;
	for(unsigned short int i=0; i < Global_Scan_List.size(); ++i)
	{
		if( Global_Scan_List[i]->HasName() )
		{
			if( Global_Scan_List[i]->GetName() == param_name )
			{
				Returnable_Param = Global_Scan_List[i];
			}
		}
	}
	if( Returnable_Param->HasSigma() )
	{
		vector<double> range = GetRange( Returnable_Param );
		Returnable_Param->SetMax( range[0] );
		Returnable_Param->SetMin( range[1] );
	}
	return Returnable_Param;
}

pair<ScanParam*, ScanParam*> OutputConfiguration::Get2DScanParams( const string param_1, const string param_2 )
{
	pair<ScanParam*, ScanParam* > Returnable_Pair;
	for( int i=int(Global_2DScan_List.size()-1); i >= 0 ; --i )
	{
		if( ( Global_2DScan_List[unsigned(i)].first->HasName() ) && ( Global_2DScan_List[unsigned(i)].second->HasName() ) )
		{
			if( ( Global_2DScan_List[unsigned(i)].first->GetName() == param_1 ) && ( Global_2DScan_List[unsigned(i)].second->GetName() == param_2 ) )
			{
				Returnable_Pair = Global_2DScan_List[unsigned(i)];
			}
		}
	}

	if( Returnable_Pair.first->HasSigma() || Returnable_Pair.second->HasSigma() )
	{
		vector<double> range_1 = GetRange( Returnable_Pair.first );
		vector<double> range_2 = GetRange( Returnable_Pair.second );

		Returnable_Pair.first->SetMax( range_1[0] );
		Returnable_Pair.first->SetMin( range_1[1] );
		Returnable_Pair.second->SetMax( range_2[0] );
		Returnable_Pair.second->SetMin( range_2[1] );
	}

	return Returnable_Pair;
}

//Return the requested Scans
vector<string> OutputConfiguration::GetScanList( )
{
	vector<string> ScanReturnList;
	for(unsigned short int i=0; i < Global_Scan_List.size(); ++i)
	{
		ScanReturnList.push_back( Global_Scan_List[i]->GetName() );
	}
	return ScanReturnList;
}

void OutputConfiguration::ClearScanList()
{
	while( !Global_Scan_List.empty() )
	{
		if( Global_Scan_List.back() != NULL ) delete Global_Scan_List.back();
		Global_Scan_List.pop_back();
	}
	return;
}

void OutputConfiguration::Clear2DScanList()
{
	while( !Global_2DScan_List.empty() )
	{
		if( Global_2DScan_List.back().first != NULL ) delete Global_2DScan_List.back().first;
		if( Global_2DScan_List.back().second != NULL ) delete Global_2DScan_List.back().second;
		Global_2DScan_List.pop_back();
	}
	return;
}

//Return the requested Scans
vector<pair<string, string> > OutputConfiguration::Get2DScanList( )
{
	vector<pair<string, string> > ScanReturnList;
	for(unsigned short int i=0; i < Global_2DScan_List.size(); ++i)
	{
		string new_first = Global_2DScan_List[i].first->GetName();
		string new_second = Global_2DScan_List[i].second->GetName();
		pair<string,string> new_pair;
		new_pair.first = new_first;
		new_pair.second = new_second;
		ScanReturnList.push_back( new_pair );
	}
	return ScanReturnList;
}

//  [0] Max, [1] Min, [2] nPoints
vector<double> OutputConfiguration::GetRange( const string wanted_param )
	//  Return a vector containing [0]=Max, [1]=Min and [2]=Points
{
	ScanParam* Wanted_Param = GetScanParam( wanted_param );
	return GetRange( Wanted_Param );
}

//  [0] Max, [1] Min, [2] nPoints
vector<double> OutputConfiguration::GetRange( ScanParam* Wanted_Param )
{
	//  [0] Max, [1] Min, [2] nPoints
	vector<double> ReturnRange;
	string wanted_param = Wanted_Param->GetName();

	if( !Stored_Fit_Results.empty() && Wanted_Param->HasSigma() )
	{

		double error = Stored_Fit_Results[0]->GetResultParameter( wanted_param )->GetError();
		error = Wanted_Param->GetSigma() * Wanted_Param->GetSigma() * error;
		double maximum = Stored_Fit_Results[0]->GetResultParameter( wanted_param )->GetValue() + error;
		double minimum = Stored_Fit_Results[0]->GetResultParameter( wanted_param )->GetValue() - error;
		ReturnRange.push_back( maximum );
		ReturnRange.push_back( minimum );

		if( Wanted_Param->HasPoints() )  {
			ReturnRange.push_back( Wanted_Param->GetPoints() );
		} else {	ReturnRange.push_back( 5 ); }

		return ReturnRange;
	} else if( Wanted_Param->HasMax() && Wanted_Param->HasMin() ) {
		ReturnRange.push_back( Wanted_Param->GetMax() );
		ReturnRange.push_back( Wanted_Param->GetMin() );
		if( Wanted_Param->HasPoints() )  {
			ReturnRange.push_back( Wanted_Param->GetPoints() );
		} else {	ReturnRange.push_back( 5 ); }

		return ReturnRange;
	} else {
		cerr << "Unable to Properly Determine the Range for : " << wanted_param << " Giving Stupid Limits that are safe" << endl;
		ReturnRange.push_back( 1.0 );
		ReturnRange.push_back( -1.0 );
		ReturnRange.push_back( 5 );
		return ReturnRange;
	}
}

//	Pair of:  [0] Max, [1] Min, [2] nPoints
pair<vector<double>, vector<double> > OutputConfiguration::Get2DRange( const string param_1, const string param_2 )
{
	pair<vector<double>, vector<double> >  Returnable_Range;

	pair<ScanParam*, ScanParam*> Param_Set = Get2DScanParams( param_1, param_2 );

	vector<double> Range1 = GetRange( Param_Set.first  );
	vector<double> Range2 = GetRange( Param_Set.second );

	Returnable_Range.first = Range1;
	Returnable_Range.second = Range2;

	return Returnable_Range;
}

//Return whether to do pull plots
bool OutputConfiguration::DoPullPlots()
{
	return !( pullType == "None" );
}

//Make the requested output from a single result
void OutputConfiguration::OutputFitResult( FitResult * TheResult )
{
	//Output information aboout the fit
	//ResultFormatter::LatexOutputFitResult(TheResult);
	//ResultFormatter::WriteOutputLatex(TheResult);
	//ResultFormatter::LatexOutputCovarianceMatrix(TheResult);

	ResultFormatter::WriteOutputLatex(TheResult);

	//Output any calculated contours
	ResultFormatter::PlotFitContours( TheResult, contourFileName );


	//	I spent long enough writing componentprojections that these will take priority over the old projection code

	if( !proj_components.empty() ) this->OutputCompProjections( TheResult );

	return;
}

void OutputConfiguration::OutputCompProjections( FitResult* TheResult )
{
	PhysicsBottle* resultBottle = TheResult->GetPhysicsBottle();

	string outputFolderName = ResultFormatter::GetOutputFolder();
	if( outputFolderName.empty() )
	{
		ResultFormatter::initOutputFolder();
		ResultFormatter::GetOutputFolder();
	}

	gSystem->cd( outputFolderName.c_str() );

	for( vector<CompPlotter_config*>::iterator projection_i = proj_components.begin(); projection_i != proj_components.end(); ++projection_i )
	{
		vector<vector<TGraph*> > all_components_for_all_results;

		vector<TGraphErrors*> all_datasets_for_all_results;

		vector<ComponentPlotter*> allComponentPlotters;

		TString folderName( "RapidFit_Component_Projection_" ); folderName.Append( (*projection_i)->observableName ); folderName.Append( "_" );
		folderName.Append( StringProcessing::TimeString() );

		gSystem->mkdir( folderName );
		gSystem->cd( folderName );

		TString filename( "RapidFit_Component_Projections_" ); filename.Append( (*projection_i)->observableName ); filename.Append( "_" );
		filename.Append( StringProcessing::TimeString() );
		filename.Append( ".root" );
		TFile* output_file = new TFile( filename, "UPDATE" );

		vector< pair<double,double> > chi2_results;

		vector< vector<double> > pullFunctionEvals;

		for( int resultIndex = 0; resultIndex < resultBottle->NumberResults(); ++resultIndex )
		{
			string thisObservable = (*projection_i)->observableName;

			vector<string> known_observables = resultBottle->GetResultPDF(resultIndex)->GetPrototypeDataPoint();
			int num = StringProcessing::VectorContains( &known_observables, &thisObservable );

			vector<string> bad_observables = resultBottle->GetResultPDF(resultIndex)->GetDoNotIntegrateList();
			int num2 = StringProcessing::VectorContains( &bad_observables, &thisObservable );

			if( num2 != -1 )
			{
				cerr << "Observable: " << thisObservable << " is on Do Not Integrate List, cannot perform projection!" << endl;
				continue;
			}
			if( num == -1 )
			{
				cerr << "Observable: " << thisObservable << " is not constrained by this PDF!" << endl;
				continue;
			}

			cout << "Projecting ToFit: " << resultIndex+1 << endl << endl;
			TString PDFStr = "PDF_";PDFStr+=resultIndex;

			if( debug != NULL )
			{
				if( debug->DebugThisClass( "OutputConfiguration" ) )
				{
					cout << "OutputConfiguration: Constructing ComponentPlotter: " << resultIndex+1 << " of " << resultBottle->NumberResults() << endl;
				}
			}

			//	ComponentPlotter requires a PDF, Dataset, output_file, Observable to project, a plot configuration object and a string for the path for where the output for this PDF belongs
			ComponentPlotter* thisPlotter = new ComponentPlotter( resultBottle->GetResultPDF(resultIndex), resultBottle->GetResultDataSet(resultIndex),

					PDFStr, output_file, thisObservable, (*projection_i), resultIndex, debug );


			thisPlotter->SetDebug( debug );

			if( weightedEventsWereUsed ) thisPlotter->SetWeightsWereUsed( weightName );

			allComponentPlotters.push_back( thisPlotter );

			if( debug != NULL )
                        {
                                if( debug->DebugThisClass( "OutputConfiguration" ) )
                                {
                                        cout << "OutputConfiguration: Requesting Projection: " << resultIndex+1 << " of " << resultBottle->NumberResults() << endl;
                                }
                        }

			//	In ComponentPlotter this still does all of the work, but each ComponentPlotter object is created for each observable to allow you to do more easily
			thisPlotter->ProjectObservable();

			cout << "Projected" << endl;

			if( (*projection_i)->CalcChi2 ) chi2_results.push_back( thisPlotter->GetChi2Numbers() );

			pullFunctionEvals.push_back( thisPlotter->GetFunctionEval() );

			all_datasets_for_all_results.push_back( thisPlotter->GetBinnedData()[0] );
			all_components_for_all_results.push_back( thisPlotter->GetComponents()[0] );
		}

		if( !chi2_results.empty() )
		{
			double total_chi2 = 0.;
			double total_N = 0.;
			for( unsigned int i=0; i< chi2_results.size(); ++i )
			{
				total_chi2 += chi2_results[i].first;
				total_N += chi2_results[i].second;
			}
			double final_corrected_chi2 = total_chi2 / ( total_N - (double)resultBottle->GetResultPDF(0)->GetPhysicsParameters()->GetAllFloatNames().size() - 1. );
			(*projection_i)->Chi2Value = final_corrected_chi2;

			cout << endl << "\tFinal chi2/ndof = " << setprecision(10) << final_corrected_chi2 << endl << endl;
		}

		if( all_components_for_all_results.empty() )
		{
			output_file->Close();
			remove( filename.Data() );
			gSystem->cd( ".." );
			rmdir( folderName.Data() );
			continue;
		}
		int num=(int) all_components_for_all_results[0].size();
		bool compatible=true;

		for( unsigned int i=0; i< all_components_for_all_results.size(); ++i )
		{
			if( (int)all_components_for_all_results[i].size() != num ) compatible = false;
		}

		TGraphErrors* Total_BinnedData = NULL;

		vector<TGraph*> Total_Components;

		//	We only want to overlay the output when all of the PDFs have the same number of components
		//	Eg it is difficult to justify trying to overlay data from JpsiPhi and Jpsif0 in one component plot, although I'm sure it's possible
		if( compatible )
		{

			//	These static functions will take the global Random function from ROOT if need be, but lets not force that issue and give it an already initialized generator
			//	This, alas is a ROOT solution to a ROOT problem, and it creates a MESS that often is left to some code outside of ROOT to cleanup!
			Total_BinnedData = ComponentPlotter::MergeBinnedData( all_datasets_for_all_results, resultBottle->GetResultPDF(0)->GetRandomFunction() );

			Total_Components = ComponentPlotter::MergeComponents( all_components_for_all_results, resultBottle->GetResultPDF(0)->GetRandomFunction() );

			output_file->cd();
			if( output_file->GetDirectory( "Total" ) == 0 ) output_file->mkdir( "Total" );
			output_file->cd( "Total" );

			ComponentPlotter::WriteData( Total_BinnedData, Total_Components, "Final_ProjectionData" );

			ComponentPlotter::OutputPlot( Total_BinnedData, Total_Components, (*projection_i)->observableName, "_All_Data",

					resultBottle->GetResultDataSet(0)->GetBoundary(), resultBottle->GetResultPDF(0)->GetRandomFunction(), (*projection_i), debug );
		}

		vector<double> finalPullEvals;

		if( (Total_BinnedData != NULL) && ((*projection_i)->DrawPull == true) )
		{
			for( unsigned int i=0; i< pullFunctionEvals[0].size(); ++i )		//	Loop over all Bins
			{
				double this_bin=0.;
				for( unsigned int j=0; j< pullFunctionEvals.size(); ++j )	//	Loop over all PDFs (all ToFits)
				{
					this_bin+= pullFunctionEvals[j][i];
				}
				finalPullEvals.push_back( this_bin );
			}

			TGraphErrors* pullGraph = ComponentPlotter::PullPlot1D( finalPullEvals, Total_BinnedData, (*projection_i)->observableName, "_PullPlot", resultBottle->GetResultPDF(0)->GetRandomFunction() );
			(void) pullGraph;

			ComponentPlotter::OutputPlot( Total_BinnedData, Total_Components, (*projection_i)->observableName, "_All_Data_wPulls", resultBottle->GetResultDataSet(0)->GetBoundary(),
								resultBottle->GetResultPDF(0)->GetRandomFunction(), (*projection_i), debug, finalPullEvals );
		}

		CompPlotter_config* datasets_config = new CompPlotter_config( *(*projection_i) );
		vector<TGraph*> allDataSubSets;
		allDataSubSets.push_back( Total_Components.back() );
		vector<string> datasetID; datasetID.push_back( "All Datasets" );
		for( unsigned int i=0; i< all_components_for_all_results.size(); ++i )
		{
			allDataSubSets.push_back( all_components_for_all_results[i].back() );
			allDataSubSets.back()->SetLineColor( (Color_t)(i+2) );
			datasetID.push_back( resultBottle->GetResultPDF( i )->GetLabel() );
		}

		datasets_config->component_names = datasetID;
		datasets_config->LegendTextSize = (Size_t)0.02;

		ComponentPlotter::OutputPlot( Total_BinnedData, allDataSubSets, (*projection_i)->observableName, "_All_SubSets",
						resultBottle->GetResultDataSet(0)->GetBoundary(), resultBottle->GetResultPDF(0)->GetRandomFunction(), datasets_config, debug );

		if( !finalPullEvals.empty() )
		{
			ComponentPlotter::OutputPlot( Total_BinnedData, allDataSubSets, (*projection_i)->observableName, "_All_SubSets_wPulls",
						resultBottle->GetResultDataSet(0)->GetBoundary(), resultBottle->GetResultPDF(0)->GetRandomFunction(), datasets_config, debug, finalPullEvals );
		}

		delete datasets_config;

		//	it is now safe to remove the instances of ComponentPlotter
		for( vector<ComponentPlotter*>::iterator compP_i = allComponentPlotters.begin(); compP_i != allComponentPlotters.end(); ++compP_i )
		{
			if( (*compP_i) != NULL ) delete (*compP_i);
		}

		//	From experience I tend not to call delete on ROT objects as it causes segfaults when ROOT starts to clean up after itself
		output_file->Close();

		gSystem->cd( ".." );

	}

	gSystem->cd( ".." );
}

//Make the requested output from a toy study
void OutputConfiguration::OutputToyResult( FitResultVector* TheResult )
{
	ResultFormatter::MakePullPlots( pullType, pullFileName, TheResult );
}

//Make all possible function projections, and save them in given location
void OutputConfiguration::MakeAllPlots( string FileName )
{
	makeAllPlots = true;
	projectionFileName = FileName;
}

//Set the location to store contour plots
void OutputConfiguration::SetContourFileName( string FileName )
{
	contourFileName = FileName;
}

//Set the location to store pull plots
void OutputConfiguration::SetPullFileName( string FileName )
{
	pullFileName = FileName;
}

//Set the location to store LLscans plots
void OutputConfiguration::SetLLscanFileName( string FileName )
{
	LLscanFileName = FileName;
}

//Set the location to store LLcontours plots
void OutputConfiguration::SetLLcontourFileName( string FileName )
{
	LLcontourFileName = FileName;
}

void OutputConfiguration::AddContour( const string X_axis, const string Y_axis )
{
	vector<string> Contour_X_Vals = StringProcessing::SplitString( X_axis, ',' );
	vector<string> Contour_Y_Vals = StringProcessing::SplitString( Y_axis, ',' );

	if( Contour_X_Vals.size() != 4 ){
		cerr << "Runtime Defined Contour Badly Defined:" << endl;
		for(unsigned short int i=0; i < Contour_X_Vals.size(); ++i )
		{
			cerr << Contour_X_Vals[i] << "\t";
		}
		cerr << endl;
		exit(13);
	}
	if( Contour_Y_Vals.size() != 4 ){
		cerr << "Runtime Defined Contour Badly Defined:" << endl;
		for(unsigned short int i=0; i < Contour_Y_Vals.size(); ++i )
		{
			cerr << Contour_Y_Vals[i] << "\t";
		}
		cerr << endl;
		exit(13);
	}

	string x_name = Contour_X_Vals[0];
	double x_max = strtod( Contour_X_Vals[1].c_str(), NULL );
	double x_min = strtod( Contour_X_Vals[2].c_str(), NULL );
	int x_points = atoi( Contour_X_Vals[3].c_str() );
	ScanParam* new_X = new ScanParam( x_name, x_max, x_min, x_points );
	string y_name = Contour_Y_Vals[0];
	double y_max = strtod( Contour_Y_Vals[1].c_str(), NULL );
	double y_min = strtod( Contour_Y_Vals[2].c_str(), NULL );
	int y_points = atoi( Contour_Y_Vals[3].c_str() );
	ScanParam* new_Y = new ScanParam( y_name, y_max, y_min, y_points );

	pair<ScanParam*, ScanParam*> new_Contour;
	new_Contour.first = new_X;
	new_Contour.second = new_Y;

	Global_2DScan_List.push_back( new_Contour );

}

void OutputConfiguration::AddScan( const string X_axis )
{
	vector<string> Contour_X_Vals = StringProcessing::SplitString( X_axis, ',' );

	if( Contour_X_Vals.size() != 4 ){
		cerr << "Runtime Defined Contour Badly Defined:" << endl;
		for(unsigned short int i=0; i < Contour_X_Vals.size(); ++i )
		{
			cerr << Contour_X_Vals[i] << "\t";
		}
		cerr << endl;
		exit(13);
	}

	string x_name = Contour_X_Vals[0];
	double x_max = strtod( Contour_X_Vals[1].c_str(), NULL );
	double x_min = strtod( Contour_X_Vals[2].c_str(), NULL );
	int x_points = atoi( Contour_X_Vals[3].c_str() );
	ScanParam* new_X = new ScanParam( x_name, x_max, x_min, x_points );

	Global_Scan_List.push_back( new_X );

}

//Setso that it knows that weighted events were used
void OutputConfiguration::SetWeightsWereUsed( string _weightName )
{
	weightedEventsWereUsed = true ;
	weightName = _weightName ;
}

void OutputConfiguration::SetInputResults( ResultParameterSet* oneResult )
{
	Stored_Fit_Results.push_back( oneResult );
}

void OutputConfiguration::SetDebug( DebugClass* input_debug )
{
	if( debug != NULL ) delete debug;
	debug = new DebugClass( *input_debug );
}

string OutputConfiguration::GetPullFileName() const
{
	return pullFileName;
}

