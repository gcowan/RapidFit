/**
  @class OutputConfiguration

  A container for all information configuring RapidFit output

  @author Benjamin Wynne bwynne@cern.ch
 */

#include "OutputConfiguration.h"
#include "ResultFormatter.h"
#include "Plotter.h"
#include "ScanParam.h"
#include <time.h>
#include "ScanParam.h"
#include "TTree.h"
#include "StringProcessing.h"
#include "stdlib.h"

//Default constructor
OutputConfiguration::OutputConfiguration() : 
   pullType("None"), 
   makeAllPlots(false), 
   pullFileName("pullPlots.root"), 
   projectionFileName("projectionPlots.root"),  
   LLscanFileName("LLscans.root"),
   contourFileName("contourPlots.root"),
   weightedEventsWereUsed(false)
{
}

//Constructor with correct arguments
OutputConfiguration::OutputConfiguration( vector< pair< string, string > > InputContours, vector<string> InputProjections, string PullPlotType, vector<ScanParam*> ScanParameters, vector<pair<ScanParam*, ScanParam*> > _2DScanParameters ) : 
    contours(InputContours), 
    projections(InputProjections),
    pullType(PullPlotType),
	makeAllPlots(false), 
    pullFileName("pullPlots.root"), 
    projectionFileName("projectionPlots.root"), 
    LLscanFileName("LLscanPlots.root"),
	LLcontourFileName("LLcontourPlots.root"),
    contourFileName("contourPlots.root"),
    weightedEventsWereUsed(false), 
    Global_Scan_List( ScanParameters ),
    Global_2DScan_List( _2DScanParameters )
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
	return projections;
}

ScanParam* OutputConfiguration::GetScanParam( string param_name )
{

	ScanParam* Returnable_Param=NULL;
	for(unsigned short int i=0; i < Global_Scan_List.size(); i++)
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

pair<ScanParam*, ScanParam*> OutputConfiguration::Get2DScanParams( string param_1, string param_2 )
{
	pair<ScanParam*, ScanParam* > Returnable_Pair;
	for(unsigned short int i=0; i < Global_2DScan_List.size(); i++)
	{
		if( ( Global_2DScan_List[i].first->HasName() ) && ( Global_2DScan_List[i].second->HasName() ) )
		{
			if( ( Global_2DScan_List[i].first->GetName() == param_1 ) && ( Global_2DScan_List[i].second->GetName() == param_2 ) )
			{
				Returnable_Pair = Global_2DScan_List[i];
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
	for(unsigned short int i=0; i < Global_Scan_List.size(); i++)
	{
		ScanReturnList.push_back( Global_Scan_List[i]->GetName() );
	}
	return ScanReturnList;
}

//Return the requested Scans
vector<pair<string, string> > OutputConfiguration::Get2DScanList( )
{
	vector<pair<string, string> > ScanReturnList;
	for(unsigned short int i=0; i < Global_2DScan_List.size(); i++)
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
vector<double> OutputConfiguration::GetRange( string wanted_param )
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
		error = Wanted_Param->GetSigma() * error;
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
pair<vector<double>, vector<double> > OutputConfiguration::Get2DRange( string param_1, string param_2 )
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
	ResultFormatter::LatexOutputFitResult(TheResult);
	ResultFormatter::LatexOutputCovarianceMatrix(TheResult);

	//Output any calculated contours
	ResultFormatter::PlotFitContours( TheResult, contourFileName );
	
	//Make any requested projections
	//for ( int projectionIndex = 0; projectionIndex < projections.size(); projectionIndex++ )
	for (unsigned int projectionIndex = 0; projectionIndex < 1; projectionIndex++ )
	{

		PhysicsBottle * resultBottle = TheResult->GetPhysicsBottle();

		//Loop over all PDFs, and plot
		if ( makeAllPlots || projections.size() > 0 )
		{
			for (int resultIndex = 0; resultIndex < resultBottle->NumberResults(); resultIndex++ )
			{
				Plotter * testPlotter = new Plotter( resultBottle->GetResultPDF(resultIndex), resultBottle->GetResultDataSet(resultIndex) );
				if( weightedEventsWereUsed ) testPlotter->SetWeightsWereUsed( weightName ) ;
				
				char fileNumber[100];
				sprintf( fileNumber, "fit%d.", resultIndex );

				if (makeAllPlots)
				{
					testPlotter->PlotAllObservables( fileNumber + projectionFileName );
				}
				else
				{
					testPlotter->PlotObservables( fileNumber + projectionFileName, projections );
				}
			}
		}
	}
}


//Make the requested output from a single result
void OutputConfiguration::OutputLLcontourResult( vector<LLscanResult2D*> scanResults )
{
	ResultFormatter::MakeLLcontourPlots( scanResults, LLcontourFileName );
}


//Make the requested output from a single result
void OutputConfiguration::OutputLLscanResult( vector<LLscanResult*> scanResults )
{
	ResultFormatter::MakeLLscanPlots( scanResults, LLscanFileName );
}


//Make the requested output from a toy study
void OutputConfiguration::OutputToyResult( ToyStudyResult * TheResult )
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

void OutputConfiguration::AddContour( string X_axis, string Y_axis )
{
	vector<string> Contour_X_Vals = StringProcessing::SplitString( X_axis, ',' );
	vector<string> Contour_Y_Vals = StringProcessing::SplitString( Y_axis, ',' );

	if( Contour_X_Vals.size() != 4 ){
		cerr << "Runtime Defined Contour Badly Defined:" << endl;
		for(unsigned short int i=0; i < Contour_X_Vals.size(); i++ )
		{
			cerr << Contour_X_Vals[i] << "\t";
		}
		cerr << endl;
		exit(13);
	}
	if( Contour_Y_Vals.size() != 4 ){
		cerr << "Runtime Defined Contour Badly Defined:" << endl;
		for(unsigned short int i=0; i < Contour_Y_Vals.size(); i++ )
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

void OutputConfiguration::AddScan( string X_axis )
{
	vector<string> Contour_X_Vals = StringProcessing::SplitString( X_axis, ',' );

	if( Contour_X_Vals.size() != 4 ){
		cerr << "Runtime Defined Contour Badly Defined:" << endl;
		for(unsigned short int i=0; i < Contour_X_Vals.size(); i++ )
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
