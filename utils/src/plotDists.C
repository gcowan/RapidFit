
#include "TError.h"
#include "TTree.h"
#include "TH1.h"
#include "TString.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TAxis.h"

#include "ROOT_File_Processing.h"
#include "TTree_Processing.h"
#include "Template_Functions.h"
#include "StringOperations.h"

#include "EdStyle.h"
#include "StringProcessing.h"

#include <sstream>
#include <iomanip>
#include <string>
#include <limits>
#include <cfloat>
#include <cstdlib>

using namespace::std;

struct PlottingConfig
{
	explicit PlottingConfig() :
		num_bins(100), inputFileName(""), inputTreeName(""), requestedParameters(),
		explicit_min(-DBL_MAX), explicit_max(DBL_MAX), cutString(""), drawOptions( "PE9" ), weightName(""),
		setLogY(false), setLogX(false), y_min(-DBL_MAX), y_max(DBL_MAX)
	{}
	int num_bins;
	string inputFileName;
	string inputTreeName;
	vector<string> requestedParameters;
	double explicit_min;
	double explicit_max;
	string cutString;
	string drawOptions;
	string weightName;
	bool setLogY;
	bool setLogX;
	double y_min;
	double y_max;
};

string ProcessThisInputOption( unsigned int& counter, int argc, char* argv[], string ErrorMsg )
{
	if( (counter+1) < (unsigned) argc )
	{
		++counter;
		return string( argv[(unsigned)counter] );
	}
	else
	{
		cerr << ErrorMsg << endl;
		exit(0);
	}
}

void ProcessInputOptions( struct PlottingConfig& configToModify, int argc, char* argv[] )
{
	if( argc == 1 )
	{
		// ...
	}
	else if( argc == 2 )
	{
		configToModify.inputFileName = argv[1];
		configToModify.inputTreeName = "";
	}
	else
	{
		for( unsigned int i=1; i< (unsigned) argc; ++i )
		{
			string thisArg = argv[i];
			//cout << thisArg << "  " << i << " of " << argc << endl;
			if( thisArg == "-f" )
			{
				configToModify.inputFileName = ProcessThisInputOption( i, argc, argv, "Missing FileName" );
			}
			else if( thisArg == "--DrawOptions" )
			{
				configToModify.drawOptions = ProcessThisInputOption( i, argc, argv, "Missing Draw Options" );
			}
			else if( thisArg == "--SetBins" )
			{
				configToModify.num_bins = atoi( ProcessThisInputOption( i, argc, argv, "Missing Bin Number" ).c_str() );
			}
			else if( thisArg == "--SetLogY" )
			{
				configToModify.setLogY = true;
			}
			else if( thisArg == "--SetLogX" )
			{
				configToModify.setLogX = true;
			}
			else if( thisArg == "--Parameters" )
			{
				configToModify.requestedParameters = StringProcessing::SplitString( ProcessThisInputOption( i, argc, argv, "Missing Parameter List" ), ':' );
			}
			else if( thisArg == "--WeightName" )
			{
				configToModify.weightName = ProcessThisInputOption( i, argc, argv, "Missing Weight Name" );
			}
			else if( thisArg == "--SetRangeMin" )
			{
				configToModify.explicit_min = atof( ProcessThisInputOption( i, argc, argv, "Missing Minimum" ).c_str() );
			}
			else if( thisArg == "--setRangeMax" )
			{
				configToModify.explicit_max = atof( ProcessThisInputOption( i, argc, argv, "Missing Maximum" ).c_str() );
			}
			else if( thisArg == "--setMax" )
			{
				configToModify.y_max = atof( ProcessThisInputOption( i, argc, argv, "Missing Minimum" ).c_str() );
			}
			else if( thisArg == "--setMin" )
			{
				configToModify.y_min = atof( ProcessThisInputOption( i, argc, argv, "Missing Maximum" ).c_str() );
			}
			else if( thisArg == "--Cut" )
			{
				configToModify.cutString = ProcessThisInputOption( i, argc, argv, "Missing Cut String" );
			}
		}
		// ....
	}
}

void SafeFill( TH1* thisTH1, double thisValue, double thisWeight, double thisMin, double thisMax )
{
	if( thisValue > thisMin && thisValue > thisMin )
	{
		thisTH1->Fill( thisValue, thisWeight );
	}
	//      Yes this code 'should' do nothing, but ROOT INSISTS on having issues with events sitting on bin-edges...
	//      this is a severe irritant when plotting things for papers so I fix it at some level here...
	else if( thisValue <= thisMin )
	{
		if( fabs( thisValue - thisMin ) < 1E-6 )
		{
			thisTH1->Fill( thisMin+1E-6, thisWeight );
		}
	}
	else if( thisValue >= thisMax )
	{
		if( fabs( thisValue - thisMax ) < 1e-6 )
		{
			thisTH1->Fill( thisMax-1E-6, thisWeight );
		}
	}
}

int main( int argc, char* argv[] )
{
	EdStyle* thisStyle = new EdStyle();

	thisStyle->SetStyle();

	gStyle->SetOptStat(0);

        //      Mute ROOT
	gErrorIgnoreLevel = kFatal;

	struct PlottingConfig thisConfig = PlottingConfig();

	ProcessInputOptions( thisConfig, argc, argv );

	if( !thisConfig.cutString.empty() && !thisConfig.weightName.empty() )
	{
		cout << "WARNING!!! You have made a Cut and are making a weighted plot!!!" << endl;
		cout << "Be Careful!!!" << endl;
	}

	TTree* inputTree = ROOT_File_Processing::GetFirstTree( thisConfig.inputFileName, "READ" );

	if( thisConfig.requestedParameters.empty() )
	{
		thisConfig.requestedParameters =  StringOperations::TString2string( TTree_Processing::get_branch_names( inputTree ) );
		if( !thisConfig.weightName.empty() )
		{
			vector<string> allParameters;
			for( unsigned int i=0; i< thisConfig.requestedParameters.size(); ++i )
			{
				if( thisConfig.requestedParameters[i] != thisConfig.weightName ) allParameters.push_back( thisConfig.requestedParameters[i] );
			}
			thisConfig.requestedParameters = allParameters;
		}
	}

	TString plotDistsOutput = "PlotDists_Output";

	cout << "Output Plots Sotred in: " << plotDistsOutput << endl;

	gSystem->mkdir( plotDistsOutput );

	vector<Double_t>* weightData = NULL;

	if( !thisConfig.weightName.empty() )
	{
		weightData = TTree_Processing::Buffer_Branch( inputTree, thisConfig.weightName, thisConfig.cutString );
	}

	for( unsigned int i=0; i< thisConfig.requestedParameters.size(); ++i )
	{
		string thisParameterName = thisConfig.requestedParameters[i];

		vector<Double_t>* thisData = TTree_Processing::Buffer_Branch( inputTree, thisParameterName, thisConfig.cutString );

		TString CanvasName = "Canvas_"; CanvasName.Append( thisParameterName );

		TCanvas* thisCanvas = EdStyle::RapidFitCanvas( CanvasName );

		thisCanvas->SetLogy( thisConfig.setLogY );
		thisCanvas->SetLogx( thisConfig.setLogX );

		TString TH1Name = "TH1_"; TH1Name.Append( thisParameterName );

		Double_t thisMin=0.;
		Double_t thisMax=0.;

		if( thisConfig.explicit_min <= -DBL_MAX )	thisMin = get_minimum( *thisData );
		else						thisMin = thisConfig.explicit_min;

		if( thisConfig.explicit_max >= DBL_MAX )	thisMax = get_maximum( *thisData );
		else						thisMax = thisConfig.explicit_max;

		TH1* thisTH1 = new TH1D( TH1Name, TH1Name, thisConfig.num_bins, thisMin, thisMax );

		for( unsigned int j=0; j< thisData->size(); ++j )
		{
			double weight=0.;
			if( weightData == NULL )	weight = 1.;
			else				weight = weightData->at(j);

			SafeFill( thisTH1, thisData->at(j), weight, thisMin, thisMax );

		}

		thisTH1->Draw( thisConfig.drawOptions.c_str() );
		thisCanvas->Update();

		thisTH1->GetXaxis()->SetRangeUser( thisMin, thisMax );
		thisCanvas->Update();

		TString XaxisLabel = EdStyle::GetParamRootName( thisParameterName );
		XaxisLabel.Append( " " );
		XaxisLabel.Append( EdStyle::GetParamRootUnit( thisParameterName ) );

		thisTH1->GetXaxis()->SetTitle( XaxisLabel );

		TString YaxisLabel = "Candidates / (";
		stringstream thisStream;
		thisStream << setprecision(2) << scientific << fabs(thisMax-thisMin)/((double)thisConfig.num_bins);
		YaxisLabel.Append( thisStream.str().c_str() );
		YaxisLabel.Append( ")" );

		thisTH1->GetYaxis()->SetTitle( YaxisLabel );

		thisCanvas->Update();

		double ymin = thisTH1->GetBinContent( thisTH1->GetMinimumBin() ) - thisTH1->GetBinError( thisTH1->GetMinimumBin() );
		double ymax = thisTH1->GetBinContent( thisTH1->GetMaximumBin() ) + thisTH1->GetBinError( thisTH1->GetMaximumBin() );

		if( thisConfig.y_min > -DBL_MAX )		ymin = thisConfig.y_min;
		if( thisConfig.y_max < DBL_MAX )		ymax = thisConfig.y_max;

		thisTH1->GetYaxis()->SetRangeUser( ymin, ymax );
		thisTH1->SetMinimum( ymin );
		thisTH1->SetMaximum( ymax );

		thisCanvas->Update();

		TString localPDFName = thisParameterName.c_str();

		localPDFName.Append( "_Output.pdf" );

		TString TotalPDFName = plotDistsOutput; TotalPDFName.Append("/");
		TotalPDFName.Append( localPDFName );

		thisCanvas->Print( TotalPDFName );

		delete thisData;
	}

}

