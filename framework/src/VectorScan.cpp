//	RapidFit Headers
#include "ScanParam.h"
#include "FitResultVector.h"
#include "FitAssembler.h"
#include "MinimiserConfiguration.h"
#include "VectorScan.h"
//	System Haeaders
#include <vector>
#include <string>

using namespace::std;


//	Simple function for assembling a 1D scan based on user input and calling the fitting Engine to perform the scan
FitResultVector* VectorScan::Scan1D( MinimiserConfiguration* theMinimiser, ScanParam* local_param )
{

	double maximum = local_param->GetMax();
	double minimum = local_param->GetMin();
	int points = local_param->GetPoints();
	double step = ( minimum - maximum )/double(points);

	vector<vector<double> > new_coordinates;

	for( int i=0; i< points; ++i )
	{
		vector<double> new_point;
		new_point.push_back( minimum + step*double(i) );
		new_coordinates.push_back( new_point );
		cout << new_coordinates.back().back() << endl;
	}

	vector<string> new_param;
	new_param.push_back( local_param->GetName() );

	FitResultVector* scan_result = VectorScanEngine( theMinimiser, new_coordinates, new_param );

	return scan_result;
}

//      Simple function for assembling a 2D scan based on user input and calling the fitting Engine to perform the scan
FitResultVector* VectorScan::Scan2D( MinimiserConfiguration* theMinimiser, pair<ScanParam*, ScanParam*> Param_Set )
{

	vector<string> new_param_names;
	new_param_names.push_back( Param_Set.first->GetName() );
	new_param_names.push_back( Param_Set.second->GetName() );

	double maximum_1 = Param_Set.first->GetMax();
	double minimum_1 = Param_Set.first->GetMin();
	int points_1 = Param_Set.first->GetPoints();
	double step_1 = ( maximum_1 - minimum_1 )/double(points_1);

	double maximum_2 = Param_Set.second->GetMax();
	double minimum_2 = Param_Set.second->GetMin();
	int points_2 = Param_Set.second->GetPoints();
	double step_2 = ( maximum_2 - minimum_2 )/double(points_2);

	vector<vector<double> > full_coordinates;

	for( int i=0; i< points_1; ++i )
	{
		double value_1 = minimum_1+double(i)*step_1;

		for( int j=0; j< points_2; ++j )
		{
			cout << step_2 << endl;
			double value_2 = minimum_2+double(j)*step_2;

			vector<double> new_coordinate;
			new_coordinate.push_back( value_1 );
			new_coordinate.push_back( value_2 );
			cout << value_1 << "\t" << value_2 << endl;
			full_coordinates.push_back( new_coordinate );
		}
	}

	FitResultVector* scan_result = VectorScanEngine( theMinimiser, full_coordinates, new_param_names );

	return scan_result;
}


//	Engine for performing fits across a consecutive vector of fixed corrdinates
FitResultVector* VectorScan::VectorScanEngine( MinimiserConfiguration * MinimiserConfig, vector<vector<double> > thePoints, vector<string> theParameters )
{
        IMinimiser* Minimiser = MinimiserConfig->GetMinimiser();
        vector<PhysicsParameter*> scan_parameters;
        vector<FitResultVector*> the_fitresults;

        vector<string> allNames = Minimiser->GetFitFunction()->GetParameterSet()->GetAllNames();

        for( unsigned int point_i=0; point_i< thePoints.size(); ++point_i )
        {
                cout << "Scan Number:\t" << point_i+1 << "\tof:\t" << thePoints.size() << endl;
                if( thePoints[point_i].size() != theParameters.size() )
                {
                        cerr << "Bad Coordinate passed to _DoScan" << endl;
                }

                for( unsigned int j=0; j< theParameters.size(); ++j )
                {
                        cout << thePoints[point_i][j] << endl;
                }

                Minimiser->FixParameters( thePoints[point_i], theParameters );

                FitAssembler::SafeMinimise( Minimiser );

                FitResultVector* new_resultvec = new FitResultVector( allNames );

                new_resultvec->StartStopwatch();

                FitResult* new_result = Minimiser->GetFitResult();

                new_resultvec->AddFitResult( new_result );

                the_fitresults.push_back( new_resultvec );
        }

        FitResultVector* new_result_vector = new FitResultVector( the_fitresults );

        return new_result_vector;
}
