/**
  @class ResultFormatter

  A collection of static methods for outputting RapidFit data objects

  @author Benjamin M Wynne bwynne@cern.ch
  @author Greig A Cowan greig.cowan@cern.ch
  @date 2009-10-02
 */

#include "TFile.h"
#include "TNtuple.h"
#include "Rtypes.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TAxis.h"

#include "ResultFormatter.h"
#include "StatisticsFunctions.h"
#include "EdStyle.h"
#include "StringProcessing.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h>

//Output data as a RootNTuple
void ResultFormatter::MakeRootDataFile( string FileName, IDataSet * OutputData )
{
	cout << "ResultFormatter writing to " << FileName << endl;
	//Make a string naming all observables
	string observableNames = "";
	vector<string> allNames = OutputData->GetBoundary()->GetAllNames();
	for (unsigned short int nameIndex = 0; nameIndex < allNames.size(); ++nameIndex)
	{
		if (nameIndex == 0)
		{
			observableNames += allNames[0];
		}
		else
		{
			observableNames += ":" + allNames[nameIndex];
		}
	}

	//Make the file and NTuple
	TFile * rootFile = new TFile( FileName.c_str(), "RECREATE" );
	TNtuple * dataNTuple = new TNtuple( "dataNTuple", "All data", observableNames.c_str() );

	//Loop over all data points and add them to the NTuple
	for ( int dataIndex = 0; dataIndex < OutputData->GetDataNumber(); ++dataIndex)
	{
		//Retrieve the values of all observables
		Float_t observables[ allNames.size() ];
		for (unsigned short int nameIndex = 0; nameIndex < allNames.size(); ++nameIndex)
		{
			DataPoint * temporaryDataPoint = OutputData->GetDataPoint(dataIndex);
			observables[nameIndex] = Float_t(temporaryDataPoint->GetObservable( allNames[nameIndex] )->GetValue());
			//delete temporaryDataPoint;
		}

		//Populate the NTuple
		dataNTuple->Fill(observables);
	}

	//Write the file
	rootFile->Write("dataNTuple");
	rootFile->Close();
}

//Display the results of a fit using cout
void ResultFormatter::DebugOutputFitResult( FitResult * OutputData )
{
	cout << "Fit status: " << OutputData->GetFitStatus() << endl;
	cout << "Minimum function value: " << OutputData->GetMinimumValue() << endl;
	cout << "Name | Value | Minimum | Maximum" << endl;

	//Ouput each parameter
	ResultParameterSet * outputParameters = OutputData->GetResultParameterSet();
	vector<string> allNames = outputParameters->GetAllNames();
	vector<string>::iterator nameIterator;
	for ( nameIterator = allNames.begin(); nameIterator != allNames.end(); ++nameIterator )
	{
		ResultParameter * outputParameter = outputParameters->GetResultParameter( *nameIterator );
		cout << *nameIterator << " | " << outputParameter->GetValue() << " | ";
		cout << outputParameter->GetMinimum() << " | " << outputParameter->GetMaximum() << endl;
	}
}

void ResultFormatter::PlotFitContours( FitResult * OutputData, string contourFileName )
{
	vector< FunctionContour* > contours = OutputData->GetContours();
	if ( contours.size() == 0 )
	{
		return;
	}

	string name = contourFileName;// + ".root";
	TFile * contourFile = new TFile( name.c_str(), "RECREATE");
	int colours[2] = {42,38};

	//Loop over all contour plots
	for (unsigned short int plotIndex = 0; plotIndex < contours.size(); ++plotIndex )
	{
		TMultiGraph * graph = new TMultiGraph();
		FunctionContour * plotContour = contours[plotIndex];

		//Make the plot canvas
		string canvasName = plotContour->GetXName() + "vs" + plotContour->GetYName() + "Contour";
		string canvasTitle = plotContour->GetXName() + " vs " + plotContour->GetYName() + " Contour";
		TCanvas * bothPlots = new TCanvas( canvasName.c_str(), canvasTitle.c_str() );

		//Plot each contour, starting at highest sigma
		for ( int sigma = plotContour->GetContourNumber(); sigma > 0; --sigma )
		{
			vector< pair< double, double > > sigmaContour = plotContour->GetPlot(sigma);

			//Retrieve each point
			double xCoordinates[ sigmaContour.size() ], yCoordinates[ sigmaContour.size() ];
			for (unsigned short int pointIndex = 0; pointIndex < sigmaContour.size(); ++pointIndex )
			{
				xCoordinates[pointIndex] = sigmaContour[pointIndex].first;
				yCoordinates[pointIndex] = sigmaContour[pointIndex].second;
			}

			//Make the graph
			TGraphErrors * contourGraph = new TGraphErrors( int(sigmaContour.size()), xCoordinates, yCoordinates );
			contourGraph->SetFillColor( Color_t(colours[sigma-1]) );
			graph->Add(contourGraph);
		}

		//Format the graph
		graph->SetTitle("1 and 2 sigma contours");
		graph->Draw( "ALF" ); //Smooth fill area drawn

		//Titles in format: ParameterName (ParameterUnit)
		string xTitle = plotContour->GetXName() + " (" + OutputData->GetResultParameterSet()->GetResultParameter( plotContour->GetXName() )->GetUnit() + ")";
		string yTitle = plotContour->GetYName() + " (" + OutputData->GetResultParameterSet()->GetResultParameter( plotContour->GetYName() )->GetUnit() + ")";
		graph->GetXaxis()->SetTitle( xTitle.c_str() );
		graph->GetYaxis()->SetTitle( yTitle.c_str() );

		//Store the graph
		bothPlots->Modified();
		bothPlots->Update();
		bothPlots->Write();
		//bothPlots->SaveAs( (contourFileName + "." + plotContour->GetXName() + "." + plotContour->GetYName() + ".png").c_str() );
	}

	contourFile->Close();
}

//Display the covariance matrix of a fit in a LaTeX table using cout
void ResultFormatter::LatexOutputCovarianceMatrix( FitResult * OutputData )
{
	vector<double> covarianceMatrix = OutputData->GetCovarianceMatrix();
	vector<string> allNames = OutputData->GetResultParameterSet()->GetAllNames();
	vector<string>::iterator nameIterator;
	int numberOfFreeParameters = 0;

	string columns = "\\begin{tabular}{|c|";
	string parameterNames = "";
	for ( nameIterator = allNames.begin(); nameIterator != allNames.end(); ++nameIterator )
	{
		bool isFree = ResultFormatter::IsParameterFree( OutputData, *nameIterator );
		if (isFree)
		{
			columns += "c|";
			//string name = FindAndReplaceString( *nameIterator );
			//string name = StringProcessing::ReplaceString( *nameIterator, "_", "\\_" );
			string name = *nameIterator;
			std::stringstream ResultStream;
			ResultStream << setw(10) << EdStyle::GetParamLatexName( name );
			parameterNames += " & " + ResultStream.str();
			numberOfFreeParameters += 1;
		}
	}
	columns += "}\n\\hline";
	parameterNames += "\\\\ \\hline \\hline";

	cout << "Correlation matrix" << endl;
	cout << "\n\\begin{center}" << endl;
	cout << columns << endl;
	cout << setw(20) << " " <<  setw(16) << parameterNames << endl;

	int row = 0;
	for ( nameIterator = allNames.begin(); nameIterator != allNames.end(); ++nameIterator )
	{
		bool isFree = ResultFormatter::IsParameterFree( OutputData, *nameIterator );
		if (!isFree) continue;

		//string name = FindAndReplaceString( *nameIterator );
		//string name = StringProcessing::ReplaceString( *nameIterator, "_", "\\_" );
		string name = *nameIterator;
		cout << setw(20) << EdStyle::GetParamLatexName( name );
		if ( covarianceMatrix.size() == 0 )
		{
			cerr << "No correlation matrix returned from fit!" << endl;
			break;
		}

		double drow = GetElementFromCovarianceMatrix( covarianceMatrix, row, row );

		for ( int col = 0; col < numberOfFreeParameters; ++col)
		{
			double dcol = GetElementFromCovarianceMatrix( covarianceMatrix, col, col );
			double covariance = GetElementFromCovarianceMatrix( covarianceMatrix, row, col );
			double correlation = covariance/sqrt(fabs(drow * dcol));
			if (col >= row)
			{
				if ( fabs(correlation) > 0.5 && ( col != row ) )
				{
					std::stringstream ResultStream;
					ResultStream << std::setprecision(2) << correlation;
					TString formatted("\\bf{"); formatted.Append( ResultStream.str() ); formatted.Append("}") ;
					cout << " & " << setw(12) << std::setprecision(2) << formatted;
				}
				else
				{
					cout << " & " << setw(12) << std::setprecision(2) << correlation;
				}
			}
			else
			{
				cout << " & " << setw(12) << " ";
			}

		}
		cout << " \\\\" << endl;
		row += 1;
	}

	cout << "\\hline \n\\end{tabular}" << endl;
	cout << "\\end{center}\n" << endl;
}

bool ResultFormatter::IsParameterFree( FitResult * OutputData, string ParameterName )
{
	bool decision = true;
	string type = OutputData->GetResultParameterSet()->GetResultParameter( ParameterName )->GetType();
	if ( type == "Fixed") decision = false;
	return decision;
}

double ResultFormatter::GetElementFromCovarianceMatrix( vector<double> matrix, int row, int col)
{
	if(row > col) return matrix[col+row*(row+1)/2];
	else return matrix[row+col*(col+1)/2];
}


//Display the results of a fit in a LaTeX table using cout
void ResultFormatter::LatexOutputFitResult( FitResult * OutputData )
{
	//.............................................
	// Standard table for MC toys with pulls
	cout << "Fit result for MC toys with pulls" << endl;
	cout << "\n\\begin{center}" << endl;
	cout << "Fit status: " << OutputData->GetFitStatus() << endl;
	cout << setprecision(8) << "Minimum function value: " << OutputData->GetMinimumValue() << endl;
	cout << "\\begin{tabular}{|c|c|c|} \n\\hline" << endl;
	cout << setw(20) << "Parameter"<< " & " << setw(25) << "Fit result and error" << setw(21) << " & " << setw(20) << "$\\sigma$ from input \\\\ \t\t\\hline \\hline\n" << endl;

	//Ouput each parameter
	ResultParameterSet * outputParameters = OutputData->GetResultParameterSet();
	vector<string> allNames = outputParameters->GetAllNames();
	vector<string>::iterator nameIterator;
	for ( nameIterator = allNames.begin(); nameIterator != allNames.end(); ++nameIterator )
	{
		ResultParameter * outputParameter = outputParameters->GetResultParameter( *nameIterator );

		double fitValue = outputParameter->GetValue();
//		double minValue = outputParameter->GetMinimum();
//		double inputValue = outputParameter->GetOriginalValue();
		double fitError = outputParameter->GetError();
		double sigmaFromInputValue = outputParameter->GetPull();
		string unit = outputParameter->GetUnit(); 
		//if (fitError > 0.0) sigmaFromInputValue = (fitValue - inputValue)/fitError;

		//boost::regex pattern ("_",boost::regex_constants::icase|boost::regex_constants::perl);
		//string replace ("\\_");
		//string newName = boost::regex_replace (*nameIterator, pattern, replace);

		//string name = FindAndReplaceString( *nameIterator );
		//string name = StringProcessing::ReplaceString( *nameIterator, "_", "\\_" );
		string name = *nameIterator;
		cout << setw(20) << EdStyle::GetParamLatexName(name) << " & "
			<< setw(12) << setprecision(5) << fitValue << " $\\pm$ "
			<< setw(10) <<  		  fitError << " " << setw(15) << EdStyle::GetParamLatexUnit(unit) << " & "
			<< setw(20) << setprecision(2) << sigmaFromInputValue << "\\\\" << endl;
	}

	cout << "\\hline \n\\end{tabular}" << endl;
	cout << "\\end{center}\n" << endl;

	//.................................................
	//longer table for MC pull fits with absolute offsets
	cout << endl ;
	cout << "Fit result - for MC toys with pulls and absolute offsets " << endl;
	cout << "\n\\begin{center}" << endl;
	cout << "Fit status: " << OutputData->GetFitStatus() << endl;
	cout << setprecision(8) << "Minimum function value: " << OutputData->GetMinimumValue() << endl;
	cout << "\\begin{tabular}{|c|c|c|c|} \n\\hline" << endl;
	cout << setw(20)<< "Parameter"<< " & " << setw(25) << "Fit result and error" << setw(21) << " & "<< setw(20) <<"$\\sigma$ from input" << " & " << setw(20) << "Abs from input \\\\ \t\t\\hline \\hline\n" << endl;

	//Ouput each parameter
	for ( nameIterator = allNames.begin(); nameIterator != allNames.end(); ++nameIterator )
	{
		ResultParameter * outputParameter = outputParameters->GetResultParameter( *nameIterator );

		double fitValue = outputParameter->GetValue();
//		double minValue = outputParameter->GetMinimum();
		double inputValue = outputParameter->GetOriginalValue();
		double fitError = outputParameter->GetError();
		double sigmaFromInputValue = outputParameter->GetPull();
		string unit = outputParameter->GetUnit();
		//if (fitError > 0.0) sigmaFromInputValue = (fitValue - inputValue)/fitError;

		//boost::regex pattern ("_",boost::regex_constants::icase|boost::regex_constants::perl);
		//string replace ("\\_");
		//string newName = boost::regex_replace (*nameIterator, pattern, replace);

		//string name = FindAndReplaceString( *nameIterator );
		//string name = StringProcessing::ReplaceString( *nameIterator, "_", "\\_" );
		string name = *nameIterator;
		cout << setw(20) << EdStyle::GetParamLatexName(name) << " & "
			<< setw(12) << setprecision(5) << fitValue << " $\\pm$ "
			<< setw(10) <<  		  fitError << " " << setw(15) << EdStyle::GetParamLatexUnit(unit) << " & "
			<< setw(20) << setprecision(2) << sigmaFromInputValue << " & "
			<< setw(15) << setprecision(5) << fitValue-inputValue << "\\\\" << endl;
	}

	cout << "\\hline \n\\end{tabular}" << endl;
	cout << "\\end{center}\n" << endl;

	//........................................
	//short table for data fits
	cout << endl ;
	cout << "\n\\begin{center}" << endl;
	cout << "Fit result - for Data fits" << endl;
	cout << "Fit status: " << OutputData->GetFitStatus() << endl;
	cout << setprecision(8) << "Minimum function value: " << OutputData->GetMinimumValue() << endl;
	cout << "\\begin{tabular}{|c|c|} \n\\hline" << endl;
	cout << setw(20) << "Parameter" << " & " << setw(21) << "Fit result and error" << setw(21) << " " << " \\\\ \\hline \\hline\n" << endl;

	//Will need to do some comparisons
//	double Rperp =0, Rzp =0, ePerp =0 , eZp=0;

	//Ouput each parameter
	for ( nameIterator = allNames.begin(); nameIterator != allNames.end(); ++nameIterator )
	{
		ResultParameter * outputParameter = outputParameters->GetResultParameter( *nameIterator );

		double fitValue = outputParameter->GetValue();
//		double minValue = outputParameter->GetMinimum();
//		double inputValue = outputParameter->GetOriginalValue();
		double fitError = outputParameter->GetError();
//		double sigmaFromInputValue = outputParameter->GetPull();
		string unit = outputParameter->GetUnit();
		//string name = StringProcessing::ReplaceString( *nameIterator, "_", "\\_" );
		string name = *nameIterator;
		cout << setw(20) << EdStyle::GetParamLatexName(name) << " & "
			<< setw(12) << setprecision(3) << fitValue << " $\\pm$ "
			<< setw(10) <<  		  fitError << " " << setw(15) << EdStyle::GetParamLatexUnit(unit)  << "\\\\" << endl;
	}
	cout << "\\hline \n\\end{tabular}" << endl;
	cout << "\\end{center}\n" << endl;	

	cout << endl << endl;
	cout << "--------------------------------------------------" <<endl;
	cout << "\nFit Review:\t\tStatus:\t" <<OutputData->GetFitStatus()<<endl<<endl;

	//Ouput each parameter
	for ( nameIterator = allNames.begin(); nameIterator != allNames.end(); ++nameIterator )
	{
		ResultParameter * outputParameter = outputParameters->GetResultParameter( *nameIterator );
		double fitValue = outputParameter->GetValue();
		double fitError = outputParameter->GetError();
		string unit = outputParameter->GetUnit();
		string name = *nameIterator;
		cout << setw(25) << name << " : "
			<< setw(13) << setprecision(5) << fitValue << "  \\pm  "
			<< setw(13) << setprecision(5) << fitError << endl;
	}
	cout << endl;
	cout << "--------------------------------------------------" <<endl;
	cout << endl <<endl;

}

/*string ResultFormatter::FindAndReplaceString( string name )
  {
// This isn't very general, and probably won't work for names
// containing lots of "_".
int pos = 0;
int newPos = 0;
int size = name.size();
while ( pos < size )
{
newPos = name.find("_", pos);
if( newPos != string::npos ) {
name.replace(newPos, 1, "\\_");
}
else break;
pos = newPos + 2;
}
return name;
}*/

//Chose which pull plot method to use
void ResultFormatter::MakePullPlots( string Type, string FileName, ToyStudyResult * ToyResult )
{
	if ( Type == "FlatNTuple" )
	{
		return FlatNTuplePullPlots( FileName, ToyResult );
	}
	else if ( Type == "SeparateParameter" )
	{
		return SeparateParameterPullPlots( FileName, ToyResult );
	}
	else
	{
		cout << "Unrecognised pull plot type \"" << Type << "\" - defaulting to SeparateParameter" << endl;
		return SeparateParameterPullPlots( FileName, ToyResult );
	}
}

void ResultFormatter::WriteFlatNtuple( string Filename, ToyStudyResult* ToyResult )
{
	ResultFormatter::FlatNTuplePullPlots( Filename, ToyResult );
}

//Make pull plots from the output of a toy study
void ResultFormatter::FlatNTuplePullPlots( string FileName, ToyStudyResult * ToyResult )
{
	TFile * rootFile = new TFile( FileName.c_str(), "RECREATE" );
	TNtuple * parameterNTuple;
	parameterNTuple = new TNtuple("RapidFitResult", "RapidFitResult", ToyResult->GetFlatResultHeader());
	Float_t * resultArr;
	for ( int resultIndex = 0; resultIndex < ToyResult->NumberResults(); ++resultIndex )
	{
		vector<double> result = ToyResult->GetFlatResult(resultIndex);
		resultArr = new Float_t [result.size()];
		copy( result.begin(), result.end(), resultArr);
		parameterNTuple->Fill(resultArr);
		delete [] resultArr;
	}
	rootFile->Write();
	rootFile->Close();
}

//Make pull plots from the output of a toy study
void ResultFormatter::SeparateParameterPullPlots( string FileName, ToyStudyResult * ToyResult )
{
	TFile * rootFile = new TFile( FileName.c_str(), "RECREATE" );
	string header = "value:error:pull";
	vector<string> allNames = ToyResult->GetAllNames();
	Float_t valueErrorPull[3];
	vector<double> parameterValues, parameterErrors, parameterPulls;
	TH1F * pullHistogram;
	TNtuple * parameterNTuple;

	//Plots for each observable
	for (unsigned short int nameIndex = 0; nameIndex < allNames.size(); ++nameIndex)
	{
		//Prepare the NTuple
		parameterNTuple = new TNtuple( allNames[nameIndex].c_str(), "Parameter fit results", header.c_str() );
		parameterValues = ToyResult->GetParameterValues( allNames[nameIndex] );
		parameterErrors = ToyResult->GetParameterErrors( allNames[nameIndex] );
		parameterPulls = ToyResult->GetParameterPulls( allNames[nameIndex] );

		//Prepare the pull histogram
		string histogramName = allNames[nameIndex] + "PullPlot";
		string histogramTitle = allNames[nameIndex] + " pull plot";
		double maximumPull = StatisticsFunctions::Maximum(parameterPulls);
		double minimumPull = StatisticsFunctions::Minimum(parameterPulls);
		bool makeHistogram = !isnan(maximumPull) && !isnan(minimumPull);
		if (makeHistogram)
		{
			pullHistogram = new TH1F( histogramName.c_str(), histogramTitle.c_str(),
					StatisticsFunctions::OptimumBinNumber(parameterPulls) + 2, minimumPull, maximumPull );
		}

		//Plot the results
		for ( int resultIndex = 0; resultIndex < ToyResult->NumberResults(); ++resultIndex )
		{
			valueErrorPull[0] = Float_t(parameterValues[resultIndex]);
			valueErrorPull[1] = Float_t(parameterErrors[resultIndex]);
			valueErrorPull[2] = Float_t(parameterPulls[resultIndex]);

			parameterNTuple->Fill(valueErrorPull);

			if (makeHistogram)
			{
				pullHistogram->Fill( parameterPulls[resultIndex] );
			}
		}

		//Fit a gaussian distribution to the distribution
		if (makeHistogram)
		{
			pullHistogram->Fit("gaus");
			EdStyle * greigFormat = new EdStyle();
			pullHistogram->UseCurrentStyle();
			delete greigFormat;
		}
	}

	//Write out the fit times as well
	vector<double> allRealTimes = ToyResult->GetAllRealTimes();
	vector<double> allCPUTimes = ToyResult->GetAllCPUTimes();
	TNtuple * fitInfoNTuple = new TNtuple( "fitInfo", "Information about fits", "realTime:cpuTime:fitStatus" );
	Float_t timeCPUStatus[3];
	for (unsigned short  int timeIndex = 0; timeIndex < allRealTimes.size(); ++timeIndex )
	{
		timeCPUStatus[0] = Float_t(allRealTimes[timeIndex]);
		timeCPUStatus[1] = Float_t(allCPUTimes[timeIndex]);
		timeCPUStatus[2] = Float_t(ToyResult->GetFitResult(timeIndex)->GetFitStatus());

		fitInfoNTuple->Fill( timeCPUStatus );
	}

	//Write the file
	rootFile->Write();
	rootFile->Close();
}


//.........................................
// New Method form PELC to plot LL scan results

void ResultFormatter::MakeLLscanPlots( vector<LLscanResult*> scanResults, string filename )
{
	// 1=BLACK    4=BLUE    3=GREEN     5=YELLOW

	TFile * LLscanFile = new TFile( filename.c_str(), "RECREATE");
	LLscanFile->SetCompressionLevel(9);


	//Set some numbers
	int nscans = int(scanResults.size());

	//Make a canvas for all the plots
	TCanvas cv ;

	for( int ii=0; ii<nscans; ++ii )
	{
		cv.SetGrid();
		cv.GetFrame()->SetFillColor(21);
		cv.GetFrame()->SetBorderSize(12);

		cv.cd(ii+1) ;
		TGraph * grnew = scanResults[ii]->GetGraph() ;
		grnew->Draw("ALP") ;
		grnew->Draw() ;
		//grnew->Write();
		cv.Update();
		cv.Write();

	}

	LLscanFile->Close();

}

//.........................................
// Send LL contour results to a file

void ResultFormatter::MakeLLcontourPlots( vector<LLscanResult2D*> scanResults, string filename )
{	
	// 1=BLACK    4=BLUE    3=GREEN     5=YELLOW

	TFile * LLcontourFile = new TFile( filename.c_str(), "RECREATE");
	LLcontourFile->SetCompressionLevel(9);

	//Set some numbers
	int nscans = int(scanResults.size());

	for( int ii=0; ii<nscans; ++ii )
	{
		TH2D * hist = scanResults[ii]->GetTH2D() ;
		hist->Draw("cont1") ;
		hist->Write();		
	}

	LLcontourFile->Close();

}


//====================================================================================================
//Do a likelihood scan
LLscanResult* ResultFormatter::LLScan( ToyStudyResult* new_results, string scanName )
{

	cout << "Constructing LLscan for parameter " << scanName << endl ;

	// Need to set up a loop , fixing the scan parameter at each point
	vector<double> scanParameterValues = new_results->GetParameterValues( scanName );
	vector<double> scanLLValues = new_results->GetAllMLL() ;

	LLscanResult * result = new LLscanResult( scanName, scanParameterValues, scanLLValues ) ;
	result->print() ; //PELC

	return result;
}

LLscanResult2D* ResultFormatter::LLScan2D( vector<ToyStudyResult*> new_results, string scanName, string scanName2 )
{

	// THIS CODE WILL NOT WORK WITH ANYTHING LESS THAN A PERFECT SQUARE OF FIT RESULTS

	cout << "Constructing LLcontour for parameters " << scanName << "\t" << scanName2 << endl ;

	vector<LLscanResult* > LLScanResults;
	vector<double> scanParameterValues;
	vector<double> scanParameterValues2 = new_results[0]->GetParameterValues( scanName2 );

	for(unsigned int si=0; si < new_results.size(); ++si) {
		vector<double> scanLLValues = new_results[si]->GetAllMLL();
		LLscanResult * _1D_temp_result = new LLscanResult( scanName2, scanParameterValues2, scanLLValues ) ;
		LLScanResults.push_back( _1D_temp_result );
		vector<double> temp_values = new_results[si]->GetParameterValues( scanName );
		scanParameterValues.push_back( temp_values[0] );
	}

	LLscanResult2D * result = new LLscanResult2D( scanName, scanParameterValues, scanName2, scanParameterValues2, LLScanResults );

	return result;
}

