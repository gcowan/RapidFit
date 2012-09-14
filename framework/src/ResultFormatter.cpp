/**
  @class ResultFormatter

  A collection of static methods for outputting RapidFit data objects

  @author Benjamin M Wynne bwynne@cern.ch
  @author Greig A Cowan greig.cowan@cern.ch
  @date 2009-10-02
  */

///	ROOT Headers
#include "TFile.h"
#include "TNtuple.h"
#include "Rtypes.h"
//	Used for Minimiser Contours
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TSystem.h"
///	RapidFit Headers
#include "ResultFormatter.h"
#include "StatisticsFunctions.h"
#include "EdStyle.h"
#include "StringProcessing.h"
#include "RapidFitMatrix.h"
///	System Headers
#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <cmath>
#include <sstream>
#include <fstream>

using namespace::std;

//Output data as a RootNTuple
void ResultFormatter::MakeRootDataFile( string FullFileName, vector<IDataSet*> OutputData )
{

	//	Remove the extention from the filename if it exists
	string ext_dot=".";
	vector<string> temp_strings = StringProcessing::SplitString( FullFileName, *(ext_dot.c_str()) );
	TString FileName_Pre_Suffix = StringProcessing::CondenseStrings( temp_strings, 0, int(temp_strings.size() -1) );

	int counter = 0;
	for( vector<IDataSet*>::iterator data_iter = OutputData.begin(); data_iter!=OutputData.end(); ++data_iter, ++counter )
	{
		TString TString_FileName = FileName_Pre_Suffix;
		//	If we have multiple datasets write it to multuple files
		if( OutputData.size() > 1 )
		{
			TString_FileName.Append("_");
			TString_FileName+=counter;
		}
		TString_FileName.Append(".root");

		string FileName = TString_FileName.Data();

		cout << "ResultFormatter writing to " << FileName << endl;
		//Make a string naming all observables
		string observableNames = "";
		vector<string> allNames = (*data_iter)->GetBoundary()->GetAllNames();
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
		for ( int dataIndex = 0; dataIndex < (*data_iter)->GetDataNumber(); ++dataIndex)
		{
			//Retrieve the values of all observables
			Float_t* observables = new Float_t[ allNames.size() ];
			for (unsigned short int nameIndex = 0; nameIndex < allNames.size(); ++nameIndex)
			{
				DataPoint * temporaryDataPoint = (*data_iter)->GetDataPoint(dataIndex);
				observables[nameIndex] = Float_t(temporaryDataPoint->GetObservable( allNames[nameIndex] )->GetValue());
				//delete temporaryDataPoint;
			}

			//Populate the NTuple
			dataNTuple->Fill(observables);
			delete[] observables;
		}

		//Write the file
		rootFile->Write("dataNTuple");
		rootFile->Close();

	}
	return;
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
	int colours[3] = {2,3,4};
	TString confs[3] = {"68.0","90.0","95.0"};


	//Loop over all contour plots
	for (unsigned short int plotIndex = 0; plotIndex < contours.size(); ++plotIndex )
	{


		TMultiGraph * graph = new TMultiGraph();
		FunctionContour * plotContour = contours[plotIndex];

		//Make the plot canvas
		string canvasName = plotContour->GetXName() + "vs" + plotContour->GetYName() + "Profile LL";
		string canvasTitle = plotContour->GetXName() + " vs " + plotContour->GetYName() + " Profile LL";
		TCanvas * bothPlots = new TCanvas( canvasName.c_str(), canvasTitle.c_str() );
		TLegend *leg = new TLegend(0.80,0.89,0.95,0.7);
		leg->SetHeader("Conf. Levels");
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);

		//Plot each contour, starting at highest sigma
		for ( int sigma = plotContour->GetContourNumber(); sigma > 0; --sigma )
		{
			TString confname = "";
			confname +=confs[sigma-1];
			confname += "\% C.L.";

			vector< pair< double, double > > sigmaContour = plotContour->GetPlot(sigma);

			//Retrieve each point
			double* xCoordinates = new double[ sigmaContour.size() ];
			double* yCoordinates = new double[ sigmaContour.size() ];
			for (unsigned short int pointIndex = 0; pointIndex < sigmaContour.size(); ++pointIndex )
			{
				xCoordinates[pointIndex] = sigmaContour[pointIndex].first;
				yCoordinates[pointIndex] = sigmaContour[pointIndex].second;
			}

			//Make the graph
			TGraphErrors * contourGraph = new TGraphErrors( int(sigmaContour.size()), xCoordinates, yCoordinates );

			contourGraph->SetLineColor(Color_t(colours[sigma-1]));
			contourGraph->SetLineWidth(2);
			leg->AddEntry(contourGraph,confname, "L");
			graph->Add(contourGraph);
		}

		//Format the graph
		graph->SetTitle("Profile LL Contours");
		graph->Draw( "AL" ); //Smooth fill area drawn //FIXME
		leg->Draw();


		//Titles in format: ParameterName (ParameterUnit)
		string xTitle = plotContour->GetXName() + " (" + OutputData->GetResultParameterSet()->GetResultParameter( plotContour->GetXName() )->GetUnit() + ")";
		string yTitle = plotContour->GetYName() + " (" + OutputData->GetResultParameterSet()->GetResultParameter( plotContour->GetYName() )->GetUnit() + ")";
		graph->GetXaxis()->SetTitle( xTitle.c_str() );
		graph->GetYaxis()->SetTitle( yTitle.c_str() );

		//Store the graph
		bothPlots->Modified();
		bothPlots->Update();
		bothPlots->Write();
		bothPlots->SaveAs( (contourFileName + "_" + plotContour->GetXName() + "_" + plotContour->GetYName() + ".pdf").c_str() );
	}

	contourFile->Close();
}

//Display the covariance matrix of a fit in a LaTeX table using cout
void ResultFormatter::LatexOutputCovarianceMatrix( FitResult * OutputData )
{
	TMatrixDSym* covarianceMatrix = OutputData->GetCovarianceMatrix()->thisMatrix;
	vector<string> freeNames = OutputData->GetCovarianceMatrix()->theseParameters;

	vector<string>::iterator nameIterator;

	int numberOfFreeParameters = 0;

	string columns = "\\begin{tabular}{|c|";
	string parameterNames = "";
	for( nameIterator = freeNames.begin(); nameIterator != freeNames.end(); ++nameIterator )
	{
		columns += "c|";
		//string name = FindAndReplaceString( *nameIterator );
		//string name = StringProcessing::ReplaceString( *nameIterator, "_", "\\_" );
		std::stringstream ResultStream;
		ResultStream << setw(10) << EdStyle::GetParamLatexName( *nameIterator );
		parameterNames += " & " + ResultStream.str();
		numberOfFreeParameters += 1;
	}
	columns += "}\n\\hline";
	parameterNames += "\\\\ \\hline \\hline";

	if( covarianceMatrix == NULL )
	{
		cerr << "No correlation matrix returned from fit!" << endl;
	}
	else
	{
		cout << "Correlation matrix" << endl;
		cout << "\n\\begin{center}" << endl;
		cout << columns << endl;
		cout << setw(20) << " " <<  setw(16) << parameterNames << endl;

		int row = 0;
		for ( nameIterator = freeNames.begin(); nameIterator != freeNames.end(); ++nameIterator )
		{

			string name = *nameIterator;
			cout << setw(20) << EdStyle::GetParamLatexName( name );

			double drow = (*covarianceMatrix)( row, row );
			for ( int col = 0; col < numberOfFreeParameters; ++col)
			{
				double dcol = (*covarianceMatrix)( col, col );

				double covariance = (*covarianceMatrix)( row, col );
				double correlation = covariance/sqrt(fabs(drow * dcol));
				if( col >= row )
				{
					//	If there is a strong Correlation Make the Number Bold
					if( fabs(correlation) > 0.5 && ( col != row ) )
					{
						std::stringstream ResultStream;
						ResultStream << std::setprecision(2) << correlation;
						TString formatted("\\bf{"); formatted.Append( ResultStream.str() ); formatted.Append("}") ;
						cout << " & " << std::setprecision(2) << fixed << formatted;
					}
					else // else just print the number
					{
						cout << " & " << setw(12) << std::setprecision(2) << fixed << correlation;
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

	cout.unsetf( ios_base::fixed );
}

bool ResultFormatter::IsParameterFree( FitResult * OutputData, string ParameterName )
{
	bool decision = true;
	string type = OutputData->GetResultParameterSet()->GetResultParameter( ParameterName )->GetType();
	if( type == "Fixed") decision = false;
	return decision;
}

void ResultFormatter::WriteOutputLatex( FitResult* OutputData )
{
	TString output_folder("OutputLatex_");
	output_folder.Append( StringProcessing::TimeString() );

	gSystem->mkdir( output_folder );

	stringstream latex;

	ResultFormatter::LatexDocHeader( latex );

	latex << "\n\\input{ ./FullTable }" << endl;

	latex << "\n\\clearpage" << endl;

	latex << "\n\\input{ ./SimpleTable }" << endl;

	latex << "\n\\clearpage" << endl;

	latex << "\n\\input{ ./MinimalTable }" << endl;

	latex << "\n\\clearpage" << endl;

	latex << "\n\\input{ covMatrix }" << endl;

	latex << "\n\\clearpage\n" << endl;

	ResultFormatter::LatexDocFooter( latex );

	ofstream outFile;
	outFile.open( output_folder+"/main.tex" );
	outFile << latex.str();
	outFile.close();

	outFile.open( output_folder+"/FullTable.tex" );
	stringstream full;
	ResultFormatter::LatexFullFitResultTable( OutputData, full );
	outFile << full.str();
	outFile.close();

	outFile.open( output_folder+"/SimpleTable.tex" );
	stringstream simple;
	ResultFormatter::LatexSimpleFitResultTable( OutputData, simple );
	outFile << simple.str();
	outFile.close();

	outFile.open( output_folder+"/MinimalTable.tex" );
	stringstream minimal;
	ResultFormatter::LatexMinimumFitResultTable( OutputData, minimal );
	outFile << minimal.str();
	outFile.close();

	outFile.open( output_folder+"/covMatrix.tex" );
	stringstream matrix;
	ResultFormatter::LatexCovMatrix( OutputData, matrix );
	outFile << matrix.str();
	outFile.close();

	gSystem->cd( output_folder );

	system("pdflatex main.tex >/dev/null");

	gSystem->cd( ".." );
}

void ResultFormatter::LatexDocHeader( stringstream& latex )
{
	latex << "\\documentclass[a4paper,10pt]{article}" << endl;
	latex << "\\nonstopmode" << endl;
	latex << "\\usepackage[utf8]{inputenc}" << endl;
	latex << "\\title{}" << endl;
	latex << "\\author{}" << endl;
	latex << "\\date{}" << endl;
	latex << endl;
	latex << "\\pdfinfo{" << endl;
	latex << "/Title    ()" << endl;
	latex << "/Author   ()" << endl;
	latex << "/Creator  ()" << endl;
	latex << "/Producer ()" << endl;
	latex << "/Subject  ()" << endl;
	latex << "/Keywords ()" << endl;
	latex << "}" << endl;
	latex << endl;
	latex << "\\usepackage{amsmath}" << endl;
	latex << "\\usepackage{caption}" << endl;
	latex << "\\usepackage{rotating}" << endl;
	//latex << "\\usepackage{subcaption}" << endl;
	latex << "\\usepackage{graphicx}" << endl;
	latex << "\\usepackage[a4paper]{geometry}" << endl;
	latex << endl;
	latex << "\\let\\oldpm\\pm" << endl;
	latex << endl;
	latex << "\\renewcommand*{\\arraystretch}{1.}" << endl;
	latex << endl;
	latex << "\\begin{document}" << endl;
	latex << endl;
}

void ResultFormatter::LatexDocFooter( stringstream& latex )
{
	latex << endl;
	latex << "\\end{document}" << endl;
	latex << endl;
}

void ResultFormatter::TableHeader( stringstream& latex, unsigned int colnum )
{
	latex << endl;
	latex << "\\renewcommand{\\pm}{\\ensuremath{\\oldpm} }" << endl;
	latex << "\\begin{table}[h]" << endl;
	latex << "\\begin{center}" << endl;
	latex << "\\begin{tabular}{@{}|l";
	for( unsigned int i=1; i< colnum; ++i )
	{
		latex << "|r";
	}
	latex << "|@{}}" << endl;
	latex << "\\hline" << endl;
}

void ResultFormatter::TableHeaderLandscape( stringstream& latex, unsigned int colnum )
{
	latex << "\\newgeometry{left=3cm,bottom=2cm}" << endl;
	latex << "\\renewcommand{\\pm}{\\ensuremath{\\oldpm} }" << endl;
	latex << "\\begin{sidewaystable}[h]" << endl;
	latex << "\\begin{center}" << endl;
	latex << "\\begin{tabular}{@{}|l";
	for( unsigned int i=1; i< colnum; ++i )
	{
		latex << "|r";
	}
	latex << "|@{}}" << endl;
	latex << "\\hline" << endl;
}

void ResultFormatter::TableFooter( stringstream& latex )
{
	latex << "\\hline" << endl;
	latex << "\\end{tabular}" << endl;
	latex << "\\caption{Some Caption}" << endl;
	latex << "\\label{thisTable}" << endl;
	latex << "\\end{center}" << endl;
	latex << "\\end{table}" << endl;
	latex << "\\renewcommand{\\pm}{\\oldpm}" << endl;
	latex << endl;
}

void ResultFormatter::TableFooterLandscape( stringstream& latex )
{
	latex << "\\hline" << endl;
	latex << "\\end{tabular}" << endl;
	latex << "\\caption{Some Caption}" << endl;
	latex << "\\label{thisTable}" << endl;
	latex << "\\end{center}" << endl;
	latex << "\\end{sidewaystable}" << endl;
	latex << "\\renewcommand{\\pm}{\\oldpm}" << endl;
	latex << "\\restoregeometry" << endl;
	latex << endl;
}

void ResultFormatter::LatexSimpleFitResultTable( FitResult * OutputData, stringstream& latex )
{
	ResultFormatter::TableHeader( latex, 3 );
	latex << "Parameter & Fit result and error & $\\sigma$ from input \\\\ \t\t\\hline \\hline\n" << endl;

	ResultParameterSet * outputParameters = OutputData->GetResultParameterSet();
	vector<string> allNames = outputParameters->GetAllNames();
	for( vector<string>::iterator nameIterator = allNames.begin(); nameIterator != allNames.end(); ++nameIterator )
	{
		ResultParameter * outputParameter = outputParameters->GetResultParameter( *nameIterator );
		double fitValue = outputParameter->GetValue();
		double fitError = outputParameter->GetError();
		double sigmaFromInputValue = outputParameter->GetPull();
		string unit = outputParameter->GetUnit();
		string name = *nameIterator;

		latex << setw(20) << EdStyle::GetParamLatexName(name) << " & "
			<< setw(12) << setprecision(5) << fitValue << " \\pm "
			<< setw(10) <<                    fitError << " " << setw(15) << EdStyle::GetParamLatexUnit(unit) << " & "
			<< setw(20) << setprecision(2) << sigmaFromInputValue << "\\\\" << endl;
	}

	ResultFormatter::TableFooter( latex );
}

void ResultFormatter::LatexFullFitResultTable( FitResult * OutputData, stringstream& latex )
{
	ResultFormatter::TableHeader( latex, 4 );
	latex << "Parameter & Fit result and error & $\\sigma$ from input & Abs from input \\\\ \t\t\\hline \\hline\n" << endl;

	ResultParameterSet * outputParameters = OutputData->GetResultParameterSet();
	vector<string> allNames = outputParameters->GetAllNames();
	for( vector<string>::iterator nameIterator = allNames.begin(); nameIterator != allNames.end(); ++nameIterator )
	{
		ResultParameter * outputParameter = outputParameters->GetResultParameter( *nameIterator );

		double fitValue = outputParameter->GetValue();
		double inputValue = outputParameter->GetOriginalValue();
		double fitError = outputParameter->GetError();
		double sigmaFromInputValue = outputParameter->GetPull();
		string unit = outputParameter->GetUnit();
		string name = *nameIterator;
		latex << setw(20) << EdStyle::GetParamLatexName(name) << " & "
			<< setw(12) << setprecision(5) << fitValue << " \\pm "
			<< setw(10) <<  		  fitError << " " << setw(15) << EdStyle::GetParamLatexUnit(unit) << " & "
			<< setw(20) << setprecision(2) << sigmaFromInputValue << " & "
			<< setw(15) << setprecision(5) << fitValue-inputValue << "\\\\" << endl;
	}

	ResultFormatter::TableFooter( latex );
}

void ResultFormatter::LatexMinimumFitResultTable( FitResult * OutputData, stringstream& latex )
{
	ResultFormatter::TableHeader( latex, 2 );
	latex << "Parameter & Fit result and error  \\\\ \\hline \\hline\n" << endl;

	ResultParameterSet * outputParameters = OutputData->GetResultParameterSet();
	vector<string> allNames = outputParameters->GetAllNames();
	for( vector<string>::iterator nameIterator = allNames.begin(); nameIterator != allNames.end(); ++nameIterator )
	{
		ResultParameter * outputParameter = outputParameters->GetResultParameter( *nameIterator );

		double fitValue = outputParameter->GetValue();
		double fitError = outputParameter->GetError();
		string unit = outputParameter->GetUnit();
		string name = *nameIterator;

		latex << setw(20) << EdStyle::GetParamLatexName(name) << " & "
			<< setw(12) << setprecision(3) << fitValue << " \\pm "
			<< setw(10) <<  		  fitError << " " << setw(15) << EdStyle::GetParamLatexUnit(unit)  << "\\\\" << endl;
	}

	ResultFormatter::TableFooter( latex );
}

//Display the covariance matrix of a fit in a LaTeX table using cout
void ResultFormatter::LatexCovMatrix( FitResult * OutputData, stringstream& latex )
{
	if( OutputData->GetCovarianceMatrix() == NULL ) return;
	TMatrixDSym* covarianceMatrix = OutputData->GetCovarianceMatrix()->thisMatrix;
	if( covarianceMatrix == NULL )
	{
		cerr << "No correlation matrix returned from fit!" << endl;
	}
	else
	{
		vector<string> freeNames = OutputData->GetCovarianceMatrix()->theseParameters;

		if( freeNames.size() > 10 && freeNames.size() < 15 )
		{
			latex << "\\small" << endl;
			latex << "\\renewcommand*{\\arraystretch}{1.}" << endl;
		}
		if( freeNames.size() > 15 && freeNames.size() < 20 )
		{
			latex << "\\footnotesize" << endl;
			latex << "\\renewcommand*{\\arraystretch}{0.95}" << endl;
			latex << "\\renewcommand{\\tabcolsep}{2.pt}" << endl;
		}
		if( freeNames.size() > 20 )
		{
			latex << "\\scriptsize" << endl;
			latex << "\\renewcommand*{\\arraystretch}{0.9}" << endl;
			latex << "\\renewcommand{\\tabcolsep}{1.pt}" << endl;
		}

		ResultFormatter::TableHeaderLandscape( latex, (unsigned)(freeNames.size()+1) );

		string parameterNames = "";
		for( vector<string>::iterator nameIterator = freeNames.begin(); nameIterator != freeNames.end(); ++nameIterator )
		{
			parameterNames += " & " + EdStyle::GetParamLatexName( *nameIterator );
		}

		parameterNames += "\\\\ \\hline \\hline\n";

		latex << parameterNames;

		int row = 0;
		for ( vector<string>::iterator nameIterator = freeNames.begin(); nameIterator != freeNames.end(); ++nameIterator )
		{
			string name = *nameIterator;
			latex << EdStyle::GetParamLatexName( name );
			double drow = (*covarianceMatrix)( row, row );
			for ( int col = 0; col < (int)freeNames.size(); ++col)
			{
				double dcol = (*covarianceMatrix)( col, col );
				double covariance = (*covarianceMatrix)( row, col );
				double correlation = covariance/sqrt(fabs(drow * dcol));
				if( col >= row )
				{
					//	If there is a strong Correlation Make the Number Bold
					if( fabs(correlation) > 0.5 && ( col != row ) )
					{
						std::stringstream ResultStream;
						ResultStream << std::setprecision(2) << correlation;
						TString formatted("\\bf{"); formatted.Append( ResultStream.str() ); formatted.Append("}") ;
						latex << " & " << std::setprecision(2) << fixed << formatted;
					}
					else // else just print the number
					{
						latex << " & " << std::setprecision(2) << fixed << correlation;
					}
				}
				else
				{
					latex << " & ";
				}

			}
			latex << " \\\\" << endl;
			row += 1;
		}

		ResultFormatter::TableFooterLandscape( latex );
	}
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

}


//Display the results of a fit in a LaTeX table using cout
void ResultFormatter::ReviewOutput( FitResult * OutputData )
{
	ResultParameterSet * outputParameters = OutputData->GetResultParameterSet();
	vector<string> allNames = outputParameters->GetAllNames();
	vector<string>::iterator nameIterator;

	cout << endl << endl;
	cout << "--------------------------------------------------" <<endl;
	cout << "\nFit Review:\t\tStatus:\t" <<OutputData->GetFitStatus()<<"\t\tNLL:\t"<<setprecision(10)<<OutputData->GetMinimumValue()<<endl<<endl;

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

//Chose which pull plot method to use
void ResultFormatter::MakePullPlots( string Type, string FileName, FitResultVector* ToyResult )
{
	if ( Type == "FlatNTuple" )
	{
		return WriteFlatNtuple( FileName, ToyResult );
	}
	else if ( Type == "SeparateParameter" )
	{
		return SeparateParameterPullPlots( FileName, ToyResult );
	}
	else
	{
		//cout << "Unrecognised pull plot type \"" << Type << "\" - defaulting to SeparateParameter" << endl;
		//return SeparateParameterPullPlots( FileName, ToyResult );
		return WriteFlatNtuple( FileName, ToyResult );
	}
}

void ResultFormatter::FlatNTuplePullPlots( string Filename, FitResultVector* ToyResult )
{
	ResultFormatter::WriteFlatNtuple( Filename, ToyResult );
}

//Make pull plots from the output of a toy study
void ResultFormatter::WriteFlatNtuple( string FileName, FitResultVector* ToyResult )
{
	TFile * rootFile = new TFile( FileName.c_str(), "RECREATE" );
	//	Ntuples are 'stupid' objects in ROOT and the basic one can only handle float type objects
	TNtuple * parameterNTuple;
	parameterNTuple = new TNtuple("RapidFitResult", "RapidFitResult", ToyResult->GetFlatResultHeader());
	Float_t * resultArr;
	for( int resultIndex = 0; resultIndex < ToyResult->NumberResults(); ++resultIndex )
	{
		vector<double> result = ToyResult->GetFlatResult(resultIndex);
		resultArr = new Float_t [result.size()];
		copy( result.begin(), result.end(), resultArr);
		parameterNTuple->Fill( resultArr );
		delete [] resultArr;
	}
	rootFile->Write();
	rootFile->Close();
	//	THIS SHOULD BE SAFE... BUT THIS IS ROOT so 'of course' it isn't...
	//delete parameterNTuple;
	delete rootFile;
}

//Make pull plots from the output of a toy study
void ResultFormatter::SeparateParameterPullPlots( string FileName, FitResultVector* ToyResult )
{
	TFile * rootFile = new TFile( FileName.c_str(), "RECREATE" );
	string header = "value:error:pull";
	vector<string> allNames = ToyResult->GetAllNames();
	Float_t valueErrorPull[3];
	vector<double> parameterValues, parameterErrors, parameterPulls;
	TH1F * pullHistogram=NULL;
	TNtuple * parameterNTuple=NULL;

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
			valueErrorPull[0] = Float_t(parameterValues[unsigned(resultIndex)]);
			valueErrorPull[1] = Float_t(parameterErrors[unsigned(resultIndex)]);
			valueErrorPull[2] = Float_t(parameterPulls[unsigned(resultIndex)]);

			parameterNTuple->Fill(valueErrorPull);

			if (makeHistogram)
			{
				pullHistogram->Fill( parameterPulls[unsigned(resultIndex)] );
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

//      Pass this the pointer to the ntuple you're handing and it will give you a list
//      of branches that it contains in an easy to handle manner
vector<TString> ResultFormatter::get_branch_names( TTree* local_tree )
{
	//      To be populated and returned to the user
	vector<TString> temp_branch_names;

	//      Get the list of branches within the TTree
	TObjArray* branch_obj_array = local_tree->GetListOfBranches();

	//      Loop over all found branch objects and request their names
	for( unsigned short int i=0; i < branch_obj_array->GetEntries() ; ++i )
	{
		TObject* branch_object = (*branch_obj_array)[i];
		temp_branch_names.push_back((const char*) branch_object->GetName());
	}

	//      Return the vector of names I have found
	return temp_branch_names;
}

