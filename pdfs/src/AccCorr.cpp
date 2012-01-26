#include "AccCorr.h"

AccCorr::AccCorr():
useAccFit(false),histo(),nxbins(0),nybins(0),nzbins(0),xmin(0),xmax(0),ymin(0),
 ymax(0),zmin(0),zmax(0),deltax(0),deltay(0),deltaz(0),total_num_entries(0)
,angNorm(0.),w1(0.0),w2(0.0),w3(0.0),w4(0.0),w5(0.0),w6(0.0)
{
}

AccCorr::~AccCorr(){
}

AccCorr::AccCorr( const AccCorr& ACclass): 
w1(ACclass.w1),w2(ACclass.w2),w3(ACclass.w3),w4(ACclass.w4),w5(ACclass.w5),w6(ACclass.w6),histo(ACclass.histo)
,useAccFit(ACclass.useAccFit),nxbins(ACclass.nxbins),nybins(ACclass.nybins),nzbins(ACclass.nzbins),xmin(ACclass.xmin),xmax(ACclass.xmax)
,ymin(ACclass.ymin),ymax(ACclass.ymax),zmin(ACclass.zmin),zmax(ACclass.zmax),deltax(ACclass.deltax),deltay(ACclass.deltay),deltaz(ACclass.deltaz)
,total_num_entries(ACclass.total_num_entries),angNorm(ACclass.angNorm)
{
}

void AccCorr::AccCorrInit(PDFConfigurator* config) {

	//Find name of histogram needed to define 3-D angular acceptance distribution
	string fileName = config->getConfigurationValue( "AccHist" ) ;
	
	//Initialise depending upon whether configuration parameter was found
	if( fileName == "" ) {
		cout << "   No acceptance histogram was found, will use fit instead" << endl ;
		useAccFit = true ;
	}
	else {
		cout << "   Acceptance histogram requested: " << fileName << endl ;
		useAccFit = false ;

		//File location
		ifstream input_file;
		input_file.open( fileName.c_str(), ifstream::in );
		input_file.close();
		
		if( !getenv("RAPIDFITROOT") && input_file.fail() )
		{
			cerr << "\n\n" << endl;
			cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
			cerr << "$RAPIDFITROOT NOT DEFINED, PLEASE DEFINE IT SO I CAN USE ACCEPTANCE DATA" << endl;
			cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
			cerr << "\n\n" << endl;
			exit(-987);
		}
		string fullFileName;
		
		if( getenv("RAPIDFITROOT") )
		{
			string path( getenv("RAPIDFITROOT") ) ;
			
			cout << "RAPIDFITROOT defined as: " << path << endl;
			
			fullFileName = path+"/pdfs/configdata/"+fileName ;
		} else if( !input_file.fail() )
		{
			fullFileName = fileName;
		} else {
			cerr << "Shouldn't end up Here in the code!" << endl;
			exit(-892);
		}
		cout << "Found file successfully: " << fullFileName.c_str() << endl;	
		
		//Read in histo
		TFile* f =  TFile::Open(fullFileName.c_str());
		histo = new TH3D(  *(    (TH3D*)f ->Get("fit3d")      )     ); //(fileName.c_str())));
		// *********************************** CONSISTENCY CHECK ********************************
		// xaxis
		TAxis* xaxis = histo->GetXaxis();
        	xmin = xaxis->GetXmin();
        	xmax = xaxis->GetXmax();
        	nxbins = histo->GetNbinsX();
        	deltax = (xmax-xmin)/nxbins;
		// yaxis
        	TAxis* yaxis = histo->GetYaxis();
        	ymin = yaxis->GetXmin();
        	ymax = yaxis->GetXmax();
        	nybins = histo->GetNbinsY();
        	deltay = (ymax-ymin)/nybins;
		// zaxis
        	TAxis* zaxis = histo->GetZaxis();
        	zmin = zaxis->GetXmin();
        	zmax = zaxis->GetXmax();
       	 	nzbins = histo->GetNbinsZ();
        	deltaz = (zmax-zmin)/nzbins;

		//method for Checking whether histogram is sensible
		total_num_entries = histo->GetEntries();
		int total_num_bins = nxbins * nybins * nzbins;
		double sum = 0.;
		vector<int> zero_bins;

		//loop over each bin in histogram and print out how many zero bins there are
		for (int i=1; i < nxbins+1; i++){
			for (int j=1; j < nybins+1; j++){
				for (int k=1; k < nzbins+1; k++){

					double bin_content = histo->GetBinContent(i,j,k);
					//cout << "Bin content: " << bin_content << endl;
					if(bin_content<=0) { zero_bins.push_back(1);}
					//cout << " Zero bins " << zero_bins.size() << endl;}
					else if (bin_content>0){
					sum +=  bin_content;}
				}
			}
		}
		
		double average_bin_content = sum / total_num_bins;

		cout << "\n\n\t\t" << "****" << "For total number of bins " << total_num_bins << " there are " << zero_bins.size() << " bins containing zero events " << "****" << endl;
		cout <<  "\t\t\t" << "***" << "Average number of entries of non-zero bins: " << average_bin_content << "***" << endl;
		cout << endl;
		cout << endl;

		if ((xmax-xmin) < 2. || (ymax-ymin) < 2. || (zmax-zmin) < 2.*TMath::Pi() ){
			cout << "In TimeIntegratedBkg_3Dangular_Helicity::TimeIntegratedBkg_3Dangular_Helicity: The full angular range is not used in this histogram - the PDF does not support this case" << endl;
			exit(1);
		}

		cout << "Finishing processing histo" << endl;
		// **************************************** END OF CONSISTENCY CHECK *************************************

		// Normalise histogram
		//double AccEntries = histo->GetEntries();
	//cout << "histo (Init)=" << histo << endl;
	}

	MakeWeights(config);
	return;
}

double AccCorr::AccEval(vector<double> x) {
	double factor= 0.0;	
	if(useAccFit==true){
		double y[3];
		y[0]=x[0];y[1]=x[1];y[2]=x[2];
		Bs2PhiPhi_AccCorr FIT = Bs2PhiPhi_AccCorr();
		factor = FIT.MDF(y);
		//cout << "factor =" << factor << endl;
	}
	else{
		// Find the bin
		int globalbin = histo->FindBin(x[0], x[1], x[2]);
		double num_entries_bin = histo->GetBinContent(globalbin);
		//Angular factor normalized with phase space of histogram and total number of entries in the histogram
    		factor = num_entries_bin ;
	}
	//cout << "factor (END) =" << factor << endl;
	return factor;
}

void AccCorr::MakeWeights(PDFConfigurator* config) {
	
	// Get the sub-XML from configdata
	vector<pair<string, string> >* XMLOverrideList = new vector<pair<string,string> >;	
	XMLConfigReader* xmlFileACC = new XMLConfigReader( config->getConfigurationValue( "ACCEPXMLname" ), XMLOverrideList);
	
	// Extra code to use weights
	FitFunctionConfiguration * theFunction=NULL;
	OutputConfiguration * makeOutput = NULL;
	makeOutput = xmlFileACC->GetOutputConfiguration();
	theFunction = xmlFileACC->GetFitFunctionConfiguration();
	
	// If weights were specified then we need to let the output plotting know
	if( theFunction->GetWeightsWereUsed() ) makeOutput->SetWeightsWereUsed( theFunction->GetWeightName() ) ;
	
	// Get the PDF and dataset from our sub-XML
	PDFWithData * pdfAndData = xmlFileACC->GetPDFsAndData()[0];
	pdfAndData->SetPhysicsParameters( xmlFileACC->GetFitParameters() );
	IDataSet * dataSet = pdfAndData->GetDataSet();
	IPDF * PDF = pdfAndData->GetPDF();	
	PhaseSpaceBoundary * boundary = dataSet->GetBoundary();
	
	// Initialise our variables
	const int numAngularTerms = 6;
	double*  f = new double[numAngularTerms]; // the angular functions
	double xi[numAngularTerms]; // the angular weights
	double ct1, phi, ct2, time; (void) time; double weightval=0.0;
	double evalPDFraw, evalPDFnorm, val;
	int numEvents = dataSet->GetDataNumber();
	for (int i = 0; i < numAngularTerms; i++) xi[i] = 0.0;

	double Sum[numAngularTerms];
	double Sum_sq[numAngularTerms][numAngularTerms];
	double cov[numAngularTerms][numAngularTerms];
	for (int i = 0; i < numAngularTerms; i++) {
		Sum[i] = 0.0;
		for (int k = 0; k < numAngularTerms; k++) {
			cov[i][k] = 0.0;
			Sum_sq[i][k] = 0.0;
		}
	}
	
	for (int e = 0; e < numEvents; e++)
	{
		if (e % 1000 == 0) cout << "Event # " << e << endl;
		
		// Get the observables for the event
		DataPoint * event = dataSet->GetDataPoint(e);
		ct1 = event->GetObservable("ctheta_1")->GetValue();
		phi      = event->GetObservable("phi")->GetValue();
		ct2   = event->GetObservable("ctheta_2")->GetValue();
		time     = event->GetObservable("time")->GetValue();
		
		// ONLY PHI PHI SPECIFIC PART
		Mathematics::getBs2PhiPhiAngularFunctions(f[0], f[1], f[2], f[3], f[4], f[5], ct1, ct2, phi);
		
		// The method sums the angular factor f_i divided by the sum_i(A_i*f_i)
		// for each accepted event. I'm not sure if dividing by Evaluate is exactly the
		// same here, particularly if you look at untagged events.
		evalPDFraw = PDF->Evaluate( event );
		// Now need to calculate the normalisation when integrated over the 3 angles
		vector<string> dontIntegrate = PDF->GetDoNotIntegrateList();
		dontIntegrate.push_back("time");
		evalPDFnorm = PDF->Integral (event,boundary);
	
		// Code to make the weights
		if(theFunction->GetWeightsWereUsed()) weightval = event->GetObservable(theFunction->GetWeightName())->GetValue();
		if( (theFunction->GetWeightsWereUsed()) && (weightval > 0.0))
			val = evalPDFraw/evalPDFnorm*weightval;
		else {
			val = evalPDFraw/evalPDFnorm;
		}

		for (int i = 0; i < numAngularTerms; i++)
		{
			Sum[i] = Sum[i] + f[i]/val;
			xi[i] += f[i]/val;
			for (int k = 0; k < numAngularTerms; k++)
			{
			Sum_sq[i][k] += f[i]/val*f[k]/val;
			}
		}
	}

	double AvNonZero = (xi[0]+xi[1]+xi[2])/3./numEvents;
	// Print out the results and assign the weights we will use in PDF
	for (int i = 0; i < numAngularTerms; i++)
	{
		for (int k = 0; k < numAngularTerms; k++)
		{
			cov[i][k] = 1./numEvents/numEvents * ( Sum_sq[i][k] - Sum[i]*Sum[k]/numEvents);
			cout << cov[i][k]/AvNonZero/AvNonZero << "\t";
		}
		cout << endl;
	}
	cout << "Weights:" << endl;
	for (int i = 0; i < numAngularTerms; i++)
	{
			cout << xi[i]/numEvents/AvNonZero << " \\pm " << sqrt(cov[i][i])/AvNonZero << endl;
	}	
	w1 = xi[0]/numEvents/AvNonZero;
	w2 = xi[1]/numEvents/AvNonZero;
	w3 = xi[2]/numEvents/AvNonZero;
	w4 = xi[3]/numEvents/AvNonZero;
	w5 = xi[4]/numEvents/AvNonZero;
	w6 = xi[5]/numEvents/AvNonZero;
	
	return;	

}

vector<double> AccCorr::ReturnWeights() {
	vector<double> retVect;
	retVect.push_back(w1);
	retVect.push_back(w2);
	retVect.push_back(w3);
	retVect.push_back(w4);
	retVect.push_back(w5);
	retVect.push_back(w6);
	return retVect;
}
