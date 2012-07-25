/**
 @class SlicedAcceptance
 
 A class for holding a sliced propertime acceptance
 
 @author Pete Clarke
 @data 2011-06-07
 */


#include "SlicedAcceptance.h"

//............................................
// Constructor for flat acceptance
SlicedAcceptance::SlicedAcceptance( double tl, double th ) :
  slices(), nullSlice(new AcceptanceSlice(0.,0.,0.)), tlow(tl), thigh(th), beta(0)
{
	
	//Reality checks
	if( tl > th ) { cout << "SlicedAcceptance::SlicedAcceptance :  tl > th  - exiting" << endl ; exit(1) ; }
	
	//Create single slice with acceptance=1.0
	slices.push_back( new AcceptanceSlice( tlow, thigh, 1.0 ) ) ;
		
	//....done.....
	
}

SlicedAcceptance::SlicedAcceptance( const SlicedAcceptance& input ) : slices(), nullSlice( new AcceptanceSlice(0.,0.,0.) ), tlow( input.tlow ), thigh( input.thigh ), beta( input.beta )
{
	for( unsigned int i=0; i< input.slices.size(); ++i )
	{
		slices.push_back( new AcceptanceSlice( *(input.slices[i]) ) );
	}
}

//............................................
// Constructor for simple upper time acceptance only
SlicedAcceptance::SlicedAcceptance( double tl, double th, double b ) :
 slices(), nullSlice(new AcceptanceSlice(0.,0.,0.)), tlow(tl), thigh(th), beta(b)
{

	//Reality checks
	if( tlow > thigh ) { cout << "SlicedAcceptance::SlicedAcceptance :  tl > th  - exiting" << endl ; exit(1) ; }
	if( beta < 0. ) { cout << "SlicedAcceptance::SlicedAcceptance :  beta < 0  - exiting" << endl ; exit(1) ; }
		
	// The upper time acceptance formula is (1 - beta*time) so determine highest and lowest acceptance 
	double accUp = 1.0 - beta*tlow ;
	double accLow = 1.0 - beta*thigh ;
	
	
	//Create N more, equispaced 
	int N = 15 ;
	double dh = (accUp-accLow) / N ;
	double dt = (thigh - tlow) / N ;

	//First slice which goes underneath the whole slope at the accLow level +dh/2
	slices.push_back( new AcceptanceSlice( tlow, thigh, (accLow+(dh/2.0)) ) ) ;

	for( int is=1; is <= (N-1) ; is++ )    //The N-1 is important as the N one is included as the +dh/2 in first slice
	{
		double th = tlow + is*dt  ;
		slices.push_back( new AcceptanceSlice( tlow, th, dh ) ) ;
	}
		
	//....done.....
	
}


//............................................
// Constructor for simple 2010 version of lower time acceptance only
SlicedAcceptance::SlicedAcceptance( string s ) :
slices(), nullSlice(new AcceptanceSlice(0.,0.,0.)), tlow(), thigh(), beta()
{
	(void)s;	
	int N = 31;
	
	double ts[31] = {0, 0.05,  0.1, 0.15, 0.2, 0.25, 0.3, 0.34, 0.38, 0.42, 0.46, 0.5, 0.54, 0.58, 0.62, 0.66, 
		0.7,  0.75,  0.8, 0.85 , 0.9, 0.95, 1.0, 1.05, 1.1,  1.2, 1.3, 1.4,1.6 , 1.8, 2};

	double ac[31] = {0.000122069 , 0.00436829 , 0.022288 , 0.0633692 , 0.132094 , 0.225401 , 0.321979 , 0.409769 , 
		0.494438 , 0.570185 , 0.637997 , 0.695043 , 0.744576 , 0.784799 , 0.816805 , 0.844606 , 0.869917 , 
		0.891761 , 0.910512 , 0.926418 , 0.937439 , 0.947105 , 0.953599 , 0.960726 , 0.968278 , 0.976149 , 
		0.980624 , 0.98685 , 0.991765 , 0.995575 , 1 };
	
		
	for( int is=0; is < N ; is++ ) 
	{
		double tlow = ts[is] ;
		double thigh = 14.0 ;
		double height = is > 0 ? ac[is] - ac[is-1] : ac[is] ;
		slices.push_back( new AcceptanceSlice( tlow, thigh, height ) ) ;
	}
	
	//....done.....
	
}


//............................................
// Constructor for accpetance from a file
SlicedAcceptance::SlicedAcceptance( string type, string fileName ) :
slices(), nullSlice(new AcceptanceSlice(0.,0.,0.)), tlow(), thigh(), beta()
{
	(void)type;	
	
	if( type != "File" ) {   }//do nothing for now  

	ifstream input_file;

	input_file.open( fileName.c_str(), ifstream::in );
	input_file.close();

	bool local_fail = input_file.fail();

	if( !getenv("RAPIDFITROOT") && local_fail )
	{
		cerr << "\n" << endl;
		//cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		cerr << "$RAPIDFITROOT NOT DEFINED, PLEASE DEFINE IT SO I CAN USE ACCEPTANCE DATA" << endl;
		//cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		cerr << "\n" << endl;
		exit(-987);
	}

	string fullFileName;

	if( getenv("RAPIDFITROOT") )
	{
		string path( getenv("RAPIDFITROOT") ) ;

		cout << "RAPIDFITROOT defined as: " << path << endl;
	
		fullFileName = path+"/pdfs/configdata/"+fileName ;

		input_file.open( fullFileName.c_str(), ifstream::in );
		input_file.close();
	}
	bool elsewhere_fail = input_file.fail();

	if( elsewhere_fail && local_fail )
	{
		cerr << "\n\tFILE NAMED:\t" << fullFileName << "\t NOT FOUND PLEASE CHECK YOUR RAPIDFITROOT" << endl;
		exit(-89);
	}

	if( fullFileName.empty() ) fullFileName = fileName; 

	ifstream in;
	in.open(fullFileName.c_str());
	if( in.fail() ) { cout << "SlicedAcceptance::SlicedAcceptance : failed to open acceptance file  '  " << fullFileName  << "  '  " << endl ; exit(1) ; }

	std::vector<double> lowEdge; std::vector<double> hiEdge; std::vector<double> binContent;
	int ok = true;
	while (ok) {  
		lowEdge.push_back(stream(in));  hiEdge.push_back(stream(in)); binContent.push_back(stream(in) );
		if (in.eof()  == true || in.good() == false) ok =false; 
	}		
	in.close();
	
	for( unsigned int is=0; is < lowEdge.size() ; is++ ) 
	{
		double tlow = lowEdge[is];
		double thigh = hiEdge[is] ;
		double height = binContent[is]  ;
		//cout << " Adding slice " << tlow << " /  " << thigh << " /  " << height << endl ;
		slices.push_back( new AcceptanceSlice( tlow, thigh, height ) ) ;
	}
	
	//....done.....
	
}



//............................................
// Return numerator for evaluate
double SlicedAcceptance::getValue( double t ) const
{
	double returnValue = 0 ;
	for( unsigned int is = 0; is < slices.size() ; ++is )
	{
		if( (t>=slices[is]->tlow()) && (t<slices[is]->thigh()) ) returnValue += slices[is]->height() ;
	}
	return returnValue;
}

double SlicedAcceptance::getValue( Observable* time, double timeOffset ) const
{
	if( time->GetBinNumber() != -1 )
	{
//#pragma GCC diagnostic ignored "-Wfloat-equal"
		if( time->GetOffSet() == timeOffset ) return time->GetAcceptance();
//#pragma GCC diagnostic pop
	}

	double t= time->GetValue() - timeOffset;
	double returnValue = 0;
	unsigned int is = 0;
       	for( ; is < slices.size() ; ++is )
	{
		if( (t>=slices[is]->tlow()) && (t<slices[is]->thigh()) ) returnValue += slices[is]->height();
       	}
	time->SetBinNumber( (int)is );
	time->SetAcceptance( returnValue );
	time->SetOffSet( timeOffset );
	return returnValue;

}

//............................................
// Return the number of slices
unsigned int SlicedAcceptance::numberOfSlices() const
{
	return (unsigned)slices.size();
}

//............................................
// Return a slice
AcceptanceSlice * SlicedAcceptance::getSlice( unsigned int s ) const
{
	if( s>=(unsigned)slices.size() ) return nullSlice;
	return slices[s];
}

//............................................
//Helpers
double SlicedAcceptance::stream(ifstream& stream)
{
	double tmpVal;
	stream >> tmpVal;
	return tmpVal;
}

