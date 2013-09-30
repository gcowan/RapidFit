//
/**
  @class MemoryDataSet

  A data set which simply stores a vector of pointers to datapoint objects

  @author Benjamin M Wynne bwynne@cern.ch
  @date 2009-10-02
  */

//	RapidFit Headers
#include "MemoryDataSet.h"
#include "IConstraint.h"
#include "ObservableDiscreteConstraint.h"
#include "StringProcessing.h"
#include "StatisticsFunctions.h"
//	System Headers
#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <iomanip>

#define DOUBLE_TOLERANCE_DATA 1E-6

using namespace::std;

MemoryDataSet::MemoryDataSet( PhaseSpaceBoundary* NewBoundary, vector<DataPoint*> inputData ) :
	allData(), dataBoundary( new PhaseSpaceBoundary(*NewBoundary) ), allSubSets(), WeightName(""), useWeights(false), alpha(1.), alphaName("uninitialized"), canDelete(false)
{
	for( unsigned int i=0; i< (unsigned)dataBoundary->GetNumberCombinations(); ++i )
	{
		allSubSets.push_back( -1 );
	}
	for( unsigned int i=0; i< inputData.size(); ++i )
	{
		this->SafeAddDataPoint( inputData[i] );
	}
}

//Constructor with correct argument
MemoryDataSet::MemoryDataSet( PhaseSpaceBoundary* NewBoundary ) :
	allData(), dataBoundary( new PhaseSpaceBoundary(*NewBoundary) ), allSubSets(), WeightName(""), useWeights(false), alpha(1.), alphaName("uninitialized"), canDelete(true)
{
	for( unsigned int i=0; i< (unsigned)dataBoundary->GetNumberCombinations(); ++i )
	{
		allSubSets.push_back( -1 );
	}
}

//Destructor
MemoryDataSet::~MemoryDataSet()
{
	delete dataBoundary;
	if( canDelete )
	{
		while( !allData.empty() )
		{
			if( allData.back() != NULL ) delete allData.back();
			allData.pop_back();
		}
	}
}

void MemoryDataSet::ClearAllPseudoObservables()
{
	for( unsigned int i=0; i< allData.size(); ++i )
	{
		allData[i]->ClearPseudoObservable();
	}
}

//Add a data point to the set
bool MemoryDataSet::AddDataPoint( DataPoint* NewDataPoint )
{
	if( dataBoundary->IsPointInBoundary(NewDataPoint) )
	{
		allData.push_back( NewDataPoint );
		allData.back()->SetPhaseSpaceBoundary( dataBoundary );
		return true;
	}
	else
	{
		delete NewDataPoint;
		//cerr << "Data point is not within data set boundary" << endl;
		return false;
	}
}

//Add a data point to the set
void MemoryDataSet::SafeAddDataPoint( DataPoint* NewDataPoint )
{
	allData.push_back( NewDataPoint );
	allData.back()->SetPhaseSpaceBoundary( dataBoundary );
}

//Retrieve the data point with the given index
DataPoint * MemoryDataSet::GetDataPoint( int Index ) const
{
	if ( Index < int(allData.size()) )
	{
		allData[unsigned(Index)]->SetPhaseSpaceBoundary( dataBoundary );
		return allData[unsigned(Index)];
	}
	else
	{
		cerr << "Index (" << Index << ") out of range in DataSet" << endl;
	}
	return NULL;
}

//Get the number of data points in the set
int MemoryDataSet::GetDataNumber( DataPoint* templateDataPoint ) const
{
	if( templateDataPoint == NULL )	return (int)allData.size();
	else
	{
		pair< vector<ObservableRef>, vector<double > > thisPointInfo = this->GetBoundary()->GetDiscreteInfo( templateDataPoint );
		unsigned int dataPointNum = (unsigned)this->GetDiscreteSubSet( thisPointInfo.first, thisPointInfo.second ).size();
		return (int) dataPointNum;
	}
}

//Get the data bound
PhaseSpaceBoundary * MemoryDataSet::GetBoundary() const
{
	return dataBoundary;
}

void MemoryDataSet::SetBoundary( const PhaseSpaceBoundary* Input )
{
	delete dataBoundary;
	dataBoundary = new PhaseSpaceBoundary( *Input );
}

//Empty the data set
void MemoryDataSet::Clear()
{
	vector<DataPoint*> empty;
	while( !allData.empty() )
	{
		if( allData.back() != NULL ) delete allData.back();
		allData.pop_back();
	}
	allData.swap(empty);
}

void MemoryDataSet::SortBy( string parameter )
{

	cout << "Sorting" << endl;
	if( allData.size() > 0 )
	{
		vector<pair<DataPoint*,ObservableRef> > allData_sort;

		for( vector<DataPoint*>::iterator data_i = allData.begin(); data_i != allData.end(); ++data_i )
		{
			allData_sort.push_back( make_pair( *data_i, ObservableRef( parameter ) ) );
		}

		//cout << "hello" << endl;
		sort( allData_sort.begin(), allData_sort.end(), DataPoint() );
		cout << "sorted" << endl;
		//	Sort the data in memory

		while( !allData.empty() ) allData.pop_back();

		for( vector<pair<DataPoint*,ObservableRef> >::iterator sort_i = allData_sort.begin(); sort_i != allData_sort.end(); ++sort_i )
		{
			allData.push_back( sort_i->first );
		}

		cout << allData.size() << endl;
	}
	cout << "Sorted" << endl;
}

IDataSet* MemoryDataSet::GetDiscreteDataSet( const vector<ObservableRef> discreteParam, const vector<double> discreteVal ) const
{
	return (IDataSet*) new MemoryDataSet( this->GetBoundary(), GetDiscreteSubSet( discreteParam, discreteVal ) );
}

vector<DataPoint*> MemoryDataSet::GetDiscreteSubSet( const vector<ObservableRef> discreteParam, const vector<double> discreteVal ) const
{
	if( discreteParam.empty() || discreteVal.empty() )
	{
		return allData;
	}

	vector<DataPoint*> returnable_subset;
	if( discreteParam.size() != discreteVal.size() )
	{
		cerr << "\n\n\t\tBadly Defined definitio`n of a subset, returning 0 events!\n\n" << endl;
		return returnable_subset;
	}

	for( unsigned int i=0; i< allData.size(); ++i )
	{
		DataPoint* data_i = allData[i];

		bool decision = true;       
		for( unsigned int j=0; j< discreteParam.size(); ++j )
		{
			if( !( fabs( data_i->GetObservable( discreteParam[j] )->GetValue() - discreteVal[j] ) < DOUBLE_TOLERANCE_DATA ) )
			{
				decision = false;
			}
		}

		if( decision ) returnable_subset.push_back( data_i );

	}

	return returnable_subset;
}


vector<DataPoint*> MemoryDataSet::GetDiscreteSubSet( const vector<string> discreteParam, const vector<double> discreteVal ) const
{
	if( discreteParam.empty() || discreteVal.empty() )
	{
		return allData;
	}

	vector<DataPoint*> returnable_subset;
	if( discreteParam.size() != discreteVal.size() )
	{
		cerr << "\n\n\t\tBadly Defined definition of a subset, returning 0 events!\n\n" << endl;
		return returnable_subset;
	}

	for( unsigned int i=0; i< allData.size(); ++i )
	{
		DataPoint* data_i = allData[i];
		vector<ObservableRef*> temp_ref;
		for( unsigned int j=0; j< discreteParam.size(); ++j )
		{
			temp_ref.push_back( new ObservableRef( discreteParam[j] ) );
		}

		bool decision = true;
		for( unsigned int j=0; j< discreteParam.size(); ++j )
		{
			if( !( fabs( data_i->GetObservable( *(temp_ref[j]) )->GetValue() - discreteVal[j] ) < DOUBLE_TOLERANCE_DATA ) )
			{
				decision = false;
			}
		}

		if( decision ) returnable_subset.push_back( data_i );

		while( !temp_ref.empty() )
		{
			if( temp_ref.back() != NULL ) delete temp_ref.back();
			temp_ref.pop_back();
		}
	}

	return returnable_subset;
}

double MemoryDataSet::Yield()
{
	if( useWeights )	return this->GetSumWeights();
	else			return this->GetDataNumber();
}

double MemoryDataSet::YieldError()
{
	if( useWeights )	return sqrt( this->GetSumWeightsSq() );
	else			return sqrt( this->GetDataNumber() );
}

string MemoryDataSet::GetWeightName() const
{
	return WeightName;
}

bool MemoryDataSet::GetWeightsWereUsed() const
{
	return useWeights;
}

void MemoryDataSet::UseEventWeights( const string Name )
{
	WeightName = Name;
	useWeights = true;
	for( unsigned int i=0; i< allData.size(); ++i )
	{
		allData[i]->SetEventWeight( allData[i]->GetObservable( WeightName )->GetValue() );
	}
}

void MemoryDataSet::NormaliseWeights()
{
	if( useWeights )
	{
		double sum_Val=0.;
		double sum_Val2=0.;
		for( unsigned int i=0; i< allData.size(); ++i )
		{
			double thisVal=allData[i]->GetEventWeight();
			sum_Val += thisVal;
			sum_Val2 += thisVal*thisVal;
		}

		this->ApplyAlpha( sum_Val, sum_Val2 );
	}
}

void MemoryDataSet::Print() const
{
	if( this->GetWeightsWereUsed() && this->GetDataNumber() > 0 )
	{
		double total=0.;
		double avr=0.;
		double avr_sq=0.;
		double val=0.;
		ObservableRef weightRef( WeightName );
		double err=0.;
		for( unsigned int i=0; i< allData.size(); ++i )
		{
			val=allData[i]->GetObservable( WeightName )->GetValue();
			avr+=val;
			avr_sq+=val*val;
			err+=val*val;
		}
		total=avr;
		avr/=(double)allData.size(); avr_sq/=(double)allData.size();
		err = sqrt(err);
		//err = sqrt( avr_sq - avr*avr );
		cout << "DataSet contains a total of:     " << total << " ± " << err << "     SIGNAL events.(" << allData.size() << " total). In " << this->GetBoundary()->GetNumberCombinations() << " Discrete DataSets." << endl;
		if( this->GetBoundary()->GetNumberCombinations() > 1 )
		{
			vector<DataPoint*> combinations = this->GetBoundary()->GetDiscreteCombinations();
			for( unsigned int i=0; i< combinations.size(); ++i )
			{
				string description = this->GetBoundary()->DiscreteDescription( combinations[i] );
				for( unsigned int j=0; j< description.size(); ++j )
				{
					if( description[j] == '\n' ) description[j] = ' ';
					if( description[j] == '\t' ) description[j] = ' ';
				}
				pair< vector<ObservableRef>, vector<double > > thisPointInfo = this->GetBoundary()->GetDiscreteInfo( combinations[i] );
				unsigned int dataPointNum = (unsigned)this->GetDiscreteSubSet( thisPointInfo.first, thisPointInfo.second ).size();
				if( dataPointNum > 0 ) cout << "Combination: " << description << " has: " << dataPointNum << " events." << endl;
			}
		}
	}
	else
	{
		cout << "DataSet contains a total of:     " << allData.size() << "     events. In " << this->GetBoundary()->GetNumberCombinations() << " Discrete DataSets." << endl;
		if( this->GetBoundary()->GetNumberCombinations() > 1 )
		{
			vector<DataPoint*> combinations = this->GetBoundary()->GetDiscreteCombinations();
			for( unsigned int i=0; i< combinations.size(); ++i )
			{
				string description = this->GetBoundary()->DiscreteDescription( combinations[i] );
				for( unsigned int j=0; j< description.size(); ++j )
				{
					if( description[j] == '\n' ) description[j] = ' ';
					if( description[j] == '\t' ) description[j] = ' ';
				}
				pair< vector<ObservableRef>, vector<double > > thisPointInfo = this->GetBoundary()->GetDiscreteInfo( combinations[i] );
				unsigned int dataPointNum = (unsigned)this->GetDiscreteSubSet( thisPointInfo.first, thisPointInfo.second ).size();
				if( dataPointNum > 0 ) cout << "Combination: " << description << " has: " << dataPointNum << " events." << endl;
			}
		}
	}
}

double MemoryDataSet::GetSumWeights() const
{
	if( this->GetWeightsWereUsed() )
	{
		double total=0.;
		ObservableRef weightRef( WeightName );
		for( unsigned int i=0; i< allData.size(); ++i )
		{
			total+=allData[i]->GetObservable( WeightName )->GetValue();
		}
		return total;
	}
	else
	{
		return (double)this->GetDataNumber();
	}
}

double MemoryDataSet::GetSumWeightsSq() const
{
	if( this->GetWeightsWereUsed() )
	{
		double total=0.;
		ObservableRef weightRef( WeightName );
		for( unsigned int i=0; i< allData.size(); ++i )
		{
			total+=(allData[i]->GetObservable( WeightName )->GetValue()*allData[i]->GetObservable( WeightName )->GetValue());
		}
		return total;
	}
	else
	{
		return (double)this->GetDataNumber();
	}
}

void MemoryDataSet::ApplyAlpha( const double total_sum, const double total_sum_sq )
{
	alpha= fabs(total_sum / total_sum_sq);
	ObservableRef WeightNameRef( WeightName );
	for( unsigned int i=0; i< allData.size(); ++i )
	{
		allData[i]->SetEventWeight( allData[i]->GetObservable( WeightNameRef )->GetValue() * alpha );
	}
	cout << "alpha = " << setprecision(10) << total_sum << "  /  " << total_sum_sq << endl;
	cout << "Correction Factor: " << setprecision(5) << fabs(alpha) << " applied to DataSet containing " << allData.size() << " events." << endl << endl;
}

double MemoryDataSet::GetAlpha() const
{
	if( alphaName != "uninitialized" )
	{
		double alphaSum=0.;
		ObservableRef alphaNameRef( alphaName );
		for( unsigned int i=0; i< allData.size(); ++i )
		{
			alphaSum+=allData[i]->GetObservable( alphaNameRef )->GetValue();
		}
		alphaSum/=(double)allData.size();
		return fabs(alphaSum);
	}
	else
	{
		return alpha;
	}
}

void MemoryDataSet::ApplyExternalAlpha( const string AlphaName )
{
	alphaName = AlphaName;
	ObservableRef alphaNameRef( alphaName );
	ObservableRef WeightNameRef( WeightName );
	double avr=0.;
	for( unsigned int i=0; i< allData.size(); ++i )
	{
		double thisalpha = allData[i]->GetObservable( alphaNameRef )->GetValue();
		allData[i]->SetEventWeight( allData[i]->GetObservable( WeightNameRef )->GetValue() * thisalpha );
		avr+=thisalpha;
	}
	avr/=(double)allData.size();
	cout << "Using Observable: " << alphaName << " to apply a per-event alpha correction to the per-event weights used. Average Weight: " << avr << endl << endl;
}


