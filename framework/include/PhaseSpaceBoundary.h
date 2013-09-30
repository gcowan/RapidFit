/**
        @class PhaseSpaceBoundary

        A collection of constraints on observables, defining the phase space in which a data point exists

        @author Benjamin M Wynne bwynne@cern.ch
	@date 2009-10-02
*/

#pragma once
#ifndef PHASE_SPACE_BOUNDARY_H
#define PHASE_SPACE_BOUNDARY_H

//	RapidFit Headers
#include "IConstraint.h"
#include "DataPoint.h"
#include "ObservableRef.h"
//#include "IPDF.h"
//	System Headers
#include <vector>
#include <string>

using namespace::std;

class IPDF;

class PhaseSpaceBoundary
{
	public:
		PhaseSpaceBoundary( const PhaseSpaceBoundary& );
		PhaseSpaceBoundary( vector<string> );
		virtual ~PhaseSpaceBoundary();

		vector<string> GetAllNames() const;
		vector<string> GetDiscreteNames() const;
		vector<string> GetContinuousNames() const;

		pair< vector<ObservableRef>, vector<double > > GetDiscreteInfo( DataPoint* input ) const;

		bool SetConstraint( string, IConstraint* );
		bool SetConstraint( string, double, double, string );
		bool SetConstraint( string, vector<double>, string );
		void AddConstraint( string, IConstraint*, bool overwrite=false );
		IConstraint * GetConstraint( ObservableRef& ) const;
		IConstraint * GetConstraint( string ) const;
		bool IsPointInBoundary( DataPoint*, bool silence=false );

		void RemoveConstraint( string Name );

		virtual void Print() const;

		virtual string DiscreteDescription( const DataPoint* ) const;

		virtual DataPoint* GetMidPoint() const;

		virtual vector<DataPoint*> GetDiscreteCombinations() const;

		virtual int GetNumberCombinations() const;

		virtual unsigned int GetDiscreteIndex( DataPoint* Input, const bool silence=false ) const;

		void CheckPhaseSpace( IPDF* toCheck ) const;

		string XML() const;

		size_t GetID() const;

	private:
		PhaseSpaceBoundary& operator=(const PhaseSpaceBoundary&);
		vector< IConstraint* > allConstraints;
		vector<string> allNames;
		mutable int DiscreteCombinationNumber;

		mutable size_t uniqueID;

};

#endif

