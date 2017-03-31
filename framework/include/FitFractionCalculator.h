#pragma once
#include <map>
#include "TFile.h"
#include "IPDF.h"

class FitFractionCalculator
{
	public:
		FitFractionCalculator(IPDF& pdf, PhaseSpaceBoundary& boundary);
		std::map<std::string, double> GetFitFractions(DataPoint* combination) const;
		double GetComponentFitFraction(DataPoint* combination, std::string ComponentName) const;
		void Print() const;
		void WriteToFile(std::string filename);
	private:
		FitFractionCalculator(const FitFractionCalculator&); // Don't copy. All the work is done in the constructor
		std::map<int, std::map<std::string, double>> ComponentFitFractions;
		int GetIndex(DataPoint* combination) const;
		PhaseSpaceBoundary storedBoundary;
};

