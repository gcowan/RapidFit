#ifndef MULTIDIMCHI2_H
#define MULTIDIMCHI2_H
// Standard C++
#include <string>
#include <vector>
#include <memory>
// ROOT
#include "THn.h"
// RapidFit
#include "ObservableRef.h"
#include "PhaseSpaceBoundary.h"
#include "PDFWithData.h"

class MultiDimChi2
{
	public:
		MultiDimChi2(const std::vector<PDFWithData*>& _allObjects, vector<std::string> wantedObservables);
		void PerformMuiltDimTest() const;
	private:
		// Stuff to store locally
		vector<ObservableRef> Observables;
		std::vector<std::unique_ptr<THnD>> BinnedData; // THnD is neither copyable nor moveable so I have to resort to this nonsense. Thanks ROOT!
		std::vector<PDFWithData*> allObjects;
		// Helper functions
		double CalculateExpected(IPDF& thisPDF, PhaseSpaceBoundary& fullPhaseSpace, const IDataSet& thisDataSet, const THnD& DataHist, const std::vector<int>& indices) const;
		std::vector<int> GetIndices(unsigned binNum, const THnD& DataHist) const;
		double CalcChi2(const std::vector<double>& expected_events, const std::vector<double>& observed_events, const std::vector<double>& errors) const;
};

#endif

