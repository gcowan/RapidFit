/*!
 *
 * @brief This Class Contains all of the information to configure the plots and Projections which should be performed
 *
 * @author Robert Currie rcurrie@cern.ch
 */

#pragma once
#ifndef COMP_PLOTTER_CONFIG_H
#define COMP_PLOTTER_CONFIG_H

///	ROOT Headers
#include "TString.h"
///	RapidFit Headers
#include "RapidFitIntegratorConfig.h"
///	System Headers
#include <vector>
#include <string>

using namespace::std;

/*!
 * This was originally a struct and I passed it around several places, but rootcint hates things like vector<struct something*> so it's a dummy class
 * This can be replaced by an interface function with get/set methods, but this class is defined to exist ONLY to provide standard objects as input,
 * it is not expected to do any work itself
 */
class CompPlotter_config
{
        public:
                /*!
                 * @brief Explicit Constructor containing the hard-coded defaults
                 */
                CompPlotter_config() :
                        data_bins(100), PDF_points(128), observableName("undefined"), logY(false), logX(false), color_key(), style_key(), width_key(), component_names(), PlotTitle(""),
                        xmin(-99999), xmax(-99999), ymin(-99999), ymax(-99999), xtitle(""), ytitle(""), CalcChi2(false), Chi2Value(-99999), OnlyZero(false), ScaleNumerical(true), combination_names(),
			DrawPull(false), LegendTextSize(0.05), addLHCb(false), TopRightLegend(true), TopLeftLegend(false), BottomRightLegend(false), BottomLeftLegend(false),
			useLegend(true), LimitPulls(false), useSpline(true), addRightLHCb(false), integratorConfig(new RapidFitIntegratorConfig()), plotAllCombinations(true), defaultCombinationValue(0.),
			XaxisTitleScale(1.), XaxisLabelScale(1.), YaxisTitleScale(1.), YaxisLabelScale(1.)
		{}

		~CompPlotter_config()
		{
			if( integratorConfig != NULL ) delete integratorConfig;
		}

		CompPlotter_config( const CompPlotter_config& input ) :
			data_bins(input.data_bins), PDF_points(input.PDF_points), observableName(input.observableName), logY(input.logY), logX(input.logX), color_key(input.color_key),
			style_key(input.style_key), width_key(input.width_key), component_names(input.component_names), PlotTitle(input.PlotTitle), xmin(input.xmin), xmax(input.xmax),
			ymin(input.ymin), ymax(input.ymax), xtitle(input.xtitle), ytitle(input.ytitle), CalcChi2(input.CalcChi2), Chi2Value(input.Chi2Value), OnlyZero(input.OnlyZero),
			ScaleNumerical(input.ScaleNumerical), DrawPull(input.DrawPull), LegendTextSize(input.LegendTextSize), addLHCb(input.addLHCb), TopRightLegend(input.TopRightLegend),
			TopLeftLegend(input.TopLeftLegend), BottomRightLegend(input.BottomRightLegend), BottomLeftLegend(input.BottomLeftLegend), useLegend(input.useLegend),
			LimitPulls(input.LimitPulls), useSpline(input.useSpline), addRightLHCb(input.addRightLHCb), integratorConfig(NULL), combination_names(input.combination_names),
			plotAllCombinations(input.plotAllCombinations), defaultCombinationValue(input.defaultCombinationValue),
			XaxisTitleScale(input.XaxisTitleScale), XaxisLabelScale(input.XaxisLabelScale), YaxisTitleScale(input.YaxisTitleScale), YaxisLabelScale(input.YaxisLabelScale)
		{
			if( input.integratorConfig != NULL ) integratorConfig = new RapidFitIntegratorConfig( *(input.integratorConfig) );
		}

                int data_bins;                  /*!     This is the total number of bins which should be used when binning the dataset                                          */
                int PDF_points;                 /*!     This is the number of points which the PDF should be interrogated for across the whole Observable range                 */
                string observableName;          /*!     This is the Name of the Observable that is being Projected                                                              */
                bool logY;                      /*!     This is the boolean governing if the Y axis should be on a log scale                                                    */
		bool logX;			/*!	This is the boolean governing if the X axis shoulw be on a log scale							*/
                vector<int> color_key;          /*!     This is a vector of numbers, each number corresponding to the Color_t which is given to the corresponding component     */
                vector<int> style_key;          /*!     This is a vector of numbers, each number corresponding to the Style_t which is given to the corresponding component     */
                vector<int> width_key;          /*!     This is a vector of numbers, each number corresponding to the Width_t which is given to the corresponding component     */
                vector<string> component_names; /*!     This ia a vector of strings which are names passed to the TLegend for all of the Components                             */
		vector<string> combination_names;/*!	This is a vector of strings which are names passed to a seperate TLegend for the different Combinations			*/
                double LegendTextSize;		/*!	This chaneges the ROOT text size in the Legend box									*/
                string PlotTitle;               /*!     This is the Name of the Plot                                                                                            */
                double xmin, xmax, ymin, ymax;  /*!     These are the ranges that should be plotted on the X and Y axis                                                         */
                TString xtitle, ytitle;         /*!     These are the titles of the X and Y axis                                                                                */
                bool CalcChi2;                  /*!     This boolean lets the class know you want the Chi2 value calculated between this DataSet and PDF                        */
                double Chi2Value;               /*!     This is where the Chi2 final value is stored after it's calculated                                                      */
                bool OnlyZero;                  /*!     If true this class will mimic the correct behaviour of the old Plotter class                                            */
		bool ScaleNumerical;		/*!	Do you scale the Numerical or Analytical Integral                                                                       */
		bool addLHCb;			/*!	Add an LHCb 'Stamp' to the plot												*/
		bool addRightLHCb;		/*!	Add an LHCb 'Stamp' to the plot on the right										*/
		bool DrawPull;			/*!	Should I draw the Pull Plot from this projection over the data?								*/
		bool TopRightLegend;		/*!	Use a Top Right Legend Box												*/
		bool TopLeftLegend;		/*!	Use a Top Left Legend Box												*/
		bool BottomRightLegend;		/*!	Use a Bottom Right Legend Box												*/
		bool BottomLeftLegend;		/*!	Use a Bottom Left Legend Box												*/
		bool useLegend;			/*!	Use a Legend in the Plot												*/
		bool LimitPulls;		/*!	Place a 5-sigma limit on the Pull Distribution to fight a few bogus pulls killing the Distribution			*/
		bool useSpline;			/*!	Use Spline to Interpolate between the projected points?									*/
		RapidFitIntegratorConfig* integratorConfig;
		bool plotAllCombinations;	/*!	Plot a separate Projection for each Discrete Combination for all PDFs							*/
		double defaultCombinationValue;	/*!	Default Value for the Discrete Values to be given if plotAllCombinations==false						*/

		double XaxisTitleScale;
		double XaxisLabelScale;
		double YaxisTitleScale;
		double YaxisLabelScale;

	private:
		CompPlotter_config& operator= ( const CompPlotter_config& input );
};

#endif

