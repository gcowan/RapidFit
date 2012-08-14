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
                        data_bins(100), PDF_points(128), observableName("undefined"), logY(false), color_key(), style_key(), width_key(), component_names(), PlotTitle(""),
                        xmin(-999), xmax(-999), ymin(-999), ymax(-999), xtitle(""), ytitle(""), CalcChi2(false), Chi2Value(-999), OnlyZero(false), ScaleNumerical(true),
			DrawPull(false), LegendTextSize(0.03), addLHCb(false)
                {}


                int data_bins;                  /*!     This is the total number of bins which should be used when binning the dataset                                          */
                int PDF_points;                 /*!     This is the number of points which the PDF should be interrogated for across the whole Observable range                 */
                string observableName;          /*!     This is the Name of the Observable that is being Projected                                                              */
                bool logY;                      /*!     This is the boolean governing if the Y axis should be on a log scale                                                    */
                vector<int> color_key;          /*!     This is a vector of numbers, each number corresponding to the Color_t which is given to the corresponding component     */
                vector<int> style_key;          /*!     This is a vector of numbers, each number corresponding to the Style_t which is given to the corresponding component     */
                vector<int> width_key;          /*!     This is a vector of numbers, each number corresponding to the Width_t which is given to the corresponding component     */
                vector<string> component_names; /*!     This ia a vector of strings which are names passed to the TLegend for all of the Components                             */
                double LegendTextSize;
                string PlotTitle;               /*!     This is the Name of the Plot                                                                                            */
                double xmin, xmax, ymin, ymax;  /*!     These are the ranges that should be plotted on the X and Y axis                                                         */
                TString xtitle, ytitle;         /*!     These are the titles of the X and Y axis                                                                                */
                bool CalcChi2;                  /*!     This boolean lets the class know you want the Chi2 value calculated between this DataSet and PDF                        */
                double Chi2Value;               /*!     This is where the Chi2 final value is stored after it's calculated                                                      */
                bool OnlyZero;                  /*!     If true this class will mimic the correct behaviour of the old Plotter class                                            */
		bool ScaleNumerical;		/*!	Do you scale the Numerical or Analytical Integral                                                                       */
		bool addLHCb;
		bool DrawPull;			/*!	Should I draw the Pull Plot from this projection over the data?								*/
};

#endif
