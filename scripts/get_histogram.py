#!/usr/bin/env python

#python get_histograms.py something.root
# creates a time acceptance file for each histogram in the root file

import copy
import ROOT
import sys
import random
import math
import array


file = ROOT.TFile.Open(sys.argv[1])
names = [k.GetName() for k in file.GetListOfKeys() if "Error" not in k.GetName()]

for name in names:
	output = open(name+'.txt', 'w')
	histogram=file.Get(name)
	bins=histogram.GetNbinsX()
	for index in range(bins):
		bin=histogram.GetBin(index)		 
		lower=histogram.GetBinLowEdge(bin)
		width=histogram.GetBinWidth(bin)
		number=histogram.GetBinContent(bin)  
		upper=lower+width
		if bin is not 0:
			output.write(str(lower)+' '+str(upper)+' '+str(number)+'\n')
	output.close()
	print name, bins





