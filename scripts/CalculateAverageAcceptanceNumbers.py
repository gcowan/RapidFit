#! /usr/env python
import ROOT
from optparse import OptionParser
from array import array
from math import pi, cos, sin, sqrt

__author__ = 'Greig A Cowan'
__date__ = '21/01/2010'

'''
Takes 3D angular acceptance surface histogram and calculates the integrals
of the 6 angular functions across this surface. These will be fed into the fit.
'''


def angularFunctions( cosTheta, phi, cosPsi ):
	sinPhi = sin(phi)
	cosPhi = cos(phi)
	sinTheta = sqrt(1. - cosTheta*cosTheta)
	sinPsi   = sqrt(1. - cosPsi*cosPsi)
	sin2Theta = 2.*sinTheta*cosTheta
	sin2Psi   = 2.*sinPsi*cosPsi
	sin2Phi   = 2.*sinPhi*cosPhi

        f1 =  2.* cosPsi*cosPsi * ( 1. - sinTheta*sinTheta * cosPhi*cosPhi )
        f2 =      sinPsi*sinPsi * ( 1. - sinTheta*sinTheta * sinPhi*sinPhi )
        f3 =      sinPsi*sinPsi * sinTheta*sinTheta
        f4 = -1.* sinPsi*sinPsi * sin2Theta * sinPhi
        f5 =      sin2Psi * sinTheta*sinTheta * sin2Phi/sqrt(2.)
        f6 =      sin2Psi * sin2Theta * cosPhi/sqrt(2.)

	return f1, f2, f3, f4, f5, f6

def junk(acceptanceName):
	acceptanceFile = ROOT.TFile.Open(acceptanceName, "RECREATE")
	acceptanceHisto = ROOT.TH3D("flat", "flat", 40, -1., 1., 40, -pi, pi, 40, -1., 1.)
       	uniform = ROOT.TF1("uniform","1")
	acceptanceHisto.FillRandom("uniform", 20000000);
    	acceptanceHisto.Write()
	acceptanceFile.Close()

def calculate(acceptanceHisto, nameOfHisto):
    xbins = acceptanceHisto.GetXaxis().GetNbins()
    ybins = acceptanceHisto.GetYaxis().GetNbins()
    zbins = acceptanceHisto.GetZaxis().GetNbins()

    fInt = [0. for i in range(6)]
    largestNumEvents = 1
    normalisation = 1000
    # ROOT starts bin numbering from 1...
    for xbin in range(1, xbins+1):
	for ybin in range(1, ybins+1):
		for zbin in range(1, zbins+1):
			bin = acceptanceHisto.GetBin(xbin, ybin, zbin)
			numEvents = acceptanceHisto.GetBinContent(bin)
			if numEvents > largestNumEvents: largestNumEvents = numEvents
			xWidth = acceptanceHisto.GetXaxis().GetBinWidth(xbin)
			yWidth = acceptanceHisto.GetYaxis().GetBinWidth(ybin)
			zWidth = acceptanceHisto.GetZaxis().GetBinWidth(zbin)
			x = acceptanceHisto.GetXaxis().GetBinCenter(xbin)
			y = acceptanceHisto.GetYaxis().GetBinCenter(ybin)
			z = acceptanceHisto.GetZaxis().GetBinCenter(zbin)
			f = angularFunctions(x, y, z)
			for i in range(6):
				if nameOfHisto == "flat":
					numEvents = 1
					largestNumEvents = 1
					normalisation = 1
				fInt[i] += f[i]*xWidth*yWidth*zWidth*numEvents

    # Use normalisation just to get numbers to same scale as 32pi/9. It's only the
    # relative sizes of the numbers that ultimately matters in the NLL minimisation
    vals = ["%.5f" % (i/largestNumEvents*normalisation) for i in fInt]
    vals = ["%.5f" % (i/fInt[0]) for i in fInt]
    print vals


def main():
    parser = OptionParser(
        usage = 'usage: %prog acceptance_histo.root <name_of_histogram>')

    (options, args) = parser.parse_args()
    acceptanceName = args[0]
    numberOfHistos = int(args[1])
    if len(args) == 2:
	nameOfHisto = "acceptance"
    else:
	nameOfHisto = args[2]

    acceptanceFile = ROOT.TFile.Open(acceptanceName)

    for i in range(numberOfHistos):
    	acceptanceHisto = acceptanceFile.Get(nameOfHisto + str(i))
    	#acceptanceHisto = acceptanceFile.Get(nameOfHisto)
    	calculate(acceptanceHisto, nameOfHisto)

    acceptanceFile.Close()

if __name__ == "__main__":
	#junk("../../data/mc09_Bs2JpsiPhi_acceptance_flat.root")
	#import sys
	#sys.exit(1)
	main()

