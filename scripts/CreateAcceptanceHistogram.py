#! /usr/env python
import ROOT
from optparse import OptionParser
from array import array
from math import cos, pi

__author__ = 'Greig A Cowan'
__date__ = '1/11/2009'

'''
Takes LHCb MC data ntuple and a toy generated data ntuple from Bs2JpsiPhi
RapidFit PDF and creates an acceptance surface. Both data sets should use
the same observable names, as defined by the PDF.
'''

def createPyl( leaves ):
    pyl = PyListOfLeaves()
    for i in range(leaves.GetEntries()):
       	leaf = leaves.At(i)
       	name = leaf.GetName()
       	pyl.__setattr__(name,leaf)
    return pyl

def fillHistogram( pyl, ntuple, histo, weight, start, end ):
    for event in range(start, end):
	ntuple.GetEntry(event)
	cosTheta = pyl.cosTheta.GetValue()
	phi = pyl.phi.GetValue()
	cosPsi = pyl.cosPsi.GetValue()
	time = pyl.time.GetValue()
	#histo.Fill( array("f", [ cosTheta, phi, cosPsi, time ]) )
	histo.Fill( cosTheta, phi, cosPsi, weight )

def fillHistogramTime( pyl, ntuple, histo, weight, start, end ):
    for event in range(start, end):
	ntuple.GetEntry(event)
	time = pyl.time.GetValue()
	histo.Fill( time, weight )

def main():
    parser = OptionParser(
        usage = 'usage: %prog MC09_data.root Toy_data.root acceptance_histo.root')

    (options, args) = parser.parse_args()
    mcDataName = args[0]
    toyDataName = args[1]
    acceptanceName = args[2]

    mcFile = ROOT.TFile.Open(mcDataName)
    toyFile = ROOT.TFile.Open(toyDataName)
    acceptanceFile = ROOT.TFile( acceptanceName, "RECREATE")

    mcNtuple = mcFile.Get("ntuple")
    toyNtuple = toyFile.Get("dataNTuple")

    try:
	mcLeaves = mcNtuple.GetListOfLeaves()
	toyLeaves = toyNtuple.GetListOfLeaves()
    except Exception, e:
	print e

    mcData = createPyl( mcLeaves )
    toyData = createPyl( toyLeaves )

    # convention is cosTheta, phi, cosPsi, time
    bins = array("i", [20, 20, 20, 10])
    xmin = array("f", [-1., -pi, -1., -0.5])
    xmax = array("f", [ 1.,  pi,  1.,  15.])
    #mcHisto = ROOT.THnSparseD("mc", "mc", 4, bins, min, max)
    #toyHisto = ROOT.THnSparseD("toy", "toy", 4, bins, min, max)

    numMCEvents = mcNtuple.GetEntries()
    numToyEvents = toyNtuple.GetEntries()

    numHistosToProduce = 1

    #toyWeight = float(numMCEvents)/numToyEvents

    #print "There are %d MC events and %d toy events. Weight will be %f" % (numMCEvents, numToyEvents, toyWeight)

    acceptanceHistos = []
    acceptanceTimeHistos = []
    mcHistos = []
    toyHistos = []
    mcHistosTime = []
    toyHistosTime = []
    startMC = 0
    startToy = 0
    for i in range(numHistosToProduce):
	num = str(i)
    	mcHistos.append(ROOT.TH3D("mc"+num, "mc"+num, 20, -1., 1., 20, -pi, pi, 20, -1., 1.))
    	toyHistos.append(ROOT.TH3D("toy"+num, "toy"+num, 20, -1., 1., 20, -pi, pi, 20, -1., 1.))
    	mcHistos[i].Sumw2()
    	toyHistos[i].Sumw2()
    	mcHistosTime.append(ROOT.TH1D("mcTime"+num, "mcTime"+num, 20, -2., 18.))
    	toyHistosTime.append(ROOT.TH1D("toyTime"+num, "toyTime"+num, 20, -2., 18.))
    	mcHistosTime[i].Sumw2()
    	toyHistosTime[i].Sumw2()

	fillHistogram( mcData, mcNtuple, mcHistos[i], 1., startMC, startMC + numMCEvents/numHistosToProduce)
    	fillHistogram( toyData, toyNtuple, toyHistos[i], 1., startToy, startToy + numToyEvents/numHistosToProduce)

    	fillHistogramTime( mcData, mcNtuple, mcHistosTime[i], 1., startMC, startMC + numMCEvents/numHistosToProduce )
    	fillHistogramTime( toyData, toyNtuple, toyHistosTime[i], 1., startToy, startToy + numToyEvents/numHistosToProduce )

	startMC += numMCEvents/numHistosToProduce
	startToy += numToyEvents/numHistosToProduce

    	#mcHisto.Scale(1./mcHisto.Integral())
    	#toyHisto.Scale(1./toyHisto.Integral())

    	# Create the 4D histogram for the acceptance surface
    	acceptanceHistos.append(mcHistos[i].Clone("acceptance"+num))
    	acceptanceHistos[i].SetTitle("acceptance"+num)
    	acceptanceHistos[i].Divide(toyHistos[i])
    	acceptanceHistos[i].Scale(1./acceptanceHistos[i].Integral())

    	acceptanceTimeHistos.append(mcHistosTime[i].Clone("acceptanceTime"+num))
    	acceptanceTimeHistos[i].SetTitle("acceptanceTime"+num)
    	acceptanceTimeHistos[i].Divide(toyHistosTime[i])
    	acceptanceTimeHistos[i].Scale(1./acceptanceTimeHistos[i].Integral())

    	# Let's see the 1D projections for each of the angles
    	accCosTheta = acceptanceHistos[i].Project3D("x")
    	accCosTheta.SetTitle("cosTheta"+num)
    	accCosTheta.SetName("cosTheta"+num)
    	accPhi = acceptanceHistos[i].Project3D("y")
    	accPhi.SetTitle("phi"+num)
    	accPhi.SetName("phi"+num)
    	accCosPsi = acceptanceHistos[i].Project3D("z")
    	accCosPsi.SetTitle("cosPsi"+num)
    	accCosPsi.SetName("cosPsi"+num)

    	# Now let's see the 2D projections
    	accCosThetaPhi = acceptanceHistos[i].Project3D("xy")
    	accCosThetaPhi.SetTitle("phi vs cosTheta"+num)
    	accCosThetaPhi.SetName("phi_cosTheta"+num)
    	accPhiCosPsi = acceptanceHistos[i].Project3D("yz")
    	accPhiCosPsi.SetTitle("cosPsi vs phi"+num)
    	accPhiCosPsi.SetName("cosPsi_phi"+num)
    	accCosPsiCosTheta = acceptanceHistos[i].Project3D("zx")
    	accCosPsiCosTheta.SetTitle("cosTheta vs cosPsi"+num)
    	accCosPsiCosTheta.SetName("cosTheta_cosPsi"+num)

    acceptanceFile.Write()


    mcFile.Close()
    toyFile.Close()
    acceptanceFile.Close()

if __name__ == "__main__":
	class PyListOfLeaves(dict) :
        	pass
	main()
