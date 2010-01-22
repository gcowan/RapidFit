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

def fillHistogram( pyl, ntuple, histo, weight ):
    for event in range(ntuple.GetEntries()):
	ntuple.GetEntry(event)
	cosTheta = pyl.cosTheta.GetValue()
	phi = pyl.phi.GetValue()
	cosPsi = pyl.cosPsi.GetValue()
	time = pyl.time.GetValue()
	#histo.Fill( array("f", [ cosTheta, phi, cosPsi, time ]) ) 
	histo.Fill( cosTheta, phi, cosPsi, weight ) 

def fillHistogramTime( pyl, ntuple, histo, weight ):
    for event in range(ntuple.GetEntries()):
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
    mcHisto = ROOT.TH3D("mc", "mc", 20, -1., 1., 20, -pi, pi, 20, -1., 1.)
    toyHisto = ROOT.TH3D("toy", "toy", 20, -1., 1., 20, -pi, pi, 20, -1., 1.)
    mcHisto.Sumw2()
    toyHisto.Sumw2()
    mcHistoTime = ROOT.TH1D("mcTime", "mcTime", 20, -2., 18.)
    toyHistoTime = ROOT.TH1D("toyTime", "toyTime", 20, -2., 18.)
    mcHistoTime.Sumw2()
    toyHistoTime.Sumw2()

    numMCEvents = mcNtuple.GetEntries()
    numToyEvents = toyNtuple.GetEntries()
 
    #toyWeight = float(numMCEvents)/numToyEvents

    #print "There are %d MC events and %d toy events. Weight will be %f" % (numMCEvents, numToyEvents, toyWeight)

    fillHistogram( mcData, mcNtuple, mcHisto, 1. )
    fillHistogram( toyData, toyNtuple, toyHisto, 1. )
    
    fillHistogramTime( mcData, mcNtuple, mcHistoTime, 1. )
    fillHistogramTime( toyData, toyNtuple, toyHistoTime, 1. )
    
    #mcHisto.Scale(1./mcHisto.Integral())
    #toyHisto.Scale(1./toyHisto.Integral())

    # Create the 4D histogram for the acceptance surface
    acceptance = mcHisto.Clone()
    acceptance.SetTitle("acceptance")
    acceptance.SetName("acceptance")
    acceptance.Divide(toyHisto)
    acceptance.Scale(1./acceptance.Integral())
    
    acceptanceTime = mcHistoTime.Clone()
    acceptanceTime.SetTitle("acceptanceTime")
    acceptanceTime.SetName("acceptanceTime")
    acceptanceTime.Divide(toyHistoTime)
    acceptanceTime.Scale(1./acceptanceTime.Integral())

    # Let's see the 1D projections for each of the angles
    accCosTheta = acceptance.Project3D("x")
    accCosTheta.SetTitle("cosTheta")
    accCosTheta.SetName("cosTheta")
    accPhi = acceptance.Project3D("y")
    accPhi.SetTitle("phi")
    accPhi.SetName("phi")
    accCosPsi = acceptance.Project3D("z")
    accCosPsi.SetTitle("cosPsi")
    accCosPsi.SetName("cosPsi")

    # Now let's see the 2D projections
    accCosThetaPhi = acceptance.Project3D("xy")
    accCosThetaPhi.SetTitle("phi vs cosTheta")
    accCosThetaPhi.SetName("phi_cosTheta")
    accPhiCosPsi = acceptance.Project3D("yz")
    accPhiCosPsi.SetTitle("cosPsi vs phi")
    accPhiCosPsi.SetName("cosPsi_phi")
    accCosPsiCosTheta = acceptance.Project3D("zx")
    accCosPsiCosTheta.SetTitle("cosTheta vs cosPsi")
    accCosPsiCosTheta.SetName("cosTheta_cosPsi")
 
    mcFile.Close()
    toyFile.Close()
    acceptanceFile.Write()

if __name__ == "__main__":
	class PyListOfLeaves(dict) :
        	pass
	main()
