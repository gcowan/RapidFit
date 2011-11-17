#! /usr/env python
import ROOT
from optparse import OptionParser
from array import array
import math

__author__ = 'Greig A Cowan'
__date__ = '1/11/2009'
__version__ = '$Id: CreateNtupleWithCorrectObservableNames.py,v 1.1 2009/11/01 22:59:59 gcowan Exp $'

'''
Takes an ntuple as input and produces an ntuple as output where the output has
the correct branch names such that they can be used as input to RapidFit.
Obviously, the second ntuple will only contain the branches that you specify
in this script.
'''

def main():
    parser = OptionParser(
        usage = 'usage: %prog ntuple.root ntuple/path')

    (options, args) = parser.parse_args()
    inputFileName = args[0]
    pathToNtuple = args[1]

    inputFile = ROOT.TFile(inputFileName)
    name = inputFileName.split(".")
    outputFileName = name[0] + "_RapidFit." + name[1]
    outputFile = ROOT.TFile(outputFileName, "RECREATE")  

    inputNtuple = inputFile.Get(pathToNtuple)

    try:
	leaves = inputNtuple.GetListOfLeaves()
    except Exception, e:
	print e

    pyl = PyListOfLeaves()
    for i in range(leaves.GetEntries()):
       	leaf = leaves.At(i)
       	name = leaf.GetName()
       	pyl.__setattr__(name,leaf)

    # This list MUST be in the correct order, as defined by the 
    # valueList lines in the loop below. It's not pretty, but it works.
    newBranches = [ "time"
		  #, "cosTheta"
		  #, "cosPsi"
		  #, "phi"
		  , "tag"
		  , "mistag"
		  , "mass"
		  , "truetime"
		  , "residual"
		  ]

    outputNtupleStructure = ":".join(newBranches) 
    outputNtuple = ROOT.TNtuple("ntuple","RapidFit", outputNtupleStructure);

    for event in range(inputNtuple.GetEntries()):
	inputNtuple.GetEntry(event)
    	valueList = []

	# ThetaTr, ThetaK, ThetaVtr are the names that the P2VVAngleCalculator
	# in DaVinci gives to the angles theta, psi and phi respectively.	
	'''
  	valueList.append( pyl.B_s_TAU.GetValue() )
	valueList.append( math.cos( pyl.B_s_ThetaTr.GetValue() ))
	valueList.append( math.cos( pyl.B_s_ThetaK.GetValue() ))
	valueList.append( pyl.B_s_ThetaVtr.GetValue() )
	valueList.append( pyl.B_s_TAGDECISION.GetValue() )
	valueList.append( pyl.B_s_TAGOMEGA.GetValue() )
	'''
  	valueList.append( pyl.B_s_TAU.GetValue()*1000 )
	valueList.append( pyl.B_s_TAGDECISION.GetValue() )
	valueList.append( pyl.B_s_TAGOMEGA.GetValue() )
	valueList.append( pyl.B_s_MM.GetValue() )
	valueList.append( pyl.B_s_TRUETAU.GetValue()*1000 )
	valueList.append( pyl.B_s_TAU.GetValue()*1000 - pyl.B_s_TRUETAU.GetValue()*1000 )
  
  	if pyl.B_s_TRUETAU.GetValue() < -0.2: continue
	data = array("f", valueList)
	outputNtuple.Fill(data)	

    inputFile.Close()
    outputFile.Write()
    outputFile.Close()

if __name__ == "__main__":
	class PyListOfLeaves(dict) :
        	pass
	main()
