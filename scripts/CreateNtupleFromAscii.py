#! /usr/env python
import ROOT
from optparse import OptionParser
from array import array

__author__ = 'Greig A Cowan'
__date__ = '1/11/2009'
__version__ = '$Id: CreateNtupleFromAscii.py,v 1.1 2009/11/01 22:59:59 gcowan Exp $'

'''
Takes an ASCII file containing data with headers at the top of the file and
creates the corresponding nTuple using the same header names for the branches.
Currently is just creates everything as floats, so I'm not sure of the
behaviour if we are trying to compare to ints...
'''

def main():
    parser = OptionParser(
        usage = 'usage: %prog ascii_file.txt ntuple.root')

    (options, args) = parser.parse_args()
    inputName = args[0]
    outputName = args[1]

    inputFile = open(inputName)
    outputFile = ROOT.TFile( outputName, "RECREATE")  
 
    data = inputFile.readlines()
    readFirstLine = False

    i = 0
    for datum in data:
	observables = datum.rstrip("\n").split()
	if not readFirstLine:
		string = ":".join(observables[1:])
		ntuple = ROOT.TNtuple("ntuple","data from ascii file", string);
		readFirstLine = True
		continue
	doubleList = []
	for item in observables:
		doubleList.append(float(item))
        event = array('f', doubleList)
	ntuple.Fill(event)	
	if i >= 9999: 
		break
	else:
		i += 1		

    inputFile.close()
    outputFile.Write()


if __name__ == "__main__":
	main()
