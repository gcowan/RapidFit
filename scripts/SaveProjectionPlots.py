import ROOT
from optparse import OptionParser

__author__ = 'Greig A Cowan'
__date__ = '7/10/2009'
__version__ = '$Id: PrintLatexTableOfPullResults.py,v 1.2 2009/10/07 11:49:17 gcowan Exp $'

'''
'''

def main():
    parser = OptionParser(
        usage = 'usage: %prog projections.root')

    (options, args) = parser.parse_args()

    ROOT.gROOT.Reset()
    ROOT.gROOT.SetStyle('Plain')
    ROOT.gStyle.SetTitleStyle(1)
    ROOT.gStyle.SetOptFit(1111)
    ROOT.gStyle.SetOptStat(1111)

    f = ROOT.TFile.Open(args[0])
    keys = f.GetListOfKeys()
    keyIterator = keys.MakeIterator()
    while True:
	key = keyIterator.Next()
	if not key: break
	name = key.GetName()
	title = key.GetTitle()
	if "Projection" in name:
		c = f.Get(name)
		c.SaveAs( name + ".pdf")

if __name__ == '__main__':
	main()
