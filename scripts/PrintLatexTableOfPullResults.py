import ROOT
import os
from xml.dom.minidom import parse
from optparse import OptionParser

__author__ = 'Greig A Cowan'
__date__ = '7/10/2009'
__version__ = '$Id: PrintLatexTableOfPullResults.py,v 1.2 2009/10/07 11:49:17 gcowan Exp $'

'''
Used to summarise the outputs from multiple RapidFit toy studies that have been 
run on a batch system. Will produce the final pull plots with fitted Gaussian curves
and a LaTeX table containing a summary of the parameter sensitivities and pull widths
and means. This can be pasted into your paper!
'''

def main():
    parser = OptionParser(
        usage = 'usage: %prog RapidFit_config.xml pullPlots.root rapidFit_output.root')

    (options, args) = parser.parse_args()

    ROOT.gROOT.Reset()
    ROOT.gROOT.SetStyle('Plain')
    ROOT.gStyle.SetTitleStyle(1)
    ROOT.gStyle.SetOptFit(1111)
    ROOT.gStyle.SetOptStat(1111)

    xmlFile  = args[0]
    rootFile = ROOT.TFile.Open( args[1])
    dir = 'RapidFitPlots/'
    dom = parse( xmlFile ) # parse an XML file by name

    physicsParameters = []
    for e in dom.getElementsByTagName('PhysicsParameter'):
	if e.getElementsByTagName('Type')[0].childNodes[0].data != 'Fixed':
		physicsParameters.append(e.getElementsByTagName('Name')[0].childNodes[0].data)
		
	
	trees = []
    histos = []
    resultsList = []

    dir = 'RapidFitPlots/'
    if not os.path.isdir(dir):
	    os.mkdir(dir)
    outputFile = ROOT.TFile.Open( dir+args[2], 'RECREATE')

    for i in range( len(physicsParameters) ):
	param = physicsParameters[i]
	tree = rootFile.Get( param )
	try:
		leaves = tree.GetListOfLeaves()
    	except Exception, e:
		continue
    	canvas = ROOT.TCanvas()
	pyl = PyListOfLeaves()

    	for i in range(leaves.GetEntries()):
        	leaf = leaves.At(i)
        	name = leaf.GetName()
        	pyl.__setattr__(name,leaf)

        tree.Draw('>>evtList' + param)
        evtList = ROOT.gDirectory.Get('evtList'+param)
        nev = evtList.GetN()
	histo_pull = ROOT.TH1F( param + 'Pull', '%s pull distribution' % param, 100, -4.0, 4.0)
	histo_error = ROOT.TH1F( param + 'Error', '%s error distribution' % param, 100, -4.0, 4.0)
	histo_pull.Sumw2()
	histo_error.Sumw2()
	for i in range(nev):
        	try:
            		tree.GetEntry(evtList.GetEntry(i))
            		pull = pyl.pull.GetValue()
            		error = pyl.error.GetValue()
        	except Exception, e:
            		print e
		histo_pull.Fill( pull )
		histo_error.Fill( error )

	histo_pull.Draw()
	gfit_pull = ROOT.TF1("GaussianP","gaus",0.,100.)
	gfit_error = ROOT.TF1("GaussianE","gaus",0.,100.)
	histo_error.Fit( gfit_error)
	histo_pull.Fit( gfit_pull )
	
	chisq = gfit_pull.GetChisquare()
	ndf = gfit_pull.GetNDF()
	#chisqdf = chisq/ndf

	pull_constant = gfit_pull.GetParameter(0)
	pull_econstant = gfit_pull.GetParError(0) 
	pull_mean = gfit_pull.GetParameter(1)
	pull_emean = gfit_pull.GetParError(1) 
	pull_sigma = gfit_pull.GetParameter(2)
	pull_esigma = gfit_pull.GetParError(2)
		
	error_mean = gfit_error.GetParameter(1)
	error_emean = gfit_error.GetParError(1)

	resultsList.append( [ param, error_mean, error_emean, pull_mean, pull_emean, pull_sigma, pull_esigma])


 	canvas.Update()
    	canvas.SaveAs(dir+'%s.png' % param)

	outputFile.Write()
    outputFile.Close()

    print '\\begin{center}'
    print '\\begin{tabular}{|c|c|c|c|}' 
    print '\hline'
    print 'Physics Parameter & Sensitivity (check I have used correct definition) & Pull mean & Pull sigma \\\\'
    for result in resultsList:
    	print "%s & %.3f \pm %.3f & %.3f \pm %.3f & %.3f \pm %.3f \\\\" % tuple(result)
    print '\hline'
    print '\end{tabular}'
    print '\end{center}'


if __name__ == '__main__':
	class PyListOfLeaves(dict) :
        	pass
	main()
