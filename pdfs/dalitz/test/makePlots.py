#!/usr/bin/python

from ROOT import TH1F, TTree, TFile, TCanvas, TLegend, TLatex
from ROOT import gDirectory, THStack, gROOT
import ROOT
import sys

myfile1 = TFile( '../AmplitudeUnderstanding/test.root')
myfile2 = TFile( './test/testJpsiK892.root')

treeEvt = myfile1.Get("tree")
treeMyAmp = myfile2.Get("tree")

c = TCanvas()
c.Draw()

pdfFile = sys.argv[1]

c.Print(pdfFile+'[')

treeEvt.Draw("mKpi>>hh1")
hh1=gROOT.FindObject("hh1")
hh2=hh1.Clone('hh2')
hh2.Sumw2()
treeMyAmp.Draw("mKpi>>hh2","w")

c.Clear()
hh2.SetLineColor(2)
hh2.SetMarkerColor(2)
hh2.SetMarkerSize(0.05)
hh1.DrawNormalized()
hh2.DrawNormalized('same')
c.Print(pdfFile)

treeEvt.Draw("mJpsipi>>hh3")
hh3=gROOT.FindObject("hh3")
hh4=hh3.Clone('hh4')
hh4.Sumw2()
treeMyAmp.Draw("mJpsipi>>hh4","w")

c.Clear()
hh4.SetLineColor(2)
hh4.SetMarkerColor(2)
hh4.SetMarkerSize(0.05)
hh3.DrawNormalized()
hh4.DrawNormalized('same')
c.Print(pdfFile)

treeEvt.Draw("mJpsiK>>hh5")
hh5=gROOT.FindObject("hh5")
hh6=hh5.Clone('hh6')
hh6.Sumw2()
treeMyAmp.Draw("mJpsiK>>hh6","w")

c.Clear()
hh6.SetLineColor(2)
hh6.SetMarkerColor(2)
hh6.SetMarkerSize(0.05)
hh5.DrawNormalized()
hh6.DrawNormalized('same')
c.Print(pdfFile)

treeEvt.Draw("theta1>>hh7")
hh7=gROOT.FindObject("hh7")
hh8=hh7.Clone('hh8')
hh8.Sumw2()
treeMyAmp.Draw("cosTheta1>>hh8","w")

c.Clear()
hh8.SetLineColor(2)
hh8.SetMarkerColor(2)
hh8.SetMarkerSize(0.05)
hh7.DrawNormalized()
hh8.DrawNormalized('same')
c.Print(pdfFile)

treeEvt.Draw("theta2>>hh7")
treeMyAmp.Draw("cosTheta2>>hh8","w")

c.Clear()
hh7.DrawNormalized()
hh8.DrawNormalized('same')
c.Print(pdfFile)

treeEvt.Draw("phi>>hh9")
hh9=gROOT.FindObject("hh9")
hh10=hh9.Clone('hh10')
hh10.Sumw2()
treeMyAmp.Draw("phi>>hh10","w")

c.Clear()
hh10.SetLineColor(2)
hh10.SetMarkerColor(2)
hh10.SetMarkerSize(0.05)
hh9.DrawNormalized()
hh10.DrawNormalized('same')
c.Print(pdfFile)


c.Print(pdfFile+']')
