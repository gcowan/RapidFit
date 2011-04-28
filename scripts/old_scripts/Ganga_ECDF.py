from math import *

def fcsplitter(XML='XML.xml',par1='x',par1min=0.,par1max=0.,par1res=1,par2='y',par2min=0.,par2max=0.,par2res=1):
        args = []
        step1 = abs(par1min-par1max)/float(par1res)
        step2 = abs(par2min-par2max)/float(par2res)
        for i in range(0,par1res+1,1):
                par1val = par1min + float(i)*step1
                for j in range(0,par2res+1,1):
                        par2val = par2min + float(j)*step2
                        #sjarg = ['-f','./XML_UB_B.xml','--doLLcontour','--defineContour']
			sjarg=[]
			XMLStr = str( XML )
			paramstr1 = str(  par1+','+str(par1val)+','+str(par1val)+',1' )
			paramstr2 = str(  par2+','+str(par2val)+','+str(par2val)+',1' )
			sjarg.append(str(XMLStr))
			sjarg.append(str(paramstr1))
			sjarg.append(str(paramstr2))
			#print sjarg
			args.append(sjarg)
	print args
	return args

j = Job( application = Root() )
j.name = 'rootrun'
j.application.script = File( name='./RapidFit.C' )
j.inputsandbox = ['./RapidFit.C','../config/betas_tagged_analysis/tagged_cfit_3_free_ECDF.xml','../lib/libRapidRun.so']
j.outputsandbox = ['LLcontourScanData.root']
j.merger = RootMerger( files=['LLcontourScanData.root'] )
j.backend=SGE()
j.splitter=ArgSplitter( args = fcsplitter('tagged_cfit_3_free_ECDF.xml','Phi_s',0,2*pi,39,'deltaGamma',-0.7,0.7,39) )
j.submit()

