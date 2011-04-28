from math import *

#	AUTOMATICALLY DERRIVE ALL JOB CONFIG BASED ON USER INPUT
def fcsplitter(XML='XML.xml',toys=50,toyspersj=50,par1='x',par1min=0.,par1max=0.,par1res=1,par2='y',par2min=0.,par2max=0.,par2res=1,SEED=12345):
#	ALL ARGUMENTS TO PASS BACK TO GANGA
        args = []
#	STEP SIZE
        step1 = abs(par1min-par1max)/float(par1res)
        step2 = abs(par2min-par2max)/float(par2res)
#	ALL FC
	for k in range(0,toys,toyspersj):
#	One FC layer of toysperj per grid-point
        	for i in range(0,par1res+1,1):
                	par1val = par1min + float(i)*step1
                	for j in range(0,par2res+1,1):
				SEED=SEED+1
                        	par2val = par2min + float(j)*step2
				sjarg=[]
#	RapidFit_FC.C takes arguments in the following Order
				XMLstr = str( XML )
				paramstr1 = str(  par1+','+str(par1val)+','+str(par1val)+',1' )
				paramstr2 = str(  par2+','+str(par2val)+','+str(par2val)+',1' )
				events_str = str( toyspersj )
				seed_str = str( SEED )
#	Put Arguments in a list
				sjarg.append(str(XML))
				sjarg.append(str(paramstr1))
				sjarg.append(str(paramstr2))
				sjarg.append(str(events_str))
				sjarg.append(str(seed_str))
				args.append(sjarg)
#	OUTPUT TO CHECK BY USER
	print args
#	RETURN
	return args



#	WANT A ROOT JOB, 5.28.00 broken on Grid
j = Job( application = Root( version = '5.26.00b' ) )

#	Name of Jobs
j.name = 'RapidFit FC Edinburgh'

#	Sctipt to Wrap RapidFit in
j.application.script = File( name='./RapidFit_FC.C' )

#	Input Required by jobs
j.inputsandbox = ['./RapidFit_FC.C','../config/betas_tagged_analysis/tagged_cfit_3_free_ECDF.xml','../lib/libRapidRun.so']
#	Output given from FC Scan
j.outputsandbox = ['FCOutput.root']

#	ROOT MERGE OUTPUT broken for large number of files but leave in anyway incase root fixes it
j.merger = RootMerger( files=['FCOutput.root'] )

#	INTERACTIVE FOR TESTING
#j.backend=Interactive()

j.backend=SGE()

#	THE CONFIGURATION OF ALL JOBS WITHIN THE SCAN
j.splitter=ArgSplitter( args = fcsplitter('tagged_cfit_3_free.xml',10,10,'Phi_s',0,2*pi,39,'deltaGamma',-0.7,0.7,39,67890) )

#	RUN
j.submit()

