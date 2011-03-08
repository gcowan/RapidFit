from math import *

#def fcsplitter(commonargs=['-f','XML_File.xml','--doFCscan'],toys=10,toyspersj=5,par1='x',par1min=0.,par1max=0.,par1res=1,par2='y',par2min=0.,par2max=0.,par2res=1):
def fcsplitter(commonargs='-f XML_File.xml --doFCscan',toys=10,toyspersj=5,par1='x',par1min=0.,par1max=0.,par1res=1,par2='y',par2min=0.,par2max=0.,par2res=1,seed=0):
	nseed = seed
        args = []
        step1 = abs(par1min-par1max)/float(par1res)
        step2 = abs(par2min-par2max)/float(par2res)
	for k in range(0,toys,toyspersj):
	        for i in range(0,par1res+1,1):
                	par1val = par1min + float(i)*step1
	                for j in range(0,par2res+1,1):
        	                par2val = par2min + float(j)*step2
				nseed = nseed + 1
                        	sjarg = []
				sjarg.append(commonargs+' -repeats '+str(toyspersj)+' --defineContour '+ par1+','+str(par1val)+','+str(par1val)+',1 '+par2+','+str(par2val)+','+str(par2val)+',1'+' --useUUID --SetSeed '+str(nseed))
				#sjarg.extend(commonargs)
				#sjarg.extend(['-repeats',str(toyspersj),'--defineContour',par1+','+str(par1val)+','+str(par1val)+',1',par2+','+str(par2val)+','+str(par2val)+',1'])
				#print sjarg
                        	args.append(sjarg)
        return args

j = Job (
name = 'Feldman_Cousins_test' ,
 outputsandbox = ['pullPlots.root','fit.log'] ,
 info = JobInfo (
    ) ,
 inputdata = None ,
 merger = RootMerger (
    files = ['FCOutput.root'] ,
    args = None ,
    ignorefailed = False ,
    overwrite = False 
    ) ,
 inputsandbox = [ ] ,
 application = Executable (
    exe = '/phys/linux/s0127440/sWave_scripts/runrapidfit.sh',
    env = {} ,
    args = [] 
    ) ,
 outputdata = None ,
 splitter = ArgSplitter (
 args = fcsplitter('-f /phys/linux/s0127440/sWave_scripts/pete-unbiased-upper.xml --doFCscan',20,100,'Phi_s',0.0,2.0*pi,24,'deltaGamma',-0.9,0.9,24,2)
    ) ,
 backend = Condor ( ) 
    ) 
j.submit()
