run_arg_0_0=['-f','XML_File.xml','--doLLcontour','--doFCscan','-repeats','50','--defineContour','Phi_s,0.,0.5,10','deltaGamma,0.,0.5,10']
j_0_0=Job(application=Executable(exe='/exports/home/s0957039/trunk/bin/fitting',args=run_arg_0_0))
j_0_0.inputsandbox=[File('./XML_File.xml')]
j_0_0.outputsandbox=['FCOutput.root']
j_0_0.backend=SGE()
j_0_0.submit()
