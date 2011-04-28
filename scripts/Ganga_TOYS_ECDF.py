def repeatjobs( number, per_core, XML ):
	args = []
	for i in range(0,number,per_core):
		temp = []
		temp.append("-f")
		temp.append(str( XML ))
		temp.append("-repeats")
		temp.append(str(per_core))
		temp.append("--doPulls")
		temp.append("PullPlots.root")
		args.append(temp)
	return args

INPUT_XML_PATH = "./"
INPUT_XML_FILE = "TOYS.xml"

j = Job (
name = 'TOYS' ,
 outputsandbox = ['PullPlots.root'] ,
 info = JobInfo (
    ) ,
 inputdata = None ,
 merger = RootMerger (
    files = ['PullPlots.root'] ,
    args = None ,
    ignorefailed = False ,
    overwrite = False 
    ) ,
 inputsandbox = [ INPUT_XML_PATH+INPUT_XML_FILE ] ,
 application = Executable (
    exe = '/exports/work/physics_ifp_ppe/sw/fitting-latest/RapidFit-latest/bin/fitting',
    env = {} ,
    args = [] 
    ) ,
 outputdata = None ,
 splitter = ArgSplitter( args = repeatjobs( 1000, 100, INPUT_XML_FILE ) ) ,
 backend = SGE()
    ) 
j.submit()
