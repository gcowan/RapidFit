#	For any constants and such, will more than likely be used as we're doing complex stuff
from math import *
#	Interfacing with the command line arguments is required
import sys, os
#       Very useful for string manipulations
import string
#	Used for tagging the date/time into the job desc
import datetime
now = datetime.datetime.now()

#	This script is intended to be run as 'ganga script.py', it will submit itself and the relative configuration as a pyROOT job
is_ganga = "Ganga.Core" in sys.modules



#	Configurables

job_name = "BIAS-TOY"

#	Starting Seed and number of toys we want in total and how many per core
SEED = 12345
toys_in_study = 5000
toys_per_job = 250

USE_UUID="true"

#	All the possible output files right now
output_file_list = [ 'PullPlots.root' ]


RapidFit_Path = "/afs/cern.ch/user/r/rcurrie/RapidFit_Source/"
RapidFit_Library=RapidFit_Path+"lib/libRapidRun.so"

ROOT_VERSION='5.26.00b'

#	Splitter for toy studies
def TOY_Splitter( XML='XML.xml', SEED='12345', total_number_of_toys=20, number_of_toys_per_job=5):
	args = []
	real_SEED = SEED
	for i in range(0,total_number_of_toys,number_of_toys_per_job):
		temp = []
		temp.append( str( XML ) )
		temp.append( str( number_of_toys_per_job ) )
		temp.append( str( real_SEED ) )
		real_SEED = real_SEED + 1
                args.append( temp )
	print args
	return args




#	GANGA JOB

#	This is the section of code which will be executed within ganga
if is_ganga:

	#	By definition of how this should be run!
	script_name = str( sys.argv[0] )
	xml = str( sys.argv[1] )

	#	i.e.	> ganga script.py some.xml

        #       Input Parameters
        script_onlyname = script_name
        script_list = string.split( script_name, "/" )
        if len( script_list ) == 1:
                script_onlyname = string.split( script_name, "\\" )
        script_onlyname = script_list[ int(len(script_list)-1) ]


	#	Create the job
	j = Job( application = Root( version = ROOT_VERSION ) )

	datetimeinfo = str( "_" + str( now.strftime("%Y-%m-%d_%H.%M") ) )
	#       Change the name of your job for records
	j.name = str(job_name + "_" + str(script_onlyname) + datetimeinfo)


	print script_name
	j.application.script = File( name=script_name )
	#	Tell the script where the RapidFit library is
	j.inputsandbox = [ script_name, xml, RapidFit_Library ]
	#	Outputs Wanted
	j.outputsandbox = output_file_list


	#	Backend to submit jobs to
	#	This changes based on the system your on

	host_name = os.popen('hostname').readline()

	if ( host_name == "frontend01"  ) or ( host_name == "frontend02" ):
		print "Running on ECDF, submitting to SGE"
		j.backend = SGE()

	elif ( string.find( host_name, "lxplus" ) != -1 ):
		choice = int( raw_input("Running on LXPLUS, submit to 1) GRID 2) lxplus Batch or 3) Interactive?\t") )
		while ( choice != 1 ) and ( choice != 2 ) and ( choice != 3 ):
			choice = int( raw_input( "try again...  " ) )
		if choice == 1:
			j.backend = Dirac()
		if choice == 2:
			j.backend = LSF()
			j.backend.queue = '8nh'		#	1nh, 8nh, 1nd, 2nd, 1nw, 2nw
		if choice == 3:
			j.backend = Interactive()

	elif ( string.find( host_name, "ppe" ) != -1 ):
		choice = int( raw_input("Running on PPE, submit to 1) CONDOR or 2) Interactive? ") )
		while ( choice != 1 ) and ( choice != 2 ):
			choice = int( raw_input( "try again... " ) )
		if choice == 1 :
			j.backend = Condor()
		if choice == 2 :
			j.backend = Interactive()

	else:
		print "Unknown system, just running the job, check this is what you want"


        #       Input Parameters
        FIT_XML = xml
        FIT_LIST = string.split( FIT_XML, "/" )
        if len( FIT_LIST ) == 1:
                FIT_LIST = string.split( FIT_XML, "\\" )

        #       just need the absolute name of the XML in order to run on the backend
        FIT_XML = FIT_LIST[ int(len(FIT_LIST)-1) ]


	#	Splitter to use for job
	j.splitter=ArgSplitter( args = TOY_Splitter( FIT_XML, SEED, toys_in_study, toys_per_job ) )
	#	submit the job
	j.submit()




#	RapidFit JOB

#	Actual pyROOT code which will be executed
if ( __name__ == '__main__' ) and ( not is_ganga ) :

	#	Just to record the input for any debugging
	for i in sys.argv:
		print i

	#	We want ROOT
	import ROOT

	FIT_XML = sys.argv[1]

	num_repeats = sys.argv[2]

	SEED = sys.argv[3]


	#	Load the RapidFit binary library
	ROOT.gSystem.Load("libRapidRun")

	#	RapidFit arguments
	args = ROOT.TList()
	#	Construct the RapidFit Arguments as you would when running the fitter binary
	args.Add( ROOT.TObjString( "fitting"      ) )
	args.Add( ROOT.TObjString( "-f"           ) )
	args.Add( ROOT.TObjString( str( FIT_XML ) ) )
	args.Add( ROOT.TObjString( "-repeats"     ) )
	args.Add( ROOT.TObjString(str(num_repeats)) )
	args.Add( ROOT.TObjString( "--doPulls"    ) )
	args.Add( ROOT.TObjString("PullPlots.root") )
	args.Add( ROOT.TObjString( "--SetSeed"    ) )
	args.Add( ROOT.TObjString( str( SEED )    ) )
	if USE_UUID:
		args.Add( ROOT.TObjString( "--useUUID" ) )

	#	Construct an instance of the Fitting algorithm
	fitter = ROOT.RapidRun( args )
	#	Run the Fit
	result = fitter.run()

	#	Exit
	ROOT.gApplication.Terminate()

