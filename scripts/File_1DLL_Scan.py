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


#	USAGE:
#
#	ganga script_name.py xml_for_job.xml file1.root file2.root	to run on most backends
#
#	ganga script_name.py xml_for_job.xml LFN1  LFN2 		to run on the GRID
#

#	Configurables

job_name = "TOY-1D"

#	Semi-Frequently changed script veriables:
TOTAL_STEPS = 2					#	Total unique subsets in file
MC_FILE_STEP = 10000				#	Size of Step to take in file
MC_FIT_STEP = 5000				#	Step size between subset selection

#	All the possible output files right now
output_file_list = [ 'MC_Study.root' ]


RapidFit_Path = "/afs/cern.ch/user/r/rcurrie/RapidFit_Source/"
RapidFit_Library=RapidFit_Path+"lib/libRapidRun.so"

ROOT_VERSION='5.26.00b'

LFN_LIST=[]
FILE_LIST=[]



#	written up here for clarity, process all possible LFNs and PFNs that have to be processed
if is_ganga:
	i=0
	for arg in sys.argv:
		if i > 1:
			if string.find( arg, "LFN:" ) != -1 :
				LFN_LIST.append( str( arg ) )
			else:
				FILE_LIST.append( str( argv ) )
		i=i+1
	print "LFNs:"
	print LFN_LIST
	print "FILEs:"
	print FILE_LIST

#	Splitter for toy studies
def MCStudy_Splitter( XML='XML.xml', ALL_Steps=30, MC_FILE_Step = 10000, MC_FIT_Step = 5000 ):
	args = []
	steps_per_jump = int( MC_FILE_Step / MC_FIT_Step )

	for i in range( 0, ALL_Steps, 1 ):
		temp = []
		temp.append( str( XML ) )
		temp.append( str( MC_FILE_Step*i ) )
		temp.append( str( MC_FIT_Step ) )
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

	#
	j.application.script = File( name=script_name )
	#	Tell the script where the RapidFit library is
	j.inputsandbox = [ script_name, xml, RapidFit_Library ]
	#	Outputs Wanted
	j.outputdata = output_file_list


	#	Backend to submit jobs to
	#	This changes based on the system your on

	host_name = os.popen('hostname').readline()

	if ( host_name == "frontend01"  ) or ( host_name == "frontend02" ):
		print "Running on ECDF, submitting to SGE"
		j.backend = SGE()

	elif ( string.find( host_name, "lxplus" ) != -1 ):
		choice = int( raw_input("Running on LXPLUS, submit to 1) GRID 2) lxplus Batch or 3) Interactive?\t") )
		while ( choice != 1 ) and ( choice != 2 ):
			choice = int( raw_input( "try again...  " ) )
		if choice == 1:
			j.backend = Dirac()
			j.inputdata = LFN_LIST			#	Point the job to the data
			j.backend.inputSandboxLFNs = LFN_LIST	#	Tell Dirac we need a local copy in order to process it
			if choice == 2:
				j.backend = LSF()
				j.backend.queue = '1nh'		#	1nh, 8nh, 1nd, 2nd, 1nw, 2nw
				j.inputdata = FILE_LIST
			if choice == 3:
				j.backend = Interactive()
				j.inputdata = FILE_LIST

	elif ( string.find( host_name, "ppe" ) != -1 ):
		choice = int( raw_input("Running on PPE, submit to 1) CONDOR or 2) Interactive? ") )
		while ( choice != 1 ) and ( choice != 2 ):
			choice = int( raw_input( "try again... " ) )
		if choice == 1 :
			j.backend = Condor()
			j.inputdata = FILE_LIST
		if choice == 2 :
			j.backend = Interactive()
			j.inputdata = FILE_LIST

	else:
		print "Unknown system, just running the job, check this is what you want"
		j.inputdata = FILE_LIST

        #       Input Parameters
        FIT_XML = xml
        FIT_LIST = string.split( FIT_XML, "/" )
        if len( FIT_LIST ) == 1:
                FIT_LIST = string.split( FIT_XML, "\\" )

        #       just need the absolute name of the XML in order to run on the backend
        FIT_XML = FIT_LIST[ int(len(FIT_LIST)-1) ]

	#	Splitter to use for job
	j.splitter=ArgSplitter( args = MCStudy_Splitter( FIT_XML, TOTAL_STEPS, MC_FILE_STEP, MC_FIT_STEP ) )
	#	submit the job
	j.submit()



#	Actual pyROOT code which will be executed
if ( __name__ == '__main__' ) and ( not is_ganga ) :

	#	Just to record the input for any debugging
	for i in sys.argv:
		print i

	#	We want ROOT
	import ROOT

	#	Input Parameters
	FIT_XML = sys.argv[1]

	MC_START = sys.argv[2]
	MC_STEP = sys.argv[3]

	#	Load the RapidFit binary library
	ROOT.gSystem.Load("libRapidRun")


	#	RapidFit arguments
	args = ROOT.TList()
	#	Construct the RapidFit Arguments as you would when running the fitter binary
	args.Add( ROOT.TObjString( "RapidFit"     ) )
	args.Add( ROOT.TObjString( "-f"           ) )
	args.Add( ROOT.TObjString( str( FIT_XML ) ) )
	args.Add( ROOT.TObjString( "--MCStudy"    ) )
	args.Add( ROOT.TObjString("--MCStartEntry") )
	args.Add( ROOT.TObjString( str( MC_START )) )
	args.Add( ROOT.TObjString( "--MCStepSize" ) )
	args.Add( ROOT.TObjString( str( MC_STEP ) ) )

	#	Print the command that is being run for reference
	#print args

	#	Construct an instance of the Fitting algorithm
	fitter = ROOT.RapidRun( args )
	#	Run the Fit
	result = fitter.run()

	#	Exit
	ROOT.gApplication.Terminate()

