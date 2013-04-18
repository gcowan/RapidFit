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

job_name = "untagged_Bc"

#	All the possible output files right now
#
#			THIS HAS TO BE CHANGED BASED ON THE OUTPUT FROM YOUR SCAN

output_file_list = [ '*.root', '*/*.root' ]


LFN_LIST=[]
FILE_LIST=[]


xml = str()
script_name = str()

#	written up here for clarity, process all possible LFNs and PFNs that have to be processed
if is_ganga:
	for arg in sys.argv:
		if string.find( arg, "LFN:" ) != -1 :
			LFN_LIST.append( str( arg ) )
		elif string.find( arg, ".xml" ) != -1 :
			xml = str( arg )
		elif string.find( arg, ".py" ) != -1 :
			script_name = str( arg )
		else:
			FILE_LIST.append( str( arg ) )
	if LFN_LIST:
		print "LFNs:"
		print LFN_LIST
	print "FILEs:"
	print FILE_LIST

#	Splitter for toy studies
def Arg_Splitter( XML='XML.xml', RapidFit_Path="" ):
	args = []

	temp = []
	temp.append( str( XML ) )
	temp.append( RapidFit_Path )
        args.append( temp )
	print args
	return args

#	GANGA JOB

#	This is the section of code which will be executed within ganga
if is_ganga:

	ROOT_VERSION = str( os.popen("root-config --version | sed -e \'s/\\//\./g' ").readline() )[:-1]
	RapidFit_Path = os.environ.get("RAPIDFITROOT")
	if not RapidFit_Path:
		print ""
		print "Couldn't find $RAPIDFITROOT, Please check this!!"
		print ""
		sys.exit(-3)

	RapidFit_Library=RapidFit_Path+"/lib/libRapidRun.so"

	if not os.path.isfile( RapidFit_Library ):
		print "Please (re) compile RapidFit for submitting to a batch system"
		print "              run:    'make lib'           not just:   'make'"
		print ""
		print "This could also mean you haven't defined '$RAPIDFITROOT, please check!!"
		print ""
		sys.exit(-42)

	#	By definition of how this should be run!
	script_name = str( sys.argv[0] )

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


	#	Backend to submit jobs to
	#	This changes based on the system your on

	host_name = os.popen('hostname').readline()

	choice = -1

	if ( string.find( host_name, "frontend" ) != -1 ):
		print "Running on ECDF, submitting to SGE"
		j.backend = SGE()
		j.outputsandbox = output_file_list
		sandbox_data = [ script_name, xml, RapidFit_Library ]
		for k in FILE_LIST:
			sandbox_data.append( k )
		j.inputsandbox = sandbox_data

	elif ( string.find( host_name, "lxplus" ) != -1 ):

		if LFN_LIST:
			choice = 1
		else:
			choice = int( raw_input("Running on LXPLUS, submit to 1) GRID 2) lxplus Batch or 3) Interactive?\t") )
			while ( choice != 1 ) and ( choice != 2 ) and ( choice != 3 ):
				choice = int( raw_input( "try again...  " ) )

		if choice == 1:
			j.backend = Dirac()
			#	Only do this if the LFN_LIST is NOT empty
			if LFN_LIST:
				j.inputdata = LFN_LIST                  #       Point the job to the data
				j.backend.inputSandboxLFNs = LFN_LIST   #       Tell Dirac we need a local copy in order to process it
			sandbox_data = [ script_name, xml, RapidFit_Library ]
			#print sandbox_data
			for k in FILE_LIST:
				sandbox_data.append( k )
			#print sandbox_data
			j.inputsandbox = sandbox_data
			#	Now only works with ganga 6
			j.outputfiles = output_file_list
			#j.outputdata = output_file_list
		if choice == 2:
			j.backend = LSF()
			j.backend.queue = '8nh'         #       1nh, 8nh, 1nd, 2nd, 1nw, 2nw
			j.inputdata = FILE_LIST
			#j.outputsandbox = j.outputdata.files
			#j.outputdata=[]
		if choice == 3:
			j.backend = Interactive()
			j.inputsandbox = FILE_LIST

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
		j.outputdata = output_file_list

	elif ( string.find( host_name, "epfl" ) != -1 ):
		choice = int( raw_input("Running on EPFL, submit to 1) PBS or 2) Interactive? ") )
		while ( choice != 1 ) and ( choice != 2 ):
			choice = int( raw_input( "try again... " ) )
		if choice == 1 :
			j.backend = PBS()
			j.backend.queue = 'long'
			j.inputdata = FILE_LIST
		if choice == 2 :
			j.backend = Interactive()
			j.inputdata = FILE_LIST
		j.outputdata = output_file_list

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
	if choice == 3 :
		j.splitter=ArgSplitter( args = Arg_Splitter( xml, RapidFit_Path ) )
	else :
		j.splitter=ArgSplitter( args = Arg_Splitter( FIT_XML, RapidFit_Path ) )
	#	submit the job
	j.submit()



#	Actual pyROOT code which will be executed
if ( __name__ == '__main__' ) and ( not is_ganga ) :

	#	Just to record the input for any debugging
	for i in sys.argv:
		print i

	sys.path.append( sys.argv[2]+"/lib" )
	sys.path.append("/exports/work/ppe/sw/fitting-latest/root-5.29/lib")

	print sys.path

	#	We want ROOT
	import ROOT

	ROOT.gSystem.AddDynamicPath( sys.argv[2]+"/lib" )

	#	Input Parameters
	FIT_XML = sys.argv[1]

	#	Load the RapidFit binary library
	ROOT.gSystem.Load("libRapidRun")


	#	RapidFit arguments
	args = ROOT.TList()
	#	Construct the RapidFit Arguments as you would when running the fitter binary
	args.Add( ROOT.TObjString( "RapidFit"     ) )
	args.Add( ROOT.TObjString( "-f"           ) )
	args.Add( ROOT.TObjString( str( FIT_XML ) ) )
	args.Add( ROOT.TObjString( "-repeats" ))
	args.Add( ROOT.TObjString( "10" ))
	#	Print the command that is being run for reference
	#print args

	#	Construct an instance of the Fitting algorithm
	fitter = ROOT.RapidRun( args )
	#	Run the Fit
	result = fitter.run()

	#	Exit
	ROOT.gApplication.Terminate()

