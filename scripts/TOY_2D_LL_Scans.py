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



job_name = "TOY-2D"

output_file_list = [ 'Global_Fit_Result.root', 'LLcontourScanData.root' ]

SEED = 12345
num_of_scans = 1

X_par="Phi_s"
Y_par="deltaGamma"

X_min=-1
X_max=5
Y_min=-0.3
Y_max=0.3

resolution=50
sqrt_jobs_per_core=5

RapidFit_Path = "afs/cern.ch/user/r/rcurrie/RapidFit_Source/"
RapidFit_Library=RapidFit_Path+"lib/libRapidRun.so"

ROOT_VERSION='5.26.00b'

#	Determine the grid points between limits at a given resolution
def grid_points( par1min=0., par1max=0., par1res=1, par2min=0., par2max=0., par2res=1 ):
	grid=[]
	step1 = abs(par1min-par1max)/float(par1res)
        step2 = abs(par2min-par2max)/float(par2res)
        for i in range(0,par1res+1,1):
                row=[]
                par1val = par1min + float(i)*step1
                for j in range(0,par2res+1,1):
                        par2val = par2min + float(j)*step2
                        point=[]
                        point.append( str(par1val) )
                        point.append( str(par2val) )
                        row.append(point)
                grid.append(row)
                #print row
        return grid

#	Find the ranges for the subjobs
def cell_maker( grid_points, cellsize=2 ):
        len_x = len(grid_points)-1
        len_y = len(grid_points[0])-1
        grid_step = float(grid_points[0][1][1])-float(grid_points[0][0][0])
        split_grid=[]
        for i in range(0,len(grid_points),cellsize):
                for j in range(0,len(grid_points[0]),cellsize):
                        min_i = i
                        min_j = j
                        max_i = i+cellsize-1
                        max_j = j+cellsize-1
                        points_x = cellsize
                        points_y = cellsize
                        if max_i > len_x :
                                points_x-= int((max_i-len_x)/grid_step)
                                max_i = len_x
                        if max_j > len_y :
                                points_y-= int((max_j-len_y)/grid_step)
                                max_j = len_y
                        grid_limits = []
                        grid_limits.append( [ grid_points[min_i][min_j][0],
						grid_points[max_i][max_j][0], str(points_x) ] )
                        grid_limits.append( [ grid_points[min_i][min_j][1],
						grid_points[max_i][max_j][1], str(points_y) ] )
                        split_grid.append( grid_limits )
        return split_grid

#	2DLL Splitter
def LL_2DSplitter(XML='XML.xml',par1='x',par1min=0.,par1max=0.,par1res=1,par2='y',par2min=0.,par2max=0.,par2res=1,cellsize=2,SEED,number_of_scans):
        args = []
        par1res -= 1
        par2res -= 1
        grid_array = grid_points( par1min, par1max, par1res, par2min, par2max, par2res )
        grid_cells = cell_maker( grid_array, cellsize )
	real_SEED = SEED
	for i in range( 0,number_of_scans,1 )
		for cell in grid_cells:
			par1str= par1+","+str(cell[0][0])+","+str(cell[0][1])+","+str(cell[0][2])
			par2str= par2+","+str(cell[1][0])+","+str(cell[1][1])+","+str(cell[1][2])
			temp_arg = [ XML, par1str, par2str, str(real_SEED) ]
			real_SEED = real_SEED+3
		args.append( temp_arg )
	print args
	return args


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
				j.backend.queue = '1nh'		#	1nh, 8nh, 1nd, 2nd, 1nw, 2nw
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
        print FIT_LIST

        #       just need the absolute name of the XML in order to run on the backend
        FIT_XML = FIT_LIST[ int(len(FIT_LIST)-1) ]

	#	Splitter to use for job
	j.splitter=ArgSplitter( args = LL_2DSplitter( FIT_XML, X_par,X_min,X_max,resolution, Y_par,Y_min,Y_max,resolution,
						      sqrt_jobs_per_core, SEED, num_of_scans ) )

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

	#	Input Parameters
	FIT_XML = sys.argv[1]
	SEED = sys.argv[4]

	x_cont = sys.argv[2]
	y_cont = sys.argv[3]


	#	Load the RapidFit binary library
	ROOT.gSystem.Load("libRapidRun")

	#	RapidFit arguments
	args = ROOT.TList()

	#	Construct the RapidFit Arguments as you would when running the fitter binary
	args.Add( ROOT.TObjString( "fitting"       ) )
	args.Add( ROOT.TObjString( "-f"            ) )
	args.Add( ROOT.TObjString( str( FIT_XML )  ) )
	args.Add( ROOT.TObjString( "--doLLcontour" ) )
	args.Add( ROOT.TObjString("--defineContour") )
	args.Add( ROOT.TObjString( str( x_cont )   ) )
	args.Add( ROOT.TObjString( str( y_cont )   ) )
	args.Add( ROOT.TObjString( "--BurnToROOT"  ) )
	args.Add( ROOT.TObjString( "--SetSeed"     ) )
	args.Add( ROOT.TObjString( str( SEED )    ) )

	#	Construct an instance of the Fitting algorithm
	fitter = ROOT.RapidRun( args )

	#	Run the Fit
	result = fitter.run()

	#	Exit
	ROOT.gApplication.Terminate()

