#!/usr/bin/env python
#	Required for interfaceing with the OS & FS
import os
#	Needed to run the actual fitter Run the 2 tests in parallel
import subprocess
#	Required for logging the output into a sensible output
import logging
#	Required for unique output location
import time
#	Required for pi...... mmmmm pie :P
from math import *

#	First before we go anywhere lets test that RapidFit is built
os.system( 'cd .. && make -j3 && cd scripts' )

INPUT_DATA = "../unittest/Xdata/goldstandard_biased_timecut3PELC.root  ../unittest/Xdata/goldstandard_unbiased_timecut3.root"

Rapid_Bin = "../bin/fitting"

cFit_XML = "../unittest/cFit/DataFile_Tagged_cFit_MOMP.xml"
cFit_TOY = "../unittest/cFit/DataFile_Tagged_cFit_MOMP_TOYS.xml"
sFit_XML = "../unittest/sFit/DataFile_Tagged_sFit_MOMP.xml"

Phi_s_str = 'Phi_s,-'+str(pi)+','+str(pi)+',40'
RapidFit_args_cFit = [ './fitting', '-f', 'DataFile_Tagged_cFit_MOMP.xml', '--doLLcontour', '--defineContour', Phi_s_str, 'deltaGamma,-0.7,0.7,40' ]
Phi_s_str = 'Phi_s,-'+str(3.2)+','+str(3.2)+',40'
RapidFit_args_sFit = [ './fitting', '-f', 'DataFile_Tagged_sFit_MOMP.xml', '--doLLcontour', '--defineContour', Phi_s_str, 'deltaGamma,-0.7,0.7,40' ]

total_toys = 1000
toys_per_job = 10
TOY_SEED = 12345
RapidFit_args_cFit_TOY = [ './fitting', '-f', 'DataFile_Tagged_cFit_MOMP.xml', '-repeats', str(toys_per_job), '--SetSEED' ]

print "Setting up Output Folders and moving around files"

#	This sort of thing is handled by Ganga in other scripts

Unit_Test_Dir = "/tmp/Unit_Test_" +  str(int(time.time()))
cFit_output = Unit_Test_Dir + "/cFit_output"
cFit_TOYS_output = Unit_Test_Dir + "/cFit_TOYS_output"
sFit_output = Unit_Test_Dir + "/sFit_output"

cFit_logfile = cFit_output + "/" + "UnitTest_cFit_output.log"
cFit_TOYS_logout = cFit_TOYS_output + "/" + "UnitTest_cFit_TOYS_output.log"
sFit_logfile = sFit_output + "/" + "UnitTest_sFit_output.log"

if not os.path.exists(Unit_Test_Dir):
	os.makedirs(sFit_output)

if not os.path.exists(cFit_output):
	os.makedirs(cFit_output)

if not os.path.exists(cFit_TOYS_output):
        os.makedirs(cFit_TOYS_output)

if not os.path.exists(sFit_output):
	os.makedirs(sFit_output)


os.popen( 'cp ' + INPUT_DATA + ' ' + cFit_output )
os.popen( 'cp ' + INPUT_DATA + ' ' + sFit_output )
os.popen( 'cp ' + INPUT_DATA + ' ' + cFit_TOYS_output )
os.popen( 'cp ' + cFit_XML + ' ' + cFit_output )
os.popen( 'cp ' + sFit_XML + ' ' + sFit_output )
os.popen( 'cp ' + cFit_TOY + ' ' + cFit_TOYS_output )
os.popen( 'cp ' + Rapid_Bin + ' ' + cFit_output )
os.popen( 'cp ' + Rapid_Bin + ' ' + sFit_output )
os.popen( 'cp ' + Rapid_Bin + ' ' + cFit_TOYS_output )


print "Setting up and Starting Fits"

logging.basicConfig(level=logging.INFO)

cFit_Log = logging.getLogger( "cLog" )
cFit_Handler = logging.FileHandler( cFit_logfile )
cFit_Log.addHandler( cFit_Handler )
sFit_Log = logging.getLogger( "sLog" )
sFit_Handler = logging.FileHandler( sFit_logfile )
sFit_Log.addHandler( sFit_Handler )
TOY_Log = logging.getLogger( "TOYS" )
TOY_Handler = logging.FileHandler( cFit_TOYS_logout )
TOY_Log.addHandler( TOY_Handler )

os.chdir( cFit_output )

print "Starting cFit UnitTest:"

os.chdir( cFit_output )

#	Open a new subprocess running the cFit UnitTest, merging the cout and cerr statements (i.e. logging what appears on screen... nobody's perfect)
cFit = subprocess.Popen( RapidFit_args_cFit, stdout = subprocess.PIPE, stderr=subprocess.STDOUT )

os.chdir( sFit_output )

print "Starting sFit UnitTest:"
#       Open a new subprocess running the sFit UnitTest, merging the cout and cerr statements (i.e. logging what appears on screen... nobody's perfect)
sFit = subprocess.Popen( RapidFit_args_sFit, stdout = subprocess.PIPE, stderr=subprocess.STDOUT )

print "Logging output"

while True:
	cFit_line = cFit.stdout.readline()
	sFit_line = sFit.stdout.readline()
	exitcode_cFit = cFit.poll()
	exitcode_sFit = sFit.poll()
	if (not cFit_line) and (exitcode_cFit is not None):
		print "cFit UnitTest Finished, Results are in: " + cFit_output
	if (not sFit_line) and (exitcode_sFit is not None):
		print "sFit UnitTest Finished, Results are in: " + sFit_output
	if ( (not sFit_line) and (exitcode_sFit is not None) ) and ( (not cFit_line) and (exitcode_cFit is not None) ):
		break
	if cFit_line and not exitcode_cFit:
		cFit_line = cFit_line[:-1]
		cFit_Log.info("%s", cFit_line)
	if sFit_line and not exitcode_sFit:
		sFit_line = sFit_line[:-1]
		sFit_Log.info("%s", sFit_line)	



print "Starting TOY UnitTest:"

#	This is unique to all of the 2D scans on a single machine as we now neeed to replicate the conditions that we have on distributed systems
for i in range( 0, total_toys, toys_per_job ):
	#	Need to loop over multiple toys
	cFit_TOY_ARGS = RapidFit_args_cFit_TOY
	cFit_TOY_ARGS.append( str( TOY_SEED ) )
	cFit_TOY = subprocess.Popen( cFit_TOY_ARGS, stdout = subprocess.PIPE, stderr=subprocess.STDOUT )
	#	Logging for this particular toy
	while True:
		cFit_TOY_line = cFit_TOY.stdout.readline()
		exitcode_cFit_TOY = cFit_TOY.poll()
		if (not cFit_TOY_line) and (exitcode_cFit_TOY is not None):
			print "Finished a cFit toy with a Unique Seed"
		if ( (not cFit_TOY_line) and (exitcode_cFit_TOY is not None) ):
			break
		if cFit_TOY_line and not exitcode_cFit_TOY:
			cFit_TOY_line = cFit_TOY_line[:-1]
			cFit_TOY_Log.info("%s", cFit_TOY_line)
	TOY_SEED = TOY_SEED + 1

print
print
print "Finished All UnitTests!"
print "Outputs are in:"
print str(cFit_output) + "\t&\t" + str(sFit_output) + "\t&\t" + str(cFit_TOYS_output)
