from DIRAC.Interfaces.API.Dirac import Dirac
from DIRAC.Interfaces.API.Job import Job
import sys,os
dirac = Dirac()
logfile = sys.argv[1]
out_folder = sys.argv[2]
#print dirac.getOutputSandbox(jobid)
#print dirac.getJobOutputData(jobid)

print os.getlogin()

if not os.path.exists(out_folder):
    os.makedirs(out_folder)

os.chdir( out_folder )

i=0
for line in open( logfile, 'r'):
	jobid = int( line )
	if not os.path.exists( out_folder+"/"+str(i) ):
		os.makedirs( out_folder+"/"+str(i) )
		os.chdir( out_folder+"/"+str(i) )
#		print str(jobid)
		print dirac.getOutputSandbox( jobid )
		print dirac.getJobOutputData( jobid )
	i=i+1
	os.chdir( out_folder )
	print str(line)+"\t"+str(i)
