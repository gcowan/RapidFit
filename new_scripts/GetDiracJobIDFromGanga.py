jobnum = xxx
name=str(jobnum)+'.log'
logfile = open( name, 'w' )
for j in jobs(jobnum).subjobs:
	diracid=j.backend.id
	diracname=str(diracid)+"\n"
	logfile.write(diracname)
logfile.close()
