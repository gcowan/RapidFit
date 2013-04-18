#
#	Provide the JOB_NUM which matches the file JOB_NUM.log which is local and contains a series of DIRAC jobID values
#	
#	This script will request the output from each of these values and will store the output in a folder /tmp/username/JOB_NUM
#
#	This will perform as many parallel requests for job outputs as you want by adjusting parallel_requests, 50 is normally excessive but fast
#

export JOB_NUM="xxx"

parallel_requests=50
username=$(whoami)

for i in $JOB_NUM; do
	for j in $(seq 1 ${parallel_requests} ); do
		sleep 0.$j && $PWD/OutputScripts.sh ${i} /tmp/${username}/${i} &
	done
done

