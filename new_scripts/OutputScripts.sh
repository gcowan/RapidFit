#!/bin/bash
#
#	This is intended to be run ./OutputScripts.sh as jobID
#
#	This uses the new Dirac API which has split out the get output and output-data into nice scripts
#
#	This script, DOES NOT test wether you have a valid proxy
#

num=0
start=$PWD

while read line; do
	#echo $line
	#dirac-wms-job-get-output-data 
	#dirac-wms-job-get-output 
	#dirac-wmd-job-get-input 

	sleep 0.5;

	if [ -d "$2/$1/$num" ]; then
		#echo "moving on from " $line
		sleep 0.5;
	else
		mkdir -p $2/$1/$num
		cd $2/$1/$num
		#echo $PWD
		dirac-wms-job-get-output-data $line
		dirac-wms-job-get-output $line
		#dirac-wms-job-get-input $line
	fi

	num=$(($num+1))
	cd $start
done < $1.log


