#  XML-Config
#  
#  YOU WANT TO EDIT THESE


wanted_X_segments=10
wanted_Y_segments=10


X_PARAM="Phi_s"
Y_PARAM="deltaGamma"

wanted_X_points=20
wanted_Y_points=20
X_max="3.14159"
X_min="-3.14159"
Y_max="1.0"
Y_min="-1.0"



#  RAPIDFIT / XML CONFIG
#
#     YOU NEED TO EDIT THESE

#  FULL EXPLICIT PATH OF RAPIDFIT EXECUTABLE
RAPIDFIT_BINARY="/exports/home/s0957039/trunk/bin/fitting"

#  THESE MUST BE STRUCTURED AS ELEMENTS IN A LIST FOR PYTHON AND MUST BE AT LEAST ONE DUMMY ELEMENT i.e. ''
RAPIDFIT_ARGS="'--doLLcontour'"
OUTPUT_FILES="'LLcontourScanData.root'"



#  THE REST OF THIS IS FULLY AUTOMATED TO CONSTRUCT THE RELEVENT XML AND GANGA EXECUTION SCRIPTS AND PRODUCE A run_me.sh WHICH YOU CAN EXECUTE AS ./run_me.sh AND IT WILL SUBMIT ALL OF THE JOBS FOR YOU!


# Some pre-calculations for the algerbra

# num of decimal places in algebra
floating_precision=5

X_Range=$(echo "scale=$floating_precision; $X_max - $X_min"| bc)
Y_Range=$(echo "scale=$floating_precision; $Y_max - $Y_min"| bc)
X_Step=$(echo "scale=$floating_precision; $X_Range / $wanted_X_segments" | bc)
Y_Step=$(echo "scale=$floating_precision; $Y_Range / $wanted_Y_segments" | bc)
X_Points_per_Segment=$(echo "scale=0; $wanted_X_points / $wanted_X_segments" | bc)
Y_Points_per_Segment=$(echo "scale=0; $wanted_Y_points / $wanted_Y_segments" | bc)
X_Point_Step=$(echo "scale=$floating_precision;  $X_Range / $wanted_X_points" | bc)
Y_Point_Step=$(echo "scale=$floating_precision;  $Y_Range / $wanted_Y_points" | bc)


#  Some XML TAGS for later
_2Dopen="\<Output\>\\n\t\<2DScan\>"
_2Dclose="\t\<\/2DScan\>\\n\<\/Output\>"
X_Param_open="\t\t\<X_Param\>"
X_Param_close="\t\t\<\/X_Param\>"
Y_Param_open="\t\t\\<Y_Param\>"
Y_Param_close="\t\t\<\/Y_Param\>"
X_Name="\t\t\t\<Name\>$X_PARAM\<\/Name\>"
Y_Name="\t\t\t\<Name\>$Y_PARAM\<\/Name\>"
X_Points="\t\t\t\<Points\>$X_Points_per_Segment\<\/Points\>"
Y_Points="\t\t\t\<Points\>$Y_Points_per_Segment\<\/Points\>"

if [ -f "run_me.sh" ]; then rm run_me.sh; fi

#  Loop over All X_Axis configurations
i=0
while [ $i -lt $wanted_X_segments ]; do
	# RESET FILENAME
	filename_x="XML_Output_"
	filename_x=$filename_x$i"_"
	X_TAGS=$X_Name"\n"$X_Points
	X_step_add=0
	#  Calculate Limits
	if [ $i -gt 0 ]; then X_step_add=$X_Point_Step; fi
	limit_X_max=$(echo "scale=$floating_precision; $X_min + ( $i + 1 ) * $X_Step" | bc)
	limit_X_min=$(echo "scale=$floating_precision; $X_min + $i* $X_Step + $X_step_add" | bc)
	# CONSRUCT XML TAGS
	X_MAX="\t\t\t\<Maximum\>$limit_X_max\<\/Maximum\>"
	X_MIN="\t\t\t\<Minimum\>$limit_X_min\<\/Minimum\>"
	X_TAGS=$X_TAGS"\\n"$X_MAX"\\n"$X_MIN
	j=0
	#  Loop over All Y_Axis xonfigurations
	while [ $j -lt $wanted_Y_segments ]; do
		filename=$filename_x$j".xml"
		Y_TAGS=$Y_Name"\n"$Y_Points
	        Y_step_add=0
	        if [ $j -gt 0 ]; then Y_step_add=$Y_Point_Step; fi
		#  Calculate Limits
	        limit_Y_max=$(echo "scale=$floating_precision; $Y_min + ( $j + 1 ) * $Y_Step" | bc)
	        limit_Y_min=$(echo "scale=$floating_precision; $Y_min + $j* $Y_Step + $Y_step_add" | bc)
		# CONSTRUCT XML TAGS
	        Y_MAX="\t\t\t\<Maximum\>$limit_Y_max\<\/Maximum\>"
	        Y_MIN="\t\t\t\<Minimum\>$limit_Y_min\<\/Minimum\>"
	        Y_TAGS=$Y_TAGS"\\n"$Y_MAX"\\n"$Y_MIN

		#  PUT IT ALL TOGETHER
		OUTPUT_XML=$_2Dopen"\\n"$X_Param_open"\\n"$X_TAGS"\\n"$X_Param_close"\\n"$Y_Param_open"\\n"$Y_TAGS"\\n"$Y_Param_close"\\n"$_2Dclose
		#  USE THE TEMPLATE TO CONSTRUCT A NEW XML FILE
		$(sed "s/___CHANGE_ME___/$OUTPUT_XML/g" < XML_TEMPLATE.xml > $filename )
		#  CONSTRUCT THE GANGA SCRIPT
		GANGA_ARGS="run_arg=['-f','$filename',$RAPIDFIT_ARGS]"
		GANGA_JOB="j=Job(application=Executable(exe='$RAPIDFIT_BINARY',args=run_arg))"
		GANGA_INPUT_XML="j.inputsandbox=[File('./$filename')]"
		GANGA_OUTPUT="j.outputsandbox=[$OUTPUT_FILES]"
		GANGA_SUBMIT="j.backend=SGE()\nj.submit()"
		GANGA_FILE="ganga_"$i"_"$j".py"
		echo -e $GANGA_ARGS"\n"$GANGA_JOB"\n"$GANGA_INPUT_XML"\n"$GANGA_OUTPUT"\n"$GANGA_SUBMIT > $GANGA_FILE
		echo "ganga $GANGA_FILE" >> run_me.sh
		chmod +x run_me.sh
	j=$(($j+1))
	done
i=$((i+1))
done

