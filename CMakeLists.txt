#  CMAKE CONFIG

#  Some Defaults to allow people to just use the file

#	Would REALLY like to make this 2.8 but 2.6.x is provided in SC_Linux 5.x
#	2.8 handles system includes in a MUCH better way
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )
MARK_AS_ADVANCED( CMAKE_BACKWARDS_COMPATABILITY )
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )



STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" INSOURCE)
IF( INSOURCE )
  MESSAGE( " " )
  MESSAGE( FATAL_ERROR "\nDO NOT BUILD RAPIDFIT IN THE SOURCE DIRECTORY!!!\n\nI want to keep the Makefile for Historic reasons and backwards compatability :D" )
ENDIF( INSOURCE )


#  COMPILER OPTIONS

#  -O3 optimizations
#SET( CMAKE_BUILD_TYPE Release )
#  -O2 optimizations and -g
SET( CMAKE_BUILD_TYPE Debug )

#  Get back the old style verbosity :) by setting this to 1
#SET( CMAKE_VERBOSE_MAKEFILE 1 )
SET( CMAKE_VERBOSE_MAKEFILE 0 )



#  PROJECT OPTIONS

#Set Project name and version

PROJECT( RapidFit )

SET( ${PROJECT_NAME}_MAJOR_VERSION 0 )
SET( ${PROJECT_NAME}_MINOR_VERSION 1 )
SET( ${PROJECT_NAME}_PATCH_LEVEL 0 )

EXECUTE_PROCESS(
       		COMMAND	svnversion -n ..
		OUTPUT_VARIABLE SVN_REV )
MESSAGE( "SVN VERSION:  ${SVN_REV}" )

EXECUTE_PROCESS(
		COMMAND date +%Y_%m_%d
		OUTPUT_VARIABLE BUILD_DATE )

#  Grab the info on ROOT and store the global prarmaters
SET( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/scripts )
FIND_PACKAGE( ROOT )

#  INCLUDE & LINKING DIRS

INCLUDE_DIRECTORIES( "${PROJECT_SOURCE_DIR}/framework/include" )
INCLUDE_DIRECTORIES( "${PROJECT_SOURCE_DIR}/pdfs/include" )
INCLUDE_DIRECTORIES( "${PROJECT_SOURCE_DIR}/utils/include" )

#	Apple is full of fomplete... WHAT THE HELL?
IF( NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
	INCLUDE_DIRECTORIES( SYSTEM "${ROOT_INCLUDE_DIR}" )
ENDIF( NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )


IF( NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
	SET( CMAKE_CXX_FLAGS "-DSVN_REV=${SVN_REV} -Weffc++ -Os -msse -msse2 -msse3 -m3dnow -ansi -fPIC -fmerge-all-constants -funroll-all-loops -g -D__ROOFIT_NOBANNER -Wconversion -Weffc++ -Wextra -Wsign-compare -Wfloat-equal -Wmissing-noreturn -Wall -Wno-non-virtual-dtor" )
ENDIF( NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )

IF( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
	SET( CMAKE_CXX_FLAGS "-DSVN_REV=${SVN_REV} -fPIC -O3 -msse -msse2 -msse3 -funroll-all-loops -m3dnow -ansi -g -Wextra -Wall -Wno-non-virtual-dtor" )
	SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ${ROOT_INCLUDE_DIR}" )
ENDIF( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )

LINK_DIRECTORIES( ${ROOT_LIBRARY_DIR} )
LINK_DIRECTORIES( ${LINK_DIRECTORIES} )

SET( LIBRARY_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/lib" )
SET( EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin" )



#  PROJECT SOURCES

#  Define the PDF Source Files to be Built in this case all source files in the directory
AUX_SOURCE_DIRECTORY( ${PROJECT_SOURCE_DIR}/pdfs/src PDF_SRCS )
AUX_SOURCE_DIRECTORY( ${PROJECT_SOURCE_DIR}/framework/src FITTER_SRCS )


#  This is here due to 'broken' sources called Roo* I have to remove these from the list of source files in framework/srcs
#  If these files weren't here I could remove this code which is less than 

SET( ALL_HEADERS )
LIST(LENGTH FITTER_SRCS SRC_LEN )
MATH(EXPR SRC_LEN ${SRC_LEN}-1 )
# This forms the equivalent of a for loop in cmake and removes all list elements containing the string "Roo"
SET( I 0 )
WHILE( ${I} LESS ${SRC_LEN} )
	MATH(EXPR PRESENT_FILENUM ${I}+1 )
	LIST( GET FITTER_SRCS ${PRESENT_FILENUM} THIS_FILE )
	IF("${THIS_FILE}" MATCHES "Roo+")
		LIST( REMOVE_AT FITTER_SRCS ${PRESENT_FILENUM} )
		MATH( EXPR I ${I}-1 )
		MATH( EXPR SRC_LEN ${SRC_LEN}-1 )
	ELSEIF(NOT "${THIS_FILE}" MATCHES "Roo+" )
		IF("${THIS_FILE}" MATCHES "${PROJECT_SOURCE_DIR}/framework/src/unused+")
                	LIST( REMOVE_AT FITTER_SRCS ${PRESENT_FILENUM} )
                	MATH( EXPR I ${I}-1 )
                	MATH( EXPR SRC_LEN ${SRC_LEN}-1 )
        	ENDIF("${THIS_FILE}" MATCHES "${PROJECT_SOURCE_DIR}/framework/src/unused+")
	ENDIF("${THIS_FILE}" MATCHES "Roo+")
        SET(HEADER_IN "${THIS_FILE}")
	STRING(REGEX REPLACE ".cpp" ".h" HEADER_OUT "${HEADER_IN}" )
	STRING(REGEX REPLACE "framework/src/" "framework/include/" FINAL_HEADER "${HEADER_OUT}")
	LIST( APPEND ALL_HEADERS ${FINAL_HEADER} )
	IF("${THIS_FILE}" MATCHES "RapidRun+")
                LIST( REMOVE_AT FITTER_SRCS ${PRESENT_FILENUM} )
                MATH( EXPR I ${I}-1 )
                MATH( EXPR SRC_LEN ${SRC_LEN}-1 )
        ELSEIF(NOT "${THIS_FILE}" MATCHES "RapidRun+" )
                IF("${THIS_FILE}" MATCHES "${PROJECT_SOURCE_DIR}/framework/src/RapidRun+")
                        LIST( REMOVE_AT FITTER_SRCS ${PRESENT_FILENUM} )
                        MATH( EXPR I ${I}-1 )
                        MATH( EXPR SRC_LEN ${SRC_LEN}-1 )
                ENDIF("${THIS_FILE}" MATCHES "${PROJECT_SOURCE_DIR}/framework/src/RapidRun+")
        ENDIF("${THIS_FILE}" MATCHES "RapidRun+")
	SET( I ${PRESENT_FILENUM}  )
ENDWHILE( ${I} LESS ${SRC_LEN} )


#  OUTPUT OPTIONS

ADD_LIBRARY( fits SHARED ${FITTER_SRCS} )
TARGET_LINK_LIBRARIES( fits ${ROOT_LIBRARIES} ${ROOT_FIT_LIBRARIES} )

ADD_LIBRARY( pdfs SHARED ${PDF_SRCS} )
TARGET_LINK_LIBRARIES( pdfs ${ROOT_LIBRARIES} )

#  Add the main program which offers the standard way into RapidFit and the usecase where most users will want it
ADD_EXECUTABLE( fitting ${PROJECT_SOURCE_DIR}/framework/src/main.cpp )
TARGET_LINK_LIBRARIES( fitting ${ROOT_LIBRARIES} fits pdfs )






#  Add Conor's and My tools to the build to be able to get these out of a make

#ADD_EXECUTABLE( rapidresults ${PROJECT_SOURCE_DIR}/utils/src/rapidresults.cc )
#TARGET_LINK_LIBRARIES( rapidresults ${ROOT_LIBRARIES} )

#	Depreciated tool
#ADD_EXECUTABLE( rapidfit_llscanresults ${PROJECT_SOURCE_DIR}/utils/src/rapidfit_llscanresults.cc )
#TARGET_LINK_LIBRARIES( rapidfit_llscanresults ${ROOT_LIBRARIES} )

#	Old tool for producing FC plots
ADD_EXECUTABLE( rapidfit_fcscanresults ${PROJECT_SOURCE_DIR}/utils/src/rapidfit_fcscanresults.cc )
TARGET_LINK_LIBRARIES( rapidfit_fcscanresults ${ROOT_LIBRARIES} )

#	Newest tool, currently broken for FC plots
ADD_EXECUTABLE( rapidfit_fcscanresults_2 ${PROJECT_SOURCE_DIR}/utils/src/rapidfit_fcscanresults_2.cc )
TARGET_LINK_LIBRARIES( rapidfit_fcscanresults_2 ${ROOT_LIBRARIES} )

#	Conors tool for analysing toy studies
#	main difference it plots the actual fit results
ADD_EXECUTABLE( rapidfit_toyresults ${PROJECT_SOURCE_DIR}/utils/src/rapidfit_toyresults.cc )
TARGET_LINK_LIBRARIES( rapidfit_toyresults ${ROOT_LIBRARIES} )

#	Tool for comparing results of 2DLL within Edinburgh
ADD_EXECUTABLE( merge_plot ${PROJECT_SOURCE_DIR}/utils/src/merge_plot.C )
TARGET_LINK_LIBRARIES( merge_plot ${ROOT_LIBRARIES} )

#       Tool that sWeights the betas ntuple
ADD_EXECUTABLE( betas_sweightfitter ${PROJECT_SOURCE_DIR}/utils/src/betas_sweightfitter.cc )
TARGET_LINK_LIBRARIES( betas_sweightfitter ${ROOT_LIBRARIES} ${ROOT_FIT_LIBRARIES} )




#	Libraries which share a LOT of code useful in the analysis of the results from RapidFit

ADD_LIBRARY( NTuple_Processing SHARED ${PROJECT_SOURCE_DIR}/utils/src/NTuple_Processing.C )
TARGET_LINK_LIBRARIES( NTuple_Processing ${ROOT_LIBRARIES} )

ADD_LIBRARY( TString_Processing SHARED ${PROJECT_SOURCE_DIR}/utils/src/TString_Processing.C ${PROJECT_SOURCE_DIR}/framework/src/StringProcessing.cpp )
TARGET_LINK_LIBRARIES( TString_Processing ${ROOT_LIBRARIES} )

ADD_LIBRARY( Histo_Processing SHARED ${PROJECT_SOURCE_DIR}/utils/src/Histo_Processing.C )
TARGET_LINK_LIBRARIES( Histo_Processing ${ROOT_LIBRARIES} )

ADD_LIBRARY( RapidFit_Style SHARED ${PROJECT_SOURCE_DIR}/framework/src/EdStyle.cpp )
TARGET_LINK_LIBRARIES( RapidFit_Style ${ROOT_LIBRARIES} )



#       My tool for analysing toy studies
#       main difference is it performs a binned analysis
ADD_EXECUTABLE( tinter ${PROJECT_SOURCE_DIR}/utils/src/tinter.C )
TARGET_LINK_LIBRARIES( tinter ${ROOT_LIBRARIES} RapidFit_Style Histo_Processing TString_Processing NTuple_Processing )

#       New tool for plotting 2DLL and FC
ADD_EXECUTABLE( RapidPlot ${PROJECT_SOURCE_DIR}/utils/src/RapidPlot.C )
TARGET_LINK_LIBRARIES( RapidPlot ${ROOT_LIBRARIES} RapidFit_Style Histo_Processing TString_Processing NTuple_Processing )

#       New tool for plotting 1D LL and overlaying multiple copies of the same
ADD_EXECUTABLE( RapidLL ${PROJECT_SOURCE_DIR}/utils/src/RapidLL.C )
TARGET_LINK_LIBRARIES( RapidLL ${ROOT_LIBRARIES} RapidFit_Style Histo_Processing TString_Processing NTuple_Processing )

#       New tool for plotting 1D LL and overlaying multiple copies of the same
ADD_EXECUTABLE( print ${PROJECT_SOURCE_DIR}/utils/src/print.C )
TARGET_LINK_LIBRARIES( print ${ROOT_LIBRARIES} RapidFit_Style Histo_Processing TString_Processing NTuple_Processing )





#  Add the custom command to 'make doc' incase this was/is still used by anyone
add_custom_target (doc doxygen RapidFit_doxygen.cfg && tar cvfz ./RapidFit_html.tgz ./html/ && scp ./RapidFit_html.tgz ph-ppe.ph.ed.ac.uk:~/WWW/RapidFit/doc/  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/doc COMMENT "Generating documentation" OUTPUT docs/RapidFit_html.tgz docs/html)
SET_DIRECTORY_PROPERTIES( PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES"${PROJECT_SOURCE_DIR}/doc/RapidFit_html.tgz;${PROJECT_SOURCE_DIR}/doc/html" )


#	Add the new files to the Headers that are required for building
#	RapidFit as a library extention to ROOT
SET( HEADERS "${PROJECT_SOURCE_DIR}/framework/include/RapidRun.h" )
SET( LINKDEF "${PROJECT_SOURCE_DIR}/framework/include/LinkDef.h" )

#	This is a custom command in cmake required to produce the psuedo-code
#	that RapidFit requires to build the fitter as a library extension to
#	ROOT to run on the grid
ADD_CUSTOM_COMMAND ( OUTPUT  ${PROJECT_BINARY_DIR}/rapidfit_dict.cpp
		    COMMAND rootcint -f ${PROJECT_BINARY_DIR}/rapidfit_dict.cpp -c ${ALL_HEADERS} ${LINKDEF}
		    DEPENDS ${ALL_HEADERS} )

#	This is what is required to build the object libRapidRun.so
#	It requires the custom command above
ADD_LIBRARY ( RapidRun SHARED ${PROJECT_BINARY_DIR}/rapidfit_dict.cpp ${PROJECT_SOURCE_DIR}/framework/src/RapidRun.cpp )
TARGET_LINK_LIBRARIES( RapidRun ${ROOT_LIBRARIES} ${ROOT_FIT_LIBRARIES} pdfs fits )
