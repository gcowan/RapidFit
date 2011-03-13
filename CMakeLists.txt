#  CMAKE CONFIG

#  Some Defaults to allow people to just use the file

CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )
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
#SET( CMAKE_BUILD_TYPE Debug )

SET( CMAKE_CXX_FLAGS_RELEASE "" )

SET( CMAKE_CXX_FLAGS "-O3 -mtune=native -ffast-math -g -ansi -fPIC -funroll-all-loops -D__ROOFIT_NOBANNER -Wconversion -Wextra -Wsign-compare -Wfloat-equal -Wmissing-noreturn -Wall -Wno-non-virtual-dtor" )

#  Get back the old style verbosity :) by setting this to 1
#SET( CMAKE_VERBOSE_MAKEFILE 1 )
SET( CMAKE_VERBOSE_MAKEFILE 0 )



#  PROJECT OPTIONS

#Set Project name and version

PROJECT( RapidFit )

SET( ${PROJECT_NAME}_MAJOR_VERSION 0 )
SET( ${PROJECT_NAME}_MINOR_VERSION 1 )
SET( ${PROJECT_NAME}_PATCH_LEVEL 0 )



#  ROOT OPTIONS

#  Grab the info on ROOT and store the global prarmaters
SET( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/scripts )
FIND_PACKAGE( ROOT )



#  INCLUDE & LINKING DIRS

INCLUDE_DIRECTORIES( "${PROJECT_SOURCE_DIR}/framework/include" )
INCLUDE_DIRECTORIES( "${PROJECT_SOURCE_DIR}/pdfs/include" )
INCLUDE_DIRECTORIES( SYSTEM ${ROOT_INCLUDE_DIR} )

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
		IF("${THIS_FILE}" MATCHES "${PROJECT_SOURCE_DIR}/framework/src/Roo+")
                	LIST( REMOVE_AT FITTER_SRCS ${PRESENT_FILENUM} )
                	MATH( EXPR I ${I}-1 )
                	MATH( EXPR SRC_LEN ${SRC_LEN}-1 )
        	ENDIF("${THIS_FILE}" MATCHES "${PROJECT_SOURCE_DIR}/framework/src/Roo+")
	ENDIF("${THIS_FILE}" MATCHES "Roo+")
        IF("${THIS_FILE}" MATCHES "main+")
                LIST( REMOVE_AT FITTER_SRCS ${PRESENT_FILENUM} )
                MATH( EXPR I ${I}-1 )
                MATH( EXPR SRC_LEN ${SRC_LEN}-1 )
	ELSEIF(NOT "${THIS_FILE}" MATCHES "main+" )
		IF("${THIS_FILE}" MATCHES "${PROJECT_SOURCE_DIR}/framework/src/main+")
                	LIST( REMOVE_AT FITTER_SRCS ${PRESENT_FILENUM} )
                	MATH( EXPR I ${I}-1 )
                	MATH( EXPR SRC_LEN ${SRC_LEN}-1 )
        	ENDIF("${THIS_FILE}" MATCHES "${PROJECT_SOURCE_DIR}/framework/src/main+")
        ENDIF("${THIS_FILE}" MATCHES "main+")
	SET( I ${PRESENT_FILENUM}  )
ENDWHILE( ${I} LESS ${SRC_LEN} )



#  OUTPUT OPTIONS

#  Add the pdfs to a Static Library to be accessed by the main program
#  (also less of a pain and more clearly removes it as a dependence when rebuilding)
ADD_LIBRARY( pdfs SHARED ${PDF_SRCS} )
TARGET_LINK_LIBRARIES( pdfs ${ROOT_LIBRARIES})

#  Add the main RapidFit architecture to a simple new fits library to allow for dynamically linking into segments of RapidFit Code
ADD_LIBRARY( fits SHARED ${FITTER_SRCS} )
TARGET_LINK_LIBRARIES( fits ${ROOT_LIBRARIES} ${ROOT_FIT_LIBRARIES} )

#  Add the main program which offers the standard way into RapidFit and the usecase where most users will want it
ADD_EXECUTABLE( fitting ${PROJECT_SOURCE_DIR}/framework/src/main.cpp )
TARGET_LINK_LIBRARIES( fitting fits pdfs ${ROOT_LIBRARIES} ${ROOT_FIT_LIBRARIES} )



#  Add Conor's and My tools to the build to be able to get these out of a make

#ADD_EXECUTABLE( rapidresults ${PROJECT_SOURCE_DIR}/utils/src/rapidresults.cc )
#TARGET_LINK_LIBRARIES( rapidresults ${ROOT_LIBRARIES} )

#ADD_EXECUTABLE( rapidfit_llscanresults ${PROJECT_SOURCE_DIR}/utils/src/rapidfit_llscanresults.cc )
#TARGET_LINK_LIBRARIES( rapidfit_llscanresults ${ROOT_LIBRARIES} )

ADD_EXECUTABLE( rapidfit_fcscanresults ${PROJECT_SOURCE_DIR}/utils/src/rapidfit_fcscanresults.cc )
TARGET_LINK_LIBRARIES( rapidfit_fcscanresults ${ROOT_LIBRARIES} )

ADD_EXECUTABLE( rapidfit_toyresults ${PROJECT_SOURCE_DIR}/utils/src/rapidfit_toyresults.cc )
TARGET_LINK_LIBRARIES( rapidfit_toyresults ${ROOT_LIBRARIES} )

ADD_EXECUTABLE( tinter ${PROJECT_SOURCE_DIR}/utils/src/tinter-static.C )
TARGET_LINK_LIBRARIES( tinter ${ROOT_LIBRARIES} )

#ADD_EXECUTABLE( merge_plot ${PROJECT_SOURCE_DIR}/utils/src/merge_plot.C )
#TARGET_LINK_LIBRARIES( merge_plot ${ROOT_LIBRARIES} )


#  Add the custom command to 'make doc' incase this was/is still used by anyone

add_custom_target (doc doxygen RapidFit_doxygen.cfg && tar cvfz ./RapidFit_html.tgz ./html/ && scp ./RapidFit_html.tgz ph-ppe.ph.ed.ac.uk:~/WWW/RapidFit/doc/  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/doc COMMENT "Generating documentation" OUTPUT docs/RapidFit_html.tgz docs/html)
SET_DIRECTORY_PROPERTIES( PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES"${PROJECT_SOURCE_DIR}/doc/RapidFit_html.tgz;${PROJECT_SOURCE_DIR}/doc/html" )
