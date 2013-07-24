CMAKE_MINIMUM_REQUIRED( VERSION 2.8.7 )
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
###########################################################################
#Project-specific Configuration
###########################################################################

#Set a project name
SET(PROJECT_NAME Water)
PROJECT(${PROJECT_NAME})

#Set Nimbus root
SET(NIMBUS_ROOT ${PROJECT_SOURCE_DIR}/../..)

include_directories(${NIMBUS_ROOT})

#Manually list your cpp and h files here - do NOT use GLOB/GLOB_RECURSE
SET(${PROJECT_NAME}_cpp
        ${PROJECT_SOURCE_DIR}/water_app.cc
	)
SET(${PROJECT_NAME}_h
        ${PROJECT_SOURCE_DIR}/water_app.h
	)
        
#List the libraries you want to link against - including a higher-level
#library automatically includes the lower-level ones
SET(LIBRARIES_TO_LINK PhysBAM_Rendering-OpenGL)

INCLUDE(${PROJECT_SOURCE_DIR}/../physbam-lib/Scripts/CMake/Nimbus_Application_Library.cmake)
