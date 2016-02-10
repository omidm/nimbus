CMAKE_MINIMUM_REQUIRED( VERSION 2.8.7 )
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
###########################################################################
#Project-specific Configuration
###########################################################################

#Set a project name
SET(PROJECT_NAME Water)
PROJECT(${PROJECT_NAME})

#Manually list your cpp and h files here - do NOT use GLOB/GLOB_RECURSE
SET(${PROJECT_NAME}_cpp
        ${PROJECT_SOURCE_DIR}/main.cpp
        ${PROJECT_SOURCE_DIR}/water-driver.cpp
        ${PROJECT_SOURCE_DIR}/water-example.cpp
        ${PROJECT_SOURCE_DIR}/advection-velocity.cpp
	)
SET(${PROJECT_NAME}_h
        ${PROJECT_SOURCE_DIR}/water-driver.h
        ${PROJECT_SOURCE_DIR}/water-example.h
        ${PROJECT_SOURCE_DIR}/myinclude.h
        ${PROJECT_SOURCE_DIR}/advection-velocity.h
        ${PROJECT_SOURCE_DIR}/global-repo.h
        ${PROJECT_SOURCE_DIR}/job-context.h
	)
        
#List the libraries you want to link against - including a higher-level
#library automatically includes the lower-level ones
SET(LIBRARIES_TO_LINK PhysBAM_Rendering-OpenGL)
        
#you can force cache variables to be on or off here
#SET(PHYSBAM_USE_CUDA ON CACHE BOOL ON FORCE)
        
INCLUDE(${PROJECT_SOURCE_DIR}/../physbam-lib/Scripts/CMake/Nimbus_Application_Library.cmake)
