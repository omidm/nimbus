CMAKE_MINIMUM_REQUIRED( VERSION 2.8.7 )
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
###########################################################################
#Project-specific Configuration
###########################################################################

#Set a project name
SET(PROJECT_NAME test)
PROJECT(${PROJECT_NAME})

#Set Nimbus root and prtobuf files directory
SET(NIMBUS_ROOT ${PROJECT_SOURCE_DIR}/../..)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${NIMBUS_ROOT})

#Manually list your cpp and h files here - do NOT use GLOB/GLOB_RECURSE
SET(${PROJECT_NAME}_cpp
        ${PROJECT_SOURCE_DIR}/main.cc)
SET(${PROJECT_NAME}_h
        ${PROJECT_SOURCE_DIR}/translator_physbam_test.h)
        
#List the libraries you want to link against - including a higher-level
#library automatically includes the lower-level ones
SET(LIBRARIES_TO_LINK PhysBAM_Rendering-OpenGL)
LIST(APPEND LIBRARIES_TO_LINK ${PROJECT_SOURCE_DIR}/../../libnimbus.so)
LIST(APPEND LIBRARIES_TO_LINK "/usr/local/lib/libboost_system.dylib")
        
#you can force cache variables to be on or off here
#SET(PHYSBAM_USE_CUDA ON CACHE BOOL ON FORCE)
        
#INCLUDE(${PROJECT_SOURCE_DIR}/../../Scripts/CMake/Project.cmake)
INCLUDE(${PROJECT_SOURCE_DIR}/../physbam-lib/Scripts/CMake/Nimbus_Application.cmake)
