CMAKE_MINIMUM_REQUIRED( VERSION 2.8.7 )
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
###########################################################################
#Project-specific Configuration
###########################################################################

#Set a project name
SET(PROJECT_NAME projection)
PROJECT(${PROJECT_NAME})

#Set Nimbus root and prtobuf files directory
SET(NIMBUS_ROOT ${PROJECT_SOURCE_DIR}/../../..)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${NIMBUS_ROOT})

#Manually list your cpp and h files here - do NOT use GLOB/GLOB_RECURSE
SET(${PROJECT_NAME}_cpp
        ${PROJECT_SOURCE_DIR}/app.cc
        ${PROJECT_SOURCE_DIR}/init_main.cc
        ${PROJECT_SOURCE_DIR}/projection_driver.cc
        ${PROJECT_SOURCE_DIR}/projection_example.cc
        ${PROJECT_SOURCE_DIR}/projection_data.cc
        ${PROJECT_SOURCE_DIR}/nimbus_pcg_sparse_mpi.cc)
SET(${PROJECT_NAME}_h
        ${PROJECT_SOURCE_DIR}/app.h
        ${PROJECT_SOURCE_DIR}/init_main.h
        ${PROJECT_SOURCE_DIR}/projection_driver.h
        ${PROJECT_SOURCE_DIR}/projection_example.h
        ${PROJECT_SOURCE_DIR}/projection_data.h
        ${PROJECT_SOURCE_DIR}/nimbus_pcg_sparse_mpi.h)
        
#List the libraries you want to link against - including a higher-level
#library automatically includes the lower-level ones
SET(LIBRARIES_TO_LINK PhysBAM_Rendering-OpenGL)
        
#you can force cache variables to be on or off here
#SET(PHYSBAM_USE_CUDA ON CACHE BOOL ON FORCE)
        
#INCLUDE(${PROJECT_SOURCE_DIR}/../../Scripts/CMake/Project.cmake)
INCLUDE(${PROJECT_SOURCE_DIR}/../../physbam-new-lib/Scripts/CMake/Nimbus_Application_Library.cmake)
