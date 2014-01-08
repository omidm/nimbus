IF(APPLE)
    SET(CMAKE_C_COMPILER /usr/local/bin/gcc-4.5)
    SET(CMAKE_CXX_COMPILER /usr/local/bin/g++-4.5)
ENDIF()

CMAKE_MINIMUM_REQUIRED( VERSION 2.8.7 )
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
###########################################################################
#Project-specific Configuration
###########################################################################

#Set a project name
SET(PROJECT_NAME water_app)
PROJECT(${PROJECT_NAME})

#Set Nimbus root and prtobuf files directory
SET(NIMBUS_ROOT ${PROJECT_SOURCE_DIR}/../..)

include_directories(${NIMBUS_ROOT})
LINK_DIRECTORIES(${NIMBUS_ROOT})

#Manually list your cpp and h files here - do NOT use GLOB/GLOB_RECURSE
SET(${PROJECT_NAME}_cpp
    ${PROJECT_SOURCE_DIR}/app_utils.cc
    ${PROJECT_SOURCE_DIR}/data_app.cc
    ${PROJECT_SOURCE_DIR}/job_initialize.cc
    ${PROJECT_SOURCE_DIR}/job_iteration.cc
    ${PROJECT_SOURCE_DIR}/job_loop.cc
    ${PROJECT_SOURCE_DIR}/job_main.cc
    ${PROJECT_SOURCE_DIR}/job_super_2.cc
    ${PROJECT_SOURCE_DIR}/water_app.cc
    ${PROJECT_SOURCE_DIR}/water_driver.cc
    ${PROJECT_SOURCE_DIR}/water_example.cc
    ${PROJECT_SOURCE_DIR}/water_sources.cc
    )
SET(${PROJECT_NAME}_h
    ${PROJECT_SOURCE_DIR}/app_utils.h
    ${PROJECT_SOURCE_DIR}/data_app.h
    ${PROJECT_SOURCE_DIR}/job_initialize.h
    ${PROJECT_SOURCE_DIR}/job_iteration.h
    ${PROJECT_SOURCE_DIR}/job_loop.h
    ${PROJECT_SOURCE_DIR}/job_main.h
    ${PROJECT_SOURCE_DIR}/job_super_2.h
    ${PROJECT_SOURCE_DIR}/physbam_include.h
    ${PROJECT_SOURCE_DIR}/water_app.h
    ${PROJECT_SOURCE_DIR}/water_driver.h
    ${PROJECT_SOURCE_DIR}/water_example.h
    ${PROJECT_SOURCE_DIR}/water_sources.h
    )

#List the libraries you want to link against - including a higher-level
#library automatically includes the lower-level ones
SET(LIBRARIES_TO_LINK PhysBAM_Rendering-OpenGL)
SET(PROTO_LIB protobuf)
SET(NIMBUS_LIB nimbus)
SET(BOOST_THREAD boost_thread-mt)
SET(BOOST_SYSTEM boost_system-mt)
LIST(APPEND LIBRARIES_TO_LINK ${PROTO_LIB})
LIST(APPEND LIBRARIES_TO_LINK ${NIMBUS_LIB})
LIST(APPEND LIBRARIES_TO_LINK ${BOOST_THREAD})
LIST(APPEND LIBRARIES_TO_LINK ${BOOST_SYSTEM})

#you can force cache variables to be on or off here
#SET(PHYSBAM_USE_CUDA ON CACHE BOOL ON FORCE)

INCLUDE(${NIMBUS_ROOT}/application/physbam-lib/Scripts/CMake/Nimbus_Application_Library.cmake)
