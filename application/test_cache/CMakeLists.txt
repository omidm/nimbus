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
SET(PROJECT_NAME cache_app)
PROJECT(${PROJECT_NAME})

#Set Nimbus root and prtobuf files directory
SET(NIMBUS_ROOT ${PROJECT_SOURCE_DIR}/../..)

include_directories(${NIMBUS_ROOT})
LINK_DIRECTORIES(${NIMBUS_ROOT})

######################################################
# Manually list your cpp and h files here - do NOT   #
# use GLOB/GLOB_RECURSE                              #
######################################################

#Manually list your cpp and h files here - do NOT use GLOB/GLOB_RECURSE
SET(${PROJECT_NAME}_cpp
    ${PROJECT_SOURCE_DIR}/cache_face_array.cc
    # ${PROJECT_SOURCE_DIR}/cache_scalar_array.cc
    ${PROJECT_SOURCE_DIR}/cache_test.cc
    ${PROJECT_SOURCE_DIR}/data_face_array.cc
    # ${PROJECT_SOURCE_DIR}/data_scalar_array.cc
    ${PROJECT_SOURCE_DIR}/physbam_tools.cc
    )
SET(${PROJECT_NAME}_h
    ${PROJECT_SOURCE_DIR}/cache_face_array.h
    # ${PROJECT_SOURCE_DIR}/cache_scalar_array.h
    ${PROJECT_SOURCE_DIR}/cache_test.h
    ${PROJECT_SOURCE_DIR}/data_face_array.h
    # ${PROJECT_SOURCE_DIR}/data_scalar_array.h
    ${PROJECT_SOURCE_DIR}/physbam_include.h
    ${PROJECT_SOURCE_DIR}/physbam_tools.h
    )

#List the libraries you want to link against - including a higher-level
#library automatically includes the lower-level ones
SET(LIBRARIES_TO_LINK PhysBAM_Rendering-OpenGL)
SET(PROTO_LIB protobuf)
SET(NIMBUS_LIB nimbus)
SET(BOOST_THREAD boost_thread)
SET(BOOST_SYSTEM boost_system)
LIST(APPEND LIBRARIES_TO_LINK ${PROTO_LIB})
LIST(APPEND LIBRARIES_TO_LINK ${NIMBUS_LIB})
LIST(APPEND LIBRARIES_TO_LINK ${BOOST_THREAD})
LIST(APPEND LIBRARIES_TO_LINK ${BOOST_SYSTEM})

#you can force cache variables to be on or off here
#SET(PHYSBAM_USE_CUDA ON CACHE BOOL ON FORCE)

INCLUDE(${NIMBUS_ROOT}/application/physbam-lib/Scripts/CMake/Nimbus_Application_Library.cmake)
