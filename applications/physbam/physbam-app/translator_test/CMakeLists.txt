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
SET(NIMBUS_ROOT ${PROJECT_SOURCE_DIR}/../../..)
SET(LEVELDB_ROOT ${NIMBUS_ROOT}/leveldb)
include_directories(${NIMBUS_ROOT})
LINK_DIRECTORIES(${NIMBUS_ROOT})
LINK_DIRECTORIES(${LEVELDB_ROOT})

######################################################
# Manually list your cpp and h files here - do NOT   #
# use GLOB/GLOB_RECURSE                              #
######################################################
#Add flags
add_definitions(-DDUMP_NDATA=FALSE )

#Manually list your cpp and h files here - do NOT use GLOB/GLOB_RECURSE
SET(${PROJECT_NAME}_cpp
    ${PROJECT_SOURCE_DIR}/data_face_array.cc
    ${PROJECT_SOURCE_DIR}/data_particle_array.cc
    ${PROJECT_SOURCE_DIR}/face_array_test.cc
    ${PROJECT_SOURCE_DIR}/main.cc
    ${PROJECT_SOURCE_DIR}/particle_test.cc
    ${PROJECT_SOURCE_DIR}/water_driver.cc
    ${PROJECT_SOURCE_DIR}/water_example.cc
    )
SET(${PROJECT_NAME}_h
    ${PROJECT_SOURCE_DIR}/data_face_array.h
    ${PROJECT_SOURCE_DIR}/data_particle_array.h
    ${PROJECT_SOURCE_DIR}/face_array_test.h
    ${PROJECT_SOURCE_DIR}/particle_test.h
    ${PROJECT_SOURCE_DIR}/physbam_include.h
    ${PROJECT_SOURCE_DIR}/water_driver.h
    ${PROJECT_SOURCE_DIR}/water_example.h
    )

#List the libraries you want to link against - including a higher-level
#library automatically includes the lower-level ones
SET(LIBRARIES_TO_LINK PhysBAM_Rendering-OpenGL)
SET(PROTO_LIB protobuf)
SET(NIMBUS_LIB nimbus)
SET(LEVELDB_LIB leveldb)
SET(BOOST_THREAD boost_thread)
SET(BOOST_SYSTEM boost_system)
LIST(APPEND LIBRARIES_TO_LINK ${PROTO_LIB})
LIST(APPEND LIBRARIES_TO_LINK ${NIMBUS_LIB})
LIST(APPEND LIBRARIES_TO_LINK ${LEVELDB_LIB})
LIST(APPEND LIBRARIES_TO_LINK ${BOOST_THREAD})
LIST(APPEND LIBRARIES_TO_LINK ${BOOST_SYSTEM})

#you can force app_data variables to be on or off here
#SET(PHYSBAM_USE_CUDA ON CACHE BOOL ON FORCE)

INCLUDE(${NIMBUS_ROOT}/application/physbam-lib/Scripts/CMake/Nimbus_Application.cmake)
include_directories(${NIMBUS_ROOT}/leveldb/include/)
