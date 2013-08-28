CMAKE_MINIMUM_REQUIRED( VERSION 2.8.7 )
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
###########################################################################
#Project-specific Configuration
###########################################################################

#Set a project name
SET(PROJECT_NAME water_test)
PROJECT(${PROJECT_NAME})

#Set Nimbus root and prtobuf files directory
SET(NIMBUS_ROOT ${PROJECT_SOURCE_DIR}/../..)
SET(PROTO_DIR ${PROJECT_SOURCE_DIR}/proto_files)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${NIMBUS_ROOT})
include_directories(${PROTO_DIR})

SET(PROTO_FILES_h
        ${PROTO_DIR}/vector_2d.pb.h
        ${PROTO_DIR}/range_2d.pb.h
        ${PROTO_DIR}/grid_2d.pb.h
        ${PROTO_DIR}/face_array_2d.pb.h
        ${PROTO_DIR}/face_array_ap_2d.pb.h
    )

SET(PROTO_FILES_cc
        ${PROTO_DIR}/vector_2d.pb.cc
        ${PROTO_DIR}/range_2d.pb.cc
        ${PROTO_DIR}/grid_2d.pb.cc
        ${PROTO_DIR}/face_array_2d.pb.cc
        ${PROTO_DIR}/face_array_ap_2d.pb.cc
    )

SET(PROTO_HELPER_cc
        ${PROTO_DIR}/serialize_data_ap.cc
    )

#Manually list your cpp and h files here - do NOT use GLOB/GLOB_RECURSE
SET(${PROJECT_NAME}_cpp
        ${PROJECT_SOURCE_DIR}/water_data_types.cc
        ${PROJECT_SOURCE_DIR}/water_helper.cc
        ${PROJECT_SOURCE_DIR}/water_driver.cc
        ${PROJECT_SOURCE_DIR}/water_app.cc
        ${PROTO_FILES_h}
	)

SET(${PROJECT_NAME}_h
        ${PROJECT_SOURCE_DIR}/physbam_include.h
        ${PROJECT_SOURCE_DIR}/water_data_types.h
        ${PROJECT_SOURCE_DIR}/water_driver.h
        ${PROJECT_SOURCE_DIR}/water_app.h
        ${PROTO_FILES_cc}
        ${PROTO_HELPER_cc}
	)
        
#List the libraries you want to link against - including a higher-level
#library automatically includes the lower-level ones
SET(LIBRARIES_TO_LINK PhysBAM_Rendering-OpenGL)

INCLUDE(${PROJECT_SOURCE_DIR}/../physbam-lib/Scripts/CMake/Nimbus_Application_Library.cmake)
