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
        ${PROTO_DIR}/adv_vel_par.pb.h
        ${PROTO_DIR}/physbam_face_array_2d.pb.h
        ${PROTO_DIR}/physbam_grid_2d.pb.h
        ${PROTO_DIR}/physbam_range_2d.pb.h
        ${PROTO_DIR}/physbam_vector_2d.pb.h
    )

SET(PROTO_FILES_cc
        ${PROTO_DIR}/adv_vel_par.pb.cc
        ${PROTO_DIR}/physbam_face_array_2d.pb.cc
        ${PROTO_DIR}/physbam_grid_2d.pb.cc
        ${PROTO_DIR}/physbam_range_2d.pb.cc
        ${PROTO_DIR}/physbam_vector_2d.pb.cc
    )

SET(PROTO_HELPER_h
        ${PROTO_DIR}/pb_include_2d.h
        ${PROTO_DIR}/physbam_deserialize_data_arrays_2d.h
        ${PROTO_DIR}/physbam_deserialize_data_common_2d.h
        ${PROTO_DIR}/physbam_serialize_data_arrays_2d.h
        ${PROTO_DIR}/physbam_serialize_data_common_2d.h
    )

SET(PROTO_HELPER_cc
        ${PROTO_DIR}/physbam_deserialize_data_arrays_2d.cc
        ${PROTO_DIR}/physbam_deserialize_data_common_2d.cc
        ${PROTO_DIR}/physbam_serialize_data_arrays_2d.cc
        ${PROTO_DIR}/physbam_serialize_data_common_2d.cc
    )

#Manually list your cpp and h files here - do NOT use GLOB/GLOB_RECURSE
SET(${PROJECT_NAME}_cpp
        ${PROJECT_SOURCE_DIR}/advection.cc
        ${PROJECT_SOURCE_DIR}/data_face_arrays.cc
        ${PROJECT_SOURCE_DIR}/data_utils.cc
        ${PROJECT_SOURCE_DIR}/water_app.cc
        ${PROJECT_SOURCE_DIR}/water_data_driver.cc
        ${PROJECT_SOURCE_DIR}/water_driver.cc
        ${PROJECT_SOURCE_DIR}/water_helper.cc
        ${PROTO_FILES_h}
        ${PROTO_HELPER_h}
	)

SET(${PROJECT_NAME}_h
        ${PROJECT_SOURCE_DIR}/advection.h
        ${PROJECT_SOURCE_DIR}/app_config.h
        ${PROJECT_SOURCE_DIR}/app_utils.h
        ${PROJECT_SOURCE_DIR}/data_face_arrays.h
        ${PROJECT_SOURCE_DIR}/data_fwd_decl.h
        ${PROJECT_SOURCE_DIR}/data_utils.h
        ${PROJECT_SOURCE_DIR}/physbam_include.h
        ${PROJECT_SOURCE_DIR}/water_app.h
        ${PROJECT_SOURCE_DIR}/water_data_driver.h
        ${PROJECT_SOURCE_DIR}/water_driver.h
        ${PROTO_FILES_cc}
        ${PROTO_HELPER_cc}
	)
        
#List the libraries you want to link against - including a higher-level
#library automatically includes the lower-level ones
SET(PROTO_LIB protobuf)
SET(LIBRARIES_TO_LINK PhysBAM_Rendering-OpenGL)
LIST(APPEND LIBRARIES_TO_LINK ${PROTO_LIB})

INCLUDE(${PROJECT_SOURCE_DIR}/../physbam-lib/Scripts/CMake/Nimbus_Application_Library.cmake)
