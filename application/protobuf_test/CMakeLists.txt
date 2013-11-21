CMAKE_MINIMUM_REQUIRED( VERSION 2.8.7 )
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
###########################################################################
#Project-specific Configuration
###########################################################################

#Set a project name
SET(PROJECT_NAME projection_test)
PROJECT(${PROJECT_NAME})

#Set Nimbus root and prtobuf files directory
SET(NIMBUS_ROOT ${PROJECT_SOURCE_DIR}/../..)
SET(PROTO_DIR ${PROJECT_SOURCE_DIR}/protocol_buffer)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${NIMBUS_ROOT})
include_directories(${PROTO_DIR})

SET(PROTO_FILES_h
        ${PROTO_DIR}/Sparse_Matrix_Float.pb.h
   )

SET(PROTO_FILES_cc
        ${PROTO_DIR}/Sparse_Matrix_Float.pb.cc
   )

SET(PROTO_HELPER_h
   )

SET(PROTO_HELPER_cc
   )

#Manually list your cpp and h files here - do NOT use GLOB/GLOB_RECURSE
SET(${PROJECT_NAME}_cpp
        ${PROJECT_SOURCE_DIR}/ProtoBuf_test.cc
        ${PROTO_FILES_h}
        ${PROTO_HELPER_h}
	)

SET(${PROJECT_NAME}_h
        ${PROJECT_SOURCE_DIR}/ProtoBuf_test.h
        ${PROTO_FILES_cc}
        ${PROTO_HELPER_cc}
	)
        
#List the libraries you want to link against - including a higher-level
#library automatically includes the lower-level ones
SET(PROTO_LIB protobuf)
SET(LIBRARIES_TO_LINK PhysBAM_Rendering-OpenGL)
LIST(APPEND LIBRARIES_TO_LINK ${PROTO_LIB})

INCLUDE(${PROJECT_SOURCE_DIR}/../physbam-lib/Scripts/CMake/Nimbus_Application_Library.cmake)
