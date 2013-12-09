CMAKE_MINIMUM_REQUIRED( VERSION 2.8.7 )
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
###########################################################################
#Project-specific Configuration
###########################################################################

#Set a project name
SET(PROJECT_NAME projection)
PROJECT(${PROJECT_NAME})

#Set Nimbus root and prtobuf files directory
SET(NIMBUS_ROOT ${PROJECT_SOURCE_DIR}/../..)
SET(PROTO_DIR ${PROJECT_SOURCE_DIR}/protocol_buffer)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${NIMBUS_ROOT})
include_directories(${PROTO_DIR})

SET(PROTO_FILES_h
        ${PROTO_DIR}/Sparse_Matrix_Float.pb.h
	    ${PROTO_DIR}/Vector_Float.pb.h
   )

SET(PROTO_FILES_cc
        ${PROTO_DIR}/Sparse_Matrix_Float.pb.cc
    	${PROTO_DIR}/Vector_Float.pb.cc
   )

#Manually list your cpp and h files here - do NOT use GLOB/GLOB_RECURSE
SET(${PROJECT_NAME}_cpp
        ${PROJECT_SOURCE_DIR}/app.cc
        ${PROJECT_SOURCE_DIR}/data_impl.cc
        ${PROJECT_SOURCE_DIR}/job_impl.cc
        ${PROJECT_SOURCE_DIR}/job_impl_internal.cc
        ${PROJECT_SOURCE_DIR}/projection_driver.cc
        ${PROJECT_SOURCE_DIR}/projection_example.cc
        ${PROJECT_SOURCE_DIR}/projection_data.cc
        ${PROJECT_SOURCE_DIR}/PCG_Sparse_Solver.cc
        ${PROJECT_SOURCE_DIR}/nimbus_pcg_sparse_mpi.cc)
SET(${PROJECT_NAME}_h
        ${PROJECT_SOURCE_DIR}/app.h
        ${PROJECT_SOURCE_DIR}/data_impl.h
        ${PROJECT_SOURCE_DIR}/job_impl.h
        ${PROJECT_SOURCE_DIR}/projection_driver.h
        ${PROJECT_SOURCE_DIR}/projection_example.h
        ${PROJECT_SOURCE_DIR}/projection_data.h
        ${PROJECT_SOURCE_DIR}/PCG_Sparse_Solver.cc
        ${PROJECT_SOURCE_DIR}/nimbus_pcg_sparse_mpi.h)
        
#List the libraries you want to link against - including a higher-level
#library automatically includes the lower-level ones
SET(PROTO_LIB protobuf)
SET(LIBRARIES_TO_LINK PhysBAM_Rendering-OpenGL)
LIST(APPEND LIBRARIES_TO_LINK ${PROTO_LIB})
        
#you can force cache variables to be on or off here
#SET(PHYSBAM_USE_CUDA ON CACHE BOOL ON FORCE)
        
#INCLUDE(${PROJECT_SOURCE_DIR}/../../Scripts/CMake/Project.cmake)
INCLUDE(${PROJECT_SOURCE_DIR}/../physbam-lib/Scripts/CMake/Nimbus_Application_Library.cmake)
