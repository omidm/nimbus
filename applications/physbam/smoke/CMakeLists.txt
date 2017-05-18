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
#SET(PROJECT_NAME water_app)
SET(PROJECT_NAME smoke_app)
PROJECT(${PROJECT_NAME})

#Set Nimbus root and prtobuf files directory
SET(NIMBUS_ROOT ${PROJECT_SOURCE_DIR}/../../..)

# Add the include path for nimbus, boost, leveldb, and protobuf.
INCLUDE(${NIMBUS_ROOT}/Makeinclude.cmake)
add_definitions(${IFLAGS})

include_directories(${NIMBUS_ROOT})
LINK_DIRECTORIES(${NIMBUS_ROOT})

#Manually list your cpp and h files here - do NOT use GLOB/GLOB_RECURSE
SET(${PROJECT_NAME}_cpp
    ${PROJECT_SOURCE_DIR}/app_utils.cc
    ${PROJECT_SOURCE_DIR}/app_data_face_array.cc
    ${PROJECT_SOURCE_DIR}/app_data_prototypes.cc
    ${PROJECT_SOURCE_DIR}/app_data_scalar_array.cc
    ${PROJECT_SOURCE_DIR}/app_data_compressed_scalar_array.cc
    ${PROJECT_SOURCE_DIR}/data_face_array.cc
    ${PROJECT_SOURCE_DIR}/data_scalar_array.cc
    ${PROJECT_SOURCE_DIR}/data_compressed_scalar_array.cc
    ${PROJECT_SOURCE_DIR}/job_initialize.cc
    ${PROJECT_SOURCE_DIR}/job_loop_frame.cc
    ${PROJECT_SOURCE_DIR}/job_loop_iteration.cc
    ${PROJECT_SOURCE_DIR}/job_loop_iteration_part_two.cc
    ${PROJECT_SOURCE_DIR}/job_main.cc
    ${PROJECT_SOURCE_DIR}/job_update_ghost_velocities.cc
    ${PROJECT_SOURCE_DIR}/job_update_ghost_densities.cc
    ${PROJECT_SOURCE_DIR}/job_write_output.cc
    ${PROJECT_SOURCE_DIR}/job_substep.cc
    ${PROJECT_SOURCE_DIR}/job_scalar_advance.cc
    ${PROJECT_SOURCE_DIR}/job_convect.cc
    ${PROJECT_SOURCE_DIR}/nimbus_thread_queue.cc
    ${PROJECT_SOURCE_DIR}/physbam_tools.cc
    ${PROJECT_SOURCE_DIR}/physbam_utils.cc
    ${PROJECT_SOURCE_DIR}/smoke_app.cc
    ${PROJECT_SOURCE_DIR}/smoke_driver.cc
    ${PROJECT_SOURCE_DIR}/smoke_driver_initializers.cc
    ${PROJECT_SOURCE_DIR}/smoke_example.cc
    ${PROJECT_SOURCE_DIR}/projection/app_data_array_m2c.cc
    ${PROJECT_SOURCE_DIR}/projection/app_data_sparse_matrix.cc
    ${PROJECT_SOURCE_DIR}/projection/app_data_vector.cc
    ${PROJECT_SOURCE_DIR}/projection/app_data_raw_grid_array.cc
    ${PROJECT_SOURCE_DIR}/projection/data_raw_grid_array.cc
    ${PROJECT_SOURCE_DIR}/projection/data_raw_array_m2c.cc
    ${PROJECT_SOURCE_DIR}/projection/data_raw_vector_nd.cc
    ${PROJECT_SOURCE_DIR}/projection/data_sparse_matrix.cc
    ${PROJECT_SOURCE_DIR}/projection/job_projection_main.cc
    ${PROJECT_SOURCE_DIR}/projection/job_projection_transform_pressure.cc
    ${PROJECT_SOURCE_DIR}/projection/job_projection_calculate_boundary_condition_part_one.cc
    ${PROJECT_SOURCE_DIR}/projection/job_projection_calculate_boundary_condition_part_two.cc
    ${PROJECT_SOURCE_DIR}/projection/job_projection_construct_matrix.cc
    ${PROJECT_SOURCE_DIR}/projection/job_projection_wrapup.cc
    ${PROJECT_SOURCE_DIR}/projection/laplace_solver_wrapper.cc
    ${PROJECT_SOURCE_DIR}/projection/projection_driver.cc
    ${PROJECT_SOURCE_DIR}/projection/projection_driver_caching.cc
    ${PROJECT_SOURCE_DIR}/projection/projection_helper.cc
    ${PROJECT_SOURCE_DIR}/projection/job_projection_local_initialize.cc
    ${PROJECT_SOURCE_DIR}/projection/job_projection_global_initialize.cc
    ${PROJECT_SOURCE_DIR}/projection/job_projection_loop_iteration.cc
    ${PROJECT_SOURCE_DIR}/projection/job_projection_step_one.cc
    ${PROJECT_SOURCE_DIR}/projection/job_projection_reduce_rho.cc
    ${PROJECT_SOURCE_DIR}/projection/job_projection_step_two.cc
    ${PROJECT_SOURCE_DIR}/projection/job_projection_step_three.cc
    ${PROJECT_SOURCE_DIR}/projection/job_projection_reduce_alpha.cc
    ${PROJECT_SOURCE_DIR}/projection/job_projection_step_four.cc
    )
SET(${PROJECT_NAME}_h
    ${PROJECT_SOURCE_DIR}/app_utils.h
    ${PROJECT_SOURCE_DIR}/app_data_include.h
    ${PROJECT_SOURCE_DIR}/app_data_face_array.h
    ${PROJECT_SOURCE_DIR}/app_data_options.h
    ${PROJECT_SOURCE_DIR}/app_data_prototypes.h
    ${PROJECT_SOURCE_DIR}/app_data_scalar_array.h
    ${PROJECT_SOURCE_DIR}/data_include.h
    ${PROJECT_SOURCE_DIR}/data_names.h
    ${PROJECT_SOURCE_DIR}/data_face_array.h
    ${PROJECT_SOURCE_DIR}/data_scalar_array.h
    ${PROJECT_SOURCE_DIR}/app_data_compressed_scalar_array.h
    ${PROJECT_SOURCE_DIR}/data_compressed_scalar_array.h
    ${PROJECT_SOURCE_DIR}/job_include.h
    ${PROJECT_SOURCE_DIR}/job_names.h
    ${PROJECT_SOURCE_DIR}/job_initialize.h
    ${PROJECT_SOURCE_DIR}/job_loop_frame.h
    ${PROJECT_SOURCE_DIR}/job_loop_iteration.h
    ${PROJECT_SOURCE_DIR}/job_loop_iteration_part_two.h
    ${PROJECT_SOURCE_DIR}/job_main.h
    ${PROJECT_SOURCE_DIR}/job_update_ghost_velocities.h
    ${PROJECT_SOURCE_DIR}/job_update_ghost_densities.h
    ${PROJECT_SOURCE_DIR}/job_write_output.h
    ${PROJECT_SOURCE_DIR}/job_substep.h
    ${PROJECT_SOURCE_DIR}/job_scalar_advance.h
    ${PROJECT_SOURCE_DIR}/job_convect.h
    ${PROJECT_SOURCE_DIR}/nimbus_thread_queue.h
    ${PROJECT_SOURCE_DIR}/options.h
    ${PROJECT_SOURCE_DIR}/parameters.h
    ${PROJECT_SOURCE_DIR}/physbam_include.h
    ${PROJECT_SOURCE_DIR}/physbam_tools.h
    ${PROJECT_SOURCE_DIR}/physbam_utils.h
    ${PROJECT_SOURCE_DIR}/smoke_app.h
    ${PROJECT_SOURCE_DIR}/smoke_driver.h
    ${PROJECT_SOURCE_DIR}/smoke_example.h
    ${PROJECT_SOURCE_DIR}/projection/app_data_array_m2c.h
    ${PROJECT_SOURCE_DIR}/projection/app_data_sparse_matrix.h
    ${PROJECT_SOURCE_DIR}/projection/app_data_vector.h
    ${PROJECT_SOURCE_DIR}/projection/app_data_raw_grid_array.h
    ${PROJECT_SOURCE_DIR}/projection/data_raw_grid_array.h
    ${PROJECT_SOURCE_DIR}/projection/data_raw_array_m2c.h
    ${PROJECT_SOURCE_DIR}/projection/data_raw_vector_nd.h
    ${PROJECT_SOURCE_DIR}/projection/data_sparse_matrix.h
    ${PROJECT_SOURCE_DIR}/projection/job_projection_main.h
    ${PROJECT_SOURCE_DIR}/projection/job_projection_transform_pressure.h
    ${PROJECT_SOURCE_DIR}/projection/job_projection_calculate_boundary_condition_part_one.h
    ${PROJECT_SOURCE_DIR}/projection/job_projection_calculate_boundary_condition_part_two.h
    ${PROJECT_SOURCE_DIR}/projection/job_projection_construct_matrix.h
    ${PROJECT_SOURCE_DIR}/projection/job_projection_wrapup.h
    ${PROJECT_SOURCE_DIR}/projection/laplace_solver_wrapper.h
    ${PROJECT_SOURCE_DIR}/projection/translator_util.h
    ${PROJECT_SOURCE_DIR}/projection/projection_driver.h
    ${PROJECT_SOURCE_DIR}/projection/projection_helper.h
    ${PROJECT_SOURCE_DIR}/projection/job_projection_local_initialize.h
    ${PROJECT_SOURCE_DIR}/projection/job_projection_global_initialize.h
    ${PROJECT_SOURCE_DIR}/projection/job_projection_loop_iteration.h
    ${PROJECT_SOURCE_DIR}/projection/job_projection_step_one.h
    ${PROJECT_SOURCE_DIR}/projection/job_projection_reduce_rho.h
    ${PROJECT_SOURCE_DIR}/projection/job_projection_step_two.h
    ${PROJECT_SOURCE_DIR}/projection/job_projection_step_three.h
    ${PROJECT_SOURCE_DIR}/projection/job_projection_reduce_alpha.h
    ${PROJECT_SOURCE_DIR}/projection/job_projection_step_four.h
    )

#List the libraries you want to link against - including a higher-level
#library automatically includes the lower-level ones
SET(LIBRARIES_TO_LINK PhysBAM_Rendering-OpenGL)
# The rest of the libraries, will be linked by the executable that uses the
# shared library generated here. -omidm
# SET(PROTO_LIB protobuf)
# SET(NIMBUS_LIB nimbus)
# SET(BOOST_THREAD boost_thread)
# SET(BOOST_SYSTEM boost_system)
# LIST(APPEND LIBRARIES_TO_LINK ${PROTO_LIB})
# LIST(APPEND LIBRARIES_TO_LINK ${NIMBUS_LIB})
# LIST(APPEND LIBRARIES_TO_LINK ${BOOST_THREAD})
# LIST(APPEND LIBRARIES_TO_LINK ${BOOST_SYSTEM})

#you can force cache variables to be on or off here
#SET(PHYSBAM_USE_CUDA ON CACHE BOOL ON FORCE)

INCLUDE(${NIMBUS_ROOT}/applications/physbam/physbam-lib/Scripts/CMake/Nimbus_Application_Library.cmake)
