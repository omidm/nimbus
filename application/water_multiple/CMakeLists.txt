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
    ${PROJECT_SOURCE_DIR}/data_def.cc
    ${PROJECT_SOURCE_DIR}/data_face_array.cc
    ${PROJECT_SOURCE_DIR}/data_particle_array.cc
    ${PROJECT_SOURCE_DIR}/data_scalar_array.cc
    ${PROJECT_SOURCE_DIR}/job_calculate_frame.cc
    ${PROJECT_SOURCE_DIR}/job_initialize.cc
    ${PROJECT_SOURCE_DIR}/job_loop_frame.cc
    ${PROJECT_SOURCE_DIR}/job_loop_iteration.cc
    ${PROJECT_SOURCE_DIR}/job_adjust_phi.cc
    ${PROJECT_SOURCE_DIR}/job_delete_particles.cc
    ${PROJECT_SOURCE_DIR}/job_main.cc
    ${PROJECT_SOURCE_DIR}/job_modify_levelset_part_one.cc
    ${PROJECT_SOURCE_DIR}/job_modify_levelset_part_two.cc
    ${PROJECT_SOURCE_DIR}/job_reincorporate_removed_particles.cc
    ${PROJECT_SOURCE_DIR}/job_update_ghost_velocities.cc
    ${PROJECT_SOURCE_DIR}/job_extrapolate_phi.cc
    ${PROJECT_SOURCE_DIR}/job_advect_phi.cc
    ${PROJECT_SOURCE_DIR}/job_step_particles.cc
    ${PROJECT_SOURCE_DIR}/job_advect_removed_particles.cc
    ${PROJECT_SOURCE_DIR}/job_advect_v.cc
    ${PROJECT_SOURCE_DIR}/job_apply_forces.cc
    ${PROJECT_SOURCE_DIR}/job_extrapolation.cc
    ${PROJECT_SOURCE_DIR}/job_synchronize_particles.cc
    ${PROJECT_SOURCE_DIR}/job_write_frame.cc
    ${PROJECT_SOURCE_DIR}/physbam_utils.cc
    ${PROJECT_SOURCE_DIR}/reg_def.cc
    ${PROJECT_SOURCE_DIR}/water_app.cc
    ${PROJECT_SOURCE_DIR}/water_driver.cc
    ${PROJECT_SOURCE_DIR}/water_example.cc
    ${PROJECT_SOURCE_DIR}/water_sources.cc
    ${PROJECT_SOURCE_DIR}/projection/data_raw_grid_array.cc
    ${PROJECT_SOURCE_DIR}/projection/data_raw_array_m2c.cc
    ${PROJECT_SOURCE_DIR}/projection/data_raw_vector_nd.cc
    ${PROJECT_SOURCE_DIR}/projection/data_sparse_matrix.cc
    ${PROJECT_SOURCE_DIR}/projection/job_projection_calculate_boundary_condition.cc
    ${PROJECT_SOURCE_DIR}/projection/job_projection_construct_matrix.cc
    ${PROJECT_SOURCE_DIR}/projection/job_projection_core.cc
    ${PROJECT_SOURCE_DIR}/projection/job_projection_wrapup.cc
    ${PROJECT_SOURCE_DIR}/projection/laplace_solver_wrapper.cc
    ${PROJECT_SOURCE_DIR}/projection/projection_driver.cc
    ${PROJECT_SOURCE_DIR}/projection/projection_helper.cc
    )
SET(${PROJECT_NAME}_h
    ${PROJECT_SOURCE_DIR}/app_utils.h
    ${PROJECT_SOURCE_DIR}/data_include.h
    ${PROJECT_SOURCE_DIR}/data_names.h
    ${PROJECT_SOURCE_DIR}/data_def.h
    ${PROJECT_SOURCE_DIR}/data_face_array.h
    ${PROJECT_SOURCE_DIR}/data_particle_array.h
    ${PROJECT_SOURCE_DIR}/data_scalar_array.h
    ${PROJECT_SOURCE_DIR}/job_include.h
    ${PROJECT_SOURCE_DIR}/job_names.h
    ${PROJECT_SOURCE_DIR}/job_calculate_frame.h
    ${PROJECT_SOURCE_DIR}/job_initialize.h
    ${PROJECT_SOURCE_DIR}/job_loop_frame.h
    ${PROJECT_SOURCE_DIR}/job_loop_iteration.h
    ${PROJECT_SOURCE_DIR}/job_adjust_phi.h
    ${PROJECT_SOURCE_DIR}/job_delete_particles.h
    ${PROJECT_SOURCE_DIR}/job_main.h
    ${PROJECT_SOURCE_DIR}/job_modify_levelset_part_one.h
    ${PROJECT_SOURCE_DIR}/job_modify_levelset_part_two.h
    ${PROJECT_SOURCE_DIR}/job_reincorporate_removed_particles.h
    ${PROJECT_SOURCE_DIR}/job_update_ghost_velocities.h
    ${PROJECT_SOURCE_DIR}/job_extrapolate_phi.h
    ${PROJECT_SOURCE_DIR}/job_advect_phi.h
    ${PROJECT_SOURCE_DIR}/job_step_particles.h
    ${PROJECT_SOURCE_DIR}/job_advect_removed_particles.h
    ${PROJECT_SOURCE_DIR}/job_advect_v.h
    ${PROJECT_SOURCE_DIR}/job_apply_forces.h
    ${PROJECT_SOURCE_DIR}/job_extrapolation.h
    ${PROJECT_SOURCE_DIR}/job_synchronize_particles.h
    ${PROJECT_SOURCE_DIR}/job_write_frame.h
    ${PROJECT_SOURCE_DIR}/options.h
    ${PROJECT_SOURCE_DIR}/physbam_include.h
    ${PROJECT_SOURCE_DIR}/physbam_utils.h
    ${PROJECT_SOURCE_DIR}/reg_def.h
    ${PROJECT_SOURCE_DIR}/water_app.h
    ${PROJECT_SOURCE_DIR}/water_driver.h
    ${PROJECT_SOURCE_DIR}/water_example.h
    ${PROJECT_SOURCE_DIR}/water_sources.h
    ${PROJECT_SOURCE_DIR}/projection/data_raw_grid_array.h
    ${PROJECT_SOURCE_DIR}/projection/data_raw_array_m2c.h
    ${PROJECT_SOURCE_DIR}/projection/data_raw_vector_nd.h
    ${PROJECT_SOURCE_DIR}/projection/data_sparse_matrix.h
    ${PROJECT_SOURCE_DIR}/projection/job_projection_calculate_boundary_condition.h
    ${PROJECT_SOURCE_DIR}/projection/job_projection_construct_matrix.h
    ${PROJECT_SOURCE_DIR}/projection/job_projection_core.h
    ${PROJECT_SOURCE_DIR}/projection/job_projection_wrapup.h
    ${PROJECT_SOURCE_DIR}/projection/laplace_solver_wrapper.h
    ${PROJECT_SOURCE_DIR}/projection/translator_util.h
    ${PROJECT_SOURCE_DIR}/projection/projection_driver.h
    ${PROJECT_SOURCE_DIR}/projection/projection_helper.h
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
