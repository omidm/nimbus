/*
 * Copyright 2013 Stanford University.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the
 *   distribution.
 *
 * - Neither the name of the copyright holders nor the names of
 *   its contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
 * THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * This file defines the name of jobs that will be used for registration and
 * spawning the jobs.
 *
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#ifndef NIMBUS_APPLICATION_WATER_MULTIPLE_DATA_NAMES_H_
#define NIMBUS_APPLICATION_WATER_MULTIPLE_DATA_NAMES_H_

#define APP_FACE_VEL "face_vel"
#define APP_FACE_VEL_GHOST "face_vel_ghost"
#define APP_PHI "phi"
#define APP_POS_PARTICLES "pos_particles"
#define APP_NEG_PARTICLES "neg_particles"
#define APP_POS_REM_PARTICLES "pos_rem_particles"
#define APP_NEG_REM_PARTICLES "neg_rem_particles"
#define APP_LAST_UNIQUE_PARTICLE_ID "last_unique_particle_id"

#define APP_PSI_D "psi_d"
#define APP_PSI_N "psi_n"
#define APP_PRESSURE "pressure"
#define APP_FILLED_REGION_COLORS "filled_region_colors"
#define APP_DIVERGENCE "divergence"
#define APP_U_INTERFACE "u_interface"

#define APP_MATRIX_A "matrix_a"
#define APP_VECTOR_B "vector_b"
#define APP_INDEX_C2M "index_c2m"
#define APP_INDEX_M2C "index_m2c"
#define APP_PROJECTION_LOCAL_N "projection_local_n"
#define APP_PROJECTION_INTERIOR_N "projection_interior_n"

#define APP_PROJECTION_LOCAL_TOLERANCE "projection_local_tolerance"
#define APP_PROJECTION_GLOBAL_TOLERANCE "projection_global_tolerance"
#define APP_PROJECTION_GLOBAL_N "projection_global_n"
#define APP_PROJECTION_DESIRED_ITERATIONS "projection_desired_iterations"

#define APP_PROJECTION_LOCAL_RESIDUAL "projection_local_residual"
#define APP_PROJECTION_LOCAL_RHO "projection_local_rho"
#define APP_PROJECTION_GLOBAL_RHO "projection_global_rho"
#define APP_PROJECTION_GLOBAL_RHO_OLD "projection_global_rho_old"
#define APP_PROJECTION_LOCAL_DOT_PRODUCT_FOR_ALPHA "projection_local_dot_product_for_alpha"
#define APP_PROJECTION_ALPHA "projection_alpha"
#define APP_PROJECTION_BETA "projection_beta"
#define APP_MATRIX_C "matrix_c"
#define APP_VECTOR_Z "vector_z"
#define APP_VECTOR_P "vector_p"
#define APP_VECTOR_TEMP "vector_temp"

#endif  // NIMBUS_APPLICATION_WATER_MULTIPLE_DATA_NAMES_H_
