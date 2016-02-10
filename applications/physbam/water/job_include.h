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
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#ifndef NIMBUS_APPLICATION_WATER_MULTIPLE_JOB_INCLUDE_H_
#define NIMBUS_APPLICATION_WATER_MULTIPLE_JOB_INCLUDE_H_

#include "applications/physbam/water//job_adjust_phi.h"
#include "applications/physbam/water//job_update_ghost_velocities.h"
#include "applications/physbam/water//job_advect_phi.h"
#include "applications/physbam/water//job_advect_removed_particles.h"
#include "applications/physbam/water//job_advect_v.h"
#include "applications/physbam/water//job_apply_forces.h"
#include "applications/physbam/water//job_calculate_dt.h"
#include "applications/physbam/water//job_delete_particles.h"
#include "applications/physbam/water//job_extrapolate_phi.h"
#include "applications/physbam/water//job_extrapolation.h"
#include "applications/physbam/water//job_initialize.h"
#include "applications/physbam/water//job_loop_frame.h"
#include "applications/physbam/water//job_loop_iteration.h"
#include "applications/physbam/water//job_loop_iteration_part_two.h"
#include "applications/physbam/water//job_main.h"
#include "applications/physbam/water//job_modify_levelset_part_one.h"
#include "applications/physbam/water//job_modify_levelset_part_two.h"
#include "applications/physbam/water//job_make_signed_distance.h"
#include "applications/physbam/water//job_barrier.h"
#include "applications/physbam/water//job_names.h"
#include "applications/physbam/water//projection/job_projection_main.h"
#include "applications/physbam/water//projection/job_projection_transform_pressure.h"
#include "applications/physbam/water//projection/job_projection_calculate_boundary_condition_part_one.h"
#include "applications/physbam/water//projection/job_projection_calculate_boundary_condition_part_two.h"
#include "applications/physbam/water//projection/job_projection_construct_matrix.h"
#include "applications/physbam/water//projection/job_projection_wrapup.h"
#include "applications/physbam/water//job_reincorporate_removed_particles.h"
#include "applications/physbam/water//job_step_particles.h"
#include "applications/physbam/water//job_synchronize_particles.h"
#include "applications/physbam/water//job_write_output.h"
#include "applications/physbam/water//job_reseed_particles.h"

#include "applications/physbam/water///projection/job_projection_global_initialize.h"
#include "applications/physbam/water///projection/job_projection_local_initialize.h"
#include "applications/physbam/water///projection/job_projection_loop_iteration.h"
#include "applications/physbam/water///projection/job_projection_loop_bottleneck.h"
#include "applications/physbam/water//projection/job_projection_step_one.h"
#include "applications/physbam/water//projection/job_projection_reduce_rho.h"
#include "applications/physbam/water//projection/job_projection_step_two.h"
#include "applications/physbam/water//projection/job_projection_step_three.h"
#include "applications/physbam/water//projection/job_projection_reduce_alpha.h"
#include "applications/physbam/water//projection/job_projection_step_four.h"

#endif  // NIMBUS_APPLICATION_WATER_MULTIPLE_JOB_INCLUDE_H_
