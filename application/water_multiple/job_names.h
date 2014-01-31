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
 * Author: Omid Mashayekhi <omidm@stanford.edu>
 */

#ifndef NIMBUS_APPLICATION_WATER_ALTERNATE_FINE_JOB_NAMES_H_
#define NIMBUS_APPLICATION_WATER_ALTERNARE_FINE_JOB_NAMES_H_


#define SUPER_1 "super_1"
#define SUPER_2 "super_2"
#define SUPER_3 "super_3"

#define MAIN "main"
#define INITIALIZE "initialize"
#define LOOP_FRAME "loop_frame"
#define LOOP_ITERATION "loop_iteration"
#define CALCULATE_FRAME "calculate_frame"
#define WRITE_FRAME "write_frame"

#define ADJUST_PHI_WITH_OBJECTS "adjust_phi_with_objects"
#define ADVECT_PHI "advect_phi"
#define STEP_PARTICLES "step_particles"
#define ADVECT_REMOVED_PARTICLES "advect_removed_particles"
#define ADVECT_V "advect_v"
#define APPLY_FORCES "apply_forces"
#define MODIFY_LEVELSET "modify_levelset"
#define ADJUST_PHI "adjust_phi"
#define DELETE_PARTICLES "delete_particles"
#define REINCORPORATE_PARTICLES "reincorporate_particles"
#define PROJECTION "projection"
#define EXTRAPOLATION "extrapolation"


#endif  // NIMBUS_APPLICATION_WATER_ALTERNATE_FINE_JOB_NAMES_H_
