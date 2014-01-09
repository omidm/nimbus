/* Copyright 2013 Stanford University.
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
 * This file contains super job 2 that : modifies phi with particles, adjusts
 * levelset with sources, deletes particles and reincorporates removed
 * particles.
 *
 * This job uses all data -- face velocities, phi, particles and removed
 * particles!
 *
 * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
 */

#include "application/water_alternate_fine/app_utils.h"
#include "application/water_alternate_fine/job_super_2.h"
#include "application/water_alternate_fine/physbam_include.h"
#include "application/water_alternate_fine/water_driver.h"
#include "application/water_alternate_fine/water_example.h"
#include "application/water_alternate_fine/water_sources.h"
#include "data/physbam/physbam_data.h"
#include "shared/dbg.h"
#include "shared/nimbus.h"
#include <sstream>
#include <string>

namespace application {

    JobSuper2::JobSuper2(nimbus::Application *app) {
        set_application(app);
    };

    nimbus::Job* JobSuper2::Clone() {
        return new JobSuper2(application());
    }

    void JobSuper2::Execute(nimbus::Parameter params, const nimbus::DataArray& da) {
        dbg(APP_LOG, "Executing 2nd super job\n");

        int frame;
        T time, dt;
        std::string params_str(params.ser_data().data_ptr_raw(),
                               params.ser_data().size());
        LoadParameter(params_str, &frame, &time, &dt);

        dbg(APP_LOG, "Frame %i in super job 2\n", frame);

        // initialize configuration and state
        PhysBAM::WATER_EXAMPLE<TV> *example =
            new PhysBAM::WATER_EXAMPLE<TV>(PhysBAM::STREAM_TYPE((RW())));

        example->Initialize_Grid(TV_INT::All_Ones_Vector()*kScale,
                                 PhysBAM::RANGE<TV>(TV(),
                                                    TV::All_Ones_Vector())
                                 );
        PhysBAM::WaterSources::Add_Source(example);
        PhysBAM::WATER_DRIVER<TV> driver(*example);
        driver.init_phase = false;
        driver.current_frame = frame;
        driver.Initialize(this, da, false);

        // face velocity for ghost + interior
        FaceArray face_velocities_ghost;
        face_velocities_ghost.Resize(example->incompressible.grid,
                                     example->number_of_ghost_cells,
                                     false);
        example->incompressible.boundary->
            Fill_Ghost_Cells_Face(example->mac_grid,
                                  example->face_velocities,
                                  face_velocities_ghost,
                                  time+dt,
                                  example->number_of_ghost_cells);

        // modify levelset
        dbg(APP_LOG, "Modify Levelset ...\n");
        example->particle_levelset_evolution.particle_levelset.
            Exchange_Overlap_Particles();
        example->particle_levelset_evolution.
            Modify_Levelset_And_Particles(&face_velocities_ghost);

        // adjust phi with sources
        dbg(APP_LOG, "Adjust Phi ...\n");
        example->Adjust_Phi_With_Sources(time+dt);

        // delete particles
        dbg(APP_LOG, "Delete Particles ...\n");
        example->particle_levelset_evolution.Delete_Particles_Outside_Grid();
        example->particle_levelset_evolution.particle_levelset.
            Delete_Particles_In_Local_Maximum_Phi_Cells(1);
        example->particle_levelset_evolution.particle_levelset.
            Delete_Particles_Far_From_Interface(); // uses visibility
        example->particle_levelset_evolution.particle_levelset.
            Identify_And_Remove_Escaped_Particles(face_velocities_ghost,
                                                  1.5,
                                                  time + dt);

        // reincorporate removed particles
        dbg(APP_LOG, "Reincorporate Removed Particles ...\n");
        if (example->particle_levelset_evolution.particle_levelset.
                use_removed_positive_particles ||
            example->particle_levelset_evolution.particle_levelset.
                use_removed_negative_particles)
            example->particle_levelset_evolution.particle_levelset.
                Reincorporate_Removed_Particles(1, 1, 0, true);

        // free resources
        delete example;

        dbg(APP_LOG, "Completed executing 2nd super job job\n");
    }

} // namespace application
