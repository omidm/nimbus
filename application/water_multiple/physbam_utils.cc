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

#include "application/water_multiple/app_utils.h"
#include "application/water_multiple/water_driver.h"
#include "application/water_multiple/water_example.h"
#include "application/water_multiple/water_sources.h"
#include "shared/nimbus.h"

#include "application/water_multiple_fine/physbam_utils.h"

namespace application {

Range GridToRange(
    const GeometricRegion& global_region,
    const GeometricRegion& local_region) {
  TV start, end;
  start(1) = (local_region.x() - 1) / global_region.dx();
  start(2) = (local_region.y() - 1) / global_region.dy();
  start(3) = (local_region.z() - 1) / global_region.dz();
  end(1) =  (local_region.x() + local_region.dx() - 1) / global_region.dx();
  end(2) =  (local_region.y() + local_region.dy() - 1) / global_region.dy();
  end(3) =  (local_region.z() + local_region.dz() - 1) / global_region.dz();
  return Range(start, end);
}

bool InitializeExampleAndDriver(
    const InitConfig& init_config,
    const DataConfig& data_config,
    const nimbus::Job* job,
    const nimbus::DataArray& da,
    PhysBAM::WATER_EXAMPLE<TV>*& example,
    PhysBAM::WATER_DRIVER<TV>*& driver) {
  dbg(APP_LOG, "HANG:Enter initialize_example_driver.\n");
  example = new PhysBAM::WATER_EXAMPLE<TV>(PhysBAM::STREAM_TYPE((RW())));
  example->Initialize_Grid(
      TV_INT(init_config.local_region.dx(),
             init_config.local_region.dy(),
             init_config.local_region.dz()),
      GridToRange(init_config.global_region, init_config.local_region));
  PhysBAM::WaterSources::Add_Source(example);
  example->data_config.Set(data_config);
  driver= new PhysBAM::WATER_DRIVER<TV>(*example);
  driver->init_phase = init_config.init_phase;
  driver->current_frame = init_config.frame;
  driver->time = init_config.time;
  dbg(APP_LOG, "Before enter driver->Initialize.\n");
  driver->Initialize(job, da, init_config.set_boundary_condition);
  dbg(APP_LOG, "HANG:Exit initialize_example_driver.\n");
  return true;
}

void DestroyExampleAndDriver(
    PhysBAM::WATER_EXAMPLE<TV>*& example,
    PhysBAM::WATER_DRIVER<TV>*& driver) {
  delete example;
  example = NULL;
  delete driver;
  driver = NULL;
}

}  // namespace application
