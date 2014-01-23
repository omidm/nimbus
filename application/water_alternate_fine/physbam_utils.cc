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

#include "shared/nimbus.h"
#include "application/water_alternate_fine/app_utils.h"
#include "application/water_alternate_fine/water_driver.h"
#include "application/water_alternate_fine/water_example.h"
#include "application/water_alternate_fine/water_sources.h"

#include "application/water_alternate_fine/physbam_utils.h"

namespace application {

bool InitializeExampleAndDriver(
    const InitConfig& init_config,
    const nimbus::Job* job,
    const nimbus::DataArray& da,
    PhysBAM::WATER_EXAMPLE<TV>*& example,
    PhysBAM::WATER_DRIVER<TV>*& driver) {
  example = new PhysBAM::WATER_EXAMPLE<TV>(PhysBAM::STREAM_TYPE((RW())));
  example->Initialize_Grid(
      init_config.grid_size,
      PhysBAM::RANGE<TV>(TV(), TV::All_Ones_Vector()));
  PhysBAM::WaterSources::Add_Source(example);
  driver= new PhysBAM::WATER_DRIVER<TV>(*example);
  driver->init_phase = init_config.init_phase;
  driver->current_frame = init_config.frame;
  driver->time = init_config.time;
  driver->Initialize(job, da, init_config.set_boundary_condition);
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

template<class TV>
bool InitializeParticleLevelsetEvolutionHelper(
    const DataConfig& data_config,
    const PhysBAM::GRID<TV>& grid_input,
    PhysBAM::PARTICLE_LEVELSET_EVOLUTION_UNIFORM<PhysBAM::GRID<TV> >*
    particle_levelset_evolution) {
  typedef float T;
  PhysBAM::PARTICLE_LEVELSET_UNIFORM<PhysBAM::GRID<TV> >* particle_levelset =
      &particle_levelset_evolution->particle_levelset;
  assert(grid_input.Is_MAC_Grid());
  particle_levelset_evolution->grid = grid_input;
  // Resizes phi here.
  if (data_config.GetFlag(DataConfig::LEVELSET)) {
    particle_levelset_evolution->phi.Resize(
        grid_input.Domain_Indices(particle_levelset->number_of_ghost_cells));
  }
  // Resizes particles.
  if (data_config.GetFlag(DataConfig::POSITIVE_PARTICLE)) {
    particle_levelset->positive_particles.Resize(
        particle_levelset->levelset.grid.Block_Indices(
            particle_levelset->number_of_ghost_cells));
  }
  if (data_config.GetFlag(DataConfig::NEGATIVE_PARTICLE)) {
    particle_levelset->negative_particles.Resize(
        particle_levelset->levelset.grid.Block_Indices(
            particle_levelset->number_of_ghost_cells));
  }
  particle_levelset->use_removed_positive_particles=true;
  particle_levelset->use_removed_negative_particles=true;
  // Resizes removed particles.
  if (data_config.GetFlag(DataConfig::REMOVED_POSITIVE_PARTICLE)) {
    particle_levelset->removed_positive_particles.Resize(
        particle_levelset->levelset.grid.Block_Indices(
            particle_levelset->number_of_ghost_cells));
  }
  if (data_config.GetFlag(DataConfig::REMOVED_NEGATIVE_PARTICLE)) {
    particle_levelset->removed_negative_particles.Resize(
        particle_levelset->levelset.grid.Block_Indices(
            particle_levelset->number_of_ghost_cells));
  }

  particle_levelset->Set_Minimum_Particle_Radius(
      (T).1*particle_levelset->levelset.grid.Minimum_Edge_Length());
  particle_levelset->Set_Maximum_Particle_Radius(
      (T).5*particle_levelset->levelset.grid.Minimum_Edge_Length());
  if (particle_levelset->half_band_width &&
      particle_levelset->levelset.grid.Minimum_Edge_Length()) {
   particle_levelset->Set_Band_Width(particle_levelset->half_band_width /
                   ((T).5*particle_levelset->levelset.grid.Minimum_Edge_Length()));
  } else {
    particle_levelset->Set_Band_Width();
  }
  particle_levelset->levelset.Initialize_Levelset_Grid_Values();
  if (particle_levelset_evolution->
      levelset_advection.semi_lagrangian_collidable) {
    particle_levelset->levelset.Initialize_Valid_Masks(grid_input);
  }
  return true;
}

template
bool InitializeParticleLevelsetEvolutionHelper<PhysBAM::VECTOR<float,2> >(
    const DataConfig& data_config,
    const PhysBAM::GRID<PhysBAM::VECTOR<float,2> >& grid_input,
    PhysBAM::PARTICLE_LEVELSET_EVOLUTION_UNIFORM
    <PhysBAM::GRID<PhysBAM::VECTOR<float,2> > >* particle_levelset_evolution);

template
bool InitializeParticleLevelsetEvolutionHelper<PhysBAM::VECTOR<float,3> >(
    const DataConfig& data_config,
    const PhysBAM::GRID<PhysBAM::VECTOR<float,3> >& grid_input,
    PhysBAM::PARTICLE_LEVELSET_EVOLUTION_UNIFORM
    <PhysBAM::GRID<PhysBAM::VECTOR<float,3> > >* particle_levelset_evolution);

}  // namespace application
