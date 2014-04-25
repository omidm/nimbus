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
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#include <string>

#include "application/water_multiple/cache_particle_levelset_evolution.h"
#include "application/water_multiple/data_names.h"
#include "application/water_multiple/physbam_include.h"
#include "application/water_multiple/physbam_tools.h"
#include "data/cache/cache_object.h"
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "worker/data.h"

namespace application {
template<class TS> CacheParticleLevelsetEvolution<TS>::
CacheParticleLevelsetEvolution(std::string type,
               const nimbus::GeometricRegion &global_region,
               const int ghost_width,
               const nimbus::GeometricRegion &app_region)
    : CacheObject(type, app_region),
      ghost_width_(ghost_width),
      global_region_(global_region),
      local_region_(app_region.NewEnlarged(-ghost_width_)) {
    shift_.x = local_region_.x() - global_region.x();
    shift_.y = local_region_.y() - global_region.y();
    shift_.z = local_region_.z() - global_region.z();
    enlarge_ =  nimbus::GeometricRegion(1-ghost_width_,
                                        1-ghost_width_,
                                        1-ghost_width_,
                                        local_region_.dx()+2*ghost_width_,
                                        local_region_.dy()+2*ghost_width_,
                                        local_region_.dz()+2*ghost_width_);
    scale_ = global_region_.dx();
    if (local_region_.dx() > 0 && local_region_.dy() > 0 && local_region_.dz() > 0) {
        Range domain = RangeFromRegions<TV>(global_region, local_region_);
        TV_INT count = CountFromRegion(local_region_);
        mac_grid.Initialize(count, domain, true);
        data_ = new PhysBAMPLE(mac_grid, ghost_width);
        {
            data_->grid = mac_grid;
            PhysBAMParticleContainer *particle_levelset =
                &data_->particle_levelset;
            particle_levelset->Set_Band_Width(6);
            // Resizes particles.
            particle_levelset->positive_particles.Resize(
                    particle_levelset->levelset.grid.Block_Indices(
                        particle_levelset->number_of_ghost_cells));
            particle_levelset->negative_particles.Resize(
                    particle_levelset->levelset.grid.Block_Indices(
                        particle_levelset->number_of_ghost_cells));
            particle_levelset->use_removed_positive_particles=true;
            particle_levelset->use_removed_negative_particles=true;
            // Resizes removed particles.
            particle_levelset->removed_positive_particles.Resize(
                    particle_levelset->levelset.grid.Block_Indices(
                        particle_levelset->number_of_ghost_cells));
            particle_levelset->removed_negative_particles.Resize(
                    particle_levelset->levelset.grid.Block_Indices(
                        particle_levelset->number_of_ghost_cells));
            }
        {
            // policies etc
            data_->Set_CFL_Number((TS).9);
            data_->Set_Number_Particles_Per_Cell(16);
            data_->Initialize_FMM_Initialization_Iterative_Solver(true);
            data_->Bias_Towards_Negative_Particles(false);
            data_->particle_levelset.Use_Removed_Positive_Particles();
            data_->particle_levelset.Use_Removed_Negative_Particles();
            data_->particle_levelset.Store_Unique_Particle_Id();
            data_->Use_Particle_Levelset(true);
            data_->particle_levelset.Set_Collision_Distance_Factors(.1,1);
        }
    }
}

template<class TS> void CacheParticleLevelsetEvolution<TS>::
ReadDiffToCache(const nimbus::DataArray &read_set,
                const nimbus::DataArray &diff,
                const nimbus::GeometricRegion &reg) {
    dbg(DBG_WARN, "\n--- Reading %i elements into particles for region %s\n", read_set.size(), reg.toString().c_str());
    bool read_outer_only = true;
    for (size_t i = 0; i < diff.size(); ++i) {
        nimbus::Data *d = diff[i];
        nimbus::GeometricRegion dr = d->region();
        if (local_region_.Covers(&dr)) {
            read_outer_only = false;
            break;
        }
    }
    nimbus::DataArray final_read;
    if (read_outer_only) {
        for (size_t i = 0; i < read_set.size(); ++i)
            final_read.push_back(read_set[i]);
    } else {
        final_read = read_set;
    }
    nimbus::DataArray pos, neg, pos_rem, neg_rem;
    for (size_t i = 0; i < final_read.size(); ++i) {
        nimbus::Data *d = final_read[i];
        if (d->name() == APP_POS_PARTICLES) {
            pos.push_back(d);
        } else if (d->name() == APP_NEG_PARTICLES) {
            neg.push_back(d);
        } else if (d->name() == APP_POS_REM_PARTICLES) {
            pos_rem.push_back(d);
        } else if (d->name() == APP_NEG_REM_PARTICLES) {
            neg_rem.push_back(d);
        }
    }
    // TODO(Chinmayee): something faster??
    PhysBAMParticleContainer *particle_levelset = &data_->particle_levelset;
    if (read_outer_only) {
        data_->Delete_Particles_Outside_Grid();
        dbg(DBG_WARN, "\n--- Finally reading %i elements into particles for region %s\n", final_read.size(), reg.toString().c_str());
    }
    Translator::ReadParticles(enlarge_, shift_, pos, particle_levelset, scale_, true, read_outer_only);
    Translator::ReadParticles(enlarge_, shift_, neg, particle_levelset, scale_, false, read_outer_only);
    Translator::ReadRemovedParticles(enlarge_, shift_, pos_rem, particle_levelset, scale_, true, read_outer_only);
    Translator::ReadRemovedParticles(enlarge_, shift_, neg_rem, particle_levelset, scale_, false, read_outer_only);
}

template<class TS> void CacheParticleLevelsetEvolution<TS>::
WriteFromCache(const nimbus::DataArray &write_set,
               const nimbus::GeometricRegion &reg) const {
    dbg(DBG_WARN, "\n Writing %i elements into particles for region %s\n", write_set.size(), reg.toString().c_str());
    nimbus::DataArray pos, neg, pos_rem, neg_rem;
    for (size_t i = 0; i < write_set.size(); ++i) {
        nimbus::Data *d = write_set[i];
        if (d->name() == APP_POS_PARTICLES) {
            pos.push_back(d);
        } else if (d->name() == APP_NEG_PARTICLES) {
            neg.push_back(d);
        } else if (d->name() == APP_POS_REM_PARTICLES) {
            pos_rem.push_back(d);
        } else if (d->name() == APP_NEG_REM_PARTICLES) {
            neg_rem.push_back(d);
        }
    }
    PhysBAMParticleContainer *particle_levelset = &data_->particle_levelset;
    Translator::WriteParticles(enlarge_, shift_, pos, particle_levelset, scale_, true);
    Translator::WriteParticles(enlarge_, shift_, neg, particle_levelset, scale_, false);
    Translator::WriteRemovedParticles(enlarge_, shift_, pos_rem, particle_levelset, scale_, true);
    Translator::WriteRemovedParticles(enlarge_, shift_, neg_rem, particle_levelset, scale_, false);
}

template<class TS> nimbus::CacheObject *CacheParticleLevelsetEvolution<TS>::
CreateNew(const nimbus::GeometricRegion &ar) const {
    return new CacheParticleLevelsetEvolution(type(),
                              global_region_,
                              ghost_width_,
                              ar);
}

template class CacheParticleLevelsetEvolution<float>;

} // namespace application
