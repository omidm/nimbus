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
#include <vector>

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
CacheParticleLevelsetEvolution(
        const nimbus::GeometricRegion &global_reg,
        int ghost_width,
        bool make_proto)
    : CacheStruct(NUM_PARTICLE_TYPES),
      global_region_(global_reg),
      ghost_width_(ghost_width) {
    if (make_proto)
        MakePrototype();
}

template<class TS> CacheParticleLevelsetEvolution<TS>::
CacheParticleLevelsetEvolution(
        const nimbus::GeometricRegion &global_reg,
        const nimbus::GeometricRegion &ob_reg,
        const int ghost_width)
    : CacheStruct(NUM_PARTICLE_TYPES, ob_reg),
      global_region_(global_reg),
      local_region_(ob_reg.NewEnlarged(-ghost_width)),
      inner_region_(ob_reg.NewEnlarged(-2*ghost_width)),
      ghost_width_(ghost_width) {
    {
        nimbus::int_dimension_t x = local_region_.x();
        nimbus::int_dimension_t y = local_region_.y();
        nimbus::int_dimension_t z = local_region_.z();
        nimbus::int_dimension_t dx = local_region_.dx();
        nimbus::int_dimension_t dy = local_region_.dy();
        nimbus::int_dimension_t dz = local_region_.dz();
        if (local_region_.x() == global_reg.x()) {
            x -= ghost_width;
            dx += ghost_width;
        }
        if (local_region_.x() + local_region_.dx() ==
            global_reg.x() + global_reg.dx()) {
            dx += ghost_width;
        }
        if (local_region_.y() == global_reg.y()) {
            y -= ghost_width;
            dy += ghost_width;
        }
        if (local_region_.y() + local_region_.dy() ==
            global_reg.y() + global_reg.dy()) {
            dy += ghost_width;
        }
        if (local_region_.z() == global_reg.z()) {
            z -= ghost_width;
            dz += ghost_width;
        }
        if (local_region_.z() + local_region_.dz() ==
            global_reg.z() + global_reg.dz()) {
            dz += ghost_width;
        }
        wgb_region_ = nimbus::GeometricRegion(x, y, z, dx, dy, dz);
    }
    shift_.x = local_region_.x() - global_reg.x();
    shift_.y = local_region_.y() - global_reg.y();
    shift_.z = local_region_.z() - global_reg.z();
    scale_ = global_reg.dx();
    if (local_region_.dx() > 0 && local_region_.dy() > 0 &&
        local_region_.dz() > 0) {
        Range domain = RangeFromRegions<TV>(global_reg, local_region_);
        TV_INT count = CountFromRegion(local_region_);
        mac_grid_.Initialize(count, domain, true);
        data_ = new PhysBAMPLE(mac_grid_, ghost_width);
        {
            data_->grid = mac_grid_;
            PhysBAMParticleContainer *particle_levelset =
                &data_->particle_levelset;
            particle_levelset->Set_Band_Width(6);
            // Resize phi
            data_->phi.Resize(mac_grid_.
                    Domain_Indices(particle_levelset->number_of_ghost_cells));
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
            particle_levelset->Set_Minimum_Particle_Radius(
                (TS).1*particle_levelset->levelset.grid.Minimum_Edge_Length());
            particle_levelset->Set_Maximum_Particle_Radius(
                (TS).5*particle_levelset->levelset.grid.Minimum_Edge_Length());
            if (particle_levelset->half_band_width &&
                particle_levelset->levelset.grid.Minimum_Edge_Length()) {
              particle_levelset->Set_Band_Width(particle_levelset->half_band_width /
                  ((TS).5*particle_levelset->levelset.grid.Minimum_Edge_Length()));
            } else {
              particle_levelset->Set_Band_Width();
            }
            particle_levelset->levelset.Initialize_Levelset_Grid_Values();
            if (data_->levelset_advection.semi_lagrangian_collidable) {
              particle_levelset->levelset.Initialize_Valid_Masks(mac_grid_);
            }
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

template<class TS> nimbus::CacheStruct *CacheParticleLevelsetEvolution<TS>::
CreateNew(const nimbus::GeometricRegion &ob_reg) const {
    return new CacheParticleLevelsetEvolution(
                              global_region_,
                              ob_reg,
                              ghost_width_);
}

template<class TS> void CacheParticleLevelsetEvolution<TS>::
ReadToCache(
        const std::vector<nimbus::cache::type_id_t> &var_type,
        const std::vector<nimbus::DataArray> &read_sets,
        const nimbus::GeometricRegion &read_reg) {
    // dbg(DBG_WARN, "\n--- Reading %i elements into particles for region %s\n", read_set.size(), reg.toString().c_str());
    bool merge = true;
    PhysBAMParticleContainer *particle_levelset = &data_->particle_levelset;
    for (size_t t = 0; t < var_type.size(); ++t) {
        const nimbus::DataArray &read_set = read_sets[t];
        std::vector<nimbus::GeometricRegion> regs;
        for (size_t i = 0; i < read_set.size(); ++i) {
            regs.push_back(read_set[i]->region());
        }
        switch (var_type[t]) {
            case POS:
                Translator::DeleteParticles(
                        shift_, regs, particle_levelset, scale_, true);
                Translator::ReadParticles(
                        object_region(), shift_, read_set,
                        particle_levelset, scale_, true, merge);
                break;
            case NEG:
                Translator::DeleteParticles(
                        shift_, regs, particle_levelset, scale_, false);
                Translator::ReadParticles(
                        object_region(), shift_, read_set,
                        particle_levelset, scale_, false, merge);
                break;
            case POS_REM:
                Translator::DeleteRemovedParticles(
                        shift_, regs, particle_levelset, scale_, true);
                Translator::ReadRemovedParticles(
                        object_region(), shift_, read_set,
                        particle_levelset, scale_, true, merge);
                break;
            case NEG_REM:
                Translator::DeleteRemovedParticles(
                        shift_, regs, particle_levelset, scale_, false);
                Translator::ReadRemovedParticles(
                        object_region(), shift_, read_set,
                        particle_levelset, scale_, false, merge);
                break;
            default:
                dbg(DBG_ERROR, "Unexpected particle type in ReadToCache function\n");
                exit(-1);
        }
    }
}

template<class TS> void CacheParticleLevelsetEvolution<TS>::
WriteFromCache(
        const std::vector<nimbus::cache::type_id_t> &var_type,
        const std::vector<nimbus::DataArray> &write_sets,
        const nimbus::GeometricRegion &write_reg) const {
    // dbg(DBG_WARN, "\n--- Writing %i elements into particles for region %s\n", write_set.size(), reg.toString().c_str());
    PhysBAMParticleContainer *particle_levelset = &data_->particle_levelset;
    for (size_t t = 0; t < var_type.size(); ++t) {
        const nimbus::DataArray &write_set = write_sets[t];
        switch (var_type[t]) {
            case POS:
                Translator::WriteParticles(
                        object_region(), shift_, write_set,
                        particle_levelset, scale_, true);
                break;
            case NEG:
                Translator::WriteParticles(
                        object_region(), shift_, write_set,
                        particle_levelset, scale_, false);
                break;
            case POS_REM:
                Translator::WriteRemovedParticles(
                        object_region(), shift_, write_set,
                        particle_levelset, scale_, true);
                break;
            case NEG_REM:
                Translator::WriteRemovedParticles(
                        object_region(), shift_, write_set,
                        particle_levelset, scale_, false);
                break;
            default:
                dbg(DBG_ERROR, "Unexpected particle type in WriteFromCache function\n");
                exit(-1);
        }
    }
}

template class CacheParticleLevelsetEvolution<float>;

} // namespace application
