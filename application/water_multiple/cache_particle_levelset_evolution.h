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

#ifndef NIMBUS_APPLICATION_WATER_MULTIPLE_CACHE_PARTICLE_LEVELSET_EVOLUTION_H_
#define NIMBUS_APPLICATION_WATER_MULTIPLE_CACHE_PARTICLE_LEVELSET_EVOLUTION_H_

#include <string>

#include "application/water_multiple/physbam_include.h"
#include "data/cache/cache_struct.h"
#include "data/physbam/translator_physbam.h"
#include "shared/geometric_region.h"
#include "worker/data.h"

namespace application {

enum Particles {
    POS,
    NEG,
    POS_REM,
    NEG_REM
};

template<class TS = float>
class CacheParticleLevelsetEvolution : public nimbus::CacheStruct {
        typedef typename PhysBAM::VECTOR<TS, 3> TV;
        typedef typename PhysBAM::VECTOR<int, 3> TV_INT;
        typedef typename PhysBAM::RANGE<TV> Range;
        typedef typename PhysBAM::GRID<TV> Grid;
        typedef typename PhysBAM::PARTICLE_LEVELSET_EVOLUTION_UNIFORM<Grid> PhysBAMPLE;
        typedef typename PhysBAM::PARTICLE_LEVELSET_UNIFORM<Grid> PhysBAMParticleContainer;
        typedef typename PhysBAM::PARTICLE_LEVELSET_PARTICLES<TV> PhysBAMParticleBucket;
        typedef typename PhysBAM::PARTICLE_LEVELSET_EVOLUTION<Grid> PhysBAMParticleArray;
        typedef typename nimbus::TranslatorPhysBAM<TS> Translator;

    public:
        explicit CacheParticleLevelsetEvolution(
                const nimbus::GeometricRegion &global_reg,
                const nimbus::GeometricRegion &ob_reg =
                nimbus::GeometricRegion(),
                const int ghost_width = 0);

        PhysBAMPLE *data() {
            return data_;
        }
        void set_data(PhysBAMPLE *d) {
            data_ = d;
        }

    protected:
        virtual nimbus::CacheStruct *CreateNew(
                const nimbus::GeometricRegion &ob_reg) const;

        virtual void ReadToCache(
                const std::vector<nimbus::cache::type_id_t> &var_type,
                const std::vector<nimbus::DataArray> &read_sets,
                const nimbus::GeometricRegion &read_reg);
        virtual void WriteFromCache(
                const std::vector<nimbus::cache::type_id_t> &var_type,
                const std::vector<nimbus::DataArray> &write_sets,
                const nimbus::GeometricRegion &write_reg) const;

    private:
        nimbus::GeometricRegion global_region_;
        nimbus::GeometricRegion local_region_;
        nimbus::GeometricRegion wgb_region_;
        nimbus::GeometricRegion inner_region_;
        int ghost_width_;
        nimbus::Coord shift_;
        PhysBAMPLE *data_;
        Grid mac_grid_;
        int scale_;
}; // class CacheParticleLevelsetEvolution
} // namespace application

#endif // NIMBUS_APPLICATION_WATER_MULTIPLE_CACHE_PARTICLE_LEVELSET_EVOLUTION_H_
