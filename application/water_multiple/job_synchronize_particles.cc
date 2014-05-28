/* Copyright 2013 Stanford University.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 vd* are met:
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

#include <string>
#include <vector>

#include "application/water_multiple/app_utils.h"
#include "application/water_multiple/cache_particle_levelset_evolution.h"
#include "application/water_multiple/data_particle_array.h"
#include "application/water_multiple/job_synchronize_particles.h"
#include "application/water_multiple/job_names.h"
#include "application/water_multiple/options.h"
#include "application/water_multiple/parameters.h"
#include "application/water_multiple/physbam_include.h"
#include "application/water_multiple/physbam_utils.h"
#include "application/water_multiple/water_driver.h"
#include "application/water_multiple/water_example.h"
#include "data/physbam/translator_physbam.h"
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "shared/nimbus.h"

namespace application {

JobSynchronizeParticles::JobSynchronizeParticles(nimbus::Application *app) {
    set_application(app);
};

nimbus::Job* JobSynchronizeParticles::Clone() {
    return new JobSynchronizeParticles(application());
}

void JobSynchronizeParticles::Execute(nimbus::Parameter params, const nimbus::DataArray& da) {
    dbg(APP_LOG, "--- Executing synchronize particles job\n");

    // get time, dt, frame from the parameters.
    InitConfig init_config;
    init_config.use_cache = true;
    T dt;
    std::string params_str(params.ser_data().data_ptr_raw(),
                           params.ser_data().size());
    LoadParameter(params_str, &init_config.frame, &init_config.time, &dt,
                  &init_config.global_region, &init_config.local_region);
    dbg(APP_LOG, " Loaded parameters (Frame=%d, Time=%f, dt=%f).\n",
        init_config.frame, init_config.time, dt);


    if (!kUseCache) {
        // initializing the example and driver with state and configuration variables
        PhysBAM::WATER_EXAMPLE<TV> *example;
        PhysBAM::WATER_DRIVER<TV> *driver;
        DataConfig data_config;
        data_config.SetFlag(DataConfig::POSITIVE_PARTICLE);
        data_config.SetFlag(DataConfig::NEGATIVE_PARTICLE);
        data_config.SetFlag(DataConfig::REMOVED_POSITIVE_PARTICLE);
        data_config.SetFlag(DataConfig::REMOVED_NEGATIVE_PARTICLE);
        InitializeExampleAndDriver(init_config, data_config, this,
                                   da, example, driver);

        example->Save_To_Nimbus(this, da, init_config.frame + 1);
        DestroyExampleAndDriver(example, driver);
        return;
    }

    nimbus::GeometricRegion array_inner(
            init_config.local_region.NewEnlarged(-kGhostNum));
    nimbus::GeometricRegion array_outer(
            init_config.local_region.NewEnlarged(kGhostNum));

    typedef std::vector<nimbus::GeometricRegion> GRegArray;
    nimbus::cache::type_id_t vars[] = {POS, NEG, POS_REM, NEG_REM};
    std::vector<nimbus::cache::type_id_t> var_type(
            vars, vars + sizeof(vars)/sizeof(nimbus::cache::type_id_t));
    std::vector<nimbus::DataArray> read(4), write(4), read_inner(4), read_outer(4);
    std::vector<GRegArray> reg(4);
    std::string dtype[] = {  APP_POS_PARTICLES,
                             APP_NEG_PARTICLES,
                             APP_POS_REM_PARTICLES,
                             APP_NEG_REM_PARTICLES
                          };
    for (size_t t = 0; t < NUM_PARTICLE_TYPES; ++t) {
        GetReadData(*this, dtype[t], da, &read[t], false);
        GetWriteData(*this, dtype[t], da, &write[t], false);
        nimbus::DataArray &read_t = read[t];
        nimbus::DataArray &read_inner_t = read_inner[t];
        nimbus::DataArray &read_outer_t = read_outer[t];
        GRegArray &reg_t = reg[t];
        for (size_t i = 0; i < read_t.size(); ++i) {
            nimbus::Data *d = read_t[i];
            nimbus::GeometricRegion dreg = d->region();
            if (array_inner.Covers(&dreg)) {
                read_inner_t.push_back(d);
            } else{
                read_outer_t.push_back(d);
                reg_t.push_back(dreg);
            }
        }
    }

    nimbus::CacheManager *cm = GetCacheManager();
    nimbus::CacheStruct *cache_struct =
      cm->GetAppStruct(
              var_type,
              read_inner, array_outer,
              write, array_outer,
              kCachePLE, array_outer,
              nimbus::cache::EXCLUSIVE);
    CacheParticleLevelsetEvolution<T> *cache_ple =
        dynamic_cast<CacheParticleLevelsetEvolution<T> *>(cache_struct);
    assert(cache_ple != NULL);

    PhysBAM::PARTICLE_LEVELSET_EVOLUTION_UNIFORM<PhysBAM::GRID<TV> >
        *particle_levelset_evolution = cache_ple->data();
    PhysBAM::PARTICLE_LEVELSET_UNIFORM<PhysBAM::GRID<TV> >
        *particle_levelset = &particle_levelset_evolution->particle_levelset;

    typedef typename nimbus::TranslatorPhysBAM<T> Translator;
    nimbus::GeometricRegion local_region = init_config.local_region;
    nimbus::GeometricRegion global_region = init_config.global_region;
    nimbus::Coord shift;
    shift.x = local_region.x() - global_region.x();
    shift.y = local_region.y() - global_region.y();
    shift.z = local_region.z() - global_region.z();
    int scale = global_region.dx();

    Translator::DeleteParticles(shift, reg[POS], particle_levelset, scale, true);
    Translator::DeleteParticles(shift, reg[NEG], particle_levelset, scale, false);
    Translator::DeleteRemovedParticles(shift, reg[POS_REM], particle_levelset, scale, true);
    Translator::DeleteRemovedParticles(shift, reg[NEG_REM], particle_levelset, scale, false);
    Translator::ReadParticles(array_outer, shift, read_outer[POS], particle_levelset, scale, true, true);
    Translator::ReadParticles(array_outer, shift, read_outer[NEG], particle_levelset, scale, false, true);
    Translator::ReadRemovedParticles(array_outer, shift, read_outer[POS_REM], particle_levelset, scale, true, true);
    Translator::ReadRemovedParticles(array_outer, shift, read_outer[NEG_REM], particle_levelset, scale, false, true);

    cache_ple->ReleaseAccess();

    dbg(APP_LOG, "Completed executing synchronize particles job\n");
}

} // namespace application
