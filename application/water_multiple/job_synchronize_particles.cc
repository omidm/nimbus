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

#include "application/water_multiple/app_utils.h"
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
        typedef application::DataConfig DataConfig;
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

    nimbus::DataArray readp, readn, readpr, readnr;
    nimbus::DataArray write;

    const std::string pp_string = std::string(APP_POS_PARTICLES);
    GetReadData(*this, pp_string, da, &readp);
    GetWriteData(*this, pp_string, da, &write, false);
    const std::string np_string = std::string(APP_NEG_PARTICLES);
    GetReadData(*this, np_string, da, &readn);
    GetWriteData(*this, np_string, da, &write, false);
    const std::string prp_string = std::string(APP_POS_REM_PARTICLES);
    GetReadData(*this, prp_string, da, &readpr);
    GetWriteData(*this, prp_string, da, &write, false);
    const std::string nrp_string = std::string(APP_NEG_REM_PARTICLES);
    GetReadData(*this, nrp_string, da, &readnr);
    GetWriteData(*this, nrp_string, da, &write, false);

    // read_inner. Simply delete particles outside.
    nimbus::DataArray read_inner, read_outer,
        read_inner_p, read_inner_n,
        read_inner_pr, read_inner_nr,
        read_outer_p, read_outer_n,
        read_outer_pr, read_outer_nr;

    for (size_t i = 0; i < readp.size(); ++i) {
        nimbus::Data *d = readp[i];
        nimbus::GeometricRegion dr = d->region();
        if (array_inner.Covers(&dr)) {
            read_inner.push_back(d);
            read_inner_p.push_back(d);
        } else {
            read_outer.push_back(d);
            read_outer_p.push_back(d);;
        }
    }
    for (size_t i = 0; i < readn.size(); ++i) {
        nimbus::Data *d = readn[i];
        nimbus::GeometricRegion dr = d->region();
        if (array_inner.Covers(&dr)) {
            read_inner.push_back(d);
            read_inner_n.push_back(d);
    } else {
            read_outer.push_back(d);
            read_outer_n.push_back(d);
        }
    }
    for (size_t i = 0; i < readpr.size(); ++i) {
        nimbus::Data *d = readpr[i];
        nimbus::GeometricRegion dr = d->region();
        if (array_inner.Covers(&dr)) {
            read_inner.push_back(d);
            read_inner_pr.push_back(d);
    } else {
            read_outer.push_back(d);
            read_outer_pr.push_back(d);
        }
    }
    for (size_t i = 0; i < readnr.size(); ++i) {
        nimbus::Data *d = readnr[i];
        nimbus::GeometricRegion dr = d->region();
        if (array_inner.Covers(&dr)) {
            read_inner.push_back(d);
            read_inner_nr.push_back(d);
    } else {
            read_outer.push_back(d);
            read_outer_nr.push_back(d);
        }
    }

    nimbus::CacheManager *cm = GetCacheManager();
    dbg(DBG_WARN, "\n--- Sync particles for region %s, requesting read %i, write %i\n",
        init_config.local_region.toString().c_str(), read_inner.size(), write.size());
    nimbus::CacheObject *cache_obj =
      cm->GetAppObject(read_inner, write,
          array_outer,
          application::kCachePLE,
          nimbus::EXCLUSIVE, false, true);

    CacheParticleLevelsetEvolution<T> *cache_ple =
        dynamic_cast<CacheParticleLevelsetEvolution<T> *>(cache_obj);
    assert(cache_ple != NULL);

    PhysBAM::PARTICLE_LEVELSET_EVOLUTION_UNIFORM<PhysBAM::GRID<TV> >
        *particle_levelset_evolution = cache_ple->data();
    PhysBAM::PARTICLE_LEVELSET_UNIFORM<PhysBAM::GRID<TV> >
        *particle_levelset = &particle_levelset_evolution->particle_levelset;

    typedef typename nimbus::TranslatorPhysBAM<T> Translator;
    nimbus::GeometricRegion local_region = init_config.local_region;
    nimbus::GeometricRegion global_region = init_config.global_region;
    nimbus::GeometricRegion enlarge(1-kGhostNum,
                                    1-kGhostNum,
                                    1-kGhostNum,
                                    local_region.dx()+2*kGhostNum,
                                    local_region.dy()+2*kGhostNum,
                                    local_region.dz()+2*kGhostNum);
    nimbus::Coord shift;
    shift.x = local_region.x() - global_region.x();
    shift.y = local_region.y() - global_region.y();
    shift.z = local_region.z() - global_region.z();
    int scale = global_region.dx();

    dbg(DBG_WARN, "--- inner_p %i, inner_n %i, outer_p %i, outer_n %i\n",
            read_inner_p.size(), read_inner_n.size(),
            read_outer_p.size(), read_outer_n.size());

    // TODO(Chinmayee): get rid of all inner statements once delete is
    // implemented
    Translator::ReadParticles(enlarge, shift, read_inner_p, particle_levelset, scale, true);
    Translator::ReadParticles(enlarge, shift, read_inner_n, particle_levelset, scale, false);
    Translator::ReadRemovedParticles(enlarge, shift, read_inner_pr, particle_levelset, scale, true);
    Translator::ReadRemovedParticles(enlarge, shift, read_inner_nr, particle_levelset, scale, false);
    Translator::ReadParticles(enlarge, shift, read_outer_p, particle_levelset, scale, true, true);
    Translator::ReadParticles(enlarge, shift, read_outer_n, particle_levelset, scale, false, true);
    Translator::ReadRemovedParticles(enlarge, shift, read_outer_pr, particle_levelset, scale, true, true);
    Translator::ReadRemovedParticles(enlarge, shift, read_outer_nr, particle_levelset, scale, false, true);

    cache_ple->Write(array_outer, true);

    dbg(APP_LOG, "Finish translating particles.\n");

    dbg(APP_LOG, "Completed executing synchronize particles job\n");
}

} // namespace application
