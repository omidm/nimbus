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
 * The class provides functions that the application writer can use to get
 * scratch data for a job, scratch data from a region for synchronization
 * job and scratch name registration utilities.
 *
 * The user should provide base nimbus type name or a list of scratch names.
 *
 * Supports only 3d.
 *
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#include <string>
#include <vector>

#include "data/scratch_data_helper.h"
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "shared/nimbus_types.h"
#include "worker/application.h"
#include "worker/data.h"
#include "worker/job.h"

namespace nimbus {

ScratchDataHelper::ScratchDataHelper() {}

ScratchDataHelper::ScratchDataHelper(const int gw[DIMENSION]) {
    set_ghost_width(gw);
}

ScratchDataHelper::ScratchDataHelper(const int gw[DIMENSION],
                                     const std::string b_name) {
    set_ghost_width(gw);
    SetScratchBaseName(b_name);
}

ScratchDataHelper::~ScratchDataHelper() {}

void ScratchDataHelper::set_ghost_width(const int gw[DIMENSION]) {
    for (size_t i = 0; i < DIMENSION; i++)
        ghost_width_[i] = (gw[i] >= 0)? gw[i] : 0;
}

void ScratchDataHelper::SetScratchNames(const ScratchType st,
                                        const std::vector<std::string> &st_names) {
    size_t nl = st_names.size();
    switch (st) {
        case VERTEX:
            if (nl != VERTEX_TYPES)
                dbg(DBG_WARN, "WARNING: scratch type name vector for VERTEX does not match required size of %i\n", VERTEX_TYPES); // NOLINT
            nl = (VERTEX_TYPES < nl)? VERTEX_TYPES : nl;
            for (size_t i = 0; i < nl; i++)
                vertex_types_[i] = st_names[i];
            break;
        case EDGE:
            if (nl != EDGE_TYPES)
                dbg(DBG_WARN, "WARNING: scratch type name vector for EDGE does not match required size of %i\n", EDGE_TYPES); // NOLINT
            nl = (EDGE_TYPES < nl)? EDGE_TYPES : nl;
            for (size_t i = 0; i < nl; i++)
                edge_types_[i] = st_names[i];
            break;
        case FACE:
            if (nl != FACE_TYPES)
                dbg(DBG_WARN, "WARNING: scratch type name vector for FACE does not match required size of %i\n", FACE_TYPES); // NOLINT
            nl = (FACE_TYPES < nl)? FACE_TYPES : nl;
            for (size_t i = 0; i < nl; i++)
                face_types_[i] = st_names[i];
            break;
        default:
            dbg(DBG_WARN, "WARNING: invalid scratch type %i, ignoring it in SetScratchType\n"); // NOLINT
            return;
    }
}

void ScratchDataHelper::SetScratchBaseName(const std::string b_name) {
    base_name_ = b_name;
    for (size_t i = 0; i < VERTEX_TYPES; i++) {
        std::stringstream ss;
        ss << b_name << "_vertex_" << i+1;
        vertex_types_[i] = ss.str();
    }
    for (size_t i = 0; i < EDGE_TYPES; i++) {
        std::stringstream ss;
        ss << b_name << "_edge_" << i+1;
        edge_types_[i] = ss.str();
    }
    for (size_t i = 0; i < FACE_TYPES; i++) {
        std::stringstream ss;
        ss << b_name << "_face_" << i+1;
        face_types_[i] = ss.str();
    }
}

void ScratchDataHelper::RegisterScratchNames(Application *app,
                                             Data *data) const {
    for (size_t i = 0; i < VERTEX_TYPES; i++) {
        app->RegisterData(vertex_types_[i], data);
    }
    for (size_t i = 0; i < EDGE_TYPES; i++) {
        app->RegisterData(edge_types_[i], data);
    }
    for (size_t i = 0; i < FACE_TYPES; i++) {
        app->RegisterData(face_types_[i], data);
    }
}

void ScratchDataHelper::GetJobScratchData(Job *job,
                                          const GeometricRegion &cr,
                                          lIDSet *ids,
                                          bool clear) const {
    if (clear)
        ids->clear();

    const int cl[DIMENSION]  = {cr.x(),  cr.y(),  cr.z()};
    int cld[DIMENSION] = {cr.dx(), cr.dy(), cr.dx()};
    int l[DIMENSION]  = {0, 0, 0};
    int ld[DIMENSION] = {0, 0, 0};
    size_t n;

    // vertex scratch regions
    for (size_t d = 0; d < DIMENSION; d++) {
        ld[d] = 2*ghost_width_[d];
    }
    n  = 0;
    for (size_t i = 0; i < 2; i++) {
        l[XCOORD] = cl[XCOORD] - ghost_width_[XCOORD] +
                    i * cld[XCOORD];
        for (size_t j = 0; j < 2; j++) {
            l[YCOORD] = cl[YCOORD] - ghost_width_[YCOORD] +
                        j * cld[YCOORD];
            for (size_t k = 0; k < 2; k++) {
                l[ZCOORD] = cl[ZCOORD] - ghost_width_[ZCOORD] +
                            k * cld[ZCOORD];
                CLdoVector ldos;
                GeometricRegion region(l[XCOORD], l[YCOORD], l[ZCOORD],
                                       ld[XCOORD], ld[YCOORD], ld[ZCOORD]);
                job->GetCoveredLogicalObjects(&ldos, vertex_types_[n], &region);
                dbg(DBG_WARN, "Getting data %s for region %s\n",
                        vertex_types_[n].c_str(),
                        region.toString().c_str());
                for (size_t s = 0; s < ldos.size(); s++)
                    ids->insert(ldos[s]->id());
                n++;
            }
        }
    }

    // edge scratch regions
    for (size_t d = 0; d < DIMENSION; d++) {
        l[d]  = cl[d] + ghost_width_[d];
        ld[d] = cld[d] - 2*ghost_width_[d];
        size_t d1 = (d+1)%DIMENSION;
        size_t d2 = (d+2)%DIMENSION;
        ld[d1] = 2*ghost_width_[d1];
        ld[d2] = 2*ghost_width_[d2];
        n = 0;
        for (int i = 0; i < 2; i++) {
            l[d1] = cl[d1] - ghost_width_[d1] +
                    i * cld[d1];
            for (int j = 0; j < 2; j++) {
                l[d2] = cl[d2] - ghost_width_[d2] +
                        j * cld[d2];
                CLdoVector ldos;
                GeometricRegion region(l[XCOORD], l[YCOORD], l[ZCOORD],
                                       ld[XCOORD], ld[YCOORD], ld[ZCOORD]);
                job->GetCoveredLogicalObjects(&ldos, edge_types_[n], &region);
                for (size_t s = 0; s < ldos.size(); s++)
                    ids->insert(ldos[s]->id());
                n++;
            }
        }
    }

    // face scratch regions
    for (size_t d = 0; d < DIMENSION; d++) {
        size_t d1 = (d+1)%DIMENSION;
        size_t d2 = (d+2)%DIMENSION;
        ld[d] = 2*ghost_width_[d];
        l[d1]  = cl[d1] + ghost_width_[d1];
        l[d2]  = cl[d2] + ghost_width_[d2];
        ld[d1] = cld[d1] - 2*ghost_width_[d1];
        ld[d2] = cld[d2] - 2*ghost_width_[d2];
        n = 0;
        for (size_t i = 0; i < 2; i++) {
            l[d]  = cl[d] - ghost_width_[d] + i * cld[d];
            CLdoVector ldos;
            GeometricRegion region(l[XCOORD], l[YCOORD], l[ZCOORD],
                                   ld[XCOORD], ld[YCOORD], ld[ZCOORD]);
            job->GetCoveredLogicalObjects(&ldos, face_types_[n], &region);
            for (size_t s = 0; s < ldos.size(); s++)
                ids->insert(ldos[s]->id());
            n++;
        }
    }
}

void ScratchDataHelper::GetMainForScratchData(Job *job,
                                              const GeometricRegion &complete,
                                              const GeometricRegion &exclude,
                                              lIDSet *ids,
                                              bool clear) const {
    if (clear)
        ids->clear();

    CLdoVector ldos_complete, ldos_exclude;
    job->GetCoveredLogicalObjects(&ldos_complete, base_name_, &complete);
    job->GetCoveredLogicalObjects(&ldos_exclude, base_name_, &exclude);
    for (size_t s = 0; s < ldos_complete.size(); s++) {
        ids->insert(ldos_complete[s]->id());
    }
    for (size_t s = 0; s < ldos_exclude.size(); s++) {
        ids->remove(ldos_exclude[s]->id());
    }
}

void ScratchDataHelper::GetAllScratchData(Job *job,
                                          const GeometricRegion &region,
                                          lIDSet *ids,
                                          bool clear) const {
    if (clear)
        ids->clear();

    {
        CLdoVector ldos_v;
        for (size_t i = 0; i < VERTEX_TYPES; i++) {
            job->GetCoveredLogicalObjects(&ldos_v, vertex_types_[i], &region);
            for (size_t s = 0; s < ldos_v.size(); s++)
                ids->insert(ldos_v[s]->id());
        }
    }
    {
        CLdoVector ldos_e;
        for (size_t i = 0; i < EDGE_TYPES; i++) {
            job->GetCoveredLogicalObjects(&ldos_e, edge_types_[i], &region);
            for (size_t s = 0; s < ldos_e.size(); s++)
                ids->insert(ldos_e[s]->id());
        }
    }
    {
        CLdoVector ldos_f;
        for (size_t i = 0; i < FACE_TYPES; i++) {
            job->GetCoveredLogicalObjects(&ldos_f, face_types_[i], &region);
            for (size_t s = 0; s < ldos_f.size(); s++)
                ids->insert(ldos_f[s]->id());
        }
    }
}

void ScratchDataHelper::GetAllScratchData(Job *job,
                                          const GeometricRegion &region,
                                          ScratchType st,
                                          lIDSet *ids,
                                          bool clear) const {
    if (clear)
        ids->clear();

    CLdoVector ldos;
    switch (st) {
        case VERTEX:
            for (size_t i = 0; i < VERTEX_TYPES; i++)
                job->GetCoveredLogicalObjects(&ldos, vertex_types_[i], &region);
            break;
        case EDGE:
            for (size_t i = 0; i < EDGE_TYPES; i++)
                job->GetCoveredLogicalObjects(&ldos, edge_types_[i], &region);
            break;
        case FACE:
            for (size_t i = 0; i < FACE_TYPES; i++)
                job->GetCoveredLogicalObjects(&ldos, face_types_[i], &region);
            break;
        default:
            dbg(DBG_WARN, "WARNING: invalid scratch type %i, ignoring it in GetAllScratchData\n", st); // NOLINT
            return;
    }
    for (size_t s = 0; s < ldos.size(); s++)
        ids->insert(ldos[s]->id());
}
}  // namespace nimbus
