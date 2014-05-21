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
 * A CacheTable contains a map from geometric region to application cache
 * objects.
 *
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#include <vector>

#include "data/cache/cache_defs.h"
#include "data/cache/cache_object.h"
#include "data/cache/cache_struct.h"
#include "data/cache/cache_table.h"
#include "data/cache/cache_var.h"
#include "shared/geometric_region.h"
#include "worker/data.h"

namespace nimbus {

/**
 * \details If ttype equals cache::VAR, tstruct_ is set to NULL. Similarly, is ttype
 * equals cache::STRUCT, tvar_ is set to NULL.
 */
CacheTable::CacheTable(cache::CacheTType ttype) {
    if (ttype == cache::VAR) {
        tvar_ = new TVar();
        tstruct_ = NULL;
    } else {
        tvar_ = NULL;
        tstruct_ = new TStruct();
    }
}

/**
 * \details If there is already a mapping for the given region key,
 * AddEntry(...) appends cv to the vector of cache vars that region maps to.
 * Otherwise AddEntry(...) creates a new mapping from region key to a vector of
 * cache vars with one element, cv.
 */
void CacheTable::AddEntry(const GeometricRegion &region,
                          CacheVar *cv) {
    CacheVars *cvs = NULL;
    if (tvar_->find(region) == tvar_->end()) {
        cvs = new CacheVars();
        (*tvar_)[region] = cvs;
    } else {
        cvs = tvar_->at(region);
    }
    cvs->push_back(cv);
}

/**
 * \details If there is already a mapping for the given region key,
 * AddEntry(...) appends cs to the vector of cache structs that region maps to.
 * Otherwise AddEntry(...) creates a new mapping from region key to a vector of
 * cache structs with one element, cs.
 */
void CacheTable::AddEntry(const GeometricRegion &region,
                          CacheStruct *cs) {
    CacheStructs *css = NULL;
    if (tstruct_->find(region) == tstruct_->end()) {
        css = new CacheStructs();
        (*tstruct_)[region] = css;
    } else {
        css = tstruct_->at(region);
    }
    css->push_back(cs);
}

/**
 * \details GetClosestAvailable(...) returns NULL if no usable cache var
 * instance is available, for the given region and access mode. Otherwise it
 * finds the closest cache var instance using the CacheVar cost function.
 */
CacheVar *CacheTable::GetClosestAvailable(const GeometricRegion &region,
                                          const DataArray &read_set,
                                          cache::CacheAccess access) {
    if (tvar_->find(region) == tvar_->end())
        return NULL;
    int min = GetMinDistanceIndex(*tvar_->at(region), read_set, access);
    if (min == -1)
        return NULL;
    else
        return tvar_->at(region)->at(min);
}

/**
 * \details GetClosestAvailable(...) returns NULL if no usable cache struct
 * instance is available, for the given region and access mode. Otherwise it
 * finds the closest cache struct instance using the CacheStruct cost function.
 */
CacheStruct *CacheTable::GetClosestAvailable(const GeometricRegion &region,
                                             const std::vector<cache::type_id_t> &var_type,
                                             const std::vector<DataArray> &read_sets,
                                             cache::CacheAccess access) {
    if (tstruct_->find(region) == tstruct_->end())
        return NULL;
    int min = GetMinDistanceIndex(*tstruct_->at(region), var_type, read_sets, access);
    if (min == -1)
        return NULL;
    else
        return tstruct_->at(region)->at(min);
}

/**
 * \details
 */
int CacheTable::GetMinDistanceIndex(const CacheVars &cvs,
                                    const DataArray &read_set,
                                    cache::CacheAccess access) const {
    int min_index = -1;
    if (cvs.size() == 0 || !cvs.at(0)->IsAvailable(access))
        return min_index;
    cache::distance_t min_distance = cvs.at(0)->GetDistance(read_set);
    cache::distance_t dist;
    for (size_t i = 1; i < cvs.size(); ++i) {
        if (!cvs.at(i)->IsAvailable(access))
            continue;
        dist = cvs.at(i)->GetDistance(read_set);
        if (dist == 0) {
            min_distance = 0;
            min_index = i;
            break;
        } else if (dist < min_distance) {
            min_distance = dist;
            min_index = i;
        }
    }
    return min_index;
}

/**
 * \details
 */
int CacheTable::GetMinDistanceIndex(const CacheStructs &css,
                                    const std::vector<cache::type_id_t> &var_type,
                                    const std::vector<DataArray> &read_sets,
                                    cache::CacheAccess access) const {
    int min_index = -1;
    if (css.size() == 0 || !css.at(0)->IsAvailable(access))
        return min_index;
    cache::distance_t min_distance = css.at(0)->GetDistance(var_type, read_sets);
    cache::distance_t dist;
    for (size_t i = 1; i < css.size(); ++i) {
        if (!css.at(i)->IsAvailable(access))
            continue;
        dist = css.at(i)->GetDistance(var_type, read_sets);
        if (dist == 0) {
            min_distance = 0;
            min_index = i;
            break;
        } else if (dist < min_distance) {
            min_distance = dist;
            min_index = i;
        }
    }
    return min_index;
}

}  // namespace nimbus
