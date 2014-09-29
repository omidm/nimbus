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
 * A CacheObject is an application object corresponding to one/ multiple nimbus
 * variables.
 *
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#include <map>
#include <set>
#include <string>
#include <vector>

#include "data/cache/cache_defs.h"
#include "data/cache/cache_object.h"
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "shared/nimbus_types.h"
#include "worker/data.h"

namespace nimbus {

// TODO(concurrency) needs to be protected.
cache::co_id_t CacheObject::ids_allocated_ = 0;

/**
 * \details
 */
CacheObject::CacheObject() : unique_id_(0), pending_flag_(false), id_(0),
    access_(cache::SHARED), users_(0) {
  name_ = "default";
}

/**
 * \details
 */
CacheObject::CacheObject(const GeometricRegion &ob_reg) : unique_id_(0),
    pending_flag_(false), id_(0),
    access_(cache::SHARED), users_(0),
    object_region_(ob_reg) {
  name_ = "default";
}

/**
 * \details MakePrototype() increases ids_allocated_ for CacheObject
 * prototypes, and allocates a new id to the prototype. A prototype is used by
 * application when requesting an application object. CacheManager uses
 * prototype id to put all instances of a prototype together - if a cached
 * instance satisfies requested region and prototype id, the CacheManager can
 * return the instance (after updating to reflect the read set), provided the
 * instance is available.
 */
void CacheObject::MakePrototype() {
// TODO(concurrency) needs to be protected.
    id_ = ++ids_allocated_;
}

/**
 * \details AcquireAccess(...) ensures that only one request in cache::EXCLUSIVE mode
 * holds the object, otherwise the object is in cache::SHARED mode. It sets object
 * access mode to requested access mode, and increases the number of users by
 * one.
 */
void CacheObject::AcquireAccess(cache::CacheAccess access) {
    if (!(users_ == 0 || (access == cache::SHARED && access_ == cache::SHARED))) {
        dbg(DBG_ERROR, "Cannot acquire cache object, it is already in use!");
        exit(-1);
    }
    access_ = access;
    users_++;
}

/**
 * \details ReleaseAccess() decreases the number of users by one. A request
 * (example, application job) must release an object when it is done reading/
 * writing.
 */
void CacheObject::ReleaseAccessInternal() {
    users_--;
    if (users_ != 0 && access_ == cache::EXCLUSIVE) {
        dbg(DBG_ERROR, "Incocistency in number of users using cache object!");
        exit(-1);
    }
}

/**
 * \details IsAvailable(...) returns true if an object is available, otherwise
 * it returns false. An object is available in cache::EXCLUSIVE mode if the number of
 * users using it is zero. An object is available in cache::SHARED mode if number of
 * users is zero, or the current access mode for the object is cache::SHARED.
 */
bool CacheObject::IsAvailable(cache::CacheAccess access) const {
    return ((access == cache::EXCLUSIVE && users_ == 0) ||
            (access == cache::SHARED && (access_ == cache::SHARED || users_ == 0)));
}

/**
 * \details
 */
cache::co_id_t CacheObject::id() const {
    return id_;
}

/**
 * \details
 */
GeometricRegion CacheObject::object_region() const {
    return object_region_;
}

/**
 * \details
 */
void CacheObject::set_object_region(const GeometricRegion &object_region) {
    object_region_ = object_region;
}

}  // namespace nimbus
