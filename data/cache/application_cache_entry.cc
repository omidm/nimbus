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
 * Cache table entry containing object type, a pointer to the cached object,
 * and a list of tables for the fields that the cached object contains.
 *
 * Author: Chinmayee Shah <chshah@stanford.edu>
 */

#include "data/cache/application_cache_entry.h"
#include "worker/data.h"
#include "worker/job.h"

namespace nimbus {

ApplicationCacheEntry::ApplicationCacheEntry()
    : object_(NULL), fields_(NULL) {}

std::string ApplicationCacheEntry::object_type() {
    return object_type_;
}

void ApplicationCacheEntry::set_object_type(std::string type) {
    object_type_ = type;
}

void* ApplicationCacheEntry::object() {
    return object_;
}

void ApplicationCacheEntry::set_object(void *object) {
    object_ = object;
}

void ApplicationCacheEntry::LockData(const Job &job,
                                     const DataArray &da) {
    IDSet<nimbus::physical_data_id_t> read_set = job.read_set();
    IDSet<nimbus::physical_data_id_t> write_set = job.write_set();
    for (size_t i = 0; i < da.size(); i++) {
        Data *d = da[i];
        if (write_set.contains(d->physical_id())) {
            ApplicationField field = fields_->at(d->name());
            field.LockWrite(d->region());
        }
    }
}

}  // namespace nimbus
