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

#ifndef NIMBUS_DATA_CACHE_APPLICATION_CACHE_ENTRY_H_
#define NIMBUS_DATA_CACHE_APPLICATION_CACHE_ENTRY_H_

#include <string>
#include <vector>

#include "data/cache/application_field.h"
#include "worker/data.h"
#include "worker/job.h"

namespace nimbus {

class ApplicationCacheEntry {
    public:
        ApplicationCacheEntry();

        /* accessors */
        std::string object_type();
        void set_object_type(std::string type);
        void* object();
        void set_object(void *object);

        /* get read/ write locks on cache data */
        void LockData(const Job &job,
                      const DataArray &da);

    private:
        std::string object_type_;
        void *object_;
        CacheObjectFieldMap *fields_;
};  // class ApplicationCacheEntry

typedef std::vector<ApplicationCacheEntry> ApplicationCacheEntries;

}  // namespace nimbus

#endif  // NIMBUS_DATA_CACHE_APPLICATION_CACHE_ENTRY_H_
