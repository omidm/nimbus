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
  * Dense Version Map class. Instead of using a map, use an array and lookup
  * with the index. Since, the assumption it that the index is dence the memory
  * usage would be a problem. The user of the class has to make sure that the
  * index is dense enough.
  *
  * NOTE: the quey result can never be NIMBUS_UNDEFINED_DATA_VERSION, so it
  * does not work for difference version maps - only absolute version maps. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SRC_SCHEDULER_DENSE_VERSION_MAP_H_
#define NIMBUS_SRC_SCHEDULER_DENSE_VERSION_MAP_H_

#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>
#include "src/shared/nimbus_types.h"
#include "src/shared/dbg.h"
#include "src/scheduler/version_map.h"

namespace nimbus {

class DenseVersionMap : public VersionMap {
  public:
    typedef std::vector<data_version_t> List;

    DenseVersionMap(logical_data_id_t min_id,
                    logical_data_id_t max_id);

    DenseVersionMap(const DenseVersionMap& other);
    virtual ~DenseVersionMap();

    virtual bool query_entry(logical_data_id_t l_id, data_version_t *version) const;

    virtual void set_entry(logical_data_id_t l_id, data_version_t version);

    virtual void Print() const;

    virtual DenseVersionMap& operator=(const DenseVersionMap& right);

  private:
    logical_data_id_t min_id_;
    logical_data_id_t max_id_;
    List list_;
};


}  // namespace nimbus
#endif  // NIMBUS_SRC_SCHEDULER_DENSE_VERSION_MAP_H_
