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
  * Version Map.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SCHEDULER_VERSION_MAP_H_
#define NIMBUS_SCHEDULER_VERSION_MAP_H_

#include <boost/shared_ptr.hpp>
#include <map>
#include "shared/nimbus_types.h"
#include "shared/dbg.h"

namespace nimbus {

class VersionMap {
  public:
    typedef std::map<logical_data_id_t, data_version_t> Map;
    typedef std::map<logical_data_id_t, data_version_t>::iterator Iter;
    typedef std::map<logical_data_id_t, data_version_t>::const_iterator ConstIter;

    VersionMap();
    virtual ~VersionMap();

    Map content() const;
    const Map* content_p() const;
    bool query_entry(logical_data_id_t l_id, data_version_t *version) const;

    void set_content(const Map& content);
    void set_entry(logical_data_id_t l_id, data_version_t version);

    void Print() const;

  private:
    Map content_;
};


}  // namespace nimbus
#endif  // NIMBUS_SCHEDULER_VERSION_MAP_H_
