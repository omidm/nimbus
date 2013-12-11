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
  * A physical data object is an piece of data on a worker node that an
  * application can access. Every physical object is an instance of a
  * logical object, with a corresponding version number. A worker can
  * have multiple physical objects for one logical object, with the same
  * or different version numbers.
  */

#ifndef NIMBUS_WORKER_PHYSICAL_DATA_OBJECT_H_
#define NIMBUS_WORKER_PHYSICAL_DATA_OBJECT_H_

#include <list>
#include <set>
#include <string>
#include <vector>
#include "shared/nimbus_types.h"
#include "shared/geometric_region.h"
#include "shared/logical_data_object.h"
#include "worker/data.h"

namespace nimbus {

  class PhysicalDataInstance {
  public:
    PhysicalDataInstance();
    PhysicalDataInstance(physical_data_id_t id,
                         LogicalDataObject* lobj,
                         Data* data,
                         data_version_t version);

    virtual ~PhysicalDataInstance();

    virtual physical_data_id_t id() const;
    virtual std::string variable() const;
    virtual GeometricRegion* region() const;
    virtual partition_id_t partition() const;
    virtual LogicalDataObject* logical_object() const;
    virtual data_version_t version() const;
    virtual void set_version(data_version_t version);
    virtual data_version_t IncrementVersion();
    virtual Data* data() const;

  private:
    physical_data_id_t id_;
    data_version_t version_;
    LogicalDataObject* logical_object_;
    Data* data_;
  };

  typedef std::set<PhysicalDataInstance*> PdiSet;
  typedef std::list<PhysicalDataInstance*> PdiList;
  typedef std::vector<PhysicalDataInstance*> PdiVector;
  typedef std::vector<const PhysicalDataInstance*> CPdiVector;

}  // namespace nimbus

#endif  // NIMBUS_WORKER_PHYSICAL_DATA_OBJECT_H_
