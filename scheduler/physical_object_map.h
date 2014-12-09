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
  * Data structure for maintaining and querying for the physical instances
  * of each logical data object.
  * 
  */

#ifndef NIMBUS_SCHEDULER_PHYSICAL_OBJECT_MAP_H_
#define NIMBUS_SCHEDULER_PHYSICAL_OBJECT_MAP_H_

#include <boost/unordered_map.hpp>
#include <map>
#include <string>
#include <vector>
#include <utility>
#include "scheduler/physical_data.h"
#include "shared/logical_data_object.h"
#include "shared/geometric_region.h"
#include "shared/nimbus_types.h"


namespace nimbus {

  typedef boost::unordered_map<logical_data_id_t, PhysicalDataList*> PhysicalObjectMapType;

  class PhysicalObjectMap {
  public:
    PhysicalObjectMap();
    virtual ~PhysicalObjectMap();

    virtual bool AddLogicalObject(LogicalDataObject* object);

    virtual bool RemoveLogicalObject(logical_data_id_t id);

    virtual bool AddPhysicalInstance(LogicalDataObject* object,
                                     const PhysicalData& instance);

    virtual bool RemovePhysicalInstance(LogicalDataObject* object,
                                        const PhysicalData& instance);

    virtual size_t RemoveAllInstanceByWorker(worker_id_t worker_id);

    virtual size_t ResetAllInstances();

    virtual bool UpdatePhysicalInstance(LogicalDataObject* object,
                                        const PhysicalData& old_instance,
                                        const PhysicalData& new_instance);

    virtual const PhysicalDataList* AllInstances(LogicalDataObject* object);

    virtual int AllInstances(LogicalDataObject* object,
                             PhysicalDataList* dest);

    virtual int InstancesByWorker(LogicalDataObject* object,
                                  worker_id_t worker,
                                  PhysicalDataList* dest);

    virtual int InstancesByVersion(LogicalDataObject* object,
                                   data_version_t version,
                                   PhysicalDataList* dest);

    virtual int InstancesByWorkerAndVersion(LogicalDataObject* object,
                                  worker_id_t worker,
                                  data_version_t version,
                                  PhysicalDataList* dest);

  private:
    PhysicalObjectMapType data_map_;
  };
}  // namespace nimbus

#endif  // NIMBUS_SCHEDULER_PHYSICAL_OBJECT_MAP_H_
