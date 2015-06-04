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
  * The Data Manager maintains the controller's metadata on what
  * logical data objects are in the system as well as their physical
  * instances.
  */

#ifndef NIMBUS_SCHEDULER_DATA_MANAGER_H_
#define NIMBUS_SCHEDULER_DATA_MANAGER_H_

#include <map>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include "shared/nimbus_types.h"
#include "shared/logical_data_object.h"
#include "shared/ldo_index_rtree.h"
#include "shared/scheduler_server.h"
#include "scheduler/physical_data.h"
#include "scheduler/physical_object_map.h"

namespace nimbus {

  class DataManager {
  public:
    DataManager();
    virtual ~DataManager();

    /* Managing geometric partitions. */
    bool AddPartition(partition_id_t id, GeometricRegion region);
    bool RemovePartition(partition_id_t id);
    bool HasPartition(partition_id_t id);
    bool FindPartition(partition_id_t id, GeometricRegion* r);

    /* Managing logical objects. */
    bool AddLogicalObject(logical_data_id_t id,
                          std::string variable,
                          GeometricRegion region,
                          partition_id_t partition = 0);
    bool AddLogicalObject(logical_data_id_t id,
                          std::string variable,
                          partition_id_t partition);
    bool RemoveLogicalObject(logical_data_id_t id);
    const LogicalDataObject* FindLogicalObject(logical_data_id_t id);
    int FindLogicalObjects(std::string variable,
                           CLdoVector* dest);
    int FindIntersectingLogicalObjects(std::string variable,
                                       GeometricRegion* region,
                                       CLdoVector* dest);
    int FindCoveredLogicalObjects(std::string variable,
                                  GeometricRegion* region,
                                  CLdoVector* dest);
    int FindAdjacentLogicalObjects(std::string variable,
                                   GeometricRegion* region,
                                   CLdoVector* dest);

    /* Managing physical instances of logical objects. */
    bool AddPhysicalInstance(LogicalDataObject* object,
                             const PhysicalData& instance);
    bool RemovePhysicalInstance(LogicalDataObject* object,
                                const PhysicalData& instance);

    virtual size_t RemoveAllInstanceByWorker(worker_id_t worker_id);

    virtual size_t ResetAllInstances();

    bool UpdatePhysicalInstance(LogicalDataObject* object,
                                const PhysicalData& old_instance,
                                const PhysicalData& new_instance);

    bool UpdateVersionAndAccessRecord(const logical_data_id_t& ldid,
                                      const physical_data_id_t& pdid,
                                      const data_version_t& version,
                                      const IDSet<job_id_t>& list_job_read,
                                      const job_id_t& last_job_write);

    const PhysicalDataList* AllInstances(LogicalDataObject* object);
    int AllInstances(LogicalDataObject* object,
                     PhysicalDataList* dest);
    int AllInstances(LogicalDataObject* object,
                     ConstPhysicalDataPList* dest);
    int InstancesByWorker(LogicalDataObject* object,
                          worker_id_t worker,
                          PhysicalDataList* dest);
    int InstancesByVersion(LogicalDataObject* object,
                           data_version_t version,
                           PhysicalDataList* dest);
    int InstancesByWorkerAndVersion(LogicalDataObject* object,
                          worker_id_t worker,
                          data_version_t version,
                          PhysicalDataList* dest);

    partition_id_t max_defined_partition();

    GeometricRegion global_bounding_region();

    bool initialized_global_bounding_region();

    const LdoMap* ldo_map_p();

  private:
    PhysicalObjectMap physical_object_map_;
    LdoIndexRtree ldo_index_;
    LdoMap ldo_map_;
    PartitionMap partition_map_;
    partition_id_t max_defined_partition_;
    GeometricRegion global_bounding_region_;
    bool initialized_global_bounding_region_;
  };
}  // namespace nimbus

#endif  // NIMBUS_SCHEDULER_DATA_MANAGER_H_
