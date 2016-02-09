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
  * The Data Manager is thread-safe.
  */

#ifndef NIMBUS_SRC_WORKER_WORKER_LDO_MAP_H_
#define NIMBUS_SRC_WORKER_WORKER_LDO_MAP_H_

#include <pthread.h>
#include <map>
#include <string>
#include <vector>
#include <utility>
#include "src/shared/nimbus_types.h"
#include "src/shared/geometric_region.h"
#include "src/shared/logical_data_object.h"
#include "src/shared/ldo_index_rtree.h"
#include "src/shared/dbg.h"

namespace nimbus {

class WorkerLdoMap {
  public:
    WorkerLdoMap();
    virtual ~WorkerLdoMap();

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
    int FindLogicalObjects(const std::string& variable,
                           CLdoVector* dest);
    int FindIntersectingLogicalObjects(const std::string& variable,
                                       const GeometricRegion* region,
                                       CLdoVector* dest);
    int FindCoveredLogicalObjects(const std::string& variable,
                                  const GeometricRegion* region,
                                  CLdoVector* dest);
    int FindAdjacentLogicalObjects(const std::string& variable,
                                   const GeometricRegion* region,
                                   CLdoVector* dest);

  private:
    // Protects ldo_index_.
    pthread_mutexattr_t lockattr_ldo_index_;
    pthread_mutex_t lock_ldo_index_;
    LdoIndexRtree ldo_index_;
    // Protects partition_map_.
    pthread_mutexattr_t lockattr_partition_map_;
    pthread_mutex_t lock_partition_map_;
    std::map<partition_id_t, GeometricRegion> partition_map_;

    class LockGuard {
     public:
      explicit LockGuard(pthread_mutex_t* lock) {
        pthread_mutex_lock(lock);
        lock_ = lock;
      }
      virtual ~LockGuard() {
        pthread_mutex_unlock(lock_);
      }
     private:
      pthread_mutex_t* lock_;
    };
  };
}  // namespace nimbus

#endif  // NIMBUS_SRC_WORKER_WORKER_LDO_MAP_H_
