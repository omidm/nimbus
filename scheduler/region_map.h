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
  * This class holds the mapping between workers and the region they cover for
  * computation. This region is not necessarily a box and could be an
  * unstructured region. It also provides lookup facilities to assign jobs to
  * workers that match the job region.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SCHEDULER_REGION_MAP_H_
#define NIMBUS_SCHEDULER_REGION_MAP_H_

#include <boost/unordered_map.hpp>
#include <list>
#include <map>
#include <vector>
#include <utility>
#include "shared/nimbus_types.h"
#include "shared/dbg.h"
#include "shared/geometric_region.h"
#include "scheduler/region_map_entry.h"

namespace nimbus {

  class RegionMap {
  public:
    typedef boost::unordered_map<worker_id_t, RegionMapEntry*> Table;
    typedef Table::iterator TableIter;
    typedef std::map<GeometricRegion, worker_id_t> Cache;
    typedef Cache::iterator CacheIter;

    RegionMap();
    explicit RegionMap(const Table& table);

    RegionMap(const RegionMap& other);

    virtual ~RegionMap();

    Table table() const;
    const Table* table_p() const;
    Table* table_p();
    size_t table_size();
    void set_table(const Table& table);

    void Initialize(const std::vector<worker_id_t>& worker_ids,
                    const GeometricRegion& global_region);

    bool QueryWorkerWithMostOverlap(const GeometricRegion *region,
                                      worker_id_t *worker_id);

    bool WorkersAreNeighbor(worker_id_t first, worker_id_t second);

    RegionMap& operator= (const RegionMap& right);

  private:
    Table table_;
    Cache cache_;

    void ClearTable();

    void SplitDimensions(size_t worker_num, size_t *num_x, size_t *num_y, size_t *num_z);

    void GenerateTable(size_t num_x, size_t num_y, size_t num_z,
                       std::vector<size_t> weight_x,
                       std::vector<size_t> weight_y,
                       std::vector<size_t> weight_z,
                       const std::vector<worker_id_t>& worker_ids,
                       const GeometricRegion& global_region);

    bool QueryCache(const GeometricRegion *region,
                    worker_id_t *worker_id);

    void InvalidateCache();

    void CacheQueryResult(const GeometricRegion *region,
                          const worker_id_t *worker_id);

    void TrackRegionCoverage(const GeometricRegion *region,
                             const worker_id_t *worker_id);

    void InvalidateRegionCoverage();
  };

}  // namespace nimbus

#endif  // NIMBUS_SCHEDULER_REGION_MAP_H_
