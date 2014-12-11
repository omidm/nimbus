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


#include "scheduler/region_map.h"

#define WEIGHT_NUM 8
#define WEIGHT_X {1, 1, 1, 1, 1, 1, 1, 1}
#define WEIGHT_Y {1, 1, 1, 1, 1, 1, 1, 1}
#define WEIGHT_Z {1, 1, 1, 1, 1, 1, 1, 1}

namespace nimbus {

RegionMap::RegionMap() {
}

RegionMap::RegionMap(const Table& table)
  : table_(table) {
}

RegionMap::RegionMap(const RegionMap& other) {
  table_ = other.table_;
  cache_ = other.cache_;
}

RegionMap::~RegionMap() {
  ClearTable();
}

RegionMap::Table RegionMap::table() const {
  return table_;
}

const RegionMap::Table* RegionMap::table_p() const {
  return &table_;
}

RegionMap::Table* RegionMap::table_p() {
  return &table_;
}

size_t RegionMap::table_size() {
  return table_.size();
}

void RegionMap::ClearTable() {
  TableIter iter = table_.begin();
  for (; iter != table_.end(); ++iter) {
    delete iter->second;
  }

  table_.clear();

  InvalidateCache();
}

void RegionMap::set_table(const Table& table) {
  table_ = table;
}

void RegionMap::TrackRegionCoverage(DataManager *data_manager,
                                    JobEntry *job,
                                    const worker_id_t *worker_id) {
  TableIter iter = table_.find(*worker_id);
  if (iter == table_.end()) {
    dbg(DBG_WARN, "WARNING: RegionMap: worker id %lu not in the table to track region coverage.\n", *worker_id); // NOLINT
    return;
  }

  GeometricRegion region;
  if (job->GetRegion(&region)) {
  } else if (job->GetWriteSetRegion(data_manager, &region)) {
    if (region.GetSurfaceArea() < (0.25 * global_region_.GetSurfaceArea())) {
      iter->second->AddCoveredRegion(&region);
    }
  }
}


bool RegionMap::FindWorkerWithMostOverlappedRegion(const GeometricRegion *region,
                                                   worker_id_t *worker_id) {
  assert(table_.size() > 0);
  TableIter iter = table_.begin();
  worker_id_t w_id = iter->first;
  int_dimension_t common_area = iter->second->CommonSurface(region);
  ++iter;
  for (; iter != table_.end(); ++iter) {
    int_dimension_t temp = iter->second->CommonSurface(region);
    if (temp > common_area) {
      common_area = temp;
      w_id = iter->first;
    }
  }

  *worker_id = w_id;
  return true;
}


bool RegionMap::QueryWorkerWithMostOverlap(DataManager *data_manager,
                                           JobEntry *job,
                                           worker_id_t *worker_id) {
  if (table_size() == 0) {
    return false;
  }

  GeometricRegion region;
  if (job->GetRegion(&region)) {
  } else if (job->GetUnionSetRegion(data_manager, &region)) {
  } else {
    dbg(DBG_ERROR, "ERROR: RegionMap: could not get any region for job %lu.\n", job->job_id()); // NOLINT
    return false;
  }

  if (region.GetSurfaceArea() < (0.25 * global_region_.GetSurfaceArea())) {
    /* 
     * If the union set is not the global bounding then you can find the
     * intersect with only union set region not each ldo one by one. If there
     * are global variables being read by jobs, then the region of a job
     * becomes the entire space.
     */

    if (QueryCache(&region, worker_id)) {
      return true;
    }

    FindWorkerWithMostOverlappedRegion(&region, worker_id);
    CacheQueryResult(&region, worker_id);
    return true;
  } else {
    /*
     * This method is weighted polling among the workers, the worker that has
     * the highest aggregate intersect volume with each of the ldos in the
     * union set of the job would win the job.
     */

    // initialize worker ranks.
    WorkerRank worker_ranks;
    TableIter tabit = table_.begin();
    for (; tabit != table_.end(); ++tabit) {
      worker_ranks[tabit->first] = 0;
    }

    IDSet<logical_data_id_t>::IDSetIter iter;
    for (iter = job->union_set_p()->begin(); iter != job->union_set_p()->end(); ++iter) {
      const LogicalDataObject* ldo;
      ldo = data_manager->FindLogicalObject(*iter);
      WorkerRank::iterator it = worker_ranks.begin();
      for (; it != worker_ranks.end(); ++it) {
        it->second += table_[it->first]->CommonSurface(ldo->region());
      }
    }

    // find the worker that wins the poll.
    WorkerRank::iterator writ = worker_ranks.begin();
    assert(writ != worker_ranks.end());
    worker_id_t w_id = writ->first;
    int_dimension_t vol = worker_ranks[writ->first];
    for (; writ != worker_ranks.end(); ++writ) {
      if (vol < worker_ranks[writ->first]) {
        vol = worker_ranks[writ->first];
        w_id = writ->first;
      }
    }

    *worker_id = w_id;
    return true;
  }
}

bool RegionMap::WorkersAreNeighbor(worker_id_t first, worker_id_t second) {
  TableIter iter_first = table_.find(first);
  if (iter_first == table_.end()) {
    dbg(DBG_ERROR, "ERROR: RegionMap: worker %lu is not defined in region map.\n", first);
    return false;
  }

  TableIter iter_second = table_.find(second);
  if (iter_second == table_.end()) {
    dbg(DBG_ERROR, "ERROR: RegionMap: worker %lu is not defined in region map.\n", second);
    return false;
  }

  return iter_first->second->AdjacentOrIntersects(iter_second->second);
}


void RegionMap::Initialize(const std::vector<worker_id_t>& worker_ids,
                           const GeometricRegion& global_region) {
  global_region_ = global_region;

  size_t num_x, num_y, num_z;
  SplitDimensions(worker_ids.size(), &num_x, &num_y, &num_z);

  static const int arr_x[] = WEIGHT_X;
  std::vector<size_t> weight_x(arr_x, arr_x + WEIGHT_NUM);
  static const int arr_y[] = WEIGHT_Y;
  std::vector<size_t> weight_y(arr_y, arr_y + WEIGHT_NUM);
  static const int arr_z[] = WEIGHT_Z;
  std::vector<size_t> weight_z(arr_z, arr_z + WEIGHT_NUM);

  GenerateTable(num_x, num_y, num_z,
                weight_x, weight_y, weight_z,
                worker_ids, global_region);

  InvalidateCache();
}


bool RegionMap::NotifyDownWorker(worker_id_t worker_id) {
  TableIter iter_down = table_.find(worker_id);
  if (iter_down == table_.end()) {
    dbg(DBG_ERROR, "ERROR: RegionMap: worker %lu is not defined in region map.\n", worker_id);
    exit(-1);
    return false;
  }

  bool found_region_to_capture = false;
  TableIter iter = table_.begin();
  for (; iter != table_.end(); ++iter) {
    if (iter == iter_down) {
      continue;
    }
    if (iter->second->AdjacentOrIntersects(iter_down->second)) {
      found_region_to_capture = true;
      break;
    }
  }

  if (!found_region_to_capture) {
    dbg(DBG_ERROR, "ERROR: RegionMap: could not find worker to capture worker %lu region.\n", worker_id); // NOLINT
    exit(-1);
    return false;
  }

  dbg(DBG_WARN, "WARNING: RegionMap: worker %lu capturing worker %lu region.\n", iter->first, worker_id); // NOLINT

  iter->second->Grow(iter_down->second);
  delete iter_down->second;
  table_.erase(iter_down);

  InvalidateCache();
  InvalidateRegionCoverage();

  return true;
}

bool RegionMap::QueryCache(const GeometricRegion *region,
                           worker_id_t *worker_id) {
  CacheIter iter = cache_.find(*region);
  if (iter != cache_.end()) {
    *worker_id = iter->second;
    return true;
  }

  return false;
}

void RegionMap::InvalidateCache() {
  cache_.clear();
}

void RegionMap::CacheQueryResult(const GeometricRegion *region,
                                 const worker_id_t *worker_id) {
  cache_[*region] = *worker_id;
}

void RegionMap::InvalidateRegionCoverage() {
  TableIter iter = table_.begin();
  for (; iter != table_.end(); ++iter) {
    iter->second->ClearCoveredRegions();
  }
}



bool RegionMap::BalanceRegions(const worker_id_t &w_grow,
                               const worker_id_t &w_shrink) {
  if (w_shrink == w_grow) {
    dbg(DBG_ERROR, "ERROR: RegionMap: workers are the same for balancing the regions.\n"); // NOLINT
    return false;
  }

  if (!WorkersAreNeighbor(w_shrink, w_grow)) {
    dbg(DBG_ERROR, "ERROR: RegionMap: workers are not neighbors for balancing the regions.\n"); // NOLINT
    return false;
  }

  TableIter iter_shrink = table_.find(w_shrink);
  if (iter_shrink == table_.end()) {
    dbg(DBG_ERROR, "ERROR: RegionMap: worker %lu is not defined in region map.\n", w_shrink);
    return false;
  }

  TableIter iter_grow = table_.find(w_grow);
  if (iter_grow == table_.end()) {
    dbg(DBG_ERROR, "ERROR: RegionMap: worker %lu is not defined in region map.\n", w_grow);
    return false;
  }

  GeometricRegion region;
  if (!iter_shrink->second->GetRegionToGiveUp(iter_grow->second, &region)) {
    return false;
  }
  iter_shrink->second->Shrink(&region);
  iter_grow->second->Grow(&region);

  InvalidateCache();
  InvalidateRegionCoverage();

  return true;
}


std::string RegionMap::Print() {
  std::string rval;
  rval += "\n+++++++ Region Map Begin +++++++\n";

  TableIter iter = table_.begin();
  for (; iter != table_.end(); ++iter) {
    std::ostringstream ss;
    rval += "worker_id: ";
    ss << iter->first;
    rval += ss.str();
    ss.str(std::string());
    rval += "\n";
    rval += iter->second->PrintRegion();
  }

  rval += "\n++++++++ Region Map End ++++++++\n";
  return rval;
}

void RegionMap::GenerateTable(size_t num_x, size_t num_y, size_t num_z,
                              std::vector<size_t> weight_x,
                              std::vector<size_t> weight_y,
                              std::vector<size_t> weight_z,
                              const std::vector<worker_id_t>& worker_ids,
                              const GeometricRegion& global_region) {
  assert(weight_x.size() >= num_x);
  assert(weight_y.size() >= num_y);
  assert(weight_z.size() >= num_z);

  std::vector<int_dimension_t> width_x;
  size_t weight_sum_x = 0;
  for (size_t i = 0; i < num_x; ++i) {
    weight_sum_x += weight_x[i];
  }
  for (size_t i = 0; i < num_x; ++i) {
    width_x.push_back(global_region.dx() * weight_x[i] / weight_sum_x);
  }
  std::vector<int_dimension_t> marker_x;
  marker_x.push_back(global_region.x());
  for (size_t i = 0; i < num_x; ++i) {
    marker_x.push_back(marker_x[i] + width_x[i]);
  }


  std::vector<int_dimension_t> width_y;
  size_t weight_sum_y = 0;
  for (size_t i = 0; i < num_y; ++i) {
    weight_sum_y += weight_y[i];
  }
  for (size_t i = 0; i < num_y; ++i) {
    width_y.push_back(global_region.dy() * weight_y[i] / weight_sum_y);
  }
  std::vector<int_dimension_t> marker_y;
  marker_y.push_back(global_region.y());
  for (size_t i = 0; i < num_y; ++i) {
    marker_y.push_back(marker_y[i] + width_y[i]);
  }

  std::vector<int_dimension_t> width_z;
  size_t weight_sum_z = 0;
  for (size_t i = 0; i < num_z; ++i) {
    weight_sum_z += weight_z[i];
  }
  for (size_t i = 0; i < num_z; ++i) {
    width_z.push_back(global_region.dz() * weight_z[i] / weight_sum_z);
  }
  std::vector<int_dimension_t> marker_z;
  marker_z.push_back(global_region.z());
  for (size_t i = 0; i < num_z; ++i) {
    marker_z.push_back(marker_z[i] + width_z[i]);
  }

  std::vector<RegionMapEntry*> domains;
  for (size_t i = 0; i < num_x; ++i) {
    for (size_t j = 0; j < num_y; ++j) {
      for (size_t k = 0; k < num_z; ++k) {
        RegionMapEntry *rme = new RegionMapEntry();
        GeometricRegion r(marker_x[i],
                          marker_y[j],
                          marker_z[k],
                          marker_x[i + 1] - marker_x[i],
                          marker_y[j + 1] - marker_y[j],
                          marker_z[k + 1] - marker_z[k]);
        rme->Grow(&r);
        domains.push_back(rme);
      }
    }
  }

  ClearTable();

  size_t index = 0;
  assert(domains.size() == worker_ids.size());
  std::vector<RegionMapEntry*>::iterator iter = domains.begin();
  for (; iter != domains.end(); ++iter) {
    table_[worker_ids[index]] = *iter;
    ++index;
  }
}

void RegionMap::SplitDimensions(size_t worker_num, size_t *num_x, size_t *num_y, size_t *num_z) {
  switch (worker_num) {
    case 1 :
      *num_x = 1;
      *num_y = 1;
      *num_z = 1;
      break;
    case 2 :
      *num_x = 2;
      *num_y = 1;
      *num_z = 1;
      break;
    case 3 :
      *num_x = 3;
      *num_y = 1;
      *num_z = 1;
      break;
    case 4 :
      *num_x = 4;
      *num_y = 1;
      *num_z = 1;
      break;
    case 5 :
      *num_x = 5;
      *num_y = 1;
      *num_z = 1;
      break;
    case 6 :
      *num_x = 6;
      *num_y = 1;
      *num_z = 1;
      break;
    case 7 :
      *num_x = 7;
      *num_y = 1;
      *num_z = 1;
      break;
    case 8 :
      *num_x = 8;
      *num_y = 1;
      *num_z = 1;
      break;
    case 100 :
      *num_x = 5;
      *num_y = 5;
      *num_z = 4;
      break;
    default:
      dbg(DBG_ERROR, "ERROR: Do not know how to split!");
      exit(-1);
  }
}



RegionMap& RegionMap::operator= (
    const RegionMap& right) {
  table_ = right.table_;
  cache_ = right.cache_;
  return *this;
}

}  // namespace nimbus
