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
  * Class used to define data objects and partitions.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include <assert.h>
#include "application_utils/data_definer.h"

using namespace nimbus; // NOLINT

DataDefiner::DataDefiner(Job *job) {
  job_ = job;
  partition_id_ = 0;
}

DataDefiner::~DataDefiner() {
}

bool DataDefiner::DefineData(const std::string& name,
                             int_dimension_t size_x,
                             int_dimension_t size_y,
                             int_dimension_t size_z,
                             int_dimension_t bw_x,
                             int_dimension_t bw_y,
                             int_dimension_t bw_z,
                             int_dimension_t part_num_x,
                             int_dimension_t part_num_y,
                             int_dimension_t part_num_z,
                             bool with_global_bw) {
  IDSet<nimbus::partition_id_t> neighbor_partitions;

  SpaceIterator iter_x(size_x, bw_x, part_num_x, with_global_bw);
  SpaceIterator iter_y(size_y, bw_y, part_num_y, with_global_bw);
  SpaceIterator iter_z(size_z, bw_z, part_num_z, with_global_bw);

  SpaceIterator::Cursor cursor_x;
  SpaceIterator::Cursor cursor_y;
  SpaceIterator::Cursor cursor_z;

  iter_x.Initialize();
  do {
    iter_x.ReadCursor(&cursor_x);
    iter_y.Initialize();
    do {
      iter_y.ReadCursor(&cursor_y);
      iter_z.Initialize();
      do {
        iter_z.ReadCursor(&cursor_z);

        GeometricRegion region(cursor_x.point_,
                               cursor_y.point_,
                               cursor_z.point_,
                               cursor_x.delta_,
                               cursor_y.delta_,
                               cursor_z.delta_);

        partition_id_t p_id = 0;
        PartitionMap::iterator iter = part_map_.find(region.ToNetworkData());
        if (iter != part_map_.end()) {
          p_id = iter->second;
        } else {
          job_->DefinePartition(nimbus::ID<partition_id_t>(partition_id_), region);
          part_map_[region.ToNetworkData()] = partition_id_;
          p_id = partition_id_;
          ++partition_id_;
        }

        std::vector<nimbus::logical_data_id_t> data_id;
        job_->GetNewLogicalDataID(&data_id, 1);

        job_->DefineData(name, data_id[0], p_id, neighbor_partitions);
      } while (iter_z.Advance());
    } while (iter_y.Advance());
  } while (iter_x.Advance());

  return true;
}

bool DataDefiner::DefineScratchData(const std::string& name,
                                    int_dimension_t size_x,
                                    int_dimension_t size_y,
                                    int_dimension_t size_z,
                                    int_dimension_t bw_x,
                                    int_dimension_t bw_y,
                                    int_dimension_t bw_z,
                                    int_dimension_t part_num_x,
                                    int_dimension_t part_num_y,
                                    int_dimension_t part_num_z,
                                    bool with_global_bw) {
  IDSet<nimbus::partition_id_t> neighbor_partitions;

  SpaceIterator iter_x(size_x, bw_x, part_num_x, with_global_bw);
  SpaceIterator iter_y(size_y, bw_y, part_num_y, with_global_bw);
  SpaceIterator iter_z(size_z, bw_z, part_num_z, with_global_bw);

  SpaceIterator::Cursor cursor_x;
  SpaceIterator::Cursor cursor_y;
  SpaceIterator::Cursor cursor_z;

  iter_x.Initialize();
  do {
    iter_x.ReadCursor(&cursor_x);
    iter_y.Initialize();
    do {
      iter_y.ReadCursor(&cursor_y);
      iter_z.Initialize();
      do {
        iter_z.ReadCursor(&cursor_z);

        GeometricRegion region(cursor_x.point_,
                               cursor_y.point_,
                               cursor_z.point_,
                               cursor_x.delta_,
                               cursor_y.delta_,
                               cursor_z.delta_);

        partition_id_t p_id = 0;
        PartitionMap::iterator iter = part_map_.find(region.ToNetworkData());
        if (iter != part_map_.end()) {
          p_id = iter->second;
        } else {
          job_->DefinePartition(nimbus::ID<partition_id_t>(partition_id_), region);
          part_map_[region.ToNetworkData()] = partition_id_;
          p_id = partition_id_;
          ++partition_id_;
        }

        if ((cursor_x.type_ + cursor_y.type_ + cursor_z.type_) == 0) {
          std::vector<nimbus::logical_data_id_t> data_ids;
          job_->GetNewLogicalDataID(&data_ids, 8);
          job_->DefineData(name + "_vertex_1", data_ids[0], p_id, neighbor_partitions);
          job_->DefineData(name + "_vertex_2", data_ids[1], p_id, neighbor_partitions);
          job_->DefineData(name + "_vertex_3", data_ids[2], p_id, neighbor_partitions);
          job_->DefineData(name + "_vertex_4", data_ids[3], p_id, neighbor_partitions);
          job_->DefineData(name + "_vertex_5", data_ids[4], p_id, neighbor_partitions);
          job_->DefineData(name + "_vertex_6", data_ids[5], p_id, neighbor_partitions);
          job_->DefineData(name + "_vertex_7", data_ids[6], p_id, neighbor_partitions);
          job_->DefineData(name + "_vertex_8", data_ids[7], p_id, neighbor_partitions);
        } else if ((cursor_x.type_ + cursor_y.type_ + cursor_z.type_) == 1) {
          std::vector<nimbus::logical_data_id_t> data_ids;
          job_->GetNewLogicalDataID(&data_ids, 4);
          job_->DefineData(name + "_edge_1", data_ids[0], p_id, neighbor_partitions);
          job_->DefineData(name + "_edge_2", data_ids[1], p_id, neighbor_partitions);
          job_->DefineData(name + "_edge_3", data_ids[2], p_id, neighbor_partitions);
          job_->DefineData(name + "_edge_4", data_ids[3], p_id, neighbor_partitions);
        } else if ((cursor_x.type_ + cursor_y.type_ + cursor_z.type_) == 2) {
          std::vector<nimbus::logical_data_id_t> data_ids;
          job_->GetNewLogicalDataID(&data_ids, 2);
          job_->DefineData(name + "_face_1", data_ids[0], p_id, neighbor_partitions);
          job_->DefineData(name + "_face_2", data_ids[1], p_id, neighbor_partitions);
        }
      } while (iter_z.Advance());
    } while (iter_y.Advance());
  } while (iter_x.Advance());

  return true;
}



