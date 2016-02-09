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
  * Class used to build and access the partitions, in the application.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include <assert.h>
#include "src/application_utils/partition_handler.h"

using namespace nimbus; // NOLINT

PartitionHandler::PartitionHandler() {
}

PartitionHandler::~PartitionHandler() {
}

void PartitionHandler::AddPartitions(const std::string& name,
                                    int_dimension_t size_x,
                                    int_dimension_t size_y,
                                    int_dimension_t size_z,
                                    int_dimension_t bw_x,
                                    int_dimension_t bw_y,
                                    int_dimension_t bw_z,
                                    int_dimension_t part_num_x,
                                    int_dimension_t part_num_y,
                                    int_dimension_t part_num_z,
                                    PartitionType type) {
  assert((size_x / part_num_x) > (2 * bw_x));
  assert((size_y / part_num_y) > (2 * bw_y));
  assert((size_z / part_num_z) > (2 * bw_z));

  std::string key = name + GetTag(type);
  Map::iterator iter = map_.find(key);
  if (iter != map_.end()) {
    assert(false);
  }

  SpaceIterator iter_x(size_x, 0, part_num_x, false);
  SpaceIterator iter_y(size_y, 0, part_num_y, false);
  SpaceIterator iter_z(size_z, 0, part_num_z, false);

  SpaceIterator::Cursor cursor_x;
  SpaceIterator::Cursor cursor_y;
  SpaceIterator::Cursor cursor_z;

  List list;

  iter_x.Initialize();
  do {
    iter_x.ReadCursor(&cursor_x);
    iter_y.Initialize();
    do {
      iter_y.ReadCursor(&cursor_y);
      iter_z.Initialize();
      do {
        iter_z.ReadCursor(&cursor_z);

        int_dimension_t x  = 0;
        int_dimension_t y  = 0;
        int_dimension_t z  = 0;
        int_dimension_t dx = 0;
        int_dimension_t dy = 0;
        int_dimension_t dz = 0;

        switch (type) {
          case INNER:
            x = cursor_x.point_ + bw_x;
            y = cursor_y.point_ + bw_y;
            z = cursor_z.point_ + bw_z;
            dx = (size_x / part_num_x) - (2 * bw_x);
            dy = (size_y / part_num_y) - (2 * bw_y);
            dz = (size_z / part_num_z) - (2 * bw_z);
            break;
          case CENTRAL:
            x = cursor_x.point_;
            y = cursor_y.point_;
            z = cursor_z.point_;
            dx = (size_x / part_num_x);
            dy = (size_y / part_num_y);
            dz = (size_z / part_num_z);
            break;
          case CENTRAL_WGB:
            if (cursor_x.point_ == SpaceIterator::START_INDEX) {
              x = cursor_x.point_ - bw_x;
              dx = (size_x / part_num_x) + bw_x;
            } else if (cursor_x.point_ == (SpaceIterator::START_INDEX + size_x - (size_x / part_num_x))) { // NOLINT
              x = cursor_x.point_;
              dx = (size_x / part_num_x) + bw_x;
            } else {
              x = cursor_x.point_;
              dx = (size_x / part_num_x);
            }

            if (cursor_y.point_ == SpaceIterator::START_INDEX) {
              y = cursor_y.point_ - bw_y;
              dy = (size_y / part_num_y) + bw_y;
            } else if (cursor_y.point_ == (SpaceIterator::START_INDEX + size_y - (size_y / part_num_y))) { // NOLINT
              y = cursor_y.point_;
              dy = (size_y / part_num_y) + bw_y;
            } else {
              y = cursor_y.point_;
              dy = (size_y / part_num_y);
            }

            if (cursor_z.point_ == SpaceIterator::START_INDEX) {
              z = cursor_z.point_ - bw_z;
              dz = (size_z / part_num_z) + bw_z;
            } else if (cursor_z.point_ == (SpaceIterator::START_INDEX + size_z - (size_z / part_num_z))) { // NOLINT
              z = cursor_z.point_;
              dz = (size_z / part_num_z) + bw_z;
            } else {
              z = cursor_z.point_;
              dz = (size_z / part_num_z);
            }
            break;
          case OUTER:
            x = cursor_x.point_ - bw_x;
            y = cursor_y.point_ - bw_y;
            z = cursor_z.point_ - bw_z;
            dx = (size_x / part_num_x) + (2 * bw_x);
            dy = (size_y / part_num_y) + (2 * bw_y);
            dz = (size_z / part_num_z) + (2 * bw_z);
            break;
          default:
            assert(false);
            break;
        }

        GeometricRegion region(x, y, z, dx, dy, dz);
        list.push_back(region);
      } while (iter_z.Advance());
    } while (iter_y.Advance());
  } while (iter_x.Advance());

  map_[key] = list;
}

PartitionHandler::Map& PartitionHandler::map() {
  return map_;
}


