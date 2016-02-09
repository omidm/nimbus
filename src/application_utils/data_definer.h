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

#ifndef NIMBUS_SRC_APPLICATION_UTILS_DATA_DEFINER_H_
#define NIMBUS_SRC_APPLICATION_UTILS_DATA_DEFINER_H_

#include <boost/unordered_map.hpp>
#include <sstream> // NOLINT
#include <iostream> // NOLINT
#include <string>
#include <vector>
#include <map>
#include <set>
#include "src/shared/nimbus.h"
#include "src/shared/nimbus_types.h"
#include "src/application_utils/space_iterator.h"



namespace nimbus {

class Job;

class DataDefiner {
 public:
  explicit DataDefiner(Job *job);

  ~DataDefiner();

  bool DefineData(const std::string& name,
                  int_dimension_t size_x,
                  int_dimension_t size_y,
                  int_dimension_t size_z,
                  int_dimension_t bw_x,
                  int_dimension_t bw_y,
                  int_dimension_t bw_z,
                  int_dimension_t part_num_x,
                  int_dimension_t part_num_y,
                  int_dimension_t part_num_z,
                  bool with_global_bw);


  bool DefineScratchData(const std::string& name,
                         int_dimension_t size_x,
                         int_dimension_t size_y,
                         int_dimension_t size_z,
                         int_dimension_t bw_x,
                         int_dimension_t bw_y,
                         int_dimension_t bw_z,
                         int_dimension_t part_num_x,
                         int_dimension_t part_num_y,
                         int_dimension_t part_num_z,
                         bool with_global_bw);

 private:
  typedef boost::unordered_map<std::string, partition_id_t> PartitionMap;

  Job *job_;
  partition_id_t partition_id_;
  PartitionMap part_map_;
};

}  // namespace nimbus

#endif  // NIMBUS_SRC_APPLICATION_UTILS_DATA_DEFINER_H_
