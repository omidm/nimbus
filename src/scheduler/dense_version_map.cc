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
  * Dense Version Map class. Instead of using a map, use an array and lookup
  * with the index. Since, the assumption it that the index is dence the memory
  * usage would be a problem. The user of the class has to make sure that the
  * index is dense enough.
  *
  * NOTE: the quey result can never be NIMBUS_UNDEFINED_DATA_VERSION, so it
  * does not work for difference version maps - only absolute version maps. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "src/scheduler/dense_version_map.h"

using namespace nimbus; // NOLINT

DenseVersionMap::DenseVersionMap(logical_data_id_t min_id,
                                 logical_data_id_t max_id) {
  assert(max_id >= min_id);
  min_id_ = min_id;
  max_id_ = max_id;
  list_ = List(max_id - min_id + 1, NIMBUS_UNDEFINED_DATA_VERSION);
}

DenseVersionMap::DenseVersionMap(const DenseVersionMap& other) {
  list_ = other.list_;
  content_ = other.content_;
  min_id_ = other.min_id_;
  max_id_ = other.max_id_;
}

DenseVersionMap::~DenseVersionMap() {
}

bool DenseVersionMap::query_entry(logical_data_id_t l_id, data_version_t *version) const {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  if ((l_id > max_id_) || (l_id < min_id_)) {
    ConstIter iter;

    iter = content_.find(l_id);
    if (iter != content_.end()) {
      *version = iter->second;
      return true;
    }

    return false;
  } else {
    logical_data_id_t idx = l_id - min_id_;
    data_version_t v = list_[idx];
    if (v != NIMBUS_UNDEFINED_DATA_VERSION) {
      *version = v;
      return true;
    }

    return false;
  }
}

void DenseVersionMap::set_entry(logical_data_id_t l_id, data_version_t version) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  if ((l_id > max_id_) || (l_id < min_id_)) {
    content_[l_id] = version;
  } else {
    list_[l_id - min_id_] = version;
  }
}

void DenseVersionMap::Print() const {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);

  std::cout << "Content:\n";
  logical_data_id_t idx = 0;
  for (; idx < list_.size(); ++idx) {
    data_version_t v = list_[idx];
    if (v != NIMBUS_UNDEFINED_DATA_VERSION) {
      std::cout << min_id_+idx << " -> " <<  v  << std::endl;
    }
  }

  ConstIter iter;
  for (iter = content_.begin(); iter != content_.end(); ++iter) {
    std::cout << iter->first << " -> " << iter->second << std::endl;
  }
}


DenseVersionMap& DenseVersionMap::operator=(const DenseVersionMap& right) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  list_ = right.list_;
  content_ = right.content_;
  min_id_ = right.min_id_;
  max_id_ = right.max_id_;
  return (*this);
}


