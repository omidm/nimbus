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
  * Version Map.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "scheduler/version_map.h"

using namespace nimbus; // NOLINT

VersionMap::VersionMap() {
}

VersionMap::VersionMap(const VersionMap& other) {
  content_ = other.content_;
}

VersionMap::~VersionMap() {
}

bool VersionMap::query_entry(logical_data_id_t l_id, data_version_t *version) const {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  ConstIter iter;

  iter = content_.find(l_id);
  if (iter != content_.end()) {
    *version = iter->second;
    return true;
  }

  return false;
}

void VersionMap::set_entry(logical_data_id_t l_id, data_version_t version) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  content_[l_id] = version;
}

void VersionMap::Print() const {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  ConstIter iter;

  std::cout << "Content:\n";
  for (iter = content_.begin(); iter != content_.end(); ++iter) {
    std::cout << iter->first << " -> " << iter->second << std::endl;
  }
}


VersionMap& VersionMap::operator=(const VersionMap& right) {
  boost::unique_lock<boost::recursive_mutex> lock(mutex_);
  content_ = right.content_;
  return (*this);
}


