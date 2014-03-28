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
  * Scheduler Version Table. It holds the meta data for logical data version
  * context of each job. It also implements methods to merge the version tables
  * to get new version tables for jobs based on dependencies.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "scheduler/version_table.h"

using namespace nimbus; // NOLINT

VersionTable::VersionTable(version_table_id_t id) {
  id_ = id;
  root_is_set_ = false;
}

VersionTable::~VersionTable() {
}

version_table_id_t VersionTable::id() const {
  return id_;
}

boost::shared_ptr<const VersionTable::Map> VersionTable::root() const {
  return root_;
}

VersionTable::Map VersionTable::content() const {
  return content_;
}

const VersionTable::Map* VersionTable::content_p() const {
  return &content_;
}

bool VersionTable::query_entry(logical_data_id_t l_id, data_version_t *version) const {
  MapConstIter iter;

  iter = content_.find(l_id);
  if (iter != content_.end()) {
    *version = iter->second;
    return true;
  }

  if (root_is_set_) {
    iter = root_->find(l_id);
    if (iter != root_->end()) {
      *version = iter->second;
      return true;
    }
  }
  return false;
}

bool VersionTable::root_is_set() const {
  return root_is_set_;
}

void VersionTable::set_id(version_table_id_t id) {
  id_ = id;
}

void VersionTable::set_root(boost::shared_ptr<const Map> root) {
  root_ = root;
  root_is_set_ = true;
}

void VersionTable::set_content(const VersionTable::Map& content) {
  content_= content;
}

void VersionTable::set_entry(logical_data_id_t l_id, data_version_t version) {
  content_[l_id] = version;
}

void VersionTable::Print() {
  MapConstIter iter;
  std::cout << "Table id: " << id_ << std::endl;
  std::cout << "Root:\n";
  if (root_is_set_) {
    for (iter = root_->begin(); iter != root_->end(); ++iter) {
      std::cout << iter->first << " -> " << iter->second << std::endl;
    }
  }
  std::cout << "Content:\n";
  for (iter = content_.begin(); iter != content_.end(); ++iter) {
    std::cout << iter->first << " -> " << iter->second << std::endl;
  }
}




