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

VersionTable::VersionTable() {
  is_root_ = false;
}

VersionTable::~VersionTable() {
}

version_table_id_t VersionTable::id() {
  return id_;
}

boost::shared_ptr<VersionTable> VersionTable::root() {
  return root_;
}

VersionTable* VersionTable::root_raw() {
  return get_pointer(root_);
}

VersionTable::Map VersionTable::content() {
  return content_;
}

const VersionTable::Map* VersionTable::content_p() {
  return &content_;
}

bool VersionTable::is_root() {
  return is_root_;
}

bool VersionTable::query_entry(logical_data_id_t l_id, data_version_t *version) {
  // TODO(omidm): implement!
  return false;
}

void VersionTable::set_id(version_table_id_t id) {
  id_ = id;
}

void VersionTable::set_root(boost::shared_ptr<VersionTable> root) {
  root_ = root;
}

void VersionTable::set_root(VersionTable* root) {
  root_ = boost::shared_ptr<VersionTable>(root);
}

void VersionTable::set_content(const VersionTable::Map& content) {
  content_= content;
}

void VersionTable::set_is_root(bool flag) {
  is_root_ = flag;
}

void VersionTable::set_entry(logical_data_id_t l_id, data_version_t version) {
  // TODO(omidm): implement!
}



