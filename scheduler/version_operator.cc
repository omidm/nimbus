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
  * Scheduler Version Operator. It is the main class that performs operations
  * over version tables including merging and making roots out of normal nodes.
  * This module provides caching to expedite the merging operations over known
  * previously calculated merges.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "scheduler/version_operator.h"

#define VERSION_OPERATOR_CACHE_SIZE 10

using namespace nimbus; // NOLINT

VersionOperator::VersionOperator() {
  cache_size_ = VERSION_OPERATOR_CACHE_SIZE;
}

VersionOperator::~VersionOperator() {
}

bool VersionOperator::MergeVersionTables(std::list<boost::shared_ptr<VersionTable> > tables,
    boost::shared_ptr<VersionTable> *result) {
  // TODO(omidm): implement!
  return false;
}

bool VersionOperator::MakeRootTable(boost::shared_ptr<VersionTable> table,
    boost::shared_ptr<VersionTable> *result) {
  // TODO(omidm): implement!
  return false;
}

bool VersionOperator::MakeVersionTableOut(boost::shared_ptr<VersionTable> tabel_in,
    IDSet<logical_data_id_t> write_set,
    boost::shared_ptr<VersionTable> *teble_out) {
  // TODO(omidm): implement!
  return false;
}

version_table_id_t VersionOperator::GetNewVersionTableId() {
  static version_table_id_t id = 0;
  return ++id;
}

