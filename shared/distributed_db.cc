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
  * Distributed DB class implements a service that is used to save the state
  * required for checkpoint creation in order to provide fault tolerance.
  * 
  * It is basically a key-value store. Upon saving a value with a specific key,
  * the service returns a handle that can be used later to retrieve the value.
  * The handle is overloaded version of the key such that any node in the
  * cluster could find the primary node that saved the data and retrieve it
  * form that node. This way even if a worker is down, other nodes could access
  * the saved state from the non-volatile memory that worker save the data on.
  *
  * Currently each node of the distributed db has a leveldb module to save the
  * value on local disk. The handle interface makes it easy to change the
  * implementation anytime in the future to use other options like AWS S3.
  * For LevelDB implementation handle has the following format:
  * 
  * handle: <node_ip_address>-<leveldb_root>-<key>
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */


#include "shared/distributed_db.h"

using namespace nimbus; // NOLINT

DistributedDB::DistributedDB() {
  initialized_ = false;
}

DistributedDB::~DistributedDB() {
}


void DistributedDB::Initialize(const std::string& ip_address,
                               const worker_id_t& worker_id) {
  return;
}

bool DistributedDB::Put(const std::string& key,
                        const std::string& value,
                        const checkpoint_id_t& checkpoint_id,
                        std::string *handle) {
  return false;
}

bool DistributedDB::Get(const std::string& handle,
                        std::string *value) {
  return false;
}

bool DistributedDB::RemoveCheckpoint(checkpoint_id_t checkpoint_id) {
  return false;
}


bool DistributedDB::RetrieveDBFromOtherNode(const std::string& ip_address,
                                            const std::string& leveldb_root) {
  return false;
}

